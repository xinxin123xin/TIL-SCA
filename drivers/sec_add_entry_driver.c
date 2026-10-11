/*
 * sec_add_tsim_driver.c — Genuine OIM/tsim.h sec_add module-level streaming TVLA
 * + zero-inclusive per-FSM-state HD and five-entry state[3] audit.
 *
 * Build: gcc -O2 -o sec_add_tsim_driver sec_add_tsim_driver.c -lm
 * Usage: ./sec_add_tsim_driver [N_per_class] [campaign_seed] [--m3-only]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "oim_bitlvl.h"
#include "tsim.h"

#define WIDTH 13
#define MAX_CYCLES 64
#define MAX_STATES 16

static inline int popcount32(uint32_t x) {
    x = x - ((x >> 1) & 0x55555555u);
    x = (x & 0x33333333u) + ((x >> 2) & 0x33333333u);
    x = (x + (x >> 4)) & 0x0F0F0F0Fu;
    return (x * 0x01010101u) >> 24;
}
static inline uint32_t rotl32(uint32_t x, int k) { return (x << k) | (x >> (32 - k)); }
typedef struct { uint32_t s[4]; } Rng;
static inline uint32_t rng_next(Rng *r) {
    uint32_t result = rotl32(r->s[0] + r->s[3], 7) + r->s[0];
    uint32_t t = r->s[1] << 9;
    r->s[2] ^= r->s[0]; r->s[3] ^= r->s[1];
    r->s[1] ^= r->s[2]; r->s[0] ^= r->s[3];
    r->s[2] ^= t; r->s[3] = rotl32(r->s[3], 11);
    return result;
}
static inline void rng_seed(Rng *r, uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
    r->s[0]=a; r->s[1]=b; r->s[2]=c; r->s[3]=d;
    for (int i=0;i<16;i++) rng_next(r);
}
static int g_campaign_seed = 0;
static inline void derive_seed(int trace_id, uint32_t *a, uint32_t *b, uint32_t *c, uint32_t *d) {
    uint32_t h0=0x811c9dc5u, h1=0x01000193u, h2=0xdeadbeefu, h3=0xcafebabeu;
    int tid = trace_id + g_campaign_seed * 1000000;
    char buf[64]; snprintf(buf, sizeof(buf), "TSIM-SECADD|trace=%d", tid);
    for (char *p=buf; *p; p++) {
        h0=(h0^(uint8_t)*p)*0x01000193u;
        h1=h1+(uint8_t)*p*2654435761u;
        h2=(h2<<5)^((uint8_t)*p*40503u);
        h3=(h3>>3)+((uint8_t)*p<<7);
    }
    *a=h0^(uint32_t)tid; *b=h1^(uint32_t)tid*0x9e3779b9u;
    *c=h2^(uint32_t)tid*0x85ebca6bu; *d=h3^(uint32_t)tid*0xc2b2ae35u;
}
static void write_wide_port(TsimState *st, int pi, const uint8_t *bytes, int nbytes) {
    if (pi < 0) return;
    const OimPort *p = &oim_ports[pi];
    const int16_t *bits = port_bit_arrays[pi];
    for (int i = 0; i < p->width; i++) {
        if (bits[i] < 0) continue;
        int byte_idx = i / 8, bit_idx = i % 8;
        if (byte_idx >= nbytes) break;
        st->li[bits[i]] = (bytes[byte_idx] >> bit_idx) & 1;
    }
}

typedef struct {
    double n, sx, sy, sxy, sx2, sy2, sx2y, sxy2, sx2y2;
} StreamStats;
static inline void stream_init(StreamStats *s) { memset(s, 0, sizeof(*s)); }
static inline void stream_update(StreamStats *s, double x, double y) {
    s->n += 1.0; s->sx += x; s->sy += y; s->sxy += x*y;
    s->sx2 += x*x; s->sy2 += y*y;
    s->sx2y += x*x*y; s->sxy2 += x*y*y; s->sx2y2 += x*x*y*y;
}
static inline void stream_zstats(const StreamStats *s, double mux, double muy,
                                  double *out_mean, double *out_var) {
    double n = s->n;
    if (n < 2.0) { *out_mean=0; *out_var=0; return; }
    double ez  = s->sxy/n - mux*s->sy/n - muy*s->sx/n + mux*muy;
    double ez2 = s->sx2y2/n
               - 2*mux*s->sxy2/n - 2*muy*s->sx2y/n
               + mux*mux*s->sy2/n + 4*mux*muy*s->sxy/n + muy*muy*s->sx2/n
               - 2*mux*mux*muy*s->sy/n - 2*mux*muy*muy*s->sx/n
               + mux*mux*muy*muy;
    double vz = ez2 - ez*ez;
    if (vz < 0) vz = 0;
    *out_mean = ez; *out_var = vz * n / (n - 1.0);
}
static inline double welch_t(double m1, double v1, double n1,
                              double m2, double v2, double n2) {
    double se = sqrt(v1/n1 + v2/n2);
    if (se == 0) return (m1==m2)?0:copysign(1e30,m1-m2);
    return (m1-m2)/se;
}

#ifndef FIXED_A
#define FIXED_A 0x1FFFu
#endif
#ifndef FIXED_B
#define FIXED_B 0x1FFFu
#endif

int main(int argc, char **argv) {
    int per_class = 250;
    int m3_only = 0;
    int positional = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--m3-only") == 0) {
            m3_only = 1;
        } else if (positional++ == 0) {
            per_class = atoi(argv[i]);
        } else {
            g_campaign_seed = atoi(argv[i]);
        }
    }
    if (per_class < 2) {
        fprintf(stderr, "N_per_class must be at least 2\n");
        return 2;
    }
    int N = 2 * per_class;
    printf("=== sec_add tsim.h (OIM) streaming TVLA + per-FSM-state HD ===\n");
    printf("WIDTH=13 N_STAGES=4, fixed_a=0x%04X fixed_b=0x%04X per_class=%d total=%d campaign_seed=%d m3_only=%d\n\n",
           FIXED_A, FIXED_B, per_class, N, g_campaign_seed, m3_only);

    StreamStats st_fixed, st_random;
    stream_init(&st_fixed); stream_init(&st_random);
    double pool_sx=0, pool_sy=0; int pool_n=0;

    double *off_l0_fixed = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *off_l1_fixed = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *off_l0_rand  = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *off_l1_rand  = m3_only ? NULL : calloc(per_class, sizeof(double));

    /* Per-FSM-state HD: for each state, accumulate HD(g0 at entry) per trace */
    double *hd0_state_fixed[MAX_STATES] = {0};
    double *hd0_state_rand[MAX_STATES] = {0};
    double *hd1_state_fixed[MAX_STATES] = {0};
    double *hd1_state_rand[MAX_STATES] = {0};
    uint8_t *state_valid_fixed[MAX_STATES] = {0};
    uint8_t *state_valid_rand[MAX_STATES] = {0};
    for (int s = 0; s < MAX_STATES; s++) {
        if (!m3_only) {
            hd0_state_fixed[s] = calloc(per_class, sizeof(double));
            hd0_state_rand[s]  = calloc(per_class, sizeof(double));
            hd1_state_fixed[s] = calloc(per_class, sizeof(double));
            hd1_state_rand[s]  = calloc(per_class, sizeof(double));
            state_valid_fixed[s] = calloc(per_class, sizeof(uint8_t));
            state_valid_rand[s]  = calloc(per_class, sizeof(uint8_t));
        }
    }

    /* Per-cycle average HD (old method, kept for comparison) */
    double *hd0_avg_fixed = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *hd1_avg_fixed = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *hd0_avg_rand  = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *hd1_avg_rand  = m3_only ? NULL : calloc(per_class, sizeof(double));

    int func_pass = 0;
    int max_cycles_used = 0;

    int pa0 = port_idx("a0_i"), pa1 = port_idx("a1_i");
    int pb0 = port_idx("b0_i"), pb1 = port_idx("b1_i");
    int prnd = port_idx("rnd_i");
    int pstart = port_idx("start_i");
    int prst = port_idx("rst_ni");
    int ps0 = port_idx("s0_o"), ps1 = port_idx("s1_o");
    int pdone = port_idx("done_o");
    int pg0 = port_idx("dbg_g0_o"), pg1 = port_idx("dbg_g1_o");
    int pstate = port_idx("dbg_state_o");
    int pstage = port_idx("dbg_stage_o");

    FILE *trace_csv = m3_only ? NULL : fopen("trace_level_sec_add.csv", "w");
    FILE *state_trace_csv = m3_only ? NULL : fopen("trace_state_hd_sec_add.csv", "w");
    FILE *entry_csv = m3_only ? NULL : fopen("trace_state3_entries_sec_add.csv", "w");
    if (!m3_only && (!trace_csv || !state_trace_csv || !entry_csv)) {
        perror("opening audit CSV");
        return 2;
    }
    if (trace_csv) {
        fprintf(trace_csv, "trace_id,group,L0,L1,hd0_avg,hd1_avg,states_visited,max_state\n");
    }
    if (state_trace_csv) {
        fprintf(state_trace_csv, "trace_id,group,state,hd_g0,hd_g1,transition_count\n");
    }
    if (entry_csv) {
        fprintf(entry_csv, "trace_id,group,entry,from_state,stage,cycle,hd_g0,hd_g1\n");
    }

    for (int trace_id = 0; trace_id < N; trace_id++) {
        int is_fixed = (trace_id % 2 == 0);
        uint32_t sa, sb, sc, sd;
        derive_seed(trace_id, &sa, &sb, &sc, &sd);
        Rng rng; rng_seed(&rng, sa, sb, sc, sd);

        uint32_t a, b;
        if (is_fixed) { a = FIXED_A; b = FIXED_B; }
        else { a = rng_next(&rng) & 0x1FFF; b = rng_next(&rng) & 0x1FFF; }

        uint32_t a0 = rng_next(&rng) & 0x1FFF;
        uint32_t b0 = rng_next(&rng) & 0x1FFF;
        uint32_t a1 = a ^ a0;
        uint32_t b1 = b ^ b0;

        uint8_t rnd_bytes[23];
        for (int i = 0; i < 23; i++) rnd_bytes[i] = (uint8_t)(rng_next(&rng) & 0xFF);

        TsimState st;
        tsim_init(&st);

        tsim_write_port(&st, prst, 0);
        tsim_write_port(&st, pstart, 0);
        tsim_cycle(&st);

        tsim_write_port(&st, prst, 1);
        tsim_write_port(&st, pa0, a0);
        tsim_write_port(&st, pa1, a1);
        tsim_write_port(&st, pb0, b0);
        tsim_write_port(&st, pb1, b1);
        write_wide_port(&st, prnd, rnd_bytes, 23);

        tsim_write_port(&st, pstart, 1);
        tsim_cycle(&st);
        tsim_write_port(&st, pstart, 0);

        /* Track per-FSM-state: record g0/g1 when entering each new state */
        uint32_t prev_g0[MAX_STATES], prev_g1[MAX_STATES];
        int state_hd_count[MAX_STATES] = {0};
        double state_hd0_sum[MAX_STATES] = {0};
        double state_hd1_sum[MAX_STATES] = {0};
        int last_state = -1;
        int state3_entries = 0;

        /* Per-cycle avg HD (old method) */
        uint32_t prev_g0_cycle = 0, prev_g1_cycle = 0;
        double hd0_cycle_sum = 0, hd1_cycle_sum = 0;
        int n_cycle_hd = 0;

        int cycles = 0;
        while (cycles < MAX_CYCLES) {
            tsim_cycle(&st);
            cycles++;

            uint32_t g0 = (uint32_t)tsim_read_port(&st, pg0) & 0x1FFF;
            uint32_t g1 = (uint32_t)tsim_read_port(&st, pg1) & 0x1FFF;
            uint32_t state = (uint32_t)tsim_read_port(&st, pstate) & 0xF;

            /* Per-cycle HD */
            if (cycles > 1) {
                hd0_cycle_sum += (double)popcount32(g0 ^ prev_g0_cycle);
                hd1_cycle_sum += (double)popcount32(g1 ^ prev_g1_cycle);
                n_cycle_hd++;
            }
            prev_g0_cycle = g0;
            prev_g1_cycle = g1;

            /* Per-FSM-state HD: when state changes, compute HD from previous state entry */
            if ((int)state != last_state && state < MAX_STATES) {
                if (last_state >= 0) {
                    int hd0_entry = popcount32(g0 ^ prev_g0[last_state]);
                    int hd1_entry = popcount32(g1 ^ prev_g1[last_state]);
                    state_hd0_sum[state] += (double)hd0_entry;
                    state_hd1_sum[state] += (double)hd1_entry;
                    state_hd_count[state]++;
                    if (state == 3) {
                        state3_entries++;
                        if (entry_csv) {
                            fprintf(entry_csv, "%d,%s,%d,%d,%u,%d,%d,%d\n", trace_id,
                                    is_fixed ? "fixed" : "random", state3_entries,
                                    last_state, (unsigned)tsim_read_port(&st, pstage),
                                    cycles, hd0_entry, hd1_entry);
                        }
                    }
                }
                prev_g0[state] = g0;
                prev_g1[state] = g1;
                last_state = state;
            }

            if (tsim_read_port(&st, pdone)) break;
        }
        if (state3_entries != 5) {
            fprintf(stderr, "FAIL trace %d: state[3] entries=%d, expected 5\n",
                    trace_id, state3_entries);
            return 1;
        }
        if (cycles > max_cycles_used) max_cycles_used = cycles;

        uint32_t s0 = (uint32_t)tsim_read_port(&st, ps0) & 0x1FFF;
        uint32_t s1 = (uint32_t)tsim_read_port(&st, ps1) & 0x1FFF;

        uint32_t expected = (a + b) & 0x1FFF;
        if ((s0 ^ s1) != expected) {
            fprintf(stderr, "FAIL trace %d: s0^s1=%04x expected=%04x\n", trace_id, s0^s1, expected);
            return 1;
        }
        func_pass++;

        double L0 = (double)popcount32(s0);
        double L1 = (double)popcount32(s1);
        double hd0_cycle_avg = (n_cycle_hd > 0) ? hd0_cycle_sum / n_cycle_hd : 0;
        double hd1_cycle_avg = (n_cycle_hd > 0) ? hd1_cycle_sum / n_cycle_hd : 0;

        /* Find max state seen this trace */
        int max_state = 0;
        for (int s = 0; s < MAX_STATES; s++)
            if (state_hd_count[s] > 0) max_state = s;

        if (trace_csv) {
            fprintf(trace_csv, "%d,%s,%.4f,%.4f,%.4f,%.4f,%d,%d\n",
                    trace_id, is_fixed?"fixed":"random", L0, L1,
                    hd0_cycle_avg, hd1_cycle_avg, cycles, max_state);
        }

        if (is_fixed) {
            if (!m3_only) {
                off_l0_fixed[trace_id/2] = L0;
                off_l1_fixed[trace_id/2] = L1;
                hd0_avg_fixed[trace_id/2] = hd0_cycle_avg;
                hd1_avg_fixed[trace_id/2] = hd1_cycle_avg;
            }
            stream_update(&st_fixed, L0, L1);
            for (int s = 0; !m3_only && s < MAX_STATES; s++)
                if (state_hd_count[s] > 0) {
                    hd0_state_fixed[s][trace_id/2] = state_hd0_sum[s] / state_hd_count[s];
                    hd1_state_fixed[s][trace_id/2] = state_hd1_sum[s] / state_hd_count[s];
                    state_valid_fixed[s][trace_id/2] = 1;
                    fprintf(state_trace_csv, "%d,fixed,%d,%.17g,%.17g,%d\n", trace_id, s,
                            hd0_state_fixed[s][trace_id/2], hd1_state_fixed[s][trace_id/2],
                            state_hd_count[s]);
                }
        } else {
            if (!m3_only) {
                off_l0_rand[trace_id/2] = L0;
                off_l1_rand[trace_id/2] = L1;
                hd0_avg_rand[trace_id/2] = hd0_cycle_avg;
                hd1_avg_rand[trace_id/2] = hd1_cycle_avg;
            }
            stream_update(&st_random, L0, L1);
            for (int s = 0; !m3_only && s < MAX_STATES; s++)
                if (state_hd_count[s] > 0) {
                    hd0_state_rand[s][trace_id/2] = state_hd0_sum[s] / state_hd_count[s];
                    hd1_state_rand[s][trace_id/2] = state_hd1_sum[s] / state_hd_count[s];
                    state_valid_rand[s][trace_id/2] = 1;
                    fprintf(state_trace_csv, "%d,random,%d,%.17g,%.17g,%d\n", trace_id, s,
                            hd0_state_rand[s][trace_id/2], hd1_state_rand[s][trace_id/2],
                            state_hd_count[s]);
                }
        }
        pool_sx += L0; pool_sy += L1; pool_n++;
    }
    if (trace_csv) fclose(trace_csv);
    if (state_trace_csv) fclose(state_trace_csv);
    if (entry_csv) fclose(entry_csv);

    printf("Functional: %d/%d passed, max cycles=%d\n\n", func_pass, N, max_cycles_used);

    double mu0 = pool_sx / pool_n;
    double mu1 = pool_sy / pool_n;
    printf("Pooled mean L0=%.4f L1=%.4f\n\n", mu0, mu1);

    double zf_m, zf_v, zr_m, zr_v;
    stream_zstats(&st_fixed, mu0, mu1, &zf_m, &zf_v);
    stream_zstats(&st_random, mu0, mu1, &zr_m, &zr_v);
    double t_stream = welch_t(zf_m, zf_v, st_fixed.n, zr_m, zr_v, st_random.n);

    if (m3_only) {
        int pass = (func_pass == N);
        char result_name[112];
        snprintf(result_name, sizeof(result_name),
                 "m3_only_sec_add_N%d_seed%d.txt", N, g_campaign_seed);
        FILE *res = fopen(result_name, "w");
        if (res) {
            fprintf(res,
                    "module=sec_add engine=OIM/tsim.h mode=M3-only N=%d per_class=%d "
                    "campaign_seed=%d second_order_t=%.12g abs_t=%.12g functional=%d/%d "
                    "raw_trace_matrix_written=0 retained_mixed_moments=18 pass=%d\n",
                    N, per_class, g_campaign_seed, t_stream, fabs(t_stream),
                    func_pass, N, pass);
            fclose(res);
        }
        printf("=== M3-only second-order centered-product Welch t ===\n");
        printf("Streaming: %.6f\n", t_stream);
        printf("raw_trace_matrix_written=0  retained_mixed_moments=18\n");
        printf("=== RESULT: %s ===\n", pass ? "PASS" : "FAIL");
        return pass ? 0 : 1;
    }

    int nf = per_class, nr = per_class;
    double *zf = malloc(nf*sizeof(double));
    double *zr = malloc(nr*sizeof(double));
    for (int i=0;i<nf;i++) zf[i] = (off_l0_fixed[i]-mu0)*(off_l1_fixed[i]-mu1);
    for (int i=0;i<nr;i++) zr[i] = (off_l0_rand[i]-mu0)*(off_l1_rand[i]-mu1);
    double mf=0, mr=0;
    for (int i=0;i<nf;i++) mf+=zf[i]; mf/=nf;
    for (int i=0;i<nr;i++) mr+=zr[i]; mr/=nr;
    double vf=0, vr=0;
    for (int i=0;i<nf;i++) vf+=(zf[i]-mf)*(zf[i]-mf);
    for (int i=0;i<nr;i++) vr+=(zr[i]-mr)*(zr[i]-mr);
    vf/=(nf-1); vr/=(nr-1);
    double t_offline = welch_t(mf, vf, nf, mr, vr, nr);

    /* First-order */
    double s0f=0,s0r=0,s1f=0,s1r=0;
    for (int i=0;i<nf;i++) { s0f+=off_l0_fixed[i]; s1f+=off_l1_fixed[i]; }
    for (int i=0;i<nr;i++) { s0r+=off_l0_rand[i]; s1r+=off_l1_rand[i]; }
    s0f/=nf; s0r/=nr; s1f/=nf; s1r/=nr;
    double vs0f=0,vs0r=0,vs1f=0,vs1r=0;
    for (int i=0;i<nf;i++) { vs0f+=(off_l0_fixed[i]-s0f)*(off_l0_fixed[i]-s0f); vs1f+=(off_l1_fixed[i]-s1f)*(off_l1_fixed[i]-s1f); }
    for (int i=0;i<nr;i++) { vs0r+=(off_l0_rand[i]-s0r)*(off_l0_rand[i]-s0r); vs1r+=(off_l1_rand[i]-s1r)*(off_l1_rand[i]-s1r); }
    vs0f/=nf-1; vs0r/=nr-1; vs1f/=nf-1; vs1r/=nr-1;
    double t_s0 = welch_t(s0f,vs0f,nf,s0r,vs0r,nr);
    double t_s1 = welch_t(s1f,vs1f,nf,s1r,vs1r,nr);

    /* Per-cycle avg HD (old) */
    double hd0f=0,hd0r=0,hd1f=0,hd1r=0;
    for (int i=0;i<nf;i++) { hd0f+=hd0_avg_fixed[i]; hd1f+=hd1_avg_fixed[i]; }
    for (int i=0;i<nr;i++) { hd0r+=hd0_avg_rand[i]; hd1r+=hd1_avg_rand[i]; }
    hd0f/=nf; hd0r/=nr; hd1f/=nf; hd1r/=nr;
    double vhd0f=0,vhd0r=0,vhd1f=0,vhd1r=0;
    for (int i=0;i<nf;i++) { vhd0f+=(hd0_avg_fixed[i]-hd0f)*(hd0_avg_fixed[i]-hd0f); vhd1f+=(hd1_avg_fixed[i]-hd1f)*(hd1_avg_fixed[i]-hd1f); }
    for (int i=0;i<nr;i++) { vhd0r+=(hd0_avg_rand[i]-hd0r)*(hd0_avg_rand[i]-hd0r); vhd1r+=(hd1_avg_rand[i]-hd1r)*(hd1_avg_rand[i]-hd1r); }
    vhd0f/=nf-1; vhd0r/=nr-1; vhd1f/=nf-1; vhd1r/=nr-1;
    double t_hd0_cycle = welch_t(hd0f,vhd0f,nf,hd0r,vhd0r,nr);
    double t_hd1_cycle = welch_t(hd1f,vhd1f,nf,hd1r,vhd1r,nr);

    /* Include zero HD whenever the state was visited; report every state. */
    double max_state_t0 = 0, max_state_t1 = 0;
    int max_state0_id = -1, max_state1_id = -1;
    FILE *state_csv = fopen("state_hd_sec_add.csv", "w");
    if (!state_csv) { perror("opening state summary CSV"); return 2; }
    fprintf(state_csv, "state,channel,n_fixed,n_random,zero_fixed,zero_random,"
                       "mean_fixed,mean_random,var_fixed,var_random,t,abs_t,"
                       "threshold_4_5,status\n");
    printf("=== Per-FSM-state serial HD |t| ===\n");
    for (int s = 0; s < MAX_STATES; s++) {
        for (int channel = 0; channel < 2; channel++) {
            double *fixed = channel ? hd1_state_fixed[s] : hd0_state_fixed[s];
            double *random = channel ? hd1_state_rand[s] : hd0_state_rand[s];
            int count_fixed = 0, count_random = 0, zero_fixed = 0, zero_random = 0;
            double sum_fixed = 0, sum_random = 0;
            for (int i = 0; i < nf; i++) if (state_valid_fixed[s][i]) {
                count_fixed++; sum_fixed += fixed[i]; zero_fixed += fixed[i] == 0.0;
            }
            for (int i = 0; i < nr; i++) if (state_valid_rand[s][i]) {
                count_random++; sum_random += random[i]; zero_random += random[i] == 0.0;
            }

            double mean_fixed = NAN, mean_random = NAN;
            double var_fixed = NAN, var_random = NAN, t = NAN;
            const char *status = count_fixed == 0 && count_random == 0
                               ? "unvisited" : "insufficient";
            if (count_fixed >= 2 && count_random >= 2) {
                mean_fixed = sum_fixed / count_fixed;
                mean_random = sum_random / count_random;
                var_fixed = 0; var_random = 0;
                for (int i = 0; i < nf; i++) if (state_valid_fixed[s][i]) {
                    double d = fixed[i] - mean_fixed; var_fixed += d*d;
                }
                for (int i = 0; i < nr; i++) if (state_valid_rand[s][i]) {
                    double d = random[i] - mean_random; var_random += d*d;
                }
                var_fixed /= count_fixed - 1;
                var_random /= count_random - 1;
                double se = sqrt(var_fixed/count_fixed + var_random/count_random);
                status = se == 0.0 ? "zero_variance" : "ok";
                if (se > 0.0) t = (mean_fixed - mean_random) / se;
            }
            int crossed = isfinite(t) && fabs(t) >= 4.5;
            fprintf(state_csv, "%d,g%d,%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,"
                               "%.17g,%.17g,%d,%s\n", s, channel, count_fixed,
                    count_random, zero_fixed, zero_random, mean_fixed, mean_random,
                    var_fixed, var_random, t, fabs(t), crossed, status);
            if (isfinite(t)) {
                printf("  state[%2d] HD(g%d) |t|=%.3f (n=%d/%d, zeros=%d/%d)\n",
                       s, channel, fabs(t), count_fixed, count_random,
                       zero_fixed, zero_random);
                if (channel == 0 && fabs(t) > max_state_t0) {
                    max_state_t0 = fabs(t); max_state0_id = s;
                }
                if (channel == 1 && fabs(t) > max_state_t1) {
                    max_state_t1 = fabs(t); max_state1_id = s;
                }
            }
        }
    }
    fclose(state_csv);
    printf("  Max HD(g0) state=%d |t|=%.3f\n", max_state0_id, max_state_t0);
    printf("  Max HD(g1) state=%d |t|=%.3f\n\n", max_state1_id, max_state_t1);

    double delta = fabs(t_stream - t_offline);
    printf("=== Second-order centered-product Welch t ===\n");
    printf("Offline:   %.6f\n", t_offline);
    printf("Streaming: %.6f\n", t_stream);
    printf("|delta|: %.2e  Match: %s\n\n", delta, (delta<1e-9)?"YES":"NO");

    printf("=== First-order share-local |t| ===\n");
    printf("s0: %.4f   s1: %.4f\n\n", fabs(t_s0), fabs(t_s1));

    printf("=== Serial HD (per-cycle average) |t| ===\n");
    printf("HD(g0): %.4f   HD(g1): %.4f\n\n", fabs(t_hd0_cycle), fabs(t_hd1_cycle));

    printf("=== Threshold (|t|>=4.5) ===\n");
    printf("2nd-order:          %s\n", (fabs(t_stream)>=4.5)?"CROSSED":"clean");
    printf("1st-order s0:        %s\n", (fabs(t_s0)>=4.5)?"CROSSED":"clean");
    printf("1st-order s1:        %s\n", (fabs(t_s1)>=4.5)?"CROSSED":"clean");
    printf("HD(g0) per-cycle:   %s\n", (fabs(t_hd0_cycle)>=4.5)?"CROSSED":"clean");
    printf("HD(g1) per-cycle:    %s\n", (fabs(t_hd1_cycle)>=4.5)?"CROSSED":"clean");
    printf("HD(g0) per-FSM-state: %s\n", (max_state_t0>=4.5)?"CROSSED":"clean");
    printf("HD(g1) per-FSM-state: %s\n", (max_state_t1>=4.5)?"CROSSED":"clean");

    int pass = (delta < 1e-9) && (func_pass == N);
    printf("\n=== RESULT: %s ===\n", pass ? "PASS" : "FAIL");
    printf("Wrote trace_level_sec_add.csv (%d rows)\n", N);

    for (int s = 0; s < MAX_STATES; s++) {
        free(hd0_state_fixed[s]); free(hd0_state_rand[s]);
        free(hd1_state_fixed[s]); free(hd1_state_rand[s]);
        free(state_valid_fixed[s]); free(state_valid_rand[s]);
    }
    free(off_l0_fixed); free(off_l1_fixed);
    free(off_l0_rand); free(off_l1_rand);
    free(hd0_avg_fixed); free(hd1_avg_fixed);
    free(hd0_avg_rand); free(hd1_avg_rand);
    free(zf); free(zr);
    return pass ? 0 : 1;
}
