/*
 * sec_add_tsim_driver.c — Genuine OIM/tsim.h sec_add module-level streaming TVLA
 * + optional trace-level CSV + per-FSM-state serial HD.
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

/* ADDED 2026-10-03: streaming 1st-order HW accumulators (n/mean/M2) */
typedef struct { double n, sum, sumsq; } MeanVar;
static inline void mv_init(MeanVar *m) { memset(m, 0, sizeof(*m)); }
static inline void mv_update(MeanVar *m, double v) { m->n++; m->sum+=v; m->sumsq+=v*v; }
static inline double mv_mean(MeanVar *m) { return m->n>0 ? m->sum/m->n : 0; }
static inline double mv_var(MeanVar *m) {
    if (m->n < 2) return 0;
    double mean = m->sum/m->n;
    return (m->sumsq/m->n - mean*mean) * m->n/(m->n-1);
}

#define FIXED_A 0x1FFFu
#define FIXED_B 0x1FFFu

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
    printf("WIDTH=13 N_STAGES=4, per_class=%d total=%d campaign_seed=%d m3_only=%d\n\n",
           per_class, N, g_campaign_seed, m3_only);

    StreamStats st_fixed, st_random;
    stream_init(&st_fixed); stream_init(&st_random);
    double pool_sx=0, pool_sy=0; int pool_n=0;

    double *off_l0_fixed = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *off_l1_fixed = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *off_l0_rand  = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *off_l1_rand  = m3_only ? NULL : calloc(per_class, sizeof(double));

    /* ADDED 2026-10-03: streaming 1st-order HW accumulators, share-local */
    MeanVar mv0_f, mv1_f, mv0_r, mv1_r;
    mv_init(&mv0_f); mv_init(&mv1_f); mv_init(&mv0_r); mv_init(&mv1_r);

    /* Per-FSM-state HD: for each state, accumulate HD(g0 at entry) per trace */
    double *hd0_state_fixed[MAX_STATES] = {0};
    double *hd0_state_rand[MAX_STATES] = {0};
    double *hd1_state_fixed[MAX_STATES] = {0};
    double *hd1_state_rand[MAX_STATES] = {0};
    for (int s = 0; s < MAX_STATES; s++) {
        if (!m3_only) {
            hd0_state_fixed[s] = calloc(per_class, sizeof(double));
            hd0_state_rand[s]  = calloc(per_class, sizeof(double));
            hd1_state_fixed[s] = calloc(per_class, sizeof(double));
            hd1_state_rand[s]  = calloc(per_class, sizeof(double));
        }
    }
    int state_seen[MAX_STATES] = {0};

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

    FILE *trace_csv = m3_only ? NULL : fopen("trace_level_sec_add.csv", "w");
    if (trace_csv) {
        fprintf(trace_csv, "trace_id,group,L0,L1,hd0_avg,hd1_avg,states_visited,max_state\n");
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
                state_seen[state] = 1;
                if (last_state >= 0) {
                    state_hd0_sum[state] += (double)popcount32(g0 ^ prev_g0[last_state]);
                    state_hd1_sum[state] += (double)popcount32(g1 ^ prev_g1[last_state]);
                    state_hd_count[state]++;
                }
                prev_g0[state] = g0;
                prev_g1[state] = g1;
                last_state = state;
            }

            if (tsim_read_port(&st, pdone)) break;
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
            mv_update(&mv0_f, L0); mv_update(&mv1_f, L1); /* ADDED 2026-10-03 */
            for (int s = 0; !m3_only && s < MAX_STATES; s++)
                if (state_hd_count[s] > 0) {
                    hd0_state_fixed[s][trace_id/2] = state_hd0_sum[s] / state_hd_count[s];
                    hd1_state_fixed[s][trace_id/2] = state_hd1_sum[s] / state_hd_count[s];
                }
        } else {
            if (!m3_only) {
                off_l0_rand[trace_id/2] = L0;
                off_l1_rand[trace_id/2] = L1;
                hd0_avg_rand[trace_id/2] = hd0_cycle_avg;
                hd1_avg_rand[trace_id/2] = hd1_cycle_avg;
            }
            stream_update(&st_random, L0, L1);
            mv_update(&mv0_r, L0); mv_update(&mv1_r, L1); /* ADDED 2026-10-03 */
            for (int s = 0; !m3_only && s < MAX_STATES; s++)
                if (state_hd_count[s] > 0) {
                    hd0_state_rand[s][trace_id/2] = state_hd0_sum[s] / state_hd_count[s];
                    hd1_state_rand[s][trace_id/2] = state_hd1_sum[s] / state_hd_count[s];
                }
        }
        pool_sx += L0; pool_sy += L1; pool_n++;
    }
    if (trace_csv) fclose(trace_csv);

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

    /* Per-FSM-state HD: report max |t| across states */
    double max_state_t0 = 0, max_state_t1 = 0;
    int max_state0_id = 0, max_state1_id = 0;
    printf("=== Per-FSM-state serial HD |t| ===\n");
    for (int s = 0; s < MAX_STATES; s++) {
        if (!state_seen[s]) continue;
        /* Count non-zero entries */
        int nf_s = 0, nr_s = 0;
        for (int i = 0; i < nf; i++) if (hd0_state_fixed[s][i] > 0) nf_s++;
        for (int i = 0; i < nr; i++) if (hd0_state_rand[s][i] > 0) nr_s++;
        if (nf_s < 5 || nr_s < 5) continue;

        double m0f=0, m0r=0, m1f=0, m1r=0;
        int cf0=0, cr0=0, cf1=0, cr1=0;
        for (int i = 0; i < nf; i++) {
            if (hd0_state_fixed[s][i] > 0) { m0f += hd0_state_fixed[s][i]; cf0++; }
            if (hd1_state_fixed[s][i] > 0) { m1f += hd1_state_fixed[s][i]; cf1++; }
        }
        for (int i = 0; i < nr; i++) {
            if (hd0_state_rand[s][i] > 0) { m0r += hd0_state_rand[s][i]; cr0++; }
            if (hd1_state_rand[s][i] > 0) { m1r += hd1_state_rand[s][i]; cr1++; }
        }
        if (cf0 < 5 || cr0 < 5) continue;
        m0f /= cf0; m0r /= cr0; m1f /= cf1; m1r /= cr1;

        double v0f=0, v0r=0, v1f=0, v1r=0;
        for (int i = 0; i < nf; i++) {
            if (hd0_state_fixed[s][i] > 0) v0f += (hd0_state_fixed[s][i]-m0f)*(hd0_state_fixed[s][i]-m0f);
            if (hd1_state_fixed[s][i] > 0) v1f += (hd1_state_fixed[s][i]-m1f)*(hd1_state_fixed[s][i]-m1f);
        }
        for (int i = 0; i < nr; i++) {
            if (hd0_state_rand[s][i] > 0) v0r += (hd0_state_rand[s][i]-m0r)*(hd0_state_rand[s][i]-m0r);
            if (hd1_state_rand[s][i] > 0) v1r += (hd1_state_rand[s][i]-m1r)*(hd1_state_rand[s][i]-m1r);
        }
        v0f /= (cf0-1); v0r /= (cr0-1); v1f /= (cf1-1); v1r /= (cr1-1);
        double t0 = welch_t(m0f, v0f, cf0, m0r, v0r, cr0);
        double t1 = welch_t(m1f, v1f, cf1, m1r, v1r, cr1);
        printf("  state[%2d]: HD(g0) |t|=%.3f (n=%d/%d)   HD(g1) |t|=%.3f (n=%d/%d)\n",
               s, fabs(t0), cf0, cr0, fabs(t1), cf1, cr1);
        if (fabs(t0) > max_state_t0) { max_state_t0 = fabs(t0); max_state0_id = s; }
        if (fabs(t1) > max_state_t1) { max_state_t1 = fabs(t1); max_state1_id = s; }
    }
    printf("  Max HD(g0) state=%d |t|=%.3f\n", max_state0_id, max_state_t0);
    printf("  Max HD(g1) state=%d |t|=%.3f\n\n", max_state1_id, max_state_t1);

    double delta = fabs(t_stream - t_offline);
    printf("=== Second-order centered-product Welch t ===\n");
    printf("Offline:   %.6f\n", t_offline);
    printf("Streaming: %.6f\n", t_stream);
    printf("|delta|: %.2e  Match: %s\n\n", delta, (delta<1e-9)?"YES":"NO");

    printf("=== First-order share-local |t| ===\n");
    printf("s0: %.4f   s1: %.4f\n\n", fabs(t_s0), fabs(t_s1));

    /* ADDED 2026-10-03: streaming vs offline 1st-order parity + CSV rows */
    double t_s0_stream = welch_t(mv_mean(&mv0_f), mv_var(&mv0_f), mv0_f.n,
                                 mv_mean(&mv0_r), mv_var(&mv0_r), mv0_r.n);
    double t_s1_stream = welch_t(mv_mean(&mv1_f), mv_var(&mv1_f), mv1_f.n,
                                 mv_mean(&mv1_r), mv_var(&mv1_r), mv1_r.n);
    double d_s0 = fabs(t_s0_stream - t_s0);
    double d_s1 = fabs(t_s1_stream - t_s1);
    printf("Streaming 1st-order s0: %.4f  s1: %.4f  |d|(s0)=%.2e |d|(s1)=%.2e\n",
           fabs(t_s0_stream), fabs(t_s1_stream), d_s0, d_s1);
    {
        FILE *fo_csv = fopen("sec_add_first_order.csv", "a");
        if (fo_csv) {
            fseek(fo_csv, 0, SEEK_END);
            if (ftell(fo_csv) == 0)
                fprintf(fo_csv, "module,seed,share,n,max|t|,threshold_crossed,delta_t\n");
            fprintf(fo_csv, "sec_add,%d,s0,%d,%.4f,%s,%.2e\n",
                    g_campaign_seed, nf, fabs(t_s0), fabs(t_s0)>=4.5?"yes":"no", d_s0);
            fprintf(fo_csv, "sec_add,%d,s1,%d,%.4f,%s,%.2e\n",
                    g_campaign_seed, nf, fabs(t_s1), fabs(t_s1)>=4.5?"yes":"no", d_s1);
            fclose(fo_csv);
        }
    }

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
    }
    free(off_l0_fixed); free(off_l1_fixed);
    free(off_l0_rand); free(off_l1_rand);
    free(hd0_avg_fixed); free(hd1_avg_fixed);
    free(hd0_avg_rand); free(hd1_avg_rand);
    free(zf); free(zr);
    return pass ? 0 : 1;
}
