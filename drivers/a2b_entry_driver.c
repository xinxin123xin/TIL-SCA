/*
 * a2b_fsm_entry_driver.c — a2b_core FSM-entry HD analysis (OIM/tsim.h path).
 *
 * Mirrors the sec_add FSM-entry observation boundary (zero-inclusive state
 * entry HD; unit = trace; Welch t between fixed/random classes; |t|>=4.5
 * screening threshold; t undefined when both group variances are zero).
 *
 * Observation: for each selected working register r, at every a2b_core FSM
 * state ENTRY (state_q change), sample r's value and compute
 *     h_k = HW( r(tau_k) xor r(tau_{k-1}) )
 * where tau_{k-1} is the previous state entry. Per trace, per state, the
 * per-entry HDs are averaged; the trace is the statistical unit.
 *
 * Registers observed (24-bit, exposed by existing dbg_* ports — NO RTL change):
 *   r0 = dbg_raw0_o  (raw0_r, pass-1 sum share 0, written at ADD1_RUN->FILL_RND2)
 *   r1 = dbg_raw1_o  (raw1_r, pass-1 sum share 1, written at ADD1_RUN->FILL_RND2)
 *   r2 = dbg_trial0_o(trial0_r, pass-2 trial share 0, written at ADD2_RUN->SELECT)
 *   r3 = dbg_trial1_o(trial1_r, pass-2 trial share 1, written at ADD2_RUN->SELECT)
 *
 * Inputs / randomness / functional check are byte-for-byte the same
 * conventions as a2b_tsim_driver.c (fixed x=1234, mod Q=3329, wide_mode=0,
 * 408-bit rnd_bypass_i from per-trace RNG, y0^y1 == x check, 5-cycle reset,
 * start pulse, max 200 cycles).
 *
 * Build:
 *   gcc -std=c11 -O2 -Wall -Wextra
 *       -Wno-unused-parameter -Wno-misleading-indentation
 *       -I <a2b_module> a2b_fsm_entry_driver.c -o a2b_fsm_entry_driver.exe -lm
 * Usage: a2b_fsm_entry_driver.exe [per_class] [campaign_seed]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "oim_bitlvl.h"
#include "tsim.h"

#define WIDTH 24
#define MAX_CYCLES 200
#define MAX_STATES 16
#define N_REGS 4

static inline int popcount32(uint32_t x) {
    x = x - ((x >> 1) & 0x55555555u);
    x = (x & 0x33333333u) + ((x >> 2) & 0x33333333u);
    x = (x + (x >> 4)) & 0x0F0F0F0Fu;
    return (x * 0x01010101u) >> 24;
}
static inline int popcount24(uint32_t x) { return popcount32(x & 0xFFFFFF); }

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
    char buf[64]; snprintf(buf, sizeof(buf), "TSIM-A2B|trace=%d", tid);
    for (char *p=buf; *p; p++) {
        h0=(h0^(uint8_t)*p)*0x01000193u;
        h1=h1+(uint8_t)*p*2654435761u;
        h2=(h2<<5)^((uint8_t)*p*40503u);
        h3=(h3>>3)+((uint8_t)*p<<7);
    }
    *a=h0^(uint32_t)tid; *b=h1^(uint32_t)tid*0x9e3779b9u;
    *c=h2^(uint32_t)tid*0x85ebca6bu; *d=h3^(uint32_t)tid*0xc2b2ae35u;
}

static void write_wide(TsimState *st, int pi, uint64_t *words) {
    const int16_t *bits = port_bit_arrays[pi];
    int width = oim_ports[pi].width;
    for (int i = 0; i < width; i++)
        if (bits[i] >= 0) st->li[bits[i]] = (words[i/64] >> (i%64)) & 1;
}

static inline double welch_t(double m1, double v1, double n1,
                             double m2, double v2, double n2) {
    double se = sqrt(v1/n1 + v2/n2);
    if (se == 0) return (m1==m2)?0:copysign(1e30,m1-m2);
    return (m1-m2)/se;
}

/* streaming accumulator (mean/var from running sum) */
typedef struct { double n, sum, sumsq; } RAcc;
static inline void racc_init(RAcc *a) { a->n=0; a->sum=0; a->sumsq=0; }
static inline void racc_update(RAcc *a, double v) { a->n += 1.0; a->sum += v; a->sumsq += v*v; }
static inline double racc_mean(const RAcc *a) { return a->n>0 ? a->sum/a->n : 0; }
static inline double racc_var(const RAcc *a) {
    if (a->n < 2) return 0;
    double m = a->sum/a->n;
    return (a->sumsq/a->n - m*m) * a->n/(a->n-1);
}

static const char *state_name(int s) {
    switch (s) {
        case 0: return "ST_IDLE";
        case 1: return "ST_FILL_RND1";
        case 2: return "ST_ADD1_START(unused)";
        case 3: return "ST_ADD1_RUN";
        case 4: return "ST_FILL_RND2";
        case 5: return "ST_ADD2_START(unused)";
        case 6: return "ST_ADD2_RUN";
        case 7: return "ST_SELECT";
        case 8: return "ST_DONE";
        default: return "ST_UNUSED";
    }
}
static const char *reg_name(int r) {
    switch (r) {
        case 0: return "raw0_r";
        case 1: return "raw1_r";
        case 2: return "trial0_r";
        case 3: return "trial1_r";
        default: return "?";
    }
}
/* predecessor of the transition INTO state s (within a2b linear FSM) */
static const char *transition_pair(int s) {
    switch (s) {
        case 1: return "ST_IDLE->ST_FILL_RND1";
        case 4: return "ST_ADD1_RUN->ST_FILL_RND2";
        case 7: return "ST_ADD2_RUN->ST_SELECT";
        default: return "-";
    }
}

#define Q 3329u
#define W_MASK 0x1FFF

int main(int argc, char **argv) {
    int per_class = 250;
    int positional = 0;
    for (int i = 1; i < argc; i++) {
        if (positional++ == 0) per_class = atoi(argv[i]);
        else g_campaign_seed = atoi(argv[i]);
    }
    if (per_class < 2) {
        fprintf(stderr, "N_per_class must be at least 2\n");
        return 2;
    }
    int N = 2 * per_class;

    printf("=== a2b_core FSM-entry HD (OIM/tsim.h, zero-inclusive) ===\n");
    printf("registers: raw0_r raw1_r trial0_r trial1_r (24-bit each, dbg_* ports)\n");
    printf("per_class=%d total=%d campaign_seed=%d\n\n", per_class, N, g_campaign_seed);

    /* offline per-trace arrays: [state][reg][class] */
    double *hd_fixed[MAX_STATES][N_REGS];
    double *hd_rand [MAX_STATES][N_REGS];
    RAcc stream_fixed[MAX_STATES][N_REGS];
    RAcc stream_rand [MAX_STATES][N_REGS];
    int state_seen[MAX_STATES] = {0};

    for (int s = 0; s < MAX_STATES; s++) {
        for (int r = 0; r < N_REGS; r++) {
            hd_fixed[s][r] = calloc(per_class, sizeof(double));
            hd_rand[s][r]  = calloc(per_class, sizeof(double));
            racc_init(&stream_fixed[s][r]);
            racc_init(&stream_rand[s][r]);
        }
    }

    int func_pass = 0;
    int max_cycles_used = 0;

    int pstart = port_idx("start_i");
    int prst = port_idx("rst_ni");
    int ps0i = port_idx("share0_i");
    int ps1i = port_idx("share1_i");
    int prndbypass = port_idx("rnd_bypass_i");
    int pwm  = port_idx("wide_mode_i");
    int ps0o = port_idx("share0_o");
    int ps1o = port_idx("share1_o");
    int pdone = port_idx("done_o");
    int pstate = port_idx("dbg_state_o");
    int preg[N_REGS];
    preg[0] = port_idx("dbg_raw0_o");
    preg[1] = port_idx("dbg_raw1_o");
    preg[2] = port_idx("dbg_trial0_o");
    preg[3] = port_idx("dbg_trial1_o");

    char trace_csv_name[96];
    snprintf(trace_csv_name, sizeof(trace_csv_name), "trace_state_hd_a2b_N%d_seed%d.csv", N, g_campaign_seed);
    FILE *trace_csv = fopen(trace_csv_name, "w");
    if (trace_csv) {
        fprintf(trace_csv, "trace_id,group,state,state_name,hd_raw0,hd_raw1,hd_tr0,hd_tr1,cycles\n");
    }

    for (int trace_id = 0; trace_id < N; trace_id++) {
        int is_fixed = (trace_id % 2 == 0);
        uint32_t sa, sb, sc, sd;
        derive_seed(trace_id, &sa, &sb, &sc, &sd);
        Rng rng; rng_seed(&rng, sa, sb, sc, sd);

        uint32_t x = is_fixed ? 1234 : (rng_next(&rng) % Q);
        uint32_t s0 = rng_next(&rng) % Q;
        uint32_t s1 = (x - s0 + Q) % Q;

        uint64_t rnd_words[7];
        for (int w = 0; w < 7; w++)
            rnd_words[w] = ((uint64_t)rng_next(&rng) << 32) | rng_next(&rng);

        TsimState st;
        tsim_init(&st);

        tsim_write_port(&st, prst, 0);
        tsim_write_port(&st, pstart, 0);
        for (int i=0;i<5;i++) tsim_cycle(&st);

        tsim_write_port(&st, prst, 1);
        tsim_write_port(&st, ps0i, s0);
        tsim_write_port(&st, ps1i, s1);
        tsim_write_port(&st, pwm, 0);
        write_wide(&st, prndbypass, rnd_words);

        tsim_write_port(&st, pstart, 1);
        tsim_cycle(&st);
        tsim_write_port(&st, pstart, 0);

        /* FSM-entry tracking: per-state per-register HD */
        uint32_t prev_val[N_REGS][MAX_STATES];
        memset(prev_val, 0, sizeof(prev_val));
        double hd_sum[MAX_STATES][N_REGS];
        memset(hd_sum, 0, sizeof(hd_sum));
        int hd_cnt[MAX_STATES];
        memset(hd_cnt, 0, sizeof(hd_cnt));
        int visited[MAX_STATES];
        memset(visited, 0, sizeof(visited));
        int last_state = -1;

        int cycles = 0;
        while (cycles < MAX_CYCLES) {
            tsim_cycle(&st);
            cycles++;

            uint32_t v[N_REGS];
            v[0] = (uint32_t)tsim_read_port(&st, preg[0]) & 0xFFFFFF;
            v[1] = (uint32_t)tsim_read_port(&st, preg[1]) & 0xFFFFFF;
            v[2] = (uint32_t)tsim_read_port(&st, preg[2]) & 0xFFFFFF;
            v[3] = (uint32_t)tsim_read_port(&st, preg[3]) & 0xFFFFFF;
            uint32_t state = (uint32_t)tsim_read_port(&st, pstate) & 0xF;

            /* State entry detected: sample ALL registers, aligned to the entry */
            if ((int)state != last_state && state < MAX_STATES) {
                visited[state] = 1;
                state_seen[state] = 1;
                if (last_state >= 0) {
                    for (int r = 0; r < N_REGS; r++)
                        hd_sum[state][r] += (double)popcount24(v[r] ^ prev_val[r][last_state]);
                    hd_cnt[state]++;
                }
                for (int r = 0; r < N_REGS; r++)
                    prev_val[r][state] = v[r];
                last_state = (int)state;
            }

            if (tsim_read_port(&st, pdone)) break;
        }
        if (cycles > max_cycles_used) max_cycles_used = cycles;

        uint32_t y0 = (uint32_t)tsim_read_port(&st, ps0o) & W_MASK;
        uint32_t y1 = (uint32_t)tsim_read_port(&st, ps1o) & W_MASK;

        if ((y0 ^ y1) != x) {
            fprintf(stderr, "FAIL trace %d: y0^y1=%04x expected=%04x\n",
                    trace_id, y0^y1, x);
            return 1;
        }
        func_pass++;

        /* per-trace per-state per-register mean HD -> offline array + streaming */
        for (int s = 0; s < MAX_STATES; s++) {
            if (!visited[s] || hd_cnt[s] <= 0) continue;
            double avg[N_REGS];
            for (int r = 0; r < N_REGS; r++) avg[r] = hd_sum[s][r] / (double)hd_cnt[s];
            if (is_fixed) {
                for (int r = 0; r < N_REGS; r++) {
                    hd_fixed[s][r][trace_id/2] = avg[r];
                    racc_update(&stream_fixed[s][r], avg[r]);
                }
            } else {
                for (int r = 0; r < N_REGS; r++) {
                    hd_rand[s][r][trace_id/2] = avg[r];
                    racc_update(&stream_rand[s][r], avg[r]);
                }
            }
            if (trace_csv) {
                fprintf(trace_csv, "%d,%s,%d,%s,%.6f,%.6f,%.6f,%.6f,%d\n",
                        trace_id, is_fixed?"fixed":"random", s, state_name(s),
                        avg[0], avg[1], avg[2], avg[3], cycles);
            }
        }
    }
    if (trace_csv) fclose(trace_csv);

    printf("Functional: %d/%d, max cycles=%d\n\n", func_pass, N, max_cycles_used);

    char state_csv_name[96];
    snprintf(state_csv_name, sizeof(state_csv_name), "state_hd_a2b_N%d_seed%d.csv", N, g_campaign_seed);
    FILE *state_csv = fopen(state_csv_name, "w");
    if (state_csv) {
        fprintf(state_csv, "state,state_name,register,width,transition_pair,"
                           "n_fixed,n_random,zero_fixed,zero_random,"
                           "mean_fixed,mean_random,var_fixed,var_random,"
                           "welch_t_offline,welch_t_streaming,delta_t,match\n");
    }

    /* Welch per (state, register); zero-inclusive; t undefined if both vars 0 */
    double max_t[N_REGS] = {0};
    int max_t_state[N_REGS];
    int max_t_defined[N_REGS] = {0};
    double max_delta = 0;
    int any_mismatch = 0;

    printf("=== Per-FSM-state-entry HD Welch |t| (fixed vs random) ===\n");
    for (int s = 0; s < MAX_STATES; s++) {
        if (!state_seen[s]) continue;
        for (int r = 0; r < N_REGS; r++) {
            int nf = per_class, nr = per_class;
            int zf = 0, zr = 0;
            double mf = 0, mr = 0, vf = 0, vr = 0;
            for (int i = 0; i < nf; i++) { mf += hd_fixed[s][r][i]; if (hd_fixed[s][r][i] == 0.0) zf++; }
            for (int i = 0; i < nr; i++) { mr += hd_rand [s][r][i]; if (hd_rand [s][r][i] == 0.0) zr++; }
            mf /= nf; mr /= nr;
            for (int i = 0; i < nf; i++) vf += (hd_fixed[s][r][i]-mf)*(hd_fixed[s][r][i]-mf);
            for (int i = 0; i < nr; i++) vr += (hd_rand [s][r][i]-mr)*(hd_rand [s][r][i]-mr);
            vf /= (nf-1); vr /= (nr-1);

            double t_off = 0, t_str = 0;
            int defined = (vf > 0.0) || (vr > 0.0);
            if (defined) {
                t_off = welch_t(mf, vf, nf, mr, vr, nr);
                t_str = welch_t(racc_mean(&stream_fixed[s][r]), racc_var(&stream_fixed[s][r]), nf,
                                racc_mean(&stream_rand [s][r]), racc_var(&stream_rand [s][r]), nr);
            }
            double delta = fabs(t_off - t_str);
            if (delta > max_delta) max_delta = delta;
            int match = (delta < 1e-9);
            if (!match) any_mismatch = 1;

            if (defined) {
                printf("  state[%d] %-8s (%s)  HD(%-8s) n=%d/%d  zero=%d/%d  "
                       "mean=%.6f/%.6f  |t_off|=%.4f  |t_str|=%.4f  dT=%.2e %s\n",
                       s, state_name(s), transition_pair(s), reg_name(r), nf, nr, zf, zr,
                       mf, mr, fabs(t_off), fabs(t_str), delta, match?"match":"MISMATCH");
                if (fabs(t_off) > max_t[r]) { max_t[r] = fabs(t_off); max_t_state[r] = s; max_t_defined[r] = 1; }
            } else {
                printf("  state[%d] %-8s (%s)  HD(%-8s) n=%d/%d  zero=%d/%d  "
                       "HD identically 0 in both classes -> t undefined\n",
                       s, state_name(s), transition_pair(s), reg_name(r), nf, nr, zf, zr);
            }
            if (state_csv) {
                char t_off_s[32], t_str_s[32];
                snprintf(t_off_s, sizeof(t_off_s), defined ? "%.9f" : "undef", t_off);
                snprintf(t_str_s, sizeof(t_str_s), defined ? "%.9f" : "undef", t_str);
                fprintf(state_csv, "%d,%s,%s,%d,%s,%d,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%s,%s,%.2e,%s\n",
                        s, state_name(s), reg_name(r), WIDTH, transition_pair(s),
                        nf, nr, zf, zr, mf, mr, vf, vr,
                        t_off_s, t_str_s, delta, match ? "match" : "MISMATCH");
            }
        }
    }
    if (state_csv) fclose(state_csv);

    printf("\n=== max|t| per register (over states with defined Welch t) ===\n");
    char summary_name[96];
    snprintf(summary_name, sizeof(summary_name), "a2b_fsm_entry_N%d_seed%d.txt", N, g_campaign_seed);
    FILE *res = fopen(summary_name, "w");
    if (res) {
        fprintf(res, "module=a2b_core engine=OIM/tsim.h observation=FSM-entry-HD "
                     "zero_inclusive=1 N=%d per_class=%d campaign_seed=%d "
                     "functional=%d/%d max_cycles=%d threshold=4.5 "
                     "randomness_interface=rnd_bypass_i\n", N, per_class, g_campaign_seed,
                     func_pass, N, max_cycles_used);
    }
    for (int r = 0; r < N_REGS; r++) {
        const char *crossed = "no";
        double mt = max_t_defined[r] ? max_t[r] : 0.0;
        if (max_t_defined[r] && mt >= 4.5) crossed = "yes";
        if (max_t_defined[r]) {
            printf("  %-8s (width %d)  max|t|=%.4f at state[%d] (%s -> %s)  crossed=%s\n",
                   reg_name(r), WIDTH, mt, max_t_state[r],
                   state_name(max_t_state[r]), transition_pair(max_t_state[r]), crossed);
        } else {
            printf("  %-8s (width %d)  max|t|=undef (HD identically 0 in both classes)  crossed=no\n",
                   reg_name(r), WIDTH);
        }
        if (res) {
            if (max_t_defined[r]) {
                fprintf(res, "register=%s width=%d transition=%s N=%d max_abs_t=%.4f "
                             "threshold_crossed=%s streaming_delta_t=%.2e streaming_match=%s\n",
                        reg_name(r), WIDTH, transition_pair(max_t_state[r]), N, mt,
                        crossed, max_delta, any_mismatch ? "no" : "yes");
            } else {
                fprintf(res, "register=%s width=%d transition=- N=%d max_abs_t=undefined "
                             "threshold_crossed=no streaming_delta_t=%.2e streaming_match=%s\n",
                        reg_name(r), WIDTH, N, max_delta, any_mismatch ? "no" : "yes");
            }
        }
    }
    printf("\nStreaming vs offline: max|delta_t| over all (state,register) = %.2e  Match: %s\n",
           max_delta, any_mismatch ? "NO" : "YES");
    printf("Wrote %s, %s, %s\n", trace_csv_name, state_csv_name, summary_name);

    if (res) {
        fprintf(res, "max_streaming_delta_t=%.2e streaming_match=%s\n",
                max_delta, any_mismatch ? "no" : "yes");
        fclose(res);
    }

    for (int s = 0; s < MAX_STATES; s++) {
        for (int r = 0; r < N_REGS; r++) {
            free(hd_fixed[s][r]);
            free(hd_rand[s][r]);
        }
    }
    int pass = (func_pass == N) && !any_mismatch;
    printf("\n=== RESULT: %s ===\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}
