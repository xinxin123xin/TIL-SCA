/*
 * a2b_tsim_driver.c — Genuine OIM/tsim.h a2b_core streaming TVLA
 * + serial HD analysis on raw0/raw1/trial0/trial1 working registers.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "oim_bitlvl.h"
#include "tsim.h"

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
typedef struct { double n, sum, sumsq; } MeanVar;
static inline void mv_init(MeanVar *m) { memset(m, 0, sizeof(*m)); }
static inline void mv_update(MeanVar *m, double v) { m->n++; m->sum+=v; m->sumsq+=v*v; }
static inline double mv_mean(MeanVar *m) { return m->n>0 ? m->sum/m->n : 0; }
static inline double mv_var(MeanVar *m) {
    if (m->n < 2) return 0;
    double mean = m->sum/m->n;
    return (m->sumsq/m->n - mean*mean) * m->n/(m->n-1);
}

#define Q 3329u
#define W_MASK 0x1FFF

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
    printf("=== a2b_core tsim.h (OIM) streaming TVLA + serial HD ===\n");
    printf("per_class=%d total=%d campaign_seed=%d m3_only=%d\n\n",
           per_class, N, g_campaign_seed, m3_only);

    StreamStats st_fixed, st_random;
    stream_init(&st_fixed); stream_init(&st_random);
    double pool_sx=0, pool_sy=0; int pool_n=0;

    double *off_l0_fixed = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *off_l1_fixed = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *off_l0_rand  = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *off_l1_rand  = m3_only ? NULL : calloc(per_class, sizeof(double));

    /* ADDED 2026-10-03: streaming 1st-order HW accumulators (n/mean/M2),
       share-local, one per class per share */
    MeanVar mv0_f, mv1_f, mv0_r, mv1_r;
    mv_init(&mv0_f); mv_init(&mv1_f); mv_init(&mv0_r); mv_init(&mv1_r);

    double *hd_raw0_fixed = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *hd_raw0_rand  = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *hd_raw1_fixed = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *hd_raw1_rand  = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *hd_tr0_fixed  = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *hd_tr0_rand   = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *hd_tr1_fixed  = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *hd_tr1_rand   = m3_only ? NULL : calloc(per_class, sizeof(double));

    int func_pass = 0;
    int max_cycles_used = 0;
    FILE *trace_csv = m3_only ? NULL : fopen("trace_level_a2b.csv", "w");
    if (trace_csv) {
        fprintf(trace_csv, "trace_id,group,L0,L1,hd_raw0,hd_raw1,hd_tr0,hd_tr1,cycles\n");
    }

    int pstart = port_idx("start_i");
    int prst = port_idx("rst_ni");
    int ps0i = port_idx("share0_i");
    int ps1i = port_idx("share1_i");
    int prndbypass = port_idx("rnd_bypass_i");
    int pwm  = port_idx("wide_mode_i");
    int ps0o = port_idx("share0_o");
    int ps1o = port_idx("share1_o");
    int pdone = port_idx("done_o");
    int praw0 = port_idx("dbg_raw0_o");
    int praw1 = port_idx("dbg_raw1_o");
    int ptr0  = port_idx("dbg_trial0_o");
    int ptr1  = port_idx("dbg_trial1_o");

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

        uint32_t prev_raw0=0, prev_raw1=0, prev_tr0=0, prev_tr1=0;
        double hd_raw0_sum=0, hd_raw1_sum=0, hd_tr0_sum=0, hd_tr1_sum=0;
        int hd_samples = 0;

        int cycles = 0;
        while (cycles < 200) {
            tsim_cycle(&st);
            cycles++;
            uint32_t cr0 = (uint32_t)tsim_read_port(&st, praw0) & 0xFFFFFF;
            uint32_t cr1 = (uint32_t)tsim_read_port(&st, praw1) & 0xFFFFFF;
            uint32_t ct0 = (uint32_t)tsim_read_port(&st, ptr0) & 0xFFFFFF;
            uint32_t ct1 = (uint32_t)tsim_read_port(&st, ptr1) & 0xFFFFFF;
            if (cycles > 1) {
                hd_raw0_sum += popcount24(cr0 ^ prev_raw0);
                hd_raw1_sum += popcount24(cr1 ^ prev_raw1);
                hd_tr0_sum  += popcount24(ct0 ^ prev_tr0);
                hd_tr1_sum  += popcount24(ct1 ^ prev_tr1);
                hd_samples++;
            }
            prev_raw0=cr0; prev_raw1=cr1; prev_tr0=ct0; prev_tr1=ct1;
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

        double L0 = (double)popcount32(y0);
        double L1 = (double)popcount32(y1);
        double avg_hd_raw0 = hd_samples > 0 ? hd_raw0_sum/hd_samples : 0;
        double avg_hd_raw1 = hd_samples > 0 ? hd_raw1_sum/hd_samples : 0;
        double avg_hd_tr0  = hd_samples > 0 ? hd_tr0_sum/hd_samples : 0;
        double avg_hd_tr1  = hd_samples > 0 ? hd_tr1_sum/hd_samples : 0;

        if (trace_csv) {
            fprintf(trace_csv, "%d,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d\n",
                    trace_id, is_fixed?"fixed":"random", L0, L1,
                    avg_hd_raw0, avg_hd_raw1, avg_hd_tr0, avg_hd_tr1, cycles);
        }

        if (is_fixed) {
            if (!m3_only) {
                off_l0_fixed[trace_id/2] = L0;
                off_l1_fixed[trace_id/2] = L1;
                hd_raw0_fixed[trace_id/2] = avg_hd_raw0;
                hd_raw1_fixed[trace_id/2] = avg_hd_raw1;
                hd_tr0_fixed[trace_id/2]  = avg_hd_tr0;
                hd_tr1_fixed[trace_id/2]  = avg_hd_tr1;
            }
            stream_update(&st_fixed, L0, L1);
            mv_update(&mv0_f, L0); mv_update(&mv1_f, L1); /* ADDED 2026-10-03 */
        } else {
            if (!m3_only) {
                off_l0_rand[trace_id/2] = L0;
                off_l1_rand[trace_id/2] = L1;
                hd_raw0_rand[trace_id/2] = avg_hd_raw0;
                hd_raw1_rand[trace_id/2] = avg_hd_raw1;
                hd_tr0_rand[trace_id/2]  = avg_hd_tr0;
                hd_tr1_rand[trace_id/2]  = avg_hd_tr1;
            }
            stream_update(&st_random, L0, L1);
            mv_update(&mv0_r, L0); mv_update(&mv1_r, L1); /* ADDED 2026-10-03 */
        }
        pool_sx += L0; pool_sy += L1; pool_n++;
    }
    if (trace_csv) fclose(trace_csv);

    printf("Functional: %d/%d, max cycles=%d\n\n", func_pass, N, max_cycles_used);

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
                 "m3_only_a2b_N%d_seed%d.txt", N, g_campaign_seed);
        FILE *res = fopen(result_name, "w");
        if (res) {
            fprintf(res,
                    "module=a2b_core engine=OIM/tsim.h mode=M3-only N=%d per_class=%d "
                    "campaign_seed=%d second_order_t=%.12g abs_t=%.12g functional=%d/%d "
                    "raw_trace_matrix_written=0 retained_mixed_moments=18 "
                    "randomness_interface=rnd_bypass_i pass=%d\n",
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

    double t_hd_raw0, t_hd_raw1, t_hd_tr0, t_hd_tr1;
    {
        MeanVar m_f, m_r;
        mv_init(&m_f); mv_init(&m_r);
        for (int i=0;i<nf;i++) mv_update(&m_f, hd_raw0_fixed[i]);
        for (int i=0;i<nr;i++) mv_update(&m_r, hd_raw0_rand[i]);
        t_hd_raw0 = welch_t(mv_mean(&m_f), mv_var(&m_f), nf, mv_mean(&m_r), mv_var(&m_r), nr);
        mv_init(&m_f); mv_init(&m_r);
        for (int i=0;i<nf;i++) mv_update(&m_f, hd_raw1_fixed[i]);
        for (int i=0;i<nr;i++) mv_update(&m_r, hd_raw1_rand[i]);
        t_hd_raw1 = welch_t(mv_mean(&m_f), mv_var(&m_f), nf, mv_mean(&m_r), mv_var(&m_r), nr);
        mv_init(&m_f); mv_init(&m_r);
        for (int i=0;i<nf;i++) mv_update(&m_f, hd_tr0_fixed[i]);
        for (int i=0;i<nr;i++) mv_update(&m_r, hd_tr0_rand[i]);
        t_hd_tr0 = welch_t(mv_mean(&m_f), mv_var(&m_f), nf, mv_mean(&m_r), mv_var(&m_r), nr);
        mv_init(&m_f); mv_init(&m_r);
        for (int i=0;i<nf;i++) mv_update(&m_f, hd_tr1_fixed[i]);
        for (int i=0;i<nr;i++) mv_update(&m_r, hd_tr1_rand[i]);
        t_hd_tr1 = welch_t(mv_mean(&m_f), mv_var(&m_f), nf, mv_mean(&m_r), mv_var(&m_r), nr);
    }

    double delta = fabs(t_stream - t_offline);
    printf("=== Second-order centered-product Welch t ===\n");
    printf("Offline:   %.6f\n", t_offline);
    printf("Streaming: %.6f\n", t_stream);
    printf("|delta|: %.2e  Match: %s\n\n", delta, (delta<1e-9)?"YES":"NO");

    printf("=== First-order share-local |t| ===\n");
    printf("y0: %.4f   y1: %.4f\n\n", fabs(t_s0), fabs(t_s1));

    /* ADDED 2026-10-03: streaming vs offline 1st-order parity + CSV rows */
    double t_s0_stream = welch_t(mv_mean(&mv0_f), mv_var(&mv0_f), mv0_f.n,
                                 mv_mean(&mv0_r), mv_var(&mv0_r), mv0_r.n);
    double t_s1_stream = welch_t(mv_mean(&mv1_f), mv_var(&mv1_f), mv1_f.n,
                                 mv_mean(&mv1_r), mv_var(&mv1_r), mv1_r.n);
    double d_s0 = fabs(t_s0_stream - t_s0);
    double d_s1 = fabs(t_s1_stream - t_s1);
    printf("Streaming 1st-order y0: %.4f  y1: %.4f  |d|(y0)=%.2e |d|(y1)=%.2e\n",
           fabs(t_s0_stream), fabs(t_s1_stream), d_s0, d_s1);
    {
        FILE *fo_csv = fopen("a2b_first_order.csv", "a");
        if (fo_csv) {
            fseek(fo_csv, 0, SEEK_END);
            if (ftell(fo_csv) == 0)
                fprintf(fo_csv, "module,seed,share,n,max|t|,threshold_crossed,delta_t\n");
            fprintf(fo_csv, "a2b_core,%d,y0,%d,%.4f,%s,%.2e\n",
                    g_campaign_seed, nf, fabs(t_s0), fabs(t_s0)>=4.5?"yes":"no", d_s0);
            fprintf(fo_csv, "a2b_core,%d,y1,%d,%.4f,%s,%.2e\n",
                    g_campaign_seed, nf, fabs(t_s1), fabs(t_s1)>=4.5?"yes":"no", d_s1);
            fclose(fo_csv);
        }
    }

    printf("=== Serial HD |t| (avg per-trace, across FSM cycles) ===\n");
    printf("HD(raw0):  %.4f\n", fabs(t_hd_raw0));
    printf("HD(raw1):  %.4f\n", fabs(t_hd_raw1));
    printf("HD(tr0):   %.4f\n", fabs(t_hd_tr0));
    printf("HD(tr1):   %.4f\n\n", fabs(t_hd_tr1));

    printf("=== Threshold (|t|>=4.5) ===\n");
    printf("2nd-order:      %s\n", (fabs(t_stream)>=4.5)?"CROSSED":"clean");
    printf("1st-order y0:   %s\n", (fabs(t_s0)>=4.5)?"CROSSED":"clean");
    printf("1st-order y1:   %s\n", (fabs(t_s1)>=4.5)?"CROSSED":"clean");
    printf("HD(raw0):       %s\n", (fabs(t_hd_raw0)>=4.5)?"CROSSED":"clean");
    printf("HD(raw1):       %s\n", (fabs(t_hd_raw1)>=4.5)?"CROSSED":"clean");
    printf("HD(trial0):     %s\n", (fabs(t_hd_tr0)>=4.5)?"CROSSED":"clean");
    printf("HD(trial1):     %s\n", (fabs(t_hd_tr1)>=4.5)?"CROSSED":"clean");

    int pass = (delta < 1e-9) && (func_pass == N);
    printf("\n=== RESULT: %s ===\n", pass ? "PASS" : "FAIL");

    FILE *csv = fopen("serial_hd_a2b.csv", "w");
    fprintf(csv, "module,register,event_definition,N,max_abs_t,threshold_crossed\n");
    fprintf(csv, "a2b,raw0_r,HD(raw0_r) across FSM cycles,%d,%.4f,%s\n", N,
            fabs(t_hd_raw0), fabs(t_hd_raw0)>=4.5?"yes":"no");
    fprintf(csv, "a2b,raw1_r,HD(raw1_r) across FSM cycles,%d,%.4f,%s\n", N,
            fabs(t_hd_raw1), fabs(t_hd_raw1)>=4.5?"yes":"no");
    fprintf(csv, "a2b,trial0_r,HD(trial0_r) across FSM cycles,%d,%.4f,%s\n", N,
            fabs(t_hd_tr0), fabs(t_hd_tr0)>=4.5?"yes":"no");
    fprintf(csv, "a2b,trial1_r,HD(trial1_r) across FSM cycles,%d,%.4f,%s\n", N,
            fabs(t_hd_tr1), fabs(t_hd_tr1)>=4.5?"yes":"no");
    fclose(csv);
    printf("\nWrote serial_hd_a2b.csv\n");

    free(off_l0_fixed); free(off_l1_fixed);
    free(off_l0_rand); free(off_l1_rand);
    free(hd_raw0_fixed); free(hd_raw0_rand);
    free(hd_raw1_fixed); free(hd_raw1_rand);
    free(hd_tr0_fixed); free(hd_tr0_rand);
    free(hd_tr1_fixed); free(hd_tr1_rand);
    free(zf); free(zr);
    return pass ? 0 : 1;
}
