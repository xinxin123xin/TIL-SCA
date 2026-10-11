/*
 * dom_and_tsim_driver.c — Genuine OIM/tsim.h DOM-AND module-level streaming TVLA
 * + optional trace-level CSV output + same-DFF HD on four share-term banks.
 *
 * Build: gcc -O2 -o dom_and_tsim_driver dom_and_tsim_driver.c -lm
 * Usage: ./dom_and_tsim_driver [N_per_class] [--m3-only] [--campaign-seed 0|1|2]
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
static inline void derive_seed(int trace_id, uint32_t campaign_seed,
                               uint32_t *a, uint32_t *b, uint32_t *c, uint32_t *d) {
    uint32_t h0=0x811c9dc5u, h1=0x01000193u, h2=0xdeadbeefu, h3=0xcafebabeu;
    char buf[80];
    if (campaign_seed == 0) snprintf(buf, sizeof(buf), "TSIM-DOM-AND|trace=%d", trace_id);
    else snprintf(buf, sizeof(buf), "TSIM-DOM-AND|campaign=%u|trace=%d", campaign_seed, trace_id);
    for (char *p=buf; *p; p++) {
        h0=(h0^(uint8_t)*p)*0x01000193u;
        h1=h1+(uint8_t)*p*2654435761u;
        h2=(h2<<5)^((uint8_t)*p*40503u);
        h3=(h3>>3)+((uint8_t)*p<<7);
    }
    *a=h0^(uint32_t)trace_id; *b=h1^(uint32_t)trace_id*0x9e3779b9u;
    *c=h2^(uint32_t)trace_id*0x85ebca6bu; *d=h3^(uint32_t)trace_id*0xc2b2ae35u;
}

static inline uint64_t hash_word(uint64_t h, uint32_t word) {
    for (int i=0; i<4; i++) {
        h ^= (uint8_t)(word >> (8*i));
        h *= UINT64_C(1099511628211);
    }
    return h;
}

static inline uint32_t read_dff_q32(const TsimState *st, int dff_index) {
    const OimDff *dff = &oim_dffs[dff_index];
    return (uint32_t)gather_bits(st->li, &oim_bit_indices[dff->q_bits_start],
                                 dff->q_bits_len, 0, 0);
}

static inline double sample_var(double sum, double sumsq, double n) {
    if (n < 2.0) return NAN;
    double v = (sumsq - sum*sum/n)/(n-1.0);
    return v < 0.0 ? 0.0 : v;
}

static inline double welch_t_if_defined(double sum_f, double sumsq_f, double nf,
                                         double sum_r, double sumsq_r, double nr) {
    double mf = sum_f/nf, mr = sum_r/nr;
    double se = sqrt(sample_var(sum_f, sumsq_f, nf)/nf +
                     sample_var(sum_r, sumsq_r, nr)/nr);
    return se == 0.0 ? NAN : (mf-mr)/se;
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

#define FIXED_X 0xFFFFFFFFu
#define FIXED_Y 0xFFFFFFFFu

int main(int argc, char **argv) {
    int per_class = 250;
    int m3_only = 0;
    uint32_t campaign_seed = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--m3-only") == 0) m3_only = 1;
        else if (strcmp(argv[i], "--campaign-seed") == 0 && i+1 < argc)
            campaign_seed = (uint32_t)strtoul(argv[++i], NULL, 10);
        else per_class = atoi(argv[i]);
    }
    if (per_class < 2) {
        fprintf(stderr, "N_per_class must be at least 2\n");
        return 2;
    }
    if (N_DFFS < 4 ||
        strcmp(oim_dff_rtl_names[0], "reg_cross10[31:0]") != 0 ||
        strcmp(oim_dff_rtl_names[1], "reg_inner1[31:0]") != 0 ||
        strcmp(oim_dff_rtl_names[2], "reg_cross01[31:0]") != 0 ||
        strcmp(oim_dff_rtl_names[3], "reg_inner0[31:0]") != 0) {
        fprintf(stderr, "Unexpected OIM DFF mapping for sec_and HD observation\n");
        return 2;
    }
    int N = 2 * per_class;
    printf("=== DOM-AND tsim.h (OIM) streaming TVLA + same-DFF HD ===\n");
    printf("N per class=%d  total=%d  m3_only=%d  campaign_seed=%u\n\n",
           per_class, N, m3_only, campaign_seed);

    StreamStats st_fixed, st_random;
    stream_init(&st_fixed); stream_init(&st_random);
    double pool_sx=0, pool_sy=0; int pool_n=0;
    double hd_sum[2]={0,0}, hd_sumsq[2]={0,0};
    int hd_nonzero[2]={0,0};
    uint64_t input_digest[2]={UINT64_C(14695981039346656037),
                              UINT64_C(14695981039346656037)};

    double *off_l0_fixed = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *off_l1_fixed = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *off_l0_rand  = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *off_l1_rand  = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *hd_fixed = m3_only ? NULL : calloc(per_class, sizeof(double));
    double *hd_rand  = m3_only ? NULL : calloc(per_class, sizeof(double));

    int func_pass = 0;
    FILE *trace_csv = m3_only ? NULL : fopen("trace_level_dom_and.csv", "w");
    if (trace_csv) fprintf(trace_csv, "trace_id,group,L0,L1,hd_serial\n");

    for (int trace_id = 0; trace_id < N; trace_id++) {
        int is_fixed = (trace_id % 2 == 0);
        uint32_t sa, sb, sc, sd;
        derive_seed(trace_id, campaign_seed, &sa, &sb, &sc, &sd);
        Rng rng; rng_seed(&rng, sa, sb, sc, sd);

        uint32_t x, y;
        if (is_fixed) { x = FIXED_X; y = FIXED_Y; }
        else { x = rng_next(&rng); y = rng_next(&rng); }

        uint32_t x0 = rng_next(&rng);
        uint32_t y0 = rng_next(&rng);
        uint32_t x1 = x ^ x0;
        uint32_t y1 = y ^ y0;
        uint32_t rnd = rng_next(&rng);
        int group = is_fixed ? 0 : 1;
        input_digest[group] = hash_word(input_digest[group], x);
        input_digest[group] = hash_word(input_digest[group], y);
        input_digest[group] = hash_word(input_digest[group], x0);
        input_digest[group] = hash_word(input_digest[group], y0);
        input_digest[group] = hash_word(input_digest[group], rnd);

        TsimState st;
        tsim_init(&st);

        tsim_write_port(&st, port_idx("rst_ni"), 0);
        tsim_write_port(&st, port_idx("en_i"), 0);
        tsim_cycle(&st);

        tsim_write_port(&st, port_idx("rst_ni"), 1);
        tsim_write_port(&st, port_idx("a0_i"), x0);
        tsim_write_port(&st, port_idx("a1_i"), x1);
        tsim_write_port(&st, port_idx("b0_i"), y0);
        tsim_write_port(&st, port_idx("b1_i"), y1);
        tsim_write_port(&st, port_idx("rnd_i"), rnd);
        tsim_write_port(&st, port_idx("en_i"), 1);
        tsim_cycle(&st);

        /* Snapshot the four actual share-term register banks after cycle 1. */
        uint32_t q_cycle1[4];
        for (int d=0; d<4; d++) q_cycle1[d] = read_dff_q32(&st, d);

        tsim_write_port(&st, port_idx("en_i"), 0);
        tsim_cycle(&st);

        uint32_t c0 = (uint32_t)tsim_read_port(&st, port_idx("c0_o"));
        uint32_t c1 = (uint32_t)tsim_read_port(&st, port_idx("c1_o"));

        /* HD of the same four DFF Q banks across the two sampled cycles. */
        double hd_reg_reuse = 0.0;
        for (int d=0; d<4; d++)
            hd_reg_reuse += (double)popcount32(q_cycle1[d] ^ read_dff_q32(&st, d));
        if ((c0 ^ c1) != (x & y)) {
            fprintf(stderr, "FAIL: trace %d c0^c1=%08x expected=%08x\n",
                    trace_id, c0^c1, x&y);
            return 1;
        }
        func_pass++;

        double L0 = (double)popcount32(c0);
        double L1 = (double)popcount32(c1);
        double hd_serial = hd_reg_reuse;
        hd_sum[group] += hd_serial;
        hd_sumsq[group] += hd_serial*hd_serial;
        hd_nonzero[group] += (hd_serial != 0.0);

        if (trace_csv) {
            fprintf(trace_csv, "%d,%s,%.4f,%.4f,%.4f\n",
                    trace_id, is_fixed?"fixed":"random", L0, L1, hd_serial);
        }

        if (is_fixed) {
            if (!m3_only) {
                off_l0_fixed[trace_id/2] = L0;
                off_l1_fixed[trace_id/2] = L1;
                hd_fixed[trace_id/2] = hd_serial;
            }
            stream_update(&st_fixed, L0, L1);
        } else {
            if (!m3_only) {
                off_l0_rand[trace_id/2] = L0;
                off_l1_rand[trace_id/2] = L1;
                hd_rand[trace_id/2] = hd_serial;
            }
            stream_update(&st_random, L0, L1);
        }
        pool_sx += L0; pool_sy += L1; pool_n++;
    }
    if (trace_csv) fclose(trace_csv);

    printf("Functional: %d/%d passed\n\n", func_pass, N);

    double mu0 = pool_sx / pool_n;
    double mu1 = pool_sy / pool_n;
    printf("Pooled mean L0: %.4f  L1: %.4f\n\n", mu0, mu1);

    double zf_m, zf_v, zr_m, zr_v;
    stream_zstats(&st_fixed, mu0, mu1, &zf_m, &zf_v);
    stream_zstats(&st_random, mu0, mu1, &zr_m, &zr_v);
    double t_stream = welch_t(zf_m, zf_v, st_fixed.n, zr_m, zr_v, st_random.n);

    if (m3_only) {
        double t_s0 = welch_t_if_defined(st_fixed.sx, st_fixed.sx2, st_fixed.n,
                                         st_random.sx, st_random.sx2, st_random.n);
        double t_s1 = welch_t_if_defined(st_fixed.sy, st_fixed.sy2, st_fixed.n,
                                         st_random.sy, st_random.sy2, st_random.n);
        double t_hd = welch_t_if_defined(hd_sum[0], hd_sumsq[0], st_fixed.n,
                                        hd_sum[1], hd_sumsq[1], st_random.n);
        double max_first = fmax(fabs(t_s0), fabs(t_s1));
        char hd_t_buf[40];
        if (isnan(t_hd)) snprintf(hd_t_buf, sizeof(hd_t_buf), "NA_zero_variance");
        else snprintf(hd_t_buf, sizeof(hd_t_buf), "%.12g", t_hd);
        int pass = (func_pass == N);
        char result_name[96];
        snprintf(result_name, sizeof(result_name), "m3_only_sec_and_N%d_seed%u.txt",
                 N, campaign_seed);
        FILE *res = fopen(result_name, "w");
        if (res) {
            fprintf(res,
                    "module=sec_and engine=OIM/tsim.h mode=M3-only campaign_seed=%u "
                    "N=%d per_class=%d first_order_share0_t=%.12g "
                    "first_order_share1_t=%.12g first_order_max_abs_t=%.12g "
                    "second_order_t=%.12g second_order_abs_t=%.12g "
                    "hd_nonzero_fixed=%d hd_nonzero_random=%d hd_t=%s "
                    "input_digest_fixed=%016llx input_digest_random=%016llx "
                    "functional=%d/%d raw_trace_matrix_written=0 "
                    "retained_mixed_moments=18 pass=%d\n",
                    campaign_seed, N, per_class, t_s0, t_s1, max_first,
                    t_stream, fabs(t_stream), hd_nonzero[0], hd_nonzero[1],
                    hd_t_buf,
                    (unsigned long long)input_digest[0],
                    (unsigned long long)input_digest[1], func_pass, N, pass);
            fclose(res);
        }
        printf("=== M3-only second-order centered-product Welch t ===\n");
        printf("Streaming: %.6f\n", t_stream);
        printf("First-order share0/share1 |t|: %.6f / %.6f\n", fabs(t_s0), fabs(t_s1));
        printf("Same-register HD nonzero fixed/random: %d/%d; Welch t: ",
               hd_nonzero[0], hd_nonzero[1]);
        if (isnan(t_hd)) printf("undefined (zero variance)\n");
        else printf("%.6f\n", t_hd);
        printf("Input digests fixed/random: %016llx / %016llx\n",
               (unsigned long long)input_digest[0],
               (unsigned long long)input_digest[1]);
        printf("raw_trace_matrix_written=0  retained_mixed_moments=18\n");
        printf("=== RESULT: %s ===\n", pass ? "PASS" : "FAIL");
        return pass ? 0 : 1;
    }

    int nf=per_class, nr=per_class;
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

    /* Serial HD Welch */
    double hf=0, hr=0;
    for (int i=0;i<nf;i++) hf+=hd_fixed[i];
    for (int i=0;i<nr;i++) hr+=hd_rand[i];
    hf/=nf; hr/=nr;
    double vhf=0, vhr=0;
    for (int i=0;i<nf;i++) vhf+=(hd_fixed[i]-hf)*(hd_fixed[i]-hf);
    for (int i=0;i<nr;i++) vhr+=(hd_rand[i]-hr)*(hd_rand[i]-hr);
    vhf/=nf-1; vhr/=nr-1;
    double t_hd = (vhf/nf + vhr/nr == 0.0) ? NAN : welch_t(hf,vhf,nf,hr,vhr,nr);

    double delta = fabs(t_stream - t_offline);
    printf("=== Second-order centered-product Welch t ===\n");
    printf("Offline:   %.6f\n", t_offline);
    printf("Streaming: %.6f\n", t_stream);
    printf("|delta|: %.2e  Match: %s\n\n", delta, (delta<1e-9)?"YES":"NO");

    printf("=== First-order |t| ===\n");
    printf("share0: %.4f   share1: %.4f\n\n", fabs(t_s0), fabs(t_s1));
    printf("=== Same-register HD (four OIM DFF banks, cycle1->cycle2) ===\n");
    if (isnan(t_hd)) printf("Welch t: undefined (zero variance)\n\n");
    else printf("|t|: %.4f  CROSSED: %s\n\n", fabs(t_hd), fabs(t_hd)>=4.5?"yes":"no");

    int pass = (delta < 1e-9) && (func_pass == N);
    printf("=== RESULT: %s ===\n", pass ? "PASS" : "FAIL");
    printf("Wrote trace_level_dom_and.csv (%d rows)\n", N);

    free(off_l0_fixed); free(off_l1_fixed);
    free(off_l0_rand); free(off_l1_rand);
    free(hd_fixed); free(hd_rand);
    free(zf); free(zr);
    return pass ? 0 : 1;
}
