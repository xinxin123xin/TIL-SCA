/* ============================================================================
 * CBD2-256 N/P scaling driver (tsim.h OIM path, independent coprocessor IP)
 *
 * DUT   : sca_pqc_coprocessor_top (independent IP, no CPU / CV-X-IF / SoC)
 * Engine: tsim.h (custom bit-level OIM tensor simulator, NOT Verilator/Icarus)
 *
 * Usage:
 *   cbd2_256_tsim N outdir R [P] [--m3-only] [--seed-offset HEX]
 *
 *   N           : traces per repeat (even), N/2 fixed + N/2 random
 *   outdir      : output directory
 *   R           : repeats (each repeat uses a different random seed space)
 *   P           : events per trace = P (default 256; must be multiple of 8)
 *   --m3-only   : do NOT write the N x P shadow HW matrix; only in-loop Welford
 *   --seed-offset : add HEX to RAN_SEED0 (for bookkeeping)
 *
 * M3 Welford (n, mean, M2) is updated inside the tsim cycle loop on every
 * result_valid event. In --m3-only mode no N x P matrix is allocated or
 * written; only the compact moments (2 groups x P x 3 doubles) and a
 * trace-level audit summary (N rows, not N x P) are retained.
 * ==========================================================================*/
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "oim_bitlvl.h"
#include "tsim.h"

#define MAX_P         256
#define TAINT         0xFFFFFFFFu
#define RESULT_DFF    14
#define GRP_CBD       3
#define CUSTOM3       0x7Bu
#define INSTR(lane)   (((uint32_t)(lane) << 25) | ((uint32_t)GRP_CBD << 12) | CUSTOM3)

typedef struct { double n, mean, m2; } Mom;
static Mom mom[2][MAX_P];

static inline void wadd(Mom *m, double x) {
    m->n += 1.0;
    double d = x - m->mean;
    m->mean += d / m->n;
    m->m2   += d * (x - m->mean);
}
static inline void put(TsimState *st, const char *name, uint64_t v) {
    tsim_write_port(st, port_idx(name), v);
}
static inline uint64_t getp(TsimState *st, const char *name) {
    return tsim_read_port(st, port_idx(name));
}
static inline uint32_t result_reg(const TsimState *st) {
    const OimDff *d = &oim_dffs[RESULT_DFF];
    return (uint32_t)gather_bits(st->li, &oim_bit_indices[d->q_bits_start],
                                 d->q_bits_len, 0, 0);
}
static inline int popcount32(uint32_t x) { return __builtin_popcount(x); }

static uint32_t ref_coeff(uint32_t word, int lane) {
    uint32_t slice = (word >> (4u * lane)) & 0xFu;
    int a = (slice & 1) + ((slice >> 1) & 1);
    int b = ((slice >> 2) & 1) + ((slice >> 3) & 1);
    return (uint32_t)(int32_t)(a - b);
}
static inline uint32_t lcg(uint32_t x) { return x * 1664525u + 1013904223u; }

static int exec_cbd2(TsimState *st, uint32_t instr, uint32_t rs1,
                     uint32_t *out_prev, uint32_t *out_curr) {
    put(st, "x_issue_req_instr_i", instr);
    put(st, "x_issue_req_rs1_i", rs1);
    put(st, "x_issue_req_rs2_i", 0u);
    put(st, "x_issue_valid_i", 1u);
    tsim_cycle(st);
    put(st, "x_issue_valid_i", 0u);
    uint32_t prev = result_reg(st);
    tsim_cycle(st);
    if (!getp(st, "x_result_valid_o")) return -1;
    uint32_t curr = result_reg(st);
    tsim_cycle(st);
    *out_prev = prev;
    *out_curr = curr;
    return 0;
}

/* Repeat seed space stride (different R values get disjoint random pools) */
static const uint32_t REP_STRIDE = 0x9E3779B9u;

int main(int argc, char **argv) {
    int N = 100, R = 1, P = 256;
    const char *outdir = "out";
    int m3_only = 0;
    uint32_t seed_offset = 0;

    /* parse: N outdir R [P] [flags...] */
    if (argc >= 2) N = atoi(argv[1]);
    if (argc >= 3) outdir = argv[2];
    if (argc >= 4) R = atoi(argv[3]);
    if (argc >= 5) P = atoi(argv[4]);
    for (int i = 5; i < argc; i++) {
        if (strcmp(argv[i], "--m3-only") == 0) m3_only = 1;
        else if (strcmp(argv[i], "--seed-offset") == 0 && i+1 < argc) {
            seed_offset = (uint32_t)strtoul(argv[++i], NULL, 16);
        }
    }
    if (N <= 0 || (N % 2) != 0) { fprintf(stderr, "N must be positive even\n"); return 2; }
    if (P <= 0 || P > MAX_P || (P % 8) != 0) { fprintf(stderr, "P must be multiple of 8, <=%d\n", MAX_P); return 2; }
    int Nfix = N/2, Nran = N/2, WORDS = P/8;

    char path[1024];
    const uint32_t FIXED_SEED = 0x12345678u;
    const uint32_t RAN_SEED0  = 0x20260919u + seed_offset;

    /* A0 N x P matrix CSV: only in non-m3-only mode */
    FILE *fa0 = NULL;
    if (!m3_only) {
        snprintf(path, sizeof(path), "%s/a0_shadow_hw.csv", outdir);
        fa0 = fopen(path, "w");
        if (!fa0) { perror("open a0"); return 3; }
        fprintf(fa0, "trace_id,group,event,word,lane,rs1,coeff,shadow_hw\n");
    }

    /* trace-level audit CSV (N rows, NOT N x P) */
    snprintf(path, sizeof(path), "%s/trace_audit.csv", outdir);
    FILE *faud = fopen(path, "w");
    if (!faud) { perror("open audit"); return 3; }
    fprintf(faud, "trace_id,group,trace_mean_hw,trace_max_hw\n");

    memset(mom, 0, sizeof(mom));

    TsimState st;
    int n_events = 0, n_fail = 0, n_valid_check_fail = 0;
    long a0_bytes = 0;  /* logical A0 bytes (N x P uint8) */
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    for (int r = 0; r < R; r++)
    for (int t = 0; t < N; t++) {
        int group = (t < Nfix) ? 0 : 1;
        int trace_id = r * N + t;
        uint32_t words[WORDS];
        if (group == 0) {
            /* fixed group: every repeat uses the same fixed seed (deterministic) */
            uint32_t s = FIXED_SEED;
            for (int i = 0; i < WORDS; i++) { s = lcg(s); words[i] = s; }
        } else {
            /* random group: disjoint pool per repeat */
            uint32_t s = RAN_SEED0 + r * REP_STRIDE + (uint32_t)(t - Nfix) * 2654435761u;
            for (int i = 0; i < WORDS; i++) { s = lcg(s); words[i] = s; }
        }

        tsim_init(&st);
        put(&st, "rst_ni", 0u);
        put(&st, "x_result_ready_i", 1u);
        put(&st, "x_commit_valid_i", 0u);
        put(&st, "x_issue_valid_i", 0u);
        for (int i = 0; i < 4; i++) tsim_cycle(&st);
        put(&st, "rst_ni", 1u);
        tsim_cycle(&st);
        if (getp(&st, "x_result_valid_o")) { fprintf(stderr, "reset stuck\n"); n_fail++; }

        double trace_sum = 0.0; int trace_max = 0;
        for (int w = 0; w < WORDS; w++) {
            for (int lane = 0; lane < 8; lane++) {
                int event = 8*w + lane;
                uint32_t prev, curr;
                if (exec_cbd2(&st, INSTR(lane), words[w], &prev, &curr) != 0) {
                    fprintf(stderr, "trace=%d event=%d no result_valid\n", t, event);
                    n_valid_check_fail++; continue;
                }
                uint32_t expected = ref_coeff(words[w], lane);
                if (curr != expected) {
                    fprintf(stderr, "FUNC FAIL t=%d e=%d got=0x%x exp=0x%x\n", t, event, curr, expected);
                    n_fail++;
                }
                uint32_t shadow = TAINT & (curr ^ prev);
                int hw = popcount32(shadow);
                trace_sum += hw; if (hw > trace_max) trace_max = hw;

                if (fa0) {
                    fprintf(fa0, "%d,%s,%d,%d,%d,0x%08x,0x%08x,%d\n",
                            trace_id, group ? "random":"fixed", event, w, lane,
                            words[w], curr, hw);
                }
                wadd(&mom[group][event], (double)hw);
                n_events++;
            }
        }
        fprintf(faud, "%d,%s,%.4f,%d\n", trace_id, group?"random":"fixed",
                trace_sum / P, trace_max);
    }
    if (fa0) { fclose(fa0); }
    fclose(faud);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double wall = (t1.tv_sec - t0.tv_sec) + 1e-9*(t1.tv_nsec - t0.tv_nsec);

    /* ---- write M3 moments ---- */
    snprintf(path, sizeof(path), "%s/m3_moments.json", outdir);
    FILE *fm = fopen(path, "w");
    if (!fm) { perror("open moments"); return 3; }
    fprintf(fm, "{\n  \"N\": %d, \"R\": %d, \"P\": %d, \"m3_only\": %d,\n", N*R, R, P, m3_only);
    fprintf(fm, "  \"groups\": {\n");
    for (int g = 0; g < 2; g++) {
        fprintf(fm, "    \"%s\": [\n", g?"random":"fixed");
        for (int e = 0; e < P; e++) {
            fprintf(fm, "      {\"event\":%3d,\"n\":%.0f,\"mean\":%.12g,\"m2\":%.12g}%s\n",
                    e, mom[g][e].n, mom[g][e].mean, mom[g][e].m2, (e==P-1)?"":",");
        }
        fprintf(fm, "    ]%s\n", (g==0)?",":"");
    }
    fprintf(fm, "  }\n}\n");
    fclose(fm);

    /* ---- write M3 Welch ---- */
    snprintf(path, sizeof(path), "%s/m3_welch.csv", outdir);
    FILE *fw = fopen(path, "w");
    if (!fw) { perror("open welch"); return 3; }
    fprintf(fw, "event,fixed_n,fixed_mean,random_n,random_mean,welch_t\n");
    double max_abs_t = 0.0;
    for (int e = 0; e < P; e++) {
        double nF=mom[0][e].n, mF=mom[0][e].mean, sF=mom[0][e].m2;
        double nR=mom[1][e].n, mR=mom[1][e].mean, sR=mom[1][e].m2;
        double vF=(nF>1)?sF/(nF-1):0, vR=(nR>1)?sR/(nR-1):0;
        double den = sqrt(vF/nF + vR/nR);
        double tv = (den==0)?0:(mF-mR)/den;
        double at=fabs(tv); if (at>max_abs_t) max_abs_t=at;
        fprintf(fw, "%d,%.0f,%.12g,%.0f,%.12g,%.12g\n", e, nF, mF, nR, mR, tv);
    }
    fclose(fw);

    /* logical A0 size = N*R*P uint8 (shadow HW per event) */
    long a0_logical = (long)N * R * P;
    /* M3 retained: 2 groups * P * 3 doubles = 48P bytes */
    long m3_bytes = 2L * P * 3 * 8;

    fprintf(stderr,
        "CBD2_TSIM N=%d R=%d P=%d events=%d wall=%.2fs m3_only=%d a0_logical=%ld m3_bytes=%ld max|t|=%.4f fail=%d valid_fail=%d\n",
        N*R, R, P, n_events, wall, m3_only, a0_logical, m3_bytes, max_abs_t, n_fail, n_valid_check_fail);
    return (n_fail==0 && n_valid_check_fail==0 && n_events==N*R*P) ? 0 : 1;
}
