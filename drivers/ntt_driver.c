/* ntt_tsim_main.c — NTT-256 streaming-embedding main experiment on tsim.h.
 *
 * Independent WPI coprocessor (no CPU / CV-X-IF / SoC). The testbench drives
 * Custom-3 OP_NTT_BFLY_CT_MUL directly into the issue ports and samples
 * result_reg[31:0] on each result_valid event.
 *
 * Feature: shadow_hw = popcount((curr ^ q_prev) & 0xffffffff); q_prev reset
 *          to 0 at the start of every trace and chained across all 896 events.
 *
 * A0 path (baseline): materialize the N x P shadow matrix in memory, then
 *                    compute Welch t from it after the run.
 * M3 path (streaming): Welford moments (n, mean, M2) are updated IN the tsim
 *                    cycle on every result_valid event; only moments survive.
 *
 * Usage:
 *   ntt_tsim_main <mode> <inputs.hex> <events.txt> <outputs.hex> <n_traces> <out.csv>
 *     mode = func  : functional gate (verify n_traces*896 results + final coeffs)
 *     mode = tvla  : run A0 + in-loop M3, emit per-event Welch-t comparison CSV
 *     mode = m3    : run in-loop M3 only, emit per-event moments and Welch t
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "oim_bitlvl.h"
#include "tsim.h"

/* Fixed CT_MUL instruction word (rd=5, rs1=6, rs2=7; funct3=2; funct7[3:0]=2; opcode 0x7b) */
#define INSTR_CT_MUL (0x07b | (5u<<7) | (2u<<12) | (6u<<15) | (7u<<20) | (2u<<25))

#define EVENTS 896
#define COEFFS 256
#define NCLASS 2  /* 0=fixed(even trace), 1=random(odd trace) */

static int pi_clk, pi_rst, pi_commit, pi_issue_ready, pi_instr, pi_rs1, pi_rs2,
           pi_accept, pi_issue_valid, pi_result_data, pi_result_ready, pi_result_valid;

static TsimState st;
static uint64_t tsim_cycles;
static long protocol_mismatch;

static void cycle_once(void) {
    tsim_cycle(&st);
    tsim_cycles++;
}

static uint32_t read_result_reg(void) {
    uint32_t v = 0;
    const int32_t *q = &oim_bit_indices[oim_dffs[3].q_bits_start];
    for (int i = 0; i < 32; i++) if (q[i] >= 0 && st.li[q[i]]) v |= (1u << i);
    return v;
}
static uint32_t read_port_u(int pi, int w) {
    uint32_t v = 0;
    const int16_t *bits = port_bit_arrays[pi];
    for (int i = 0; i < w; i++) if (bits[i] >= 0 && st.li[bits[i]]) v |= (1u << i);
    return v;
}

/* Issue one CT_MUL butterfly. 3 tsim cycles. Returns result_reg. */
static uint32_t issue_butterfly(uint16_t a, uint16_t b, uint16_t zeta) {
    if (read_port_u(pi_issue_ready, 1) != 1) protocol_mismatch++;
    tsim_write_port(&st, pi_issue_valid, 1);
    tsim_write_port(&st, pi_instr, INSTR_CT_MUL);
    tsim_write_port(&st, pi_rs1, a);
    tsim_write_port(&st, pi_rs2, ((uint32_t)zeta << 16) | b);
    cycle_once();                                      /* IDLE -> ISE_RESULT, operands captured */
    tsim_write_port(&st, pi_issue_valid, 0);
    cycle_once();                                      /* result_reg captured -> RESULT_READY, valid=1 */
    if (read_port_u(pi_result_valid, 1) != 1) protocol_mismatch++;
    uint32_t res = read_result_reg();
    /* sanity: result_data_o must equal result_reg when valid */
    uint32_t dv = read_port_u(pi_result_data, 32);
    if (dv != res) {
        protocol_mismatch++;
        if (protocol_mismatch <= 5)
            fprintf(stderr, "PROTOCOL MISMATCH result_data_o=%08x result_reg=%08x\n", dv, res);
    }
    cycle_once();                                      /* RESULT_READY -> IDLE */
    return res;
}

typedef struct { double n, mean, M2; } Welford;

int main(int argc, char **argv) {
    if (argc < 7) {
        fprintf(stderr, "usage: %s <func|tvla|m3> <inputs.hex> <events.txt> <outputs.hex> <n_traces> <out.csv>\n", argv[0]);
        return 1;
    }
    const char *mode = argv[1];
    const char *inp_path = argv[2];
    const char *evt_path = argv[3];
    const char *out_path = argv[4];
    int n_traces = atoi(argv[5]);
    const char *csv_path = argv[6];
    int mode_func = strcmp(mode, "func") == 0;
    int mode_tvla = strcmp(mode, "tvla") == 0;
    int mode_m3 = strcmp(mode, "m3") == 0;
    if ((!mode_func && !mode_tvla && !mode_m3) || n_traces < 1) {
        fprintf(stderr, "invalid mode or trace count: mode=%s n_traces=%d\n", mode, n_traces);
        return 1;
    }

    pi_clk    = port_idx("clk_i");
    pi_rst    = port_idx("rst_ni");
    pi_commit = port_idx("x_commit_valid_i");
    pi_issue_ready = port_idx("x_issue_ready_o");
    pi_instr  = port_idx("x_issue_req_instr_i");
    pi_rs1    = port_idx("x_issue_req_rs1_i");
    pi_rs2    = port_idx("x_issue_req_rs2_i");
    pi_accept = port_idx("x_issue_resp_accept_o");
    pi_issue_valid = port_idx("x_issue_valid_i");
    pi_result_data = port_idx("x_result_data_o");
    pi_result_ready = port_idx("x_result_ready_i");
    pi_result_valid = port_idx("x_result_valid_o");

    size_t coeff_count = (size_t)n_traces * COEFFS;
    size_t event_count = (size_t)n_traces * EVENTS;

    /* Load input vectors: n_traces*256 uint16 */
    uint16_t *inputs = calloc(coeff_count, sizeof *inputs);
    uint16_t *outputs = calloc(coeff_count, sizeof *outputs);
    int *ev_a = calloc(event_count, sizeof *ev_a);
    int *ev_b = calloc(event_count, sizeof *ev_b);
    uint16_t *ev_zeta = calloc(event_count, sizeof *ev_zeta);
    uint32_t *ev_exp = calloc(event_count, sizeof *ev_exp);
    if (!inputs || !outputs || !ev_a || !ev_b || !ev_zeta || !ev_exp) {
        fprintf(stderr, "allocation failed for n_traces=%d\n", n_traces);
        return 1;
    }
    FILE *fi = fopen(inp_path, "r");
    if (!fi) { perror("inputs"); return 1; }
    size_t nread = 0;
    while (nread < coeff_count && fscanf(fi, "%hx", &inputs[nread]) == 1) nread++;
    fclose(fi);
    fprintf(stderr, "loaded %zu input coeffs (want %zu)\n", nread, coeff_count);
    if (nread != coeff_count) { fprintf(stderr, "INPUT COUNT MISMATCH\n"); return 1; }

    /* Load output vectors (for functional final-array check) */
    FILE *fo = fopen(out_path, "r");
    if (!fo) { perror("outputs"); return 1; }
    nread = 0;
    while (nread < coeff_count && fscanf(fo, "%hx", &outputs[nread]) == 1) nread++;
    fclose(fo);
    if (nread != coeff_count) { fprintf(stderr, "OUTPUT COUNT MISMATCH\n"); return 1; }

    /* Parse schedule into arrays for the requested n_traces.
     * We stream events.txt; it is ordered trace 0..99, 896 each. */
    FILE *fe = fopen(evt_path, "r");
    if (!fe) { perror("events"); return 1; }
    char hdr[512]; if (!fgets(hdr, sizeof hdr, fe)) { fprintf(stderr,"no hdr\n"); return 1; }

    /* Schedule storage */
    size_t total_events = event_count;
    size_t got = 0;
    int tr, gi, stage, len, start, j, zi; unsigned zval; int ai, bi; unsigned exp;
    while (got < total_events &&
           fscanf(fe, "%d %d %d %d %d %d %d %x %d %d %x\n",
                  &tr,&gi,&stage,&len,&start,&j,&zi,&zval,&ai,&bi,&exp) == 11) {
        if (tr >= n_traces) break;
        size_t idx = (size_t)tr * EVENTS + gi;
        ev_a[idx] = ai; ev_b[idx] = bi; ev_zeta[idx] = (uint16_t)zval; ev_exp[idx] = exp;
        got++;
    }
    fclose(fe);
    fprintf(stderr, "parsed %zu scheduled events (want %zu)\n", got, total_events);
    if (got != total_events) { fprintf(stderr, "EVENT COUNT MISMATCH\n"); return 1; }

    /* Reset DUT */
    tsim_init(&st);
    tsim_write_port(&st, pi_rst, 0);
    tsim_write_port(&st, pi_commit, 0);
    tsim_write_port(&st, pi_result_ready, 1);
    tsim_write_port(&st, pi_issue_valid, 0);
    tsim_cycles = 0;
    protocol_mismatch = 0;
    for (int i = 0; i < 5; i++) cycle_once();
    tsim_write_port(&st, pi_rst, 1);
    for (int i = 0; i < 3; i++) cycle_once();

    /* Welford moments M3 (in-loop), per event x class */
    static Welford wm[NCLASS][EVENTS];
    /* A0 is allocated only for the explicit N=100 equivalence mode. */
    int class_capacity = (n_traces + 1) / 2;
    double *a0 = mode_tvla
        ? calloc((size_t)NCLASS * EVENTS * class_capacity, sizeof *a0)
        : NULL;
    static int a0_cnt[NCLASS][EVENTS];
    if (mode_tvla && !a0) { fprintf(stderr, "A0 allocation failed\n"); return 1; }
    memset(wm, 0, sizeof wm);
    memset(a0_cnt, 0, sizeof a0_cnt);

    uint16_t state[COEFFS];
    long total_results = 0, mismatch = 0;
    clock_t t0 = clock();

    for (int tr = 0; tr < n_traces; tr++) {
        int cls = (tr & 1) ? 1 : 0;   /* even=fixed(0), odd=random(1) */
        for (int i = 0; i < COEFFS; i++) state[i] = inputs[tr * COEFFS + i];
        uint32_t q_prev = 0;
        for (int e = 0; e < EVENTS; e++) {
            size_t idx = (size_t)tr * EVENTS + e;
            uint32_t res = issue_butterfly(state[ev_a[idx]], state[ev_b[idx]], ev_zeta[idx]);
            total_results++;
            if (res != ev_exp[idx]) {
                mismatch++;
                if (mismatch <= 5)
                    fprintf(stderr, "FUNC MISMATCH tr=%d e=%d got=%08x exp=%08x\n", tr, e, res, ev_exp[idx]);
            }
            /* shadow feature */
            uint32_t shadow = __builtin_popcount(res ^ q_prev);
            q_prev = res;
            /* M3: in-loop Welford */
            Welford *w = &wm[cls][e];
            w->n += 1.0;
            double delta = (double)shadow - w->mean;
            w->mean += delta / w->n;
            w->M2 += delta * ((double)shadow - w->mean);
            /* A0: materialize only in the explicit equivalence mode. */
            if (mode_tvla) {
                int k = a0_cnt[cls][e]++;
                a0[((size_t)cls * EVENTS + e) * class_capacity + k] = (double)shadow;
            }
            /* update state for next event */
            state[ev_a[idx]] = (uint16_t)(res & 0xffff);
            state[ev_b[idx]] = (uint16_t)(res >> 16);
        }
        /* Final-array checking is mandatory in every mode. */
        for (int i = 0; i < COEFFS; i++) {
            uint16_t want = outputs[(size_t)tr * COEFFS + i];
            if (state[i] != want) {
                mismatch++;
                if (mismatch <= 10)
                    fprintf(stderr, "FINAL MISMATCH tr=%d coeff=%d got=%04x want=%04x\n", tr, i, state[i], want);
            }
        }
    }
    clock_t t1 = clock();
    double t_sim = (double)(t1 - t0) / CLOCKS_PER_SEC;
    fprintf(stderr, "done: traces=%d results=%ld mismatch=%ld protocol_mismatch=%ld cycles=%llu wall_sim=%.3fs\n",
            n_traces, total_results, mismatch, protocol_mismatch,
            (unsigned long long)tsim_cycles, t_sim);

    if (mode_func) {
        int pass = mismatch == 0 && protocol_mismatch == 0 && total_results == (long)total_events;
        if (pass)
            printf("TSIM_NTT_FUNCTIONAL_PASS traces=%d results=%ld cycles=%llu\n",
                   n_traces, total_results, (unsigned long long)tsim_cycles);
        else
            printf("TSIM_NTT_FUNCTIONAL_FAIL traces=%d results=%ld mismatch=%ld protocol_mismatch=%ld cycles=%llu\n",
                   n_traces, total_results, mismatch, protocol_mismatch,
                   (unsigned long long)tsim_cycles);
        return pass ? 0 : 2;
    }

    if (mismatch != 0 || protocol_mismatch != 0 || total_results != (long)total_events) {
        fprintf(stderr, "TSIM_NTT_STATISTICS_FAIL results=%ld mismatch=%ld protocol_mismatch=%ld\n",
                total_results, mismatch, protocol_mismatch);
        return 2;
    }

    FILE *csv = fopen(csv_path, "w");
    if (!csv) { perror("csv"); return 1; }

    if (mode_m3) {
        int leaking = 0;
        double max_abs_t = 0.0;
        int peak_event = -1;
        fprintf(csv, "event,n_fixed,mean_fixed,M2_fixed,n_random,mean_random,M2_random,t,abs_t,threshold_4_5,status\n");
        for (int e = 0; e < EVENTS; e++) {
            Welford *w0 = &wm[0][e], *w1 = &wm[1][e];
            double v0 = w0->M2 / (w0->n - 1.0);
            double v1 = w1->M2 / (w1->n - 1.0);
            double denom = sqrt(v0 / w0->n + v1 / w1->n);
            double t = denom == 0.0
                ? (w0->mean == w1->mean ? 0.0 : copysign(INFINITY, w0->mean - w1->mean))
                : (w0->mean - w1->mean) / denom;
            double at = fabs(t);
            int threshold = at >= 4.5;
            const char *status = isfinite(t) ? "finite" : "zero_variance_unequal";
            if (threshold) leaking++;
            if (at > max_abs_t) { max_abs_t = at; peak_event = e; }
            fprintf(csv, "%d,%.0f,%.17g,%.17g,%.0f,%.17g,%.17g,%.17g,%.17g,%d,%s\n",
                    e, w0->n, w0->mean, w0->M2, w1->n, w1->mean, w1->M2,
                    t, at, threshold, status);
        }
        fclose(csv);
        printf("TSIM_NTT_M3_PASS traces=%d results=%ld cycles=%llu leaking=%d peak_event=%d max_abs_t=%.12g wall_sim=%.3fs\n",
               n_traces, total_results, (unsigned long long)tsim_cycles,
               leaking, peak_event, max_abs_t, t_sim);
        return 0;
    }

    /* TVLA: compute Welch t per event from A0 and from M3 moments. */
    fprintf(csv, "event,t_a0,t_m3,dt,peak_fixed,peak_random\n");
    double max_abs_dt = 0.0;
    double max_abs_t = 0.0; int peak_event = -1;
    for (int e = 0; e < EVENTS; e++) {
        /* A0: classical two-sample from stored samples */
        double m0 = 0, m1 = 0; int n0 = a0_cnt[0][e], n1 = a0_cnt[1][e];
        for (int k = 0; k < n0; k++)
            m0 += a0[((size_t)0 * EVENTS + e) * class_capacity + k];
        m0 /= n0;
        for (int k = 0; k < n1; k++)
            m1 += a0[((size_t)1 * EVENTS + e) * class_capacity + k];
        m1 /= n1;
        double v0 = 0, v1 = 0;
        for (int k = 0; k < n0; k++) { double d = a0[((size_t)0 * EVENTS + e) * class_capacity + k] - m0; v0 += d*d; } v0 /= (n0 - 1);
        for (int k = 0; k < n1; k++) { double d = a0[((size_t)1 * EVENTS + e) * class_capacity + k] - m1; v1 += d*d; } v1 /= (n1 - 1);
        double t_a0 = (m0 - m1) / sqrt(v0/n0 + v1/n1);
        /* M3: from in-loop moments */
        double mm0 = wm[0][e].mean, mm1 = wm[1][e].mean;
        double vv0 = wm[0][e].M2 / (wm[0][e].n - 1.0);
        double vv1 = wm[1][e].M2 / (wm[1][e].n - 1.0);
        double t_m3 = (mm0 - mm1) / sqrt(vv0/n0 + vv1/n1);
        double dt = fabs(t_a0 - t_m3);
        if (dt > max_abs_dt) max_abs_dt = dt;
        double at = fabs(t_m3);
        if (at > max_abs_t) { max_abs_t = at; peak_event = e; }
        fprintf(csv, "%d,%.12g,%.12g,%.3e,%d,%d\n", e, t_a0, t_m3, dt,
                (int)(fabs(t_a0) >= 4.5), (int)(fabs(t_m3) >= 4.5));
    }
    fclose(csv);
    printf("TVLA_DONE n_traces=%d max_abs_dt=%.3e max_abs_t=%.3f peak_event=%d wall_sim=%.3fs\n",
           n_traces, max_abs_dt, max_abs_t, peak_event, t_sim);
    return 0;
}
