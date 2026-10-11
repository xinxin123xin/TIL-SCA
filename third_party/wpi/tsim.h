/*
 * tsim.h — Genuine RTeAAL tensor simulator for H_v1 PQC coprocessor
 * Bit-level LI buffer: one uint8_t per net-bit (7536 entries)
 * Per-cycle: evaluate leveled layers → DFF clock edge
 */
#ifndef TSIM_H
#define TSIM_H

/* OIM header: included by caller if OIM_BITLVL_H is already defined,
 * otherwise fall back to the canonical pqc_only relative path. */
#include <stdint.h>
#ifndef OIM_BITLVL_H
  #include "../../phase1_extract/oim_bitlvl.h"
#endif
#include <string.h>
#include <stdio.h>

/* Toggle counting: per-cell-type, per-FF (cca.pc), per-bit, and per-cell-bit-vector */
typedef struct {
    uint64_t total;                   /* total bit flips this cycle */
    uint64_t per_type[N_CELL_TYPES];  /* bit flips attributed to each cell type */
    uint64_t per_dff[N_DFFS];         /* per-FF Q-bit flip COUNT this cycle (cca.pc) */
    uint32_t per_dff_qmask[N_DFFS];   /* bit-vector: bit p set if Q[p] flipped this cycle */
    uint64_t cum_total;               /* cumulative total across all cycles */
    uint64_t cum_per_type[N_CELL_TYPES];
    uint64_t cum_per_dff[N_DFFS];     /* cumulative per-FF */
} TsimToggleStats;

typedef struct {
    uint8_t  li[LI_SIZE];
    uint8_t  prev_li[LI_SIZE];   /* previous cycle LI for toggle counting */
    uint64_t cycle;
    /* Toggle counting */
    int      toggle_enable;       /* 0=disabled, 1=enabled */
    TsimToggleStats toggle;
    /* Activity-aware eval */
    uint8_t  changed[LI_SIZE];   /* 1 if bit changed this cycle */
    int      activity_enable;     /* 0=disabled, 1=enabled */
    uint64_t activity_skipped;    /* ops skipped due to no-change inputs */
    uint64_t activity_evaluated;  /* ops actually evaluated */
    /* VCD output */
    FILE    *vcd_fp;
    int      vcd_enable;
    /* Combinational-settling stats (RTeAAL FIX 2026-05-16: iterate-to-fixpoint) */
    int      comb_passes_last_cycle;  /* passes used to converge in last cycle */
    int      comb_passes_max;         /* max passes seen across all cycles */
} TsimState;

/* =====================================================================
 * Memory model — $mem_v2 cells kept abstract (not lowered to FFs by Yosys)
 *
 * Storage lives in TsimMem[N_MEMS] side-table indexed by mem cell index.
 * Each row is stored LSB-first across ceil(width/8) bytes. tsim_mem_cycle()
 * is called every clock edge AFTER comb eval and BEFORE the DFF update —
 * write happens on the rising edge (clocked), reads drive RD_DATA bits in
 * LI for the next-cycle comb eval. This is the simplest sane model for
 * single-port-read + write-first $mem_v2 instances that yosys emits for
 * cv32e40x's ROM/SRAM.  More sophisticated multi-port + write-mask handling
 * follows the OimMem fields exactly.
 *
 * Gated behind #ifdef N_MEMS so pqc_only flow with no memories is unaffected.
 * ===================================================================== */
#ifdef N_MEMS
#include <stdlib.h>
typedef struct {
    uint8_t *data;
    uint32_t depth;
    uint16_t width_bits;
    uint16_t row_bytes;
} TsimMem;
static TsimMem _tsim_mems[N_MEMS > 0 ? N_MEMS : 1];
static int _tsim_mems_inited = 0;

static inline void tsim_mems_init_storage(void) {
    if (_tsim_mems_inited) return;
    for (int m = 0; m < N_MEMS; m++) {
        const OimMem *om = &oim_mems[m];
        _tsim_mems[m].depth      = om->size;
        _tsim_mems[m].width_bits = om->width;
        _tsim_mems[m].row_bytes  = (om->width + 7) / 8;
        size_t bytes = (size_t)om->size * _tsim_mems[m].row_bytes;
        _tsim_mems[m].data = (uint8_t*)calloc(1, bytes);
        if (om->init) memcpy(_tsim_mems[m].data, om->init, bytes);
    }
    _tsim_mems_inited = 1;
}

/* Load a .hex file (one 32-bit word per line) into a specific memory.
 * Used to override the OIM's zero-init with the actual firmware image
 * (Yosys captures $readmemh init only when it sees synthesizable initial
 * blocks; the CW305 top uses sim-time $readmemh which Yosys ignores). */
static inline int tsim_mem_load_hex(int mem_idx, const char *path) {
    if (mem_idx < 0 || mem_idx >= N_MEMS) return -1;
    if (!_tsim_mems_inited) tsim_mems_init_storage();
    TsimMem *m = &_tsim_mems[mem_idx];
    FILE *fp = fopen(path, "r");
    if (!fp) return -2;
    char line[128]; uint32_t row = 0;
    while (fgets(line, sizeof(line), fp) && row < m->depth) {
        if (line[0] == '\n' || line[0] == '/' || line[0] == '@') continue; /* skip comments / addr markers */
        uint32_t v = (uint32_t)strtoul(line, NULL, 16);
        int row_bytes = m->row_bytes < 4 ? m->row_bytes : 4;
        for (int b = 0; b < row_bytes; b++)
            m->data[row * m->row_bytes + b] = (v >> (b * 8)) & 0xFF;
        row++;
    }
    fclose(fp);
    return (int)row;
}

/* Per-clock memory eval.
 *
 * RTeAAL FIX 2026-05-16: honour Yosys $mem_v2 RD_CLK semantics.
 *   rd_clk_len == 0 → asynchronous read (transparent, this cycle)
 *   rd_clk_len  > 0 → synchronous read (data appears NEXT cycle)
 *
 * cv32e40x's instruction-fetch ROM and core SRAM are emitted as sync-read
 * by Yosys (rd_clk_len=1).  Driving RD_DATA combinationally caused the OBI
 * handshake to see data 1 cycle early → outstnd_cnt runaway → prefetcher
 * stuck at BootAddr.  Shadow registers below provide the 1-cycle delay.
 *
 * Writes are synchronous (committed at end of cycle), unchanged. */
#define TSIM_MAX_RD_PORTS 4
static uint32_t _tsim_rd_data_reg[N_MEMS > 0 ? N_MEMS : 1][TSIM_MAX_RD_PORTS];

static inline void tsim_mem_cycle(uint8_t *li) {
    for (int m = 0; m < N_MEMS; m++) {
        const OimMem *om = &oim_mems[m];
        TsimMem *ms = &_tsim_mems[m];
        int is_sync_read = (om->rd_clk_len > 0);
        /* Read ports */
        for (int p = 0; p < om->n_rd_ports; p++) {
            if (is_sync_read) {
                /* SYNC READ: drive LI from the SHADOW (latched last cycle), then
                 * sample the CURRENT addr → shadow for next cycle.  This realises
                 * the 1-cycle read latency that cv32e40x's OBI handshake assumes. */
                uint32_t prev_val = _tsim_rd_data_reg[m][p];
                for (int w = 0; w < om->width && w < 32; w++) {
                    int yi = oim_bit_indices[om->rd_data_start + p * om->width + w];
                    if (yi >= 0) li[yi] = (prev_val >> w) & 1;
                }
                /* Sample addr now (gated by rd_en if present) → shadow for next cycle */
                int rd_en_ok = 1;
                if (om->rd_en_len > 0) {
                    int ei = oim_bit_indices[om->rd_en_start + p];
                    if (ei >= 0) rd_en_ok = li[ei] ? 1 : 0;
                }
                if (rd_en_ok) {
                    uint32_t addr = 0;
                    for (int a = 0; a < om->abits; a++) {
                        int ai = oim_bit_indices[om->rd_addr_start + p * om->abits + a];
                        if (ai >= 0 && li[ai]) addr |= (1u << a);
                    }
                    uint32_t row_addr = (addr >= ms->depth) ? 0 : addr;
                    uint32_t row_val = 0;
                    int row_bytes = ms->row_bytes < 4 ? ms->row_bytes : 4;
                    for (int b = 0; b < row_bytes; b++)
                        row_val |= ((uint32_t)ms->data[row_addr * ms->row_bytes + b]) << (b * 8);
                    _tsim_rd_data_reg[m][p] = row_val;
                }
                continue;  /* skip the async-read path below */
            }
            /* ASYNC READ (rd_clk_len == 0): combinational transparent read */
            uint32_t addr = 0;
            for (int a = 0; a < om->abits; a++) {
                int ai = oim_bit_indices[om->rd_addr_start + p * om->abits + a];
                if (ai >= 0 && li[ai]) addr |= (1u << a);
            }
            uint32_t row_addr = addr;
            if (row_addr >= ms->depth) row_addr = 0;
            uint32_t row_val = 0;
            int row_bytes = ms->row_bytes < 4 ? ms->row_bytes : 4;
            for (int b = 0; b < row_bytes; b++)
                row_val |= ((uint32_t)ms->data[row_addr * ms->row_bytes + b]) << (b * 8);
            /* Drive RD_DATA bits */
            for (int w = 0; w < om->width && w < 32; w++) {
                int yi = oim_bit_indices[om->rd_data_start + p * om->width + w];
                if (yi >= 0) li[yi] = (row_val >> w) & 1;
            }
        }
        /* Write ports (synchronous, posedge committed at end of comb eval) */
        for (int p = 0; p < om->n_wr_ports; p++) {
            uint32_t wr_en_mask = 0;
            int wr_en_any = 0;
            int en_bits = om->width;
            for (int w = 0; w < en_bits && w < 32; w++) {
                int ei = oim_bit_indices[om->wr_en_start + p * en_bits + w];
                if (ei >= 0 && li[ei]) { wr_en_any = 1; wr_en_mask |= (1u << w); }
            }
            if (!wr_en_any) continue;
            uint32_t addr = 0;
            for (int a = 0; a < om->abits; a++) {
                int ai = oim_bit_indices[om->wr_addr_start + p * om->abits + a];
                if (ai >= 0 && li[ai]) addr |= (1u << a);
            }
            if (addr >= ms->depth) continue;
            uint32_t wr_val = 0;
            for (int w = 0; w < om->width && w < 32; w++) {
                int di = oim_bit_indices[om->wr_data_start + p * om->width + w];
                if (di >= 0 && li[di]) wr_val |= (1u << w);
            }
            uint32_t cur = 0;
            int row_bytes = ms->row_bytes < 4 ? ms->row_bytes : 4;
            for (int b = 0; b < row_bytes; b++)
                cur |= ((uint32_t)ms->data[addr * ms->row_bytes + b]) << (b * 8);
            uint32_t new_val = (cur & ~wr_en_mask) | (wr_val & wr_en_mask);
            for (int b = 0; b < row_bytes; b++)
                ms->data[addr * ms->row_bytes + b] = (new_val >> (b * 8)) & 0xFF;
        }
    }
}
#endif /* N_MEMS */

/* Gather bits from LI into uint64_t (LSB-first, respecting constants).
 * bit_indices: >= 0 = LI index, -1 = constant 0, -2 = constant 1.
 * cmask/cval provide fast-path for bits 0-63 (legacy). */
static inline uint64_t gather_bits(const uint8_t *li, const int32_t *idx,
                                    int len, uint64_t cmask, uint64_t cval) {
    uint64_t val = 0;
    for (int i = 0; i < len && i < 64; i++) {
        if (idx[i] >= 0) {
            if (li[idx[i]]) val |= (1ULL << i);
        } else if (idx[i] == -2) {
            val |= (1ULL << i);  /* constant 1 */
        }
        /* idx[i] == -1 → constant 0, nothing to do */
    }
    return val;
}

/* Scatter bits from uint64_t to LI */
static inline void scatter_bits(uint8_t *li, const int32_t *idx,
                                 int len, uint64_t val) {
    for (int i = 0; i < len && i < 64; i++) {
        if (idx[i] >= 0)
            li[idx[i]] = (val >> i) & 1;
    }
}

/* Scatter bits with change tracking for activity-aware eval */
static inline void scatter_bits_tracked(uint8_t *li, uint8_t *changed,
                                         const int32_t *idx,
                                         int len, uint64_t val) {
    for (int i = 0; i < len && i < 64; i++) {
        if (idx[i] >= 0) {
            uint8_t nv = (val >> i) & 1;
            if (li[idx[i]] != nv) {
                changed[idx[i]] = 1;
                li[idx[i]] = nv;
            }
        }
    }
}

/* Check if any input bit has changed (for activity-aware skip) */
static inline int any_input_changed(const uint8_t *changed, const int32_t *idx, int len) {
    for (int i = 0; i < len; i++) {
        if (idx[i] >= 0 && changed[idx[i]])
            return 1;
    }
    return 0;
}

static inline int64_t sext(uint64_t val, int w) {
    if (w <= 0 || w >= 64) return (int64_t)val;
    if (val & (1ULL << (w - 1)))
        return (int64_t)(val | ~((1ULL << w) - 1));
    return (int64_t)val;
}

static inline uint64_t wmask(int w) {
    return (w >= 64) ? ~0ULL : (1ULL << w) - 1;
}

/* Core op_r[n] — evaluate one Yosys cell */
static inline uint64_t op_eval(uint8_t ct, uint64_t a, uint64_t b, uint64_t s,
                                int aw, int bw, int yw,
                                int as, int bs, int sw) {
    int64_t sa, sb;
    switch (ct) {
    case CT_ADD:
        sa = as ? sext(a, aw) : (int64_t)a;
        sb = bs ? sext(b, bw) : (int64_t)b;
        return (uint64_t)(sa + sb) & wmask(yw);
    case CT_SUB:
        sa = as ? sext(a, aw) : (int64_t)a;
        sb = bs ? sext(b, bw) : (int64_t)b;
        return (uint64_t)(sa - sb) & wmask(yw);
    case CT_MUL:
        sa = as ? sext(a, aw) : (int64_t)a;
        sb = bs ? sext(b, bw) : (int64_t)b;
        if (aw + bw > 64) {
            __int128 p = (__int128)sa * (__int128)sb;
            return (uint64_t)p & wmask(yw);
        }
        return (uint64_t)(sa * sb) & wmask(yw);
    case CT_AND:  return (a & b) & wmask(yw);
    case CT_OR:   return (a | b) & wmask(yw);
    case CT_XOR:  return (a ^ b) & wmask(yw);
    case CT_NOT:  return (~a) & wmask(yw);
    case CT_MUX:  return (s & 1) ? b : a;
    case CT_PMUX:
        if (s == 0) return a;
        for (int k = 0; k < sw && k < 64; k++)
            if (s & (1ULL << k))
                return (b >> (k * yw)) & wmask(yw);
        return a;
    case CT_EQ:   return (a == b) ? 1 : 0;
    case CT_NE:   return (a != b) ? 1 : 0;
    case CT_LT:
        if (as && bs) return (sext(a, aw) < sext(b, bw)) ? 1 : 0;
        return (a < b) ? 1 : 0;
    case CT_GE:
        if (as && bs) return (sext(a, aw) >= sext(b, bw)) ? 1 : 0;
        return (a >= b) ? 1 : 0;
    case CT_SHL:
        /* tsim fix 2026-05-17: shift amount >= 64 must zero the result, not
         * wrap mod 64.  The old `b & 63` was a C-UB workaround but produced
         * wrong values when b >= 64 (e.g. shr-by-65 returned shr-by-1).
         * Found via tests/cells/shr_by_65. */
        if (b >= 64) return 0;
        return (a << b) & wmask(yw);
    case CT_SHR:
        if (b >= 64) return 0;
        if (as) return (uint64_t)(sext(a, aw) >> b) & wmask(yw);
        return (a >> b) & wmask(yw);
    case CT_LOGIC_AND: return ((a != 0) && (b != 0)) ? 1 : 0;
    case CT_LOGIC_NOT: return (a == 0) ? 1 : 0;
    case CT_REDUCE_OR: return (a != 0) ? 1 : 0;
    /* Extended types for CROC SoC */
    case CT_LOGIC_OR:  return ((a != 0) || (b != 0)) ? 1 : 0;
    case CT_GT:
        if (as && bs) return (sext(a, aw) > sext(b, bw)) ? 1 : 0;
        return (a > b) ? 1 : 0;
    case CT_SSHR:
        /* Arithmetic right shift: vacated bits = sign of A.  Shift >= aw
         * collapses to sign-extension (all 1s if negative, 0 if positive). */
        if (b >= (uint64_t)aw) {
            int64_t sa2 = sext(a, aw);
            return (sa2 < 0 ? ~0ULL : 0ULL) & wmask(yw);
        }
        return (uint64_t)(sext(a, aw) >> b) & wmask(yw);
    case CT_REDUCE_AND:return (a == wmask(aw)) ? 1 : 0;
    case CT_REDUCE_BOOL: return (a != 0) ? 1 : 0;
    /* RTeAAL note 2026-05-16 (Bug A from audit): $reduce_xor / $reduce_xnor
     * not defined as CT_* enum constants because gen_oim_fullchip.py doesn't
     * emit them on this design (parity gates absent in cv32e40x + PQC coproc).
     * If a future design uses them, add enum CT_REDUCE_XOR/_XNOR and:
     *   uint64_t x = a & wmask(aw);
     *   x ^= x>>32; x ^= x>>16; x ^= x>>8; x ^= x>>4; x ^= x>>2; x ^= x>>1;
     *   return x & 1;  // (or ~x & 1 for xnor)
     */
    case CT_BMUX: {
        /* Binary-addressed mux: A has 2^S_WIDTH * WIDTH bits, S selects which WIDTH-wide slice */
        int idx = (int)(s & wmask(sw));
        return (a >> (idx * yw)) & wmask(yw);
    }
    case CT_DEMUX: {
        /* Demux: single input A, S selects which output slice gets the value */
        /* For simulation, output is 0 everywhere except at position S */
        /* Y has 2^S * WIDTH bits; Y[S*WIDTH +: WIDTH] = A, rest = 0 */
        /* But since we store Y as one word, just output A shifted to position S */
        int idx = (int)(s & wmask(sw));
        return (a & wmask(yw / (1 << sw))) << (idx * (yw / (1 << sw)));
    }
    case CT_BWMUX:
        /* Bitwise mux: y[i] = s[i] ? b[i] : a[i].
         * If sw==1, semantics are scalar select (s==1 → b, s==0 → a).
         * Otherwise treat S as a per-bit select mask (Yosys semantics for
         * $bwmux: Y = (S & B) | (~S & A)).  Without this case PRNG-state
         * D-input ends up always 0 in OIM eval -> masking gadget appears
         * functionally unmasked. */
        if (sw == 1) return (s & 1) ? b : a;
        return ((b & s) | (a & ~s)) & wmask(yw);
    case CT_SHIFTX:
    case CT_SHIFT:
        /* tsim fix 2026-05-17: per Yosys simlib.v, $shift/$shiftx semantics are:
         *   B_SIGNED && $signed(B) < 0: Y = A << -B   (LEFT shift by abs(B))
         *   else:                       Y = A >> B    (RIGHT shift by B)
         * Previous code had the directions INVERTED, which caused every
         * variable-index slice write `y[s] = a` (synthesised by Yosys as
         * $neg + $shift) to evaluate as right-shift instead of left-shift —
         * silently producing 0 for any non-zero index, e.g. obi_demux's
         * mgr_ports_req_o[sel].req writes always landed in slot 0.
         * Found via tests/cells/demux_1bit_8slots — see rtl_change.md.
         * Also: shifts >= 64 must zero the result (not wrap), matching the
         * fix above for CT_SHL/CT_SHR. */
        if (bs) {
            int64_t sb = sext(b, bw);
            if (sb < 0) {
                uint64_t amt = (uint64_t)(-sb);
                if (amt >= 64) return 0;
                return (a << amt) & wmask(yw);
            }
            if ((uint64_t)sb >= 64) return 0;
            return (a >> sb) & wmask(yw);
        }
        if (b >= 64) return 0;
        return (a >> b) & wmask(yw);
    case CT_LE:
        if (as && bs) return (sext(a, aw) <= sext(b, bw)) ? 1 : 0;
        return (a <= b) ? 1 : 0;
    case CT_NEG:
        sa = as ? sext(a, aw) : (int64_t)a;
        return (uint64_t)(-sa) & wmask(yw);
    case CT_TRIBUF:
        /* Tristate: Y = S ? A : 'z' — treat 'z' as 0 in our toggle model. */
        return (s & 1) ? (a & wmask(yw)) : 0ULL;
    default: return 0;
    }
}

/* Pre-computed at startup: per-DFF-bit alias flag.
 * `_tsim_li_has_comb`   = LI-bit set if any comb op writes it
 * `_tsim_li_is_dff_q`   = LI-bit set if any DFF Q output uses it
 * `_tsim_dff_alias_bitmask[d]` has bit p set iff bit p of DFF d is Y/Q-aliased.
 * `_tsim_dff_force_en[d]` = 1 if DFF d's en_bit is dead (no driver and not a
 * DFF Q output) → Yosys-flatten $sdff→$sdffe artifact: treat as always-enabled. */
static uint8_t _tsim_li_has_comb[LI_SIZE];
static uint8_t _tsim_li_is_dff_q[LI_SIZE];
static uint8_t _tsim_dff_force_en[N_DFFS];
static uint64_t _tsim_dff_alias_bitmask[N_DFFS];
static int _tsim_dff_alias_inited = 0;

static inline void tsim_compute_dff_alias_flags(void) {
    if (_tsim_dff_alias_inited) return;
    memset(_tsim_li_has_comb, 0, sizeof(_tsim_li_has_comb));
    for (int lay = 0; lay < N_LAYERS; lay++) {
        int nops = oim_layer_sizes[lay];
        const OimOp *ops = oim_layers[lay];
        for (int s = 0; s < nops; s++) {
            const OimOp *op = &ops[s];
            for (int i = 0; i < op->y_len && i < 4096; i++) {
                int li = oim_bit_indices[op->y_start + i];
                if (li >= 0 && li < LI_SIZE) _tsim_li_has_comb[li] = 1;
            }
        }
    }
    /* Build _tsim_li_is_dff_q (which LI bits are DFF Q outputs) */
    memset(_tsim_li_is_dff_q, 0, sizeof(_tsim_li_is_dff_q));
    for (int d = 0; d < N_DFFS; d++) {
        const OimDff *dff = &oim_dffs[d];
        for (int i = 0; i < dff->q_bits_len && i < 4096; i++) {
            int li = oim_bit_indices[dff->q_bits_start + i];
            if (li >= 0 && li < LI_SIZE) _tsim_li_is_dff_q[li] = 1;
        }
    }
    /* tsim fix 2026-05-17 (audit item #2): also mark primary-INPUT-port LI
     * bits as "live" — the previous heuristic only checked comb-output and
     * DFF-Q drivers, so a FF whose enable signal was a primary input port
     * (e.g. `latch_strobe_i` driving `r_latched_q`'s EN) was wrongly
     * declared "dead enable" and force-enabled every cycle.  Found by
     * Verilator FF bit-match on the M1 masking cone: tsim's r_latched_q
     * updated every cycle, Verilator's only on latch_strobe_i=1. */
    {
        for (int pi = 0; pi < N_PORTS; pi++) {
            if (oim_ports[pi].dir != 0) continue;  /* 0 = input, 1 = output */
            const int16_t *bits = port_bit_arrays[pi];
            int w = oim_ports[pi].width;
            for (int b = 0; b < w; b++) {
                if (bits[b] >= 0 && bits[b] < LI_SIZE)
                    _tsim_li_has_comb[bits[b]] = 1;  /* treat as "has driver" */
            }
        }
    }
    /* Per-DFF: alias bitmask + force_en flag */
    int n_force_en = 0, n_alias = 0;
    for (int d = 0; d < N_DFFS; d++) {
        _tsim_dff_alias_bitmask[d] = 0;
        _tsim_dff_force_en[d] = 0;
        const OimDff *dff = &oim_dffs[d];
        int n = (dff->q_bits_len < dff->d_bits_len) ? dff->q_bits_len : dff->d_bits_len;
        if (n > 64) n = 64;
        for (int i = 0; i < n; i++) {
            int q_li = oim_bit_indices[dff->q_bits_start + i];
            int d_li = oim_bit_indices[dff->d_bits_start + i];
            if (q_li < 0 || d_li < 0 || q_li == d_li) continue;
            if (_tsim_li_has_comb[q_li] && !_tsim_li_has_comb[d_li]) {
                _tsim_dff_alias_bitmask[d] |= (1ULL << i);
            }
        }
        if (_tsim_dff_alias_bitmask[d]) n_alias++;
        /* Detect dead enable: en_bit references LI that no comb op or DFF Q drives */
        if (dff->en_bit >= 0 && dff->en_bit < LI_SIZE) {
            if (!_tsim_li_has_comb[dff->en_bit] && !_tsim_li_is_dff_q[dff->en_bit]) {
                _tsim_dff_force_en[d] = 1;
                n_force_en++;
            }
        }
    }
    fprintf(stderr, "  [tsim] OIM analysis: %d/%d DFFs with Y/Q alias bits, %d/%d DFFs with dead-enable (forced enabled)\n",
            n_alias, N_DFFS, n_force_en, N_DFFS);
#ifdef DEBUG_FF
    /* Diagnostic: dump status of a specific FF */
    {
        int d = DEBUG_FF;
        if (d >= 0 && d < N_DFFS) {
            const OimDff *dff = &oim_dffs[d];
            fprintf(stderr, "  [tsim] DEBUG ff_%d: q_bits_start=%u, q_bits_len=%u, d_bits_start=%u, d_bits_len=%u, en_bit=%d, en_pol=%d, has_arst=%d\n",
                    d, dff->q_bits_start, dff->q_bits_len, dff->d_bits_start, dff->d_bits_len,
                    dff->en_bit, dff->en_pol, dff->has_arst);
            for (int i = 0; i < dff->q_bits_len && i < 4; i++) {
                int q_li = oim_bit_indices[dff->q_bits_start + i];
                int d_li = oim_bit_indices[dff->d_bits_start + i];
                fprintf(stderr, "  [tsim]   bit %d: q_li=%d (comb=%d, dff_q=%d) d_li=%d (comb=%d, dff_q=%d)\n",
                        i, q_li,
                        (q_li >= 0 && q_li < LI_SIZE) ? _tsim_li_has_comb[q_li] : -1,
                        (q_li >= 0 && q_li < LI_SIZE) ? _tsim_li_is_dff_q[q_li] : -1,
                        d_li,
                        (d_li >= 0 && d_li < LI_SIZE) ? _tsim_li_has_comb[d_li] : -1,
                        (d_li >= 0 && d_li < LI_SIZE) ? _tsim_li_is_dff_q[d_li] : -1);
            }
            fprintf(stderr, "  [tsim]   alias_bitmask=0x%lx force_en=%d\n",
                    (unsigned long)_tsim_dff_alias_bitmask[d], _tsim_dff_force_en[d]);
            if (dff->en_bit >= 0 && dff->en_bit < LI_SIZE) {
                fprintf(stderr, "  [tsim]   en_bit LI[%d]: comb_driver=%d, dff_q=%d\n",
                        dff->en_bit, _tsim_li_has_comb[dff->en_bit], _tsim_li_is_dff_q[dff->en_bit]);
            }
        }
    }
#endif
    _tsim_dff_alias_inited = 1;
}

static inline void tsim_init(TsimState *st) {
    memset(st, 0, sizeof(TsimState));
    for (int d = 0; d < N_DFFS; d++) {
        const OimDff *dff = &oim_dffs[d];
        scatter_bits(st->li, &oim_bit_indices[dff->q_bits_start],
                     dff->q_bits_len, dff->reset_val);
    }
    memcpy(st->prev_li, st->li, sizeof(st->li));
    tsim_compute_dff_alias_flags();
#ifdef N_MEMS
    tsim_mems_init_storage();
#endif
}

/* popcount for uint8_t XOR result */
static inline int popcount8(uint8_t x) {
    x = x - ((x >> 1) & 0x55);
    x = (x & 0x33) + ((x >> 2) & 0x33);
    return (x + (x >> 4)) & 0x0F;
}

/* Count toggles after a cycle: compare li[] vs prev_li[] */
static inline void tsim_count_toggles(TsimState *st) {
    st->toggle.total = 0;
    memset(st->toggle.per_type, 0, sizeof(st->toggle.per_type));
    memset(st->toggle.per_dff,  0, sizeof(st->toggle.per_dff));
    memset(st->toggle.per_dff_qmask, 0, sizeof(st->toggle.per_dff_qmask));

    /* Count total bit flips */
    for (int i = 0; i < LI_SIZE; i++) {
        if (st->li[i] != st->prev_li[i])
            st->toggle.total++;
    }

    /* Count per-cell-type toggles: for each cell, check its Y bits */
    for (int lay = 0; lay < N_LAYERS; lay++) {
        int nops = oim_layer_sizes[lay];
        const OimOp *ops = oim_layers[lay];
        for (int s = 0; s < nops; s++) {
            const OimOp *op = &ops[s];
            int flips = 0;
            const int32_t *y_idx = &oim_bit_indices[op->y_start];
            for (int i = 0; i < op->y_len && i < 64; i++) {
                if (y_idx[i] >= 0) {
                    if (st->li[y_idx[i]] != st->prev_li[y_idx[i]])
                        flips++;
                }
            }
            st->toggle.per_type[op->cell_type] += flips;
        }
    }

    /* Per-FF (cca.pc) Q-bit flip counting + per-bit qmask + add to total */
    for (int d = 0; d < N_DFFS; d++) {
        const OimDff *dff = &oim_dffs[d];
        const int32_t *q_idx = &oim_bit_indices[dff->q_bits_start];
        uint64_t ff_flips = 0;
        uint32_t qmask = 0;
        int hi = (dff->q_bits_len < 32) ? dff->q_bits_len : 32;
        for (int i = 0; i < hi; i++) {
            if (q_idx[i] >= 0 && st->li[q_idx[i]] != st->prev_li[q_idx[i]]) {
                ff_flips++;
                qmask |= (1u << i);
            }
        }
        st->toggle.per_dff[d] = ff_flips;
        st->toggle.per_dff_qmask[d] = qmask;
        st->toggle.total      += ff_flips;
        st->toggle.cum_per_dff[d] += ff_flips;
    }

    st->toggle.cum_total += st->toggle.total;
    for (int ct = 0; ct < N_CELL_TYPES; ct++)
        st->toggle.cum_per_type[ct] += st->toggle.per_type[ct];
}

/* Per-FF cca.pc CSV emit (one row per cycle, columns: cycle,ff_0,ff_1,...,ff_{N_DFFS-1}) */
static inline void tsim_per_ff_csv_header(FILE *fp) {
    fprintf(fp, "cycle");
    for (int d = 0; d < N_DFFS; d++)
        fprintf(fp, ",ff_%d", d);
    fprintf(fp, "\n");
}
static inline void tsim_per_ff_csv_row(FILE *fp, TsimState *st) {
    fprintf(fp, "%lu", (unsigned long)st->cycle);
    for (int d = 0; d < N_DFFS; d++)
        fprintf(fp, ",%lu", (unsigned long)st->toggle.per_dff[d]);
    fprintf(fp, "\n");
}

/* Per-bit qmask CSV emit (one row per cycle: cycle, then N_DFFS hex uint32 columns).
   Each column is the bit-vector of which Q bits flipped this cycle for that cell.
   Python expands to per-bit boolean tensor for per-bit Welch t. */
static inline void tsim_per_ff_qmask_csv_header(FILE *fp) {
    fprintf(fp, "cycle");
    for (int d = 0; d < N_DFFS; d++)
        fprintf(fp, ",ff_%d_q", d);
    fprintf(fp, "\n");
}
static inline void tsim_per_ff_qmask_csv_row(FILE *fp, TsimState *st) {
    fprintf(fp, "%lu", (unsigned long)st->cycle);
    for (int d = 0; d < N_DFFS; d++)
        fprintf(fp, ",%08x", st->toggle.per_dff_qmask[d]);
    fprintf(fp, "\n");
}

/* Write toggle CSV header */
static inline void tsim_toggle_csv_header(FILE *fp) {
    fprintf(fp, "cycle,total_toggles");
    static const char *ct_names[] = {
        "add","sub","mul","and","or","xor","not","mux","pmux",
        "eq","ne","lt","ge","shl","shr","logic_and","logic_not","reduce_or",
        "logic_or","gt","sshr","reduce_and","reduce_bool","bmux","demux"
    };
    int n = sizeof(ct_names)/sizeof(ct_names[0]);
    for (int i = 0; i < N_CELL_TYPES && i < n; i++)
        fprintf(fp, ",%s", ct_names[i]);
    fprintf(fp, "\n");
}

/* Write one toggle CSV row */
static inline void tsim_toggle_csv_row(FILE *fp, TsimState *st) {
    fprintf(fp, "%lu,%lu", (unsigned long)st->cycle, (unsigned long)st->toggle.total);
    for (int i = 0; i < N_CELL_TYPES; i++)
        fprintf(fp, ",%lu", (unsigned long)st->toggle.per_type[i]);
    fprintf(fp, "\n");
}

/* Special $pmux evaluation: read wide B, select slice based on S */
static inline uint64_t eval_pmux(const uint8_t *li, const OimOp *op) {
    /* A = default (y_len bits) */
    uint64_t a = gather_bits(li, &oim_bit_indices[op->a_start],
                              op->a_len, op->a_const_mask, op->a_const_val);
    /* S = select bits */
    uint64_t s = gather_bits(li, &oim_bit_indices[op->s_start],
                              op->s_len, op->s_const_mask, op->s_const_val);
    if (s == 0) return a;

    /* Find which S bit is set, then gather the corresponding WIDTH-wide slice from B */
    int width = op->y_len;  /* WIDTH = Y width */
    const int32_t *b_idx = &oim_bit_indices[op->b_start];
    for (int k = 0; k < op->s_width && k < 64; k++) {
        if (s & (1ULL << k)) {
            /* Gather bits [k*width .. (k+1)*width-1] from B */
            int offset = k * width;
            uint64_t slice = 0;
            for (int i = 0; i < width && i < 64; i++) {
                int bi = offset + i;
                if (bi < op->b_len) {
                    int32_t bidx = b_idx[bi];
                    if (bidx >= 0) {
                        if (li[bidx]) slice |= (1ULL << i);
                    } else if (bidx == -2) {
                        slice |= (1ULL << i);  /* constant 1 */
                    }
                }
            }
            return slice;
        }
    }
    return a;
}

/* Special $bmux evaluation: binary-addressed mux with wide A
 * A has 2^S_WIDTH * WIDTH bits, S (binary) selects which WIDTH-wide slice.
 * Used for register file reads: A=1024b (32 regs x 32b), S=5b, Y=32b.
 * Constant bits beyond index 63 cannot use a_const_mask (only 64 bits).
 * For those, the bit_indices entry is -1 and we treat it as 0. */
static inline uint64_t eval_bmux(const uint8_t *li, const OimOp *op) {
    /* S = binary select */
    uint64_t s = gather_bits(li, &oim_bit_indices[op->s_start],
                              op->s_len, op->s_const_mask, op->s_const_val);
    int idx = (int)(s & ((1ULL << op->s_len) - 1));
    int width = op->y_len;
    int offset = idx * width;

    /* Gather bits [offset .. offset+width-1] from A */
    const int32_t *a_idx = &oim_bit_indices[op->a_start];
    uint64_t slice = 0;
    for (int i = 0; i < width && i < 64; i++) {
        int ai = offset + i;
        if (ai < op->a_len) {
            int32_t bidx = a_idx[ai];
            if (bidx >= 0) {
                if (li[bidx]) slice |= (1ULL << i);
            } else if (bidx == -2) {
                slice |= (1ULL << i);  /* constant 1 */
            }
        }
    }
    return slice;
}

/* Special $demux evaluation: A placed at position S in wide output
 * Y has 2^S_WIDTH * WIDTH bits, Y[S*WIDTH +: WIDTH] = A, rest = 0 */
static inline uint64_t eval_demux(const uint8_t *li, const OimOp *op) {
    /* For bit-level sim, we only return the 64-bit portion.
     * But demux output goes to Y bits which are scattered.
     * Since Y can be wider than 64, we handle it specially. */
    uint64_t a = gather_bits(li, &oim_bit_indices[op->a_start],
                              op->a_len, op->a_const_mask, op->a_const_val);
    uint64_t s = gather_bits(li, &oim_bit_indices[op->s_start],
                              op->s_len, op->s_const_mask, op->s_const_val);
    /* Write 0 to all Y bits first, then set the selected slice */
    /* This is handled by the caller via scatter — here we return 0,
     * and do a special scatter in tsim_cycle() */
    (void)a; (void)s;
    return 0; /* placeholder — real work done in eval+scatter below */
}

/* ===== VCD Waveform Output (#7) ===== */

/* VCD identifier characters: use printable ASCII starting from '!' */
/* VCD identifier: use safe printable chars (a-z, A-Z, 0-9 range) */
static inline void vcd_id_str(int idx, char *buf) {
    /* Generate 2-char ID: first char a-z, second char 0-9 or a-z */
    buf[0] = 'a' + (idx / 26) % 26;
    buf[1] = 'a' + (idx % 26);
    buf[2] = '\0';
}

/* Write VCD header with port signal names */
static inline void tsim_vcd_init(TsimState *st, FILE *fp) {
    st->vcd_fp = fp;
    st->vcd_enable = 1;

    fprintf(fp, "$date\n  tsim VCD output\n$end\n");
    fprintf(fp, "$version\n  tsim v0.2 RTeAAL RTL Tensor Simulator\n$end\n");
    fprintf(fp, "$timescale 1ns $end\n");
    fprintf(fp, "$scope module h_v1_coprocessor_top $end\n");

    /* Emit clock signal (synthesized — toggles every half-cycle) */
    fprintf(fp, "$var wire 1 ck clk_i $end\n");

    /* Emit port signals (skip clk_i — handled separately) */
    for (int pi = 0; pi < N_PORTS; pi++) {
        int w = oim_ports[pi].width;
        if (strcmp(oim_ports[pi].name, "clk_i") == 0) continue;
        char id[4]; vcd_id_str(pi, id);
        if (w == 1)
            fprintf(fp, "$var wire 1 %s %s $end\n", id, oim_ports[pi].name);
        else
            fprintf(fp, "$var wire %d %s %s [%d:0] $end\n",
                    w, id, oim_ports[pi].name, w - 1);
    }

    fprintf(fp, "$upscope $end\n");
    fprintf(fp, "$enddefinitions $end\n");

    /* Initial values */
    fprintf(fp, "#0\n$dumpvars\n");
    fprintf(fp, "0ck\n");  /* clock starts low */
    for (int pi = 0; pi < N_PORTS; pi++) {
        if (strcmp(oim_ports[pi].name, "clk_i") == 0) continue;
        char id[4]; vcd_id_str(pi, id);
        int w = oim_ports[pi].width;
        const int16_t *bits = port_bit_arrays[pi];
        if (w == 1) {
            int v = (bits[0] >= 0) ? st->li[bits[0]] : 0;
            fprintf(fp, "%d%s\n", v, id);
        } else {
            fprintf(fp, "b");
            for (int b = w - 1; b >= 0; b--) {
                int v = (bits[b] >= 0) ? st->li[bits[b]] : 0;
                fprintf(fp, "%d", v);
            }
            fprintf(fp, " %s\n", id);
        }
    }
    fprintf(fp, "$end\n");

    /* Save initial state as prev for first cycle comparison */
    memcpy(st->prev_li, st->li, sizeof(st->li));
}

/* Dump VCD changes for one cycle with clock toggling.
   Timeline per cycle:
     #(2*cycle)     — clock rises (posedge), signals update
     #(2*cycle + 1) — clock falls (negedge)
*/
static inline void tsim_vcd_cycle(TsimState *st, FILE *fp) {
    /* Offset by 1 cycle so first transition is at #10, not #0 (avoids $dumpvars collision) */
    unsigned long t_rise = (st->cycle + 1) * 10;
    unsigned long t_fall = t_rise + 5;

    /* Rising edge — clock goes high, dump ALL port signals (full dump for correctness) */
    fprintf(fp, "#%lu\n", t_rise);
    fprintf(fp, "1ck\n");

    for (int pi = 0; pi < N_PORTS; pi++) {
        if (strcmp(oim_ports[pi].name, "clk_i") == 0) continue;
        char id[4]; vcd_id_str(pi, id);
        int w = oim_ports[pi].width;
        const int16_t *bits = port_bit_arrays[pi];

        if (w == 1) {
            int v = (bits[0] >= 0) ? st->li[bits[0]] : 0;
            fprintf(fp, "%d%s\n", v, id);
        } else {
            fprintf(fp, "b");
            for (int b = w - 1; b >= 0; b--) {
                int v = (bits[b] >= 0) ? st->li[bits[b]] : 0;
                fprintf(fp, "%d", v);
            }
            fprintf(fp, " %s\n", id);
        }
    }

    /* Falling edge — clock goes low */
    fprintf(fp, "#%lu\n", t_fall);
    fprintf(fp, "0ck\n");
}

static inline void tsim_vcd_close(TsimState *st) {
    if (st->vcd_fp) {
        fclose(st->vcd_fp);
        st->vcd_fp = NULL;
        st->vcd_enable = 0;
    }
}

/* Evaluate combinational layers ONLY — no DFF capture, no toggle/VCD.
 * Used to resolve combinational feedback loops (e.g., gnt = req). */
static inline void tsim_eval_comb(TsimState *st) {
    for (int i = 0; i < N_LAYERS; i++) {
        int nops = oim_layer_sizes[i];
        const OimOp *ops = oim_layers[i];
        for (int s = 0; s < nops; s++) {
            const OimOp *op = &ops[s];
            uint64_t r;
            if (op->is_identity) {
                const int32_t *src_idx = (op->id_src == 0)
                    ? &oim_bit_indices[op->a_start]
                    : &oim_bit_indices[op->b_start];
                uint64_t src_cmask = (op->id_src == 0) ? op->a_const_mask : op->b_const_mask;
                uint64_t src_cval  = (op->id_src == 0) ? op->a_const_val  : op->b_const_val;
                int src_len = (op->id_src == 0) ? op->a_len : op->b_len;
                r = gather_bits(st->li, src_idx, src_len, src_cmask, src_cval);
            } else if (op->cell_type == CT_PMUX) {
                r = eval_pmux(st->li, op);
            } else if (op->cell_type == CT_BMUX && op->a_len > 64) {
                r = eval_bmux(st->li, op);
            } else {
                uint64_t a = gather_bits(st->li, &oim_bit_indices[op->a_start],
                                          op->a_len, op->a_const_mask, op->a_const_val);
                uint64_t b = gather_bits(st->li, &oim_bit_indices[op->b_start],
                                          op->b_len, op->b_const_mask, op->b_const_val);
                uint64_t sv = gather_bits(st->li, &oim_bit_indices[op->s_start],
                                           op->s_len, op->s_const_mask, op->s_const_val);
                r = op_eval(op->cell_type, a, b, sv,
                            op->a_len, op->b_len, op->y_len,
                            op->a_signed, op->b_signed, op->s_width);
            }
            scatter_bits(st->li, &oim_bit_indices[op->y_start], op->y_len, r);
        }
    }
}

static inline void tsim_cycle(TsimState *st) {
    /* Activity-aware: seed changed[] by comparing current li (after port writes)
     * vs prev_li (saved BEFORE DFF update of last cycle).
     * This captures both port-write changes AND DFF Q changes. */
    if (st->activity_enable) {
        for (int i = 0; i < LI_SIZE; i++)
            st->changed[i] = (st->li[i] != st->prev_li[i]) ? 1 : 0;
    }

    /* Save snapshot for toggle counting / VCD (before we modify li in eval) */
    uint8_t pre_eval_snap[LI_SIZE];
    if (st->toggle_enable || st->vcd_enable)
        memcpy(pre_eval_snap, st->li, sizeof(pre_eval_snap));

    /* Evaluate combinational layers (RTeAAL Algorithm 3)
     * NU: ops sorted by cell_type within each layer for branch prediction.
     * Identity elision: const-S $mux ops skip eval, just copy src->Y.
     * Activity-aware: skip ops whose inputs haven't changed.
     *
     * RTeAAL FIX 2026-05-16: iterate the cascade to a fixed point each cycle
     * (Verilator-style settling).  We snapshot LI before each pass and stop
     * when no LI bit changes.  This is the principled form of the earlier
     * hardcoded 2-pass workaround — it handles arbitrary-length cross-layer
     * dependency chains, not just 1-edge ones.  Convergence is asserted: if
     * MAX_COMB_PASSES is exceeded the cycle is considered non-converged and
     * a warning is emitted (does not abort, so TVLA campaigns continue and
     * the user can audit the affected cycles). */
    #define MAX_COMB_PASSES 32
    uint8_t li_before_pass[LI_SIZE];
    int converged = 0;
    int passes_used = 0;
    for (int pass = 0; pass < MAX_COMB_PASSES; pass++) {
        memcpy(li_before_pass, st->li, sizeof(li_before_pass));
        for (int i = 0; i < N_LAYERS; i++) {
        int nops = oim_layer_sizes[i];
        const OimOp *ops = oim_layers[i];
        for (int s = 0; s < nops; s++) {
            const OimOp *op = &ops[s];

            /* Activity-aware skip: check if any input bit changed.
             * First cycle after port changes: evaluate everything.
             * Layer 0 always evaluated (primary inputs).
             * Identity/PMUX/MUX always evaluated (mux is in result path). */
            if (st->activity_enable && st->cycle >= 3 && i > 0
                && !op->is_identity
                && op->cell_type != CT_PMUX
                && op->cell_type != CT_MUX) {
                int inp_changed = any_input_changed(st->changed,
                    &oim_bit_indices[op->a_start], op->a_len);
                if (!inp_changed && op->b_len > 0)
                    inp_changed = any_input_changed(st->changed,
                        &oim_bit_indices[op->b_start], op->b_len);
                if (!inp_changed && op->s_len > 0)
                    inp_changed = any_input_changed(st->changed,
                        &oim_bit_indices[op->s_start], op->s_len);
                if (!inp_changed) {
                    st->activity_skipped++;
                    continue;
                }
            }
            if (st->activity_enable)
                st->activity_evaluated++;

            uint64_t r;

            if (op->is_identity) {
                const int32_t *src_idx = (op->id_src == 0)
                    ? &oim_bit_indices[op->a_start]
                    : &oim_bit_indices[op->b_start];
                uint64_t src_cmask = (op->id_src == 0) ? op->a_const_mask : op->b_const_mask;
                uint64_t src_cval  = (op->id_src == 0) ? op->a_const_val  : op->b_const_val;
                int src_len = (op->id_src == 0) ? op->a_len : op->b_len;
                r = gather_bits(st->li, src_idx, src_len, src_cmask, src_cval);
            } else if (op->cell_type == CT_PMUX) {
                r = eval_pmux(st->li, op);
            } else if (op->cell_type == CT_BMUX && op->a_len > 64) {
                r = eval_bmux(st->li, op);
            } else if (op->cell_type == CT_BWMUX && (op->a_len > 64 || op->s_len > 64)) {
                /* Wide bitwise mux: y[i] = s[i] ? b[i] : a[i].  S, A, B all share
                 * y_len bits.  Without this path the OIM evaluator gathers only
                 * the low 64 bits, defaults the rest of S to 0, and ends up
                 * routing A everywhere — which silently kills the PRNG seed_data
                 * path through the S-gated mux at the prng-state input. */
                const int32_t *a_idx = &oim_bit_indices[op->a_start];
                const int32_t *b_idx = &oim_bit_indices[op->b_start];
                const int32_t *s_idx = &oim_bit_indices[op->s_start];
                const int32_t *y_idx = &oim_bit_indices[op->y_start];
                for (int i = 0; i < op->y_len; i++) {
                    if (y_idx[i] < 0) continue;
                    int s_pos = (op->s_len == 1) ? 0 : i;  /* scalar S broadcasts */
                    int sv = (s_pos < op->s_len && s_idx[s_pos] >= 0)
                             ? st->li[s_idx[s_pos]] : 0;
                    int av = (i < op->a_len && a_idx[i] >= 0) ? st->li[a_idx[i]] : 0;
                    int bv = (i < op->b_len && b_idx[i] >= 0) ? st->li[b_idx[i]] : 0;
                    uint8_t nv = (uint8_t)(sv ? bv : av);
                    if (st->activity_enable && st->li[y_idx[i]] != nv)
                        st->changed[y_idx[i]] = 1;
                    st->li[y_idx[i]] = nv;
                }
                goto next_op;
            } else if (op->cell_type == CT_MUX && op->y_len > 64) {
                /* Wide scalar $mux: Y = S ? B : A, where S is 1 bit and A/B/Y
                 * are y_len bits wide.  The standard op_eval / gather_bits /
                 * scatter_bits path caps at 64 bits, so anything in bits 64..
                 * is silently zero — which kills e.g. the 128-bit PRNG-state
                 * mux that drives u_prng.s D-input bits 64..127. */
                const int32_t *a_idx = &oim_bit_indices[op->a_start];
                const int32_t *b_idx = &oim_bit_indices[op->b_start];
                const int32_t *s_idx = &oim_bit_indices[op->s_start];
                const int32_t *y_idx = &oim_bit_indices[op->y_start];
                int sv = (op->s_len > 0 && s_idx[0] >= 0) ? st->li[s_idx[0]]
                       : (op->s_len > 0 && s_idx[0] == -2) ? 1 : 0;
                for (int i = 0; i < op->y_len; i++) {
                    if (y_idx[i] < 0) continue;
                    int av, bv;
                    if (i < op->a_len) {
                        av = (a_idx[i] >= 0) ? st->li[a_idx[i]]
                           : (a_idx[i] == -2) ? 1 : 0;
                    } else av = 0;
                    if (i < op->b_len) {
                        bv = (b_idx[i] >= 0) ? st->li[b_idx[i]]
                           : (b_idx[i] == -2) ? 1 : 0;
                    } else bv = 0;
                    uint8_t nv = (uint8_t)(sv ? bv : av);
                    if (st->activity_enable && st->li[y_idx[i]] != nv)
                        st->changed[y_idx[i]] = 1;
                    st->li[y_idx[i]] = nv;
                }
                goto next_op;
            } else if (op->cell_type == CT_PMUX && op->y_len > 64) {
                /* Wide $pmux: Y has y_len bits, B has y_len*s_len bits as
                 * concatenated slices, S is a one-hot select (or 0 → A).
                 * eval_pmux() caps each slice at 64 bits, so any per-slice
                 * width > 64 silently drops the high bits.  Loop bit-by-bit. */
                const int32_t *a_idx = &oim_bit_indices[op->a_start];
                const int32_t *b_idx = &oim_bit_indices[op->b_start];
                const int32_t *s_idx = &oim_bit_indices[op->s_start];
                const int32_t *y_idx = &oim_bit_indices[op->y_start];
                /* Pick which source: A (if S==0) or one slice of B */
                int sel = -1;  /* -1 = use A */
                for (int k = 0; k < op->s_len && k < op->s_width; k++) {
                    int sv = (s_idx[k] >= 0) ? st->li[s_idx[k]]
                           : (s_idx[k] == -2) ? 1 : 0;
                    if (sv) { sel = k; break; }
                }
                int width = op->y_len;
                int offset = (sel >= 0) ? sel * width : 0;
                /* RTeAAL FIX 2026-05-16 (Bug B): src must be int32_t to match
                 * a_idx/b_idx stride; the prior int16_t* declaration read
                 * garbage indices for any wide $pmux (y_len>64) and silently
                 * corrupted PRNG-state pmux outputs. */
                const int32_t *src = (sel >= 0) ? b_idx + offset : a_idx;
                int src_len = (sel >= 0)
                              ? ((op->b_len > offset) ? op->b_len - offset : 0)
                              : op->a_len;
                for (int i = 0; i < op->y_len; i++) {
                    if (y_idx[i] < 0) continue;
                    int v = 0;
                    if (i < src_len) {
                        v = (src[i] >= 0) ? st->li[src[i]]
                          : (src[i] == -2) ? 1 : 0;
                    }
                    uint8_t nv = (uint8_t)v;
                    if (st->activity_enable && st->li[y_idx[i]] != nv)
                        st->changed[y_idx[i]] = 1;
                    st->li[y_idx[i]] = nv;
                }
                goto next_op;
            } else if (op->cell_type == CT_DEMUX && op->y_len > 64) {
                /* Wide demux: clear all Y bits, then set slice at position S */
                uint64_t a = gather_bits(st->li, &oim_bit_indices[op->a_start],
                                          op->a_len, op->a_const_mask, op->a_const_val);
                uint64_t sv = gather_bits(st->li, &oim_bit_indices[op->s_start],
                                           op->s_len, op->s_const_mask, op->s_const_val);
                int sel = (int)(sv & ((1ULL << op->s_len) - 1));
                int width = op->a_len;  /* each slice is A_WIDTH wide */
                const int32_t *y_idx = &oim_bit_indices[op->y_start];
                /* Clear all Y bits */
                for (int i = 0; i < op->y_len; i++)
                    if (y_idx[i] >= 0) st->li[y_idx[i]] = 0;
                /* Set the selected slice */
                int offset = sel * width;
                for (int i = 0; i < width && i < 64; i++) {
                    int yi = offset + i;
                    if (yi < op->y_len && y_idx[yi] >= 0)
                        st->li[y_idx[yi]] = (a >> i) & 1;
                }
                goto next_op;  /* skip normal scatter */
            } else if (op->a_len > 64 &&
                       (op->cell_type == CT_REDUCE_OR ||
                        op->cell_type == CT_REDUCE_BOOL ||
                        op->cell_type == CT_REDUCE_AND)) {
                /* Wide reduce: scan all A bits directly from LI */
                const int32_t *a_idx = &oim_bit_indices[op->a_start];
                if (op->cell_type == CT_REDUCE_AND) {
                    r = 1;
                    for (int i = 0; i < op->a_len; i++) {
                        int v = (a_idx[i] >= 0) ? st->li[a_idx[i]] :
                                (a_idx[i] == -2) ? 1 : 0;
                        if (!v) { r = 0; break; }
                    }
                } else { /* REDUCE_OR / REDUCE_BOOL */
                    r = 0;
                    for (int i = 0; i < op->a_len; i++) {
                        int v = (a_idx[i] >= 0) ? st->li[a_idx[i]] :
                                (a_idx[i] == -2) ? 1 : 0;
                        if (v) { r = 1; break; }
                    }
                }
            } else {
                uint64_t a = gather_bits(st->li, &oim_bit_indices[op->a_start],
                                          op->a_len, op->a_const_mask, op->a_const_val);
                uint64_t b = gather_bits(st->li, &oim_bit_indices[op->b_start],
                                          op->b_len, op->b_const_mask, op->b_const_val);
                uint64_t sv = gather_bits(st->li, &oim_bit_indices[op->s_start],
                                           op->s_len, op->s_const_mask, op->s_const_val);
                r = op_eval(op->cell_type, a, b, sv,
                            op->a_len, op->b_len, op->y_len,
                            op->a_signed, op->b_signed, op->s_width);
            }

            if (st->activity_enable)
                scatter_bits_tracked(st->li, st->changed,
                    &oim_bit_indices[op->y_start], op->y_len, r);
            else
                scatter_bits(st->li, &oim_bit_indices[op->y_start], op->y_len, r);
            next_op:;
        }
        }  /* end layer loop */
        passes_used = pass + 1;
        /* Convergence check: did this pass change any LI bit? */
        if (memcmp(li_before_pass, st->li, sizeof(li_before_pass)) == 0) {
            converged = 1;
            break;
        }
    }  /* end pass loop */
    if (!converged) {
        fprintf(stderr, "[tsim WARN] cycle %llu: combinational logic did not converge in %d passes\n",
                (unsigned long long)st->cycle, MAX_COMB_PASSES);
    }
    st->comb_passes_last_cycle = passes_used;
    if (passes_used > st->comb_passes_max) st->comb_passes_max = passes_used;

    /* OIM Y/Q alias mitigation: PER-BIT decision. For each aliased Q-bit,
     * save the comb-eval-time Q value (which IS the next-state computed by
     * comb ops aliased onto Q), then restore the Q LI to pre-eval. */
    uint64_t alias_q_combval[N_DFFS];
    for (int d = 0; d < N_DFFS; d++) {
        alias_q_combval[d] = 0;
        uint64_t mask = _tsim_dff_alias_bitmask[d];
        if (!mask) continue;
        const OimDff *dff = &oim_dffs[d];
        int n = dff->q_bits_len < 64 ? dff->q_bits_len : 64;
        for (int i = 0; i < n; i++) {
            if (!(mask & (1ULL << i))) continue;
            int q_li = oim_bit_indices[dff->q_bits_start + i];
            if (q_li >= 0 && st->li[q_li])
                alias_q_combval[d] |= (1ULL << i);
            if (q_li >= 0)
                st->li[q_li] = pre_eval_snap[q_li];
        }
    }

    /* Memory cycle: drive RD_DATA combinationally, commit WR_DATA. */
#ifdef N_MEMS
    tsim_mem_cycle(st->li);
#endif

    /* Clock edge: capture all D, then update all Q */
#ifndef RST_NI_BIT
#define RST_NI_BIT 3
#endif
    int rst_active = (st->li[RST_NI_BIT] == 0);

    if (rst_active) {
        /* RTeAAL FIX 2026-05-16: only $adff/$adffe (has_arst=1) reset to
         * reset_val on async rst_ni.  Pure $dff/$dffe cells (no reset port in
         * RTL) hold across reset and are processed by the else-branch below.
         * Pre-fix clobbered ALL DFFs to 0, latent because every campaign holds
         * rst_n low for the first 4 cycles only and never asserts it mid-run. */
        for (int d = 0; d < N_DFFS; d++) {
            const OimDff *dff = &oim_dffs[d];
            if (!dff->has_arst) continue;
            if (st->activity_enable)
                scatter_bits_tracked(st->li, st->changed,
                    &oim_bit_indices[dff->q_bits_start],
                    dff->q_bits_len, dff->reset_val);
            else
                scatter_bits(st->li, &oim_bit_indices[dff->q_bits_start],
                             dff->q_bits_len, dff->reset_val);
        }
        /* Non-arst DFFs still get a normal D-capture; we do them in the else
         * branch but only for !has_arst this cycle */
    }
    {  /* D-capture for all DFFs EXCEPT has_arst-during-rst (which got reset above) */
        uint64_t nq[N_DFFS];
        for (int d = 0; d < N_DFFS; d++)
            {
                /* Skip has_arst DFFs during reset — they were set to reset_val above */
                if (rst_active && oim_dffs[d].has_arst) {
                    nq[d] = oim_dffs[d].reset_val;
                    continue;
                }
                /* Check enable for $adffe: only capture D if enabled.
                 * If the enable bit is "dead" (no comb driver, no DFF Q),
                 * treat as always-enabled — Yosys-flatten $sdff→$sdffe artifact. */
                int en_ok = 1;
                if (oim_dffs[d].en_bit >= 0 && !_tsim_dff_force_en[d]) {
                    int ev = st->li[oim_dffs[d].en_bit];
                    en_ok = (oim_dffs[d].en_pol) ? ev : !ev;
                }
                if (en_ok) {
                    uint64_t mask = _tsim_dff_alias_bitmask[d];
                    if (mask) {
                        /* Mixed: aliased bits come from saved alias_q_combval,
                         * non-aliased bits come from gather of D LI bits. */
                        uint64_t d_val = gather_bits(st->li,
                            &oim_bit_indices[oim_dffs[d].d_bits_start],
                            oim_dffs[d].d_bits_len, 0, 0);
                        nq[d] = (alias_q_combval[d] & mask) | (d_val & ~mask);
                    } else {
                        nq[d] = gather_bits(st->li, &oim_bit_indices[oim_dffs[d].d_bits_start],
                                            oim_dffs[d].d_bits_len, 0, 0);
                    }
                } else {
                    nq[d] = gather_bits(st->li, &oim_bit_indices[oim_dffs[d].q_bits_start],
                                        oim_dffs[d].q_bits_len, 0, 0);
                }
            }
        for (int d = 0; d < N_DFFS; d++) {
            if (st->activity_enable)
                scatter_bits_tracked(st->li, st->changed,
                    &oim_bit_indices[oim_dffs[d].q_bits_start],
                    oim_dffs[d].q_bits_len, nq[d]);
            else
                scatter_bits(st->li, &oim_bit_indices[oim_dffs[d].q_bits_start],
                             oim_dffs[d].q_bits_len, nq[d]);
        }
    }

    /* tsim fix 2026-05-17 (audit item #4): post-clock combinational re-eval.
     *
     * Without this pass, downstream observers (ports, LI probes) read
     * combinational signals that were computed BEFORE the FFs latched —
     * so any comb signal derived from a just-updated Q (e.g. `prng_ready_o
     * = &seed_done`) appears one cycle behind silicon and behind Verilator.
     *
     * Verilator's tick semantics: posedge clk → FF capture → comb re-eval.
     * Observers between posedges then read the post-clock comb.  Our prior
     * tsim_cycle skipped this re-eval — found by the Verilator oracle on
     * the masking cone (cycle 5: tsim=ready=0 vs Verilator=ready=1).
     *
     * Single-pass `tsim_eval_comb` is sufficient for any acyclic cone (no
     * non-trivial SCCs); for a design with comb feedback the outer
     * fixed-point loop above already settled the cycle's main comb fabric,
     * so a single re-eval after FF capture re-propagates the new Q values
     * through the (assumed-acyclic) consumer fanout in one pass.  If a
     * future design has an SCC whose Q-fan-in feeds back to the SCC
     * itself, this pass may need to be wrapped in another fixed-point
     * loop — but the audit shows that's not the case in the masking cone
     * (M1 has 0 non-trivial SCCs) which is the camera-ready scope.
     *
     * Update 2026-05-17 (validator gating item #3): wrap in fixed-point
     * loop so the post-clock pass also converges for SoCs with comb
     * feedback. On the masking cone (0 SCCs) the loop body runs once
     * and exits — free.  On the full SoC the loop iterates until LI
     * stabilises, matching Verilator's settling semantics for SCCs
     * whose fan-in includes a just-updated Q. */
    {
        uint8_t li_before_postclock[LI_SIZE];
        for (int pcpass = 0; pcpass < MAX_COMB_PASSES; pcpass++) {
            memcpy(li_before_postclock, st->li, sizeof(li_before_postclock));
            tsim_eval_comb(st);
            if (memcmp(li_before_postclock, st->li, sizeof(li_before_postclock)) == 0)
                break;
        }
    }

    /* Toggle counting: compare final li vs pre-eval snapshot */
    if (st->toggle_enable) {
        memcpy(st->prev_li, pre_eval_snap, sizeof(st->prev_li));
        tsim_count_toggles(st);
    }

    /* VCD output: compare li vs pre-eval snapshot for changed signals */
    if (st->vcd_enable && st->vcd_fp) {
        memcpy(st->prev_li, pre_eval_snap, sizeof(st->prev_li));
        tsim_vcd_cycle(st, st->vcd_fp);
    }

    /* For activity-aware eval: save prev_li BEFORE DFF update.
     * This way, next cycle's change detection will see DFF Q changes
     * as "changed" inputs, enabling proper propagation through comb logic.
     * Note: DFFs already updated above, so we use pre_eval_snap which
     * has the state BEFORE this cycle's comb eval + DFF update. */
    if (st->activity_enable)
        memcpy(st->prev_li, pre_eval_snap, sizeof(st->prev_li));
    else
        memcpy(st->prev_li, st->li, sizeof(st->prev_li));

    st->cycle++;
}

/* Port helpers */
static inline int port_idx(const char *name) {
    for (int i = 0; i < N_PORTS; i++)
        if (strcmp(oim_ports[i].name, name) == 0) return i;
    return -1;
}

static inline void tsim_write_port(TsimState *st, int pi, uint64_t val) {
    if (pi < 0) return;
    const int16_t *bits = port_bit_arrays[pi];
    for (int i = 0; i < oim_ports[pi].width; i++) {
        if (bits[i] >= 0) {
            uint8_t nv = (val >> i) & 1;
            if (st->activity_enable && st->li[bits[i]] != nv)
                st->changed[bits[i]] = 1;
            st->li[bits[i]] = nv;
        }
    }
}

static inline uint64_t tsim_read_port(TsimState *st, int pi) {
    if (pi < 0) return 0;
    const int16_t *bits = port_bit_arrays[pi];
    uint64_t val = 0;
    for (int i = 0; i < oim_ports[pi].width; i++)
        if (bits[i] >= 0 && st->li[bits[i]]) val |= (1ULL << i);
    return val;
}

#endif /* TSIM_H */
