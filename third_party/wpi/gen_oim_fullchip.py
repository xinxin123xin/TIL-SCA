#!/usr/bin/env python3
"""
gen_oim_bitlvl.py — Bit-level OIM Generator for genuine RTeAAL simulator.

Unlike the word-level gen_oim.py, this operates at individual net-bit
granularity — matching the actual RTeAAL reference implementation where
each LI entry is one signal value.

For H_v1 (7536 bits, 283 cells), bit-level is:
  - Simple and correct (no mixed-source ambiguity)
  - Fast enough (microseconds per cycle for this design size)
  - Matches the RTeAAL paper's formulation exactly

Each cell becomes ONE operation that reads N input bits and writes M output bits.
The evaluator reconstructs multi-bit values from individual LI bits, computes,
then scatters output bits back to LI.

Output: oim_bitlvl.h — C header with all OIM arrays
"""

import json
import sys
from collections import defaultdict, deque
from pathlib import Path

# Cell type enum — extended for full SoC (CROC+H_v1)
CELL_TYPES = {
    "$add": 0, "$sub": 1, "$mul": 2, "$and": 3, "$or": 4, "$xor": 5,
    "$not": 6, "$mux": 7, "$pmux": 8, "$eq": 9, "$ne": 10, "$lt": 11,
    "$ge": 12, "$shl": 13, "$shr": 14, "$logic_and": 15, "$logic_not": 16,
    "$reduce_or": 17,
    # Extended types for CROC SoC
    "$logic_or": 18, "$gt": 19, "$sshr": 20,
    "$reduce_and": 21, "$reduce_bool": 22,
    "$bmux": 23, "$demux": 24,
    # Added for full_soc (cv32e40x produces these): bitwise mux + shift-with-x
    "$bwmux": 25, "$shiftx": 26,
    # Further cv32e40x cells observed in full_soc flatten
    "$shift": 27, "$le": 28, "$neg": 29, "$tribuf": 30,
}
# All DFF flavors Yosys can produce after opt_dff (without -sat -nodffe -nosdff).
# $dffe/$sdff added for cv32e40x pipeline registers.
DFF_TYPES = {"$adff", "$dff", "$adffe", "$aldff", "$dffe", "$sdff"}
# Memory cells kept intact when memory_map is skipped.
MEM_TYPES = {"$mem_v2", "$mem"}
CELL_TYPE_NAMES = {v: k for k, v in CELL_TYPES.items()}


def parse_param_int(p, default=0):
    """Yosys binary-string param → int (handles '', 'x', plain int)."""
    if p is None or p == "":
        return default
    if isinstance(p, int):
        return p
    if p == "x":
        return default
    try:
        return int(p, 2)
    except (TypeError, ValueError):
        return default


def load_json(path):
    with open(path) as f:
        return json.load(f)


def bits_to_const(bits):
    """Convert a list of bit chars ('0','1') to an integer value."""
    val = 0
    for i, b in enumerate(bits):
        if b == "1":
            val |= (1 << i)
    return val


class BitLevelOIM:
    def __init__(self, mod_data):
        self.ports = mod_data["ports"]
        self.cells = mod_data["cells"]
        self._raw_module_netnames = mod_data.get("netnames", {})

        # Separate DFFs, memories, and comb. Cells without port_directions
        # (Verilog blackboxes like EF_UART_AHBL) are skipped silently — they
        # appear in the bit_producer scan below as warnings.
        self.dffs = {}
        self.combs = {}
        self.mems = {}
        for cn, cd in self.cells.items():
            t = cd["type"]
            if not cd.get("port_directions"):
                continue  # blackbox; will be reported by bit_producer scan
            if t in DFF_TYPES:
                self.dffs[cn] = cd
            elif t in MEM_TYPES:
                self.mems[cn] = cd
            elif t in CELL_TYPES:
                self.combs[cn] = cd
            else:
                print(f"  WARN: unknown cell type {t}: {cn}")

        # Net-bit producer map: bit_id -> cell_name
        # Blackbox cells (e.g. EF_PSRAM_CTRL_AHBL) and some Yosys internal cells
        # may lack port_directions; skip them gracefully.
        self.bit_producer = {}
        self.blackbox_cells = []
        for cn, cd in self.cells.items():
            pdirs = cd.get("port_directions")
            if not pdirs:
                self.blackbox_cells.append((cn, cd.get("type", "?")))
                continue
            for pn, dr in pdirs.items():
                if dr == "output":
                    for b in cd["connections"].get(pn, []):
                        if isinstance(b, int):
                            self.bit_producer[b] = cn
        if self.blackbox_cells:
            print(f"  Blackbox/unknown cells (skipped, {len(self.blackbox_cells)}):")
            for cn, t in self.blackbox_cells[:6]:
                print(f"    {t:30s}  {cn}")
            if len(self.blackbox_cells) > 6:
                print(f"    ... ({len(self.blackbox_cells)-6} more)")

        # Port bit map
        self.port_bits = {}  # port_name -> list of bit IDs
        for pn, pd in self.ports.items():
            self.port_bits[pn] = pd["bits"]

        # Find max bit ID for LI sizing
        self.max_bit = 0
        for cd in self.cells.values():
            for bits in cd["connections"].values():
                for b in bits:
                    if isinstance(b, int):
                        self.max_bit = max(self.max_bit, b)
        for pd in self.ports.values():
            for b in pd["bits"]:
                if isinstance(b, int):
                    self.max_bit = max(self.max_bit, b)

        self.li_size = self.max_bit + 1  # bit-level: one entry per net bit
        print(f"  LI size (bit-level): {self.li_size}")
        print(f"  Comb cells: {len(self.combs)}, DFFs: {len(self.dffs)}")

        # ---------------------------------------------------------------
        # Y=Y latch-replacement detector
        # ---------------------------------------------------------------
        # When a SystemVerilog `always_comb` block omits a default for a
        # signal that is assigned only in some case arms, Yosys's `proc_dlatch`
        # pass infers a $dlatch.  Subsequent opt_dff lowers that $dlatch into
        # a $mux whose A input is electrically the same wire as Y (self-feedback,
        # Y = S ? B : Y_prev).  In tsim this evaluates to a Jacobi fixed point
        # that depends on initial LI state — typically locking the signal at 0.
        #
        # SV-2012 §9.2.2.2 forbids latches in `always_comb`.  Slang's frontend
        # attaches an `always_comb` attribute to the proc, but the attribute
        # does not survive to the final cell — we have to detect the Y=Y
        # shape structurally and key it to the source-line attribute.
        #
        # Self-feedback at the *direct* level: cd['connections']['A'] == ['Y'].
        # More commonly the feedback runs through chains of buffers and
        # constant-folded muxes, which this simple check will miss.  Treat as
        # a first-line guard: catches the obi2ahbm_adapter.sv class of bug
        # that bit us; future chained patterns surface via tests/cells.
        latch_violations = []
        latch_warnings = []
        for cn, cd in self.combs.items():
            if cd["type"] != "$mux":
                continue
            a_bits = cd["connections"].get("A", [])
            y_bits = cd["connections"].get("Y", [])
            if not a_bits or not y_bits:
                continue
            # Only flag if A and Y are integer-bit lists (named net references)
            # and at least one bit overlaps exactly between them.
            a_ints = [b for b in a_bits if isinstance(b, int)]
            y_ints = [b for b in y_bits if isinstance(b, int)]
            if not a_ints or not y_ints:
                continue
            shared = set(a_ints) & set(y_ints)
            if not shared:
                continue
            src = cd.get("attributes", {}).get("src", "?")
            # always_comb-region heuristic: slang frontend embeds the source
            # location; if it points inside an always_comb body the violation
            # is hard.  Otherwise (always @* / always_latch) it's a warning.
            entry = (cn, src, len(shared), len(y_ints))
            # We don't have a perfect "is this always_comb?" detector without
            # parsing the source file.  Conservative rule: flag every Y=Y
            # $mux as an error; explicit always_latch is rare in this codebase.
            latch_violations.append(entry)

        if latch_violations:
            print(f"  Y=Y self-feedback $mux cells (latch-replacement signature): "
                  f"{len(latch_violations)}")
            for cn, src, n_shared, n_y in latch_violations[:8]:
                print(f"    {cn[:50]:50s}  bits Y∩A = {n_shared}/{n_y}  src={src[-50:]}")
            if len(latch_violations) > 8:
                print(f"    ... ({len(latch_violations)-8} more)")
            print(f"  HINT: these are typically caused by missing defaults in "
                  f"always_comb blocks (SV-2012 §9.2.2.2 violation).  Add a "
                  f"default assignment at the top of the always_comb block, "
                  f"or convert to always_latch if intentional.")
        else:
            print(f"  Y=Y self-feedback $mux cells: 0 (clean)")

    def get_cell_dependencies(self, cn):
        """Return set of comb cell names that this cell depends on."""
        cd = self.combs[cn]
        deps = set()
        for pn, dr in cd["port_directions"].items():
            if dr == "input":
                for b in cd["connections"][pn]:
                    if isinstance(b, int) and b in self.bit_producer:
                        producer = self.bit_producer[b]
                        if producer in self.combs:
                            deps.add(producer)
        return deps

    def levelize(self):
        """Topological sort comb cells into layers.

        Uses Tarjan SCC + condensation-DAG + longest-path layering so that:
          1. Cells in a true combinational SCC share ONE layer (the simulator's
             outer fixed-point loop in tsim_cycle settles them in 2-4 passes).
          2. Cells *downstream* of an SCC get correct layer numbers — Kahn's
             alone gets stuck on them because their in-degree never decrements
             past the SCC predecessor.
          3. SCC-internal cells are emitted in reverse-postorder DFS so the
             fixed-point converges in O(depth_of_SCC) iters, not O(|SCC|).

        Diagnostic output:
          - Total SCCs, non-trivial SCC count, largest SCC size.
          - One "WARN: SCC of N cells at layer L" line per non-trivial SCC.
        """
        # 1. Predecessors per cell
        preds = {cn: self.get_cell_dependencies(cn) for cn in self.combs}
        succs = defaultdict(set)
        for cn, ps in preds.items():
            for p in ps:
                succs[p].add(cn)

        # 2. Tarjan SCC (iterative to avoid recursion blow-up on deep graphs)
        index_ctr = [0]
        stack = []
        on_stack = set()
        indices = {}
        lowlinks = {}
        sccs = []  # list of lists; each inner list is one SCC (cell names)

        def strongconnect(v0):
            work = [(v0, iter(preds[v0]))]
            indices[v0] = index_ctr[0]
            lowlinks[v0] = index_ctr[0]
            index_ctr[0] += 1
            stack.append(v0)
            on_stack.add(v0)
            while work:
                v, it = work[-1]
                try:
                    w = next(it)
                    if w not in indices:
                        indices[w] = index_ctr[0]
                        lowlinks[w] = index_ctr[0]
                        index_ctr[0] += 1
                        stack.append(w)
                        on_stack.add(w)
                        work.append((w, iter(preds[w])))
                    elif w in on_stack:
                        lowlinks[v] = min(lowlinks[v], indices[w])
                except StopIteration:
                    if lowlinks[v] == indices[v]:
                        comp = []
                        while True:
                            x = stack.pop()
                            on_stack.discard(x)
                            comp.append(x)
                            if x == v:
                                break
                        sccs.append(comp)
                    work.pop()
                    if work:
                        parent = work[-1][0]
                        lowlinks[parent] = min(lowlinks[parent], lowlinks[v])

        for v in self.combs:
            if v not in indices:
                strongconnect(v)

        # 3. SCC id per cell + condensation edges
        scc_id = {}
        for i, comp in enumerate(sccs):
            for cn in comp:
                scc_id[cn] = i
        n_sccs = len(sccs)
        scc_succs = defaultdict(set)
        scc_preds = defaultdict(set)
        for cn, ps in preds.items():
            sid = scc_id[cn]
            for p in ps:
                pid = scc_id[p]
                if pid != sid:
                    scc_preds[sid].add(pid)
                    scc_succs[pid].add(sid)

        # 4. Layer = longest path from a root in the condensation DAG.
        #    Tarjan was run on `preds` (i.e. DFS on the reversed dep graph),
        #    so sccs come out in TOPOLOGICAL order of the original graph
        #    (sources first, sinks last).  Iterate range(n_sccs) directly.
        scc_layer = [0] * n_sccs
        for sid in range(n_sccs):
            if scc_preds[sid]:
                scc_layer[sid] = max(scc_layer[p] for p in scc_preds[sid]) + 1
            else:
                scc_layer[sid] = 0

        # 5. Report non-trivial SCCs (real combinational feedback loops)
        nontriv = []
        for i, comp in enumerate(sccs):
            is_self_loop = (len(comp) == 1 and comp[0] in preds[comp[0]])
            if len(comp) > 1 or is_self_loop:
                nontriv.append((i, comp))
        if nontriv:
            print(f"  Combinational SCCs (real feedback loops): {len(nontriv)}")
            for sid, comp in sorted(nontriv, key=lambda x: -len(x[1])):
                print(f"    WARN: SCC of {len(comp):4d} cells at layer {scc_layer[sid]:3d}  "
                      f"(sample: {comp[0][:60]})")
        else:
            print(f"  Combinational SCCs (real feedback loops): 0  (pure DAG)")

        # 6. Within each non-trivial SCC, order cells by reverse-postorder DFS
        #    so the simulator's fixed-point loop converges in O(SCC_depth)
        #    iterations rather than O(|SCC|).  For trivial SCCs the single
        #    cell is its own order.
        def rpo_within_scc(comp_list):
            comp_set = set(comp_list)
            visited = set()
            order = []  # post-order
            for root in comp_list:
                if root in visited:
                    continue
                # iterative DFS
                wk = [(root, iter([s for s in succs[root] if s in comp_set]))]
                visited.add(root)
                while wk:
                    node, it = wk[-1]
                    try:
                        nxt = next(it)
                        if nxt not in visited:
                            visited.add(nxt)
                            wk.append((nxt, iter([s for s in succs[nxt] if s in comp_set])))
                    except StopIteration:
                        order.append(node)
                        wk.pop()
            order.reverse()  # reverse-postorder
            return order

        # 7. Assemble per-layer cell lists; SCCs ordered internally by RPO.
        n_layers = (max(scc_layer) + 1) if scc_layer else 1
        layers = [[] for _ in range(n_layers)]
        for i, comp in enumerate(sccs):
            ordered = rpo_within_scc(comp) if len(comp) > 1 else comp
            layers[scc_layer[i]].extend(ordered)

        print(f"  Tarjan: {n_sccs} SCCs ({len(nontriv)} non-trivial), "
              f"{n_layers} layers, "
              f"{sum(len(c) for _,c in nontriv)} cells inside non-trivial SCCs")
        return layers

    def resolve_bits(self, bits):
        """Resolve a connection bit list to (bit_indices, const_mask, const_val).

        bit_indices[i] = LI index for bit i, or:
          -1 = constant 0 (or x/unknown)
          -2 = constant 1
        const_mask/const_val are still used for bits 0-63 (backward compat).
        For bits >= 64, the -1/-2 encoding in bit_indices is authoritative.
        """
        indices = []
        const_mask = 0
        const_val = 0
        for i, b in enumerate(bits):
            if isinstance(b, int):
                indices.append(b)
            elif b == "1":
                indices.append(-2)  # constant 1
                if i < 64:
                    const_mask |= (1 << i)
                    const_val |= (1 << i)
            elif b == "0":
                indices.append(-1)  # constant 0
                if i < 64:
                    const_mask |= (1 << i)
            else:  # "x" or other
                indices.append(-1)
                if i < 64:
                    const_mask |= (1 << i)
        return indices, const_mask, const_val

    def detect_identity(self, cd):
        """Detect identity ops: $mux with constant S that always selects one input.
        Returns (is_identity, id_src) where id_src=0 means A, id_src=1 means B."""
        if cd["type"] != "$mux":
            return 0, 0
        conns = cd["connections"]
        s_bits = conns.get("S", [])
        if len(s_bits) != 1:
            return 0, 0
        # S must be a constant (not a net)
        b = s_bits[0]
        if isinstance(b, int):
            return 0, 0  # S is a net, not constant
        # S is constant "0" or "1"
        if b == "0":
            return 1, 0  # always selects A
        elif b == "1":
            return 1, 1  # always selects B
        return 0, 0

    def generate_c_header(self, layers, path):
        """Generate oim_bitlvl.h with all OIM data.
        Compressed format:
        - Packed flags byte (is_identity, id_src, a_signed, b_signed)
        - NuRange arrays retained for forward-compatibility
        - Const masks only emitted when non-zero
        """
        L = []
        L.append("/* oim_bitlvl.h — Auto-generated bit-level OIM for RTeAAL tsim */")
        L.append("/* NU (N-rank Unrolling) + Identity Elision enabled */")
        L.append("/* Compressed format: packed flags, optimized layout */")
        L.append("#ifndef OIM_BITLVL_H")
        L.append("#define OIM_BITLVL_H")
        L.append("#include <stdint.h>")
        L.append("")

        L.append(f"#define LI_SIZE     {self.li_size}")
        L.append(f"#define N_LAYERS    {len(layers)}")
        n_dffs = len(self.dffs)
        L.append(f"#define N_DFFS      {n_dffs}")
        n_mems = len(self.mems)
        L.append(f"#define N_MEMS      {n_mems}")
        total_ops = sum(len(layer) for layer in layers)
        L.append(f"#define TOTAL_OPS   {total_ops}")
        L.append(f"#define N_CELL_TYPES {max(CELL_TYPES.values()) + 1}")

        # Robust reset-bit lookup: any top port whose name looks like a
        # synchronous- or asynchronous-low reset (rst_ni, rstn, rst_n, ...).
        rst_bit = -1
        for pn, pd in self.ports.items():
            lname = pn.lower()
            if lname in ("rst_ni", "rstn_i", "rst_n", "rstn", "rst_n_i", "aresetn"):
                bits = pd["bits"]
                if bits and isinstance(bits[0], int):
                    rst_bit = bits[0]
                    break
        if rst_bit < 0:
            # Fallback: second port bit (matches legacy hardcoded li[3])
            rst_bit = 3
        L.append(f"#define RST_NI_BIT  {rst_bit}")
        L.append("")

        # Cell type enum
        L.append("enum CellType {")
        for name, val in sorted(CELL_TYPES.items(), key=lambda x: x[1]):
            L.append(f"    CT_{name[1:].upper():16s} = {val},")
        L.append("};")
        L.append("")

        # NuRange struct for per-layer per-type ranges
        L.append("typedef struct { uint16_t start; uint16_t end; } NuRange;")
        L.append("")

        # Operation struct — compressed layout
        L.append("typedef struct {")
        L.append("    uint8_t  cell_type;")
        L.append("    uint32_t a_start; uint16_t a_len;   /* indices into bit_indices[] for port A */")
        L.append("    uint32_t b_start; uint16_t b_len;   /* port B */")
        L.append("    uint32_t s_start; uint16_t s_len;   /* port S (for mux/pmux) */")
        L.append("    uint32_t y_start; uint16_t y_len;   /* output port Y */")
        L.append("    uint64_t a_const_mask, a_const_val;  /* constant bits for A */")
        L.append("    uint64_t b_const_mask, b_const_val;  /* constant bits for B */")
        L.append("    uint64_t s_const_mask, s_const_val;  /* constant bits for S */")
        L.append("    uint8_t  a_signed, b_signed;")
        L.append("    uint16_t s_width;  /* for $pmux */")
        L.append("    uint8_t  is_identity;  /* 1 if this op is a pass-through */")
        L.append("    uint8_t  id_src;       /* 0=A, 1=B — which input to copy for identity */")
        L.append("} OimOp;")
        L.append("")

        L.append("typedef struct {")
        L.append("    uint32_t d_bits_start; uint16_t d_bits_len;  /* D input bit indices */")
        L.append("    uint32_t q_bits_start; uint16_t q_bits_len;  /* Q output bit indices */")
        L.append("    uint64_t reset_val;")
        L.append("    int16_t  en_bit;    /* Enable bit in LI (-1 = always enabled) */")
        L.append("    uint8_t  en_pol;    /* Enable polarity (1=active-high) */")
        L.append("    uint8_t  has_arst;  /* 1 if $adff/$adffe (async reset honors global RST_NI_BIT) */")
        L.append("} OimDff;")
        L.append("")
        L.append("/* OimMem — Yosys $mem_v2 cell kept intact (not lowered to FFs).")
        L.append(" * Read/Write port slices live in oim_bit_indices[]. Storage is held in")
        L.append(" * a parallel TsimMemState side-table, NOT in LI[]. RD_DATA bits in LI[]")
        L.append(" * are driven by tsim_mem_cycle() each clock edge. */")
        L.append("typedef struct {")
        L.append("    uint32_t size;            /* number of rows (depth)             */")
        L.append("    uint16_t width;           /* bits per row                       */")
        L.append("    uint8_t  abits;           /* address width (bits)               */")
        L.append("    uint8_t  n_rd_ports;")
        L.append("    uint8_t  n_wr_ports;")
        L.append("    /* Single-port slice indices (this design's memories are 1R1W or 1R0W).")
        L.append("     * For multi-port memories the slices concatenate ports along the wide axis. */")
        L.append("    uint32_t rd_clk_start;    /* RD_CLK[0..n_rd_ports-1] bit indices */")
        L.append("    uint16_t rd_clk_len;")
        L.append("    uint32_t rd_en_start;     /* RD_EN[0..n_rd_ports-1]              */")
        L.append("    uint16_t rd_en_len;")
        L.append("    uint32_t rd_addr_start;   /* RD_ADDR concatenated (abits * n_rd_ports) */")
        L.append("    uint16_t rd_addr_len;")
        L.append("    uint32_t rd_data_start;   /* RD_DATA concatenated (width * n_rd_ports) */")
        L.append("    uint16_t rd_data_len;")
        L.append("    uint32_t wr_clk_start;")
        L.append("    uint16_t wr_clk_len;")
        L.append("    uint32_t wr_en_start;     /* WR_EN: width bits per port (byte-enable mask) */")
        L.append("    uint16_t wr_en_len;")
        L.append("    uint32_t wr_addr_start;")
        L.append("    uint16_t wr_addr_len;")
        L.append("    uint32_t wr_data_start;")
        L.append("    uint16_t wr_data_len;")
        L.append("    const uint8_t *init;      /* depth * ceil(width/8) bytes, LSB-first; NULL = no init */")
        L.append("    const char    *name;")
        L.append("} OimMem;")
        L.append("")

        # Build flat bit-index array and ops
        all_bit_indices = []  # flat array of int16_t (LI indices, -1 for const)
        all_ops = []  # list of OimOp dicts per layer
        layer_sizes = []

        total_identity = 0
        identity_by_type = defaultdict(int)

        for layer in layers:
            layer_ops = []
            for cn in layer:
                cd = self.combs[cn]
                ct = CELL_TYPES[cd["type"]]
                conns = cd["connections"]
                params = cd.get("parameters", {})

                # Resolve each port
                def resolve_port(port_name):
                    if port_name in conns:
                        bits = conns[port_name]
                        indices, cmask, cval = self.resolve_bits(bits)
                        start = len(all_bit_indices)
                        all_bit_indices.extend(indices)
                        return start, len(indices), cmask, cval
                    return 0, 0, 0, 0

                a_s, a_l, a_cm, a_cv = resolve_port("A")
                b_s, b_l, b_cm, b_cv = resolve_port("B")
                s_s, s_l, s_cm, s_cv = resolve_port("S")
                y_s, y_l, y_cm, y_cv = resolve_port("Y")

                a_signed = params.get("A_SIGNED", 0)
                b_signed = params.get("B_SIGNED", 0)
                s_width = params.get("S_WIDTH", 0)
                # Parse binary string params from Yosys JSON
                if isinstance(a_signed, str): a_signed = int(a_signed, 2)
                if isinstance(b_signed, str): b_signed = int(b_signed, 2)
                if isinstance(s_width, str): s_width = int(s_width, 2)

                # Detect identity ops
                is_id, id_src = self.detect_identity(cd)
                if is_id:
                    total_identity += 1
                    identity_by_type[cd["type"]] += 1

                layer_ops.append({
                    "ct": ct, "cn": cn,
                    "a_s": a_s, "a_l": a_l, "a_cm": a_cm, "a_cv": a_cv,
                    "b_s": b_s, "b_l": b_l, "b_cm": b_cm, "b_cv": b_cv,
                    "s_s": s_s, "s_l": s_l, "s_cm": s_cm, "s_cv": s_cv,
                    "y_s": y_s, "y_l": y_l,
                    "a_signed": a_signed, "b_signed": b_signed,
                    "s_width": s_width,
                    "is_identity": is_id, "id_src": id_src,
                })

            # NU: Sort ops within each layer by cell_type (stable sort)
            layer_ops.sort(key=lambda op: op["ct"])
            all_ops.append(layer_ops)
            layer_sizes.append(len(layer_ops))

        # Print identity statistics
        print(f"  Identity ops (const-S $mux): {total_identity} total")
        for tname, count in sorted(identity_by_type.items()):
            print(f"    {tname}: {count}")

        # DFFs
        dff_list = []
        for cn, cd in self.dffs.items():
            d_bits = cd["connections"].get("D", [])
            q_bits = cd["connections"].get("Q", [])
            d_indices, _, _ = self.resolve_bits(d_bits)
            q_indices, _, _ = self.resolve_bits(q_bits)

            d_start = len(all_bit_indices)
            all_bit_indices.extend(d_indices)
            q_start = len(all_bit_indices)
            all_bit_indices.extend(q_indices)

            # Reset value:
            #   $adff/$adffe : ARST_VALUE (async, honors global RST_NI_BIT)
            #   $sdff        : SRST_VALUE (sync; treated as async here — close enough
            #                  for our boot smoke test, since the synth flow uses
            #                  a single global rst_ni domain)
            #   $aldff       : AD port (async data) which may be constants
            #   $dff/$dffe   : no reset port (has_arst=0; ignored on global reset)
            ctype = cd["type"]
            params = cd.get("parameters", {})
            if ctype == "$aldff":
                ad_bits = cd["connections"].get("AD", [])
                arst_val = 0
                for i, b in enumerate(ad_bits):
                    if b == "1":
                        arst_val |= (1 << i)
            elif ctype == "$sdff":
                arst_val = parse_param_int(params.get("SRST_VALUE", "0"))
            else:
                arst_val = parse_param_int(params.get("ARST_VALUE", "0"))

            has_arst = 1 if ctype in ("$adff", "$adffe", "$sdff", "$aldff") else 0

            # Enable bit for $adffe / $dffe
            en_bits = cd["connections"].get("EN", [])
            en_bit = -1  # -1 = always enabled ($adff, $aldff, $dff, $sdff)
            en_pol = 1
            if en_bits and isinstance(en_bits[0], int):
                en_bit = en_bits[0]
            en_pol_param = params.get("EN_POLARITY", "1")
            if isinstance(en_pol_param, str):
                en_pol = 1 if en_pol_param.endswith("1") else 0
            elif isinstance(en_pol_param, int):
                en_pol = en_pol_param

            dff_list.append({
                "d_start": d_start, "d_len": len(d_indices),
                "q_start": q_start, "q_len": len(q_indices),
                "reset_val": arst_val, "cn": cn,
                "en_bit": en_bit, "en_pol": en_pol,
                "has_arst": has_arst,
            })

        # ---------------------------------------------------------------
        # Memory cells ($mem_v2): parse first so port slices land in
        # all_bit_indices BEFORE the flat array is written. Emission of
        # the OimMem records and init data happens later (after DFFs).
        # ---------------------------------------------------------------
        mem_list = []
        for mi, (cn, cd) in enumerate(self.mems.items()):
            params = cd.get("parameters", {})
            conns  = cd.get("connections", {})
            size   = parse_param_int(params.get("SIZE"))
            abits  = parse_param_int(params.get("ABITS"))
            width  = parse_param_int(params.get("WIDTH"))
            nrd    = parse_param_int(params.get("RD_PORTS"))
            nwr    = parse_param_int(params.get("WR_PORTS"))

            def resolve_conn(name):
                bits = conns.get(name, [])
                if not bits:
                    return 0, 0
                indices, _, _ = self.resolve_bits(bits)
                start = len(all_bit_indices)
                all_bit_indices.extend(indices)
                return start, len(indices)

            rd_clk_s,  rd_clk_l  = resolve_conn("RD_CLK")
            rd_en_s,   rd_en_l   = resolve_conn("RD_EN")
            rd_addr_s, rd_addr_l = resolve_conn("RD_ADDR")
            rd_data_s, rd_data_l = resolve_conn("RD_DATA")
            wr_clk_s,  wr_clk_l  = resolve_conn("WR_CLK")
            wr_en_s,   wr_en_l   = resolve_conn("WR_EN")
            wr_addr_s, wr_addr_l = resolve_conn("WR_ADDR")
            wr_data_s, wr_data_l = resolve_conn("WR_DATA")

            # INIT: SIZE*WIDTH MSB-first binary string or "x".
            # Pack to a row-major, LSB-first byte array.
            init_param = params.get("INIT", "x")
            init_bytes = None
            if isinstance(init_param, str) and init_param not in ("", "x"):
                bits = init_param[::-1]              # → LSB-first across whole string
                bytes_per_row = (width + 7) // 8
                ba = bytearray(size * bytes_per_row)
                for row in range(size):
                    base = row * width
                    for bit in range(width):
                        idx = base + bit
                        if idx < len(bits) and bits[idx] == "1":
                            ba[row * bytes_per_row + (bit // 8)] |= 1 << (bit % 8)
                init_bytes = bytes(ba)

            mem_list.append({
                "idx": mi, "name": cn, "size": size, "width": width, "abits": abits,
                "nrd": nrd, "nwr": nwr,
                "rd_clk_s": rd_clk_s, "rd_clk_l": rd_clk_l,
                "rd_en_s":  rd_en_s,  "rd_en_l":  rd_en_l,
                "rd_addr_s": rd_addr_s, "rd_addr_l": rd_addr_l,
                "rd_data_s": rd_data_s, "rd_data_l": rd_data_l,
                "wr_clk_s": wr_clk_s, "wr_clk_l": wr_clk_l,
                "wr_en_s":  wr_en_s,  "wr_en_l":  wr_en_l,
                "wr_addr_s": wr_addr_s, "wr_addr_l": wr_addr_l,
                "wr_data_s": wr_data_s, "wr_data_l": wr_data_l,
                "init": init_bytes,
            })

        print(f"  Memories: {len(mem_list)}")
        for m in mem_list:
            init_kib = (len(m["init"]) / 1024.0) if m["init"] else 0.0
            print(f"    {m['name'][:60]:60s} size={m['size']} W={m['width']} "
                  f"rd={m['nrd']} wr={m['nwr']} init={init_kib:.1f}KiB")

        print(f"  Flat bit-index array: {len(all_bit_indices)} entries")
        print(f"  Total ops: {total_ops} across {len(layers)} layers")

        # Emit flat bit-index array (compressed: 40 entries per line)
        L.append(f"/* Flat bit-index array ({len(all_bit_indices)} entries) */")
        L.append(f"static const int32_t oim_bit_indices[{len(all_bit_indices)}] = {{")
        for i in range(0, len(all_bit_indices), 40):
            chunk = all_bit_indices[i:i+40]
            L.append(",".join(str(b) for b in chunk) + ",")
        L.append("};")
        L.append("")

        # Emit layer sizes
        L.append(f"static const uint16_t oim_layer_sizes[N_LAYERS] = {{")
        L.append("    " + ", ".join(str(s) for s in layer_sizes))
        L.append("};")
        L.append("")

        # Emit ops per layer
        for li, layer_ops in enumerate(all_ops):
            L.append(f"static const OimOp oim_layer_{li}[{len(layer_ops)}] = {{")
            for op in layer_ops:
                L.append(f"    {{ {op['ct']}, "
                         f"{op['a_s']},{op['a_l']}, "
                         f"{op['b_s']},{op['b_l']}, "
                         f"{op['s_s']},{op['s_l']}, "
                         f"{op['y_s']},{op['y_l']}, "
                         f"0x{op['a_cm']:x},0x{op['a_cv']:x}, "
                         f"0x{op['b_cm']:x},0x{op['b_cv']:x}, "
                         f"0x{op['s_cm']:x},0x{op['s_cv']:x}, "
                         f"{op['a_signed']},{op['b_signed']}, "
                         f"{op['s_width']}, "
                         f"{op['is_identity']},{op['id_src']} }},")
            L.append("};")

        L.append("")
        L.append(f"static const OimOp* const oim_layers[N_LAYERS] = {{")
        for li in range(len(all_ops)):
            L.append(f"    oim_layer_{li},")
        L.append("};")
        L.append("")

        # NuRange arrays eliminated — ops already sorted by cell_type within layers.
        # The evaluator iterates sequentially, benefiting from NU sort order
        # for branch prediction without needing explicit range lookups.
        L.append("")

        # DFFs
        L.append(f"static const OimDff oim_dffs[N_DFFS] = {{")
        for dff in dff_list:
            L.append(f"    {{ {dff['d_start']},{dff['d_len']}, "
                     f"{dff['q_start']},{dff['q_len']}, "
                     f"0x{dff['reset_val']:x}, "
                     f"{dff['en_bit']},{dff['en_pol']},{dff['has_arst']} }},  "
                     f"/* {dff['cn'][:50]} */")
        L.append("};")
        L.append("")
        # Full hierarchical DFF instance names (for per-FF Welch t RTL traceback)
        L.append(f"static const char* const oim_dff_names[N_DFFS] = {{")
        for dff in dff_list:
            safe = dff['cn'].replace('\\', '_').replace('"', '_')
            L.append(f'    "{safe}",')
        L.append("};")
        L.append("")
        # Per-DFF Q widths (bits per cell) — used by per-bit CSV expansion
        L.append(f"static const uint8_t oim_dff_widths[N_DFFS] = {{")
        for dff in dff_list:
            L.append(f"    {dff['q_len']},")
        L.append("};")
        L.append("")
        # Total DFF Q bits across all cells (for per-bit Welch t)
        total_dff_bits = sum(dff['q_len'] for dff in dff_list)
        L.append(f"#define N_DFF_BITS  {total_dff_bits}")
        L.append("")
        # RTL-named DFFs: invert netname.bits to get RTL signal name + bit offset
        # for each Q bit, then compress per cell into "<name>[hi:lo]" form.
        # Uses yosys's preserved netnames (which carry the original SV signal name).
        bit_to_nets = {}
        for nn, nd in self.ports.items():
            pass  # ports handled separately below
        # Re-walk netnames from the raw module data to build bit_to_nets
        # (self has self.cells but not self.netnames; pull from outer loader)
        # We piggyback on the JSON we already loaded — but BitLevelOIM doesn't
        # store netnames. Re-fetch via the module dict held during init.
        netnames = self._raw_module_netnames if hasattr(self, '_raw_module_netnames') else {}
        for nn, nd in netnames.items():
            bits = nd.get('bits', [])
            for off, b in enumerate(bits):
                if isinstance(b, int):
                    bit_to_nets.setdefault(b, []).append((nn, off))
        L.append(f"static const char* const oim_dff_rtl_names[N_DFFS] = {{")
        per_bit_names = []  # accumulate per-Q-bit names for the per-bit array
        for dff in dff_list:
            cn = dff['cn']
            cd = self.cells.get(cn, self.dffs.get(cn, {}))
            qb = cd.get('connections', {}).get('Q', [])
            # Resolve each Q bit to its non-internal netname
            base_offs = {}  # base_name -> list of (bit_pos_within_cell, off_within_net)
            per_bit_for_cell = []  # length = len(qb)
            for pos, b in enumerate(qb):
                if not isinstance(b, int):
                    per_bit_for_cell.append("(unresolved)")
                    continue
                cands = [(n, off) for n, off in bit_to_nets.get(b, []) if not n.startswith('$')]
                if cands:
                    n, off = cands[0]
                    base_offs.setdefault(n, []).append((pos, off))
                    per_bit_for_cell.append(f"{n}[{off}]")
                else:
                    per_bit_for_cell.append("(unresolved)")
            if base_offs:
                base = max(base_offs, key=lambda k: len(base_offs[k]))
                offs = [o for _, o in base_offs[base]]
                if offs:
                    rng = f"[{max(offs)}:{min(offs)}]" if max(offs) != min(offs) else f"[{offs[0]}]"
                else:
                    rng = ""
                rtl = f"{base}{rng}"
            else:
                rtl = "(unresolved)"
            safe = rtl.replace('\\', '_').replace('"', '_')
            L.append(f'    "{safe}",')
            per_bit_names.extend(per_bit_for_cell)
        L.append("};")
        L.append("")
        # Per-Q-bit RTL names — flat across all DFFs in oim_dffs order
        L.append(f"static const char* const oim_dff_bit_rtl_names[N_DFF_BITS] = {{")
        for nm in per_bit_names:
            safe = nm.replace('\\', '_').replace('"', '_')
            L.append(f'    "{safe}",')
        L.append("};")
        L.append("")

        # ---------------------------------------------------------------
        # Memory cells: emit init data + OimMem records (parsing already
        # done above so the flat bit-index array is complete).
        # ---------------------------------------------------------------
        # Emit init byte arrays (only for memories that have INIT data).
        for m in mem_list:
            if m["init"] is None:
                continue
            L.append(f"/* {m['name'][:60]} init: {len(m['init'])} bytes */")
            L.append(f"static const uint8_t mem_init_{m['idx']}[{len(m['init'])}] = {{")
            chunk = []
            for i, b in enumerate(m["init"]):
                chunk.append(f"0x{b:02x}")
                if (i + 1) % 16 == 0:
                    L.append("    " + ",".join(chunk) + ",")
                    chunk = []
            if chunk:
                L.append("    " + ",".join(chunk) + ",")
            L.append("};")
            L.append("")

        # Emit OimMem array (always emit, even if empty, so #if checks compile).
        if mem_list:
            L.append(f"static const OimMem oim_mems[N_MEMS] = {{")
            for m in mem_list:
                init_ptr = f"mem_init_{m['idx']}" if m["init"] is not None else "0"
                # Truncate / escape the name for the C string.
                safe_name = m["name"].replace("\\", "_").replace('"', "_")[:80]
                L.append(f"    {{ {m['size']}, {m['width']}, {m['abits']}, "
                         f"{m['nrd']}, {m['nwr']}, "
                         f"{m['rd_clk_s']}, {m['rd_clk_l']}, "
                         f"{m['rd_en_s']}, {m['rd_en_l']}, "
                         f"{m['rd_addr_s']}, {m['rd_addr_l']}, "
                         f"{m['rd_data_s']}, {m['rd_data_l']}, "
                         f"{m['wr_clk_s']}, {m['wr_clk_l']}, "
                         f"{m['wr_en_s']}, {m['wr_en_l']}, "
                         f"{m['wr_addr_s']}, {m['wr_addr_l']}, "
                         f"{m['wr_data_s']}, {m['wr_data_l']}, "
                         f'{init_ptr}, "{safe_name}" }},')
            L.append("};")
            L.append("")
        else:
            L.append("/* No $mem_v2 cells in this design */")
            L.append("static const OimMem oim_mems[1] = {{0}};")
            L.append("")

        # Port map — for testbench to drive/read ports
        L.append("typedef struct { const char *name; uint8_t dir; uint16_t width; uint16_t bits_start; } OimPort;")
        port_list = []
        for pn, pd in sorted(self.ports.items()):
            bits = pd["bits"]
            indices, _, _ = self.resolve_bits(bits)
            start = len(all_bit_indices)
            # We already have bit indices in the main array for cells,
            # but ports need their own entries. Store them separately.
            port_list.append({
                "name": pn, "dir": 0 if pd["direction"] == "input" else 1,
                "width": len(bits), "bit_ids": [b if isinstance(b, int) else -1 for b in bits],
            })

        L.append(f"#define N_PORTS {len(port_list)}")
        # Emit port bit arrays
        for pi, p in enumerate(port_list):
            L.append(f"static const int16_t port_{pi}_bits[{len(p['bit_ids'])}] = "
                     f"{{ {', '.join(str(b) for b in p['bit_ids'])} }};")

        L.append(f"static const OimPort oim_ports[N_PORTS] = {{")
        for pi, p in enumerate(port_list):
            L.append(f'    {{ "{p["name"]}", {p["dir"]}, {p["width"]}, 0 }},')
        L.append("};")
        L.append("")

        # Helper: port bit arrays indexed by port index
        L.append("static const int16_t* const port_bit_arrays[N_PORTS] = {")
        for pi in range(len(port_list)):
            L.append(f"    port_{pi}_bits,")
        L.append("};")
        L.append("")

        L.append("#endif /* OIM_BITLVL_H */")

        with open(path, "w") as f:
            f.write("\n".join(L) + "\n")
        print(f"  Wrote: {path} ({len(L)} lines)")

    def generate_debug_json(self, layers, path):
        """Write summary JSON."""
        summary = {
            "li_size": self.li_size,
            "n_layers": len(layers),
            "layer_sizes": [len(l) for l in layers],
            "n_dffs": len(self.dffs),
            "n_combs": len(self.combs),
            "ports": {pn: {"dir": pd["direction"], "width": len(pd["bits"])}
                      for pn, pd in self.ports.items()},
        }
        with open(path, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"  Wrote: {path}")


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Generate bit-level OIM header from Yosys flat JSON.")
    ap.add_argument("--json", default="phase1_extract/tsim_pqc_only_flat.json",
                    help="Path to Yosys-flattened JSON (input)")
    ap.add_argument("--out-h", default="phase1_extract/oim_bitlvl.h",
                    help="Output C header path")
    ap.add_argument("--out-json", default="phase1_extract/oim_bitlvl.json",
                    help="Output debug JSON summary path")
    args = ap.parse_args()

    json_path = Path(args.json)
    out_h     = Path(args.out_h)
    out_json  = Path(args.out_json)

    if not json_path.exists():
        print(f"Error: {json_path} not found")
        sys.exit(1)

    print("=== tsim_mlkem_mldsa  Bit-Level OIM Generator ===")
    print(f"  Input JSON: {json_path}")
    data = load_json(json_path)
    mod = list(data["modules"].values())[0]

    builder = BitLevelOIM(mod)
    layers = builder.levelize()
    print(f"  Layers: {len(layers)}")
    for i, layer in enumerate(layers[:5]):
        print(f"    Layer {i}: {len(layer)} cells")
    if len(layers) > 5:
        print(f"    ... ({len(layers) - 5} more layers)")

    builder.generate_c_header(layers, out_h)
    builder.generate_debug_json(layers, out_json)
    print("=== Done ===")


if __name__ == "__main__":
    main()
