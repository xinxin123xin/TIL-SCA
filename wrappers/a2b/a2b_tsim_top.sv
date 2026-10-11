// a2b_tsim_top.sv — Wrapper for a2b_core module-level tsim.h experiment
`include "a2b_core.sv"

module a2b_tsim_top (
  input  wire        clk_i,
  input  wire        rst_ni,
  input  wire        start_i,
  input  wire [31:0] share0_i,
  input  wire [31:0] share1_i,
  input  wire [31:0] random_i,
  input  wire        wide_mode_i,
  input  wire [407:0] rnd_bypass_i,
  output wire [31:0] share0_o,
  output wire [31:0] share1_o,
  output wire        done_o,
  output wire [3:0]  dbg_state_o,
  output wire [23:0] dbg_raw0_o,
  output wire [23:0] dbg_raw1_o,
  output wire [23:0] dbg_trial0_o,
  output wire [23:0] dbg_trial1_o
);

  a2b_core u_dut (
    .clk_i      (clk_i),
    .rst_ni     (rst_ni),
    .start_i    (start_i),
    .share0_i   (share0_i),
    .share1_i   (share1_i),
    .random_i   (random_i),
    .wide_mode_i(wide_mode_i),
    .rnd_bypass_i(rnd_bypass_i),
    .share0_o   (share0_o),
    .share1_o   (share1_o),
    .done_o     (done_o),
    .dbg_state_o(dbg_state_o),
    .dbg_raw0_o (dbg_raw0_o),
    .dbg_raw1_o (dbg_raw1_o),
    .dbg_trial0_o(dbg_trial0_o),
    .dbg_trial1_o(dbg_trial1_o)
  );

endmodule
