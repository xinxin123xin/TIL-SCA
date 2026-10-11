// sec_add_tsim_top.sv — Wrapper for sec_add module-level tsim.h experiment
// WIDTH=13, N_STAGES=4 (ML-KEM Kogge-Stone)

`include "sec_add.sv"

module sec_add_tsim_top #(
  parameter WIDTH    = 13,
  parameter N_STAGES = 4
)(
  input  wire                          clk_i,
  input  wire                          rst_ni,
  input  wire                          start_i,

  input  wire [WIDTH-1:0]              a0_i, a1_i,
  input  wire [WIDTH-1:0]              b0_i, b1_i,

  input  wire [(2+3*N_STAGES)*WIDTH-1:0] rnd_i,

  output wire [WIDTH-1:0]             s0_o, s1_o,
  output wire                         done_o,

  // Debug observation
  output wire [WIDTH-1:0]             dbg_g0_o, dbg_g1_o,
  output wire [WIDTH-1:0]             dbg_p0_o, dbg_p1_o,
  output wire [3:0]                   dbg_state_o,
  output wire [2:0]                   dbg_stage_o,
  output wire                         dbg_and_done_o,
  output wire                         dbg_or_done_o,
  output wire                         dbg_and_en_o,
  output wire [WIDTH-1:0]             dbg_and_c0_o,
  output wire [WIDTH-1:0]             dbg_and_c1_o,
  output wire [WIDTH-1:0]             dbg_and_a0_o,
  output wire [WIDTH-1:0]             dbg_and_b0_o
);

  sec_add #(.WIDTH(WIDTH), .N_STAGES(N_STAGES)) u_dut (
    .clk_i     (clk_i),
    .rst_ni    (rst_ni),
    .start_i   (start_i),
    .a0_i      (a0_i),
    .a1_i      (a1_i),
    .b0_i      (b0_i),
    .b1_i      (b1_i),
    .rnd_i     (rnd_i),
    .s0_o      (s0_o),
    .s1_o      (s1_o),
    .done_o    (done_o),
    .dbg_g0_o  (dbg_g0_o),
    .dbg_g1_o  (dbg_g1_o),
    .dbg_p0_o  (dbg_p0_o),
    .dbg_p1_o  (dbg_p1_o),
    .dbg_state_o (dbg_state_o),
    .dbg_stage_o (dbg_stage_o),
    .dbg_and_done_o (dbg_and_done_o),
    .dbg_or_done_o  (dbg_or_done_o),
    .dbg_and_en_o   (dbg_and_en_o),
    .dbg_and_c0_o   (dbg_and_c0_o),
    .dbg_and_c1_o   (dbg_and_c1_o),
    .dbg_and_a0_o   (dbg_and_a0_o),
    .dbg_and_b0_o   (dbg_and_b0_o)
  );

endmodule
