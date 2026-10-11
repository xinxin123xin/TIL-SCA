// dom_and_tsim_top.sv — Minimal DOM-AND module-level top for tsim.h
//
// Instantiates sec_and (WPI Boolean DOM-AND gadget, WIDTH=32) and a
// serial result_reg that outputs share0 then share1 on consecutive cycles.
// This models the "serial result register reuse" observation path.
//
// Port list (all driven/read by C driver via tsim_write_port/tsim_read_port):
//   clk_i, rst_ni, en_i
//   a0_i[31:0], a1_i[31:0]   — shares of operand a
//   b0_i[31:0], b1_i[31:0]   — shares of operand b
//   rnd_i[31:0]              — fresh random word (same as AND_LATCH/REUSE)
//   result_o[31:0]           — serial output: cycle N = c0, cycle N+1 = c1
//   valid_o                  — 1 when result_o holds a valid share
//   share_sel_o              — 0 = outputting c0, 1 = outputting c1

`include "sec_and.sv"

module dom_and_tsim_top #(
  parameter WIDTH = 32
)(
  input  wire              clk_i,
  input  wire              rst_ni,
  input  wire              en_i,

  input  wire [WIDTH-1:0]  a0_i,
  input  wire [WIDTH-1:0]  a1_i,
  input  wire [WIDTH-1:0]  b0_i,
  input  wire [WIDTH-1:0]  b1_i,
  input  wire [WIDTH-1:0]  rnd_i,

  output wire [WIDTH-1:0]  result_o,
  output reg               valid_o,
  output reg               share_sel_o
);

  wire [WIDTH-1:0] c0, c1;
  wire             and_done;

  sec_and #(.WIDTH(WIDTH)) u_dom_and (
    .clk_i    (clk_i),
    .rst_ni   (rst_ni),
    .en_i     (en_i),
    .a0_i     (a0_i),
    .a1_i     (a1_i),
    .b0_i     (b0_i),
    .b1_i     (b1_i),
    .rnd_i    (rnd_i),
    .c0_o     (c0),
    .c1_o     (c1),
    .done_o   (and_done)
  );

  // Serial result register: on done, present c0 first cycle, c1 second.
  reg [WIDTH-1:0] result_reg;
  reg             phase;  // 0 = output c0, 1 = output c1

  always @(posedge clk_i or negedge rst_ni) begin
    if (!rst_ni) begin
      result_reg <= {WIDTH{1'b0}};
      valid_o    <= 1'b0;
      share_sel_o <= 1'b0;
      phase      <= 1'b0;
    end else begin
      if (and_done && phase == 1'b0) begin
        result_reg  <= c0;
        valid_o     <= 1'b1;
        share_sel_o <= 1'b0;
        phase       <= 1'b1;
      end else if (phase == 1'b1) begin
        result_reg  <= c1;
        valid_o     <= 1'b1;
        share_sel_o <= 1'b1;
        phase       <= 1'b0;
      end else begin
        valid_o     <= 1'b0;
      end
    end
  end

  assign result_o = result_reg;

endmodule
