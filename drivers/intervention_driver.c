/* Real OIM/tsim execution; identical stimulus and protocol in both RTL variants. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "oim_bitlvl.h"
#include "tsim.h"

typedef struct { uint32_t s[4]; } Rng;
static uint32_t rol(uint32_t x,int k){return (x<<k)|(x>>(32-k));}
static uint32_t next(Rng*r){uint32_t z=rol(r->s[0]+r->s[3],7)+r->s[0],t=r->s[1]<<9;
 r->s[2]^=r->s[0];r->s[3]^=r->s[1];r->s[1]^=r->s[2];r->s[0]^=r->s[3];r->s[2]^=t;r->s[3]=rol(r->s[3],11);return z;}
static void seed(Rng*r,uint32_t a,uint32_t b,uint32_t c,uint32_t d){r->s[0]=a;r->s[1]=b;r->s[2]=c;r->s[3]=d;for(int i=0;i<16;i++)next(r);}
static void derive(int tid,int campaign,uint32_t*w){
 uint32_t h0=0x811c9dc5u,h1=0x01000193u,h2=0xdeadbeefu,h3=0xcafebabeu;
 char buf[80];
 if(campaign==0)snprintf(buf,sizeof(buf),"TSIM-DOM-SOC|trace=%d",tid);
 else snprintf(buf,sizeof(buf),"TSIM-DOM-SOC|campaign=%d|trace=%d",campaign,tid);
 for(char*p=buf;*p;p++){h0=(h0^(uint8_t)*p)*0x01000193u;h1+=(uint8_t)*p*2654435761u;h2=(h2<<5)^((uint8_t)*p*40503u);h3=(h3>>3)+((uint8_t)*p<<7);}
 w[0]=h0^(uint32_t)tid;w[1]=h1^(uint32_t)tid*0x9e3779b9u;w[2]=h2^(uint32_t)tid*0x85ebca6bu;w[3]=h3^(uint32_t)tid*0xc2b2ae35u;
}
static int hw(uint32_t x){int n=0;while(x){x&=x-1;n++;}return n;}
static uint64_t hash32(uint64_t h,uint32_t x){for(int k=0;k<4;k++){h^=(x>>(8*k))&255;h*=1099511628211ULL;}return h;}
static uint32_t instr(int g,int f){return 0x7bu|((uint32_t)g<<12)|((uint32_t)f<<25);}
typedef struct{double n,m,m2;} Stats;
static void add(Stats*s,double x){double d=x-s->m;s->n++;s->m+=d/s->n;s->m2+=d*(x-s->m);}
typedef struct{double n,sx,sy,sx2,sy2,sxy,sx2y,sxy2,sx2y2;} Mixed;
static void mixed_add(Mixed*m,double x,double y){m->n++;m->sx+=x;m->sy+=y;m->sx2+=x*x;m->sy2+=y*y;m->sxy+=x*y;m->sx2y+=x*x*y;m->sxy2+=x*y*y;m->sx2y2+=x*x*y*y;}
static void mixed_z(Mixed*m,double mx,double my,double*mean,double*var){
 double sz=m->sxy-my*m->sx-mx*m->sy+m->n*mx*my;
 double sz2=m->sx2y2-2*my*m->sx2y-2*mx*m->sxy2+my*my*m->sx2+mx*mx*m->sy2+4*mx*my*m->sxy-2*mx*my*my*m->sx-2*mx*mx*my*m->sy+m->n*mx*mx*my*my;
 *mean=sz/m->n;*var=(sz2-sz*sz/m->n)/(m->n-1);
 if(*var<0&&*var>-1e-9)*var=0;
}
static double test(double mf,double vf,double nf,double mr,double vr,double nr){double se=sqrt(vf/nf+vr/nr);return se>0?(mf-mr)/se:NAN;}
static int dff_find(const char*name){int found=-1;size_t l=strlen(name);for(int i=0;i<N_DFFS;i++)if(strncmp(oim_dff_rtl_names[i],name,l)==0&&oim_dff_rtl_names[i][l]=='['){if(found!=-1||oim_dffs[i].q_bits_len!=32){fprintf(stderr,"Ambiguous DFF %s\n",name);exit(2);}found=i;}return found;}
static uint32_t read_dff(TsimState*s,int d){if(d<0)return 0;uint32_t v=0;for(int k=0;k<32;k++)v|=(uint32_t)s->li[oim_bit_indices[oim_dffs[d].q_bits_start+k]]<<k;return v;}

enum{F_SHARED,F_B0_WRITE,F_B1_WRITE,F_B0_INTERVAL,F_B1_INTERVAL,F_BUS_ENDPOINT,F_BUS_SUM,F_S0_HW,F_S1_HW,NF};
static const char*names[NF]={"result_reg_between_returns","share_return0_write","share_return1_write","share_return0_between_returns","share_return1_between_returns","return_bus_valid_endpoint","return_bus_adjacent_sum","share0_hw","share1_hw"};
typedef struct{
 TsimState st;int valid,instr,rs1,rs2,ready,rvalid,data,issue_ready,accept;
 int trace,capture;uint32_t prev_bus;int sum_bus;FILE*cycles;uint64_t protocol_hash;int issue_id;
 int dr,d0,d1;uint32_t prev_reg[3];int reg_hd_sum[3];
} Run;
static void tick(Run*r){tsim_cycle(&r->st);if(r->capture){uint32_t b=(uint32_t)tsim_read_port(&r->st,r->data);int h=hw(b^r->prev_bus);r->sum_bus+=h;
 uint32_t q[3]={read_dff(&r->st,r->dr),read_dff(&r->st,r->d0),read_dff(&r->st,r->d1)};int hd[3];for(int k=0;k<3;k++){hd[k]=hw(q[k]^r->prev_reg[k]);r->reg_hd_sum[k]+=hd[k];r->prev_reg[k]=q[k];}
 fprintf(r->cycles,"%d,%llu,%d,0x%08x,%d,0x%08x,%d,0x%08x,%d,0x%08x,%d\n",r->trace,(unsigned long long)r->st.cycle,(int)tsim_read_port(&r->st,r->rvalid),b,h,q[0],hd[0],q[1],hd[1],q[2],hd[2]);r->prev_bus=b;}}
static uint32_t issue(Run*r,uint32_t ins,uint32_t a,uint32_t b,int start_capture){
 if(!tsim_read_port(&r->st,r->issue_ready)){fprintf(stderr,"Not issue-ready trace=%d op=%d\n",r->trace,r->issue_id);exit(3);}
 uint64_t begin=r->st.cycle;
 tsim_write_port(&r->st,r->instr,ins);tsim_write_port(&r->st,r->rs1,a);tsim_write_port(&r->st,r->rs2,b);tsim_write_port(&r->st,r->valid,1);tsim_write_port(&r->st,r->ready,0);
 tsim_eval_comb(&r->st); /* Settle newly driven instruction before checking combinational accept. */
 if(!tsim_read_port(&r->st,r->accept)){fprintf(stderr,"Not accepted trace=%d\n",r->trace);exit(3);}tick(r);tsim_write_port(&r->st,r->valid,0);
 int guard=0;while(!tsim_read_port(&r->st,r->rvalid)){tick(r);if(++guard>5000){fprintf(stderr,"Timeout\n");exit(3);}}
 uint32_t result=(uint32_t)tsim_read_port(&r->st,r->data);
 r->protocol_hash=hash32(r->protocol_hash,(uint32_t)begin);r->protocol_hash=hash32(r->protocol_hash,(uint32_t)r->st.cycle);r->protocol_hash=hash32(r->protocol_hash,(uint32_t)r->issue_id++);
 if(start_capture){r->capture=1;r->prev_bus=result;r->sum_bus=0;r->prev_reg[0]=read_dff(&r->st,r->dr);r->prev_reg[1]=read_dff(&r->st,r->d0);r->prev_reg[2]=read_dff(&r->st,r->d1);fprintf(r->cycles,"%d,%llu,1,0x%08x,0,0x%08x,0,0x%08x,0,0x%08x,0\n",r->trace,(unsigned long long)r->st.cycle,result,r->prev_reg[0],r->prev_reg[1],r->prev_reg[2]);}
 if(r->issue_id==8)r->capture=0; /* Endpoint has been captured by the final tick. */
 tsim_write_port(&r->st,r->ready,1);tick(r);tsim_write_port(&r->st,r->ready,0);
 if(tsim_read_port(&r->st,r->rvalid)||!tsim_read_port(&r->st,r->issue_ready)){fprintf(stderr,"Retirement protocol error\n");exit(3);}
 return result;
}

int main(int argc,char**argv){
 int per=argc>1?atoi(argv[1]):10,campaign=argc>2?atoi(argv[2]):0,N=2*per;
 if(per<2||campaign<0||campaign>2)return 2;
 int dr=dff_find("result_reg"),d0=dff_find("share_return0_reg"),d1=dff_find("share_return1_reg");int split=d0>=0&&d1>=0;if(dr<0||(d0>=0)!=(d1>=0))return 2;
 printf("variant=%s N=%d per_class=%d seed=%d DFF_shared=%d DFF_bank0=%d DFF_bank1=%d\n",split?"split_storage":"baseline",N,per,campaign,dr,d0,d1);
 printf("cycle_snapshots: port returns zero outside ST_RESULT_READY; bus endpoint HD is distinct from adjacent-cycle sum\n");
 Stats stats[NF][2]={0};Mixed moments[2]={0};double*raw=calloc((size_t)N*NF,sizeof(double));
 FILE*csv=fopen("trace_audit.csv","w"),*cycles=fopen("return_bus_cycles.csv","w");if(!raw||!csv||!cycles)return 2;
 fprintf(csv,"trace_id,group,input_hash,output_hash,protocol_hash,cycles,x,y,share0,share1,identity_ok,functional_ok,protocol_ok");for(int f=0;f<NF;f++)fprintf(csv,",%s",names[f]);fprintf(csv,"\n");
 fprintf(cycles,"trace_id,cycle,result_valid,result_data,adjacent_hd,result_reg,result_reg_adjacent_hd,share_return0_reg,bank0_adjacent_hd,share_return1_reg,bank1_adjacent_hd\n");
 int functional=0,identity=0;uint64_t poolhash=14695981039346656037ULL;
 for(int tid=0;tid<N;tid++){
  int group=tid%2;uint32_t sw[4];derive(tid,campaign,sw);Rng rng;seed(&rng,sw[0],sw[1],sw[2],sw[3]);
  uint32_t x=0xffffffffu,y=0xffffffffu;if(group){x=next(&rng);y=next(&rng);}uint32_t x0=next(&rng),y0=next(&rng),x1=x^x0,y1=y^y0;
  uint32_t seedw[4]={next(&rng),next(&rng),next(&rng),next(&rng)};
  uint64_t ih=14695981039346656037ULL;uint32_t input[10]={x,y,x0,x1,y0,y1,seedw[0],seedw[1],seedw[2],seedw[3]};for(int k=0;k<10;k++)ih=hash32(ih,input[k]);poolhash=hash32(hash32(poolhash,(uint32_t)ih),(uint32_t)(ih>>32));
  Run r={0};r.trace=tid;r.cycles=cycles;r.protocol_hash=14695981039346656037ULL;r.dr=dr;r.d0=d0;r.d1=d1;
  r.valid=port_idx("x_issue_valid_i");r.instr=port_idx("x_issue_req_instr_i");r.rs1=port_idx("x_issue_req_rs1_i");r.rs2=port_idx("x_issue_req_rs2_i");r.ready=port_idx("x_result_ready_i");r.rvalid=port_idx("x_result_valid_o");r.data=port_idx("x_result_data_o");r.issue_ready=port_idx("x_issue_ready_o");r.accept=port_idx("x_issue_resp_accept_o");
  tsim_init(&r.st);int rst=port_idx("rst_ni");tsim_write_port(&r.st,rst,0);tsim_write_port(&r.st,port_idx("x_commit_valid_i"),0);tick(&r);tsim_write_port(&r.st,rst,1);tick(&r);
  uint32_t out[8];for(int k=0;k<4;k++)out[k]=issue(&r,instr(7,0),seedw[k],k,0);
  out[4]=issue(&r,instr(6,6),x0,y1,0);out[5]=issue(&r,instr(6,7),x1,y0,0);
  uint32_t before0=read_dff(&r.st,d0),before1=read_dff(&r.st,d1);
  out[6]=issue(&r,instr(6,3),x0&y0,out[4],1);
  uint32_t sr0=read_dff(&r.st,dr),b00=read_dff(&r.st,d0),b10=read_dff(&r.st,d1);
  if(split&&(b00!=out[6]||b10!=before1)){fprintf(stderr,"Bank0 write/hold mismatch\n");return 3;}
  out[7]=issue(&r,instr(6,3),x1&y1,out[5],0);
  uint32_t sr1=read_dff(&r.st,dr),b01=read_dff(&r.st,d0),b11=read_dff(&r.st,d1);
  int ok=(out[6]^out[7])==(x&y)&&out[6]==((x0&y0)^out[4])&&out[7]==((x1&y1)^out[5]);
  if(!ok||(split&&(b01!=b00||b11!=out[7]))||(!split&&(sr0!=out[6]||sr1!=out[7]))){fprintf(stderr,"Functional/register mismatch trace=%d\n",tid);return 1;}functional++;
  double values[NF]={(double)hw(sr0^sr1),NAN,NAN,NAN,NAN,(double)hw(out[6]^out[7]),(double)r.sum_bus,(double)hw(out[6]),(double)hw(out[7])};
  if(split){values[F_B0_WRITE]=hw(before0^b00);values[F_B1_WRITE]=hw(before1^b11);values[F_B0_INTERVAL]=hw(b00^b01);values[F_B1_INTERVAL]=hw(b10^b11);}
  int id=(int)values[F_BUS_ENDPOINT]==hw(x&y);if(!id)return 1;identity++;
  /* The return bus is zero between transactions: retire share0 then present share1. */
  if(r.sum_bus!=hw(out[6])+hw(out[7])){fprintf(stderr,"Unexpected bus window trace=%d sum=%d\n",tid,r.sum_bus);return 3;}
  if(r.reg_hd_sum[0]!=hw(sr0^sr1)||r.reg_hd_sum[1]!=hw(b00^b01)||r.reg_hd_sum[2]!=hw(b10^b11)){fprintf(stderr,"Intervening register changes trace=%d\n",tid);return 3;}
  uint64_t oh=14695981039346656037ULL;for(int k=0;k<8;k++)oh=hash32(oh,out[k]);
  fprintf(csv,"%d,%s,%016llx,%016llx,%016llx,%llu,0x%08x,0x%08x,0x%08x,0x%08x,1,1,1",tid,group?"random":"fixed",(unsigned long long)ih,(unsigned long long)oh,(unsigned long long)r.protocol_hash,(unsigned long long)r.st.cycle,x,y,out[6],out[7]);
  for(int f=0;f<NF;f++){raw[(size_t)tid*NF+f]=values[f];if(isfinite(values[f]))add(&stats[f][group],values[f]);fprintf(csv,",%.17g",values[f]);}fprintf(csv,"\n");
  mixed_add(&moments[group],values[F_S0_HW],values[F_S1_HW]);
 }
 fclose(csv);fclose(cycles);
 FILE*summary=fopen("statistics.csv","w");fprintf(summary,"variant,seed,N,feature,n_fixed,n_random,mean_fixed,mean_random,var_fixed,var_random,t_online,t_offline,abs_delta,status,threshold_crossed\n");double max_delta=0;
 for(int f=0;f<NF;f++){
  if(stats[f][0].n==0){fprintf(summary,"%s,%d,%d,%s,0,0,nan,nan,nan,nan,nan,nan,nan,not_present,0\n",split?"split_storage":"baseline",campaign,N,names[f]);continue;}
  double mean[2]={0},var[2]={0};for(int tid=0;tid<N;tid++)mean[tid%2]+=raw[(size_t)tid*NF+f];for(int g=0;g<2;g++)mean[g]/=per;
  for(int tid=0;tid<N;tid++){double d=raw[(size_t)tid*NF+f]-mean[tid%2];var[tid%2]+=d*d;}for(int g=0;g<2;g++)var[g]/=(per-1);
  double vf=stats[f][0].m2/(per-1),vr=stats[f][1].m2/(per-1);
  double ton=test(stats[f][0].m,vf,per,stats[f][1].m,vr,per),toff=test(mean[0],var[0],per,mean[1],var[1],per),delta=fabs(ton-toff);
  if(isfinite(ton)!=isfinite(toff)||(isfinite(delta)&&delta>1e-9)){fprintf(stderr,"Paired mismatch %s\n",names[f]);return 1;}if(isfinite(delta)&&delta>max_delta)max_delta=delta;
  fprintf(summary,"%s,%d,%d,%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%s,%d\n",split?"split_storage":"baseline",campaign,N,names[f],per,per,stats[f][0].m,stats[f][1].m,vf,vr,ton,toff,delta,isfinite(ton)?"finite":"undefined_zero_denominator",isfinite(ton)&&fabs(ton)>=4.5);
  printf("feature=%s meanF=%.6f meanR=%.6f t=%.9g paired_delta=%.3g\n",names[f],stats[f][0].m,stats[f][1].m,ton,delta);
 }
 /* Predefined output HW pair, Q=1: online nine-moment CP versus materialized offline. */
 double mx=(moments[0].sx+moments[1].sx)/N,my=(moments[0].sy+moments[1].sy)/N,zm[2],zv[2],om[2]={0},ov[2]={0};
 for(int g=0;g<2;g++)mixed_z(&moments[g],mx,my,&zm[g],&zv[g]);
 for(int tid=0;tid<N;tid++)om[tid%2]+=(raw[(size_t)tid*NF+F_S0_HW]-mx)*(raw[(size_t)tid*NF+F_S1_HW]-my);
 for(int g=0;g<2;g++)om[g]/=per;
 for(int tid=0;tid<N;tid++){double z=(raw[(size_t)tid*NF+F_S0_HW]-mx)*(raw[(size_t)tid*NF+F_S1_HW]-my),d=z-om[tid%2];ov[tid%2]+=d*d;}for(int g=0;g<2;g++)ov[g]/=(per-1);
 double cp_on=test(zm[0],zv[0],per,zm[1],zv[1],per),cp_off=test(om[0],ov[0],per,om[1],ov[1],per),cp_delta=fabs(cp_on-cp_off);
 if(isfinite(cp_on)!=isfinite(cp_off)||(isfinite(cp_delta)&&cp_delta>1e-9)){fprintf(stderr,"CP paired mismatch %.17g\n",cp_delta);return 1;}if(isfinite(cp_delta)&&cp_delta>max_delta)max_delta=cp_delta;
 fprintf(summary,"%s,%d,%d,output_share_hw_centered_product,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%s,%d\n",split?"split_storage":"baseline",campaign,N,per,per,zm[0],zm[1],zv[0],zv[1],cp_on,cp_off,cp_delta,isfinite(cp_on)?"finite":"undefined_zero_denominator",isfinite(cp_on)&&fabs(cp_on)>=4.5);
 printf("feature=output_share_hw_centered_product Q=1 moments=18 meanF=%.6f meanR=%.6f t=%.9g paired_delta=%.3g\n",zm[0],zm[1],cp_on,cp_delta);
 fclose(summary);free(raw);
 printf("PASS functional=%d/%d protocol=%d/%d endpoint_identity=%d/%d input_pool_hash=%016llx max_paired_delta=%.3g\n",functional,N,functional,N,identity,N,(unsigned long long)poolhash,max_delta);return 0;
}
