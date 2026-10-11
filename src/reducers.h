/* Extracted from the archived intervention driver. See provenance/sources.csv. */
#ifndef TIL_SCA_REDUCERS_H
#define TIL_SCA_REDUCERS_H
#include <math.h>
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

#endif
