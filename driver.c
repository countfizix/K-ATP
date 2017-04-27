#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "atp.h"
#include "atp.c"
#include "cblock.h"

/* global variable declarations */

/* double F[N];
double Y[N];
long int *np;
double *xp; */
 
double current[C],vset1;
double gsyn,iapp,epsilon;
double kact;

typedef struct stim_c stim_c;

struct stim_c {
   int time;
   double value;
};
 
int main(int argc, char *argv[])
{
double state[N];
double fstate[N];
double rtol[N],atol[N];

int itol,ijac,mljac,mujac;
int imas,mlmas,mumas;
int iout,lwork,liwork,lrcont;
int iwork[LIWORK];
double work[LWORK];
extern int out_(), dummy_out_();
int n,idid;

double x,xend,h;
FILE *fp;
int i,j;
fp=fopen("interval.data","w");
fclose(fp);

int sz = 0;
int ch = 0;


fp = fopen(STIM_FILE,"r");
while(!feof(fp)){
   ch = fgetc(fp);
   if(ch=='\n') sz++;  
}
  
//sz = ftell(fp)+1; //return line number

//printf("%d\n", sz);
rewind(fp);//

//const int stim_len = sz;

struct stim_c stim_s[sz];  

for(i=0;i<sz; i++){
   fscanf(fp, "%d", &stim_s[i].time);
   fscanf(fp, "%lf", &stim_s[i].value); //read stimulus info
   //printf("%d %e\n", stim_s[i].time, stim_s[i].value);
}

fclose(fp);


n=N;
//iapp=0.1*I_APP/(3.1415926538979*D_S*L_S);
epsilon = EPSILON;
for(i=0;i<7;i++) {iwork[i]=0;
                  work[i]=0.0;}
for(i=1;i<N;i++) {rtol[i]=1.0e-6;//9
                  atol[i]=1.0e-8;}//15
atol[V1]=1.0e-6;
rtol[V1]=1.0e-4;
//atol[A]=1.0e-8;
//rtol[A]=1.0e-6;
atol[Ca0]=1.0e-12;
rtol[Ca0]=1.0e-8;
//atol[Ca1]=1.0e-12;
//rtol[Ca1]=1.0e-8;
atol[CB0]=1.0e-12;
rtol[CB0]=1.0e-8;
//atol[CB1]=1.0e-8;
//rtol[CB1]=1.0e-6;
//atol[CA]=1.0e-8;
//rtol[CA]=1.0e-6;
//atol[NA]=1.0e-8;
//rtol[NA]=1.0e-6;

iwork[1] = 1000000000;
iwork[2] = 1000000;
iwork[3] = 1;
itol=1;
ijac=0;
mljac=n;
mujac=0;
imas=0;
mlmas=n;
mumas=0;
lwork=LWORK;
liwork=LIWORK;
lrcont=LRCONT;
iout=1;

h=1e-4;
x=START_TIME;
scan_(state);
deriv_(&n,&x,state,fstate);
xend = x;

/*xend = WARMUP;
iapp = 0.0;
radau5_(&n,deriv_,&x,state,&xend,&h,
   rtol,atol,&itol,
   dummy,&ijac,&mljac,&mujac,
   mas,&imas,&mlmas,&mumas,
   dummy_out_,0,
   work,&lwork,iwork,&liwork,&lrcont,&idid);

x = START_TIME;*/

if (stim_s[0].time > 0){
   xend=stim_s[0].time;
   iapp = 0.0;
   radau5_(&n,deriv_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        out_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);
}

for(i=1;i<sz;i++){
   if(x<ENDTIME){
      if (stim_s[i].time < ENDTIME) xend = stim_s[i].time;
      else xend = ENDTIME;
      iapp = stim_s[i-1].value;
      radau5_(&n,deriv_,&x,state,&xend,&h,
           rtol,atol,&itol,
           dummy,&ijac,&mljac,&mujac,
           mas,&imas,&mlmas,&mumas,
           out_,&iout,
           work,&lwork,iwork,&liwork,&lrcont,&idid);
   }
}//end for

if(x<ENDTIME){
   xend = ENDTIME;
   iapp = stim_s[sz-1].value;
   radau5_(&n,deriv_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        out_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);
   j = 0;
   /*while(idid < 0){
     if(j>500) break;
     j++;
     for(i=0;i<LWORK;i++) {
                    work[i]=0.0;}
     for(i=0;i<LIWORK;i++) {
                    iwork[i]=0.0;}

     radau5_(&n,deriv_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        out_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);
   }*/
}

dump_(state);
fflush(stdout);
}

