//4 CALCIUM compartments with dynamic buffering

#ifndef ATP_C_
#define ATP_C_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "atp.h"

#ifndef EFUN_
#define EFUN_
double efun(z) 
double z;
{
	if (fabs(z) < 1e-4) {
	return( 1 - z/2);
	}else{
		return( z/(exp(z) - 1));
	}
}
#endif

float rand_gauss (void) {
  float v1,v2,s;

  do {
    v1 = 2.0 * ((float) rand()/RAND_MAX) - 1;
    v2 = 2.0 * ((float) rand()/RAND_MAX) - 1;

    s = v1*v1 + v2*v2;
  } while ( s >= 1.0 );

  if (s == 0.0)
    return 0.0;
  else
    return (v1*sqrt(-2.0 * log(s) / s));
}
 
 
double  gaussian(v,a,b,c,d) 
double  v,a,b,c,d ;
{
double arg;
arg = pow(((v-c)/b),2.0);
if ((arg>50.0) || (arg<-50.0))
{if (arg>50.0) return(d);
else  return(d+a);}
else return  (d + a*exp(-arg));
}

#ifndef BOLTZ_
#define BOLTZ_
double boltz(v,half,slope)
double v,half,slope;
{
double arg;
arg = -(v-half)/slope;
if ((arg>50.0) || (arg<-50.0))
{if (arg>50.0) return(0.0);
else  return(1.0);}
else return(1.0/(1.0 + exp(arg)));}
#endif

double tanhsig(v,half,slope) //sigmoidal tanh function
double v, half, slope;
{
double arg;
arg = (v-half)/slope;
if(arg > 50.0) return 1.0;
else if(arg < -50.0) return 0.0;
else return 0.5*(1.0+tanh(arg));
}



int deriv_(np,xp,Y,F)
double *F,*Y;
int *np;
double *xp;
{
#ifndef G_NMDA
extern double G_NMDA;
#endif

#ifndef atau
extern double atau;
#endif

extern double current[C]; 
int el;
double time_;
double iapp; //=I_APP;
double anm,bnm;
double minf,hinf,htau,ntau,mtau,ninf,arg,dlinf,dltau,nminf;
//atau
time_ = *xp;
double fadp,cadp,sinf;
el = *np;


//optional squarewave pulses
if (time_<2000 || time_>2200){ 
  iapp = -G_GABA*(Y[V1]+80.0);
}

if (time_>7000 && time_<7200){
  iapp = 0;
}

//optional tonic gabaergic stimulus
iapp = -G_GABA*(Y[V1]+80.0);

//HH channels
minf= boltz(Y[V1],-28.0907,8.7264); //sodium activation -28 9
hinf= boltz(Y[V1],-54.0,-8.7665); //sodium inactivation -54
ninf= boltz(Y[V1],-20.0,7.0); //fast potassium


htau = 56.0*(boltz(Y[V1],-21.0,-4.5) - boltz(Y[V1],-41.0,-2.0))+1.0;
mtau= (0.01 + 1.0/((15.6504+0.4043*Y[V1])/(1.0-exp(-19.565-0.50542*Y[V1])) +3.0212*exp(-7.463e-3*Y[V1])));
if(Y[V1]>-60.0 ) ntau= 1.0 + 19.0*exp(-pow(20*log(1.0 +0.05*(Y[V1]+40.0)),2.0)/300.0); //1 19 40 the 20 does all the work
else ntau=1.0; //1.0


// L-type
dlinf= boltz(Y[V1],dloff,dlslope); //-45 7.5 //L type Calcium
dltau= 1.0/(-0.020876*(Y[V1]+39.726)/(exp(-(Y[V1]+39.726)/4.711)-1) + 0.19444*exp(-(Y[V1]+15.338)/224.21));

//NMDA

anm= 1.0/(1.0+(1.4/3.57)*exp(-0.062*Y[V1])); //NMDA activation
bnm= NMC+(1.0-NMC)/(1.0+(1.2/NMOFF)*exp(-Y[V1]/NMSLOPE)); //NMDA activation
nminf = NMDATYPE*anm+(1-NMDATYPE)*bnm; //NMDA type switches between broad and sharp activation types, for paper NMDATYPE = 0


current[I_NMDA_CA] = FNC*G_NMDA*nminf*(Y[V1]- E_CA); //calcium component of NMDA


F[Ca0] = -2.0e6*fsca*(current[I_LCa] + current[I_CaL] +current[I_CAP] + current[I_NMDA_CA])/(D_S*0.0001*FARADAY); //calcium currents


F[CB0] = -KU*Y[CB0] + KB*Y[Ca0]*(1.0-Y[CB0]); // Y[CB0] fraction of buffering stuff bound to calcium


F[Ca0] = F[Ca0] - TBUFF*F[CB0];  //dynmaic buffering added to calcium dynamics



current[I_NA1] = G_NA*pow(Y[M1],3.0)*Y[H1]*(Y[V1] - E_NA);
current[I_K1] = G_K*pow(Y[N1],3.0)*(Y[V1] - E_K);

current[I_L1] =  G_L*(Y[V1] - E_L);
current[I_LCa] =G_LCa*(Y[V1]- E_CA);

current[I_CaL] =G_CaL*Y[DL]*(Y[V1]- E_CA);

current[I_NMDA] = G_NMDA*nminf*(Y[V1]);



current[I_CAP] = pow(Y[Ca0],CSLOPE)/(pow(Y[Ca0],CSLOPE)+pow(CHALF,CSLOPE)); //calcium pump activation (non electrogenic)


//K_ATP activation and ADP dynamics

F[A] = BEE*current[I_CAP] - KNST*Y[A]/(BASE+Y[A]); //unscaled I_CAPMAX without loss of generality


//The MM component of the ATP dynamics (KNST*Y[A]...) is important.  At high enough levels of calcium pump/sodium activity, the mitochondria simply cannot keep up

sinf = Y[A]/(S_HALF+Y[A]);  //sensitivity of K-ATP channel to ADP (Y[A])  

//sinf = Y[A]/(S_HALF+Y[A]); //this needs to be changed such that it has a positive slope

//fadp = 1.0/pow(1.0+S_HALF/Y[A],fpower);
//cadp = 1.0 - fadp;
//sinf = 1.0 - pow(cadp,S_SLOPE);

//ATP variations assumed to be unimportant (for now)s

current[I_CAP] = current[I_CAP]*I_CAPMAX; //scaled calcium pump current


current[I_KATP] = G_KATP*pow(sinf,spower)*(Y[V1]- E_K);

F[V1] = 1000*(iapp - current[I_NA1] - current[I_K1] - current[I_L1] - current[I_LCa] -current[I_CaL] -current[I_NMDA]- current[I_KATP])/CM +SCALE*rand_gauss();



F[M1] = (minf-Y[M1])/mtau;
F[H1] = (hinf-Y[H1])/htau;
F[N1] = (ninf-Y[N1])/ntau;
F[DL] = (dlinf-Y[DL])/dltau;
//F[CA] = -2e6*fsca*(current[I_LCa] + current[I_CaL] +current[I_CAP] + current[I_NMDA_CA] + current[I_CaN])/(D_S*0.0001*FARADAY);

//printf("%e %e\n", F[V1], F[R1]);

return 0;
   }


void scan_(Y) 
double Y[N];
{FILE *fopen(),*sp;
int i;
#ifdef BIFURCATION
sp = fopen(STATEFILE,"r");
#else
sp = fopen("state.data","r");
#endif
for(i=0;i<N;i++) fscanf(sp,"%lf\n",&Y[i]);
fclose(sp);}

void dump_(Y) 
double Y[N];
{FILE *fopen(),*sp;
int i;
if(Y[0]==Y[0]){ //sanity check
  #ifdef BIFURCATION //used for tests not in manuscript
  sp = fopen("end.data","w");
  #else
  sp = fopen("end.data","w");
  #endif
  //sp = fopen("end_b.data","w");
  for(i=0;i<N;i++) fprintf(sp,"%.16f\n",Y[i]);
  fclose(sp);
  }
}

int mas(n,amas,l)
        int *n;
        double *amas;
        int *l;
{return 0;}

int dummy(n,t,y,ydot)
        int *n;
        double *t;
        double *y;
        double *ydot;
{return 0;}

#endif

