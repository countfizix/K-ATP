
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

int main(){

double current[C];
double STEMP;
double CATEMP;
double ATEMP; 

int check;
double tempcheck;
double step;
int dir;

int i,j;
//int el;
//double time_;
double V0;
double iapp = IAPP;
//double iapp = -1.0e-4;
double start;
double anm,bnm;
double sinf,minf,hinf,htau,ntau,mtau,ninf,arg,dlinf,dltau,nminf,ainf;
double dminf, dhinf;
double alpm,alph,betm,beth;

double sinf2, diff;

double valtemp;


for(V0=-90;V0<30;V0+=0.1){
  //minf= boltz(V0,-28.0907+KOFF,9.7264); //sodium activation -28 9
  //hinf= boltz(V0,HHALF,-10.7665); //sodium inactivation -54
  //ninf= boltz(V0,-25.0,12.0); //fast potassium
  dlinf= boltz(V0,dloff,dlslope); //-45 7.5 //L type Calcium
  anm= 1.0/(1.0+(1.4/3.57)*exp(-0.062*V0)); //NMDA activation
  bnm= NMC+(1.0-NMC)/(1.0+(1.2/NMOFF)*exp(-V0/NMSLOPE)); //NMDA activation
  nminf = NMDATYPE*anm+(1-NMDATYPE)*bnm;

  minf= boltz(V0,-28.0907,8.7264); //sodium activation -28 9
  hinf= boltz(V0,-54.0,-8.7665); //sodium inactivation -54
  ninf= boltz(V0,-20.0,7.0); //fast potassium


  alph = 1.6e-4*exp(-V0/48.4);
  beth = boltz(V0,39.0,10.0);

  alpm = 1.967*efun((-1.0*V0+19.8)/10.0);
  betm = 0.046*exp(-V0/20.73);


  dminf = alpm/(alpm+betm);
  dhinf = alph/(alph+beth);

  //ainf= boltz(CATEMP,AHALF,ASLOPE); //worked at 120 2
  //sinf= boltz(ainf,S_HALF,S_SLOPE);

  current[I_NA1] = G_NA*pow(minf,3.0)*hinf*(V0 - E_NA);
  current[I_K1] = G_K*pow(ninf,3.0)*(V0 - E_K);
  current[I_L1] =  G_L*(V0 - E_L);
  current[I_LCa] =G_LCa*(V0- E_CA);
  current[I_CaL] =G_CaL*dlinf*(V0- E_CA);
  current[I_NMDA] = G_NMDA*nminf*(V0);

  current[I_NMDA_CA] = FNC*G_NMDA*nminf*(V0- E_CA);
  //current[I_CaN] = G_CaN*pow(dminf,2.0)*dhinf*(V0- E_CA);

  iapp = -2.0e-6*(V0+80.0);


//solve for K-ATP activation
  STEMP = iapp  - current[I_NA1] - current[I_K1] - current[I_L1] - current[I_LCa] -current[I_CaL] -current[I_NMDA];
  valtemp = STEMP +current[I_NMDA]+current[I_LCa]+ current[I_L1];

  sinf = pow(STEMP/(G_KATP*(V0-E_K)),1.0/spower); 

  fflush(stdout);

  if(!(sinf > 1.0 || sinf < 0.0)){
    //ATEMP = 1.0 - pow(1.0-sinf,S_SLOPE); //fadp
    //ATEMP = pow(ATEMP,-1.0/fpower);  //1 + SHALF/A
    //ATEMP = S_HALF/(ATEMP-1.0);//ADP 
    //ATEMP = pow(sinf,0.5);
    ATEMP =  sinf*S_HALF/(1.0-sinf);
    //ATEMP = ATEMP - NN414;

    sinf2 = ATEMP/(S_HALF+ATEMP+2.5);
    diff = G_KATP*(pow(sinf,spower)-pow(sinf2,spower))*(E_K-V0); // additional current from NN414 (steady state)
    CATEMP = ATEMP*KNST/BEE/(ATEMP+BASE); //Solve for I_CAP
    CATEMP = CHALF/pow(1.0/CATEMP-1.0,1.0/CSLOPE); //solve for Y[CA]
    //CATEMP = AHALF - ASLOPE*log(1.0/(S_HALF-S_SLOPE*log(1.0/sinf-1.0))-1.0);
    if(CATEMP==CATEMP&&CATEMP>0&&CATEMP<300.0){
      printf("%e  %e  %e  %e  %e\n",CATEMP, V0, ATEMP,sinf,diff+current[I_NMDA]);
      fflush(stdout);
    }

  }//endif*/

  /*CATEMP = current[I_LCa] + current[I_CaL] + current[I_NMDA_CA];
  ATEMP = -CATEMP*BEE/(KNST*I_CAPMAX);
  if(ATEMP < 1){
     ATEMP = BASE*ATEMP/(1.0-ATEMP);
     sinf = ATEMP/(ATEMP+S_HALF);
     sinf = pow(sinf,spower);
  }
  else{
     ATEMP = -1.0;
     sinf = 1.0;
  }
  valtemp = G_KATP*pow(sinf,spower)*(V0-E_K) +current[I_NMDA]+current[I_LCa]+ current[I_L1];
  CATEMP = -CATEMP/I_CAPMAX;
  CATEMP = CHALF/pow(1.0/CATEMP-1.0,1.0/CSLOPE);
  if(CATEMP==CATEMP&&CATEMP>0){ 
    printf("%e  %e  %e  %e",CATEMP, V0, ATEMP, sinf, -valtemp);
    //for(j=0;j<C;j++){
    //  printf("%e  ", -1e6*current[j]);
    //}
    printf("\n");
  }
  /*else{
      printf("%e  %e  %e  %e\n",V0, -1.0, -1.0, -1.0); 
  }*/
}


//fclose(fp);



return 0;
}
