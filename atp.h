#ifndef ATP_H_ //only defines once - important as we can use undefine on a parameter to vary it.
#define ATP_H_

#define THRESHOLD -14.0 /* -14.0 */
#define INCREMENT 100 
#define SAMPLE 0


#define V1 0  /* mV */
#define H1 1  
#define N1 2 
#define M1 3  
#define Ca0 4  /* nM */
#define CB0 5  /* nM */
#define A 6  /* uM? */
#define DL 7  


#define N  8/* number of state variables*/

#define PI 3.14159

#define EPSILON 0.0

#define I_NA1 0  
#define I_K1 1 
#define I_L1 2
#define I_LCa 3 //calcium leak 
#define I_CaL 4 //l-type calcium
#define I_NMDA 5 
#define I_NMDA_CA 6 
#define I_KATP 7 
#define I_CAP 8

#define C 9 /* number of currents*/

#define CM 1.00 /*uF/cm2*/
#define PHI 5.0

#define I_APP 0.0e-6  /*uA/cm2*/

#define E_NA  60.0
#define E_CA  60.0 //based on voltage curve fit

#define E_HCN -40.0

#define E_K  -90.0
#define E_L  -53.0 /* -60 mV*/
#define G_NA 1500.0e-6 //0 for TTX
#define G_K   250.0e-6 
#define G_L   5.0e-6  

#define G_HCN 0.0e-6

#ifdef BIFURCATION
//extern double params[3];

#define PMIN_0 2.5
#define PMAX_0 5.0
#define PSTEP_0 2.5

#define PMIN_1 0.0e-5
#define PMAX_1 5.0e-5
#define PSTEP_1 1.0e-6

#define PMIN_2 0.0
#define PMAX_2 0.0
#define PSTEP_2 1e-6
#endif


#define G_LCa  0.1e-6 //old 0.2 
#define G_CaL  1.0e-6  //new 1.8 old  1.5!!! //2.0


#define CHALF 500.0
#define CSLOPE 1.0


#define G_KATP  40.0e-6 //for all

#define BASE 1.0 //

#define TBUFF 3.0e3 //2.0e3 for Calbindin negative

#define KB 1.0e-5 //1e-5 -its largely insensitive to this - acts as a (biased) low pass filter
#define KU KB*(TBUFF)  //unbuffering 1/ms
  //buffering 1/(nM ms)

#define DCA 0.006


#define STIM_FILE "stim.dat"
#define SCALE 0.0e-4


#define RATIO 10.0
#define FARADAY  96485.0 /*old value 96520 was work 96485 need to redo with this value*/
#define D_S  10.0 //big compartment means reduced calcium influx 7.5
#define L_S  25.0
   
#define ALPHA 6.25 /* ms */
#define TAUSYN 1.00 /*  1.0 3.0 2.0   */
#define SS1 0.0
#define SS2 0.0
#define SS1X 0.0
#define SS2X 0.0

#define CSLOPE 1.0

#define fsca 1.0 // buffering is dynamic rather than instantaneous in this model
#define START_TIME  0
#define FNC 0.00 // no NMDA calcium
//varied parameters

#define NMC 0.02

#ifndef BIFURCATION
#define G_NMDA 2.500000e-05 //0 for control
#define S_HALF 5.000000e+00 //2.5 for NN414
#endif
#define IAPP 0.0

#define S_SLOPE 0.5
#define spower 4.0

#define I_CAPMAX 0.00016 // 0.001 !!!!



#define BEE 1.5 //ADP accumulation 0.25


//model is bistable
//#define BEE  2.0e-2
//#define KNST 1.0e-3
#define KNST 0.25*BEE //0.2*BEE for BK acid

#define HHALF -54.0289 //-54

#define NMSLOPE 9.0 //9
#define ELECT 0
#define NMOFF 102.0 //102 larger moves left 250
#define dloff -35.0 //-45
#define dlslope 7.5 //7.5

#define NMDATYPE 0 // 1 = original 0 = sharp

#define KOFF 0.0
#define SOFF 0.0

#define TIMING_FILE "derp.dat"
#define PRINT_STATES 1
#define ENDTIME 10000
#define WARMUP 40000

//############################

#define PRINT  0
#define DEBUG 0


/* INTEGRATION PARAMETERS */

/* these only apply to C code */
#define LWORK N*(4*N+8)+7 /* N(4*N + 8)+7   679 for NN = 3   327 for NN=2*/ 
#define LIWORK 3*N+7 /*3*N+7 43 for NN = 3  31 for NN= 2*/
#define LRCONT 4*N+4 /*4*N+4  52 for NN = 3  36 for NN=2*/
#define DURATION 5
#define TIMEON 0
#define TIMEOFF 0
#define PERIOD  250

#endif //end ifndef
//endall
