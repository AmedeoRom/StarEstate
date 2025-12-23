/* software by chris belczynski: last modified on Nov21, 2022 */
/* single star evolution, after Jarrod Hurley 2000 */
/* MODIFIED VERSION FOR MY BINARY.C PROGRAM */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sinbin.h"


                                              /* for AGB stars: */
/* Mcbur1 defined in sinbin.h */              /* He core mass at BAGB: for M_he(t=t_bagb) < Mcbur1 -> */ 
/* 1.6 Msun Hurley value (Mzams=6.3) */       /* formation of degenerate CO core (ONe core never forms). */ 
/* 1.83 standard used here (Mzams=7) */       /* So for M_he(t=t_bagb) < Mcbur1 we have either: */
                                              /* 1) if M_co reaches M_Ch ->SN IIa (K=15): CO ignition under */
                                              /*    degenerate conditions -- total disruption, no remnanat */ 
                                              /* 2) if CO core gets exposed with M_co<M_Ch; CO WD formation */ 
                                              /*    a) wind mass loss */ 
                                              /*    b) RLOF mass loss */
                                              /*    c) CO grows all the way (to surface) in He star (K=8,9) */ 


double Mcocrit=1.08;                          /* Mcbur1 <= M_he(t=t_bagb) <= Mcbur2: formation of */
double Mecs=1.38;                             /* semi-degenerate CO core, if Mco<Mcocrit and core exposed */
double Mnsecs=1.26;                           /* then CO WD formation, otherwise: */ 
                                              /* C ignites off-center when M_co reaches Mcocrit=1.08 Msun */ 
                                              /* that forms degenerate ONe core (1.08 Msun), then: */  
                                              /* 1) ONe core reches Mecs=1.38 Msun and explodes in ECS SNa */ 
                                              /*    forming NS of Mnsecs=1.26Msun: Mbaryon-Mgrav=0.075*Mgrav^2 */
                                              /*    also possibly lowered kick */
                                              /* 2) ONe core gets exposed with mass<Mecs; ONe WD formation */
                                              /*    a) wind mass loss */
                                              /*    b) RLOF mass loss */
                                              /*    c) ONe grows all the way (to surface) in He star (K=8,9) */
                                              /* small cheat here: i treat ONe core as CO core, and i use */ 
                                              /* its mass and M_he(t=t_bagb) to know what type of core it is */ 


/* Mcbur2=2.25; (standard) */                 /* He core mass at BAGB: for M_he(t=t_bagb) > Mcbur2 -> */
/* Mcbur2 defined in sinbin.h */              /* formation of nondegenerate CO core and then heavier elements */
                                              /* ignite nondegenerately; leading to formation of FeNi core */
                                              /* and then SN II and NS/BH formation */


double Mgrow=1.05;                            /* Hurley2000 eq.(75): Mc,SN=max(M_Ch,0.773*M_he(t=t_bagb)-0.35) */
                                              /* is in his code (and StarTrack as well) different: */ 
                                              /* Mc,SN=max(M_Ch,0.773*M_he(t=t_bagb)-0.35,Mgrow*M_co(t=t_bagb)) */
                                              /* the third element: Mgrow*M_co(t=t_bagb) works only for AGB stars */
                                              /* otherwise (e.g. He stars) it is 0 */    
                                              /* M_co(t=t_bagb) -- starting CO core mass at Base of AGB */

                                          
                                              /* for He stars */ 
                                              /* the same procedure is follwed as above for AGB stars, */ 
                                              /* but i adopt M_he(t=t_bagb)=initial mass of He star: M_he,i  */
                                              /* see last paragraf of 6.1 Hurley et al. 2000, rest deduced */ 
                                              /* from his paper and code: */   
                                              /* 1) M_he,i=mass of helium star at the end of helium MS if */ 
                                              /*    He star is evolved through helium MS */
                                              /* 2) M_he,i=mass of helium core of exposed AGB star, if He star */
                                              /*    is formed out of AGB star */ 
                                             

double Mzams;                                 /* starting ZAMS mass [Msun] -- always constant */
double M0;                                    /* effective initial mass [Msun] -- changed in evolution */

int mflag=1;                                  /* 1 -- activates mass loss, 0 -- no mass loss */     
double ETA=0.5;                               /* Reimers coefficient to increase or decrease mass loss */

/* TB,TVIR,TE:  control times [Myr] for calculations: last step was for TB at TVIR and now return values for TE */
/* MC:  He core mass from last step only for HGf() in single() and used throughtout binary.c */ 
/* FLAG:  0-M0 changes on HG, 1-M0 stoped changing on HG */ 
/* MHER, MCOR -- always keeps mass of Helium core and CO core -- needed for mass transfer in binary.c */
/* DT -- sets the next time step taken from single() and passed to binary.c */
/* MPRE -- keeps the mass of star just before SN explosion */
/* KP -- if type K was changed before Remnantf() it keeps the previous one */
/* TSTART: passes just value of tstart for HSMSf and HSGBf from binary.c to single() during mass transfer */
/* and when star changes from K=7 to K=8 */
/* FRAC -- if there was supernova it gives info on amount of fall back [0-1 of fall back] */
/* MOUT1 -- [Msun] sub Chandrasekhar SN Ia outburst mass for H accreting He WD */
/* MOUT2 -- [Msun] transformation of He WD to He MS mass (He accretion at high rate) */
/* ECSSN -- 1: marks ECS SN of regular star (AGB or He giant) */
/*          2: marks ECS/AIC of ONe WD -> NS */ 
/*          0: all other cases (regular SN and no SN) */ 
/* SNR: for SN type; 0-no SN, 1-SN took place */
/* FMS: fraction of star mass accreted during RLOF on MS: for compactf() functions */


double TB,TVIR,TE,DT,MPRE;
double LRR,RRR,MCRR,FRAC;
double DMMT;                                                                       /* not changed here!!! */
double MOUT1,MOUT2;                                                                /* not changed here!!! */
int DERIV,SN;                                                                      /* not changed here!!! */ 
int FLAG,ECSSN,FB;
double Fms;                                                                        /* not changed here!!! */

void singl(double *MzamsR, double *M0R, double *MR, int *KR, double *TBR, double *TVIRR, double *TER, 
            double *LR, double *RR, double *MCR, double *MHER, double *MCOR, int *FLAGR, double *DTR, 
            double *MPRER, int *KPR, double *TSTART, double *FRACR, double DMMTR, int DERIVR, 
            double MOUTR1, double MOUTR2, int *ECSSNR, int *FBR, int SNR, double FMS)
{ /* function to be called by my BINARY program */
  /* M_hookf(),M_HeFf(),M_FGBf() must be caaled before single() in binary.c */
  /* coeff_aa() and coeff_bb() must be caaled before single() in binary.c */
  
 double t,tstart;
 double M,dM,Mc,Mhe,Mco;
 double L,R;
 int K,Kp;
 int stop=0;                                                 /* control when to stop calculations */
 int mark=0;                                                 /* to mark if it is first calculation step or not */

 Mzams=(*MzamsR);                                              /* ZAMS mass */
 M0=(*M0R);                                                    /* effective initial mass from last step */ 
 M=(*MR);                                                      /* current mass from last step */
 K=(*KR);                                                      /* type from last step */
 TB=(*TBR);
 TVIR=(*TVIRR);
 TE=(*TER);
 MCRR=(*MCR);
 Mhe=(*MHER);
 Mco=(*MCOR);
 FLAG=(*FLAGR);
 DT=(*DTR);
 MPRE=(*MPRER);
 Kp=(*KPR);
 tstart=(*TSTART);
 FRAC=(*FRACR);                                   
 DMMT=1.0e-06*DMMTR;                /* [Msun/Myr] --> [Msun/yr], positive: mass loss, negative: mass gain */
 DERIV=DERIVR;
 MOUT1=MOUTR1;
 MOUT2=MOUTR2;
 ECSSN=(*ECSSNR);
 FB=(*FBR);
 SN=SNR;
 Fms=FMS;                           /* set as global variable; so it is visible in compactf() functions */ 
 
 if(fabs(TE)<acc) {
   t=0.0;
   mark=1;
   ZAMSf(&M,&t,&Mhe,&Mco,&stop);                            /* first call -- initialize the star and stop */
 }  

 if((K==0 || K==1) && stop==0) {
   t=TB;
   mark=1;
   K=MSf(&M,&t,&Mhe,&Mco,&stop);                            /* change of M0 if wind is present */  
 }
 
 if(K==2 && stop==0) {
   if(mark==0) t=TB;
   mark=1;
   K=HGf(&M,&t,&tstart,&Mhe,&Mco,&Kp,&stop);             /* change of M0 if wind is present */
 }

 if(K==3 && stop==0) {                                       /* change of M0 if wind is present */
   if(mark==0) t=TB;
   mark=1;
   K=RGf(&M,&t,&tstart,&Mhe,&Mco,&Kp,&stop); 
 }

 if(K==4 && stop==0) {                                       /* change of M0 if wind is present and if M0<=M_HeF */
   if(mark==0) t=TB;
   mark=1;
   K=HEf(&M,&t,&tstart,&Mhe,&Mco,&stop);             
 }

 if((K==5 || K==6) && stop==0) {                             /* evolves also K==6 */
   if(mark==0) t=TB;
   mark=1;
   K=AGBf(&M,&t,&tstart,&Mhe,&Mco,&Kp,&stop,K);            
 }

 if(K==7 && stop==0) {                                       /* change of M0: by formation of HS and/or by wind */
   if(mark==0) t=TB;
   mark=1;  
   K=HSMSf(&M,&t,&tstart,&Mhe,&Mco,&stop);                              
 }

 if((K==8 || K==9) && stop==0) {                             /* change of M0 if HSMSf() skipped in evolution */                             
   if(mark==0) t=TB;
   mark=1;
   K=HSGBf(&M,&t,tstart,&Mhe,&Mco,&Kp,&stop);
 }

 if(K==-1 && stop==0) {                                      /* change of M0: by formation of remnant */
   t=TE;
   K=REMNANTf(&M,&t,&Mhe,&Mco,Kp);
   stop=1;
 }

 if((K>=10 && K<=17) && stop==0) {
   L=1.0e-50;                                 /* so small, as to avoid getting significant Teff for NS/BH */
   Mc=0.0;                                    /* with get_T() func., for WDs use Lwdf() to calculate L */
   t=TE;
   if(K==10 || K==11 || K==12 || K==16 || K==17) {
     dM=DMMT*1.0e+6*(TE-TB);                  /* DMMT -- MT rate in [Msun/yr]: [Msun/Myr] * [Myr] = [Msun] */
     M-=dM;                                   /* possible both accretion and mass loss */
     M0=M; 
     if(M>MOUT2 && (K==10 || K==17)) {        /* He shell flash, He WD -> He ZAMS transformation */
       K=7;
       M0=M;
       TVIR=tstart=0.0;
       DT=delms*thsmsf(M);
     }
     else if((M>MCh || M>MOUT1) && (K==10 || K==11 || K==17)) {    /* SN Ia */
       K=15;
       R=0.0;
       stop=1;
     }  
     else if((M>Mecs || M>MOUT1) && K==12) {    /* M>Mecs: AIC: ONeMg WD->NS, M>MOUT1: SN Ia */
       if(M>Mecs && ECS==1) {   /* AIC, now takes place over Mecs and not over MCh */
         K=13;
         M=Mnsecs;
         R=Rnsf();
         FRAC=0.0;                        /* writes how much fallback [0:1] took place into global var. */
         FB=0;
         ECSSN=2; 
       }
       else if(M>MOUT1) {                       /* SN Ia */
         K=15;
         R=0.0;
         stop=1;
       }  
       else if(M>MCh) {                    /* SN Ia -- can get here only if ECS not allowed */
         K=15;                             /* if M:  Mecs<M<MCh do nothing: accumulate */
         R=0.0;
         stop=1;
       }                  
       else                                /* if M:  Mecs<M<MCh do nothing: accumulate */
         R=Rwdf(M,K);                      /* can get here only if ECS/AIC not allowed */
     } 
     else
       R=Rwdf(M,K);
   }  
   else if(K==13) {
     dM=DMMT*1.0e+6*(TE-TB);                  /* only accretion possible */
     M-=dM;
     M0=M;
     if(M>Mmaxns) {                           /* AIC: NS->BH, no SN assumed */
       K=14;
       R=Rbhf(M);
     }
     else  
       R=Rnsf();
   }
   else if(K==14) {
     dM=DMMT*1.0e+6*(TE-TB);                  /* only accretion possible */  
     M-=dM;
     M0=M;
     R=Rbhf(M);
   }
   else
     R=acc;    
   (LRR)=L; (RRR)=R; MCRR=Mc; TVIR=0.0; MPRE=0.0; 
   DT=1.0e+50;                                           /* big number to avoid its usage in binary.c */     
 } 


 (*MzamsR)=Mzams;
 (*M0R)=M0;                                                    
 (*MR)=M;                                                     
 (*KR)=K;                                                     
 (*TBR)=TB;
 (*TVIRR)=TVIR;
 (*TER)=TE;
 (*LR)=LRR;
 (*RR)=RRR;
 (*MCR)=MCRR;
 (*MHER)=Mhe;
 (*MCOR)=Mco;
 (*FLAGR)=FLAG;
 (*DTR)=DT;
 (*MPRER)=MPRE; 
 (*KPR)=Kp;
 (*TSTART)=tstart;
 (*FRACR)=FRAC;
 (*ECSSNR)=ECSSN;
 (*FBR)=FB;
}  



/* -------------------------------------- EVOLUTIONARY FUNCTIONS -------------------------------------- */

void ZAMSf(double *M, double *t, double *Mhe, double *Mco, int *stop)
{ /* writes ZAMS values for star of M at time t which must be 0.0 */
  /* pointers to M,t passed only for compatibility with other evolutionary functions */ 
 double Mche,Mcco,dM,dt,L,R; 
 int K;
 
 if((*M)<0.7) K=0;
 else K=1;
 Mche=Mcco=0.0;
 dM=0.0;
 dt=0.0;
 L=Lzamsf(*M);
 R=Rzamsf(*M);
 (LRR)=L; (RRR)=R; MCRR=Mche; TVIR=0.0;  DT=0.0; MPRE=0.0;
 (*Mhe)=Mche;
 (*Mco)=Mcco;
 (*stop)=1;
}


int MSf(double *M, double *t, double *Mhe, double *Mco, int *stop)
{ /* evolve star of M through Main Sequance and writes the data, supplied *t shall be 0.0 */
  /* first step is for MS and not for ZAMS (so first files output is for t!=0.0) */
  /* 1-means value from previous step, 2 from current step; p-star of deminished mass */
  /* there are descrepancies in L,R, and t with Jarrod code, I think them numerical noise */   
  /* t2p is virtual time here */
 double trest1,trest2,dt,t2p,dtms2p,tms2p,t1,dtms1,t2,tms2;
 double Mt,Mche,Mcco,dM,Mold,rat;
 double T,L,R,L2p,L1,R2p,R1;
 int K;
 int i=0,mark2=0;           /* to mark when to stop decreasing time step because of jump in R */
 int mm=0;                  /* make maximum of 5 time step decreaments because of jump in R */

 Mt=(*M);
 if(Mt<0.7) K=0;
 else K=1;
 Mche=Mcco=0.0;
 
 t2p=TVIR;
 L2p=Lmsf(Mt,t2p);
 R2p=Rmsf(Mt,t2p);
 tms2p=tmsf(Mt);
 if(DERIV==1)                               /* derivative estimation run */
   dtms2p=TE-(*t);
 else
   dtms2p=min(delms*tms2p,TE-(*t));         /* choose FIRST time step as miniumum from evolution of this star */
                                            /* or as given from binary.c taken from other consideration */
 while(t2p<tms2p || fabs(tms2p-t2p)<acc) {
   t1=t2p; 
   L1=L2p;
   R1=R2p;
   dtms1=dtms2p;
   if(WIND==1)
     dM=windf1(Mt,L1,R1,Mche,K)*1.0e+6*dtms1;
   else if(WIND==2)
     dM=windf2(Mt,L1,R1,Mche,K)*1.0e+6*dtms1;
   else if(WIND==3)
     dM=windf3(Mt,L1,R1,Mche,K)*1.0e+6*dtms1;   
   
   if(dM>(maxdM*Mt)) {                       /* check that mass loss do not go over 1% of star mass */
     dtms2p=0.01*dtms2p*Mt/dM;
     dtms1=dtms2p;
     dM=0.01*Mt;
     if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;    /* if DMMT=0, no MT and do what requested. during MT dt forced */
   }                                                  /* during derivative calls do not stop here */                                              

   t2=t1+dtms1;
   tms2=tms2p;
   Mold=Mt;
   Mt-=dM;
   M0=Mt;
   
   tms2p=tmsf(Mt);
   if(dM>acc)                /* mass loss from the star: rat not applicable: fresh fuel not taken away from the core */
     rat=1.0;
   else {                    /* mass gain: rat may apply, as fresh fuel added to the core */  
     if(Mt>0.35 && Mt<1.25) rat=1.0;         /* star has radiative core: Hurley et al 2002, Sec.2.6.6.1) */ 
     else rat=Mold/Mt;                       /* convective core: Tout et al. 1997, MNRAS 291, 732 Sec. 5.1, eq. 41 */
   }
   t2p=rat*(tms2p/tms2)*t2;

   L2p=Lmsf(Mt,t2p);
   R2p=Rmsf(Mt,t2p);
  
   if(fabs(R2p-R1)>(maxdR*R1) && mark2==0) {              /* check that radius change no more then 10% */ 
     Mt+=dM;
     M0=Mt;
     t2p=t1;                                /* undo above part of loop */
     tms2p=tmsf(Mt);                        /* t1,L1,R1,tms2 unchanged */
     L2p=L1;                                /* dtms1,dM,t2 changed via others changed param. */ 
     R2p=R1;
     dtms2p*=0.5;                           /* decrease by half time step */
     i++;
     if(i>5)                                /* make maximum 5 time step decreaments */
       mark2=1;
     if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
     continue;                              /* go directly from here back to begining of while() loop */
   }
   i=0;
   mark2=0;
                                    
   L=L2p;
   R=R2p;
   T=get_T(L,R);
   (*t)+=dtms2p;
   dt=dtms2p;
   (*M)=Mt;

   tms2p=tmsf(Mt);
   dtms2p=delms*tms2p;
   trest1=max(tms2p-t2p,0.1*acc);
   dtms2p=min(dtms2p,trest1);

   if(fabs(TE-(*t))<acc || mm==1) {                               /* to stop calculation at TE */
     if(mm==1) TE=(*t);
     (LRR)=L; (RRR)=R; MCRR=Mche; TVIR=t2p; MPRE=0.0; DT=dtms2p; 
     (*stop)=1;
     (*Mhe)=Mche;
     (*Mco)=Mcco;
     break;                               
   }   
   trest2=max(TE-(*t),0.1*acc);
   dtms2p=min(dtms2p,trest2);
   if(fabs(tms2p-t2p)<acc) 
     break;                                               /* to stop the loop after last step */ 
 }                      
 
 if((*stop)!=1 || fabs(tms2p-t2p)<acc) {
   K=2;                                       /* HG star */      
   TVIR=0.0;
   DT=delhg*thgf(Mt);                         /* sets step for next evolutionary phase to be supplied to binary.c */ 
 }  
 
 return K; 
} 


int HGf(double *M, double *t, double *tstart, double *Mhe, double *Mco, int *Kp, int *stop)
{ /* evolve star of M through Hartzprung Gap and writes the data, supplied current time *t from end of MS */
  /* returns type after HG phase K: 3-red giant, 4-HB star, 7-helium MS star, 10-helium WD */   
  /* first output is already somewhat (in time) into HG */
  /* t2p is virtual time here */
 double trest1,trest2,dt,t2p,dthg2p,thg2p,t1,dthg1,t2,thg2,tend2p;
 double Mt,Mche1,Mche2p,Mcco,dM;
 double T,L,R,L2p,L1,R2p,R1;
 double Mcbgb;
 int K;
 int mark=0;                                /* to mark the moment when M0 stops changing */
 int mark1=0;                               /* to mark when core mass exceeds star mass */
 int i=0,mark2=0;
 int mm=0;
 
 Mt=(*M); 
 K=2;
 Mcco=0.0;

 if(fabs(TVIR)<acc) {
   t2p=tmsf(Mt);
   thg2p=thgf(Mt);
   Mche2p=Mchgf(Mt,t2p);
   L2p=Lhgf(Mt,t2p);
   R2p=Rhgf(Mt,M0,t2p);
   perturb(Mt,Mche2p,&L2p,&R2p,t2p,K);
   FLAG=0;
 }  
 else { 
   t2p=TVIR;
   thg2p=thgf(Mt);
   Mche2p=MCRR;
   L2p=Lhgf(Mt,t2p);
   R2p=Rhgf(Mt,M0,t2p);
   perturb(Mt,Mche2p,&L2p,&R2p,t2p,K);
   if(FLAG==1) {
     t1=t2p;
     thg2p=thgf(M0);
     dthg1=min(delhg*thgf(M0),TE-(*t));
     L2p=Lhgf(M0,t2p);
     R2p=Rhgf(Mt,M0,t2p);
     perturb(Mt,Mche2p,&L2p,&R2p,t2p,K); 
     mark=1;   
   }
 }  
 
 tend2p=tmsf(Mt)+thgf(Mt);
 if(DERIV==1)                                /* derivative estimation run */
   dthg2p=TE-(*t);  
 else    
   dthg2p=min(delhg*thgf(Mt),TE-(*t));
 while((t2p<tend2p || fabs(tend2p-t2p)<acc) && mark==0) {
   t1=t2p; 
   L1=L2p;
   R1=R2p;
   Mche1=Mche2p;
   dthg1=dthg2p;
   if(WIND==1)
     dM=dthg1*1.0e+6*windf1(Mt,L1,R1,Mche1,K);
   else if(WIND==2)
     dM=dthg1*1.0e+6*windf2(Mt,L1,R1,Mche1,K);
   else if(WIND==3)
     dM=dthg1*1.0e+6*windf3(Mt,L1,R1,Mche1,K);
   
   if(dM>(maxdM*Mt)) {                       /* check that mass loss do not go over 1% of star mass */
     dthg2p=0.01*dthg2p*Mt/dM;
     dthg1=dthg2p;
     dM=0.01*Mt;
     if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
   }
 
   t2=t1+dthg1;
   thg2=thg2p;
   Mt-=dM;
   M0=Mt;
   thg2p=thgf(Mt);
   tend2p=tmsf(Mt)+thgf(Mt);
   t2p=tmsf(Mt)+(t2-tmsf(Mt+dM))*thg2p/thg2;

   Mche2p=Mchgf(Mt,t2p);    
   Mche2p=max(Mche1,Mche2p);                /* stop of core mass grow: see end of sec. 5.1.2 of paper0 */ 
   if(Mche2p>Mt) {                                 /* check that He core mass will not exceed star mass */
     Mt+=dM;                                       /* undo above part of loop */
     M0=Mt;
     t2p=t1;
     Mche2p=Mt;                                    /* approximation: I set that Mc=Mt from previous and last step */
     mark1=1;
     break;
   }

   Mcbgb=Mcbgbf1a(Mt);
   if(Mcbgb<Mche2p) 
     if(mark==0) {
       Mt+=dM;
       M0=Mt;
       mark=1;
       FLAG=1;
       break;
     }  

   L2p=Lhgf(Mt,t2p);
   R2p=Rhgf(Mt,M0,t2p);
   perturb(Mt,Mche2p,&L2p,&R2p,t2p,K);

   if(fabs(R2p-R1)>(maxdR*R1) && mark2==0) {              /* check that radius change no more then 10% */
     Mt+=dM;
     M0=Mt;
     t2p=t1;                                /* undo above part of loop */
     tend2p=tmsf(Mt)+thgf(Mt);
     thg2p=thgf(Mt);                        /* Mche1,t1,L1,R1,thg2 unchanged */
     L2p=L1;                                /* dthg1,dM,t2 changed via others changed param. */
     R2p=R1;
     Mche2p=Mche1;
     dthg2p*=0.5;                           /* decrease by half time step */
     i++;
     if(i>5)                                /* make maximum 5 time step decreaments */
       mark2=1;
     if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
     continue;                              /* go directly from here back to begining of while() loop */
   }
   i=0;
   mark2=0;

   L=L2p;
   R=R2p;
   T=get_T(L,R);
   dt=dthg2p;
   (*t)+=dthg2p;
   dt=dthg2p;
   (*M)=Mt;
   
   thg2p=thgf(Mt);
   dthg2p=delhg*thg2p;
   trest1=max(tend2p-t2p,0.1*acc);
   dthg2p=min(dthg2p,trest1);
   
   if(fabs(TE-(*t))<acc || mm==1) {                                         /* to stop calculation at TE */
     if(mm==1) TE=(*t);
     (LRR)=L; (RRR)=R; MCRR=Mche2p; TVIR=t2p; MPRE=0.0; DT=dthg2p;
     (*stop)=1;
     (*Mhe)=Mche2p;
     (*Mco)=0.0;
     break;                               
   } 
   trest2=max(TE-(*t),0.1*acc);
   dthg2p=min(dthg2p,trest2);
   if(fabs(tend2p-t2p)<acc) 
     break;                                 /* to stop the loop after last step */ 
 }                                            

 if(mark==1) {                                /* M0 stopped decreasing */
   t2p=t1;
   tend2p=tmsf(M0)+thgf(M0);
   dthg2p=dthg1;
   mark2=i=0;
   while(t2p<tend2p || fabs(tend2p-t2p)<acc) {
     L1=L2p;
     R1=R2p;
     Mche1=Mche2p;    
     if(WIND==1)
       dM=dthg2p*1.0e+6*windf1(Mt,L1,R1,Mche1,K);
     else if(WIND==2)
       dM=dthg2p*1.0e+6*windf2(Mt,L1,R1,Mche1,K);
     else if(WIND==3)
       dM=dthg2p*1.0e+6*windf3(Mt,L1,R1,Mche1,K);
    
    
     if(dM>(maxdM*Mt)) {                       /* check that mass loss do not go over 1% of star mass */
       dthg2p=0.01*dthg2p*Mt/dM;
       dM=0.01*Mt;
       if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
     }
    
     Mt-=dM;
     t2p+=dthg2p;
     
     Mche2p=Mchgf(M0,t2p);
     Mche2p=max(Mche1,Mche2p);                /* stop of core mass grow: see end of sec. 5.1.2 of paper0 */
     if(Mche2p>Mt) {                                 /* check that He core mass will not exceed star mass */
       Mt+=dM;                                       /* undo above part of loop */
       t2p-=dthg2p;
       Mche2p=Mt;                                    /* approximation: I set that Mc=Mt from previous and last step */
       mark1=1;
       break;
      }
           
     L2p=Lhgf(M0,t2p);
     R2p=Rhgf(Mt,M0,t2p);
     perturb(Mt,Mche2p,&L2p,&R2p,t2p,K);

     if(fabs(R2p-R1)>(maxdR*R1) && mark2==0) {              /* check that radius change no more then 10% */
       L2p=L1;
       R2p=R1;
       Mche2p=Mche1;
       Mt+=dM;
       t2p-=dthg2p;
       dthg2p*=0.5;                           /* decrease by half time step */
       i++;
       if(i>5)                                /* make maximum 5 time step decreaments */
         mark2=1;
       if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
       continue;                              /* go directly from here back to begining of while() loop */
     }
     i=0;
     mark2=0;

     L=L2p;
     R=R2p;
     T=get_T(L,R);
     (*t)+=dthg2p;         
     dt=dthg2p;
     (*M)=Mt;
     
     thg2p=thgf(Mt);
     dthg2p=delhg*thg2p;
     trest1=max(tend2p-t2p,0.1*acc);
     dthg2p=min(dthg2p,trest1);
     
     if(fabs(TE-(*t))<acc || mm==1) {                                         /* to stop calculation at TE */
       if(mm==1) TE=(*t);
       (LRR)=L; (RRR)=R; MCRR=Mche2p; TVIR=t2p; MPRE=0.0; DT=dthg2p; 
       (*stop)=1;
       (*Mhe)=Mche2p;
       (*Mco)=0.0;
       break;                               
     } 
     trest2=max(TE-(*t),0.1*acc);
     dthg2p=min(dthg2p,trest2);
     if(fabs(tend2p-t2p)<acc) 
       break;                                 /* to stop the loop after last step */     
   }
 }
                                             /* mark1==0 means that core mass is smaller then star mass */ 
                                             /* mark1==1 core mass has reached star mass */             
 (*tstart)=0.0;                              /* if Helium star is formed, it is on its ZAMS */
 
 if((*stop)!=1 || fabs(tend2p-t2p)<acc) {
   if(mark1==0 && M0<=M_FGB) {                 /* red giant */
     K=3;
     TVIR=0.0;
     if(tbgbf(M0)<=txf(M0))                                           
       DT=delrg*(tinf1f(M0)-tbgbf(M0));
     else
       DT=delrg*(tinf2f(M0)-tbgbf(M0));  
     DT=min(DT,the1f(M0)-tbgbf(M0));
   }  
   else if(mark1==0 && M0>M_FGB) {             /* CHeB star */ 
     K=4;
     TVIR=0.0;
     DT=delhe*thef(M0);
   }
   else if(mark1==1 && M0>M_HeF) {               /* ZAMS Helium star */
     K=7; 
     TVIR=0.0;
   }  
   else if(mark1==1 && M0<=M_HeF) {            /* formation of remnant */
     (*Mhe)=Mche2p;
     (*Mco)=0.0;
     (*Kp)=K;
     K=-1;
     TVIR=0.0;
   }
   else 
     fprintf(fp0,"error: in HGf() star type K was not set\n");        
 }

 return K;
}


int RGf(double *M, double *t, double *tstart, double *Mhe, double *Mco, int *Kp, int *stop)
{ /* evolve star of M through Red Giant Branch and writes the data, supplied current time *t */
  /* returns type after RG phase K: 4-HB star, 7-helium MS star, 10-helium WD */
  /* first output is already somewhat (in time) into RG */
 double Mt,Mche1,Mche2,Mcd,Mcco,dM;
 double tbgb,the1,trg,tx,tinf1,tinf2,dtrg,trest1,trest2,dt,tvir;
 double L,R,T,L1,L2,R1,R2;
 int K;
 int mark1=0;                                             /* to mark when core mass exceeds star mass */
 int i=0,mark2=0;
 int mm=0;

 Mt=(*M);
 K=3;
 Mcco=0.0;

 tbgb=tbgbf(M0);
 the1=the1f(M0);
 trg=the1-tbgb;
 if(fabs(TVIR)<acc) 
   tvir=tbgb;
 else  
   tvir=TVIR;

 if(trg<acc) {                                                   /* ensures that star has RG branch */
   fprintf(fp0,"error: I should not be in RGf()\n");
   K=4;
   return K; 
 }

 tx=txf(M0);
 tinf1=tinf1f(M0);
 tinf2=tinf2f(M0);

 if(DERIV==1)                                /* derivative estimation run */
   dtrg=TE-(*t); 
 else if(tvir<=tx)                                           /* set the first time step */
   dtrg=min(delrg*(tinf1-tvir),TE-(*t));
 else
   dtrg=min(delrg*(tinf2-tvir),TE-(*t));  
 trest1=max(the1-tvir,0.1*acc);                                       /* time left on first GB */                                        
 dtrg=min(dtrg,trest1);
 
 Mcd=Mcgbf1b(M0,tvir);
 L2=Lgbf2(M0,Mcd);
 R2=Rgbf(Mt,L2);
 Mche2=Mcgbf1a(M0,tvir);
 perturb(Mt,Mche2,&L2,&R2,tvir,K);

 while(tvir<the1 || fabs(the1-tvir)<acc) {
   L1=L2;
   R1=R2;
   Mche1=Mche2;
   if(WIND==1)
     dM=dtrg*1.0e+6*windf1(Mt,L1,R1,Mche1,K);
   else if(WIND==2)
     dM=dtrg*1.0e+6*windf2(Mt,L1,R1,Mche1,K); 
   else if(WIND==3)
     dM=dtrg*1.0e+6*windf3(Mt,L1,R1,Mche1,K); 


   if(dM>(maxdM*Mt)) {                       /* check that mass loss do not go over 1% of star mass */
     dtrg=0.01*dtrg*Mt/dM;
     dM=0.01*Mt;
     if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
   }

   tvir+=dtrg;
   Mt-=dM;
   Mcd=Mcgbf1b(M0,tvir);
   L2=Lgbf2(M0,Mcd);
   R2=Rgbf(Mt,L2);
   Mche2=Mcgbf1a(M0,tvir); 
  
   if(Mche2>Mt) {                                  /* check that He core mass will not exceed star mass */
     Mt+=dM;                                       /* undo above part of loop */
     tvir-=dtrg;
     Mche2=Mt;                                     /* approximation: I set that Mc=Mt from previous and last step */
     mark1=1;
     break;
   }
   perturb(Mt,Mche2,&L2,&R2,tvir,K);

   if(fabs(R2-R1)>(maxdR*R1) && mark2==0) {              /* check that radius change no more then 10% */
     Mt+=dM;
     tvir-=dtrg;                           /* undo above part of loop */
     L2=L1;                                /* dM,Mcd changed via others changed param. */
     R2=R1;
     Mche2=Mche1;
     dtrg*=0.5;                            /* decrease by half time step */
     i++;
     if(i>5)                                /* make maximum 5 time step decreaments */
       mark2=1;
     if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
     continue;                              /* go directly from here back to begining of while() loop */
   }
   i=0;
   mark2=0;

   L=L2;
   R=R2;
   T=get_T(L,R);
   dt=dtrg;
   (*t)+=dtrg;
   (*M)=Mt;
   
   if(tvir<=tx)
     dtrg=delrg*(tinf1-tvir);
   else
     dtrg=delrg*(tinf2-tvir);
   dtrg=max(dtrg,0.1*acc);    
   trest1=max(the1-tvir,0.1*acc);                                    /* time left on GB */                                        
   dtrg=min(dtrg,trest1);   
   
   if(fabs(TE-(*t))<acc || mm==1) {                                         /* to stop calculation at TE */
     if(mm==1) TE=(*t); 
     (LRR)=L; (RRR)=R; MCRR=Mche2; TVIR=tvir; MPRE=0.0; DT=dtrg;
     (*stop)=1;
     (*Mhe)=Mche2;
     (*Mco)=0.0;
     break;                               
   }     
   trest2=max(TE-(*t),0.1*acc);
   dtrg=min(dtrg,trest2);
   if(fabs(the1-tvir)<acc) 
     break;                      /* to stop the loop after last step */
 } 
                                                       /* mark1==0 means that core mass is smaller then star mass */ 
                                                       /* mark1==1 core mass has reached star mass */             
 (*tstart)=0.0;                                        /* if Helium star is formed, it is on its ZAMS */ 

 if((*stop)!=1 || fabs(the1-tvir)<acc) {
   if(mark1==0) {                                        /* CHeB star */
     K=4;
     TVIR=0.0;
     if(M0<=M_HeF)
       DT=delhe*thef(Mt);
     else
       DT=delhe*thef(M0);  
   }
   else if(mark1==1 && M0>M_HeF) {                       /* ZAMS Helium star */
     K=7;
     TVIR=0.0;
   }  
   else if(mark1==1 && M0<=M_HeF) {                      /* formation of remnant */
     (*Mhe)=Mche2;
     (*Mco)=0.0;
     (*Kp)=K;
     K=-1;
     TVIR=0.0;
   }  
   else
     fprintf(fp0,"error: in HGf() star type K was not set\n");        
 }

 return K;
}


int HEf(double *M, double *t, double *tstart, double *Mhe, double *Mco, int *stop)
{ /* evolve star of M through Core Helium Burning Phase and writes the data, supplied current time *t */
  /* returns type after CHeB phase K: 7-helium MS star, 5-EAGB star */
  /* first output is already somewhat (in time) into CHeB */
  /* M0 is changed here for LM stars */
 double Mt,dM,Mche1,Mche2,Mcco;
 double the1,the,tvir,tend,dthe,trest1,trest2,dt; 
 double L,R,T,L1,L2,R1,R2;
 int K;
 int mark1=0;                                     /* to mark when core mass exceeds star mass */
 int i=0,mark2=0;
 int mm=0;
 
 Mt=(*M);    
 K=4;  
 Mcco=0.0;

 if(fabs(TVIR)<acc && M0<=M_HeF)                  /* change of M0 for LM star after He flash */    
   M0=(*M);                                       /* and only if present run does not come directly here */

 the1=the1f(M0);                                  /* virtual start of CHeB phase */
 the=thef(M0);                                    /* virtual and real lifetime of CHeB phase */
 tend=the1+the;                                   /* virtual end of CHeB phase */
 if(fabs(TVIR)<acc)
   tvir=the1;
 else
   tvir=TVIR;
 

 if(the<acc) {                                                /* ensures that star has CHeB phase */
   fprintf(fp0,"error: I should not be in HEf()\n");
   K=-10;
   return K;
 }
                                                              /* mark1=1 if Mcorehe>Mstar */
 mark1=CHeBf(Mt,tvir,&R2,&L2,&Mche2);                         /* fills L2,Mche2 for M0, Rd--dummy */
 if(mark1==0) {
   perturb(Mt,Mche2,&L2,&R2,tvir,K);
   if(DERIV==1)                                               /* derivative estimation run */
     dthe=TE-(*t); 
   else {
     dthe=min(delhe*the,TE-(*t));
     trest1=max(tend-tvir,0.1*acc);                           /* time left on CHeB */
     dthe=min(dthe,trest1);
   }
 }   
 
 while((tvir<tend || fabs(tend-tvir)<acc) && mark1==0) {
   L1=L2;
   R1=R2;
   Mche1=Mche2;
   if(WIND==1)
     dM=dthe*1.0e+6*windf1(Mt,L1,R1,Mche1,K);
   else if(WIND==2)
     dM=dthe*1.0e+6*windf2(Mt,L1,R1,Mche1,K); 
   else if(WIND==3)
     dM=dthe*1.0e+6*windf3(Mt,L1,R1,Mche1,K); 

   if(dM>(maxdM*Mt)) {                       /* check that mass loss do not go over 1% of star mass */
     dthe=0.01*dthe*Mt/dM;
     dM=0.01*Mt;
     if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
   }
                
   tvir+=dthe;
   Mt-=dM;

   mark1=CHeBf(Mt,tvir,&R2,&L2,&Mche2);                        /* fills L2,Mche2 for M0, Rd for Mt */
   if(Mche2>Mt || mark1==1) {                                  /* check that He core mass will not exceed star mass */
     Mt+=dM;                                                   /* undo above part of loop */  
     tvir-=dthe;
     Mche2=Mt;                                                 /* approximation: I set that Mc=Mt from previous and last step */ 
     mark1=1;
     break;
   }    
   perturb(Mt,Mche2,&L2,&R2,tvir,K);

   if(fabs(R2-R1)>(maxdR*R1) && mark2==0) {                    /* check that radius change no more then 10% */
     Mt+=dM;
     tvir-=dthe;                                               /* undo above part of loop */
     L2=L1;                                                    /* dM,Mcd changed via others changed param. */
     R2=R1;
     Mche2=Mche1;
     dthe*=0.5;                                                /* decrease by half time step */
     i++;
     if(i>5)                                                   /* make maximum 5 time step decreaments */
       mark2=1;
     if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;    
     continue;                                                 /* go directly from here back to begining of while() loop */
   }
   i=0;
   mark2=0;

   L=L2;
   R=R2;
   T=get_T(L,R);
   dt=dthe;
   (*t)+=dthe;
   (*M)=Mt;
   
   dthe=delhe*the;                                                      
   trest1=max(tend-tvir,0.1*acc);                                       /* time left on CHeB */
   dthe=min(dthe,trest1);
   if(fabs(TE-(*t))<acc || mm==1) {                                         /* to stop calculation at TE */
     if(mm==1) TE=(*t);
     (LRR)=L; (RRR)=R; MCRR=Mche2; TVIR=tvir; MPRE=0.0; DT=dthe;
     (*stop)=1;
     (*Mhe)=Mche2;
     (*Mco)=Mcco;
     break;                               
   } 
   trest2=max(TE-(*t),0.1*acc);
   dthe=min(dthe,trest2);
   if(fabs(tend-tvir)<acc) 
     break;                        /* to stop the loop after last step */
 }
                                                       /* mark1==0 means that core mass is smaller then star mass */ 
                                                       /* mark1==1 core mass has reached star mass */             
 if((*stop)!=1 || fabs(tend-tvir)<acc) {
   if(mark1==0) {                                        /* EAGB star */
     K=5;
     (*tstart)=0.0;
     TVIR=0.0;
     if(tbagbf(M0)<=txef(M0))
       DT=delag*(tinf1ef(M0)-tbagbf(M0));
     else
       DT=delag*(tinf2ef(M0)-tbagbf(M0));
   }
   else if(mark1==1) {                                  /* evolved MS Helium star */
     K=7;
     (*tstart)=thsmsf(Mt)*(tvir-the1)/the;
     TVIR=0.0;
   }
   else
     fprintf(fp0,"error: in HGf() star type K was not set\n");
 }

 return K;
}


int AGBf(double *M, double *t, double *tstart, double *Mhe, double *Mco, int *Kp, int *stop, int Klast)
{ /* evolve star of M through EAGB branch and TP AGB phase and writes the data, supplied current time *t */
  /* returns type after these phases K: 8-helium giant, 11-CO WD, 12-ONe WD, 13-NS, 14-BH, 15-NoRem */
  /* first output is already somewhat (in time) into EAGB */
 double Mt,Mcbagb,Mche,Mcco,Mcsn,dM,Mcdu;
 double tbagb,tvir,txe,tinf1e,tinf2e,dtagb,dt,tdu,trest1,trest2;
 double Mc1,Mc2,Mcd,Mcco0;
 double tinf1t,tinf2t,txt;
 double T,L,R,L1,R1,L2,R2;
 int K;
 int mark=0;                                                  /* to mark SN */  
 int mark1=0;                                                 /* to mark when core mass exceeds star mass */
 int mark3=0;                                                 /* to set tvir for K=6 */
 int i=0,mark2=0;
 int mm=0;

 Mt=(*M);
 tbagb=tbagbf(M0);                                            /* time when star gets to Base of AGB */ 
 Mcbagb=Mcbagbf(M0);                                          /* mass of He core at Base of AGB */
 Mche=Mcbagb;                                                 /* mass of He core througout EAGB */   
 Mcco0=Mceagbf1(M0,tbagb);                                    /* CO core mass at BAGB */
 Mcsn=Mcsnf(Mcbagb,Mcco0);                                    /* critical mass of CO/One core when SN occurs */ 

if(Klast==5) {

 K=5;
 if(fabs(TVIR)<acc)
   tvir=tbagb;
 else
   tvir=TVIR;

 txe=txef(M0);                                                  
 tinf1e=tinf1ef(M0);
 tinf2e=tinf2ef(M0);
 if(DERIV==1)                                /* derivative estimation run */
   dtagb=TE-(*t);
 else if(tvir<=txe)                                              
   dtagb=min(delag*(tinf1e-tvir),TE-(*t));
 else
   dtagb=min(delag*(tinf2e-tvir),TE-(*t));

 Mcco=Mceagbf1(M0,tvir);                                      /* starting CO core mass */
 L2=Leagbf2(M0,Mcco);
 R2=Ragbf(Mt,L2);
 perturb(Mt,Mche,&L2,&R2,tvir,K);                             /* Mche as He core to the H envelope */ 


 if(Mcbagb<=0.8) {                                            /* EAGB stops when Mcco reaches Mche */
   while(Mcco<Mche && Mcco<Mcsn) {                            /* Mcco grows, Mche constant: no 2DU, TPAGB */
     if(ECS==1 && Mcbagb>=Mcbur1 && Mcbagb<=Mcbur2 && Mcco>Mecs) {   /* ECS SNa -> NS formation */ 
       mark=2;
       break;
     }
     L1=L2;                                                  
     R1=R2;
     if(WIND==1)
       dM=dtagb*1.0e+6*windf1(Mt,L1,R1,Mche,K);              /* Mche as He core may blow wind through H envelope */
     else if(WIND==2)
       dM=dtagb*1.0e+6*windf2(Mt,L1,R1,Mche,K);
     else if(WIND==3)
       dM=dtagb*1.0e+6*windf3(Mt,L1,R1,Mche,K);

     if(dM>(maxdM*Mt)) {    
       dtagb=0.01*dtagb*Mt/dM;
       dM=0.01*Mt;
       if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
     }
    
     tvir+=dtagb;
     Mt-=dM;
     if(Mche>Mt) {                                   /* check that He core mass will not exceed star mass */
       Mt+=dM;                                       /* undo above part of loop */  
       tvir-=dtagb;
       Mche=Mt;                                      /* approximation: I set that Mc=Mt from previous and last step */ 
       mark1=1;
       break;
     }    

     Mcco=Mceagbf1(M0,tvir);
     L2=Leagbf2(M0,Mcco);
     R2=Ragbf(Mt,L2);
     perturb(Mt,Mche,&L2,&R2,tvir,K);        /* Mche as He core to the H envelope */ 

     if(fabs(R2-R1)>(maxdR*R1) && mark2==0) {              /* check that radius change no more then 10% */
       Mt+=dM;
       tvir-=dtagb;                          /* undo above part of loop */  
       Mcco=Mceagbf1(M0,tvir);               /* added, as it is used in a while loop above, and continue sends code there */
       L2=L1;                                /* dM changed via others changed param. */
       R2=R1;
       dtagb*=0.5;                           /* decrease by half time step */
       i++;
       if(i>5)                                /* make maximum 5 time step decreaments */
         mark2=1;
       if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;         /* used to be: if(mm==1) TE=(*t); but for: */
                                        /* some cases it would lkead to dt=0.0, and program would halt! */
       continue;                              /* go directly from here back to begining of while() loop */
     }
     i=0;
     mark2=0;
           
     L=L2;
     R=R2;
     T=get_T(L,R);
     dt=dtagb;
     (*t)+=dtagb;
     (*M)=Mt;
     if(tvir<=txe)                                              
       dtagb=max(delag*(tinf1e-tvir),0.1*acc);
     else
       dtagb=max(delag*(tinf2e-tvir),0.1*acc);     
     dtagb=max(dtagb,0.1*acc);
     if(fabs(TE-(*t))<acc || mm==1) {                                         /* to stop calculation at TE */
       if(mm==1) TE=(*t);
       (LRR)=L; (RRR)=R; MCRR=Mche; TVIR=tvir; MPRE=0.0; DT=dtagb;
       (*stop)=1;
       (*Mhe)=Mche;
       (*Mco)=Mcco;
       return K;
     }
     trest2=max(TE-(*t),0.1*acc);
     dtagb=min(dtagb,trest2);
   }
 }
 else if(Mcbagb<2.25) {                                            /* EAGB stops at tdu, when Mcco get Mdu */
   tdu=tduf(M0);                                                   /* tduf() valid only for 0.8<Mcbagb<2.25 */
   trest1=max(tdu-tvir,0.1*acc);                                                 /* time left on EAGB */
   dtagb=min(dtagb,trest1);
   Mcdu=Mcduf(M0);                                                 /* Mcdu is CO core mass at the time of 2DU */
   while((tvir<tdu || fabs(tdu-tvir)<acc) && Mcco<Mcsn) { /* Mcco grows but don't reach Mche-const: 2DU,TPAGB */        
     if(ECS==1 && Mcbagb>=Mcbur1 && Mcbagb<=Mcbur2 && Mcco>Mecs) {   /* ECS SNa -> NS formation */ 
       mark=2;
       break;
     }
     L1=L2;                                                         
     R1=R2;
     if(WIND==1)
       dM=dtagb*1.0e+6*windf1(Mt,L1,R1,Mche,K);              /* Mche as He core may blow wind through H envelope */
     else if(WIND==2)
       dM=dtagb*1.0e+6*windf2(Mt,L1,R1,Mche,K);
     else if(WIND==3)
       dM=dtagb*1.0e+6*windf3(Mt,L1,R1,Mche,K);
       
     if(dM>(maxdM*Mt)) {    
       dtagb=0.01*dtagb*Mt/dM;
       dM=0.01*Mt;
       if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
     }
    
     tvir+=dtagb;
     Mt-=dM;
     if(Mche>Mt) {                                   /* check that He core mass will not exceed star mass */
       Mt+=dM;                                       /* undo above part of loop */  
       tvir-=dtagb;
       Mche=Mt;                                      /* approximation: I set that Mc=Mt from previous and last step */ 
       mark1=1;
       break;
     }    

     Mcco=Mceagbf1(M0,tvir);
     L2=Leagbf2(M0,Mcco);
     R2=Ragbf(Mt,L2);
     perturb(Mt,Mche,&L2,&R2,tvir,K);                         /* Mche as He core to the H envelope */ 

     if(fabs(R2-R1)>(maxdR*R1) && mark2==0) {              /* check that radius change no more then 10% */
       Mt+=dM;
       tvir-=dtagb;                          /* undo above part of loop */
       Mcco=Mceagbf1(M0,tvir);
       L2=L1;                                /* dM changed via others changed param. */
       R2=R1;
       dtagb*=0.5;                           /* decrease by half time step */
       i++;
       if(i>5)                                /* make maximum 5 time step decreaments */
         mark2=1;
       if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
       continue;                              /* go directly from here back to begining of while() loop */
     }
     i=0;
     mark2=0;
           
     L=L2;
     R=R2;
     T=get_T(L,R);
     dt=dtagb;
     (*t)+=dtagb;
     (*M)=Mt;
     if(tvir<=txe)                                              
       dtagb=delag*(tinf1e-tvir);
     else
       dtagb=delag*(tinf2e-tvir);
     dtagb=max(dtagb,0.1*acc);
     trest1=max(tdu-tvir,0.1*acc);
     dtagb=min(dtagb,trest1);      
     if(fabs(TE-(*t))<acc || mm==1) {                                         /* to stop calculation at TE */
       if(mm==1) TE=(*t); 
       (LRR)=L; (RRR)=R; MCRR=Mche; TVIR=tvir; MPRE=0.0; DT=dtagb;
       (*stop)=1;
       (*Mhe)=Mche;
       (*Mco)=Mcco;
       return K;
     }
     trest2=max(TE-(*t),0.1*acc);
     dtagb=min(dtagb,trest2);
     if(fabs(tdu-tvir)<acc) 
       break;                            /* to stop the loop after last step */     
   }
   if(mark1!=1)                                       /* if mark==1 then star lost envelope before 2DU -so no 2DU */
     Mche=Mcdu;                                                    /* 2DU: He core reduced to Mcdu */ 
 }               
 else {                                                            /* EAGB stops when Mcco get to Msn: SN on EAGB */
   while(Mcco<Mche && Mcco<Mcsn) {      /* CO core ignites non-degenerately on EAGB */
                                                    /* CO core grows but do not rech Mche-const, no 2DU, no TPAGB */
     if(ECS==1 && Mcbagb>=Mcbur1 && Mcbagb<=Mcbur2 && Mcco>Mecs) {   /* ECS SNa -> NS formation */ 
       mark=2;
       break;
     }
     L1=L2;    
     R1=R2; 
     if(WIND==1)
       dM=dtagb*1.0e+6*windf1(Mt,L1,R1,Mche,K);              /* Mche as He core may blow wind through H envelope */
     else if(WIND==2)
       dM=dtagb*1.0e+6*windf2(Mt,L1,R1,Mche,K); 
     else if(WIND==3)
       dM=dtagb*1.0e+6*windf3(Mt,L1,R1,Mche,K); 

     if(dM>(maxdM*Mt)) {    
       dtagb=0.01*dtagb*Mt/dM;
       dM=0.01*Mt;
       if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
     }
    
     tvir+=dtagb;
     Mt-=dM;
     if(Mche>Mt) {                                   /* check that He core mass will not exceed star mass */
       Mt+=dM;                                       /* undo above part of loop */  
       tvir-=dtagb;
       Mche=Mt;                                      /* approximation: I set that Mc=Mt from previous and last step */ 
       mark1=1;
       break;
     }    

     Mcco=Mceagbf1(M0,tvir);
     L2=Leagbf2(M0,Mcco);
     R2=Ragbf(Mt,L2);
     perturb(Mt,Mche,&L2,&R2,tvir,K);                              /* Mche as He core to the H envelope */ 

     if(fabs(R2-R1)>(maxdR*R1) && mark2==0) {                                    /* check that radius change no more then 10% */
       Mt+=dM;
       tvir-=dtagb;                                                 /* undo above part of loop */
       Mcco=Mceagbf1(M0,tvir);
       L2=L1;                                                      /* dM changed via others changed param. */
       R2=R1;
       dtagb*=0.5;                                                 /* decrease by half time step */
       i++;
       if(i>5)                                /* make maximum 5 time step decreaments */
         mark2=1;
       if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
       continue;                              /* go directly from here back to begining of while() loop */
     }
     i=0;
     mark2=0;
           
     L=L2;
     R=R2;
     T=get_T(L,R);
     dt=dtagb;
     (*t)+=dtagb;
     (*M)=Mt;
     if(tvir<=txe)                                              
       dtagb=max(delag*(tinf1e-tvir),0.1*acc);
     else
       dtagb=max(delag*(tinf2e-tvir),0.1*acc);
     dtagb=max(dtagb,0.1*acc);  
     if(fabs(TE-(*t))<acc || mm==1) {                                         /* to stop calculation at TE */
       if(mm==1) TE=(*t); 
       (LRR)=L; (RRR)=R; MCRR=Mche; TVIR=tvir; MPRE=0.0; DT=dtagb;
       (*stop)=1;
       (*Mhe)=Mche;
       (*Mco)=Mcco;
       return K;
     }
     trest2=max(TE-(*t),0.1*acc);
     dtagb=min(dtagb,trest2);
   }
 }

 if(Mcco>Mt) mark1=1;                       /* CO core mass has exceeded star mass */
                                            /* the same condition is checked fo He core mass in above loops */ 

 if(Mcco>=Mcsn) mark=1;                     /* 1-regular SN */
                                            /* mark1==0 means that He core mass is smaller then star mass */
                                            /* mark1==1 He core mass has reached star mass */
 if(mark1==1) {                             /* Helium giant formed (K=8 or 9), but I assume always 8 */                                          
   K=8;                                     /* to pass control to apropriate function */
   (*tstart)=thsevolf(Mche,Mcco);           /* calculate how evolved is HS Giant */ 
   TVIR=0.0;
   return K;
 }
 else if(mark!=0) {                         /* SN occured: 0-no SN, 1-regular SN, 2-ECS SN */
   (*Mhe)=Mche;
   (*Mco)=Mcco;
   (*Kp)=K;
   K=-1;                                    /* formation of remnant */
   TVIR=0.0;
   return K;
 }  
 else 
   mark3=1;                                 /* TP AGB star formed */                                                                    

}


 if(mark3==0)                           
   tvir=TVIR;
        
 K=6;                                       /* calculates and writes TP AGB evolution */ 
 tinf1t=tinf1tf(M0);                        /* now there is one Mc=Mche=Mcco  =Mc2 */                          
 tinf2t=tinf2tf(M0);                        /* both He and CO cores grow together */ 
 txt=txtf(M0);
 if(DERIV==1)                                /* derivative estimation run */
   dtagb=TE-(*t);
 else if(tvir<=txt)
   dtagb=min(delag*(tinf1t-tvir),TE-(*t));
 else
   dtagb=min(delag*(tinf2t-tvir),TE-(*t));
 
 Mc2=Mctagbf1a(M0,tvir);                               
 Mcd=Mctagbf1b(M0,tvir);
 L2=Ltagbf2(M0,Mcd);
 R2=Ragbf(Mt,L2);
 perturb(Mt,Mc2,&L2,&R2,tvir,K);

 i=mark2=0;                                                 /* Mc2=Mcco=Mche grows until reaching Msn -> SN */ 
 while(Mc2<Mcsn) {              /* or until envelope is lost */

   if(Mzams>30.0 && (*t)>8.0) {   /* force to form remnant: BH */
     (*stop)=2;  
     break;
   }  

 
   if(ECS==1 && Mcbagb>=Mcbur1 && Mcbagb<=Mcbur2 && Mc2>Mecs) {   /* ECS SNa -> NS formation */ 
     mark=2;
     break;
   }
   L1=L2;                                                   
   R1=R2;
   Mc1=Mc2;
   if(WIND==1)
     dM=dtagb*1.0e+6*windf1(Mt,L1,R1,Mc1,K); 
   else if(WIND==2)
     dM=dtagb*1.0e+6*windf2(Mt,L1,R1,Mc1,K);
   else if(WIND==3)
     dM=dtagb*1.0e+6*windf3(Mt,L1,R1,Mc1,K);

   if(dM>(maxdM*Mt)) {                                          
     dtagb=0.01*dtagb*Mt/dM;
     dM=0.01*Mt;
     if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
   }
                           
   tvir+=dtagb;
   Mt-=dM;
   Mc2=Mctagbf1a(M0,tvir);
   if(Mc2>Mt) {                                    /* check that He/CO core mass will not exceed star mass */
     Mt+=dM;                                       /* undo above part of loop */
     tvir-=dtagb;
     Mc2=Mt;                                       /* approximation: I set that Mc=Mt from previous and last step */
     mark1=1;
     break;
   }
   
   Mcd=Mctagbf1b(M0,tvir);
   L2=Ltagbf2(M0,Mcd);
   R2=Ragbf(Mt,L2);
   perturb(Mt,Mc2,&L2,&R2,tvir,K);

   if(fabs(R2-R1)>(maxdR*R1) && mark2==0) {              /* check that radius change no more then 10% */
     Mt+=dM;
     tvir-=dtagb;                          /* undo above part of loop */
     L2=L1;                                /* dM,Mcd changed via others changed param. */
     R2=R1;
     Mc2=Mc1;
     dtagb*=0.5;                           /* decrease by half time step */
     i++;
     if(i>5)                                /* make maximum 5 time step decreaments */
       mark2=1;
     if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
     continue;                              /* go directly from here back to begining of while() loop */
   }
   i=0;
   mark2=0;

   L=L2;
   R=R2;
   T=get_T(L,R);
   dt=dtagb;
   (*t)+=dtagb;
   (*M)=Mt;
   if(tvir<=txt)
     dtagb=max(delag*(tinf1t-tvir),0.1*acc);
   else
     dtagb=max(delag*(tinf2t-tvir),0.1*acc);   
   dtagb=max(dtagb,0.1*acc);  
   if(fabs(TE-(*t))<acc || mm==1) {                                         /* to stop calculation at TE */
     if(mm==1) TE=(*t);
     (LRR)=L; (RRR)=R; MCRR=Mc2; TVIR=tvir; MPRE=0.0; DT=dtagb;
     (*stop)=1;
     (*Mhe)=Mc2;
     (*Mco)=Mc2;
     return K;
   }
   trest2=max(TE-(*t),0.1*acc);
   dtagb=min(dtagb,trest2);
 }
                                 /* mark1==0 means that He core mass=CO core mass is smaller then star mass */
                                 /* mark1==1 He core mass=CO core mass has reached star mass */
 (*tstart)=0.0;
 (*Mhe)=Mc2;
 (*Mco)=Mc2;
 (*Kp)=K;
 
 if((*stop)!=1) {
   K=-1;                                    /* always formation of remnant */
   TVIR=0.0;
 }   
    
 return K;           
} 


int HSMSf(double *M, double *t, double *tstart, double *Mhe, double *Mco, int *stop)
{ /* evolve Helium Star of M through Main Sequance and writes the data, supplied current *t */
  /* and virtual start time if star is already evolved Helium MS, first step is already into  HSMS */
  /* 1-means value from previous step, 2 from current step; p-star of deminished mass */
  /* tstart is the relative age of HS on its MS at this function call */
  /* t2p is virtual time here */
 double trest1,trest2,dt,t2p,dtms2p,tms2p,t1,dtms1,t2,tms2;
 double Mt,Mche,Mcco,dM,Mold,rat;
 double T,L,R,L2p,L1,R2p,R1;
 int K;
 int i=0,mark2=0;
 int mm=0;

 Mt=(*M);
 K=7;
 Mche=Mcco=0.0;                  /* although Mche has no meaning for Helium Star it is retained for output lines */
 
 if(fabs(TVIR)<acc) {
   t2p=(*tstart); 
   M0=Mt;                        /* change of M0 at begining of HS evolution, only if present run does not */
 }                               /* come directly here */  
 else
   t2p=TVIR;
        
 L2p=Lhsmsf(Mt,t2p);
 R2p=Rhsmsf(Mt,t2p);
 tms2p=thsmsf(Mt);
 if(DERIV==1)                                /* derivative estimation run */
   dtms2p=TE-(*t);
 else
   dtms2p=min(delms*tms2p,TE-(*t));

 while(t2p<tms2p || fabs(tms2p-t2p)<acc) {
   t1=t2p; 
   L1=L2p;
   R1=R2p;
   dtms1=dtms2p;
   if(WIND==1)
     dM=dtms1*1.0e+6*windf1(Mt,L1,R1,Mcco,K);
   else if(WIND==2)
     dM=dtms1*1.0e+6*windf2(Mt,L1,R1,Mcco,K);
   else if(WIND==3)
     dM=dtms1*1.0e+6*windf3(Mt,L1,R1,Mcco,K);
   
   if(dM>(maxdM*Mt)) {                       /* check that mass loss do not go over 1% of star mass */
     dtms2p=0.01*dtms2p*Mt/dM;
     dtms1=dtms2p;
     dM=0.01*Mt;
     if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1; 
   }  

   t2=t1+dtms1;
   tms2=tms2p;
   Mold=Mt;
   Mt-=dM;
   M0=Mt;
   
   tms2p=thsmsf(Mt);
   if(dM>acc)                /* mass loss from the star: rat not applicable: fresh fuel not taken away from the core */
     rat=1.0;
   else                      /* mass gain: rat applies, as fresh fuel added to the core */
     rat=Mold/Mt;                  /* He MS stars have convective core: Tout et al. 1997, MNRAS 291, 732 Sec. 5.1 */
   t2p=rat*(tms2p/tms2)*t2;      /* see also  Hurley et al 2002, Sec.2.6.6.1): rejuvanation included */

   L2p=Lhsmsf(Mt,t2p);
   R2p=Rhsmsf(Mt,t2p);
  
   if(fabs(R2p-R1)>(maxdR*R1) && mark2==0) {              /* check that radius change no more then 10% */ 
     Mt+=dM;
     M0=Mt;
     t2p=t1;                                /* undo above part of loop */
     tms2p=thsmsf(Mt);                      /* t1,L1,R1,tms2 unchanged */
     L2p=L1;                                /* dtms1,dM,t2 changed via others changed param. */ 
     R2p=R1;
     dtms2p*=0.5;                           /* decrease by half time step */
     i++;
     if(i>5)                                /* make maximum 5 time step decreaments */
       mark2=1;
     if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1; 
     continue;                              /* go directly from here back to begining of while() loop */
   }
   i=0;
   mark2=0;
     
   L=L2p;
   R=R2p;
   T=get_T(L,R);
   (*t)+=dtms2p;
   dt=dtms2p;
   (*M)=Mt;
   
   tms2p=thsmsf(Mt);
   dtms2p=delms*tms2p; 
   trest1=max(tms2p-t2p,0.1*acc);
   dtms2p=min(dtms2p,trest1);   
   if(fabs(TE-(*t))<acc || mm==1) {                                         /* to stop calculation at TE */
     if(mm==1) TE=(*t); 
     (LRR)=L; (RRR)=R; MCRR=Mcco; TVIR=t2p; MPRE=0.0; DT=dtms2p;
     (*stop)=1;
     (*Mhe)=Mche;
     (*Mco)=Mcco;
     break;
   }   
   trest2=max(TE-(*t),0.1*acc);
   dtms2p=min(dtms2p,trest2);
   if(fabs(tms2p-t2p)<acc) 
     break;                                 /* to stop the loop after last step */ 
 }                      
 
 if((*stop)!=1 || fabs(tms2p-t2p)<acc) {
   K=8;                                              /* always: unevolved Helium Giant */
   (*tstart)=thsmsf(M0);                             /* start at the end of Helium MS */
   TVIR=0.0;
   if(thsmsf(M0)<=txhsf(M0))                                           /* set the first time step */
     DT=delrg*(tinf1hsf(M0)-thsmsf(M0));
   else
     DT=delrg*(tinf2hsf(M0)-thsmsf(M0));  
 }
 
 return K;
} 


int HSGBf(double *M, double *t, double tstart, double *Mhe, double *Mco, int *Kp, int *stop)
{ /* evolve Helium Star of M through Giant Branches and writes the data, supplied current *t */
  /* first step is already into  HSGB */
  /* 1-means value from previous step, 2 from current step; p-star of deminished mass */
  /* tstart is the relative age of HS on its GB at this function call */
 double Mt,Mcco1,Mcco2,dM,Mcmax,Mcsn;      
 double tvir,tx,tinf1,tinf2,dtrg,dt,trest2;
 double T,L,R,L2,L1,R2,R1,Lths,Rzhs;       
 int Kold,K;                                           /* set here in function Rhsgbf(): either 8 or 9 */
 int i=0,mark2=0;
 int mm=0;       
        
 Mt=(*M);

 if(fabs(TVIR)<acc) {
   tvir=tstart;
   M0=Mt;                       /* if after HSMS ok, because Mt=*M=M0 at the end of HSMS */
 }                              /* if from any evolutionary functions, it is also ok, as Mt=*M=Mche and thus =M0 */
 else
   tvir=TVIR;                   /* here M0 from previous run passed in global variable here */

 tx=txhsf(M0);
 tinf1=tinf1hsf(M0);
 tinf2=tinf2hsf(M0);
  
 if(DERIV==1)                                /* derivative estimation run */     
   dtrg=TE-(*t);
 else if(tvir<=tx)                                           /* set the first time step */
   dtrg=min(delrg*(tinf1-tvir),TE-(*t));
 else
   dtrg=min(delrg*(tinf2-tvir),TE-(*t));  

 if(dtrg<0.0) {   /* star is depleted by RLOF to such a state that is has already degenerate structure */
   (*Mhe)=Mt;      /* e.g., M=0.6 and Mhe=Mco=0.6 for K=5 and should be send to REMNANTf to make CO WD */
   (*Mco)=Mt;          /* for such a star, the function used here do not work: e.g., timestep negative */
   (*Kp)=8;                                                   /* so it is treated as K=8 in REMNANTF() */
   K=-1;                                                               /* control passed to REMNANTf() */
   TVIR=0.0;
   return K;
 }

 Mcco2=Mchsgbf1(M0,tvir); 
 Mcmax=Mcmaxf(Mt);
 L2=Lhsgbf(M0,Mcco2);
 Lths=Lthsf(M0);                                 /* M0 constant, so no need to calculate this again */
 Rzhs=Rzhsf(Mt);
 R2=Rhsgbf(Mt,L2,Lths,Rzhs,&K);

 if((Mcco2+acc)>=Mcmax) {           /* formation of low mass remnant: CO WD */
   (*Mhe)=Mt;                       /* if Mcco2>Mcmax I can't pass this to perturb(): */
   (*Mco)=Mt;                       /* as then mi=negative and I get "Floating exception" error! */
   (*Kp)=K;                         /* no change in star mass, stop=0, no output as K=8 */
   K=-1;                            /* control passed to REMNANTf() */
   TVIR=0.0;                        /* setting here (*Mco)=Mt is approximation as He shell burning */
   return K;                        /* haven't converted all env. to CO! see. end of Sec. 6.1 */
 }
 perturb(Mt,Mcco2,&L2,&R2,tvir,K);

 while(1) {                                      /* stop is incorporated in the body of this loop */
   L1=L2;
   R1=R2;
   Mcco1=Mcco2;
   if(WIND==1)
     dM=dtrg*1.0e+6*windf1(Mt,L1,R1,Mcco1,K);
   else if(WIND==2)
     dM=dtrg*1.0e+6*windf2(Mt,L1,R1,Mcco1,K);
   else if(WIND==3)
     dM=dtrg*1.0e+6*windf3(Mt,L1,R1,Mcco1,K);

   if(dM>(maxdM*Mt)) {                            /* check that mass loss do not go over 1% of star mass */
     dtrg=0.01*dtrg*Mt/dM;
     dM=0.01*Mt;
     if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
   }

   tvir+=dtrg;
   Mt-=dM;
   Mcco2=Mchsgbf1(M0,tvir);
   L2=Lhsgbf(M0,Mcco2);   
   Kold=K;
   Rzhs=Rzhsf(Mt);
   R2=Rhsgbf(Mt,L2,Lths,Rzhs,&K);

   Mcmax=Mcmaxf(Mt);                               /* stop of He shell burning, for low mass stars */
   Mcsn=Mcsnf(M0,0.0);                           
   if(ECS==1 && M0>=Mcbur1 && M0<=Mcbur2 && Mcco2>Mecs) {   /* ECS SNa -> NS formation */ 
     break;
   }     
   if(Mcco2>Mt || Mcco2>Mcmax || Mcco2>Mcsn) {     /* check that CO core mass will not exceed star mass */
     Mt+=dM;                                       /* or Mcmax or Mcsn */
     tvir-=dtrg;                                   /* undo above part of loop */
     K=Kold;
     if(Mcco2>Mt)
       Mcco2=Mt;                                     /* approximation: I set that Mc=Mt from previous and last step */
     else if(Mcco2>Mcmax)
       Mcco2=Mt;
     else
       Mcco2=Mcco1;
     break;
   }
   perturb(Mt,Mcco2,&L2,&R2,tvir,K);

   if(fabs(R2-R1)>(maxdR*R1) && mark2==0) {              /* check that radius change no more then 10% */
     Mt+=dM;
     tvir-=dtrg;                           /* undo above part of loop */
     K=Kold;
     L2=L1;                                /* dM,Mcd changed via others changed param. */
     R2=R1;
     Mcco2=Mcco1;
     dtrg*=0.5;                            /* decrease by half time step */
     i++;
     if(i>5)                                /* make maximum 5 time step decreaments */
       mark2=1;
     if(MM==1 && fabs(DMMT)<acc && DERIV!=1) mm=1;
     continue;                              /* go directly from here back to begining of while() loop */
   }
   i=0;
   mark2=0;

   L=L2;
   R=R2;
   T=get_T(L,R);
   dt=dtrg;
   (*t)+=dtrg;
   (*M)=Mt;
   if(tvir<=tx)
     dtrg=max(delrg*(tinf1-tvir),0.1*acc);
   else
     dtrg=max(delrg*(tinf2-tvir),0.1*acc);     
   dtrg=max(dtrg,0.1*acc);  
   if(fabs(TE-(*t))<acc || mm==1) {                                         /* to stop calculation at TE */
     if(mm==1) TE=(*t);
     (LRR)=L; (RRR)=R; MCRR=Mcco2; TVIR=tvir; MPRE=0.0; DT=dtrg;
     (*stop)=1;
     (*Mhe)=Mt;
     (*Mco)=Mcco2;
     break;
   }     
   trest2=max(TE-(*t),0.1*acc);
   dtrg=min(dtrg,trest2);
 } 

 (*Mhe)=Mt;
 (*Mco)=Mcco2;
 (*Kp)=K;

 if((*stop)!=1) {
   K=-1;                      /* formation of remnant */
   TVIR=0.0;
 }  

 return K;
}


int REMNANTf(double *M, double *t, double *Mhe, double *Mco, int Kp)
{ /* calculates mass and type of Remnant and writes the data, supplied current time *t and mass *M */
  /* given He core mass Mhe, and CO core mass Mco, Kp must be type of star from which remnant forms */
 double Mt,Mc,dM,Mcbagb,L,R; 
 double Mbar,a,b,c,del;
 int K;

 dM=0.0;                /* no mass loss from any remnant assumed */
 Mc=0.0;                /* no core for remnants */
 L=acc;                 /* L set to eps (bigger then 0 to get its log), it may not be case for WD */

 if(Kp==0 || Kp==1 || Kp==4 || Kp==7 || Kp>9 || Kp<0)
   fprintf(fp0,"error: in REMNANT() remnant from star for which formation of remnant isn't allowed\n"); 

 else if(Kp==2 || Kp==3) {                        /* Helium White Dwarf */
   K=10; 
   Mt=(*Mhe);
   R=Rwdf(Mt,K);
 }
 
 else if(Kp==5 || Kp==6) {
   Mcbagb=Mcbagbf(M0);                                    /* M0 remains constant for K=5,6 */
   if(Mcbagb<Mcbur1 && (*Mco)<MCh) {                      /* Carbon Oxygen White Dwarf */
     K=11;
     Mt=(*Mhe);
     R=Rwdf(Mt,K);
   }
   else if(Mcbagb<Mcbur1 && (*Mco)>=MCh) {                /* Massless Supernova Remnant: SN IIa */
     K=15;
     Mt=0.0;
     R=acc;
     if(SN==1) sn_type(*M,*Mhe,*Mco,*t,K);                          /* sets SN type */
   }
   else if(ECS==1 && Mcbagb>=Mcbur1 && Mcbagb<=Mcbur2 && (*Mco)>=Mecs) {        /* ECS SNa -> NS */
     Mt=Mnsecs;
     K=13;
     R=Rnsf(); 
     FB=0;
     FRAC=0.0;                        /* writes how much fallback [0:1] took place into global var. */
     ECSSN=1;
     if(SN==1) sn_type(*M,*Mhe,*Mco,*t,K);                          /* sets SN type */
   }  
   else if(Mcbagb>=Mcbur1 && Mcbagb<=Mcbur2 && (*Mco)<Mcocrit) {                /* CO WD */
     K=11;
     Mt=(*Mhe);
     R=Rwdf(Mt,K);
   }
   else if(Mcbagb>=Mcbur1 && Mcbagb<=Mcbur2 && (*Mco)>=Mcocrit && (*Mco)<MCh) { /* ONe WD */
     K=12;
     Mt=(*Mhe);
     R=Rwdf(Mt,K);
   }
   else if((Mcbagb>=Mcbur1 && Mcbagb<=Mcbur2 && (*Mco)>=MCh) || Mcbagb>Mcbur2 ) {  /* SN II */
     if(REMNANT==0)
       Mt=compactf0(*M,*Mco);                      /* returns gravitational mass of compact object */
     else if(REMNANT==1)
       Mt=compactf1(*M,*Mco);
     else if(REMNANT==2)
       Mt=compactf2(*M,*Mco);
     else if(REMNANT==3)
       Mt=compactf3(*M,Mco,Mhe);
     else if(REMNANT==4)
       Mt=compactf4(*M,Mco,Mhe);
     else if(REMNANT==5)
       Mt=compactf5(*M,Mco,Mhe);
     else if(REMNANT==6)
       Mt=compactf6(*M,Mco,Mhe);
     else if(REMNANT==7)
       Mt=compactf7(*M,Mco,Mhe);
     else if(REMNANT==8)
       Mt=compactf8(*M,Mco,Mhe);
     else if(REMNANT==9)
       Mt=compactf9(*M,Mco,Mhe);
     else if(REMNANT==10)
       Mt=compactf10(*M,Mco,Mhe);
     else if(REMNANT==12)
       Mt=compactf12(*M,Mco,Mhe);
     else if(REMNANT==13)
       Mt=compactf13(*M,Mco,Mhe);
     else if(REMNANT==14)
       Mt=compactf14(*M,Mco,Mhe);
     else 
       fprintf(fp0,"error: unknown choice of remnant mass function in singl.c: REMNANT: %d\n",REMNANT);

     if(Mt<0.001) {                                    /* pair instability SN: no BH */
         R=0.0;    
         K=18;
     }   
     else if(Mt<=Mmaxns) {                             /* Neutron Star */
       K=13; 
       R=Rnsf();
     }                  
     else {                                            /* Black Hole */                          
       R=Rbhf(Mt);   
       K=14;
     }     
     if(SN==1) sn_type(*M,*Mhe,*Mco,*t,K);                    /* sets SN type */ 
   }
   else {
     fprintf(fp0,"unexpected remnant case for K=5-6: %d\n",iidd_old);
     fflush(fp0);
   }
 }

 else {                                                  /* Kp=8 or Kp=9 */
   if(M0<Mcbur1 && (*Mco)<MCh) {                      /* Carbon Oxygen White Dwarf */
     K=11;
     Mt=(*Mco);
     R=Rwdf(Mt,K);
   }
   else if(M0<Mcbur1 && (*Mco)>=MCh) {                /* Massless Supernova Remnant: SN IIa */
     K=15;
     Mt=0.0;
     R=acc;
     if(SN==1) sn_type(*M,*Mhe,*Mco,*t,K);                       /* sets SN type */
   }
   else if(ECS==1 && M0>=Mcbur1 && M0<=Mcbur2 && (*Mco)>=Mecs) {        /* ECS SNa -> NS */
     Mt=Mnsecs;
     K=13;
     R=Rnsf();
     FB=0; 
     FRAC=0.0;                        /* writes how much fallback [0:1] took place into global var. */
     ECSSN=1;
     if(SN==1) sn_type(*M,*Mhe,*Mco,*t,K);                       /* sets SN type */
   }  
   else if(M0>=Mcbur1 && M0<=Mcbur2 && (*Mco)<Mcocrit) {                /* CO WD */
     K=11;
     Mt=(*Mco);
     R=Rwdf(Mt,K);
   }
   else if(M0>=Mcbur1 && M0<=Mcbur2 && (*Mco)>=Mcocrit && (*Mco)<MCh) { /* ONe WD */
     K=12;
     Mt=(*Mco);
     R=Rwdf(Mt,K);
   }
   else if((M0>=Mcbur1 && M0<=Mcbur2 && (*Mco)>=MCh) || M0>Mcbur2 ) {  /* SN II */
     if(REMNANT==0)
       Mt=compactf0(*M,*Mco);                      /* returns gravitational mass of compact object */
     else if(REMNANT==1)
       Mt=compactf1(*M,*Mco);
     else if(REMNANT==2)
       Mt=compactf2(*M,*Mco);
     else if(REMNANT==3)
       Mt=compactf3(*M,Mco,Mhe);
     else if(REMNANT==4)
       Mt=compactf4(*M,Mco,Mhe);
     else if(REMNANT==5)
       Mt=compactf5(*M,Mco,Mhe);
     else if(REMNANT==6)
       Mt=compactf6(*M,Mco,Mhe);
     else if(REMNANT==7)
       Mt=compactf7(*M,Mco,Mhe);
     else if(REMNANT==8)
       Mt=compactf8(*M,Mco,Mhe);     
     else if(REMNANT==9)
       Mt=compactf9(*M,Mco,Mhe);
     else if(REMNANT==10)
       Mt=compactf10(*M,Mco,Mhe);     
     else if(REMNANT==12)
       Mt=compactf12(*M,Mco,Mhe);
     else if(REMNANT==13)
       Mt=compactf13(*M,Mco,Mhe);
     else if(REMNANT==14)
       Mt=compactf14(*M,Mco,Mhe);
     else 
       fprintf(fp0,"error: unknown choice of remnant mass function in singl.c: REMNANT: %d\n",REMNANT);

     if(Mt<0.001) {                                    /* pair instability SN: no BH */
         R=0.0;    
         K=18;
     }   
     else if(Mt<=Mmaxns) {                             /* Neutron Star */
       K=13; 
       R=Rnsf();
     }                  
     else {                                            /* Black Hole */                          
       R=Rbhf(Mt);   
       K=14;
     }
     if(SN==1) sn_type(*M,*Mhe,*Mco,*t,K);                    /* sets SN type */ 
   }
   else {
     fprintf(fp0,"unexpected remnant case for K=8-9: %d\n",iidd_old);
     fflush(fp0);
   }
 }  


 if(newCONe==1 && (K==11 || K==12)) {
   if(K==11 && Mt>MCONe) K=12;
   if(K==12 && Mt<MCONe) K=11; 
 }

 (LRR)=L; (RRR)=R; MCRR=Mc; TVIR=0.0; MPRE=(*M);  
 DT=1.0e+50;                                             /* big number to avoid its usage in binary.c */   

 (*M)=M0=Mt;                                             /* change of M0 for every remnant */ 

 
 return K;
}


void sn_type(double M, double Mhe, double Mco, double t, int K)
{ /* increases global counters for different SN types */
  /* increases global counters for fraction of NS formed in different SNs */ 
  /* counts only SNs which occured in t_hubble time */    
  /* FB passed to here from compactf() through REMNANT() function */ 
  /* FB=0 -no fall back, FB=1 -partial fall back, FB=2 -complete fall back */
  /* IMPORTANT: call this function once K type of remnant is already set in REMNANT() */     
 int i;     
     
 i=(int)floor((Tst+t)/1000.0);    
      
 if(K!=18 && DERIV==0) {          /* do not count when singl() is used for evaluation of derivs */
                                  /* do not count pair instability suprnovae K=18 */  
      
   if(K==15) {                    /* SN IIa: CO core degenerate flash, no compact object formed */
     ss0++;                       
   }     
   else if(fabs(M-Mco)<acc) {     /* SN of type Ic: CO-rich env */
     if(FB==0) { 
       if(K==13) ns1a++;          /* fraction of NS formed in all clean Ic explosions */
       ss1a++;
       sIbc[i]++;
     }
     else if(FB==1) {
       if(K==13) ns1b++;                  
       ss1b++;
       sIbc[i]++;
     }  
     else if(FB==2) {
       if(K==13) ns1c++;
       ss1c++;
       sFa[i]++;
     }  
   }
   else if(fabs(M-Mhe)<acc) {     /* SN of type Ib: He-rich env */
     if(FB==0) {
       if(K==13) ns2a++;
       ss2a++;
       sIbc[i]++;
     }  
     else if(FB==1) {
       if(K==13) ns2b++;                  
       ss2b++;
       sIbc[i]++;
     }  
     else if(FB==2) {
       if(K==13) ns2c++;
       ss2c++;
       sFa[i]++;
     }   
   }
   else if((M-Mhe)<5.0) {         /* SN of type II-L: H-rich env, Menv<5Msun */
     if(FB==0) {
       if(K==13) ns3a++;
       ss3a++;
       sII[i]++;
     }  
     else if(FB==1) {
       if(K==13) ns3b++;                  
       ss3b++;
       sII[i]++;
     }  
     else if(FB==2) {
       if(K==13) ns3c++;
       ss3c++;
       sFa[i]++;
     }  
   } 
   else if((M-Mhe)>=5.0) {        /* SN of type II-P: H-rich env, Menv>=5Msun */
     if(FB==0) {
       if(K==13) ns4a++; 
       ss4a++;
       sII[i]++;
     }  
     else if(FB==1) {
       if(K==13) ns4b++;                  
       ss4b++;
       sII[i]++;
     }  
     else if(FB==2) {
       if(K==13) ns4c++;
       ss4c++;
       sFa[i]++;
     }  
   }
   else
     fprintf(fp0,"error: SN of unknown type!\n");          
 } 
}  


double windf1(double M, double L, double R, double Mc, int K)
{ /* old winds: 2000 (jarrod hurley) */
  /* calculates and returns mass loss dMs[Msun] for star of type K with (M,L,R,Mc) per one year */ 
  /* returned dMs is POSITIVE and per year, not per Myr!!!, for stars with K>9 dMs=0 always */
  /* if global varaible mflag=0 then returns always 0.0 -- no mass loss */
  /* Reimers coefficient ETA is global variable */
  /* calculated wind mass loss (Hurley) is multiplied by WIND1/2 parameters set in sinbin.h */ 
 double dMs,dMl,dMt,P0,x,T;
 double E,L0,kappa,mi,fz;

 if(DERIV==1 && fabs(DMMT)<acc)                       /* run for evaluation of dlnR/dt at constant mass */ 
   return 0.0;   
 else if(DERIV==1 && fabs(DMMT)>acc)                  /* run for evaluation of dlnR/dlnM at constant comp. */
   return DMMT;                                       /* do not include wind mass loss!!! */
 else if(mflag==0)                                    /* STOPS HERE FOR NO MASS LOSS, use only MT rate */ 
   return DMMT;
                                                 /* MS stars if get mass loss, they get it only here as dMs */
 if(L>4000.0 || K==1) {                          /* mass loss for massive stars, over entire HR diagram */
   if(K==1) x=1.0;
   else x=min(1.0,(L-4000.0)/500.0);
   dMs=9.6e-15*x*pow(R,0.81)*pow(L,1.24)*pow(M,0.16);
   dMs=dMs*pow(ZZ/0.02,0.5);
 }
 else 
   dMs=0.0;
   
 if(K>=2 && K<=9) {
   E=4.0e-13;
   dMl=ETA*E*L*R/M;
   if(K==5 || K==6) {                            /* enhanced mass loss for AGB stars */ 
     P0=-2.07-0.9*log10(M)+1.94*log10(R);
     P0=pow(10.0,P0);
     P0=min(P0,2000.0);
     dMt=-11.4+0.0125*(P0-100.0*max(M-2.5,0.0));
     dMt=pow(10.0,dMt);
     dMt=min(dMt,1.36e-9*L);
     dMl=max(dMl,dMt);
   }
 
   if(K>6) {                                       /* mass loss for naked helium stars (WR stars) */
     if(ZZ/ZZsun>=10) 
       fz=pow(10.0,Z11);
     else if(ZZ/ZZsun>=0.001 && ZZ/ZZsun<10)
       fz=pow(ZZ/ZZsun,Z11);                             /* Z-dependence from Vink & de Koter 2005 */
     else   
       fz=pow(0.001,Z11);
     dMs=max(dMl,1.0e-13*pow(L,1.5)*fz);
   }
   else {                                        
     dMs=max(dMl,dMs);                          
     L0=7.0e+4;
     kappa=-0.5;
     mi=((M-Mc)/M)*min(5.0,max(1.2,pow(L/L0,kappa))); 
     if(mi<1.0) {                                /* mass loss for stars with small H-envelope, WR shines through */
       if(ZZ/ZZsun>=10) 
         fz=pow(10.0,Z11);
       else if(ZZ/ZZsun>=0.001 && ZZ/ZZsun<10)
         fz=pow(ZZ/ZZsun,Z11);                             /* Z-dependence from Vink & de Koter 2005 */
       else   
         fz=pow(0.001,Z11);
       dMl=1.0e-13*pow(L,1.5)*fz*(1.0-mi);
       dMs=max(dMl,dMs);
     } 
     x=1.0e-5*R*sqrt(L);                         
     if(L>6.0e+5 && x>1.0) {                          /* LBV mass loss for stars beyond Humphreys-Davidson limit */
       dMl=0.1*pow(x-1.0,3.0)*(L/6.0e+5-1.0);
       dMs=dMs+dMl;
     } 
   }
 } 
   
 if(K==7 || K==8 || K==9) dMs*=WIND2;
 else dMs*=WIND1;    

 return dMs+DMMT;              /* DMMT -- MT rate/gain in [Msun/yr] for binaries from binary.c */
}


double windf2(double M, double L, double R, double Mc, int K)
{ /* updated winds: 2009 (Jorric Vink): FIRST ATTEMPT */ 
  /* returns mass loss in [Msun/yr] */
 double dM,dMl,dMs,dMt,mi,E,P0,L0,kappa,T,fz,x;


 if(DERIV==1 && fabs(DMMT)<acc)                       /* run for evaluation of dlnR/dt at constant mass */ 
   return 0.0;   
 else if(DERIV==1 && fabs(DMMT)>acc)                  /* run for evaluation of dlnR/dlnM at constant comp. */
   return DMMT;                                       /* do not include wind mass loss!!! */
 else if(mflag==0)                                    /* STOPS HERE FOR NO MASS LOSS, use only MT rate */ 
   return DMMT;

 T=get_T(L,R);
 x=1.0e-5*R*sqrt(L);                                   /* LBV star: over Humphreys-Davidson limit: Hurley et al. 2000 */  
 if(K>=0 && K<7 && L>6.0e+5 && x>1.0) {                /* LBV: only for H-rich stars: Vink & de Koter 2002 */ 
   dM=Flbv*1.0e-04;                                    /* range to check: 1.0e-03 -- 1.0e-05 */
 }
 else if(K>=0 && K<7 && T<=12500.0) {                  /* H rich stars (low mass): Hurley et al.  2000 winds */
   if(L>4000.0 || K==1) {                              /* MS stars if get mass loss, they get it only here as dMs */
     if(K==1) x=1.0;                                   /* mass loss for massive stars, over entire HR diagram */
     else x=min(1.0,(L-4000.0)/500.0);
     dMs=9.6e-15*x*pow(R,0.81)*pow(L,1.24)*pow(M,0.16);
     dMs=dMs*pow(ZZ/0.02,0.5);
   }
   else 
     dMs=0.0;
   
   if(K>=2) {
     E=4.0e-13;
     dMl=ETA*E*L*R/M;
     if(K==5 || K==6) {                            /* enhanced mass loss for AGB stars */ 
       P0=-2.07-0.9*log10(M)+1.94*log10(R);
       P0=pow(10.0,P0);
       P0=min(P0,2000.0);
       dMt=-11.4+0.0125*(P0-100.0*max(M-2.5,0.0));
       dMt=pow(10.0,dMt);
       dMt=min(dMt,1.36e-9*L);
       dMl=max(dMl,dMt);
     }
     dMs=max(dMl,dMs);                          
     L0=7.0e+4;
     kappa=-0.5;
     mi=((M-Mc)/M)*min(5.0,max(1.2,pow(L/L0,kappa))); 
     if(mi<1.0) {                                /* mass loss for stars with small H-envelope, WR shines through */
       if(ZZ/ZZsun>=10) 
         fz=pow(10.0,Z11);
       else if(ZZ/ZZsun>=0.001 && ZZ/ZZsun<10)
         fz=pow(ZZ/ZZsun,Z11);                             /* Z-dependence from Vink & de Koter 2005 */
       else   
         fz=pow(0.001,Z11);
       dMl=1.0e-13*pow(L,1.5)*fz*(1.0-mi);
       dMs=max(dMl,dMs);
     } 
   } 
   dM=dMs;
 }
 else if(K>=0 && K<7) {                                /* H-rich stars (high mass): Vink et al. 2001 */
          
   if(T>22500.0 && T<25000.0) T=22500.0;               /* bistability jump: ignored! */
   else if(T>25000.0 && T<27500.0) T=27500.0;          /* star Temp. shifted to the edge of the jump */
   else if(T>50000.0) T=50000.0;

   if(T>=12500.0 && T<=22500.0)
     dM=-6.688+2.210*log10(L/1.0e+05)-1.339*log10(M/30.0)-1.601*log10(1.3/2.0)+1.07*log10(T/20000.0)+Z22*log10(ZZ/ZZsun);
   else if(T>=27500.0 && T<=50500.0)
     dM=-6.697+2.194*log10(L/1.0e+05)-1.313*log10(M/30.0)-1.226*log10(2.6/2.0)+0.933*log10(T/40000.0)-10.92*pow(log10(T/40000.0),2.0)+Z33*log10(ZZ/ZZsun);
   dM=pow(10.0,dM);
 }
 else if(K==7 || K==8 || K==9) {                       /* naked helium stars: */
   if(ZZ/ZZsun>=10) 
     fz=pow(10.0,Z11);
   else if(ZZ/ZZsun>=0.001 && ZZ/ZZsun<10.0)
     fz=pow(ZZ/ZZsun,Z11);                             /* Z-dependence from Vink & de Koter 2005 */
   else   
     fz=pow(0.001,Z11);
   if(WRWINDS==1)
     dM=1.0e-13*pow(L,1.5)*fz;                         /* Hamman & Koesterke 1998 + Vink, de Koter 2005 */
   else if(WRWINDS==2)
     dM=1.9e-05*pow(M/14.7,2.93)*fz;                   /* Zdziarski, Mikolajewska, Belczynski 2012 + Vink, de Koter 2005 */ 
   else if(WRWINDS==3)
     dM=1.9e-05*pow(M/14.7,2.93);                      /* Zdziarski, Mikolajewska, Belczynski 2012 */
   else 
     fprintf(fp0,"error: unknown wind type for WR stars\n");        
 }  
 else 
   fprintf(fp0,"error: you are out of stellar type range: windf2()\n");   
   

 if(K==7 || K==8 || K==9) dM*=WIND2;
 else dM*=WIND1;    

 return dM+DMMT;                  /* DMMT -- MT rate/gain in [Msun/yr] for binaries from binary.c */
}


double windf3(double M, double L, double R, double Mc, int K)
{ /* updated winds: Nov 21, 2022 (Alex Gormaz): FIRST ATTEMPT */ 
  /* returns mass loss in [Msun/yr] */
 double dM,dMl,dMs,dMt,mi,E,P0,L0,kappa,T,Tkk,fz,x;
 double gs,lgs,ww,xx,yy,zz;
 double dM10,L10,A10,dMhe1,dMhe2;

 if(DERIV==1 && fabs(DMMT)<acc)                        /* run for evaluation of dlnR/dt at constant mass */ 
   return 0.0;   
 else if(DERIV==1 && fabs(DMMT)>acc)                   /* run for evaluation of dlnR/dlnM at constant comp. */
   return DMMT;                                        /* do not include wind mass loss!!! */
 else if(mflag==0)                                     /* STOPS HERE FOR NO MASS LOSS, use only MT rate */ 
   return DMMT;

 T=get_T(L,R);                                         /* star effective temperature [K] */
 Tkk=T/1000.0;                                         /* star effective temperature [kK] in kilo  Kelvins */
 x=1.0e-5*R*sqrt(L);                                   /* LBV star: over Humphreys-Davidson limit: Hurley et al. 2000 */  
 gs=(GGGG*M*Msun)/(R*R*Rsun*Rsun);                     /* star surface gravity accelaration in cgs units [cm s^(-2)] */
 lgs=log10(gs); 

 if(K>=0 && K<7 && L>6.0e+5 && x>1.0 && T<28000.0) {   /* LBV: only for H-rich stars: Vink & de Koter 2002 */ 
   dM=Flbv*1.0e-04;                                    /* range to check: 1.0e-03 -- 1.0e-05 */
 }
 else if(K>=0 && K<7 && T<=12500.0 && lgs<3.2) {       /* cool H-rich stars: low temp. and low gravity (small size): Hurley2000 winds */
   if(L>4000.0 || K==1) {                              /* MS stars if get mass loss, they get it only here as dMs */
     if(K==1) x=1.0;                                   /* mass loss for massive stars, over entire HR diagram */
     else x=min(1.0,(L-4000.0)/500.0);
     dMs=9.6e-15*x*pow(R,0.81)*pow(L,1.24)*pow(M,0.16);
     dMs=dMs*pow(ZZ/0.02,0.5);
   }
   else 
     dMs=0.0;
   
   if(K>=2) {
     E=4.0e-13;
     dMl=ETA*E*L*R/M;
     if(K==5 || K==6) {                            /* enhanced mass loss for AGB stars */ 
       P0=-2.07-0.9*log10(M)+1.94*log10(R);
       P0=pow(10.0,P0);
       P0=min(P0,2000.0);
       dMt=-11.4+0.0125*(P0-100.0*max(M-2.5,0.0));
       dMt=pow(10.0,dMt);
       dMt=min(dMt,1.36e-9*L);
       dMl=max(dMl,dMt);
     }
     dMs=max(dMl,dMs);                          
     L0=7.0e+4;
     kappa=-0.5;
     mi=((M-Mc)/M)*min(5.0,max(1.2,pow(L/L0,kappa))); 
     if(mi<1.0) {                                /* mass loss for stars with small H-envelope, WR shines through */
       if(ZZ/ZZsun>=10) 
         fz=pow(10.0,Z11);
       else if(ZZ/ZZsun>=0.001 && ZZ/ZZsun<10)
         fz=pow(ZZ/ZZsun,Z11);                             /* Z-dependence from Vink & de Koter 2005 */
       else   
         fz=pow(0.001,Z11);
       dMl=1.0e-13*pow(L,1.5)*fz*(1.0-mi);
       dMs=max(dMl,dMs);
     } 
   } 
   dM=dMs;
 }
 else if(K>=0 && K<7) {                 /* hot H-rich stars: high mass with small size (T>12500K and lgs>=3.2): Alex formulae */
   if(T<25000.0) {                      /* 12500K<T<25000k */
      if(T>22500.0 && T<25000.0) T=22500.0;               /* bistability jump ignored: star Temp. shifted to the edge of the jump */
     dM=-6.688+2.210*log10(L/1.0e+05)-1.339*log10(M/30.0)-1.601*log10(1.3/2.0)+1.07*log10(T/20000.0)+Z22*log10(ZZ/ZZsun);
     dM=pow(10.0,dM);                   /* [Msun/yr] */
   } 
   else {                               /* T>25000K */
     ww=log10(Tkk); 
     xx=1.0/lgs;
     yy=R;
     zz=log10(ZZ/ZZsun);
     dM=-40.314+15.438*ww+45.838*xx-8.284*ww*xx+1.0564*yy-(ww*yy)/2.36-1.1967*xx*yy+11.6*zz-4.223*ww*zz-16.377*xx*zz+(yy*zz)/81.735;
     dM=pow(10.0,dM);                   /* [Msun/yr] Gormaz-Matamala et al. 2022 */
   }
 }
 else if(K==7 || K==8 || K==9) {        /* all naked helium stars: any mass, any luminosity, any temperature: Alex formulae */
   zz=log10(ZZ/ZZsun);
   dM10=-0.75*zz-4.06;
   dM10=pow(10.0,dM10);                 /* constant [Msun/yr] */
   L10=-0.87*zz+5.06;
   L10=pow(10.0,L10);                   /* constant [Lsun] */
   A10=0.32*zz+1.4;
   dMhe1=pow(10.0,-13.3)*pow(L,1.36)*pow(ZZ/ZZsun,0.61);         /* [Msun/yr] Vink 2017: low mass He stars */
   dMhe2=dM10*pow(log10(L/L10),A10)*pow(L/(10.0*L10),0.75);      /* [Msun/yr] Sander et al. 2020: high mass He stars */
   if(L<L10)   
     dM=dMhe1;                                                   /* low-luminosity naked helium stars */
   else 
     dM=max(dMhe1,dMhe2);                                        /* high-luminosity helium stars */
 }  
 else 
   fprintf(fp0,"error: you are out of stellar type range: windf3()\n");   
   

 if(K==7 || K==8 || K==9) dM*=WIND2;
 else dM*=WIND1;    

 return dM+DMMT;                        /* DMMT -- MT rate/gain in [Msun/yr] for binaries from binary.c */
}


int perturb(double M, double Mc, double *L, double *R, double t, int K)
{ /* change luminosity *L and radius *R due to small envelope mass, sec.6.3 */
  /* uses hrdiag.f l.626,628 which defines what shall happen for K=5 if remnant is K=8 */
  /* changes done for star of M, Mc, type K at time t (time needed only for K=4!!!) */ 
  /* returns 0 if *L,*R not changed, otherwise returns 1 */
  /* DIRECTLY USES M0 -- FOR MASS LOSS!!! */
 double tbagb,the1,the,tau,L0,Lc,Lzhs,Lths,Rc,Rzhs,Mcmax,Mche,Mcco;
 double alfa,beta,mi,kappa,b,c,q,s,r;
 double Mdum;
 double lRmax,Rmax,lRmax100,Rmax100;
 int returnT;


 if(K==2 || K==3 || K==4 || K==5 || K==6) {
   L0=7.0e+4;
   kappa=-0.5;
   mi=((M-Mc)/M)*min(5.0,max(1.2,pow((*L)/L0,kappa))); 
 }
 else if(K==8 || K==9) {
   Mcmax=Mcmaxf(M);
   mi=5.0*((Mcmax-Mc)/Mcmax);
 }  
 else
   return 0;                        /* no perturbation needed: no clear core envelope division */
     
 
 if(mi>=1.0)
   returnT=0;                        /* no perturbation needed: envelope two big */
 else {
   if(K==2 || K==3) {
     if(M0>M_HeF) {
       Lc=Lzhsf(Mc);
       Rc=Rzhsf(Mc);
     }
     else {
       Lc=Lwdf(Mc,0.0,10);
       Rc=Rwdf(Mc,10);
     }
   }
   else if(K==4) {
     the1=the1f(M0);
     the=thef(M0);     
     tau=(t-the1)/the;
     alfa=max(0.0,0.85-0.08*Mc);
     beta=max(0.0,0.4-0.22*log10(Mc));
     Lzhs=Lzhsf(Mc);
     Rzhs=Rzhsf(Mc);
     Lc=Lzhs*(1.0+0.45*tau+alfa*tau*tau);
     Rc=Rzhs*(1.0+beta*tau-beta*pow(tau,6.0));            
   }
   else if(K==5) {
     tbagb=tbagbf(M0);
     tau=3.0*(t-tbagb)/(tnf(M0,K)-tbagb);
     Mche=Mcbagbf(M0);
     Mcco=Mceagbf2(M0,(*L));
     Lc=Lhsgbf(Mche,Mcco);                       /* from paper0, luminosity of Helium Giant */
     if(tau<1.0) {                               /* from hrdiag.f, looks like Luminosity of Helium Hartprung Gap */
       Lths=Lthsf(Mche); 
       Lc=Lths*pow(Lc/Lths,tau);
     }     
     Lths=Lthsf(Mche);                           /* Mche is M0 for the hipotetical Helium Star here */ 
     Rzhs=Rzhsf(Mche);
     Rc=Rhsgbf(Mche,Lc,Lths,Rzhs,&K);            /* K dummy variable here, changed to 8 or 9, */
   }                                             /* but only in the scope of this func */   
   else {                                        /* they assume that for K=6,8,9, remnant will be CO WD: K=11 */
     if(Mc>=MCh) {                               /* WATCH!!!: perturbation as applied is only my approximation! */
       Mdum=1.0;                                 /* I use Mdum=1 Msun: it is arbitrary value to result in a star */ 
       Lc=Lwdf(Mdum,0.0,11);                     /* radius of 0.008058 Rsun (5500 km). i one inputs Mdum=MCh -> */
       Rc=Rwdf(Mdum,11);                         /* then radius of a star gets to NS radius: 0.000014 Rsun (10km) */
     }                                           /* as Rwdf=max(Rns,form1(M)) and form1=0.0 for M=MCh */
     else {                                      /* Jarrod don't define pert. for K=6,8,9 if its core Mc>MCh */ 
       Lc=Lwdf(Mc,0.0,11);                       
       Rc=Rwdf(Mc,11);
     }
   }  
   b=0.002*max(1.0,2.5/M);
   c=0.006*max(1.0,2.5/M);
   q=log((*R)/Rc);
   s=(1.0+b*b*b)*pow(mi/b,3.0)/(1.0+pow(mi/b,3.0));
   if((*R)<=(Rc+acc))                                       /* taken form hrdiag.f l.651 and expanded by acc  */
     r=0.0;
   else  
     r=(1.0+c*c*c)*pow(mi/c,3.0)*pow(mi,0.1/q)/(1.0+pow(mi/c,3.0));
   (*L)=Lc*pow((*L)/Lc,s);
   (*R)=Rc*pow((*R)/Rc,r);  
 
   returnT=1;
 }
   
   
   
 if(RMAX==1 && (K==2 || K==3 || K==4 || K==5 || K==6)) {                 /* impose maximum radius on massive H-rich giants */                                                         
   
   lRmax=1.35799280*log10(2.71249545e-13*Mzams)-6.53597054e-03*Mzams+18.4428119;  /* Z-independent log10 maximum radius fit */ 
   Rmax=pow(10.0,lRmax);                                                          /* max radius [Rsun] (Amedeo fit) */

   lRmax100=1.35799280*log10(2.71249545e-13*100.0)-6.53597054e-03*100.0+18.4428119; 
   Rmax100=pow(10.0,lRmax100);

   if(Mzams>=18.0 && Mzams<100.0 && (*R)>Rmax)                           /* horizontal shift on H-R diagram */
     (*R)=Rmax;                                                          /* L is not changed so T decreases automatically */ 
   else if(Mzams>=100.0 && (*R)>Rmax100)                                 /* through get_T(L,R), e.g. for wind estimate */
     (*R)=Rmax100;                                  
 }       


 return returnT;
}





double tnf(double M, int K)
{ /* calculates the nuclear timescale, from Jarod code: star.f l.220-284 */
  /* uses tmax -estimate of time when Mcco reaches Mcmax on AGB: star.f l.193-218 */
 double tn,tmax,tau,Mt,Mx,Mcbagb,mc1,mc2,mcmax;
 double lambda,D,B,p,q,Ahp,Ahe,Ahhe;

 Mt=M;
 Mcbagb=Mcbagbf(M0);
 Mx=Mxf(M0);
 D=Df(M0);
 B=Bf(M0);
 p=pf(M0);
 q=qf(M0);
 Ahp=Ahpf(M0);
 Ahe=Ahef(); 
 Ahhe=Ahhef();
 lambda=min(0.9,0.3+0.001*pow(M0,5.0));      


 if(K<5)                              /* set mc1 as in star.f */
   mc1=Mcgbf2(M0,Lhe1f(M0));
 else {
   mc1=Mcbagb;                        /* mc1 -He core mas on EAGB and at start of TPAGB, star.f l.163 */
   if(mc1>=0.8 && mc1<=2.25)
     mc1=0.44*mc1+0.448;
 }


 tau=the1f(M0)+thef(M0);                           /* this is time not a timescale!, name from star.f */
 mc2=Mceagbf1(M0,tau);              
 mcmax=max(max(MCh,0.773*Mcbagb-0.35),1.05*mc2);   /* similar to my Mcsn */ 
 if(mcmax<=mc1) { 
   if(mcmax<=Mx) 
     tmax=tinf1ef(M0)-pow(mcmax,1.0-p)/((p-1.0)*Ahe*D);
   else
     tmax=tinf2ef(M0)-pow(mcmax,1.0-q)/((q-1.0)*Ahe*B);  
 }
 else {
   mcmax=(mcmax-lambda*mc1)/(1.0-lambda);
   if(mcmax<=Mx)
     tmax=tinf1tf(M0)-pow(mcmax,1.0-p)/((p-1.0)*Ahhe*D);
   else 
     tmax=tinf2tf(M0)-pow(mcmax,1.0-q)/((q-1.0)*Ahhe*B);
 }      
 tmax=max(tbagbf(M0),tmax);


 if(K<5 && fabs(Mt-Mcbagb)<acc)      
   tn=tbagbf(M0);
 else {
   if(Mt>Mcbagb || (Mt>=mc1 && K>4)) {
     if(K==6) 
       mc1=(Mt-lambda*mc1)/(1.0-lambda); 
     else 
       mc1=Mt;
     if(mc1<=Mx)
       tn=tinf1tf(M0)-pow(mc1,1.0-p)/((p-1.0)*Ahhe*D); 
     else
       tn=tinf2tf(M0)-pow(mc1,1.0-q)/((q-1.0)*Ahhe*B);  
   }   
   else {
     if(M0>M_FGB) {
       mc1=Mche1f(M0);
       if(Mt<=mc1)
         tn=the1f(M0);
       else
         tn=the1f(M0)+thef(M0)*((Mt-mc1)/(Mcbagb-mc1)); 
     } 
     else if(M0<=M_HeF) {
       if(K<8) {                                 /* for EAGB and TPAGB the same functions as for RG apply */
         mc1=Mcgbf2(M0,Lbgbf(M0));
         mc2=Mcgbf2(M0,Lhe1f(M0)); 
       }
       else {
         mc1=Mchsgbf2(M0,Lbgbf(M0));
         mc2=Mchsgbf2(M0,Lhe1f(M0));
       }
       if(Mt<=mc1)
         tn=tbgbf(M0);
       else if(Mt<=mc2) {
         if(Mt<=Mx)
           tn=tinf1f(M0)-pow(Mt,1.0-p)/((p-1.0)*Ahp*D);
         else
           tn=tinf2f(M0)-pow(Mt,1.0-q)/((q-1.0)*Ahp*B);                       
       }
       else
         tn=the1f(M0)+thef(M0)*((Mt-mc2)/(Mcbagb-mc2));
     }
     else {
       mc1=Mcbgbf1a(M0);
       mc2=Mche1f(M0);  
       if(Mt<=mc1)
         tn=tbgbf(M0);
       else if(Mt<=mc2)  
         tn=tbgbf(M0)+(the1f(M0)-tbgbf(M0))*((Mt-mc1)/(mc2-mc1));
       else
         tn=the1f(M0)+thef(M0)*((Mt-mc2)/(Mcbagb-mc2));  
     }   
   }    
 }
 tn=min(tn,tmax);                        

 return tn;
} 




/* ---------------------------------- MAIN SEQUENCE FORMULAE ---------------------------------------- */

double Lzamsf(double M)
{ /* for star M [M_sun] and metallicity ZZ return Lzams [L_sun] -- luminosity at ZAMS  */
  /* eq.(1,3)+table1: Tout et al. 1996, MNRAS 281,257; needs to know global variable ZZ */
 double alfaL,betaL,gammaL,deltaL,epsilonL,dzetaL,etaL;
 
 double a[]={0.3970417,8.527626,0.00025546,5.432889,5.563579,0.7886606,0.00586685};
 double b[]={-0.32913574,-24.41225973,-0.00123461,-8.62157806,-10.32345224,-2.90870942,-0.01704237};
 double c[]={0.34776688,56.43597107,-0.00023246,13.44202049,19.44322980,6.54713531,0.03872348};
 double d[]={0.37470851,37.06152575,0.00045519,14.51584135,18.97361347,4.05606657,0.02570041};
 double e[]={0.09011915,5.45624060,0.00016176,3.39793084,4.16903097,0.53287322,0.00383376};
 
 alfaL=a[0]+b[0]*log10(ZZ/ZZsun)+c[0]*pow(log10(ZZ/ZZsun),2.0)+d[0]*pow(log10(ZZ/ZZsun),3.0)+e[0]*pow(log10(ZZ/ZZsun),4.0);
 betaL=a[1]+b[1]*log10(ZZ/ZZsun)+c[1]*pow(log10(ZZ/ZZsun),2.0)+d[1]*pow(log10(ZZ/ZZsun),3.0)+e[1]*pow(log10(ZZ/ZZsun),4.0);
 gammaL=a[2]+b[2]*log10(ZZ/ZZsun)+c[2]*pow(log10(ZZ/ZZsun),2.0)+d[2]*pow(log10(ZZ/ZZsun),3.0)+e[2]*pow(log10(ZZ/ZZsun),4.0);
 deltaL=a[3]+b[3]*log10(ZZ/ZZsun)+c[3]*pow(log10(ZZ/ZZsun),2.0)+d[3]*pow(log10(ZZ/ZZsun),3.0)+e[3]*pow(log10(ZZ/ZZsun),4.0);
 epsilonL=a[4]+b[4]*log10(ZZ/ZZsun)+c[4]*pow(log10(ZZ/ZZsun),2.0)+d[4]*pow(log10(ZZ/ZZsun),3.0)+e[4]*pow(log10(ZZ/ZZsun),4.0);
 dzetaL=a[5]+b[5]*log10(ZZ/ZZsun)+c[5]*pow(log10(ZZ/ZZsun),2.0)+d[5]*pow(log10(ZZ/ZZsun),3.0)+e[5]*pow(log10(ZZ/ZZsun),4.0);
 etaL=a[6]+b[6]*log10(ZZ/ZZsun)+c[6]*pow(log10(ZZ/ZZsun),2.0)+d[6]*pow(log10(ZZ/ZZsun),3.0)+e[6]*pow(log10(ZZ/ZZsun),4.0);

 return (alfaL*pow(M,5.5)+betaL*pow(M,11.0))/(gammaL+pow(M,3.0)+deltaL*pow(M,5.0)+
         epsilonL*pow(M,7.0)+dzetaL*pow(M,8.0)+etaL*pow(M,9.5));
}         

          
double Rzamsf(double M)
{ /* for star M [M_sun] and metallicity ZZ return Rzams [R_sun] -- radius at ZAMS  */
  /* eq.(2,4)+table2: Tout et al. 1996, MNRAS 281,257; needs to know global variable ZZ */
 double tetaR,jotaR,kappaR,lambdaR,miR,niR,ksiR,omikronR,piR;
 
 double a[]={1.715359,6.597788,10.08855,1.012495,0.07490166,0.01077422,3.082234,17.84778,0.00022582};
 double b[]={0.62246212,-0.42450044,-7.11727086,0.32699690,0.02410413,0.0,0.94472050,-7.45345690,-0.00186899};
 double c[]={-0.92557761,-12.13339427,-31.67119479,-0.00923418,0.07233664,0.0,-2.15200882,-48.96066856,0.00388783};
 double d[]={-1.16996966,-10.73509484,-24.24848322,-0.03876858,0.03040467,0.0,-2.49219496,-40.05386135,0.00142402};
 double e[]={-0.30631491,-2.51487077,-5.33608972,-0.00412750,0.00197741,0.0,-0.63848738,-9.09331816,-0.00007671};
  
 tetaR=a[0]+b[0]*log10(ZZ/ZZsun)+c[0]*pow(log10(ZZ/ZZsun),2.0)+d[0]*pow(log10(ZZ/ZZsun),3.0)+e[0]*pow(log10(ZZ/ZZsun),4.0);
 jotaR=a[1]+b[1]*log10(ZZ/ZZsun)+c[1]*pow(log10(ZZ/ZZsun),2.0)+d[1]*pow(log10(ZZ/ZZsun),3.0)+e[1]*pow(log10(ZZ/ZZsun),4.0);
 kappaR=a[2]+b[2]*log10(ZZ/ZZsun)+c[2]*pow(log10(ZZ/ZZsun),2.0)+d[2]*pow(log10(ZZ/ZZsun),3.0)+e[2]*pow(log10(ZZ/ZZsun),4.0);
 lambdaR=a[3]+b[3]*log10(ZZ/ZZsun)+c[3]*pow(log10(ZZ/ZZsun),2.0)+d[3]*pow(log10(ZZ/ZZsun),3.0)+e[3]*pow(log10(ZZ/ZZsun),4.0);
 miR=a[4]+b[4]*log10(ZZ/ZZsun)+c[4]*pow(log10(ZZ/ZZsun),2.0)+d[4]*pow(log10(ZZ/ZZsun),3.0)+e[4]*pow(log10(ZZ/ZZsun),4.0);
 niR=a[5]+b[5]*log10(ZZ/ZZsun)+c[5]*pow(log10(ZZ/ZZsun),2.0)+d[5]*pow(log10(ZZ/ZZsun),3.0)+e[5]*pow(log10(ZZ/ZZsun),4.0);
 ksiR=a[6]+b[6]*log10(ZZ/ZZsun)+c[6]*pow(log10(ZZ/ZZsun),2.0)+d[6]*pow(log10(ZZ/ZZsun),3.0)+e[6]*pow(log10(ZZ/ZZsun),4.0);
 omikronR=a[7]+b[7]*log10(ZZ/ZZsun)+c[7]*pow(log10(ZZ/ZZsun),2.0)+d[7]*pow(log10(ZZ/ZZsun),3.0)+e[7]*pow(log10(ZZ/ZZsun),4.0);
 piR=a[8]+b[8]*log10(ZZ/ZZsun)+c[8]*pow(log10(ZZ/ZZsun),2.0)+d[8]*pow(log10(ZZ/ZZsun),3.0)+e[8]*pow(log10(ZZ/ZZsun),4.0);
          
 return (tetaR*pow(M,2.5)+jotaR*pow(M,6.5)+kappaR*pow(M,11.0)+lambdaR*pow(M,19.0)+miR*pow(M,19.5))/
        (niR+ksiR*pow(M,2.0)+omikronR*pow(M,8.5)+pow(M,18.5)+piR*pow(M,19.5));
}

double thookf(double M)
{ /* calculates time [10^6yrs] of Main Sequence hook of star M, eq.(7,5) */
 double mi,tbgb;
 
 tbgb=tbgbf(M);
 mi=max(0.5,1.0-0.01*max(aa6/pow(M,aa7),aa8+aa9/pow(M,aa10)));
 
 return mi*tbgb;
}
 
double tmsf(double M)
{ /* calculate lifetime [10^6yrs] of star M at Main Sequence, eq.(5,6) */
 double thook,tbgb,x;
 double dzeta=log10(ZZ/0.02);
 
 tbgb=tbgbf(M);
 thook=thookf(M);
 x=max(0.95,min(0.95-0.03*(dzeta+0.30103),0.99));
    
 return max(thook,x*tbgb);
} 

double Lmsf(double M, double t)
{ /* calculates Luminosity [Lsun] of star M during MS at time t [10^6yrs], eq.(11,12,14,15,16,18,19,19a,20) */ 
 double delL,alfaL,betaL,etaL,B,Mtmp,Ltmp;
 double tau1,tau2,epsilon;
 double Lzams,Ltms,thook,ttms,tau;                
 
 Lzams=Lzamsf(M);
 Ltms=Ltmsf(M); 
 thook=thookf(M);
 ttms=tmsf(M);
 tau=t/ttms;                  /* fractional timescale on MS */
 
 Mtmp=aa33; 
 B=min(aa34/pow(Mtmp,aa35),aa36/pow(Mtmp,aa37));
 if(M<=M_hook)
   delL=0.0;
 else if(M<aa33)
   delL=B*pow((M-M_hook)/(aa33-M_hook),0.4);  
 else
   delL=min(aa34/pow(M,aa35),aa36/pow(M,aa37));
 
 Mtmp=2.0;
 B=(aa45+aa46*pow(Mtmp,aa48))/(pow(Mtmp,0.4)+aa47*pow(Mtmp,1.9));
 if(M<0.5)
   alfaL=aa49;
 else if(M<0.7)
   alfaL=aa49+5.0*(0.3-aa49)*(M-0.5);
 else if(M<aa52)
   alfaL=0.3+(aa50-0.3)*(M-0.7)/(aa52-0.7);
 else if(M<aa53)
   alfaL=aa50+(aa51-aa50)*(M-aa52)/(aa53-aa52); 
 else if(M<2.0)
   alfaL=aa51+(B-aa51)*(M-aa53)/(2.0-aa53);  
 else
   alfaL=(aa45+aa46*pow(M,aa48))/(pow(M,0.4)+aa47*pow(M,1.9));   

 Mtmp=aa57;
 B=max(0.0,aa54-aa55*pow(Mtmp,aa56));
 betaL=max(0.0,aa54-aa55*pow(M,aa56));
 if(M>aa57 && betaL>0.0)
   betaL=max(0.0,B-10.0*(M-aa57)*B);
 
 if(ZZ<=0.0009)
   if(M<=1.0) 
     etaL=10.0;
   else if (M>=1.1) 
     etaL=20.0;
   else 
     etaL=inter_line(1.0,10.0,1.1,20.0,M); 
 else etaL=10.0;
 
 epsilon=0.01;
 tau1=min(1.0,t/thook);
 tau2=max(0.0,min(1.0,(t-(1.0-epsilon)*thook)/(epsilon*thook)));

 Ltmp=alfaL*tau+betaL*pow(tau,etaL)+(log10(Ltms/Lzams)-alfaL-betaL)*tau*tau-delL*(tau1*tau1-tau2*tau2);
 Ltmp=pow(10.0,Ltmp);
 Ltmp*=Lzams;

 return Ltmp;
}


double Rmsf(double M, double t) 
{ /* calculates Radius [Rsun] of star M during MS at time t [10^6yrs] */ 
  /* eq.(11,13,14,15,17,18,21,21a,22,22a,23,24) */
 double delR,alfaR,betaR,betaRp,gamma,B,C,X,Mtmp,Rtmp;
 double tau1,tau2,epsilon;
 double Rzams,Rtms,thook,ttms,tau;

 Rzams=Rzamsf(M);
 Rtms=Rtmsf(M);
 thook=thookf(M);
 ttms=tmsf(M);
 tau=t/ttms;                  /* fractional timescale on MS */
    
 Mtmp=2.0;   
 B=(aa38+aa39*pow(Mtmp,3.5))/(aa40*Mtmp*Mtmp*Mtmp+pow(Mtmp,aa41))-1.0;   
 if(M<=M_hook)
   delR=0.0;
 else if(M<=aa42)
   delR=aa43*pow((M-M_hook)/(aa42-M_hook),0.5);      
 else if(M<2.0)
   delR=aa43+(B-aa43)*pow((M-aa42)/(2.0-aa42),aa44);
 else
   delR=(aa38+aa39*pow(M,3.5))/(aa40*M*M*M+pow(M,aa41))-1.0;

 Mtmp=aa66;
 B=aa58*pow(Mtmp,aa60)/(aa59+pow(Mtmp,aa61));
 Mtmp=aa67;
 C=aa58*pow(Mtmp,aa60)/(aa59+pow(Mtmp,aa61));  
 if(M<0.5)
   alfaR=aa62;
 else if(M<0.65)
   alfaR=aa62+(aa63-aa62)*(M-0.5)/0.15;
 else if(M<aa68)
   alfaR=aa63+(aa64-aa63)*(M-0.65)/(aa68-0.65);
 else if(M<aa66)
   alfaR=aa64+(B-aa64)*(M-aa68)/(aa66-aa68);
 else if(M<=aa67)
   alfaR=aa58*pow(M,aa60)/(aa59+pow(M,aa61));  
 else
   alfaR=C+aa65*(M-aa67);     
  
 Mtmp=2.0; 
 B=aa69*pow(Mtmp,3.5)/(aa70+pow(Mtmp,aa71));
 Mtmp=16.0;
 C=aa69*pow(Mtmp,3.5)/(aa70+pow(Mtmp,aa71));
 if(M<=1.0)
   betaRp=1.06; 
 else if(M<aa74)
   betaRp=1.06+(aa72-1.06)*(M-1.0)/(aa74-1.06);
 else if(M<2.0)
   betaRp=aa72+(B-aa72)*(M-aa74)/(2.0-aa74); 
 else if(M<=16.0)
   betaRp=aa69*pow(M,3.5)/(aa70+pow(M,aa71));
 else
   betaRp=C+aa73*(M-16.0);     
  
 betaR=betaRp-1.0;   
   
   
 Mtmp=1.0;  
 B=aa76+aa77*pow(Mtmp-aa78,aa79);
 if(aa75>1.0)
   C=aa80;  
 else
   C=B;
 if(M<=1.0)
   gamma=aa76+aa77*pow(M-aa78,aa79);
 else if(M<=aa75)
   gamma=B+(aa80-B)*pow((M-1.0)/(aa75-1.0),aa81);
 else if(M<=(aa75+0.1))                      /* corrected: instead of 1.0 I have put 0.1 as it shall be in paper0 */                      
   gamma=C-10.0*(M-aa75)*C;
 else
   gamma=0.0;        
   
 epsilon=0.01;
 tau1=min(1.0,t/thook);
 tau2=max(0.0,min(1.0,(t-(1.0-epsilon)*thook)/(epsilon*thook)));
     
 if(tau<1.0e-04)               /* to avoid (tau,10.0) and (tau,40.0) being exceedingly small: numerics */
   Rtmp=alfaR*tau+(log10(Rtms/Rzams)-alfaR-betaR-gamma)*tau*tau*tau-delR*(tau1*tau1*tau1-tau2*tau2*tau2);    
 else    
   Rtmp=alfaR*tau+betaR*pow(tau,10.0)+gamma*pow(tau,40.0)+(log10(Rtms/Rzams)-alfaR-betaR-gamma)
        *tau*tau*tau-delR*(tau1*tau1*tau1-tau2*tau2*tau2);  
 Rtmp=pow(10.0,Rtmp);
 Rtmp*=Rzams;  
  
 if(M<0.1) {
   X=0.76-3.0*ZZ;               /* hydrogen abundance from Tout et al. 1996, MNRAS 281, 257 (p.257) */
   Rtmp=max(Rtmp,0.0258*pow(1.0+X,5.0/3.0)*pow(M,-1.0/3.0));  
 }

 return Rtmp;
}


double Ltmsf(double M)
{ /* calculate luminosity [Lsun] of star M at Terminal Main Sequence, eq.(8) */ 

 return (aa11*M*M*M+aa12*M*M*M*M+aa13*pow(M,aa16+1.8))/(aa14+aa15*M*M*M*M*M+pow(M,aa16)); 
}


double Rtmsf(double M)
{ /* calculate radius [Rsun] of star M at Terminal Main Sequence, eq.(9,9a) */
 double Rtmp,y1,y2;
 double Mstar=aa17+0.1;
 double c1=-8.672073e-2;

 if(M<=aa17)
   Rtmp=(aa18+aa19*pow(M,aa21))/(aa20+pow(M,aa22));
 else if(M>=Mstar)
   Rtmp=(c1*M*M*M+aa23*pow(M,aa26)+aa24*pow(M,aa26+1.5))/(aa25+pow(M,5.0));
 else {
   y1=(aa18+aa19*pow(aa17,aa21))/(aa20+pow(aa17,aa22));
   y2=(c1*Mstar*Mstar*Mstar+aa23*pow(Mstar,aa26)+aa24*pow(Mstar,aa26+1.5))/(aa25+pow(Mstar,5.0));
   Rtmp=inter_line(aa17,y1,Mstar,y2,M);
 }  
   
 if(M<0.5) 
   Rtmp=max(Rtmp,1.5*Rzamsf(M));  
  
 return Rtmp;
}



/* -------------------------------- HARTZPRUNG GAP PHASE --------------------------------------- */

double thgf(double M)
{ /* calculate lifetime [10^6yrs] of star M at Hartzprung Gap, eq.(5,6) */
 
 return tbgbf(M)-tmsf(M);
}


double Mchgf(double M, double t)
{ /* calculate core mass [Msun] of star M at time t during Hartzprung Gap, eq.(29,30) */ 
 double ttms,tbgb,tau,ro,Mcehg;
 
 ttms=tmsf(M);
 tbgb=tbgbf(M);
 tau=(t-ttms)/(tbgb-ttms);
 ro=(1.586+pow(M,5.25))/(2.434+1.02*pow(M,5.25));
 Mcehg=Mcehgf(M);
 
 return ((1.0-tau)*ro+tau)*Mcehg;
}


double Lhgf(double M, double t)
{ /* calculate luminosity [Lsun] of star M during Hartzprung Gap at time t [10^6yrs], eq.(25,26) */
 double tau,ttms,tbgb;
 double Ltms,Lehg;
 
 ttms=tmsf(M);
 tbgb=tbgbf(M);
 Ltms=Ltmsf(M);
 Lehg=Lehgf(M);
 
 tau=(t-ttms)/(tbgb-ttms);

 return Ltms*pow(Lehg/Ltms,tau);
}


double Rhgf(double M, double Mtmp, double t)
{ /* calculate Radius [Rsun] of star M during Hartzprung Gap at time t [10^6yrs], eq.(25,27) */
 double tau,ttms,tbgb;
 double Rtms,Rehg;
   
 ttms=tmsf(Mtmp);
 tbgb=tbgbf(Mtmp);
 tau=(t-ttms)/(tbgb-ttms);
 Rtms=Rtmsf(Mtmp);
 Rehg=Rehgf(M,Mtmp);

 return Rtms*pow(Rehg/Rtms,tau);
}
           

double Mcehgf(double M)
{ /* calculate core mass [Msun] of star M at End Hartzprung Gap, eq.(28) */
  /* inequalities are taken from hrdiag.f l.143 and that differ from eq.(28) */

 if(M<=M_HeF) 
   return Mcgbf2(M,Lbgbf(M));
 else if(M<=M_FGB)
   return Mcbgbf1a(M);  
 else
   return Mche1f(M);  
} 


double Lehgf(double M)
{ /* calculate luminosity [Lsun] of star M at End Hartzprung Gap, sec.5.1 */

 if(M<M_FGB)
   return Lbgbf(M);
 else
   return Lhe1f(M);
} 


double Rehgf(double M, double Mtmp)
{ /* calculate radius [Rsun] of star M at End Hartzprung Gap as in hrdiag.f l.177-190 */ 
  /* in paper0 in sec.5.1 when I use my functions for Rgb(),Rhe1() I won't get the wanted result */ 
 double rx,ry,rmin,mi,tau2,Lhe1,Lbgb;
 
 if(Mtmp<=M_FGB) {
   Lbgb=Lbgbf(Mtmp);
   rx=Rgbf(M,Lbgb);
 }  
 else {
   Lhe1=Lhe1f(Mtmp);
   rmin=Rmhef(Mtmp);
   ry=Ragbf(M,Lhe1);   
   rx=min(rmin,ry);
   if(Mtmp<=12.0) {
     mi=log(Mtmp/12.0)/log(M_FGB/12.0);
     rx=Rgbf(M,Lhe1);
     rx=rmin*pow(rx/rmin,mi);
   }
   tau2=taublf(Mtmp);
   if(fabs(tau2)<acc)
     rx=ry;
 }

 return rx;
}   



/* ------------------------------------ GIANT BRANCHES ---------------------------------------- */
/* -------------------------------------------------------------------------------------------- */


/* -------------------------------- BASE OF RED GIANT BRANCH ---------------------------------- */

double tbgbf(double M)
{ /* calculate time [10^6yrs] at which star of mass M reaches Base of Giant Branch, eq.(4) */
 
 return (aa1+aa2*pow(M,4.0)+aa3*pow(M,5.5)+pow(M,7.0))/(aa4*M*M+aa5*pow(M,7.0));
} 


double Mcbgbf1a(double M)
{ /* calculate REAL core mass [Msun] for star of mass M at Base of Giant Branch, eq.(39,44) */
  /* formula (44) works only for M>=M_HeF, so for M<M_HeF I use (39) with t=tbgb */
 double p,q,Ahp,B,D,C,c1,c2,Ltmp,Mtmp;
 double tbgb,tinf1,tinf2,tx,Mcbagb;

 if(M<M_HeF) {
   tbgb=tbgbf(M);
   tx=txf(M);
   Ahp=Ahpf(M);
   if(tbgb<=tx) {
     p=pf(M);
     D=Df(M);
     tinf1=tinf1f(M);
     return pow((p-1)*Ahp*D*(tinf1-tbgb),1.0/(1.0-p));
   }
   else {
     q=qf(M);           
     B=Bf(M);
     tinf2=tinf2f(M);   
     return pow((q-1)*Ahp*B*(tinf2-tbgb),1.0/(1.0-q));
   }
 }
 else {
   c1=9.20925e-5;
   c2=5.402216;
   Ltmp=Lbgbf(M_HeF);
   Mtmp=Mcgbf2(M_HeF,Ltmp);
   C=pow(Mtmp,4.0)-c1*pow(M_HeF,c2);
   Mcbagb=Mcbagbf(M);
   return min(0.95*Mcbagb,pow(C+c1*pow(M,c2),0.25));
 }
}


double Mcbgbf1b(double M)
{ /* calculate DUMMY core mass [Msun] for star of mass M at Base of Giant Branch, eq.(39) */
  /* used to get Luminosity during Red Giant Branch */
 double p,q,Ahp,B,D,tbgb,tinf1,tinf2,tx;

 tbgb=tbgbf(M);
 tx=txf(M);
 Ahp=Ahpf(M);
 if(tbgb<=tx) {
   p=pf(M);         
   D=Df(M);
   tinf1=tinf1f(M);
   return pow((p-1)*Ahp*D*(tinf1-tbgb),1.0/(1.0-p));
 }
 else {
   q=qf(M); 
   B=Bf(M);
   tinf2=tinf2f(M);
   return pow((q-1)*Ahp*B*(tinf2-tbgb),1.0/(1.0-q));
 }
}


double Lbgbf(double M)
{ /* calculate Luminosity [Lsun] for star of mass M at Base of Giant Branch, eq.(10) */
 double c2=9.301992,c3=4.637345;
 
 return (aa27*pow(M,aa31)+aa28*pow(M,c2))/(aa29+aa30*pow(M,c3)+pow(M,aa32)); 
} 


/* ------------------------------------ RED GIANT BRANCH ------------------------------------- */

double Ahpf(double M)
{ /* calculates Ahp parameter, sec.5.2 */
 double tmp;

 tmp=max(-4.8,min(-5.7+0.8*M,-4.1+0.14*M));

 return pow(10.0,tmp);
}


double Bf(double M)
{ /* calculates B parameter, sec.5.2 */

 return max(3.0e+4,500.0+1.75e+4*pow(M,0.6));
} 


double Df(double M)
{ /* calculates D parameter, sec.5.2, but instead of M_HeF 2.0 is used as in code! */
 double D0,D1,Mtmp,tmp;
 double dzeta=log10(ZZ/0.02);

 D0=5.37+0.135*dzeta;
 if(M<=2.0) tmp=D0;
 else if(M>=2.5) tmp=max(-1.0,max(0.975*D0-0.18*M,0.5*D0-0.06*M));
 else {
   Mtmp=2.5;
   D1=max(-1.0,max(0.975*D0-0.18*Mtmp,0.5*D0-0.06*Mtmp));
   tmp=inter_line(2.0,D0,2.5,D1,M);
 }

 return pow(10.0,tmp);
}


double pf(double M)
{ /* returns parameter p, sec.5.2, but instead of M_HeF 2.0 is used as in code! */
 double p;

 if(M<=2.0) p=6.0;
 else if(M>=2.5) p=5.0;
 else p=inter_line(2.0,6.0,2.5,5.0,M);

 return p;
}


double qf(double M)
{ /* returns parameter q, sec.5.2, but instead of M_HeF 2.0 is used as in code! */
 double q;

 if(M<=2.0) q=3.0;
 else if(M>=2.5) q=2.0;
 else q=inter_line(2.0,3.0,2.5,2.0,M);

 return q;
}      


double tinf1f(double M)
{ /* calculates integration constant tinf1[Myr], eq. 40 */    
 double tbgb,Lbgb,p,Ahp,D;
 
 tbgb=tbgbf(M);
 Lbgb=Lbgbf(M);
 p=pf(M);
 Ahp=Ahpf(M);
 D=Df(M);
 
 return tbgb+pow(D/Lbgb,(p-1.0)/p)/((p-1.0)*Ahp*D);
}


double tinf2f(double M)
{ /* calculates integration constant tinf2[Myr], eq. 42 */    
 double tx,Lx,q,Ahp,B;
 
 tx=txf(M);
 Lx=Lxf(M);
 q=qf(M);
 Ahp=Ahpf(M);
 B=Bf(M);
 
 return tx+pow(B/Lx,(q-1.0)/q)/((q-1.0)*Ahp*B);
}


double txf(double M)
{ /* calculates integration constant tx[Myr], eq. 42 */
 double tinf1,tbgb,Lbgb,Lx,p;
 
 tinf1=tinf1f(M);
 tbgb=tbgbf(M);
 Lbgb=Lbgbf(M);
 Lx=Lxf(M);
 p=pf(M);
 
 return tinf1-pow(Lbgb/Lx,(p-1.0)/p)*(tinf1-tbgb);
} 


double Mxf(double M)
{ /* calculates core mass Mx[Msun] for given M, eq.(38) */
 double p,q,B,D;
 
 p=pf(M);
 q=qf(M);
 B=Bf(M);
 D=Df(M);
 
 return pow(B/D,1.0/(p-q));
} 


double Lxf(double M)
{ /* calculates luminosity Lx[Lsun] for given M, from eq.(37) */
  /* B*pow(Mx,q)=D*pow(Mx,p) so it doesn't matter which one I return */
 double D,p,Mx;

 D=Df(M);
 p=pf(M);
 Mx=Mxf(M);
     
 return D*pow(Mx,p);
}


double Mcgbf1a(double M, double t)
{ /* calculates REAL core mass Mcgb[Msun] for red giant of M at time t[Myr], eq.(39,44,45) */
 double p,q,Ahp,B,D,Mcbgb,Mche1;
 double tbgb,the1,tinf1,tinf2,tx,tau;

 if(M<M_HeF) {
   tx=txf(M);
   Ahp=Ahpf(M);
   if(t<=tx) {
     p=pf(M);
     D=Df(M);
     tinf1=tinf1f(M);
     return pow((p-1)*Ahp*D*(tinf1-t),1.0/(1.0-p));
   }
   else {
     q=qf(M);
     B=Bf(M);
     tinf2=tinf2f(M);
     return pow((q-1)*Ahp*B*(tinf2-t),1.0/(1.0-q));  
   }
 }
 else {
   tbgb=tbgbf(M);
   the1=the1f(M);
   Mcbgb=Mcbgbf1a(M);
   Mche1=Mche1f(M);     
   tau=(t-tbgb)/(the1-tbgb);     
   return Mcbgb+(Mche1-Mcbgb)*tau;
 }  
}


double Mcgbf1b(double M, double t)
{ /* calculates DUMMY core mass Mcgb[Msun] for red giant of M at time t[Myr], eq.(39) */
 double p,q,Ahp,B,D,tinf1,tinf2,tx;

 tx=txf(M);
 Ahp=Ahpf(M);
 if(t<=tx) {
   p=pf(M);
   D=Df(M);
   tinf1=tinf1f(M);
   return pow((p-1)*Ahp*D*(tinf1-t),1.0/(1.0-p));
 }
 else {
   q=qf(M);
   B=Bf(M);
   tinf2=tinf2f(M);  
   return pow((q-1)*Ahp*B*(tinf2-t),1.0/(1.0-q));  
 }
}


double Mcgbf2(double M, double L)
{ /* calculates core mass Mcgb[Msun] for red giant of M and luminosity L[Lsun], inverted eq.(37) */
 double Lx,p,q,B,D;
 
 Lx=Lxf(M);
 if(L<=Lx) {
   p=pf(M);
   D=Df(M);
   return pow(L/D,1.0/p);
 }
 else {
   q=qf(M);
   B=Bf(M);
   return pow(L/B,1.0/q);
 }
}   


double Lgbf1(double M, double t)
{ /* calculates luminosity[Lsun], using DUMMY Mc, for red giant of M at time t[Myr], eq.(37) */
 double p,q,B,D,Mcgb; 
 
 p=pf(M);
 q=qf(M);
 B=Bf(M);
 D=Df(M);
 Mcgb=Mcgbf1b(M,t);

 return min(B*pow(Mcgb,q),D*pow(Mcgb,p));
}


double Lgbf2(double M, double Mc)
{ /* calculates luminosity[Lsun], for DUMMY Mc, for red giant of M, eq.(37) */
 double p,q,B,D;
 
 p=pf(M);
 q=qf(M);
 B=Bf(M);
 D=Df(M);
      
 return min(B*pow(Mc,q),D*pow(Mc,p));
}

       
double Rgbf(double M, double L)
{ /* calculate Radius [Rsun] of either red or asymptotic giant of M and L, eq.(46) */
 double A;

 A=min(bb4*pow(M,-bb5),bb6*pow(M,-bb7));
 
 return A*(pow(L,bb1)+bb2*pow(L,bb3));
}



/* -------------------------------- CORE HELIUM BURNING PHASE ------------------------------------ */
/* ----------------------------------------------------------------------------------------------- */

double Rhe1f(double M) 
{ /* calculate Radius [Rsun] of core helium ignition as it is in 'hrdiag.f' incorporating Mass Loss */
  /* for LM and IM stars I use Rhe1=rmin from the hrdiag.f (l.265,l.325) for part des. CHeB phase */
  /* for HM stars I use Rhe1=rx from the hrdiag.f(l.190) -- part of it describing HG phase */
  /* this is guessed from comparison of hrdiag.f(l.190) part and eqs.(27,R_EHG) */
  /* DIRECTLY USES M0 -- FOR MASS LOSS!!! */
 double rx,ry,rmin,mi,tau2;
 double Mche1,Lhe1,lx;
 
 Mche1=Mche1f(M0);                 /* w hrdiag.f for t=the1 for which we calculate Rhe1, we have: */
                                   /* tau=0.0 => mc=mcx=Mche1f(M0) */
 if(M0<=M_HeF) {
   rx=Rzahbf(M,Mche1);
   lx=Lzahbf(M0,Mche1);
   rmin=Rmhefa(M,lx);
   if(rmin>=rx)
     rmin=rx;
 }    
 else if(M0>M_FGB) {
   rmin=Rmhef(M0);
   Lhe1=Lhe1f(M0);
   ry=Ragbf(M,Lhe1);   
   rx=min(rmin,ry);
   if(M0<=12.0) {
     mi=log(M0/12.0)/log(M_FGB/12.0);
     rx=Rgbf(M,Lhe1);
     rx=rmin*pow(rx/rmin,mi);
   }
   tau2=taublf(M0);
   if(fabs(tau2)<acc)
     rx=ry;
   rmin=rx;    
 }
 else {
   lx=Lzahbf(M0,Mche1);
   rx=Rgbf(M,lx);
   rmin=Rmhef(M);
   if(rmin>=rx) 
     rmin=rx;
 }    

 return rmin;
}


double Mche1f(double M)
{ /* calculate core mass Mche1[Msun] at core helium ignition for star of M, sec.5.3 - just after eq.(65) */
  /* formula (44) transformed in sec.5.3 from BGB to He1 works only for M>=M_HeF, */
  /* so for M<M_HeF I use (39) with t=the1 */
 double p,q,Ahp,B,D,c1,c2,C,Ltmp,Mtmp;
 double Mche1,Mcbagb,the1,tx,tinf1,tinf2; 

 if(M<M_HeF) {
   the1=the1f(M);
   tx=txf(M);
   Ahp=Ahpf(M);
   if(the1<=tx) {
     p=pf(M);
     D=Df(M);
     tinf1=tinf1f(M);
     Mche1=pow((p-1)*Ahp*D*(tinf1-the1),1.0/(1.0-p));
   }
   else {
     q=qf(M);
     B=Bf(M);
     tinf2=tinf2f(M);
     Mche1=pow((q-1)*Ahp*B*(tinf2-the1),1.0/(1.0-q));  
   }
 }
 else {
   c1=9.20925e-5;
   c2=5.402216;
   Ltmp=Lhe1f(M_HeF);
   Mtmp=Mcgbf2(M_HeF,Ltmp);
   C=pow(Mtmp,4.0)-c1*pow(M_HeF,c2);
   Mcbagb=Mcbagbf(M);
   Mche1=min(0.95*Mcbagb,pow(C+c1*pow(M,c2),0.25));
 }  

 return Mche1;
} 


double Lhe1f(double M)
{ /* calculate Luminosity [Lsun] of core helium ignition, eq.(49) */
 double alfa1,Ltmp,Mtmp;
 
 if(M<M_HeF) {
   Mtmp=M_HeF;
   Ltmp=(bb11+bb12*pow(Mtmp,3.8))/(bb13+Mtmp*Mtmp);
   alfa1=(bb9*pow(M_HeF,bb10)-Ltmp)/Ltmp;
   return bb9*pow(M,bb10)/(1.0+alfa1*exp(15.0*(M-M_HeF)));
 }
 else
   return (bb11+bb12*pow(M,3.8))/(bb13+M*M);  
}


double the1f(double M)
{ /* returns time of helium ignition the1[Myr] for star of mass M, eq.(4,43) */
  /* if star has M>M_FGB, then I use the1=tbgb from eq.(4), and this deduced from text */
  /* orginal equation for the1 is included in eq.(43), and gives too big times for M>M_FGB */ 
 double Lhe1,Lx,tinf1,tinf2,tbgb;
 double p,q,Ahp,B,D;
 
 if(M>M_FGB) {
   tbgb=tbgbf(M);
   return tbgb;
 }  
 else {
   Lhe1=Lhe1f(M);
   Lx=Lxf(M);
   Ahp=Ahpf(M);
   if(Lhe1<=Lx) {
     p=pf(M);
     D=Df(M);
     tinf1=tinf1f(M);
     return tinf1-(1.0/((p-1.0)*Ahp*D))*pow(D/Lhe1,(p-1.0)/p);
   }
   else {
     q=qf(M);
     B=Bf(M);
     tinf2=tinf2f(M);
     return tinf2-(1.0/((q-1.0)*Ahp*B))*pow(B/Lhe1,(q-1.0)/q);  
   }
 } 
}


double thef(double M)
{ /* calculate TOTAL lifetime [10^6yrs] of star M at Core Helium Burning, eq.(57) */
  /* from He ignition to the AGB phase */
 double alfa4,Mc,mi,the,thsms,tmp;

 if(M<M_HeF) {
   Mc=Mche1f(M);
   mi=(M-Mc)/(M_HeF-Mc);
   mi=max(mi,1.0e-12);                      /* added: as in Jarrod code */
   thsms=thsmsf(Mc);
   tmp=tbgbf(M_HeF)*(bb41*pow(M_HeF,bb42)+bb43*pow(M_HeF,5.0))/(bb44+pow(M_HeF,5.0));
   alfa4=(tmp-bb39)/bb39;
   the=(bb39+(thsms-bb39)*pow(1.0-mi,bb40))*(1.0+alfa4*exp(15.0*(M-M_HeF)));
 }
 else 
   the=tbgbf(M)*(bb41*pow(M,bb42)+bb43*pow(M,5.0))/(bb44+pow(M,5.0));

 return the; 
}


double Mchef(double M, double t)
{ /* calculates core mass Mche[Msun] during core helium burning for star of M, eq.(67) */
 double tau,the1,the,Mche1,Mcbagb;
 
 the1=the1f(M);
 the=thef(M);
 tau=(t-the1)/the;
 Mche1=Mche1f(M);
 Mcbagb=Mcbagbf(M);

 return (1.0-tau)*Mche1+tau*Mcbagb;
} 


double Lminhef(double M)
{ /* calculates minimum luminosity[Lsun] for IM star of M during CHeB, eq.(51) */
  /* for LM and HM stars Lminhe should not be used */
 double c;

 c=bb17/pow(M_FGB,0.1)+(bb16*bb17-bb14)/pow(M_FGB,bb15+0.1);
 
 return Lhe1f(M)*(bb14+c*pow(M,bb15+0.1))/(bb16+pow(M,bb15));
} 


double Lzahbf(double M, double Mc)
{ /* calculates Zero Age Horizontal Branch luminosity[Lsun] for LM and IM star of M and Mc */
  /* Mc is not needed for IM stars, eq.(51,52,53), and as deduced from star.f l.125-128 */
  /* Lzahb is the minimum luminosity for CHeB phase and is not adequate for HM stars */
  /* I use here 1.647903 (as in Jarrod code) and not 1.6479 as in paper -- makes difference */ 
 double Lzahb,Lzhs,Lminhe,alfa2,mi;
 
 if(M<=M_HeF) {
   Lzhs=Lzhsf(Mc);
   Lminhe=Lminhef(M_HeF);
   alfa2=(bb18+Lzhs-Lminhe)/(Lminhe-Lzhs);
   mi=(M-Mc)/(M_HeF-Mc);
   Lzahb=Lzhs+(1.0+bb20)*bb18*pow(mi,bb19)/((1.0+bb20*pow(mi,1.647903))*(1.0+alfa2*exp(15.0*(M-M_HeF))));
 }
 else if(M<=(M_FGB+acc)) 
   Lzahb=Lminhef(M);
 else
   fprintf(fp0,"error: function Lzahbf(M,Mc) should not be called for HM stars\n");    

 return Lzahb;
}


double Rzahbf(double M, double Mc)
{ /* calculates Zero Age Horizontal Branch radius[Lsun] for star of M, eq.(54) */
 double Rzhs,Lzahb,mi,f;
 
 mi=(M-Mc)/(M_HeF-Mc);
 f=(1.0+bb21)*pow(mi,bb22)/(1.0+bb21*pow(mi,bb23));
 Rzhs=Rzhsf(Mc);
 Lzahb=Lzahbf(M,Mc);
 
 return (1.0-f)*Rzhs+f*Rgbf(M,Lzahb);
}


double Rmhef(double M)
{ /* calculates minimum radius[Rsun] of star of M on Blue Loop (during CHeB), eq.(55) */
  /* this is for IM and HM stars */
 
 return (bb24*M+pow(bb25*M,bb26)*pow(M,bb28))/(bb27+pow(M,bb28));
}


double Rmhefa(double M, double Lzahb)
{ /* calculates minimum radius[Rsun] of star of M on Blue Loop (CHeB), eq.(55a): just after eq.(55) */
  /* using Lzahb. this is for LM stars. mi in eq.(55a) is wrong in paper0 */
  /* DIRECTLY USES M0 -- FOR MASS LOSS!!! */
 double Rgb1,Rgb2,Rminhe,Mc1,mi;
   
 Rgb1=Rgbf(M,Lzahb);
 Mc1=Mcgbf2(M_HeF,Lbgbf(M_HeF));
 Rgb2=Rgbf(M_HeF,Lzahbf(M_HeF,Mc1));
 Rminhe=Rmhef(M_HeF);
 mi=M0/M_HeF;
        
 return Rgb1*pow(Rminhe/Rgb2,mi);
}


double fblf(double M)
{ /* calculates variable fbl needed in eq.(58) for taubl lifetime at Blue Phase, sec.5.3 -after eq.(58) */
 double Rmhe,Ragb,Lhe1;
 
 Rmhe=Rmhef(M);
 Lhe1=Lhe1f(M);
 Ragb=Ragbf(M,Lhe1);

 return pow(M,bb48)*pow(1.0-Rmhe/Ragb,bb49);
}


double taublf(double M)
{ /* calculate lifetime [dimenshionless] of star M at Blue Phase relative to the, */
  /* where the is lifetime [10^6yrs] of star M at Core Helium Burning */ 
  /* not from eq.(58), but from modified eq.(58) as described in  prob6.ans and also using */
  /* exponent 0.4805428 as in zdata.f l.207 and in zcnsts.f l.326 instead of 0.414 as in eq.(58) */
 double tbl,b,c;

 if(M<M_HeF) 
   tbl=1.0;
 else if(M<=M_FGB) {
   c=0.4805428;
   b=log10(M/M_FGB)/log10(M_HeF/M_FGB);
   b=max(b,1.0e-12);                     /* Jarrod's Safe Side */ 
   tbl=bb45*pow(M/M_FGB,c)+(1.0-bb45*pow(M_HeF/M_FGB,c))*pow(b,bb46);
 }
 else 
   tbl=(1.0-bb47)*fblf(M)/fblf(M_FGB);

 if(tbl<0.0)
   tbl=0.0;
 else if(tbl>1.0)
   tbl=1.0;  

 return tbl;
}


int CHeBf(double M, double t, double *R, double *L, double *Mc)
{ /* this function follows closely hrdiag.f */
  /* DIRECTLY USES M0 -- FOR MASS LOSS!!! */
 double the1,the,tau,taul,tauh,tau2,tloop;
 double Lhe1,Lzahb,lx,ly,rx,ry,rmin,ksi,mi,Mche1;

 the1=the1f(M0);  
 the=thef(M0);
 tau=(t-the1)/the;
 
 Mche1=Mche1f(M0);                     /* starting mass: in hrdiag.f it is mcx */
 (*Mc)=Mchef(M0,t);                    /* current mass: in hrdiag.f it is mc */

 if(Mche1>M0 || Mche1>M || (*Mc)>M0 || (*Mc)>M)    /* for star to be here in Core he burning: */
   return 1;                                       /* core would need to be larger than a star */
                                                   /* abort: pass the star on */
 if(M0<=M_HeF) {
   Lhe1=Lhe1f(M0);
   Mche1=Mcgbf2(M0,Lhe1);              
   lx=Lzahbf(M0,Mche1);
   ly=Lbagbf(M0); 
   rx=Rzahbf(M,*Mc);
   rmin=Rmhefa(M,lx);
   ksi=min(max(0.4,rmin/rx),2.5);
   ry=Ragbf(M,ly);
   if(rmin<rx)
     taul=pow(log(rx/rmin),1.0/3.0);
   else {
     rmin=rx;
     taul=0.0;
   }
   tauh=pow(log(ry/rmin),1.0/3.0);  
   tau2=taul*(tau-1.0)+tauh*tau;  
   (*R)=rmin*exp(pow(fabs(tau2),3.0));  
   (*L)=lx*pow(ly/lx,pow(tau,ksi));  
 }    
 else if(M0>M_FGB) {
   tau2=taublf(M0);       
   tloop=the1+tau2*the;
   rmin=Rmhef(M0);  
   Lhe1=Lhe1f(M0);   
   rx=Ragbf(M,Lhe1);
   rmin=min(rmin,rx);
   if(M0<=12.0) {
     mi=log(M0/12.0)/log(M_FGB/12.0);
     rx=Rgbf(M,Lhe1);
     rx=rmin*pow(rx/rmin,mi);
   }      
   else 
     rx=rmin;
   ksi=min(max(0.4,rmin/rx),2.5);
   (*L)=Lhe1*pow(Lbagbf(M0)/Lhe1,pow(tau,ksi));
   if(t<tloop) {
     ly=Lhe1*pow(Lbagbf(M0)/Lhe1,pow(tau2,ksi));
     ry=Ragbf(M,ly);
     taul=0.0;
     if(rmin!=rx) 
       taul=pow(log(rx/rmin),1.0/3.0);
     tauh=pow(log(ry/rmin),1.0/3.0);
     tau=(t-the1)/(tau2*the);
     tau2=taul*(tau-1.0)+tauh*tau;
     (*R)=rmin*exp(pow(fabs(tau2),3.0));
   }
   else
     (*R)=Ragbf(M,*L);
 }
 else {
   tau2=1.0-taublf(M0);
   tloop=the1+tau2*the;
   Lzahb=Lzahbf(M0,*Mc);
   if(t<tloop) {
     tau=(tloop-t)/(tau2*the);
     (*L)=Lzahb*pow(Lhe1f(M0)/Lzahb,pow(tau,3.0));
     (*R)=Rgbf(M,*L);
   }
   else {
     lx=Lzahb;
     ly=Lbagbf(M0);
     rx=Rgbf(M,lx);
     rmin=Rmhef(M);
     ksi=min(max(0.4,rmin/rx),2.5);
     ry=Ragbf(M,ly);
     if(rmin<rx) 
       taul=pow(log(rx/rmin),1.0/3.0);
     else {
       rmin=rx;
       taul=0.0;
     }
     tauh=pow(log(ry/rmin),1.0/3.0);   
     tau=(t-tloop)/(the-(tloop-the1));
     tau2=taul*(tau-1.0)+tauh*tau;
     (*R)=rmin*exp(pow(fabs(tau2),3.0));
     (*L)=lx*pow(ly/lx,pow(tau,ksi));
   }
 }    

 return 0;
}


                      

/* ----------------------------------- ASYMPTOTIC GIANT BRANCH --------------------------------------- */
/* --------------------------------------------------------------------------------------------------- */


/* ----------------------------------------- BASE OF AGB --------------------------------------------- */


double tbagbf(double M)
{ /* calculate time[10^6yrs] when star of M reaches Base of AGB, sec.5.4 -- before eq.(68) */

 return the1f(M)+thef(M);
}


double Mcbagbf(double M)
{ /* calculates Mcbagb[Msun] -- He core mass at the Base of AGB, eq.(66) */

 return pow(bb36*pow(M,bb37)+bb38,0.25);
}


double Lbagbf(double M)
{ /* calculates Lbagb[Msun] -- luminosity at the Base of AGB for star of M, eq.(56) */
 double Ltmp,alfa3;

 if(M<M_HeF) {
   Ltmp=(bb31+bb32*pow(M_HeF,bb33+1.8))/(bb34+pow(M_HeF,bb33)); 
   alfa3=(bb29*pow(M_HeF,bb30)-Ltmp)/Ltmp;
   return bb29*pow(M,bb30)/(1.0+alfa3*exp(15.0*(M-M_HeF))); 
 }  
 else
   return (bb31+bb32*pow(M,bb33+1.8))/(bb34+pow(M,bb33)); 
}


double Rbagbf(double M)
{ /* calculates Rbagb[Rsun] -- radius at the Base of AGB for star of M, sec.5.3 - just after eq.(56) */
 
 return Ragbf(M,Lbagbf(M));
} 



/* ----------------------------------------- EARLY AGB --------------------------------------------- */

double Ahef(void)
{ /* calculates Ahe [Msun Lsun^-1 Myr^-1] constant, I use value from Jarrod code and not from eq.(68) */

 return 8.0e-5; 
}


double tinf1ef(double M)
{ /* calculates integration constant tinf1a[Myr] for EAGB */    
 double tbagb,Lbagb,p,Ahe,D;
 
 tbagb=tbagbf(M);
 Lbagb=Lbagbf(M);
 p=pf(M);
 Ahe=Ahef();
 D=Df(M);

 return tbagb+pow(D/Lbagb,(p-1.0)/p)/((p-1.0)*Ahe*D);
}


double tinf2ef(double M)
{ /* calculates integration constant tinf2a[Myr] for EAGB */    
 double txe,Lx,q,Ahe,B;
 
 txe=txef(M);
 Lx=Lxf(M);
 q=qf(M);
 Ahe=Ahef();
 B=Bf(M);
 
 return txe+pow(B/Lx,(q-1.0)/q)/((q-1.0)*Ahe*B);
}


double txef(double M)
{ /* calculates integration constant txa[Myr] for EAGB */
 double tinf1e,tbagb,Lbagb,Lx,p;
 
 tinf1e=tinf1ef(M);
 tbagb=tbagbf(M);
 Lbagb=Lbagbf(M);
 Lx=Lxf(M);
 p=pf(M);
 
 return tinf1e-pow(Lbagb/Lx,(p-1.0)/p)*(tinf1e-tbagb);
} 


double Mceagbf1(double M, double t)
{ /* calculate the mass of CO Core[Msun] for Early AGB star at time t[Myr], sec. 5.4 */
 double p,q,Ahe,B,D,tinf1e,tinf2e,txe;
 
 txe=txef(M);
 Ahe=Ahef(); 
 
 if(t<=txe) {
   p=pf(M);
   D=Df(M);
   tinf1e=tinf1ef(M);   
   return pow((p-1)*Ahe*D*(tinf1e-t),1.0/(1.0-p));
 }
 else {
   q=qf(M);
   B=Bf(M);
   tinf2e=tinf2ef(M);
   return pow((q-1)*Ahe*B*(tinf2e-t),1.0/(1.0-q));  
 }
}


double Mceagbf2(double M, double L)
{ /* calculate the mass of CO Core[Msun] for Early AGB star of luminosity L[Lsun] */
  /* exactly same formulae as for Red Giant, inverted eq.(37) */
 
 return Mcgbf2(M,L);
}   


double Leagbf1(double M, double t)
{ /* calculates luminosity[Lsun] for Early AGB star of M at time t[Myr] */
  /* !!! something may be wrong with this function, was returning nan for some values */
 double p,q,B,D,Mceagb; 
 
 p=pf(M);
 q=qf(M);
 B=Bf(M);
 D=Df(M);
 Mceagb=Mceagbf1(M,t);

 return min(B*pow(Mceagb,q),D*pow(Mceagb,p));
}


double Leagbf2(double M, double Mcco)
{ /* calculates luminosity[Lsun] for Early AGB star of M with CO core mass Mcco[Msun] */
 double p,q,B,D; 
 
 p=pf(M);
 q=qf(M);
 B=Bf(M);
 D=Df(M);

 return min(B*pow(Mcco,q),D*pow(Mcco,p));
}



/* ---------------------------------- Thermally Pulsing  AGB -------------------------------------------- */ 

double Mcduf(double M)
{ /* calculate the core mass[msun] reduced by second dredge-up on AGB, eq.(69) */
  /* second dredge up valid only for stars with 0.8<Mcbagb<=2.25 */ 
  /* but for Mcbagb<0.8 Mcdu is assumed to be Mbagb, although for these stars there is no 2DU */
 double Mcbagb;
 
 Mcbagb=Mcbagbf(M);
 
 if(Mcbagb<=0.8)
   return Mcbagb;
 else  
   return 0.44*Mcbagb+0.448;
}


double Lduf(double M)
{ /* calculate luminosity[Lsun] at second dredge-up on AGB */
  /* second dredge up valid only for stars with 0.8<Mcbagb<=2.25 */
  double p,q,B,D,Mcdu;
  
 p=pf(M);
 q=qf(M);
 B=Bf(M);
 D=Df(M);
 Mcdu=Mcduf(M);
       
 return min(B*pow(Mcdu,q),D*pow(Mcdu,p));
}        


double tduf(double M)
{ /* calculate time[Myr] of second dredge-up on AGB, eq.(70) */
  /* time of 2DU valid only for stars with 0.8<Mcbagb<=2.25 */
 double p,q,Ahe,B,D,tinf1e,tinf2e,Ldu,Lx;

 Ahe=Ahef();
 Ldu=Lduf(M);
 Lx=Lxf(M);

 if(Ldu<=Lx) {
   p=pf(M);
   D=Df(M);
   tinf1e=tinf1ef(M);
   return tinf1e-pow(D/Ldu,(p-1.0)/p)/((p-1.0)*Ahe*D);
 }
 else {
   q=qf(M);
   B=Bf(M);
   tinf2e=tinf2ef(M);   
   return tinf2e-pow(B/Ldu,(q-1.0)/q)/((q-1.0)*Ahe*B);
 }
}


double Ahhef(void)
{ /* calculates Ahhe [Msun Lsun^-1 Myr^-1] constant, eq.(71) */

 return 1.27e-5;
}


double tinf1tf(double M)
{ /* calculates integration constant tinf1t[Myr] for TPAGB */    
 double tdu,Ldu,p,Ahhe,D;
 
 tdu=tduf(M);
 Ldu=Lduf(M);
 p=pf(M);
 Ahhe=Ahhef();
 D=Df(M);
 
 return tdu+pow(D/Ldu,(p-1.0)/p)/((p-1.0)*Ahhe*D);
}


double tinf2tf(double M)
{ /* calculates integration constant tinf2t[Myr] for TPAGB, eq.(72) */    
 double txt,tdu,Lx,Ldu,q,Ahhe,B;
 
 Lx=Lxf(M);
 Ldu=Lduf(M);
 q=qf(M);
 Ahhe=Ahhef();
 B=Bf(M);
  
 if(Ldu<=Lx) {
   txt=txtf(M);
   return txt+pow(B/Lx,(q-1.0)/q)/((q-1.0)*Ahhe*B);
 }
 else {
   tdu=tduf(M); 
   return tdu+pow(B/Ldu,(q-1.0)/q)/((q-1.0)*Ahhe*B);
 }
}


double txtf(double M)
{ /* calculates integration constant txt[Myr] for TPAGB */
 double tinf1t,tdu,Ldu,Lx,p;
 
 tinf1t=tinf1tf(M);
 tdu=tduf(M);
 Ldu=Lduf(M);
 Lx=Lxf(M);
 p=pf(M);
 
 return tinf1t-pow(Ldu/Lx,(p-1.0)/p)*(tinf1t-tdu);
} 


double Mctagbf1a(double M, double t)
{ /* calculate REAL mass of CO Core=He Core[Msun] for TP AGB star at time t[Myr], sec. 5.4 */
 double Mcdu,Mcd;

 Mcdu=Mcduf(M);
 Mcd=Mctagbf1b(M,t);
  
 return Mcdu+(1.0-lambdaf(M))*(Mcd-Mcdu);
}


double Mctagbf1b(double M, double t)
{ /* calculate DUMMY mass of CO Core=He Core[Msun] for TP AGB star at time t[Myr], sec. 5.4 */
 double p,q,Ahhe,B,D,tinf1t,tinf2t,txt;
 
 txt=txtf(M);
 Ahhe=Ahhef(); 

 if(t<=txt) {
   p=pf(M);
   D=Df(M);
   tinf1t=tinf1tf(M);
   return pow((p-1)*Ahhe*D*(tinf1t-t),1.0/(1.0-p));
 }
 else {
   q=qf(M);
   B=Bf(M);  
   tinf2t=tinf2tf(M);
   return pow((q-1)*Ahhe*B*(tinf2t-t),1.0/(1.0-q));  
 }
}

 
double Mctagbf2(double M, double L)
{ /* calculate the mass of CO Core=He Core[Msun] for TP AGB star of luminosity L[Lsun] */
  /* exactly same formulae as for Red Giant, inverted eq.(37) */
 
 return Mcgbf2(M,L);
}   


double Ltagbf1(double M, double t)
{ /* calculates luminosity[Lsun] for TP AGB star of M at time t[Myr] */
 double p,q,B,D,Mcd; 
 
 p=pf(M);
 q=qf(M);
 B=Bf(M);
 D=Df(M);
 Mcd=Mctagbf1b(M,t);

 return min(B*pow(Mcd,q),D*pow(Mcd,p));
}


double Ltagbf2(double M, double Mc)
{ /* calculates luminosity[Lsun] for TP AGB star of M with core mass Mc[Msun]=Mcco=Mche */
  /* DUMMY Mc shall be supplied, e.g. Mc=Mctagbf1b(M,t) */
 double p,q,B,D; 
 
 p=pf(M);
 q=qf(M);
 B=Bf(M);
 D=Df(M);

 return min(B*pow(Mc,q),D*pow(Mc,p));
}


double lambdaf(double M)
{ /* calculates the mass fraction eaten by He core in teh interpulse period, eq.(73) */   
 
 return min(0.9,0.3+0.001*pow(M,5.0));
} 


double Ragbf(double M, double L)
{ /* calculates radius[Rsun] for AGB star of M and L, eq.(74) */
 double A,tmp1,tmp2;
 
 if(M<=(M_HeF-0.2)) {  
   bb50=bb3;
   A=bb56+bb57*M;
 }   
 else if(M<M_HeF) {
   bb50=inter_line(M_HeF-0.2,bb3,M_HeF,bb55*bb3,M);
   tmp1=min(bb51*pow(M_HeF,-bb52),bb53*pow(M_HeF,-bb54));
   tmp2=bb56+bb57*(M_HeF-0.2);
   A=inter_line(M_HeF-0.2,tmp2,M_HeF,tmp1,M); 
 }
 else {
   bb50=bb55*bb3;
   A=min(bb51*pow(M,-bb52),bb53*pow(M,-bb54));
 } 
       
 return A*(pow(L,bb1)+bb2*pow(L,bb50));
} 



/* ------------------------------------- NAKED HELIUM STARS ----------------------------------------- */

double Lzhsf(double M)
{ /* calculates luminosity[Lsun] for Naked Helium Star at its ZAMS, eq.(77) */
 
 return 15262.0*pow(M,10.25)/(pow(M,9.0)+29.54*pow(M,7.5)+31.18*pow(M,6.0)+0.0469);
}  


double Rzhsf(double M)
{ /* calculates radius[Rsun] for Naked Helium Star at its ZAMS, eq.(78) */
 
 return 0.2391*pow(M,4.6)/(pow(M,4.0)+0.162*pow(M,3.0)+0.0065);
}


double thsmsf(double M)
{ /* calculates lifetime[Myr] of Naked Helium Star on its MS, eq.(79) */
                                                                 
 return (0.4129+18.81*pow(M,4.0)+1.853*pow(M,6.0))/pow(M,6.5);
}


double Lhsmsf(double M, double t)
{ /* calculates luminosity[Lsun] for Naked Helium Star at time t[Myr] during its MS, eq.(80) */
  /* time t calculated from Helum Star ZAMS */
 double Lzhs,thsms,tau,alfa;
 
 Lzhs=Lzhsf(M);
 thsms=thsmsf(M);
 tau=t/thsms;
 alfa=max(0.0,0.85-0.08*M);

 return Lzhs*(1.0+0.45*tau+alfa*tau*tau);
}


double Rhsmsf(double M, double t)
{ /* calculates radius[Rsun] for Naked Helium Star at time t[Myr] during its MS, eq.(81) */
  /* time t calculated from Helum Star ZAMS */
 double Rzhs,thsms,tau,beta;
 
 Rzhs=Rzhsf(M);
 thsms=thsmsf(M);
 tau=t/thsms;
 beta=max(0.0,0.4-0.22*log10(M));

 return Rzhs*(1.0+beta*tau-beta*pow(tau,6.0));
}


double Lthsf(double M)
{ /* calculates luminosity[Lsun] for Naked Helium Star at terminal of its MS, sec.6.1 */

 return Lhsmsf(M,thsmsf(M));
}
        
        
double Bhsf(void)
{ /* calculates Bhs parameter for Naked Helium Giants, sec.6.1 */

 return 4.1e+4;
}

         
double Dhsf(double M)
{ /* calculates Dhs parameter for Naked Helium Giant of M, sec.6.1 */

 return 5.5e+4/(1.0+0.4*pow(M,4.0));                              
} 


double phsf(void)
{ /* returns parameter phs for Naked Helium Giants, sec.6.1 from eq.(84) */
 
 return 5.0;
} 


double qhsf(void)
{ /* returns parameter qhs for Naked Helium Giants, sec.6.1 from eq.(84) */
 
 return 3.0;
} 


double tinf1hsf(double M)
{ /* calculates integration constant tinf1hs[Myr] for Naked Helium Giant of M, sec.6.1 */
 double thsms,Lths,phs,Ahe,Dhs;
 
 thsms=thsmsf(M);
 Lths=Lthsf(M);
 phs=phsf();
 Ahe=Ahef();
 Dhs=Dhsf(M);
      
 return thsms+pow(Dhs/Lths,(phs-1.0)/phs)/((phs-1.0)*Ahe*Dhs);
}


double tinf2hsf(double M)
{ /* calculates integration constant tinf2hs[Myr] for Naked Helium Giant of M, sec.6.1 */
 double thsms,txhs,Lxhs,qhs,Ahe,Bhs;
 
 thsms=thsmsf(M);
 txhs=txhsf(M);
 Lxhs=Lxhsf(M);
 qhs=qhsf();
 Ahe=Ahef();
 Bhs=Bhsf();
      
 return txhs+pow(Bhs/Lxhs,(qhs-1.0)/qhs)/((qhs-1.0)*Ahe*Bhs);
}
       

double txhsf(double M)
{ /* calculates integration constant txhs[Myr] for Naked Helium Giant of M, sec.6.1 */
 double tinf1hs,thsms,Lths,Lxhs,phs;
  
 tinf1hs=tinf1hsf(M);
 thsms=thsmsf(M);
 Lths=Lthsf(M);
 Lxhs=Lxhsf(M);
 phs=phsf();
       
 return tinf1hs-pow(Lths/Lxhs,(phs-1.0)/phs)*(tinf1hs-thsms);
}
               
       
double Lxhsf(double M) 
{ /* calculates luminosity Lxhs[Lsun] for Naked Helium Giant of M, sec.6.1 */
 double phs,qhs,Bhs,Dhs,Mxhs;
  
 phs=phsf();
 qhs=qhsf();
 Bhs=Bhsf();
 Dhs=Dhsf(M);
 Mxhs=Mxhsf(M);
            
 return min(Bhs*pow(Mxhs,qhs),Dhs*pow(Mxhs,phs));
}
                    
                    
double Mxhsf(double M)
{ /* calculates core mass Mxhs[Msun] for Naked Helium Giant of M, sec.6.1 */
 double phs,qhs,Bhs,Dhs;
  
 phs=phsf();
 qhs=qhsf();
 Bhs=Bhsf();
 Dhs=Dhsf(M);
     
 return pow(Bhs/Dhs,1.0/(phs-qhs));
}
        
        
double thsevolf(double M, double Mc)
{ /* calculate for Helium giant of mass M and CO core mass Mc how evolved it is, inversed eq.(39) */
  /* used for HS Giants -sec.6.1, so gives back virtual evolutionary time for Helium Giant of M,Mc */
  /* M0 for Helium Giant shall be supplied as M */
 double Bhs,Dhs,Ahe,phs,qhs,Mxhs,tinf1hs,tinf2hs;          

 Ahe=Ahef();
 Mxhs=Mxhsf(M);

 if(Mc<=Mxhs) {
   Dhs=Dhsf(M);
   phs=phsf();
   tinf1hs=tinf1hsf(M);
   return tinf1hs-pow(Mc,1.0-phs)/((phs-1)*Ahe*Dhs); 
 }
 else {
   Bhs=Bhsf();
   qhs=qhsf();
   tinf2hs=tinf2hsf(M);
   return tinf2hs-pow(Mc,1.0-qhs)/((qhs-1)*Ahe*Bhs); 
 }
}
        

double Mchsgbf1(double M, double t)
{ /* calculates CO core mass[Msun] for Naked Helium Giant at time t[Myr], sec.6.1 */         
 double Bhs,Dhs,Ahe,phs,qhs,txhs,tinf1hs,tinf2hs;
 
 Ahe=Ahef();
 txhs=txhsf(M);
 
 if(t<=txhs) {
   Dhs=Dhsf(M);
   phs=phsf();
   tinf1hs=tinf1hsf(M);
   return pow((phs-1)*Ahe*Dhs*(tinf1hs-t),1.0/(1.0-phs));
 }
 else {
   Bhs=Bhsf();
   qhs=qhsf();
   tinf2hs=tinf2hsf(M);
   return pow((qhs-1)*Ahe*Bhs*(tinf2hs-t),1.0/(1.0-qhs));
 }
}


double Mchsgbf2(double M, double L)
{ /* calculates CO core mass[Msun] for Naked Helium Giant of luminosity L[Lsun], sec.6.1 */         
 
 if(L<=Lxhsf(M))
   return pow(L/Dhsf(M),1.0/phsf());
 else
   return pow(L/Bhsf(),1.0/qhsf());
}


double Lhsgbf(double M, double Mc)
{ /* calculates luminosity[Lsun] for Naked Helium Giant of M and Mcco[Msun], eq.(84) */
  /* time t calculated from Helum Star ZAMS */
 double Bhs,Dhs,phs,qhs;
 
 Bhs=Bhsf();
 Dhs=Dhsf(M); 
 phs=phsf();
 qhs=qhsf();
 
 return min(Bhs*pow(Mc,qhs),Dhs*pow(Mc,phs)); 
}  


double Rhsgbf(double M, double L, double Lths, double Rzhs, int *K)
{ /* calculates radius[Rsun] for Naked Helium Giant of M and L[Lsun], eqs.(85,86,87,88) */
  /* time t calculated from Helum Star ZAMS, sets stellar type K, see after eq.(88) */
  /* Lths supplied from outside as it depends on diffrent masses at diffrent Rhsgbf() calls: */
  /* in HSGBf() it is M0, in perturb() it is not M0, Rzhs from outside for compability with Lths */
 double lambda,R1,R2;
 
 lambda=500.0*(2.0+pow(M,5.0))/pow(M,2.5);
 R1=Rzhs*pow(L/Lths,0.2)+0.02*(exp(L/lambda)-exp(Lths/lambda));
 R2=0.08*pow(L,0.75);  
 
 if(R1<R2) { 
   (*K)=8;
   return R1;
 }  
 else {
   (*K)=9;  
   return R2;
 }
}


double Mcmaxf(double M)
{ /* calculate the core mass of Helium Giant (K=8,9) at which C shell burning stops, eq.(89) */ 
  /* Helium Giant of M turns into CO WD, as M current mass Mt shall be supplied */

 return min(1.45*M-0.31,M);
}



/* ------------------------------------- COMPACT OBJECTS -------------------------------------------- */

double compactf0(double M, double Mco)
{ /* returns gravitational mass of compact object [Msun]: Timmes et al. 1996, ApJ 457, 834 */
  /* supplied current mass M and CO core mass Mcco */
  /* fall back limits (Lim1=5,Lim2=7.6)-correspond to final Mco of full wind star of Mzams=18-20 and 40-42 Msun */ 
  /* Mzams, ZAMS mass of a star is not needed here */
  /* if winds with ETA!=0.5 or different (from Jarrod) wind prescriptions are used */
  /* then new values of Mlim1,Mlim2,a2,b2 are needed */ 
  /* in FB I put information on fall back -- needed for sn_type() and get_Vkick2,5,6() */
 double a1,b1,a2,b2,a3,b3,c3,del3,Mcfn,Msn,Mgrav,frac;
 double Mfb,Mfinal;
 double Mlim1=4.999889;                            /* above this CO core mass of progenitor: fall back occurs  */
 double Mlim2=7.647335;                      /* above this CO core mass of progenitor: prompt black hole forms */
                                            /* Mlim1 corresponds to CO core mass of Mzams=20Msun at time of SN */
                                            /* Mlim2 corresponds to CO core mass of Mzams=42Msun at time of SN */ 
                                                                /* see Fryer & Kalogera 1999, astro-ph/9911312 */ 
 if((Mco-acc)>M)
   fprintf(fp0,"error: Mco>M in compactf0()\n");
 
 if(Mco<4.280)                    /* Mzams<18 Msun */ 
   Mcfn=1.496;
 else if(Mco<=6.312)              /* 18<Mzams<25 Msun */
   Mcfn=2.113;  
 else if(Mco<=6.745) {
   a1=0.692;
   b1=-2.257;
   Mcfn=a1*Mco+b1;
 }  
 else {                           /* Mzams>25 Msun */ 
   a1=0.368;                      /* coeffs for FeNi core calculations */
   b1=-0.071;
   Mcfn=a1*Mco+b1;                      
 }
 
 if(Mco<=Mlim1) {
   frac=0.0;                      /* no fall back */
   Msn=Mcfn;                      /* Mcompact is equal to mass of collapsing FeNi core, full SN explosion occurs */ 
   Mfb=0.0;                       /* fall back mass */
 }
 else if(Mco<Mlim2) {
   a2=0.377723;                   /* coeffs for fall back calculations */
   b2=-1.888571;                  /* partial fall back */
   frac=a2*Mco+b2;                /* Mcompact is equal to mass of collapsing FeNi core + mass of fall back matter */
   Msn=Mcfn+frac*(M-Mcfn);        /* fall back is 0 for Mco=Mlim1 and complete for Mco=Mlim2, smaller SN explosion */
   Mfb=frac*(M-Mcfn);
 }
 else {
   frac=1.0;                      /* direct collapse */ 
   Msn=M;                         /* Mcompact is equal to total mass of star during collapse -prompt BH formation */
   Mfb=(M-Mcfn);
 }

 if(frac<=0.001) FB=0;                    /* no fall back */
 else if(frac>0.001 && frac<0.999) FB=1;  /* partial fall back */
 else FB=2;                               /* full fall back */
 
 if(kickback==1 && Mfb<kickback0) FB=0;     /* parm. study: full kick (fall back override) */


 if(0.8*Msn<Mmaxns) {           /* for NS progenitors */
   a3=0.075;                    /* for Msn=3.0 (0.8*Msn=2.4 Msun) -> Mgrav=2.523 Msun */
   b3=1.0;         
   c3=-(Msn);        
   del3=b3*b3-4.0*a3*c3;
   Mgrav=(-b3+sqrt(del3))/(2.0*a3);   
   if(Mgrav>Mmaxns)                  /* light BH */
     Mgrav=0.9*Msn;
 }
 else                                /* for BH progenitors) */
   Mgrav=0.9*Msn;  


 FRAC=frac;                           /* writes how much fallback [0:1] took place into global var. */
 Mfinal=min(M,Mgrav);                 /* to make sure that NS or BH mass do not exceed mass M of outbursting star */

 if(SimpleNS==1 && Mfinal<Mmaxns) {          /* simple final mass prescription for NS (Mzams=8.275-20.88 Msun): */
   Mfinal=0.35*Mco+0.596;                    /* for Mco=1.44 (Mzams=8.275) Mns=1.1 Msun */
   if(Mfinal>Mmaxns) Mfinal=Mmaxns-0.0001;   /* for Mco=5.44 (Mzams=20.88) Mns=2.5 Msun */
   if(Mfinal<1.1) Mfinal=1.1;                /* linear interpolation inbetween */ 
 }

 return Mfinal;
}


double compactf1(double M, double Mco)
{ /* returns gravitational mass of compact object [Msun]: Chris F./Chris B./Vicky K. 2010 */
  /* based on RAPID SN mechanism: Chris F. compiled data from Heger&Woosley solar metallcity models */
  /* supplied current mass M and CO core mass Mco */
  /* in fb I put information on fall back -- needed for sn_type() */
  /* the actual fraction of fall back [0-1] is written into global variable FRAC */
 double a1,b1,c1,del1;
 double Mex,Menv,Mfb,Mbar,Mgrav,Mfinal,frac;

 if((Mco-acc)>M)
   fprintf(fp0,"error: Mco>M in compactf1()\n");
    
 Mex=1.0;            /* Mex [Msun] -- baryon mass of the proto-NS at the time of bounce/explosion */
                     /* the same Mex mass for all models */

 Menv=M-Mex;         /* Menv [Msun] -- baryon mass of the "envelope": mass over proto-NS */

 if(Mex>M) {         /* Mex>M: should not happen here: but look at comments at compactf2(): same place */
   Mex=M;            /* remnant can not be more massive than the mass of exploding star */
   Menv=0.0;         /* it is either a very low mass NS with almost no SN ejecta (as treated here) */
   Mfb=0.0;          /* or it could potentially never explode and form a CO/ONeMg WD */
   frac=0.0;
 }
 else if(Mco<2.5) {
   Mfb=0.20;                     /* Mfb [Msun] -- baryon mass of the fall back material accreted onto proto-NS */ 
   Mfb=min(Mfb,Menv);            /* one can not have fall back larger than mass contained in the envelope */ 
   frac=Mfb/Menv;                /* frac [0-1] -- fallback parameter: fraction of envelope accreted onto proto-NS */
 } 
 else if(Mco<6.0) {              /* Mfb=0.2Msun (Mco=2.5Msun) -- 1.2Msun (Mco=6.0): linear interpolation */
   a1=0.286; 
   b1=-0.514;; 
   Mfb=a1*Mco+b1;
   Mfb=min(Mfb,Menv);
   frac=Mfb/Menv;
 }
 else if(Mco<7.0) {              /* Mfb=1.0 Mco=6.0--7.0 Msun: full fallback: unexplained feature in Heger&Woosley models */
   frac=1.0;                     /* but A.Heger claims that it is real and has something to do with structure of star in */ 
   Mfb=Menv;                     /* nparticular mass range */
 }
 else if(Mco<11.0) {             /* frac is chosen such that at Mco=7.0 frac results in Mfb=5.1 Msun */
   a1=(1.0-(5.1/(M-Mex)))/4.0;   /* frac is chosen such that at Mco=11.0 frac results in Mfb=Menv (frac=1.0) */
   b1=1.0-11.0*a1;               /* and frac: linear interpolation for 7.0<Mco<11 */
   frac=a1*Mco+b1; 
   Mfb=frac*Menv;
   Mfb=min(Mfb,Menv);
 } 
 else {                          /* full fall back: entire star goes into remnant */
   frac=1.0; 
   Mfb=Menv; 
 }

 if(frac<=0.001) FB=0;                    /* no fall back */
 else if(frac>0.001 && frac<0.999) FB=1;  /* partial fall back */
 else FB=2;                               /* full fall back */
 
 if(kickback==0 && Mfb<0.2) FB=0;                /* stdandard: full kick (fall back override) */
 else if(kickback==1 && Mfb<kickback2) FB=0;     /* parm. study: full kick (fall back override) */
 
 Mbar=Mex+Mfb;                   /* Mbar [Msun] -- final baryon mass of the remnant */
 if(Mbar>M+acc) 
   fprintf(fp0,"compact3f(): error -- something is wrong with the mass of the remnant\n");
 
 if(0.8*Mbar<Mmaxns) {           /* for NS progenitors */
   a1=0.075;                     /* for Mbar=3.0 (0.8*Mbar=2.4 Msun) -> Mgrav=2.523 Msun */
   b1=1.0;         
   c1=-(Mbar);        
   del1=b1*b1-4.0*a1*c1;
   Mgrav=(-b1+sqrt(del1))/(2.0*a1);   
   if(Mgrav>Mmaxns)                  /* light BH */
     Mgrav=(1.0-Fneut)*Mbar;
 }
 else                                /* for BH progenitors */
   Mgrav=(1.0-Fneut)*Mbar;  

 FRAC=frac;                           /* writes how much fallback [0:1] took place into global var. */
 Mfinal=min(M,Mgrav);                 /* to make sure that NS or BH mass do not exceed mass M of exploding star */

 return Mfinal;
}


double compactf2(double M, double Mco)
{ /* returns gravitational mass of compact object [Msun]: Chris F./Chris B./Vicky K. 2010 */
  /* based on DELAYED (SASI) SN mechanism: Chris F. compiled his and other ppl data */
  /* supplied current mass M and CO core mass Mco */
  /* in FB I put information on fall back -- needed for sn_type() and get_Vkick2,5,6() */
  /* the actual fraction of fall back [0-1] is written into global variable FRAC */
 double a1,b1,c1,del1;
 double Mex,Menv,Mfb,Mbar,Mgrav,Mfinal,frac;

 if((Mco-acc)>M)
   fprintf(fp0,"error: Mco>M in compactf2()\n");
    
 if(Mco<3.5) Mex=1.20;            /* Mex [Msun] -- baryon mass of the proto-NS at the time of bounce/explosion */
 else if(Mco<6.0) Mex=1.30; 
 else if(Mco<11.0) Mex=1.40;
 else Mex=1.60; 

 Menv=M-Mex;                      /* Menv [Msun] -- baryon mass of the "envelope": mass over proto-NS */

 if(Mex>M) {                      /* Mex>M: can happen for lowest mass stars K=8 M=2.0 with removed envelope: Mco=1.15 */
   Mex=M;                         /* remnant can not be more massive than the mass of exploding star */
   Menv=0.0;                      /* it is either a very low mass NS with almost no SN ejecta (as treated here) */
   Mfb=0.0;                       /* or it could potentially never explode and form a CO/ONeMg WD */
   frac=0.0;
 }
 else if(Mco<2.5) {
   Mfb=0.20;                     /* Mfb [Msun] -- baryon mass of the fall back material accreted onto proto-NS */ 
   Mfb=min(Mfb,Menv);            /* one can not have fall back larger than mass contained in the envelope */ 
   frac=Mfb/Menv;                /* frac [0-1] -- fallback parameter: fraction of envelope accreted onto proto-NS */
 } 
 else if(Mco<3.5) {              /* Mfb=0.2 (Mco=2.5Msun) -- 0.7 (Mco=3.5) Msun: linear interpolation */
   a1=0.5; 
   b1=-1.05; 
   Mfb=a1*Mco+b1;
   Mfb=min(Mfb,Menv);
   frac=Mfb/Menv;
 }
 else if(Mco<11.0) {             /* frac is chosen such that at Mco=7.0 frac results in Mfb=0.7 Msun */
   a1=(1.0-(0.7/(M-Mex)))/7.5;   /* frac is chosen such that at Mco=11.0 frac results in Mfb=Menv (frac=1.0) */
   b1=1.0-11.0*a1;               /* and frac: linear interpolation for 3.5<Mco<11 */
   frac=a1*Mco+b1; 
   Mfb=frac*Menv;
   Mfb=min(Mfb,Menv);
 } 
 else {                          /* full fall back: entire star goes into remnant */
   frac=1.0; 
   Mfb=Menv; 
 }

 if(frac<=0.001) FB=0;                    /* no fall back */
 else if(frac>0.001 && frac<0.999) FB=1;  /* partial fall back */
 else FB=2;                               /* full fall back */
 
 if(kickback==0 && Mfb<0.2) FB=0;                /* stdandard: full kick (fall back override) */
 else if(kickback==1 && Mfb<kickback1) FB=0;     /* parm. study: full kick (fall back override) */
 
 Mbar=Mex+Mfb;                   /* Mbar [Msun] -- final baryon mass of the remnant */
 if(Mbar>M+acc) 
   fprintf(fp0,"compact2f(): error -- something is wrong with the mass of the remnant\n");
 
 if(0.8*Mbar<Mmaxns) {           /* for NS progenitors */
   a1=0.075;                     /* for Mbar=3.0 (0.8*Mbar=2.4 Msun) -> Mgrav=2.523 Msun */
   b1=1.0;         
   c1=-(Mbar);        
   del1=b1*b1-4.0*a1*c1;
   Mgrav=(-b1+sqrt(del1))/(2.0*a1);   
   if(Mgrav>Mmaxns)                  /* light BH */
     Mgrav=(1.0-Fneut)*Mbar;
 }
 else                                /* for BH progenitors */
   Mgrav=(1.0-Fneut)*Mbar;  

 FRAC=frac;                           /* writes how much fallback [0:1] took place into global var. */
 Mfinal=min(M,Mgrav);                 /* to make sure that NS or BH mass do not exceed mass M of exploding star */

 return Mfinal;
}


double compactf3(double M, double *Mco, double *Mhe)
{ /* RAPID SN + strong PPSN/PSN + adjustable BH neutrino emission + adjustable increase of CO size */
  /* returns gravitational mass of compact object [Msun]: June 13, 2019 */
 double a1,b1,c1,del1;
 double Mex,Menv,Mfb,Mbar,Mgrav,Mfinal,frac;
 double Zcrit11,delMsn1c;

 
 Zcrit11=Zcrit1+0.00001;      /* ZZ<=Zcrit1+0.00001; to make sure that model with ZZ=Zcrit1 is affected by this: */
 
                                               /* only for accretion above Fcrit1 and only for low-Z stars */
 if(Spin1==1 && Fms>Fcrit1 && ZZ<=Zcrit11) {   /* increase of CO core mass due to MS progenitor spin up from accretion */
   (*Mhe)=min(Fincr1*(*Mhe),M);                /* do not increase over entire star mass */ 
   (*Mco)=min(Fincr1*(*Mco),M);                        
 }


 if(PPSN==1 && (*Mhe)>=45.0 && (*Mhe)<65.00) {         /* PPSN regime: remove mass from the star */ 
   M=(*Mhe)=45.0;                                      /* before calculating mass of the remnant */ 
   if((*Mco)>45.0) (*Mco)=45.0;
 }
 else if(PSN==1 && (*Mhe)>=65.0 && (*Mhe)<135.0) {     /* PSN regime: explode the star and no remnant */
   Mfinal=M=0.0;
   return Mfinal;                                      /* exit function here */
 }

 if(((*Mco)-acc)>M)
   fprintf(fp0,"error: Mco>M in compactf3()\n");
    
 Mex=1.0;            /* Mex [Msun] -- baryon mass of the proto-NS at the time of bounce/explosion */
                     /* the same Mex mass for all models */

 Menv=M-Mex;         /* Menv [Msun] -- baryon mass of the "envelope": mass over proto-NS */

 if(Mex>M) {         /* Mex>M: should not happen here: but look at comments at compactf2(): same place */
   Mex=M;            /* remnant can not be more massive than the mass of exploding star */
   Menv=0.0;         /* it is either a very low mass NS with almost no SN ejecta (as treated here) */
   Mfb=0.0;          /* or it could potentially never explode and form a CO/ONeMg WD */
   frac=0.0;
 }
 else if((*Mco)<2.5) {
   Mfb=0.20;                     /* Mfb [Msun] -- baryon mass of the fall back material accreted onto proto-NS */ 
   Mfb=min(Mfb,Menv);            /* one can not have fall back larger than mass contained in the envelope */ 
   frac=Mfb/Menv;                /* frac [0-1] -- fallback parameter: fraction of envelope accreted onto proto-NS */
 } 
 else if((*Mco)<6.0) {           /* Mfb=0.2Msun (Mco=2.5Msun) -- 1.2Msun (Mco=6.0): linear interpolation */
   a1=0.286; 
   b1=-0.514;; 
   Mfb=a1*(*Mco)+b1;
   Mfb=min(Mfb,Menv);
   frac=Mfb/Menv;
 }
 else if((*Mco)<7.0) {           /* Mfb=1.0 Mco=6.0--7.0 Msun: full fallback: unexplained feature in Heger&Woosley models */
   frac=1.0;                     /* but A.Heger claims that it is real and has something to do with structure of star in */ 
   Mfb=Menv;                     /* nparticular mass range */
 }
 else if((*Mco)<11.0) {          /* frac is chosen such that at Mco=7.0 frac results in Mfb=5.1 Msun */
   a1=(1.0-(5.1/(M-Mex)))/4.0;   /* frac is chosen such that at Mco=11.0 frac results in Mfb=Menv (frac=1.0) */
   b1=1.0-11.0*a1;               /* and frac: linear interpolation for 7.0<Mco<11 */
   frac=a1*(*Mco)+b1; 
   Mfb=frac*Menv;
   Mfb=min(Mfb,Menv);
 } 
 else {                          /* full fall back: entire star goes into remnant */
   frac=1.0; 
   Mfb=Menv; 
 }

 if(frac<=0.001) FB=0;                    /* no fall back */
 else if(frac>0.001 && frac<0.999) FB=1;  /* partial fall back */
 else FB=2;                               /* full fall back */
 
 if(kickback==0 && Mfb<0.2) FB=0;                /* stdandard: full kick (fall back override) */
 else if(kickback==1 && Mfb<kickback2) FB=0;     /* parm. study: full kick (fall back override) */
 
 Mbar=Mex+Mfb;                   /* Mbar [Msun] -- final baryon mass of the remnant */
 if(Mbar>M+acc) 
   fprintf(fp0,"compact4f(): error -- something is wrong with the mass of the remnant\n");
 
 if(0.8*Mbar<Mmaxns) {           /* for NS progenitors */
   a1=0.075;                     /* for Mbar=3.0 (0.8*Mbar=2.4 Msun) -> Mgrav=2.523 Msun */
   b1=1.0;         
   c1=-(Mbar);        
   del1=b1*b1-4.0*a1*c1;
   Mgrav=(-b1+sqrt(del1))/(2.0*a1);   
   if(Mgrav>Mmaxns)                  /* light BH */
     Mgrav=(1.0-Fneut1)*Mbar;
 }
 else                                /* for BH progenitors */
   Mgrav=(1.0-Fneut2)*Mbar;  

 FRAC=frac;                           /* writes how much fallback [0:1] took place into global var. */
 Mfinal=min(M,Mgrav);                 /* to make sure that NS or BH mass do not exceed mass M of exploding star */

 if(GRB==1 && FB==2 && Mfinal>Mbhmin && Mfinal<Mbhmax) {   /* for massive BHs with full fall back */
   delMsn1c=get_t(mgrbmin,mgrbmax);     /* mass loss range [Msun] from Ashall et al. 2019: arXiv:1702.04339 */ 
   Mfinal=Mfinal-delMsn1c;              /* apply addtional mass loss due to potential SNIc/GRB */
 }

 return Mfinal;
}


double compactf4(double M, double *Mco, double *Mhe)
{ /* DELAYED SN + strong PPSN/PSN + adjustable BH neutrino emission + adjustable increase of CO size */
  /* returns gravitational mass of compact object [Msun]: Jun 13, 2019 */
 double a1,b1,c1,del1;
 double Mex,Menv,Mfb,Mbar,Mgrav,Mfinal,frac;
 double Zcrit11,delMsn1c;
 
 Zcrit11=Zcrit1+0.00001;      /* ZZ<=Zcrit1+0.00001; to make sure that model with ZZ=Zcrit1 is affected by this: */

                              /* only for accretion above Fcrit1 and only for low-Z stars */
 if(Spin1==1 && Fms>Fcrit1 && ZZ<=Zcrit11) {   /* increase of CO core mass due to MS progenitor spin up from accretion */
   (*Mhe)=min(Fincr1*(*Mhe),M);                /* do not increase over entire star mass */
   (*Mco)=min(Fincr1*(*Mco),M);
 }


 if(PPSN==1 && (*Mhe)>=45.0 && (*Mhe)<65.00) {      /* PPSN regime: remove mass from the star */ 
   M=(*Mhe)=45.0;                                /* before calculating mass of the remnant */ 
   if((*Mco)>45.0) (*Mco)=45.0;
 }
 else if(PSN==1 && (*Mhe)>=65.0 && (*Mhe)<135.0) {  /* PSN regime: explode the star and no remnant */
   Mfinal=M=0.0;
   return Mfinal;                             /* exit function here */
 }

 if(((*Mco)-acc)>M)
   fprintf(fp0,"error: Mco>M in compactf4()\n");
 if((*Mco)<3.5) Mex=1.20;            /* Mex [Msun] -- baryon mass of the proto-NS at the time of bounce/explosion */
 else if((*Mco)<6.0) Mex=1.30; 
 else if((*Mco)<11.0) Mex=1.40;
 else Mex=1.60; 

 Menv=M-Mex;                      /* Menv [Msun] -- baryon mass of the "envelope": mass over proto-NS */

 if(Mex>M) {                      /* Mex>M: can happen for lowest mass stars K=8 M=2.0 with removed envelope: Mco=1.15 */
   Mex=M;                         /* remnant can not be more massive than the mass of exploding star */
   Menv=0.0;                      /* it is either a very low mass NS with almost no SN ejecta (as treated here) */
   Mfb=0.0;                       /* or it could potentially never explode and form a CO/ONeMg WD */
   frac=0.0;
 }
 else if((*Mco)<2.5) {
   Mfb=0.20;                     /* Mfb [Msun] -- baryon mass of the fall back material accreted onto proto-NS */ 
   Mfb=min(Mfb,Menv);            /* one can not have fall back larger than mass contained in the envelope */ 
   frac=Mfb/Menv;                /* frac [0-1] -- fallback parameter: fraction of envelope accreted onto proto-NS */
 } 
 else if((*Mco)<3.5) {              /* Mfb=0.2 (Mco=2.5Msun) -- 0.7 (Mco=3.5) Msun: linear interpolation */
   a1=0.5; 
   b1=-1.05; 
   Mfb=a1*(*Mco)+b1;
   Mfb=min(Mfb,Menv);
   frac=Mfb/Menv;
 }
 else if((*Mco)<11.0) {             /* frac is chosen such that at Mco=7.0 frac results in Mfb=0.7 Msun */
   a1=(1.0-(0.7/(M-Mex)))/7.5;   /* frac is chosen such that at Mco=11.0 frac results in Mfb=Menv (frac=1.0) */
   b1=1.0-11.0*a1;               /* and frac: linear interpolation for 3.5<Mco<11 */
   frac=a1*(*Mco)+b1; 
   Mfb=frac*Menv;
   Mfb=min(Mfb,Menv);
 } 
 else {                          /* full fall back: entire star goes into remnant */
   frac=1.0; 
   Mfb=Menv; 
 }

 if(frac<=0.001) FB=0;                    /* no fall back */
 else if(frac>0.001 && frac<0.999) FB=1;  /* partial fall back */
 else FB=2;                               /* full fall back */
 
 if(kickback==0 && Mfb<0.2) FB=0;                /* stdandard: full kick (fall back override) */
 else if(kickback==1 && Mfb<kickback1) FB=0;     /* parm. study: full kick (fall back override) */
 
 Mbar=Mex+Mfb;                   /* Mbar [Msun] -- final baryon mass of the remnant */
 if(Mbar>M+acc) 
   fprintf(fp0,"compact2f(): error -- something is wrong with the mass of the remnant\n");
 
 if(0.8*Mbar<Mmaxns) {           /* for NS progenitors */
   a1=0.075;                     /* for Mbar=3.0 (0.8*Mbar=2.4 Msun) -> Mgrav=2.523 Msun */
   b1=1.0;         
   c1=-(Mbar);        
   del1=b1*b1-4.0*a1*c1;
   Mgrav=(-b1+sqrt(del1))/(2.0*a1);   
   if(Mgrav>Mmaxns)                  /* light BH */
     Mgrav=(1.0-Fneut1)*Mbar;
 }
 else                                /* for BH progenitors */
   Mgrav=(1.0-Fneut2)*Mbar;  

 FRAC=frac;                           /* writes how much fallback [0:1] took place into global var. */
 Mfinal=min(M,Mgrav);                 /* to make sure that NS or BH mass do not exceed mass M of exploding star */

 if(GRB==1 && FB==2 && Mfinal>Mbhmin && Mfinal<Mbhmax) {   /* for massive BHs with full fall back */
   delMsn1c=get_t(mgrbmin,mgrbmax);     /* mass loss range [Msun] from Ashall et al. 2019: arXiv:1702.04339 */ 
   Mfinal=Mfinal-delMsn1c;              /* apply addtional mass loss due to potential SNIc/GRB */
 }

 return Mfinal;
}


double compactf5(double M, double *Mco, double *Mhe)
{ /* RAPID SN + moderate PPSN/PSN + adjustable BH neutrino emission + adjustable increase of CO size */
  /* returns gravitational mass of compact object [Msun]: Jun 13, 2019 */
 double a1,b1,c1,del1;
 double Mcc,Mex,Menv,Mfb,Mbar,Mgrav,Mfinal,frac;
 double Zcrit11,delMsn1c;
 
 Zcrit11=Zcrit1+0.00001;      /* ZZ<=Zcrit1+0.00001; to make sure that model with ZZ=Zcrit1 is affected by this: */
 
                                               /* only for accretion above Fcrit1 and only for low-Z stars */
 if(Spin1==1 && Fms>Fcrit1 && ZZ<=Zcrit11) {   /* increase of CO core mass due to MS progenitor spin up from accretion */
   (*Mhe)=min(Fincr1*(*Mhe),M);                /* do not increase over entire star mass */ 
   (*Mco)=min(Fincr1*(*Mco),M);                        
 }

 if(PPSN==1 && (*Mhe)>=40.0 && (*Mhe)<65.00) {           /* PPSN regime: remove mass from the star */ 
   if((*Mhe)<60.00) Mcc=0.65*(*Mhe)+12.2;                 /* Mcc is mass of the colapsed helium core - pulse mass */
   else if((*Mhe)>=60.0 && (*Mhe)<62.5) Mcc=51.2;
   else Mcc=-14.3*(*Mhe)+938.0;
   M=(*Mhe)=Mcc;
   if((*Mco)>Mcc) (*Mco)=Mcc; 
 }
 else if(PSN==1 && (*Mhe)>=65.0 && (*Mhe)<135.0) {       /* PSN regime: explode the star and no remnant */
   Mfinal=M=0.0;
   return Mfinal;                                        /* exit function here */
 }

 if(((*Mco)-acc)>M)
   fprintf(fp0,"error: Mco>M in compactf5()\n");
    
 Mex=1.0;            /* Mex [Msun] -- baryon mass of the proto-NS at the time of bounce/explosion */
                     /* the same Mex mass for all models */

 Menv=M-Mex;         /* Menv [Msun] -- baryon mass of the "envelope": mass over proto-NS */

 if(Mex>M) {         /* Mex>M: should not happen here: but look at comments at compactf2(): same place */
   Mex=M;            /* remnant can not be more massive than the mass of exploding star */
   Menv=0.0;         /* it is either a very low mass NS with almost no SN ejecta (as treated here) */
   Mfb=0.0;          /* or it could potentially never explode and form a CO/ONeMg WD */
   frac=0.0;
 }
 else if((*Mco)<2.5) {
   Mfb=0.20;                     /* Mfb [Msun] -- baryon mass of the fall back material accreted onto proto-NS */ 
   Mfb=min(Mfb,Menv);            /* one can not have fall back larger than mass contained in the envelope */ 
   frac=Mfb/Menv;                /* frac [0-1] -- fallback parameter: fraction of envelope accreted onto proto-NS */
 } 
 else if((*Mco)<6.0) {           /* Mfb=0.2Msun (Mco=2.5Msun) -- 1.2Msun (Mco=6.0): linear interpolation */
   a1=0.286; 
   b1=-0.514;; 
   Mfb=a1*(*Mco)+b1;
   Mfb=min(Mfb,Menv);
   frac=Mfb/Menv;
 }
 else if((*Mco)<7.0) {           /* Mfb=1.0 Mco=6.0--7.0 Msun: full fallback: unexplained feature in Heger&Woosley models */
   frac=1.0;                     /* but A.Heger claims that it is real and has something to do with structure of star in */ 
   Mfb=Menv;                     /* nparticular mass range */
 }
 else if((*Mco)<11.0) {          /* frac is chosen such that at Mco=7.0 frac results in Mfb=5.1 Msun */
   a1=(1.0-(5.1/(M-Mex)))/4.0;   /* frac is chosen such that at Mco=11.0 frac results in Mfb=Menv (frac=1.0) */
   b1=1.0-11.0*a1;               /* and frac: linear interpolation for 7.0<Mco<11 */
   frac=a1*(*Mco)+b1; 
   Mfb=frac*Menv;
   Mfb=min(Mfb,Menv);
 } 
 else {                          /* full fall back: entire star goes into remnant */
   frac=1.0; 
   Mfb=Menv; 
 }

 if(frac<=0.001) FB=0;                    /* no fall back */
 else if(frac>0.001 && frac<0.999) FB=1;  /* partial fall back */
 else FB=2;                               /* full fall back */
 
 if(kickback==0 && Mfb<0.2) FB=0;                /* stdandard: full kick (fall back override) */
 else if(kickback==1 && Mfb<kickback2) FB=0;     /* parm. study: full kick (fall back override) */
 
 Mbar=Mex+Mfb;                   /* Mbar [Msun] -- final baryon mass of the remnant */
 if(Mbar>M+acc) 
   fprintf(fp0,"compactf5f(): error -- something is wrong with the mass of the remnant\n");
 
 if(0.8*Mbar<Mmaxns) {           /* for NS progenitors */
   a1=0.075;                     /* for Mbar=3.0 (0.8*Mbar=2.4 Msun) -> Mgrav=2.523 Msun */
   b1=1.0;         
   c1=-(Mbar);        
   del1=b1*b1-4.0*a1*c1;
   Mgrav=(-b1+sqrt(del1))/(2.0*a1);   
   if(Mgrav>Mmaxns)                    /* light BH */
     Mgrav=(1.0-Fneut1)*Mbar;
 }
 else                                  /* for BH progenitors */
   Mgrav=(1.0-Fneut2)*Mbar;  

 FRAC=frac;                            /* writes how much fallback [0:1] took place into global var. */
 Mfinal=min(M,Mgrav);                  /* to make sure that NS or BH mass do not exceed mass M of exploding star */

 if(GRB==1 && FB==2 && Mfinal>Mbhmin && Mfinal<Mbhmax) {   /* for massive BHs with full fall back */
   delMsn1c=get_t(mgrbmin,mgrbmax);     /* mass loss range [Msun] from Ashall et al. 2019: arXiv:1702.04339 */ 
   Mfinal=Mfinal-delMsn1c;              /* apply addtional mass loss due to potential SNIc/GRB */
 }

 return Mfinal;
}


double compactf6(double M, double *Mco, double *Mhe)
{ /* DELAYED SN + moderate PPSN/PSN + adjustable BH neutrino emission + adjustable increase of CO size */
  /* returns gravitational mass of compact object [Msun]: Jun 13, 2019 */ 
 double a1,b1,c1,del1;
 double Mcc,Mex,Menv,Mfb,Mbar,Mgrav,Mfinal,frac;
 double Zcrit11,delMsn1c;
 
 Zcrit11=Zcrit1+0.00001;      /* ZZ<=Zcrit1+0.00001; to make sure that model with ZZ=Zcrit1 is affected by this: */

                              /* only for accretion above Fcrit1 and only for low-Z stars */
 if(Spin1==1 && Fms>Fcrit1 && ZZ<=Zcrit11) {   /* increase of CO core mass due to MS progenitor spin up from accretion */
   (*Mhe)=min(Fincr1*(*Mhe),M);                /* do not increase over entire star mass */
   (*Mco)=min(Fincr1*(*Mco),M);
 }

 if(PPSN==1 && (*Mhe)>=40.0 && (*Mhe)<65.00) {           /* PPSN regime: remove mass from the star */ 
   if((*Mhe)<60.00) Mcc=0.65*(*Mhe)+12.2;                 /* Mcc is mass of the colapsed helium core - pulse mass */
   else if((*Mhe)>=60.0 && (*Mhe)<62.5) Mcc=51.2;
   else Mcc=-14.3*(*Mhe)+938.0;
   M=(*Mhe)=Mcc;
   if((*Mco)>Mcc) (*Mco)=Mcc; 
 }
 else if(PSN==1 && (*Mhe)>=65.0 && (*Mhe)<135.0) {       /* PSN regime: explode the star and no remnant */
   Mfinal=M=0.0;
   return Mfinal;                                        /* exit function here */
 }

 if(((*Mco)-acc)>M)
   fprintf(fp0,"error: Mco>M in compactf6()\n");
 if((*Mco)<3.5) Mex=1.20;            /* Mex [Msun] -- baryon mass of the proto-NS at the time of bounce/explosion */
 else if((*Mco)<6.0) Mex=1.30; 
 else if((*Mco)<11.0) Mex=1.40;
 else Mex=1.60; 

 Menv=M-Mex;                      /* Menv [Msun] -- baryon mass of the "envelope": mass over proto-NS */

 if(Mex>M) {                      /* Mex>M: can happen for lowest mass stars K=8 M=2.0 with removed envelope: Mco=1.15 */
   Mex=M;                         /* remnant can not be more massive than the mass of exploding star */
   Menv=0.0;                      /* it is either a very low mass NS with almost no SN ejecta (as treated here) */
   Mfb=0.0;                       /* or it could potentially never explode and form a CO/ONeMg WD */
   frac=0.0;
 }
 else if((*Mco)<2.5) {
   Mfb=0.20;                     /* Mfb [Msun] -- baryon mass of the fall back material accreted onto proto-NS */ 
   Mfb=min(Mfb,Menv);            /* one can not have fall back larger than mass contained in the envelope */ 
   frac=Mfb/Menv;                /* frac [0-1] -- fallback parameter: fraction of envelope accreted onto proto-NS */
 } 
 else if((*Mco)<3.5) {              /* Mfb=0.2 (Mco=2.5Msun) -- 0.7 (Mco=3.5) Msun: linear interpolation */
   a1=0.5; 
   b1=-1.05; 
   Mfb=a1*(*Mco)+b1;
   Mfb=min(Mfb,Menv);
   frac=Mfb/Menv;
 }
 else if((*Mco)<11.0) {             /* frac is chosen such that at Mco=7.0 frac results in Mfb=0.7 Msun */
   a1=(1.0-(0.7/(M-Mex)))/7.5;   /* frac is chosen such that at Mco=11.0 frac results in Mfb=Menv (frac=1.0) */
   b1=1.0-11.0*a1;               /* and frac: linear interpolation for 3.5<Mco<11 */
   frac=a1*(*Mco)+b1; 
   Mfb=frac*Menv;
   Mfb=min(Mfb,Menv);
 } 
 else {                          /* full fall back: entire star goes into remnant */
   frac=1.0; 
   Mfb=Menv; 
 }

 if(frac<=0.001) FB=0;                    /* no fall back */
 else if(frac>0.001 && frac<0.999) FB=1;  /* partial fall back */
 else FB=2;                               /* full fall back */
 
 if(kickback==0 && Mfb<0.2) FB=0;                /* stdandard: full kick (fall back override) */
 else if(kickback==1 && Mfb<kickback1) FB=0;     /* parm. study: full kick (fall back override) */
 
 Mbar=Mex+Mfb;                   /* Mbar [Msun] -- final baryon mass of the remnant */
 if(Mbar>M+acc) 
   fprintf(fp0,"compact2f(): error -- something is wrong with the mass of the remnant\n");
 
 if(0.8*Mbar<Mmaxns) {           /* for NS progenitors */
   a1=0.075;                     /* for Mbar=3.0 (0.8*Mbar=2.4 Msun) -> Mgrav=2.523 Msun */
   b1=1.0;         
   c1=-(Mbar);        
   del1=b1*b1-4.0*a1*c1;
   Mgrav=(-b1+sqrt(del1))/(2.0*a1);   
   if(Mgrav>Mmaxns)                  /* light BH */
     Mgrav=(1.0-Fneut1)*Mbar;
 }
 else                                /* for BH progenitors */
   Mgrav=(1.0-Fneut2)*Mbar;  

 FRAC=frac;                           /* writes how much fallback [0:1] took place into global var. */
 Mfinal=min(M,Mgrav);                 /* to make sure that NS or BH mass do not exceed mass M of exploding star */

 if(GRB==1 && FB==2 && Mfinal>Mbhmin && Mfinal<Mbhmax) {   /* for massive BHs with full fall back  */
   delMsn1c=get_t(mgrbmin,mgrbmax);     /* mass loss range [Msun] from Ashall et al. 2019: arXiv:1702.04339 */ 
   Mfinal=Mfinal-delMsn1c;              /* apply addtional mass loss due to potential SNIc/GRB */
 }

 return Mfinal;
}


double compactf7(double M, double *Mco, double *Mhe)
{ /* RAPID SN + weak PPSN/PSN + adjustable BH neutrino emission + adjustable increase of CO size */
  /* returns gravitational mass of compact object [Msun]: Jun 13, 2019 */
 double a1,b1,c1,del1;
 double Mcc,Mex,Menv,Mfb,Mbar,Mgrav,Mfinal,frac;
 double Zcrit11,delMsn1c;
 
 Zcrit11=Zcrit1+0.00001;      /* ZZ<=Zcrit1+0.00001; to make sure that model with ZZ=Zcrit1 is affected by this: */
 
                                               /* only for accretion above Fcrit1 and only for low-Z stars */
 if(Spin1==1 && Fms>Fcrit1 && ZZ<=Zcrit11) {   /* increase of CO core mass due to MS progenitor spin up from accretion */
   (*Mhe)=min(Fincr1*(*Mhe),M);                /* do not increase over entire star mass */ 
   (*Mco)=min(Fincr1*(*Mco),M);                        
 }

 if(PPSN==1 && (*Mhe)>=40.0 && (*Mhe)<65.00) {           /* PPSN regime: remove mass from the star */ 
   if((*Mhe)<60.00) Mcc=0.83*(*Mhe)+6.0;                 /* Mcc is mass of the colapsed helium core - pulse mass */
   else if((*Mhe)>=60.0 && (*Mhe)<62.5) Mcc=55.6;
   else Mcc=-14.3*(*Mhe)+938.1;
   M=(*Mhe)=Mcc;
   if((*Mco)>Mcc) (*Mco)=Mcc; 
 }
 else if(PSN==1 && (*Mhe)>=65.0 && (*Mhe)<135.0) {       /* PSN regime: explode the star and no remnant */
   Mfinal=M=0.0;
   return Mfinal;                                        /* exit function here */
 }

 if(((*Mco)-acc)>M)
   fprintf(fp0,"error: Mco>M in compactf7()\n");
    
 Mex=1.0;            /* Mex [Msun] -- baryon mass of the proto-NS at the time of bounce/explosion */
                     /* the same Mex mass for all models */

 Menv=M-Mex;         /* Menv [Msun] -- baryon mass of the "envelope": mass over proto-NS */

 if(Mex>M) {         /* Mex>M: should not happen here: but look at comments at compactf2(): same place */
   Mex=M;            /* remnant can not be more massive than the mass of exploding star */
   Menv=0.0;         /* it is either a very low mass NS with almost no SN ejecta (as treated here) */
   Mfb=0.0;          /* or it could potentially never explode and form a CO/ONeMg WD */
   frac=0.0;
 }
 else if((*Mco)<2.5) {
   Mfb=0.20;                     /* Mfb [Msun] -- baryon mass of the fall back material accreted onto proto-NS */ 
   Mfb=min(Mfb,Menv);            /* one can not have fall back larger than mass contained in the envelope */ 
   frac=Mfb/Menv;                /* frac [0-1] -- fallback parameter: fraction of envelope accreted onto proto-NS */
 } 
 else if((*Mco)<6.0) {           /* Mfb=0.2Msun (Mco=2.5Msun) -- 1.2Msun (Mco=6.0): linear interpolation */
   a1=0.286; 
   b1=-0.514;; 
   Mfb=a1*(*Mco)+b1;
   Mfb=min(Mfb,Menv);
   frac=Mfb/Menv;
 }
 else if((*Mco)<7.0) {           /* Mfb=1.0 Mco=6.0--7.0 Msun: full fallback: unexplained feature in Heger&Woosley models */
   frac=1.0;                     /* but A.Heger claims that it is real and has something to do with structure of star in */ 
   Mfb=Menv;                     /* nparticular mass range */
 }
 else if((*Mco)<11.0) {          /* frac is chosen such that at Mco=7.0 frac results in Mfb=5.1 Msun */
   a1=(1.0-(5.1/(M-Mex)))/4.0;   /* frac is chosen such that at Mco=11.0 frac results in Mfb=Menv (frac=1.0) */
   b1=1.0-11.0*a1;               /* and frac: linear interpolation for 7.0<Mco<11 */
   frac=a1*(*Mco)+b1; 
   Mfb=frac*Menv;
   Mfb=min(Mfb,Menv);
 } 
 else {                          /* full fall back: entire star goes into remnant */
   frac=1.0; 
   Mfb=Menv; 
 }

 if(frac<=0.001) FB=0;                    /* no fall back */
 else if(frac>0.001 && frac<0.999) FB=1;  /* partial fall back */
 else FB=2;                               /* full fall back */
 
 if(kickback==0 && Mfb<0.2) FB=0;                /* stdandard: full kick (fall back override) */
 else if(kickback==1 && Mfb<kickback2) FB=0;     /* parm. study: full kick (fall back override) */
 
 Mbar=Mex+Mfb;                   /* Mbar [Msun] -- final baryon mass of the remnant */
 if(Mbar>M+acc) 
   fprintf(fp0,"compactf7(): error -- something is wrong with the mass of the remnant\n");
 
 if(0.8*Mbar<Mmaxns) {           /* for NS progenitors */
   a1=0.075;                     /* for Mbar=3.0 (0.8*Mbar=2.4 Msun) -> Mgrav=2.523 Msun */
   b1=1.0;         
   c1=-(Mbar);        
   del1=b1*b1-4.0*a1*c1;
   Mgrav=(-b1+sqrt(del1))/(2.0*a1);   
   if(Mgrav>Mmaxns)                  /* light BH */
     Mgrav=(1.0-Fneut1)*Mbar;
 }
 else                                /* for BH progenitors */
   Mgrav=(1.0-Fneut2)*Mbar;  

 FRAC=frac;                           /* writes how much fallback [0:1] took place into global var. */
 Mfinal=min(M,Mgrav);                 /* to make sure that NS or BH mass do not exceed mass M of exploding star */

 if(GRB==1 && FB==2 && Mfinal>Mbhmin && Mfinal<Mbhmax) {   /* for massive BHs with full fall back */
   delMsn1c=get_t(mgrbmin,mgrbmax);     /* mass loss range [Msun] from Ashall et al. 2019: arXiv:1702.04339 */ 
   Mfinal=Mfinal-delMsn1c;              /* apply addtional mass loss due to potential SNIc/GRB */
 }

 return Mfinal;
}


double compactf8(double M, double *Mco, double *Mhe)
{ /* DELAYED SN + weak PPSN/PSN + adjustable BH neutrino emission + adjustable increase of CO size */
  /* returns gravitational mass of compact object [Msun]: Jun 13, 2019 */
 double a1,b1,c1,del1;
 double Mcc,Mex,Menv,Mfb,Mbar,Mgrav,Mfinal,frac;
 double Zcrit11,delMsn1c;
 
 Zcrit11=Zcrit1+0.00001;      /* ZZ<=Zcrit1+0.00001; to make sure that model with ZZ=Zcrit1 is affected by this: */

                              /* only for accretion above Fcrit1 and only for low-Z stars */
 if(Spin1==1 && Fms>Fcrit1 && ZZ<=Zcrit11) {   /* increase of CO core mass due to MS progenitor spin up from accretion */
   (*Mhe)=min(Fincr1*(*Mhe),M);                /* do not increase over entire star mass */
   (*Mco)=min(Fincr1*(*Mco),M);
 }

 if(PPSN==1 && (*Mhe)>=40.0 && (*Mhe)<65.00) {           /* PPSN regime: remove mass from the star */ 
   if((*Mhe)<60.00) Mcc=0.83*(*Mhe)+6.0;                 /* Mcc is mass of the colapsed helium core - pulse mass */
   else if((*Mhe)>=60.0 && (*Mhe)<62.5) Mcc=55.6;
   else Mcc=-14.3*(*Mhe)+938.1;
   M=(*Mhe)=Mcc;
   if((*Mco)>Mcc) (*Mco)=Mcc; 
 }
 else if(PSN==1 && (*Mhe)>=65.0 && (*Mhe)<135.0) {       /* PSN regime: explode the star and no remnant */
   Mfinal=M=0.0;
   return Mfinal;                                        /* exit function here */
 }

 if(((*Mco)-acc)>M)
   fprintf(fp0,"error: Mco>M in compactf8()\n");
 if((*Mco)<3.5) Mex=1.20;            /* Mex [Msun] -- baryon mass of the proto-NS at the time of bounce/explosion */
 else if((*Mco)<6.0) Mex=1.30; 
 else if((*Mco)<11.0) Mex=1.40;
 else Mex=1.60; 

 Menv=M-Mex;                      /* Menv [Msun] -- baryon mass of the "envelope": mass over proto-NS */

 if(Mex>M) {                      /* Mex>M: can happen for lowest mass stars K=8 M=2.0 with removed envelope: Mco=1.15 */
   Mex=M;                         /* remnant can not be more massive than the mass of exploding star */
   Menv=0.0;                      /* it is either a very low mass NS with almost no SN ejecta (as treated here) */
   Mfb=0.0;                       /* or it could potentially never explode and form a CO/ONeMg WD */
   frac=0.0;
 }
 else if((*Mco)<2.5) {
   Mfb=0.20;                     /* Mfb [Msun] -- baryon mass of the fall back material accreted onto proto-NS */ 
   Mfb=min(Mfb,Menv);            /* one can not have fall back larger than mass contained in the envelope */ 
   frac=Mfb/Menv;                /* frac [0-1] -- fallback parameter: fraction of envelope accreted onto proto-NS */
 } 
 else if((*Mco)<3.5) {              /* Mfb=0.2 (Mco=2.5Msun) -- 0.7 (Mco=3.5) Msun: linear interpolation */
   a1=0.5; 
   b1=-1.05; 
   Mfb=a1*(*Mco)+b1;
   Mfb=min(Mfb,Menv);
   frac=Mfb/Menv;
 }
 else if((*Mco)<11.0) {             /* frac is chosen such that at Mco=7.0 frac results in Mfb=0.7 Msun */
   a1=(1.0-(0.7/(M-Mex)))/7.5;   /* frac is chosen such that at Mco=11.0 frac results in Mfb=Menv (frac=1.0) */
   b1=1.0-11.0*a1;               /* and frac: linear interpolation for 3.5<Mco<11 */
   frac=a1*(*Mco)+b1; 
   Mfb=frac*Menv;
   Mfb=min(Mfb,Menv);
 } 
 else {                          /* full fall back: entire star goes into remnant */
   frac=1.0; 
   Mfb=Menv; 
 }

 if(frac<=0.001) FB=0;                    /* no fall back */
 else if(frac>0.001 && frac<0.999) FB=1;  /* partial fall back */
 else FB=2;                               /* full fall back */
 
 if(kickback==0 && Mfb<0.2) FB=0;                /* stdandard: full kick (fall back override) */
 else if(kickback==1 && Mfb<kickback1) FB=0;     /* parm. study: full kick (fall back override) */
 
 Mbar=Mex+Mfb;                   /* Mbar [Msun] -- final baryon mass of the remnant */
 if(Mbar>M+acc) 
   fprintf(fp0,"compactf8(): error -- something is wrong with the mass of the remnant\n");
 
 if(0.8*Mbar<Mmaxns) {           /* for NS progenitors */
   a1=0.075;                     /* for Mbar=3.0 (0.8*Mbar=2.4 Msun) -> Mgrav=2.523 Msun */
   b1=1.0;         
   c1=-(Mbar);        
   del1=b1*b1-4.0*a1*c1;
   Mgrav=(-b1+sqrt(del1))/(2.0*a1);   
   if(Mgrav>Mmaxns)                  /* light BH */
     Mgrav=(1.0-Fneut1)*Mbar;
 }
 else                                /* for BH progenitors */
   Mgrav=(1.0-Fneut2)*Mbar;  

 FRAC=frac;                           /* writes how much fallback [0:1] took place into global var. */
 Mfinal=min(M,Mgrav);                 /* to make sure that NS or BH mass do not exceed mass M of exploding star */

 if(GRB==1 && FB==2 && Mfinal>Mbhmin && Mfinal<Mbhmax) {   /* for massive BHs with full fall back */
   delMsn1c=get_t(mgrbmin,mgrbmax);     /* mass loss range [Msun] from Ashall et al. 2019: arXiv:1702.04339 */ 
   Mfinal=Mfinal-delMsn1c;              /* apply addtional mass loss due to potential SNIc/GRB */
 }

 return Mfinal;
}


double compactf9(double M, double *Mco, double *Mhe)
{ /* RAPID SN + no PPSN + high PSN + adjustable BH neutrino emission + adjustable increase of CO size */
  /* returns gravitational mass of compact object [Msun]: Sep 17, 2020 */
 double a1,b1,c1,del1;
 double Mcc,Mex,Menv,Mfb,Mbar,Mgrav,Mfinal,frac;
 double Zcrit11,delMsn1c;
 
 Zcrit11=Zcrit1+0.00001;      /* ZZ<=Zcrit1+0.00001; to make sure that model with ZZ=Zcrit1 is affected by this: */
 
                                               /* only for accretion above Fcrit1 and only for low-Z stars */
 if(Spin1==1 && Fms>Fcrit1 && ZZ<=Zcrit11) {   /* increase of CO core mass due to MS progenitor spin up from accretion */
   (*Mhe)=min(Fincr1*(*Mhe),M);                /* do not increase over entire star mass */ 
   (*Mco)=min(Fincr1*(*Mco),M);                        
 }

                                           /* 40<Mhe<90Msun: PPSN regime: becomes CC SN regime: no need to do anything */ 
 if(PSN==1 && (*Mhe)>=65.0 && (*Mhe)<135.0) {       /* PSN regime: explode the star and no remnant */
   Mfinal=M=0.0;
   return Mfinal;                                        /* exit function here */
 }

 if(((*Mco)-acc)>M)
   fprintf(fp0,"error: Mco>M in compactf9()\n");
    
 Mex=1.0;            /* Mex [Msun] -- baryon mass of the proto-NS at the time of bounce/explosion */
                     /* the same Mex mass for all models */

 Menv=M-Mex;         /* Menv [Msun] -- baryon mass of the "envelope": mass over proto-NS */

 if(Mex>M) {         /* Mex>M: should not happen here: but look at comments at compactf2(): same place */
   Mex=M;            /* remnant can not be more massive than the mass of exploding star */
   Menv=0.0;         /* it is either a very low mass NS with almost no SN ejecta (as treated here) */
   Mfb=0.0;          /* or it could potentially never explode and form a CO/ONeMg WD */
   frac=0.0;
 }
 else if((*Mco)<2.5) {
   Mfb=0.20;                     /* Mfb [Msun] -- baryon mass of the fall back material accreted onto proto-NS */ 
   Mfb=min(Mfb,Menv);            /* one can not have fall back larger than mass contained in the envelope */ 
   frac=Mfb/Menv;                /* frac [0-1] -- fallback parameter: fraction of envelope accreted onto proto-NS */
 } 
 else if((*Mco)<6.0) {           /* Mfb=0.2Msun (Mco=2.5Msun) -- 1.2Msun (Mco=6.0): linear interpolation */
   a1=0.286; 
   b1=-0.514;; 
   Mfb=a1*(*Mco)+b1;
   Mfb=min(Mfb,Menv);
   frac=Mfb/Menv;
 }
 else if((*Mco)<7.0) {           /* Mfb=1.0 Mco=6.0--7.0 Msun: full fallback: unexplained feature in Heger&Woosley models */
   frac=1.0;                     /* but A.Heger claims that it is real and has something to do with structure of star in */ 
   Mfb=Menv;                     /* nparticular mass range */
 }
 else if((*Mco)<11.0) {          /* frac is chosen such that at Mco=7.0 frac results in Mfb=5.1 Msun */
   a1=(1.0-(5.1/(M-Mex)))/4.0;   /* frac is chosen such that at Mco=11.0 frac results in Mfb=Menv (frac=1.0) */
   b1=1.0-11.0*a1;               /* and frac: linear interpolation for 7.0<Mco<11 */
   frac=a1*(*Mco)+b1; 
   Mfb=frac*Menv;
   Mfb=min(Mfb,Menv);
 } 
 else {                          /* full fall back: entire star goes into remnant */
   frac=1.0; 
   Mfb=Menv; 
 }

 if(frac<=0.001) FB=0;                    /* no fall back */
 else if(frac>0.001 && frac<0.999) FB=1;  /* partial fall back */
 else FB=2;                               /* full fall back */
 
 if(kickback==0 && Mfb<0.2) FB=0;                /* stdandard: full kick (fall back override) */
 else if(kickback==1 && Mfb<kickback2) FB=0;     /* parm. study: full kick (fall back override) */
 
 Mbar=Mex+Mfb;                   /* Mbar [Msun] -- final baryon mass of the remnant */
 if(Mbar>M+acc) 
   fprintf(fp0,"compactf9(): error -- something is wrong with the mass of the remnant\n");
 
 if(0.8*Mbar<Mmaxns) {           /* for NS progenitors */
   a1=0.075;                     /* for Mbar=3.0 (0.8*Mbar=2.4 Msun) -> Mgrav=2.523 Msun */
   b1=1.0;         
   c1=-(Mbar);        
   del1=b1*b1-4.0*a1*c1;
   Mgrav=(-b1+sqrt(del1))/(2.0*a1);   
   if(Mgrav>Mmaxns)                  /* light BH */
     Mgrav=(1.0-Fneut1)*Mbar;
 }
 else                                /* for BH progenitors */
   Mgrav=(1.0-Fneut2)*Mbar;  

 FRAC=frac;                           /* writes how much fallback [0:1] took place into global var. */
 Mfinal=min(M,Mgrav);                 /* to make sure that NS or BH mass do not exceed mass M of exploding star */

 if(GRB==1 && FB==2 && Mfinal>Mbhmin && Mfinal<Mbhmax) {   /* for massive BHs with full fall back */
   delMsn1c=get_t(mgrbmin,mgrbmax);     /* mass loss range [Msun] from Ashall et al. 2019: arXiv:1702.04339 */ 
   Mfinal=Mfinal-delMsn1c;              /* apply addtional mass loss due to potential SNIc/GRB */
 }

 return Mfinal;
}


double compactf10(double M, double *Mco, double *Mhe)
{ /* DELAYED SN + no PPSN + high PSN + adjustable BH neutrino emission + adjustable increase of CO size */
  /* returns gravitational mass of compact object [Msun]: Sep 17, 2020 */
 double a1,b1,c1,del1;
 double Mcc,Mex,Menv,Mfb,Mbar,Mgrav,Mfinal,frac;
 double Zcrit11,delMsn1c;
 double a3,b3,frac35;
 
 Zcrit11=Zcrit1+0.00001;      /* ZZ<=Zcrit1+0.00001; to make sure that model with ZZ=Zcrit1 is affected by this: */

                              /* only for accretion above Fcrit1 and only for low-Z stars */
 if(Spin1==1 && Fms>Fcrit1 && ZZ<=Zcrit11) {   /* increase of CO core mass due to MS progenitor spin up from accretion */
   (*Mhe)=min(Fincr1*(*Mhe),M);                /* do not increase over entire star mass */
   (*Mco)=min(Fincr1*(*Mco),M);
 }

                                           /* 40<Mhe<90Msun: PPSN regime: becomes CC SN regime: no need to do anything */
 if(PSN==1 && (*Mhe)>=90.0 && (*Mhe)<175.0) {       /* PSN regime: explode the star and no remnant */
   Mfinal=M=0.0;
   return Mfinal;                                        /* exit function here */
 }

 if(((*Mco)-acc)>M)
   fprintf(fp0,"error: Mco>M in compactf10()\n");
 if((*Mco)<3.5) Mex=1.20;            /* Mex [Msun] -- baryon mass of the proto-NS at the time of bounce/explosion */
 else if((*Mco)<6.0) Mex=1.30; 
 else if((*Mco)<11.0) Mex=1.40;
 else Mex=1.60; 

 Menv=M-Mex;                      /* Menv [Msun] -- baryon mass of the "envelope": mass over proto-NS */

 if(Mex>M) {                      /* Mex>M: can happen for lowest mass stars K=8 M=2.0 with removed envelope: Mco=1.15 */
   Mex=M;                         /* remnant can not be more massive than the mass of exploding star */
   Menv=0.0;                      /* it is either a very low mass NS with almost no SN ejecta (as treated here) */
   Mfb=0.0;                       /* or it could potentially never explode and form a CO/ONeMg WD */
   frac=0.0;
 }
 else if((*Mco)<2.5) {
   Mfb=0.20;                     /* Mfb [Msun] -- baryon mass of the fall back material accreted onto proto-NS */ 
   Mfb=min(Mfb,Menv);            /* one can not have fall back larger than mass contained in the envelope */ 
   frac=Mfb/Menv;                /* frac [0-1] -- fallback parameter: fraction of envelope accreted onto proto-NS */
 } 
 else if((*Mco)<3.5) {              /* Mfb=0.2 (Mco=2.5Msun) -- 0.7 (Mco=3.5) Msun: linear interpolation */
   a1=0.5; 
   b1=-1.05; 
   Mfb=a1*(*Mco)+b1;
   Mfb=min(Mfb,Menv);
   frac=Mfb/Menv;
 }
 else if((*Mco)<McoExt) {        /* frac is chosen such that at Mco=7.0 frac results in Mfb=0.7 Msun */
   a1=(1.0-(0.7/(M-Mex)))/7.5;   /* frac is chosen such that at Mco=McoExt frac results in Mfb=Menv (frac=1.0) */
   b1=1.0-11.0*a1;               /* and frac: linear interpolation for 3.5<Mco<McoExt */
   frac35=a1*3.5+b1;             /* frac at Mco=3.5Msun; same for all relations no matter how extended */ 
   a3=(frac35-1.0)/(3.5-McoExt);
   b3=1.0-a3*McoExt;
   frac=a3*(*Mco)+b3; 
   Mfb=frac*Menv;
   Mfb=min(Mfb,Menv);
 } 
 else {                          /* full fall back: entire star goes into remnant */
   frac=1.0; 
   Mfb=Menv; 
 }

 if(frac<=0.001) FB=0;                    /* no fall back */
 else if(frac>0.001 && frac<0.999) FB=1;  /* partial fall back */
 else FB=2;                               /* full fall back */
 
 if(kickback==0 && Mfb<0.2) FB=0;                /* stdandard: full kick (fall back override) */
 else if(kickback==1 && Mfb<kickback1) FB=0;     /* parm. study: full kick (fall back override) */
 
 Mbar=Mex+Mfb;                   /* Mbar [Msun] -- final baryon mass of the remnant */
 if(Mbar>M+acc) 
   fprintf(fp0,"compactf10(): error -- something is wrong with the mass of the remnant\n");
 
 if(0.8*Mbar<Mmaxns) {           /* for NS progenitors */
   a1=0.075;                     /* for Mbar=3.0 (0.8*Mbar=2.4 Msun) -> Mgrav=2.523 Msun */
   b1=1.0;         
   c1=-(Mbar);        
   del1=b1*b1-4.0*a1*c1;
   Mgrav=(-b1+sqrt(del1))/(2.0*a1);   
   if(Mgrav>Mmaxns)                  /* light BH */
     Mgrav=(1.0-Fneut1)*Mbar;
 }
 else                                /* for BH progenitors */
   Mgrav=(1.0-Fneut2)*Mbar;  

 FRAC=frac;                           /* writes how much fallback [0:1] took place into global var. */
 Mfinal=min(M,Mgrav);                 /* to make sure that NS or BH mass do not exceed mass M of exploding star */

 if(GRB==1 && FB==2 && Mfinal>Mbhmin && Mfinal<Mbhmax) {   /* for massive BHs with full fall back */
   delMsn1c=get_t(mgrbmin,mgrbmax);     /* mass loss range [Msun] from Ashall et al. 2019: arXiv:1702.04339 */ 
   Mfinal=Mfinal-delMsn1c;              /* apply addtional mass loss due to potential SNIc/GRB */
 }

 return Mfinal;
}


double compactf12(double M, double *Mco, double *Mhe)
{ /* DELAYED SN + weak PPSN/PSN + adjustable BH neutrino emission + adjustable increase of CO size */
  /* and extended partial fall back regime, from 11Msun to McoExt: sinbin.h */
  /* returns gravitational mass of compact object [Msun]: Mar 06, 2021 */
 double a1,b1,c1,del1;
 double Mcc,Mex,Menv,Mfb,Mbar,Mgrav,Mfinal,frac;
 double Zcrit11,delMsn1c;
 double a3,b3,frac35;
 
 Zcrit11=Zcrit1+0.00001;      /* ZZ<=Zcrit1+0.00001; to make sure that model with ZZ=Zcrit1 is affected by this: */

                              /* only for accretion above Fcrit1 and only for low-Z stars */
 if(Spin1==1 && Fms>Fcrit1 && ZZ<=Zcrit11) {   /* increase of CO core mass due to MS progenitor spin up from accretion */
   (*Mhe)=min(Fincr1*(*Mhe),M);                /* do not increase over entire star mass */
   (*Mco)=min(Fincr1*(*Mco),M);
 }

 if(PPSN==1 && (*Mhe)>=40.0 && (*Mhe)<65.00) {           /* PPSN regime: remove mass from the star */ 
   if((*Mhe)<60.00) Mcc=0.83*(*Mhe)+6.0;                 /* Mcc is mass of the colapsed helium core - pulse mass */
   else if((*Mhe)>=60.0 && (*Mhe)<62.5) Mcc=55.6;
   else Mcc=-14.3*(*Mhe)+938.1;
   M=(*Mhe)=Mcc;
   if((*Mco)>Mcc) (*Mco)=Mcc; 
 }
 else if(PSN==1 && (*Mhe)>=65.0 && (*Mhe)<135.0) {       /* PSN regime: explode the star and no remnant */
   Mfinal=M=0.0;
   return Mfinal;                                        /* exit function here */
 }

 if(((*Mco)-acc)>M)
   fprintf(fp0,"error: Mco>M in compactf8()\n");
 if((*Mco)<3.5) Mex=1.20;            /* Mex [Msun] -- baryon mass of the proto-NS at the time of bounce/explosion */
 else if((*Mco)<6.0) Mex=1.30; 
 else if((*Mco)<11.0) Mex=1.40;
 else Mex=1.60; 

 Menv=M-Mex;                     /* Menv [Msun] -- baryon mass of the "envelope": mass over proto-NS */

 if(Mex>M) {                     /* Mex>M: can happen for lowest mass stars K=8 M=2.0 with removed envelope: Mco=1.15 */
   Mex=M;                        /* remnant can not be more massive than the mass of exploding star */
   Menv=0.0;                     /* it is either a very low mass NS with almost no SN ejecta (as treated here) */
   Mfb=0.0;                      /* or it could potentially never explode and form a CO/ONeMg WD */
   frac=0.0;
 }
 else if((*Mco)<2.5) {
   Mfb=0.20;                     /* Mfb [Msun] -- baryon mass of the fall back material accreted onto proto-NS */ 
   Mfb=min(Mfb,Menv);            /* one can not have fall back larger than mass contained in the envelope */ 
   frac=Mfb/Menv;                /* frac [0-1] -- fallback parameter: fraction of envelope accreted onto proto-NS */
 } 
 else if((*Mco)<3.5) {           /* Mfb=0.2 (Mco=2.5Msun) -- 0.7 (Mco=3.5) Msun: linear interpolation */
   a1=0.5; 
   b1=-1.05; 
   Mfb=a1*(*Mco)+b1;
   Mfb=min(Mfb,Menv);
   frac=Mfb/Menv;
 }
 else if((*Mco)<McoExt) {        /* frac is chosen such that at Mco=7.0 frac results in Mfb=0.7 Msun */
   a1=(1.0-(0.7/(M-Mex)))/7.5;   /* frac is chosen such that at Mco=McoExt frac results in Mfb=Menv (frac=1.0) */
   b1=1.0-11.0*a1;               /* and frac: linear interpolation for 3.5<Mco<McoExt */
   frac35=a1*3.5+b1;             /* frac at Mco=3.5Msun; same for all relations no matter how extended */ 
   a3=(frac35-1.0)/(3.5-McoExt);
   b3=1.0-a3*McoExt;
   frac=a3*(*Mco)+b3; 
   Mfb=frac*Menv;
   Mfb=min(Mfb,Menv);
 } 
 else {                          /* full fall back: entire star goes into remnant */
   frac=1.0; 
   Mfb=Menv; 
 }

 if(frac<=0.001) FB=0;                    /* no fall back */
 else if(frac>0.001 && frac<0.999) FB=1;  /* partial fall back */
 else FB=2;                               /* full fall back */
 
 if(kickback==0 && Mfb<0.2) FB=0;                /* stdandard: full kick (fall back override) */
 else if(kickback==1 && Mfb<kickback1) FB=0;     /* parm. study: full kick (fall back override) */
 
 Mbar=Mex+Mfb;                   /* Mbar [Msun] -- final baryon mass of the remnant */
 if(Mbar>M+acc) 
   fprintf(fp0,"compactf8(): error -- something is wrong with the mass of the remnant\n");
 
 if(0.8*Mbar<Mmaxns) {           /* for NS progenitors */
   a1=0.075;                     /* for Mbar=3.0 (0.8*Mbar=2.4 Msun) -> Mgrav=2.523 Msun */
   b1=1.0;         
   c1=-(Mbar);        
   del1=b1*b1-4.0*a1*c1;
   Mgrav=(-b1+sqrt(del1))/(2.0*a1);   
   if(Mgrav>Mmaxns)                  /* light BH */
     Mgrav=(1.0-Fneut1)*Mbar;
 }
 else                                /* for BH progenitors */
   Mgrav=(1.0-Fneut2)*Mbar;  

 FRAC=frac;                           /* writes how much fallback [0:1] took place into global var. */
 Mfinal=min(M,Mgrav);                 /* to make sure that NS or BH mass do not exceed mass M of exploding star */

 if(GRB==1 && FB==2 && Mfinal>Mbhmin && Mfinal<Mbhmax) {   /* for massive BHs with full fall back */
   delMsn1c=get_t(mgrbmin,mgrbmax);     /* mass loss range [Msun] from Ashall et al. 2019: arXiv:1702.04339 */ 
   Mfinal=Mfinal-delMsn1c;              /* apply addtional mass loss due to potential SNIc/GRB */
 }

 return Mfinal;
}


double compactf13(double M, double *Mco, double *Mhe)
{ /* new REMNANT mass function: Fryer/Olejak 2022 */
  /* DELAYED SN + weak PPSN/PSN + adjustable BH neutrino emission + adjustable increase of CO size */
  /* and extended partial fall back regime, from 11Msun to McoExt: sinbin.h */
  /* returns gravitational mass of compact object [Msun]: Mar 06, 2021 */
 double a1,b1,c1,del1;
 double Mcc,Mex,Menv,Mfb,Mbar,Mgrav,Mfinal,frac;
 double Zcrit11,delMsn1c;
 double a3,b3,frac35;
 
 Zcrit11=Zcrit1+0.00001;      /* ZZ<=Zcrit1+0.00001; to make sure that model with ZZ=Zcrit1 is affected by this: */

                              /* only for accretion above Fcrit1 and only for low-Z stars */
 if(Spin1==1 && Fms>Fcrit1 && ZZ<=Zcrit11) {   /* increase of CO core mass due to MS progenitor spin up from accretion */
   (*Mhe)=min(Fincr1*(*Mhe),M);                /* do not increase over entire star mass */
   (*Mco)=min(Fincr1*(*Mco),M);
 }

 if(XXPSN==1) {
   if(PPSN==1 && (*Mhe)>=45.0 && (*Mhe)<65.00) {            /* STRONG PPSN regime: remove mass from the star */ 
     M=(*Mhe)=45.0;                                         /* before calculating mass of the remnant */ 
     if((*Mco)>45.0) (*Mco)=45.0;
   }
   else if(PSN==1 && (*Mhe)>=65.0 && (*Mhe)<135.0) {        /* NORMAL PSN regime: explode the star and no remnant */
     Mfinal=M=0.0;
     return Mfinal;                                         /* exit function here */
   }
 }
 else if(XXPSN==2) {
   if(PPSN==1 && (*Mhe)>=40.0 && (*Mhe)<65.00) {           /* MODERATE PPSN regime: remove mass from the star */ 
     if((*Mhe)<60.00) Mcc=0.65*(*Mhe)+12.2;                /* Mcc is mass of the colapsed helium core - pulse mass */
     else if((*Mhe)>=60.0 && (*Mhe)<62.5) Mcc=51.2;
     else Mcc=-14.3*(*Mhe)+938.0;
     M=(*Mhe)=Mcc;
     if((*Mco)>Mcc) (*Mco)=Mcc; 
   }
   else if(PSN==1 && (*Mhe)>=65.0 && (*Mhe)<135.0) {       /* NORMAL PSN regime: explode the star and no remnant */
     Mfinal=M=0.0;
     return Mfinal;                                        /* exit function here */
   }
 } 
 else if(XXPSN==3) {
   if(PPSN==1 && (*Mhe)>=40.0 && (*Mhe)<65.00) {           /* WEAK PPSN regime: remove mass from the star */ 
     if((*Mhe)<60.00) Mcc=0.83*(*Mhe)+6.0;                 /* Mcc is mass of the colapsed helium core - pulse mass */
     else if((*Mhe)>=60.0 && (*Mhe)<62.5) Mcc=55.6;
     else Mcc=-14.3*(*Mhe)+938.1;
     M=(*Mhe)=Mcc;
     if((*Mco)>Mcc) (*Mco)=Mcc; 
   }
   else if(PSN==1 && (*Mhe)>=65.0 && (*Mhe)<135.0) {       /* NORMAL PSN regime: explode the star and no remnant */
     Mfinal=M=0.0;
     return Mfinal;                                        /* exit function here */
   }
 }
 else if(XXPSN==4) {                                       /* NO PPSN regime + ULTRA HIGH PSN regime */ 
   if(PSN==1 && (*Mhe)>=90.0 && (*Mhe)<175.0) {            /* 40<Mhe<90Msun: normal PPSN regime: becomes CC SN: no need to do anything */ 
     Mfinal=M=0.0;                                         /* PSN regime (40<Mhe<90Msun): explode the star and no remnant */
     return Mfinal;                                        /* exit function here */
   }
 }
 else if(XXPSN==5) {                                       /* adjustable PPSN + adjustable PSN: stars in [Mppsn2:Mpsn1] are CC SNe */  
   if(PPSN==1 && (*Mhe)>=Mppsn1 && (*Mhe)<Mppsn2) {        /* PPSN regime [Mppsn1:Mppsn2]: remove mass from the star down to Mppsn1 */ 
     M=(*Mhe)=Mppsn1;                                      /* before calculating mass of the remnant */ 
     if((*Mco)>Mppsn1) (*Mco)=Mppsn1;                      /* for example Mpsn1=45.0, Mpsn2=65.0 for strong PPSN */
   }
   if(PSN==1 && (*Mhe)>=Mpsn1 && (*Mhe)<Mpsn2) {           /* PSN regime [Mpsn1:Mpsn2], for example Mpsn1=65.0, Mpsn2=135.0 for norml PSN */
     Mfinal=M=0.0;                                         /* PSN regime: explode the star and no remnant */
     return Mfinal;                                        /* exit function here */
   }                                                       /* stars in [Mppsn2:Mpsn1] are CC SNe, so make Mppsn2=Mpsn1 for normal flow */
 }                                                         /* lower mass stars are subject to PPSN and then smoothly transition to PSN */ 
 else 
   fprintf(fp0,"error: unknown choice on PPSN/PSN in compactf13() XXPSN: %d\n",XXPSN);


 if(((*Mco)-acc)>M)
   fprintf(fp0,"error: Mco>M in compactf13()\n");
 if((*Mco)<3.5) Mex=1.20;            /* Mex [Msun] -- baryon mass of the proto-NS at the time of bounce/explosion */
 else if((*Mco)<6.0) Mex=1.30; 
 else if((*Mco)<11.0) Mex=1.40;
 else Mex=1.60; 

 Menv=M-Mex;                     /* Menv [Msun] -- baryon mass of the "envelope": mass over proto-NS */

 if(Mex>M) {                     /* Mex>M: can happen for lowest mass stars K=8 M=2.0 with removed envelope: Mco=1.15 */
   Mex=M;                        /* remnant can not be more massive than the mass of exploding star */
   Menv=0.0;                     /* it is either a very low mass NS with almost no SN ejecta (as treated here) */
   Mfb=0.0;                      /* or it could potentially never explode and form a CO/ONeMg WD */
   frac=0.0;
 }
 else if((*Mco)<2.5) {
   Mfb=0.20;                     /* Mfb [Msun] -- baryon mass of the fall back material accreted onto proto-NS */ 
   Mfb=min(Mfb,Menv);            /* one can not have fall back larger than mass contained in the envelope */ 
   frac=Mfb/Menv;                /* frac [0-1] -- fallback parameter: fraction of envelope accreted onto proto-NS */
 } 
 else if((*Mco)<3.5) {           /* Mfb=0.2 (Mco=2.5Msun) -- 0.7 (Mco=3.5) Msun: linear interpolation */
   a1=0.5; 
   b1=-1.05; 
   Mfb=a1*(*Mco)+b1;
   Mfb=min(Mfb,Menv);
   frac=Mfb/Menv;
 }
 else if((*Mco)<McoExt) {        /* frac is chosen such that at Mco=7.0 frac results in Mfb=0.7 Msun */
   a1=(1.0-(0.7/(M-Mex)))/7.5;   /* frac is chosen such that at Mco=McoExt frac results in Mfb=Menv (frac=1.0) */
   b1=1.0-11.0*a1;               /* and frac: linear interpolation for 3.5<Mco<McoExt */
   frac35=a1*3.5+b1;             /* frac at Mco=3.5Msun; same for all relations no matter how extended */ 
   a3=(frac35-1.0)/(3.5-McoExt);
   b3=1.0-a3*McoExt;
   frac=a3*(*Mco)+b3; 
   Mfb=frac*Menv;
   Mfb=min(Mfb,Menv);
 } 
 else {                          /* full fall back: entire star goes into remnant */
   frac=1.0; 
   Mfb=Menv; 
 }

 if(frac<=0.001) FB=0;                    /* no fall back */
 else if(frac>0.001 && frac<0.999) FB=1;  /* partial fall back */
 else FB=2;                               /* full fall back */
 
 if(kickback==0 && Mfb<0.2) FB=0;                /* stdandard: full kick (fall back override) */
 else if(kickback==1 && Mfb<kickback1) FB=0;     /* parm. study: full kick (fall back override) */
 
                                                                             /* Mbar [Msun]: final baryon mass of remnant */
 Mbar=1.2+0.05*fmix+0.01*pow((*Mco/fmix),2.0)+exp(fmix*(*Mco-McritNSBH));    /* remnant mass from Fryer et al. 2022 */
 Mbar=min(M,Mbar);             /* make sure that remnant mass is not larger than collapsing star mass  */
 
 if(0.8*Mbar<Mmaxns) {           /* for NS progenitors */
   a1=0.075;                     /* for Mbar=3.0 (0.8*Mbar=2.4 Msun) -> Mgrav=2.523 Msun */
   b1=1.0;         
   c1=-(Mbar);        
   del1=b1*b1-4.0*a1*c1;
   Mgrav=(-b1+sqrt(del1))/(2.0*a1);   
   if(Mgrav>Mmaxns)                  /* light BH */
     Mgrav=(1.0-Fneut1)*Mbar;
 }
 else                                /* for BH progenitors */
   Mgrav=(1.0-Fneut2)*Mbar;  

 FRAC=frac;                           /* writes how much fallback [0:1] took place into global var. */
 Mfinal=min(M,Mgrav);                 /* to make sure that NS or BH mass do not exceed mass M of exploding star */

 if(GRB==1 && FB==2 && Mfinal>Mbhmin && Mfinal<Mbhmax) {   /* for massive BHs with full fall back */
   delMsn1c=get_t(mgrbmin,mgrbmax);     /* mass loss range [Msun] from Ashall et al. 2019: arXiv:1702.04339 */ 
   Mfinal=Mfinal-delMsn1c;              /* apply addtional mass loss due to potential SNIc/GRB */
 }

 return Mfinal;
}


double compactf14(double M, double *Mco, double *Mhe)
{ /* new REMNANT mass function: Fryer/Olejak 2022 as updated by Chris B.: fallback treatment */
  /* DELAYED SN + weak PPSN/PSN + adjustable BH neutrino emission + adjustable increase of CO size */
  /* and extended partial fall back regime, from 11Msun to McoExt: sinbin.h */
  /* returns gravitational mass of compact object [Msun]: Mar 06, 2021 */
 double a1,b1,c1,del1;
 double Mcc,Mex,Menv,Mfb,Mbar,Mgrav,Mfinal,frac;
 double Zcrit11,delMsn1c;
 double a3,b3,frac35;
 
 Zcrit11=Zcrit1+0.00001;      /* ZZ<=Zcrit1+0.00001; to make sure that model with ZZ=Zcrit1 is affected by this: */

                              /* only for accretion above Fcrit1 and only for low-Z stars */
 if(Spin1==1 && Fms>Fcrit1 && ZZ<=Zcrit11) {   /* increase of CO core mass due to MS progenitor spin up from accretion */
   (*Mhe)=min(Fincr1*(*Mhe),M);                /* do not increase over entire star mass */
   (*Mco)=min(Fincr1*(*Mco),M);
 }

 if(XXPSN==1) {
   if(PPSN==1 && (*Mhe)>=45.0 && (*Mhe)<65.00) {            /* STRONG PPSN regime: remove mass from the star */ 
     M=(*Mhe)=45.0;                                         /* before calculating mass of the remnant */ 
     if((*Mco)>45.0) (*Mco)=45.0;
   }
   else if(PSN==1 && (*Mhe)>=65.0 && (*Mhe)<135.0) {        /* NORMAL PSN regime: explode the star and no remnant */
     Mfinal=M=0.0;
     return Mfinal;                                         /* exit function here */
   }
 }
 else if(XXPSN==2) {
   if(PPSN==1 && (*Mhe)>=40.0 && (*Mhe)<65.00) {           /* MODERATE PPSN regime: remove mass from the star */ 
     if((*Mhe)<60.00) Mcc=0.65*(*Mhe)+12.2;                /* Mcc is mass of the colapsed helium core - pulse mass */
     else if((*Mhe)>=60.0 && (*Mhe)<62.5) Mcc=51.2;
     else Mcc=-14.3*(*Mhe)+938.0;
     M=(*Mhe)=Mcc;
     if((*Mco)>Mcc) (*Mco)=Mcc; 
   }
   else if(PSN==1 && (*Mhe)>=65.0 && (*Mhe)<135.0) {       /* NORMAL PSN regime: explode the star and no remnant */
     Mfinal=M=0.0;
     return Mfinal;                                        /* exit function here */
   }
 } 
 else if(XXPSN==3) {
   if(PPSN==1 && (*Mhe)>=40.0 && (*Mhe)<65.00) {           /* WEAK PPSN regime: remove mass from the star */ 
     if((*Mhe)<60.00) Mcc=0.83*(*Mhe)+6.0;                 /* Mcc is mass of the colapsed helium core - pulse mass */
     else if((*Mhe)>=60.0 && (*Mhe)<62.5) Mcc=55.6;
     else Mcc=-14.3*(*Mhe)+938.1;
     M=(*Mhe)=Mcc;
     if((*Mco)>Mcc) (*Mco)=Mcc; 
   }
   else if(PSN==1 && (*Mhe)>=65.0 && (*Mhe)<135.0) {       /* NORMAL PSN regime: explode the star and no remnant */
     Mfinal=M=0.0;
     return Mfinal;                                        /* exit function here */
   }
 }
 else if(XXPSN==4) {                                       /* NO PPSN regime + ULTRA HIGH PSN regime */ 
   if(PSN==1 && (*Mhe)>=90.0 && (*Mhe)<175.0) {            /* 40<Mhe<90Msun: normal PPSN regime: becomes CC SN: no need to do anything */ 
     Mfinal=M=0.0;                                         /* PSN regime (40<Mhe<90Msun): explode the star and no remnant */
     return Mfinal;                                        /* exit function here */
   }
 }
 else if(XXPSN==5) {                                       /* adjustable PPSN + adjustable PSN: stars in [Mppsn2:Mpsn1] are CC SNe */  
   if(PPSN==1 && (*Mhe)>=Mppsn1 && (*Mhe)<Mppsn2) {        /* PPSN regime [Mppsn1:Mppsn2]: remove mass from the star down to Mppsn1 */ 
     M=(*Mhe)=Mppsn1;                                      /* before calculating mass of the remnant */ 
     if((*Mco)>Mppsn1) (*Mco)=Mppsn1;                      /* for example Mpsn1=45.0, Mpsn2=65.0 for strong PPSN */
   }
   if(PSN==1 && (*Mhe)>=Mpsn1 && (*Mhe)<Mpsn2) {           /* PSN regime [Mpsn1:Mpsn2], for example Mpsn1=65.0, Mpsn2=135.0 for norml PSN */
     Mfinal=M=0.0;                                         /* PSN regime: explode the star and no remnant */
     return Mfinal;                                        /* exit function here */
   }                                                       /* stars in [Mppsn2:Mpsn1] are CC SNe, so make Mppsn2=Mpsn1 for normal flow */
 }                                                         /* lower mass stars are subject to PPSN and then smoothly transition to PSN */ 
 else 
   fprintf(fp0,"error: unknown choice on PPSN/PSN in compactf14() XXPSN: %d\n",XXPSN);


 if(((*Mco)-acc)>M)
   fprintf(fp0,"error: Mco>M in compactf14()\n");
   
   
 if((*Mco)<3.5) Mex=1.20;            /* Mex [Msun] -- baryon mass of the proto-NS at the time of bounce/explosion */
 else if((*Mco)<6.0) Mex=1.30; 
 else if((*Mco)<11.0) Mex=1.40;
 else Mex=1.60; 

                                                                             /* Mbar [Msun]: final baryon mass of remnant */
 Mbar=1.2+0.05*fmix+0.01*pow((*Mco/fmix),2.0)+exp(fmix*(*Mco-McritNSBH));    /* remnant mass from Fryer et al. 2022 */
 Mbar=min(M,Mbar);                                /* make sure that remnant mass is not larger than collapsing star mass  */

 Menv=M-Mex;                       /* [Msun] mass of star that is above proto-NS: so called "envelope mass" */ 
 Menv=max(0.0,Menv);               /* to make sure that envelope mass is not smaller than zero */

 Mfb=Mbar-Mex;                     /* Mfb [Msun] -- baryon mass of the fall back material accreted onto proto-NS */ 
 Mfb=max(0.0,Mfb);                 /* to make sure it is not negative value */
 Mfb=min(Mfb,Menv);                /* one can not have fall back larger than mass contained in the envelope */ 

 if(Menv<1.0e-06)         
   frac=0.0;
 else
   frac=Mfb/Menv;                  /* frac [0-1] -- fallback parameter: fraction of envelope accreted onto proto-NS */
 

 if(frac<=0.001) FB=0;                    /* no fall back */
 else if(frac>0.001 && frac<0.999) FB=1;  /* partial fall back */
 else FB=2;                               /* full fall back */
 
 if(kickback==0 && Mfb<0.2) FB=0;                /* stdandard: full kick (fall back override) */
 else if(kickback==1 && Mfb<kickback1) FB=0;     /* parm. study: full kick (fall back override) */

 if(0.8*Mbar<Mmaxns) {           /* for NS progenitors */
   a1=0.075;                     /* for Mbar=3.0 (0.8*Mbar=2.4 Msun) -> Mgrav=2.523 Msun */
   b1=1.0;         
   c1=-(Mbar);        
   del1=b1*b1-4.0*a1*c1;
   Mgrav=(-b1+sqrt(del1))/(2.0*a1);   
   if(Mgrav>Mmaxns)                  /* light BH */
     Mgrav=(1.0-Fneut1)*Mbar;
 }
 else                                /* for BH progenitors */
   Mgrav=(1.0-Fneut2)*Mbar;  

 FRAC=frac;                           /* writes how much fallback [0:1] took place into global var. */
 Mfinal=min(M,Mgrav);                 /* to make sure that NS or BH mass do not exceed mass M of exploding star */

 if(GRB==1 && FB==2 && Mfinal>Mbhmin && Mfinal<Mbhmax) {   /* for massive BHs with full fall back */
   delMsn1c=get_t(mgrbmin,mgrbmax);     /* mass loss range [Msun] from Ashall et al. 2019: arXiv:1702.04339 */ 
   Mfinal=Mfinal-delMsn1c;              /* apply addtional mass loss due to potential SNIc/GRB */
 }

 return Mfinal;
}


double Lwdf(double M, double t, int K)
{ /* calculate luminosity[Lsun] of WD of M and K at time t[Myr] since its formation, eq.(90) with change: */
  /* A for: K=10: He WD A=4.0, for K=11: CO WD and K=12: ONe WD A=16.0 -- this value(16) is average of */
  /* what is given in paper0 for CO WD A=15 and ONe WD A=17, but in hridag.f l.38 only one A=16 is used */
  /* acctulally the A^1.4 is given: aco=48.503d0 => A=16 */
 double A;       

 if(K==16) A=1.0;
 else if(K==10 || K==17) A=4.0;
 else A=16.0;

 return 635.0*M*pow(ZZ,0.4)/pow(A*(t+0.1),1.4); 
} 


double Rwdf(double M, int K)
{ /* calculates radius[Rsun] for WD of mass M, eq.(91) */
 double const1,X;
 
 if(K==16) {
   X=0.7;                       /* composition parameter: H-content */
   const1=pow(1.0+X,5.0/3.0);   /* correction to calculate radius of H-rich WD K=16 */
 }
 else
   const1=1.0;                  /* for He and CO WDs, X=0.0 */                  
 
 return max(Rnsf(),const1*0.0115*sqrt(pow(MCh/M,2.0/3.0)-pow(M/MCh,2.0/3.0)));
}


double Rnsf(void)
{ /* returns radius [Rsun] of neutron star -- value is about 10km */

 return 1.4e-05;
}


double Rnsf1(double M)
{ /* returns radius [Rsun] of neutron star with mass M: Fortin et al. 2016, Phys.Rev.C 94, 035804 */
 double Rns;


 if(M>2.17) 
   Rns=10.36;                                                /* [km] */
 else if(M<1.39)
   Rns=11.78;                                                /* [km] */
 else  
   Rns=47.4488-64.7709*M+39.1216*M*M-7.90307*M*M*M;          /* [km] */
 
 Rns/=(Rsun*1.0e-05);                                        /* [Rsun] */

 return Rns;
}


double Rbhf(double M)
{ /* returns Schwarzschild radius [Rsun] of black hole */
 
 return 4.24e-06*M;
}


double Riscof(double M, double aspin) 
{ /* returns radius [Rsun] of innermost stable circular orbit for BH accretor: */
  /* with mass M and spin aspin */ 
 double Risco,Z1,Z2;
 
 Z1=1.0+pow(1.0-aspin*aspin,0.333333)*(pow(1.0+aspin*aspin,0.333333)+pow(1.0-aspin*aspin,0.333333));
 Z2=sqrt(3.0*aspin*aspin+Z1*Z1);
 Risco=GGG*M*(3.0+Z2-sqrt((3.0-Z1)*(3.0+Z1+2.0*Z2)))/(CCC*CCC);   /* [Rsun] */ 
 
 return Risco;
}


double Eiscof(double R, double M, double aspin) 
{ /* returns specific Keplerian energy [unitless] at radius R from BH accretor: */
  /* with mass M and spin aspin.  (it is to be used for R=Risco) */ 
 double tmp,Eisco;
 
 tmp=(GGG*M)/(CCC*CCC);
 Eisco= (R*R-2.0*tmp*R+aspin*tmp*sqrt(tmp*R))/(R*sqrt(R*R-3.0*tmp*R+2.0*aspin*tmp*sqrt(tmp*R)));
 
 return Eisco;
}


double Mcsnf(double Mcbagb, double Mcco0)
{ /* calculates critical mass of CO core on AGB, for which SN occurs, eq.(75) + correction from: */
  /* hrdiag.f l.359 where they take max of eq.(75) and 1.05*Mcco0 -which is starting CO core mass at BAGB */
  /* global variable Mgrow is here this 1.05 */
  /* if the value of Mcco0 is not relevant (as for Helium Stars) just supply 0.0 */  
 double Mcsn;
 
 Mcsn=max(MCh,0.773*Mcbagb-0.35);     
 Mcsn=max(Mcsn,Mgrow*Mcco0);
   
 return Mcsn;
} 


double Mupf(void)
{ /* calculates initial star critical mass Mup, sec.6 from inversion of eq.(66) */
 double Mcbagb=1.6;

 return pow((pow(Mcbagb,4.0)-bb38)/bb36,1.0/bb37);
}


double Mecf(void)         
{ /* calculates initial star critical mass Mec, sec.6 from inversion of eq.(66) */
 double Mcbagb=2.25;                                                                                 
 
 return pow((pow(Mcbagb,4.0)-bb38)/bb36,1.0/bb37);
}



/* ----------------------------------- INITIAL FUNCTIONS  ---------------------------------------- */

double M_hookf(void)
{ /* calculate M_hook, eq.(1) */
 double dzeta=log10(ZZ/0.02);

 return 1.0185+0.16015*dzeta+0.0892*dzeta*dzeta; 
} 
 

double M_HeFf(void)
{ /* calculate M_HeF, eq.(2) */
 double dzeta=log10(ZZ/0.02);

 return 1.995+0.25*dzeta+0.087*dzeta*dzeta;
}
 

double M_FGBf(void)
{ /* calculate for global variable ZZ M_FGB from what is in Jarrod code: in zcnsts.f */
  /* not from eq.(3) in paper 0, which differs -by small amount- from what is in zcnsts.f */
 
 return 16.5*pow(ZZ,0.06)/(1.0+pow(1.0e-4/ZZ,1.27));
} 


void coef_aa(void)
{ /* sets global coefficients aa using formulae from Appendix */
  /* needs global variable ZZ */  
 double dzeta=log10(ZZ/0.02); 
 double sigma=log10(ZZ);
 double aa11p,aa12p,aa18p,aa19p,aa29p;
 double tmp;

 aa1=1.593890e+3+2.053038e+3*dzeta+1.231226e+3*dzeta*dzeta+2.327785e+2*dzeta*dzeta*dzeta;
 aa2=2.706708e+3+1.483131e+3*dzeta+5.772723e+2*dzeta*dzeta+7.411230e+1*dzeta*dzeta*dzeta;
 aa3=1.466143e+2-1.048442e+2*dzeta-6.795374e+1*dzeta*dzeta-1.391127e+1*dzeta*dzeta*dzeta;
 aa4=4.141960e-2+4.564888e-2*dzeta+2.958542e-2*dzeta*dzeta+5.571483e-3*dzeta*dzeta*dzeta;
 aa5=3.426349e-1;
 
 aa6=1.949814e+1+1.758178*dzeta-6.008212*dzeta*dzeta-4.470533*dzeta*dzeta*dzeta;
 aa7=4.903830;
 aa8=5.212154e-2+3.166411e-2*dzeta-2.750074e-3*dzeta*dzeta-2.271549e-3*dzeta*dzeta*dzeta;
 aa9=1.312179-3.294936e-1*dzeta+9.231860e-2*dzeta*dzeta+2.610989e-2*dzeta*dzeta*dzeta;
 aa10=8.073972e-1;
 
 aa11p=1.031538-2.434480e-1*dzeta+7.732821*dzeta*dzeta+6.460705*dzeta*dzeta*dzeta+1.374484*dzeta*dzeta*dzeta*dzeta;
 aa12p=1.043715-1.577474*dzeta-5.168234*dzeta*dzeta-5.596506*dzeta*dzeta*dzeta-1.299394*dzeta*dzeta*dzeta*dzeta;
 aa13=7.859573e+2-8.542048*dzeta-2.642511e+1*dzeta*dzeta-9.585707*dzeta*dzeta*dzeta;
 aa14=3.858911e+3+2.459681e+3*dzeta-7.630093e+1*dzeta*dzeta-3.486057e+2*dzeta*dzeta*dzeta-4.861703e+1*dzeta*dzeta*dzeta*dzeta;
 aa15=2.888720e+2+2.952979e+2*dzeta+1.850341e+2*dzeta*dzeta+3.797254e+1*dzeta*dzeta*dzeta;
 aa16=7.196580+5.613746e-1*dzeta+3.805871e-1*dzeta*dzeta+8.398728e-2*dzeta*dzeta*dzeta;
 aa11=aa11p*aa14;
 aa12=aa12p*aa14;
 
 aa18p=2.187715e-1-2.154437*dzeta-3.768678*dzeta*dzeta-1.975518*dzeta*dzeta*dzeta-3.021475e-1*dzeta*dzeta*dzeta*dzeta;
 aa19p=1.466440+1.839725*dzeta+6.442199*dzeta*dzeta+4.023635*dzeta*dzeta*dzeta+6.957529e-1*dzeta*dzeta*dzeta*dzeta;
 aa20=2.652091e+1+8.178458e+1*dzeta+1.156058e+2*dzeta*dzeta+7.633811e+1*dzeta*dzeta*dzeta+1.950698e+1*dzeta*dzeta*dzeta*dzeta;
 aa21=1.472103-2.947609*dzeta-3.312828*dzeta*dzeta-9.945065e-1*dzeta*dzeta*dzeta;
 aa22=3.071048-5.679941*dzeta-9.745523*dzeta*dzeta-3.594543*dzeta*dzeta*dzeta;
 aa23=2.617890+1.019135*dzeta-3.292551e-2*dzeta*dzeta-7.445123e-2*dzeta*dzeta*dzeta;
 aa24=1.075567e-2+1.773287e-2*dzeta+9.610479e-3*dzeta*dzeta+1.732469e-3*dzeta*dzeta*dzeta;
 aa25=1.476246+1.899331*dzeta+1.195010*dzeta*dzeta+3.035051e-1*dzeta*dzeta*dzeta;
 aa26=5.502535-6.601663e-2*dzeta+9.968707e-2*dzeta*dzeta+3.599801e-2*dzeta*dzeta*dzeta;
 tmp=max(0.097-0.1072*(sigma+3.0),max(0.097,min(0.1461,0.1461+0.1237*(sigma+2.0))));
 aa17=pow(10.0,tmp);
 aa18=aa18p*aa20;
 aa19=aa19p*aa20;
 
 aa27=9.511033e+1+6.819618e+1*dzeta-1.045625e+1*dzeta*dzeta-1.474939e+1*dzeta*dzeta*dzeta;
 aa28=3.113458e+1+1.012033e+1*dzeta-4.650511*dzeta*dzeta-2.463185*dzeta*dzeta*dzeta;
 aa29p=1.413057+4.578814e-1*dzeta-6.850581e-2*dzeta*dzeta-5.588658e-2*dzeta*dzeta*dzeta;
 aa30=3.910862e+1+5.196646e+1*dzeta+2.264970e+1*dzeta*dzeta+2.873680*dzeta*dzeta*dzeta;
 aa31=4.597479-2.855179e-1*dzeta+2.709724e-1*dzeta*dzeta;
 aa32=6.682518+2.827718e-1*dzeta-7.294429e-2*dzeta*dzeta;
 aa29=pow(aa29p,aa32); 
 
 aa34=1.910302e-1+1.158624e-1*dzeta+3.348990e-2*dzeta*dzeta+2.599706e-3*dzeta*dzeta*dzeta;
 aa35=3.931056e-1+7.277637e-2*dzeta-1.366593e-1*dzeta*dzeta-4.508946e-2*dzeta*dzeta*dzeta;
 aa36=3.267776e-1+1.204424e-1*dzeta+9.988332e-2*dzeta*dzeta+2.455361e-2*dzeta*dzeta*dzeta;
 aa37=5.990212e-1+5.570264e-2*dzeta+6.207626e-2*dzeta*dzeta+1.777283e-2*dzeta*dzeta*dzeta;
 aa33=min(1.4,1.5135+0.3769*dzeta);
 aa33=max(0.6355-0.4192*dzeta,max(1.25,aa33));
 
 aa38=7.330122e-1+5.192827e-1*dzeta+2.316416e-1*dzeta*dzeta+8.346941e-3*dzeta*dzeta*dzeta;
 aa39=1.172768-1.209262e-1*dzeta-1.193023e-1*dzeta*dzeta-2.859837e-2*dzeta*dzeta*dzeta;
 aa40=3.982622e-1-2.296279e-1*dzeta-2.262539e-1*dzeta*dzeta-5.219837e-2*dzeta*dzeta*dzeta;
 aa41=3.571038-2.223625e-2*dzeta-2.611794e-2*dzeta*dzeta-6.359648e-3*dzeta*dzeta*dzeta;
 aa42=1.9848+1.1386*dzeta+3.5640e-1*dzeta*dzeta;
 aa43=6.300e-2+4.810e-2*dzeta+9.840e-3*dzeta*dzeta;
 aa44=1.200+2.450*dzeta;
 aa42=min(1.25,max(1.1,aa42));
 aa44=min(1.3,max(0.45,aa44));
 
 aa45=2.321400e-1+1.828075e-3*dzeta-2.232007e-2*dzeta*dzeta-3.378734e-3*dzeta*dzeta*dzeta;
 aa46=1.163659e-2+3.427682e-3*dzeta+1.421393e-3*dzeta*dzeta-3.710666e-3*dzeta*dzeta*dzeta;
 aa47=1.048020e-2-1.231921e-2*dzeta-1.686860e-2*dzeta*dzeta-4.234354e-3*dzeta*dzeta*dzeta;
 aa48=1.555590-3.223927e-1*dzeta-5.197429e-1*dzeta*dzeta-1.066441e-1*dzeta*dzeta*dzeta;
 aa49=9.7700e-2-2.3100e-1*dzeta-7.5300e-2*dzeta*dzeta;
 aa50=2.4000e-1+1.8000e-1*dzeta+5.9500e-1*dzeta*dzeta;
 aa51=3.3000e-1+1.3200e-1*dzeta+2.1800e-1*dzeta*dzeta;
 aa52=1.1064+4.1500e-1*dzeta+1.8000e-1*dzeta*dzeta;
 aa53=1.1900+3.7700e-1*dzeta+1.7600e-1*dzeta*dzeta;
 aa49=max(aa49,0.145);
 aa50=min(aa50,0.306+0.053*dzeta);
 aa51=min(aa51,0.3625+0.062*dzeta);
 aa52=max(aa52,0.9);
 if(ZZ>0.01) aa52=min(aa52,1.0);
 aa53=max(aa53,1.0);
 if(ZZ>0.01) aa53=min(aa53,1.1);
         
 aa54=3.855707e-1-6.104166e-1*dzeta+5.676742*dzeta*dzeta+1.060894e+1*dzeta*dzeta*dzeta+5.284014*dzeta*dzeta*dzeta*dzeta;
 aa55=3.579064e-1-6.442936e-1*dzeta+5.494644*dzeta*dzeta+1.054952e+1*dzeta*dzeta*dzeta+5.280991*dzeta*dzeta*dzeta*dzeta;
 aa56=9.587587e-1+8.777464e-1*dzeta+2.017321e-1*dzeta*dzeta;
 aa57=1.5135+3.7690e-1*dzeta;
 aa57=min(1.4,aa57);
 aa57=max(0.6355-0.4192*dzeta,max(1.25,aa57));         
         
 aa58=4.907546e-1-1.683928e-1*dzeta-3.108742e-1*dzeta*dzeta-7.202918e-2*dzeta*dzeta*dzeta;
 aa59=4.537070-4.465455*dzeta-1.612690*dzeta*dzeta-1.623246*dzeta*dzeta*dzeta;
 aa60=1.796220+2.814020e-1*dzeta+1.423325*dzeta*dzeta+3.421036e-1*dzeta*dzeta*dzeta;
 aa61=2.256216+3.773400e-1*dzeta+1.537867*dzeta*dzeta+4.396373e-1*dzeta*dzeta*dzeta;
 aa62=8.4300e-2-4.7500e-2*dzeta-3.5200e-2*dzeta*dzeta;
 aa63=7.3600e-2+7.4900e-2*dzeta+4.4260e-2*dzeta*dzeta;
 aa64=1.3600e-1+3.5200e-2*dzeta;
 aa65=1.564231e-3+1.653042e-3*dzeta-4.439786e-3*dzeta*dzeta-4.951011e-3*dzeta*dzeta*dzeta-1.216530e-3*dzeta*dzeta*dzeta*dzeta;
 aa66=1.4770+2.9600e-1*dzeta;
 aa67=5.210157-4.143695*dzeta-2.120870*dzeta*dzeta; 
 aa68=1.1160+1.6600e-1*dzeta;
 aa62=max(0.065,aa62);
 if(ZZ<0.004) aa63=min(0.055,aa63);
 aa64=max(0.091,min(0.121,aa64));
 aa66=max(aa66,min(1.6,-0.308-1.046*dzeta));
 aa66=max(0.8,min(0.8-2.0*dzeta,aa66));
 aa68=max(0.9,min(aa68,1.0));
 if(aa68>aa66) aa64=aa58*pow(aa66,aa60)/(aa59*pow(aa66,aa61));
 aa68=min(aa68,aa66);    

 aa69=1.071489-1.164852e-1*dzeta-8.623831e-2*dzeta*dzeta-1.582349e-2*dzeta*dzeta*dzeta;                   
 aa70=7.108492e-1+7.935927e-1*dzeta+3.926983e-1*dzeta*dzeta+3.622146e-2*dzeta*dzeta*dzeta;
 aa71=3.478514-2.585474e-2*dzeta-1.512955e-2*dzeta*dzeta-2.833691e-3*dzeta*dzeta*dzeta;
 aa72=9.132108e-1-1.653695e-1*dzeta+3.636784e-2*dzeta*dzeta*dzeta;
 aa73=3.969331e-3+4.539076e-3*dzeta+1.720906e-3*dzeta*dzeta+1.897857e-4*dzeta*dzeta*dzeta;
 aa74=1.600+7.640e-1*dzeta+3.322e-1*dzeta*dzeta;
 if(ZZ>0.01) aa72=max(aa72,0.95);       
 aa74=max(1.4,min(aa74,1.6));        
 
 aa75=8.109e-1-6.282e-1*dzeta;
 aa76=1.192334e-2+1.083057e-2*dzeta+1.230969*dzeta*dzeta+1.551656*dzeta*dzeta*dzeta;
 aa77=-1.668868e-1+5.818123e-1*dzeta-1.105027e+1*dzeta*dzeta-1.668070e+1*dzeta*dzeta*dzeta;
 aa78=7.615495e-1+1.068243e-1*dzeta-2.011333e-1*dzeta*dzeta-9.371415e-2*dzeta*dzeta*dzeta;
 aa79=9.409838+1.522928*dzeta;
 aa80=-2.7110e-1-5.7560e-1*dzeta-8.3800e-2*dzeta*dzeta;
 aa81=2.4930+1.1475*dzeta;
 aa75=max(1.0,min(aa75,1.27));
 aa75=max(aa75,0.6355-0.4192*dzeta);
 aa76=max(aa76,-0.1015564-0.2161264*dzeta-0.05182516*dzeta*dzeta);
 aa77=max(-0.3868776-0.5457078*dzeta-0.1463472*dzeta*dzeta,min(0.0,aa77));
 aa78=max(0.0,min(aa78,7.454+9.046*dzeta));
 aa79=min(aa79,max(2.0,-13.3-18.6*dzeta));
 aa80=max(0.0585542,aa80);
 aa81=min(1.5,max(0.4,aa81));
} 


void coef_bb(void)
{ /* sets global coefficients bb using formulae from Appendix */
  /* needs global variables ZZ,M_HeF,M_FGB */
 double dzeta=log10(ZZ/0.02); 
 double sigma=log10(ZZ);
 double ro;
 double bb3p,bb11p,bb13p,bb14p,bb16p,bb24p,bb27p,bb31p,bb34p,bb36p,bb37p,bb38p,bb41p,bb44p,bb51p,bb53p;
 double bb56p,bb57p;  
 
 ro=dzeta+1.0;
   
 bb1=3.9700e-1+2.8826e-1*dzeta+5.2930e-1*dzeta*dzeta;
 bb4=9.960283e-1+8.164393e-1*dzeta+2.383830*dzeta*dzeta+2.223436*dzeta*dzeta*dzeta+8.638115e-1*dzeta*dzeta*dzeta*dzeta;
 bb5=2.561062e-1+7.072646e-2*dzeta-5.444596e-2*dzeta*dzeta-5.798167e-2*dzeta*dzeta*dzeta-1.349129e-2*dzeta*dzeta*dzeta*dzeta;
 bb6=1.157338+1.467883*dzeta+4.299661*dzeta*dzeta+3.130500*dzeta*dzeta*dzeta+6.992080e-1*dzeta*dzeta*dzeta*dzeta;
 bb7=4.022765e-1+3.050010e-1*dzeta+9.962137e-1*dzeta*dzeta+7.914079e-1*dzeta*dzeta*dzeta+1.728098e-1*dzeta*dzeta*dzeta*dzeta;
 bb1=min(0.54,bb1);
 bb2=pow(10.0,-4.6739-0.9394*sigma);
 bb2=min(max(bb2,-0.04167+55.67*ZZ),0.4771-9329.21*pow(ZZ,2.94));
 bb3p=max(-0.1451,-2.2794-1.5175*sigma-0.254*sigma*sigma);
 bb3=pow(10.0,bb3p);
 if(ZZ>0.004) bb3=max(bb3,0.7307+14265.1*pow(ZZ,3.395));
 bb4=bb4+0.1231572*pow(dzeta,5.0);
 bb6=bb6+0.01640687*pow(dzeta,5.0); 

 bb9=2.751631e+3+3.557098e+2*dzeta;
 bb10=-3.820831e-2+5.872664e-2*dzeta;
 bb11p=1.071738e+2-8.970339e+1*dzeta-3.949739e+1*dzeta*dzeta;
 bb12=7.348793e+2-1.531020e+2*dzeta-3.793700e+1*dzeta*dzeta;
 bb13p=9.219293-2.005865*dzeta-5.561309e-1*dzeta*dzeta;
 bb11=bb11p*bb11p;
 bb13=bb13p*bb13p;

 bb14p=2.917412+1.575290*dzeta+5.751814e-1*dzeta*dzeta;
 bb15= 3.629118-9.112722e-1*dzeta+1.042291*dzeta*dzeta;
 bb16p=4.916389+2.862149*dzeta+7.844850e-1*dzeta*dzeta;
 bb14=pow(bb14p,bb15);
 bb16=pow(bb16p,bb15);
 bb17=1.0;
 if(dzeta>-1.0) bb17=1.0-0.3880523*pow(dzeta+1.0,2.862149);

 bb18=5.496045e+1-1.289968e+1*dzeta+6.385758*dzeta*dzeta;
 bb19=1.832694-5.766608e-2*dzeta+5.696128e-2*dzeta*dzeta;
 bb20=1.211104e+2;
 bb21=2.214088e+2+2.187113e+2*dzeta+1.170177e+1*dzeta*dzeta-2.635340e+1*dzeta*dzeta*dzeta;
 bb22=2.063983+7.363827e-1*dzeta+2.654323e-1*dzeta*dzeta-6.140719e-2*dzeta*dzeta*dzeta;
 bb23=2.003160+9.388871e-1*dzeta+9.656450e-1*dzeta*dzeta+2.362266e-1*dzeta*dzeta*dzeta;
 bb24p=1.609901e+1+7.391573*dzeta+2.277010e+1*dzeta*dzeta+8.334227*dzeta*dzeta*dzeta;
 bb25=1.747500e-1+6.271202e-2*dzeta-2.324229e-2*dzeta*dzeta-1.844559e-2*dzeta*dzeta*dzeta;
 bb27p=2.752869+2.729201e-2*dzeta+4.996927e-1*dzeta*dzeta+2.496551e-1*dzeta*dzeta*dzeta;
 bb28=3.518506+1.112440*dzeta-4.556216e-1*dzeta*dzeta-2.179426e-1*dzeta*dzeta*dzeta;
 bb24=pow(bb24p,bb28);
 bb26=5.0-0.09138012*pow(ZZ,-0.3671407);
 bb27=pow(bb27p,2.0*bb28);

 bb29=1.626062e+2-1.168838e+1*dzeta-5.498343*dzeta*dzeta;
 bb30=3.336833e-1-1.458043e-1*dzeta-2.011751e-2*dzeta*dzeta;
 bb31p=7.425137e+1+1.790236e+1*dzeta+3.033910e+1*dzeta*dzeta+1.018259e+1*dzeta*dzeta*dzeta;
 bb32=9.268325e+2-9.739859e+1*dzeta-7.702152e+1*dzeta*dzeta-3.158268e+1*dzeta*dzeta*dzeta;
 bb33=2.474401+3.892972e-1*dzeta;
 bb34p=1.127018e+1+1.622158*dzeta-1.443664*dzeta*dzeta-9.474699e-1*dzeta*dzeta*dzeta;
 bb31=pow(bb31p,bb33);
 bb34=pow(bb34p,bb33);

 bb36p=1.445216e-1-6.180219e-2*dzeta+3.093878e-2*dzeta*dzeta+1.567090e-2*dzeta*dzeta*dzeta;
 bb37p=1.304129+1.395919e-1*dzeta+4.142455e-3*dzeta*dzeta-9.732503e-3*dzeta*dzeta*dzeta;
 bb38p=5.114149e-1-1.160850e-2*dzeta;
 bb36=pow(bb36p,4.0);
 bb37=4.0*bb37p;
 bb38=pow(bb38p,4.0);

 bb39=1.314955e+2+2.009258e+1*dzeta-5.143082e-1*dzeta*dzeta-1.379140*dzeta*dzeta*dzeta;
 bb40=1.823973e+1-3.074559*dzeta-4.307878*dzeta*dzeta;
 bb41p=2.327037+2.403445*dzeta+1.208407*dzeta*dzeta+2.087263e-1*dzeta*dzeta*dzeta;
 bb42=1.997378-8.126205e-1*dzeta;
 bb43=1.079113e-1+1.762409e-2*dzeta+1.096601e-2*dzeta*dzeta+3.058818e-3*dzeta*dzeta*dzeta;
 bb44p=2.327409+6.901582e-1*dzeta-2.158431e-1*dzeta*dzeta-1.084117e-1*dzeta*dzeta*dzeta;
 bb40=max(bb40,1.0);
 bb41=pow(bb41p,bb42);
 bb44=pow(bb44p,5.0);

 bb46=2.214315-1.975747*dzeta;
 bb48=5.072525+1.146189e+1*dzeta+6.961724*dzeta*dzeta+1.316965*dzeta*dzeta*dzeta;
 bb49=5.139740;
 bb45=1.0-(2.47162*ro-5.401682*ro*ro+3.247361*ro*ro*ro);
 if(ro<=0.0) bb45=1.0;
 bb46=-1.0*bb46*log10(M_HeF/M_FGB);
 bb47=1.127733*ro+0.2344416*ro*ro-0.3793726*ro*ro*ro;

 bb51p=1.125124+1.306486*dzeta+3.622359*dzeta*dzeta+2.601976*dzeta*dzeta*dzeta+3.031270e-1*dzeta*dzeta*dzeta*dzeta;
 bb52=3.349489e-1+4.531269e-3*dzeta+1.131793e-1*dzeta*dzeta+2.300156e-1*dzeta*dzeta*dzeta+7.632745e-2*dzeta*dzeta*dzeta*dzeta;
 bb53p=1.467794+2.798142*dzeta+9.455580*dzeta*dzeta+8.963904*dzeta*dzeta*dzeta+3.339719*dzeta*dzeta*dzeta*dzeta;
 bb54=4.658512e-1+2.597451e-1*dzeta+9.048179e-1*dzeta*dzeta+7.394505e-1*dzeta*dzeta*dzeta+1.607092e-1*dzeta*dzeta*dzeta*dzeta;
 bb55=1.0422+1.3156e-1*dzeta+4.5000e-2*dzeta*dzeta;
 bb56p=1.110866+9.623856e-1*dzeta+2.735487*dzeta*dzeta+2.445602*dzeta*dzeta*dzeta+8.826352e-1*dzeta*dzeta*dzeta*dzeta;
 bb57p=-1.584333e-1-1.728865e-1*dzeta-4.461431e-1*dzeta*dzeta-3.925259e-1*dzeta*dzeta*dzeta-1.276203e-1*dzeta*dzeta*dzeta*dzeta;
 bb51=bb51p-0.1343798*pow(dzeta,5.0);
 bb53=bb53p+0.4426929*pow(dzeta,5.0);
 bb55=min(0.99164-743.123*pow(ZZ,2.83),bb55);
 bb56=bb56p+0.1140142*pow(dzeta,5.0);
 bb57=bb57p-0.01308728*pow(dzeta,5.0);
}


void init_SN(void)
{ /* initialize SN counting tables: for rate callibration */
  /* sII[2] stores number of SNII happand at time t=[2Gyr:3Gyr] */
  /* sII[9] stores number of SNIbc happand at time t=[9Gyr:10Gyr] */ 
 int i;
 
 for(i=0;i<20;i++) 
   sII[i]=sIbc[i]=sFa[i]=0;
}


/* --------------------------- OTHER FUNCTIONS -------------------------------------------- */

double get_T(double L, double R)
{ /* for star of L[L_sun], R[R_sun] returns effective temperature in [K] */
 double sigm=5.6724e-05*Rsun*Rsun/Lsun;     /* Stefan constant [L_sun R_sun^{-2} K^{-4}] */
  
 return pow(L/(4.0*Pi*R*R*sigm),0.25);
}

   
double inter_line(double x1, double y1, double x2, double y2, double x)
{ /* interpolates straight line between two points: (x1,y1) i (x2,y2) */
  /* returns value y of this line at point x */ 
 double a,b;
 
 a=(y1-y2)/(x1-x2); 
 b=y1-a*x1;
  
 return a*x+b; 
}


/*--------------------------- funtions for binary.c -------------------------------------------- */

double tnfI(double M, double m0, int K)
{ /* Feb 11, 2002: calculates the nuclear timescale, from Jarod code: star.f l.220-284 */
  /* uses tmax -estimate of time when Mcco reaches Mcmax on AGB: star.f l.193-218 */
  /* m0 should be filled with M0 in binary.c */ 
 double tn,tmax,tau,Mt,Mx,Mcbagb,mc1,mc2,mcmax;
 double lambda,D,B,p,q,Ahp,Ahe,Ahhe;

 Mt=M;
 Mcbagb=Mcbagbf(m0);
 Mx=Mxf(m0);
 D=Df(m0);
 B=Bf(m0);
 p=pf(m0);
 q=qf(m0);
 Ahp=Ahpf(m0);
 Ahe=Ahef(); 
 Ahhe=Ahhef();
 lambda=min(0.9,0.3+0.001*pow(m0,5.0));      


 if(K<5)                              /* set mc1 as in star.f */
   mc1=Mcgbf2(m0,Lhe1f(m0));
 else {
   mc1=Mcbagb;                        /* mc1 -He core mas on EAGB and at start of TPAGB, star.f l.163 */
   if(mc1>=0.8 && mc1<=2.25)
     mc1=0.44*mc1+0.448;
 }


 tau=the1f(m0)+thef(m0);                           /* this is time not a timescale!, name from star.f */
 mc2=Mceagbf1(m0,tau);              
 mcmax=max(max(MCh,0.773*Mcbagb-0.35),1.05*mc2);   /* similar to my Mcsn */ 
 if(mcmax<=mc1) { 
   if(mcmax<=Mx) 
     tmax=tinf1ef(m0)-pow(mcmax,1.0-p)/((p-1.0)*Ahe*D);
   else
     tmax=tinf2ef(m0)-pow(mcmax,1.0-q)/((q-1.0)*Ahe*B);  
 }
 else {
   mcmax=(mcmax-lambda*mc1)/(1.0-lambda);
   if(mcmax<=Mx)
     tmax=tinf1tf(m0)-pow(mcmax,1.0-p)/((p-1.0)*Ahhe*D);
   else 
     tmax=tinf2tf(m0)-pow(mcmax,1.0-q)/((q-1.0)*Ahhe*B);
 }      
 tmax=max(tbagbf(m0),tmax);


 if(K<5 && fabs(Mt-Mcbagb)<acc)      
   tn=tbagbf(m0);
 else {
   if(Mt>Mcbagb || (Mt>=mc1 && K>4)) {
     if(K==6) 
       mc1=(Mt-lambda*mc1)/(1.0-lambda); 
     else 
       mc1=Mt;
     if(mc1<=Mx)
       tn=tinf1tf(m0)-pow(mc1,1.0-p)/((p-1.0)*Ahhe*D); 
     else
       tn=tinf2tf(m0)-pow(mc1,1.0-q)/((q-1.0)*Ahhe*B);  
   }   
   else {
     if(m0>M_FGB) {
       mc1=Mche1f(m0);
       if(Mt<=mc1)
         tn=the1f(m0);
       else
         tn=the1f(m0)+thef(m0)*((Mt-mc1)/(Mcbagb-mc1)); 
     } 
     else if(m0<=M_HeF) {
       if(K<8) {                                 /* for EAGB and TPAGB the same functions as for RG apply */
         mc1=Mcgbf2(m0,Lbgbf(m0));
         mc2=Mcgbf2(m0,Lhe1f(m0)); 
       }
       else {
         mc1=Mchsgbf2(m0,Lbgbf(m0));
         mc2=Mchsgbf2(m0,Lhe1f(m0));
       }
       if(Mt<=mc1)
         tn=tbgbf(m0);
       else if(Mt<=mc2) {
         if(Mt<=Mx)
           tn=tinf1f(m0)-pow(Mt,1.0-p)/((p-1.0)*Ahp*D);
         else
           tn=tinf2f(m0)-pow(Mt,1.0-q)/((q-1.0)*Ahp*B);                       
       }
       else
         tn=the1f(m0)+thef(m0)*((Mt-mc2)/(Mcbagb-mc2));
     }
     else {
       mc1=Mcbgbf1a(m0);
       mc2=Mche1f(m0);  
       if(Mt<=mc1)
         tn=tbgbf(m0);
       else if(Mt<=mc2)  
         tn=tbgbf(m0)+(the1f(m0)-tbgbf(m0))*((Mt-mc1)/(mc2-mc1));
       else
         tn=the1f(m0)+thef(m0)*((Mt-mc2)/(Mcbagb-mc2));  
     }   
   }    
 }
 tn=min(tn,tmax);                        

 return tn;
} 

/*---------------------------------------- THE END -------------------------------------------- */
