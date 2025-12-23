/* software by chris belczynski: last modified on Nov21, 2022 */
/* header for singl.c and binary.c */
/* parameters set as for model M30 with delayed engine and qmin and correct dJf() */

#ifdef External
#define St extern
#else
#define St
#endif



#define num_tested 50                                               /* number of tested systems */
#define NOUT 1000                          /* frequency of output into info.dat (every NOUT systems) */


#define BINARY 1                          /* 1 -- runs binaries only: num_tested binaries (STANDARD) */
                                         /* 0 -- runs single stars only: 2 * num_tested single stars */
                                         /* use IMF=1, Mmina,Mmaxa and Sal sets IMF for single stars */

#define BINOUT 0                          /* 0/1 -- no/output (general) for binary stars: bstars.dat */
                              /* general output for single stars (BINARY=0) is always on: sstars.dat */

#define COMPACT 0                /* 0/1 -- no/output for compact object binaries active: compact.dat */

#define WDWD 0                              /* 0/1 -- no/output for WD-WD at the formation: wdwd.dat */
                      /* no output for binaries with NS, BH or SN Ia is valid once this option is ON */

#define Xout1 0                     /* 0/1 -- no/output for RLOF X-ray and NS/BH binaries: DetMT.dat */
#define Lxcrit1 (1.0e+38)          /* output only if source beamed Lx is larger than Lxcrit1 [erg/s] */

#define Xout2 0                   /* 0/1 -- no/output for WIND X-ray and NS/BH binaries: DetWIND.dat */
#define Lxcrit2 (1.0e+38)          /* output only if source beamed Lx is larger than Lxcrit2 [erg/s] */



#define SINOUT 0                /* 0/1 -- no/output for single stars from disrupted binaries: sb.dat */


#define BHCAT1 0                                         /* 1-BH catalog generated, 0- not generated */
#define TIMEOUT1 (8.673997)                       /* corresponds to death of 25 Msun star at Z=0.001 */
#define TIMEOUT2 (11.017251)                      /* corresponds to death of 20 Msun star at Z=0.001 */
#define TIMEOUT3 (15.779993)                      /* corresponds to death of 15 Msun star at Z=0.001 */
#define TIMEOUT4 (41.7)                            /* corresponds to death of 8 Msun star at Z=0.001 */
#define TIMEOUT5 (103.8)                           /* corresponds to death of 5 Msun star at Z=0.001 */
#define TIMEOUT 10000.0                        /* [Myrs] center of output time for Symbiotic systsms */
#define DTOUT 1000.0                                                 /* [Myrs] time range for output */

#define LONGGRB 0                           /* 0/1 -- no/output for potential long GRBs: longgrb.dat */


#define SN1A 0                                      /* 1- output for SN 1a active, 0 -- disactivated */

#define AMCVN 0                                    /* 1- output for AM CVn active, 0 -- disactivated */

#define WDWDout 0                                   /* 1 -- output for LISA binaries, 0 -- no output */
#define Twdwd 10000.0                                        /* output time for LISA binaries [Myrs] */

#define CVout 0                                       /* 1- output for CVs active, 0 -- disactivated */
#define Tcv 10000.0                                                        /* output time for cv.dat */


#define OBCAT 0                                  /* 1- output for OB stars active, 0 -- disactivated */
#define TIMEOB 10000.0                                                  /* output time for OBbin.dat */

#define BHCAT2 0                            /* 1- output to bhsys.dat file activated, 0-disactivated */

#define CLUSTER 0                 /* 1-active, 0-disactivatied: for cluster simulations: cluster.dat */

#define SSout 0                                    /* 0/1 -- no/output for symbiotic stars: symb.dat */

#define RRwd 0                                                /* 0/1 -- no/output for RR Lyrae stars */

#define Merger 0                         /* 0/1 -- no/output for mergers of NS/BH with regular stars */

#define TranCE 0                      /* 0/1 -- no/output for CE events: they can produce transients */
                                                        /* WARNING: MAKE IT 0 FOR ALL OTHER PROJECTS */
                                                 /* 1 TERMINATES EVOLUTION OF ANY BINARY AT FIRST CE */

#define BHNSWD 0                      /* 1 -- output for BH-WD and NS-WD binaries, 0 -- disactivated */




#define SS 1                                 /* 0 -- old (standard) initial distributions: upto 2014 */
                                        /* 1 -- new Sana et al. 2012 initial parameter distributions */
                                               /* 2 -- POP III initial distributions: Taeho Ryu 2016 */

#define ZZ 1.03e-04                    /* metallicity: allowed range: 0.0001--0.03 (solar: ZZsun=0.014) */


#define REMNANT 12                      /* 0: StarTrack remnant mass calculation: Timmes et al. 1996 */
                                              /* 1: 2012 Chris F. remmant mass claculation: RAPID SN */
                                            /* 2: 2012 Chris F. remmant mass claculation: DELAYED SN */
                       /* 3:   RAPID SN + strong PPSN/PSN   + adjustable neutrino emis./CO core size */
                       /* 4: DELAYED SN + strong PPSN/PSN   + adjustable neutrino emis./CO core size */
                       /* 5:   RAPID SN + moderate PPSN/PSN + adjustable neutrino emis./CO core size */
                       /* 6: DELAYED SN + moderate PPSN/PSN + adjustable neutrino emis./CO core size */
                       /* 7:   RAPID SN + weak PPSN/PSN     + adjustable neutrino emis./CO core size */
                       /* 8: DELAYED SN + weak PPSN/PSN     + adjustable neutrino emis./CO core size */
                       /* 9:  RAPID SN + no PPSN + high PSN + adjustable neutrino emis./CO core size */
                                         /* 10: DELAYED SN + no PPSN + high PSN  + extended fallback */
                                           /* 12: DELAYED SN + weak PPSN/PSN     + extended fallback */
         /* 13: SMOOTH SN engine with FLEXIBLE PPSN/PSN: set by fmix, XXPSN, PPSN, PSN: Olka version */
           /* 14: SMOOTH SN engine with FLEXIBLE PPSN/PSN: set by fmix, XXPSN, PPSN, PSN: KB version */



                                                      /* fmix and McritNSBH work only for REMNANT=13 */
#define fmix 4.0                            /* mixing parameter in range 0.5-4.0 (Fryer et al. 2021) */
                                          /* fmix=0.5: quasi delayed engine (shallow lower mass gap) */
                                               /* fmix=4.0: quasi rapid engine (deep lower mass gap) */
                                            /* fmix=2.0: STANDARD -- partially filled lower mass gap */
#define McritNSBH 4.75               /* [Msun] mass of CO core above which BH forms (below NS forms) */
                                                          /* McritNSBH=5.75: gives mass gap ~2-7Msun */
                                                          /* McritNSBH=4.75: gives mass gap ~2-5Msun */

#define XXPSN 3                      /* controls workings of PPSN and PSN: works only for REMNANT=13 */
                                                                      /* 1: STRONG PPSN + NORMAL PSN */
                                                                    /* 2: MODERATE PPSN + NORMAL PSN */
                                                                        /* 3: WEAK PPSN + NORMAL PSN */
                                                                      /* 4: NO PPSN + ULTRA HIGH PSN */
                         /* 5: adjustable PPSN + adjustable PSN; set by Mppsn1, Mppsn2, Mpsn1, Mpsn2 */

#define PPSN 1                        /* Pair-instability pulsation supernova mass loss: 1-ON, 0-OFF */
#define PSN 1                             /* Pair-instability supernova star disruption: 1-ON, 0-OFF */

                                                         /* these parameters work only with XXPSN=5: */
#define Mppsn1 40.0                    /* lower (He core) mass bound on PPSN [Msun], 40Msun standard */
#define Mppsn2 65.0                    /* upper (He core) mass bound on PPSN [Msun], 65Msun standard */
#define Mpsn1 65.0                      /* lower (He core) mass bound on PSN [Msun], 65Msun standard */
#define Mpsn2 135.0                    /* upper (He core) mass bound on PSN [Msun], 135Msun standard */
                                                            /* for standard physics set Mppsn2=Mpsn1 */


#define McoExt 11.0                    /* 11.0 Msun: standard for rapid and delayed SN 2012 SN model */
                                   /* CO core mass above which there is direct BH formation (fb=1.0) */


#define RMAX 1             /* 0: standard radius evolution from Hurley et al. 2000 for PopI/II stars */
                                                           /* radii upto 6000 Rsun for massive stars */
                                        /* 1: maximum radius imposed based on Mzams from Paper 2023 */
                                  /* radii upto 300-3000Rsun depending on Mzams and independent of Z */



#define GRB 0           /* 0 -- no extra mass loss during direct BH formation: compactfN() functions */
                  /* 1 -- extra mass loss (from GRB/SNIc for direct BHs that did not go through PPSN */
                          /* mass loss is taken randomly (uniform distr.) in range: mgrbmin--mgrbmax */
#define mgrbmin 0.0                      /*  range adopted from Ashall et al. 2019: arXiv:1702.04339 */
#define mgrbmax 10.0
#define Mbhmin 30.0         /* above only for very massive BHs (SNIc/GRB objects): Mbhmin<Mbh<Mbhmax */
#define Mbhmax 50.0



                             /* increase of CO/BH mass due to RLOF accretion spinup of MS progenitor */
#define Spin1 0                  /* 0 -- off (standard?), 1 -- on (tests?): works only for REMNANT=3 */
#define Zcrit1 0.002        /* only stars with low Z<=Zcrit1=0.002 (10%) undergo CO/BH mass increase */
#define Fincr1 1.2           /* Mco increases by Fincr1: 1.2 (20%): Geneva models with 40% crit rot. */
#define Fcrit1 0.1         /* fraction of progenitor mass accreted during MS needs to be higher than */
                            /* few percent (we assume Fcrit1=0.1: 10%) to spin-up the accreting star */
                  /* all this works for core-collapse NS/BH, it does not work for ECS NS nor for WDs */

#define Mprotons 0.1                /* 0.1 Msun -- standard: increase of initial protoNS baryon mass */
                    /* correction to formulae presented in Fryer, Belcznski et al. 2012, ApJ 749, 91 */
                              /* works with REMNANT=1 and 2, increases/decreases protoNS baryon mass */
                                 /* Chris F. finds protoNS mass Mex=1.0 Msun, but then our NS masses */
                                         /* are about 0.1 Msun too low as compared with obsrevations */
                      /* Chris F error on Mex is at least +/- 0.1 Msun, so I increase it by 0.1 Msun */

#define Fneut 0.1                       /* old neutrino mass loss for NS and BH formation, 0.1 (10%) */
                                                                           /* used for REMNANT=0,1,2 */

#define Fneut1 0.1              /* new neutrino mass loss for light (<3Msun) BH formation, 0.1 (10%) */
#define Fneut2 0.01       /* new neutrino mass loss for normal mass (>3Msun) BH formation, 0.01 (1%) */
                                                                /* both used for REMNANT=3,4,5,6,7,8 */


#define IMF 3                 /* 1-Kroupa 1993 broken power law: -1.3/-2.2/Sal; Mzamsa=[Mmina-Mmaxa] */
                                     /* Mzamsb: [Mminb-Mmaxa] from q distr. q=[0-1]: Mzamsb=q*Mzamsa */
                /* in needs to be: Mmaxa,Mmaxb>1.0 */
                /* 2-BOTH component masses drawn independently from the Kroupa: -1.3/-2.2/Sal distr. */
                               /* Mzamsa=[Mmina,Mmaxa], Mzamsb=[Mminb,Mmaxb] and no q-distr. used!!! */
                     /* 3-BOTH Mzamsa from linear distribution and Mzamsb from linear q-distribution */
#define Sal (-2.3)        /* IMF exponent for massive stars M>1 Msun: -2.35 (clusters), -2.7 (field) */
#define Mmina 0.1              /* Z=0.02: Mmina=6, Mminb=4 -- DCO, 0.08 -- entire stellar population */
#define Mmaxa 150.0                                                               /* 150 - std value */
#define Mminb 0.1                                              /* Z=0.0001: Mmina=5, Mminb=3 -- DCO */
#define Mmaxb 150.0                                                         /* works only with IMF=2 */

#define SanaOriginal 0       /* sets range of distribution for intial mass ratio: q0 (IMF=1 && SS=1) */
                            /* 0 -- standard: q0 range: qmin-1.0; qmin allows Mb as low as 0.08 Msun */
                            /* 1 -- original Sana: q0 range 0.1--1.0 (no extreme mass ratio binaries */



#define MarkRat 1                                                                        /* MarkRat: */
#define Rat (0.0)                                 /* 1 with Rat=0.0 -- flat q distribution: STANDARD */
                                                          /* 2 -- q^(Rat) distribution, Rat=positive */
                               /* 3 -- q=const + q^(Rat) distribution (limit at q=0.2), Rat=negative */

                                                         /* 4 -- twin binaries (Rat=0), qdivide~0.95 */
#define qdivide 0.95         /* qratio (fraction): of TWINS with flat q distr in range (qdivide-1.0) */
#define qratio 0.5            /* 1-qratio: REGULAR systems with flat q distr in range (qmin-qdivide) */
                                                          /* qratio: Stanek: max 0.5, but rather 0.2 */
                                                /* qratio:  Fryer: max 0.25, but rather like 0.2-0.1 */

#define famin 2.0              /* radius multiplication factor in setting minimal initial orbit size */
                                     /* amin=(famin*Ra0+famin*Rb0)/(1.0-e0); a0=get_a(amin,1.0e+05); */
                        /* 2.0 -- standard, used exclusively before 2014; 1.0 means contact possible */

#define amax (1.0e+05)                        /* 1.0e+05 standard: maximum orbital separation [Rsun] */



double St Tst;                                                 /* formation time of the system [Myr] */
double St t_hubble;                        /* time from formation of system to allowed end (hub_val) */
#define hub_val 1000.0                             /* time marking desired stop of calculation [Myr] */
#define dtSFR  0.0            /* time marking end of Star Formation [Myr] which starts at time = 0.0 */
                            /* set dtSFR=0 -for outburst SF, and dtSFR=[0:hub_val] for continuous SF */
                                     /* hub_val=200.0 for Z=0.02 DCO; hub_val=300.0 for Z=0.0001 DCO */


#define WIND 2                               /* 1 -- (OLD) Hurley et al. (2000) wind mass loss rates */
                     /* 2 -- (STANDARD) new rates: Vink for massive stars; Hurley for low mass stars */
                    /* Hamann&Koesterke for helium stars with Z-dependence from Vink & de Koter 2005 */
                           /* 3 -- new winds (for massive H and He-rich stars from Alex Gormaz: 2022 */


#define Z11 (0.86)       /* wind-metallicity dependence for Helium stars: 0.86 -- Vink&deKotter 2005 */
                            /* set Z11=0.0 for no metallicity dependence. works for both WIND=1 or 2 */

#define Z22 (0.85)           /* wind-metallicity dependence for H-rich stars with Teff: 12500-25000K */
                             /* 0.85 -- Vink et al. 2001, 0.85 -- Vink et al. 2021 (no update neded) */

#define Z33 (0.85)         /* wind-metallicity dependence for H-rich stars with Teff: 25000-maxTempK */
                                /* 0.85 -- Vink et al. 2001, 0.42 -- Vink et al. 2021 (update neded) */

#define Flbv (1.5)           /* multiplication factor for LBV winds, LBV winds= Flbv*10^(-4) Msun/yr */
                                                           /* 1.0-standard;  USED ONLY WITH WIND=2,3 */
                                    /* range in Vink & de Koter 2002: 1.5*10^-4 -- 2.5*10^-6 Msun/yr */


#define AtmoRLOF 0                                                       /* 1 -- atmospheric RLOF on */
                                                                              /* 0 -- off (standard) */

#define WindRLOF 0                                      /* 1 -- wind RLOF (gravitationaly lensed) on */
                                                                              /* 0 -- off (standard) */



                                     /* calculated wind mass loss (Hurley or Vink) is multiplied by: */
#define WIND1 1.0               /* WIND1 for all H-rich stars (K=0,1,2,3,4,5,6) WIND1=1.0 - standard */
#define WIND2 1.0                     /* WIND2 for all He-rich stars (K=7,8,9) WIND2=1.0 - standard  */
                                 /* combination of the two allows for change of wind mass loss rates */

#define WRWINDS 1                    /* 1 -- standard: Hamman & Koesterke 1998 + Vink, de Koter 2005 */
                              /* 2 -- Zdziarski, Mikolajewska, Belczynski 2012 + Vink, de Koter 2005 */
                                                    /* 3 -- Zdziarski, Mikolajewska, Belczynski 2012 */
                            /* Hamman & Koesterke 1998: clumped winds based on models and 3 WR stars */
                             /* Zdziarski et al. 2012: clumped winds based on models and 34 WR stars */
                                           /* Vink, de Koter 2005: metalicty dependence for WR stars */
                                                                         /* works only for WIND=2!!! */

#define polar 0                                        /* 0 -- isotropic (direction) kicks: STANDARD */
                                                                           /* 1 -- fully polar kicks */



#define kick 6                          /* 1-Cordes&Chernoff kick distribution: NS and BH same kicks */
                      /* 2-BH kicks decreased, no direct collapse kicks otherwise like 1: (STANDARD) */
                                               /* 3-no BH kicks, otherwise like 2, 4-no kicks at all */
                                                 /* 5-Paczynski kicks distribution, otherwise like 2 */
                   /* 6-Hobbs 2005 single Maxwellian (Sigma3) distr: NS/BH kicks LOWERED by fallback */
               /* 7-Hobbs 2005 single Maxwellian (Sigma3) distr: NS/BH kicks NOT LOWERED by fallback */
             /* 8-Hobbs 2005 single Maxwellian (Sigma3) distr for NS, and lowered by BH mass for BHs */
                     /* nothing is LOWERED by fallback: Rodriguez 2016 uses this with Sigma3=265 kms */
                /* 9-no Sigma/Maxwellian used: based on Mns/Mbh and ejecta mass; random kick orient. */
                       /* used with alpha1 and beta1; based on Bray & Eldridge 2016, MNRAS 461, 3747 */
              /* 10-fall-back (mass) decreased kick with Sigma4 and full (neutrino) kick with Sigma5 */


#define Divide 0.4                                              /* Cordes&Chernoff98: 0.86,175.0,700 */
#define Sigma1 90.0                                      /* Arzoumanian01: 0.4,90.0,500.0 (STANDARD) */
#define Sigma2 500.0                                /* Ph.D.01: 0.8,175.0,700 (Divide,Sigma1,Sigma2) */
#define Sigma3 265.0                                        /* used only for Hobbs ditr.: 265 [km/s] */

#define alpha1 100.0                                                            /* [km/s] for kick=8 */
#define beta1 -170.0                                      /* in 2015 used to be alpha1=70, beta1=110 */
                                                                   /* in 2018: alpha=100, beta1=-170 */

#define Sigma4 265.0                   /* Maxwellian for mass (fallback-decreased) kicks: 265 [km/s] */
#define Sigma5 100.0                                 /* Maxwellian for neutrino kicks: 100 [km/s]??? */


#define kickback 0              /* 0 - (standard) fall back lowers kicks in entire range Mfb>0.0Msun */
                                              /*  for REMNANT=0 and over Mfb=0.2Msun for REMNANT=1,2 */
                               /* 1 - (optional) extends full kicks for BHs by use of the following: */
                             /* kickbackN: gives fall back mass Mfb over which kicks are lowered due */
                          /* to fall back in compactfN() functions. over that mass kicks are lowered */
                             /* proportionally to the frac (fractional fall back). so it means that: */
                         /* full kicks are applied below kickbackN and then decreased over that mass */
                              /* since frac may be already significant over kickbackN, kick decrease */
                                                         /* may start significant right away as well */

                                                                          /* used only if kickback=1 */
#define kickback0 0.6       /* [Msun] for REMNANT=0 kicks are lowered only after that fall back mass */
                         /* in standard mode=0.0, and should be in range ~ 0.0-1.0 for parm. studies */
#define kickback1 0.6       /* [Msun] for REMNANT=1 kicks are lowered only after that fall back mass */
                         /* in standard mode=0.2, and should be in range ~ 0.2-1.0 for parm. studies */
#define kickback2 0.6       /* [Msun] for REMNANT=2 kicks are lowered only after that fall back mass */
                         /* in standard mode=0.2, and should be in range ~ 0.2-1.0 for parm. studies */


#define newAcc1 1    /* 1 - (STANDARD): new He accretion onto CO/ONeMg WD from Piersanti et al. 2014 */
                                                                           /* Damian Kwiatkowski fit */
                   /* 0 - (OLD) scheme from Kato&Hachisu2004 (as described in Belczynski et al. 2008 */

#define newCONe 0        /* 0 - (STANDARD) original Hurley2000 formaule to get CO and ONeMg WD types */
                                                                  /* MASS LIMIT FOR 1.2MSUN FOR ZSUN */
                             /* 1 - uses WD mass to decide on WD type: limiting mass MCONe set below */

#define MCONe 1.0        /* [Msun], broad range includes 0.9-1.3Msun, 1.07Msun from Umeda et al.1999 */

#define newOut 0    /* 0 - (STANDARD) uses set mass (Msh) for double detonation SN1a of CO/ONeMg WD: */
                          /* 1 - mass dependent scheme which employs CONSTANT (onset of RLOF) system */
                                       /* parameters to calculate He shell mass at which WD explodes */
                         /* 2 - mass dependent scheme which employs INSTANTANEOUS mass transfer rate */
                                                  /* to calculate He shell mass at which WD explodes */
                      /* 3 - mass dependent scheme from Shen & BIldsten 2009 (ApJ 699, 1365), Fig. 5 */



#define Msh 0.1                      /* He shell mass required for double detonation SN1a on K=11/12 */
                                             /* typical values 0.05-0.2 Msun, with 0.1 most accepted */


#define MCh 1.44                                      /* Chandrasekhar mass [Msun] (WD/NS formation) */
#define Mmaxns 2.5                                         /* maximum mass of NS, over it we have BH */

#define SimpleNS 0                                        /* 0 -- standard: full NS mass calculation */
                            /* 1 -- simple calculation: NS mass increases proportionaly with CO core */
                                                     /* mass (or ZAMS mass) from Mns=1.1 to 2.5 Msun */
                                                              /* ECS NS masses are not affected !!!! */

#define ECS 1                                      /* 1- electron capture SN allowed; 0- not allowed */
                                 /* 1-AIC allowed (ONe WD->NS), 0-not allowed (ONe WD->SN Ia (K=15)) */
                            /* basically it is the same physics, so ECS and AIC are handled together */
                             /* although the range for ECS only may be changed with Mcbur1 parameter */

#define ECSlower 0.0                         /* kick as for regular SN, but mulitplied by that value */
                                                      /* make ECSlower=0 to have no kicks for ECS SN */
                     /* orbit change: mass loss from envelope + change of grav. mass + kick (if any) */

#define AIClower 0.0                         /* kick as for regular SN, but mulitplied by that value */
                                                    /* make AIClower=0.0 to have no kicks for ECS SN */
                                               /* orbit change: change of grav. mass + kick (if any) */
                                       /* BH AIC formation not affected, always no kicks for AIC BHs */

#define Mcbur1 1.83          /* He core mass at BAGB (or He star mass) over which ECS can take place */
                         /* 1.6 Msun is original Hurley2000 value, it corresponds to Mzams=6.31 Msun */
                                                                         /* Mzams=6.5 -> Mcbur1=1.66 */
                            /* Standard (to make ECS over Mzams>7.65 Msun): Mzams=7.0 -> Mcbur1=1.83 */
                                               /* Mzams=7.5 -> Mcbur1=2.0, Mzams=8.03 -> Mcbur1=2.17 */

#define Mcbur2 2.25         /* He core mass at BAGB (or He star mass) below which ECS can take place */
                                            /* Standard Mcbur2=2.25: corresponds to Mzams=8.275 Msun */
                                      /* Mzams->Mcbur2: 11.0->3.24, 10.0->2.87, 9.0->2.51, 8.0->2.16 */


#define BHSPIN 4                                        /* 1,2,3 - evolution of BH spins ON; 0 - OFF */
#define aspin_init 0.0                /* BHSPIN=1: initial BH spin independent of BH mass=aspin_init */
                                                /* BHSPIN=2: initial BH spin depends on CO core mass */
                                /* BHSPIN=3: initial BH spin from Geneva code: with 40% initial Vrot */
                  /* BHSPIN=4: initial BH spin from MESA code: magnetic dynamo with 40% initial Vrot */
                             /* BHSPIN=10: initial BH spin from Jim Fuller: a_natal=0.01 for all BHs */


#define BHSOUT 0                            /* 0/1 -- no/output for spins of BH-BH: bhbh.spin.NN.dat */
                    /* basic output for spins in compact.dat, more detailed turned on by BHSOUT here */


                                                                           /* IP paremeter settings */
#define IPfrac 0.2        /* fraction of RLOF systems with WD accretors which are IPs, rest are CVs */
#define etaIP 0.95       /* frac. of transferred material reaching WD [0.9-1.0]; rest is blown away */
#define etaBE 1.0          /* frac. of gravitational binding energy converted into radiative energy */
#define etaGEO 1.0           /* frac. of total IP X-ray luminosity that is emitted towards observer */
#define etaXIP 0.09                 /* fraction of total bolometric luminosity that is converted to */
                             /* Chandra X-ray luminosity (2-8 keV) for IPs: cf Norton & Watson 1989 */


                                                                           /* CV parameter settings */
#define etaXCV 0.5                  /* fraction of total bolometric luminosity that is converted to */
                        /* X-ray lum. for CVs: 0.5: Chris B. guess: supply something realistic here */



#define Mhecon 3.0           /* below this mass Helium giant (K=9) has convec. env: MB, Conv. Tides */
                                                    /*  otherwise it is radiative: No MB, Rad Tides */
                                            /* K=7,8 -- always radiative envelope: No MB, Rad Tides */

#define RotF 1                                           /* initial rotation of stars: choice of law */
                                      /* 1- STANDARD (A.Zezas formula), 2-Hurley et al. 2000 formula */

#define Xdtstep 1       /* 1 -- standard: timestep for Xray output dt=1 Myr (DetWIND.dat, DetMT.dat) */
                                          /* 0 -- output at every timestep the code (binary.c) takes */


#define XEdd 2                       /* assumed upper limit on MT accretion rate onto compact object */
              /* 2: (standard) -- weak BH/NS RLOF/wind accretion limited to ln() over-Eddington rate */
                                                  /* WD limited to Eddington rate: Samaresh+JP model */
                /* 1: (standard till 2018) BH accretion limited to Ohsuga 2006 critical rate (Mcrit) */
                                            /* 1: WD/NS limitied to Eddington rate (Medd ~ 10*Mcrit) */
                                           /* 0: BH accretes all mass, NS accretes all mass if TRA=1 */
                                                         /* X-ray Luminosity is limited to Eddington */


#define PT 1              /* 1 -- system either persistent or transient, judged on system parameters */
                                                                      /* 0 -- all systems persistent */
                                              /* PT: used only in old Xbin(), but not in new Xbin1() */

#define TRA 1   /* 0 -- transient system: NS do not accumulate any material for given accretion rate */
                                   /* 1 -- persistent system: NS accumulates all meterial (standard) */
                                                                   /* BH always accrets all material */


#define CR1 2.0                               /* MT diagnostic diagram: critical value between q1-q2 */
#define CR2 1.2                         	     /* CR2 -- H-rich giant donors, CR1-- all others */


#define MB 3                    /* 0-- no magnetic breaking, 1-- new MB (Sills: saturation included) */
                    /* 2-- old MB Rapapport, 3-new MB: Ivanova & Taam 2003 (STANDARD, sat. included) */
#define gamMB 2                                      /* only for MB=2, exponent in the old law [0-4] */


#define Tides 1                                                            /* 1-full tidal evolution */
                     /* 0-no tidal evolution, no pre MS synch/circ, but spins followed to include MB */

#define PMStid 4.3       /* if period of ZAMS binary is shorter than PMStid [days], then system has: */
                                          /* circularized and possibly comp. synchronized, 4.3 days: */
                             /* Mathieu et al. 1992, Binaries as Tracers of Stellar Formation, p.278 */

#define Caltid 50.0                   /* calibration of tidal/radial dump. synch/circ constant (k/T) */
#define Calrad 1.0                             /* 1--Hurley 2002, 10--standard, 100--parameter study */
                              /* !!! CHEAT: if e<0.001 Caltid=Calrad=1 to allow ODE's solve dorbdt() */

#define EE2 2                   /* 1 -- OLD: tidial damping constant E2 from Hurley et al 2002 eq.43 */
                        /* 2 -- NEW (STANDARD): the constant smaller than for OLD (so smaller tides) */
                                    /* new E2: simple fit by eye to Claret 2007 data (A&A 467, 1389) */
              /* radiative tides very small in each case: no circ/sync expected for massive binaries */

#define SymdM 10.0                       /* enhancement factor in wind-accretion rate for symbiotics */


#define HCE 3                                 /* 0- NS in any CE do accrete with hyper critical rate */
                                           /* 1-NS in CE with H-rich stars do not accrete ANY matter */
                                                      /* 2-NS/BH do not accrete ANY matter in ANY CE */
                                              /* 3 -- NS and BH accrete as indicated below: STANDARD */
                        /* orb. sep. changes through std. energy balance (not through B-H formalism) */
#define NSACC 1                /* NSACC=1: (STANDARD) NS accrets fraction=NSBOND of Bondi-Hoyle rate */
                                                   /* BH accrets fraction=BHBOND of Bondi-Hoyle rate */
#define BHBOND 0.05                                 /* for BH accretion: 0.05 standard, range: [0-1] */
#define NSBOND 0.05                                 /* for NS accretion: 0.05 standard, range: [0-1] */

#define ADDM1 0.05                                                 /* used only if HCE=3 and NSACC=0 */
#define ADDM2 0.1                 /* 0.05-0.1 Msun the lowest possible to recycle accreting NS in CE */
                                                 /*  NS accrets from range ADDM1-ADDM2 (unif. dist.) */
                                             /* BH still accrets fraction=BHBOND of Bondi-Hoyle rate */


#define TRwd 1             /* 1-use traping radius condition to check on merger before going into CE */
                  /* 0-do not use the condition for systems with WD donors (tight WD-WD/NS/BH pairs) */
                              /* for all the other type of systems, the condition is always applied! */

#define CE2 0                                             /* 0-CE FROM K=2 donor enabled, 1-disabled */
#define CE3 0                                             /* 0-CE FROM K=3 donor enabled, 1-disabled */
#define CE4 0                                             /* 0-CE FROM K=4 donor enabled, 1-disabled */
                                                                    /* (1 -- always means merger!!!) */


#define CE 3                                         /* 1- always standard energy alpha prescription */
                                                    /* 2- always Nelemans ang.momentum gamma prescr. */
                                                              /* 3 - mixed alpha/gamma prescription: */
                                           /* do gamma for non-remnant donors: K=0,1,2,3,4,5,6,7,8,9 */
                                                /* do alpha for remnant donors: 10,11,12,13,14,16,17 */
                                                /* CE=3: use as standard for SNIa and WD simulations */
                                                     /* CE=1,2,3 -- old CE/RLOF criteria: year <2021 */
                                           /* CE=4 -- new CE/RLOF criteria: Pavlovski2017+Olejak2021 */
                                           /* CE always done with standard energy alpha prescription */

#define Zrlof 0.01          /* used only with CE=4, division between high- and low-metallcity donors */
                                               /* [absolute metallicity] change in range: 0.01-0.005 */


#define Gamma 1.5        /* Nelemans CE efficiency, used if CE=2: no hyper-crit acc. incorporated!!! */
#define Alfa 1.0                  /* efficiency of orbital energy loss for CE ejection, used if CE=1 */
#define Beta 1.0        /* Podsiadlowski et al. 1992 -specific angular momentum of matter [2Pia^2/P] */
#define Fa 0.5            /* Meurs,Heuvel 1989 -fraction of mass lost by donor attached to companion */
#define Cd 6.0                   /* drag coefficent, see Bethe and Brown 1998 ApJ 506, 780 eq. (5.7) */

                         /* lambda in CE energy prescription: sets binding energy of donor ewnvelope */
#define chlambda 1            /* 0-(old standard) use single value of lambda=Lambda for H-rich stars */
                             /* 1-(standard) average value with and without internal energy, chinese */
                                           /* 2-(optional) lambda with full internal energy, chinese */
                                             /* 3-(optional) lambda with no internal energy, chinese */
                              /* 1-3: Chinese: X.-J. Xu and X.-D. Li arXiv:1004.4957 (v1, 28Apr2010) */
                                                    /* 0-3: for all He-rich stars always Natasha fit */

#define Lambda 1.0                        /* this is used as single value of lambda for H-rich stars */

#define golambda 1.0    /* increase fudge factor on lambda based on enthalpy arguments: golambda=2-5 */
                          /* Ivanova & Chaichenets 2010, ApJ 731, L36; for typical calculation use 3 */


                                              /* Delayed dynamical instability critical mass ratios: */
#define DDIcrit1 3.0             /* if Mdonor(K=0,1,2,3,4,5,6)>=DDIcrit1*Maccretor -- DDI will start */
#define DDIcrit2 1.7                         /* if Mdonor(K=7)>=DDIcrit2*Maccretor -- DDI will start */
#define DDIcrit3 3.5                       /* if Mdonor(K=8,9)>=DDIcrit3*Maccretor -- DDI will start */


/* -------------------------------------- technical settings --------------------------------------- */

#define acc 1.0e-10                                             /* accuracy for variable comparisons */

                                      /* time steps for calculation at different evolutionary stages */
#define delms 0.01                                                /* Main Sequence, HS Main Sequence */
#define delhg 0.01                                                                 /* Hartzprung Gap */
#define delrg 0.01                           /* Red Giant Branch, HS Hartzprung Gap, HS Giant Branch */
#define delhe 0.01                                                            /* Core Heluim Burning */
#define delag 0.01                                                        /* Asymptotic Giant Branch */

#define maxdM 0.01                     /* 0.01: Maximum allowed change of mass (0.01: max change 1%) */
#define maxdR 0.01               /* 0.1-0.01: Maximum allowed change of radius (0.1: max change 10%) */
                                        /* at some cases radius may change more than that: not often */

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))


#define MM 1            /* 1 --singl.c fills in tend, and may return at diffrent time than requested */
                                          /* by external driver (binary.c), 0--always return at tend */
                          /* even if MM=1, during MT the single will return always at requested time */

#define Random 1                   /* 1-- uses ran3.c function to generate random numbers (standard) */
                                              /* 2-- uses ran2.c function to generate random numbers */


/* -------------------------------------- global variables/parameters ------------------------------- */

FILE St *fp0;                       /* file with errors and warnings, common for singl.c and binary.c */
FILE St *fp108;                                                                 /* data on symbiotics */
FILE St *fp130,*fp131;                                          /* detailed output for X-ray binaries */
FILE St *fp150;                                                              /* data on WDWD binaries */
FILE St *fp151;                                                                     /* output for CVs */
FILE St *fp160;                                              /* detailed output for systems with a BH */
FILE St *fp170;                                                /* output for NS-WD and BH-WD binaries */
FILE St *fp180;                                                                   /* output for SN Ia */
FILE St *fp181;                                                                  /* output for AM CVn */
FILE St *fp190;                                                                /* output for OB stars */
FILE St *fp200;                                                        /* output for compact binaries */
FILE St *fp300;                                                        /* output for cluster binaries */
FILE St *fp400;                                                               /* output for long grbs */
FILE St *fp500;                                   /* output (SN/compact) for single stars: sstars.dat */
FILE St *fp510;                                   /* output (SN/compact) for binary stars: bstars.dat */
FILE St *fp600;                                                          /* output for WD-WD binaries */
FILE St *fp700;                                                          /* output for RR Lyrae stars */
FILE St *fp800;                                                        /* output for NS/BH CE mergers */
FILE St *fp900;                                           /* output for CE events: transient project  */
FILE St *fp1000;


#define tr1 (3.06596e+14)                         /* tr1*Lx -change Lx[Rsun^2 Msun/Myr^3] to Lx[erg/s] */
#define tr2 (0.12414)                                         /* tr2*p -change p[km/s] to p[Rsun/day] */
#define tr4 (1.2414e-06)                                      /* tr4*p -change p[cm/s] to p[Rsun/day] */
#define Pi 3.1415926536
#define Rsun 6.9599e+10                                              /* [cm], sun radius in cgs units */
#define Msun 1.989e+33                                                  /* [g], sun mass in cgs units */
#define Lsun 3.826e+33                                                                     /* [erg/s] */
#define G 2938.408721                                 /* [Rsun^(3) Msun^(-1) day^(-2)] grav. constant */
#define GG (2938.408721*pow(365.25,2.0))                              /* [Rsun^(3) Msun^(-1) yr^(-2)] */
#define GGG (2938.408721*pow(1000000.0*365.25,2.0))                  /* [Rsun^(3) Msun^(-1) Myr^(-2)] */
#define GGGG (6.6743e-08)                                          /* cgs units: cm^(3) g^(-1) s^(-2) */
#define CCC (13593198857139.902344)                                /* [Rsun Myr^(-1)]: speed of light */


double St xPacz[2000],yPacz[2000];                        /* for Paczynski kicks filled in get_Vkick5 */


double St *xp,**yp,dxsav;                                                        /* for comm_env3d1() */
int St kmax,kount;                                                               /* for comm_env3d1() */

double St minMB;                              /* min. mass of MS star over which conv. env. developes */
double St maxMB;            /* max. mass of MS star below which  conv. env. developes (depends on ZZ) */
                                    /* MS star: below minMB fully convective, over maxMB no conv env. */


long int St idum;                                         /* mark random number start of given system */
int St iidd_old;

char St evroute[128];                        /* string for storing the evolutionary route of a system */
int St nevroute;                                                      /* number of entries to evroute */



/* -------------------------------------------------------------------------------------------------- */

#define ZZsun 0.0142                             /* solar metallicity: e.g., used in wind2(), Rzamsf() */
                                               /* do not touch this until you know what you are doing */

double St M_hook,M_HeF,M_FGB;                                         /* critical evolutionary masses */

double St aa1,aa2,aa3,aa4,aa5,aa6,aa7,aa8,aa9,aa10;                        /* ZZ dependent parameters */
double St aa11,aa12,aa13,aa14,aa15,aa16,aa17,aa18,aa19,aa20;
double St aa21,aa22,aa23,aa24,aa25,aa26,aa27,aa28,aa29,aa30;
double St aa31,aa32,aa33,aa34,aa35,aa36,aa37,aa38,aa39,aa40;
double St aa41,aa42,aa43,aa44,aa45,aa46,aa47,aa48,aa49,aa50;
double St aa51,aa52,aa53,aa54,aa55,aa56,aa57,aa58,aa59,aa60;
double St aa61,aa62,aa63,aa64,aa65,aa66,aa67,aa68,aa69,aa70;
double St aa71,aa72,aa73,aa74,aa75,aa76,aa77,aa78,aa79,aa80,aa81;

double St bb1,bb2,bb3,bb4,bb5,bb6,bb7,bb9,bb10;
double St bb11,bb12,bb13,bb14,bb15,bb16,bb17,bb18,bb19,bb20;
double St bb21,bb22,bb23,bb24,bb25,bb26,bb27,bb28,bb29,bb30;
double St bb31,bb32,bb33,bb34,bb36,bb37,bb38,bb39,bb40;
double St bb41,bb42,bb43,bb44,bb45,bb46,bb47,bb48,bb49;
double St bb50,bb51,bb52,bb53,bb54,bb55,bb56,bb57;
/* bb50 set in Ragbf() called by single() not in coeff_bb() */


/* counters for different SN types and for fraction of NS formed in different SNs */
int St ss0,ss1a,ss1b,ss1c,ss2a,ss2b,ss2c,ss3a,ss3b,ss3c,ss4a,ss4b,ss4c;
int St ns1a,ns1b,ns1c,ns2a,ns2b,ns2c,ns3a,ns3b,ns3c,ns4a,ns4b,ns4c;

int St st0,st1a,st1b,st1c,st2a,st2b,st2c,st3a,st3b,st3c,st4a,st4b,st4c;
int St nt1a,nt1b,nt1c,nt2a,nt2b,nt2c,nt3a,nt3b,nt3c,nt4a,nt4b,nt4c;

int St sII[20];                                                 /* number of SN II in given Gyr */
int St sIbc[20];                                               /* number of SN Ibc in given Gyr */
int St sFa[20];                           /* number of direct core collapse events in given Gyr */


int St hce1,hce2,hce3;                                /* number of Hyper Common Envelope events */
int St hehe1;           /* number of MT of two He giants when both are more massive then Mhecon */
int St hehe2;           /* number of MT of two He giants when both are less massive then Mhecon */
int St hehe3;   /* number of MT of two He giants when donor is not massive and acceptor massive */
int St hehe4;   /* number of MT of two He giants when donor is massive and acceptor not massive */




void singl(double *MzamsR, double *M0R, double *MR, int *KR, double *TBR, double *TVIRR, double *TER,
           double *LR, double *RR, double *MCR, double *MHER, double *MC0R, int *FLAGR, double *DTR,
           double *MPRER, int *KPR, double *TSTART, double *FRACR, double DMMTR, int DERIVR,
           double MOUTR1, double MOUTR2, int *ECSSNR, int *FBR, int SNR, double FMS);



void coef_aa(void);                                              /* Initializing functions */
void coef_bb(void);
double M_hookf(void);
double M_HeFf(void);
double M_FGBf(void);
void init_SN(void);


void ZAMSf(double *M, double *t, double *Mhe, double *Mco, int *stop);          /* Evolutionary functions */
int MSf(double *M, double *t, double *Mhe, double *Mco, int *stop);
int HGf(double *M, double *t, double *tstart, double *Mhe, double *Mco, int *Kp, int *stop);
int RGf(double *M, double *t, double *tstart, double *Mhe, double *Mco, int *Kp, int *stop);
int HEf(double *M, double *t, double *tstart, double *Mhe, double *Mco, int *stop);
int AGBf(double *M, double *t, double *tstart, double *Mhe, double *Mco, int *Kp, int *stop, int Klast);
int HSMSf(double *M, double *t, double *tstart, double *Mhe, double *Mco, int *stop);
int HSGBf(double *M, double *t, double tstart, double *Mhe, double *Mco, int *Kp, int *stop);
int REMNANTf(double *M, double *t, double *Mhe, double *Mco, int Kp);
void sn_type(double M, double Mhe, double Mco, double t, int K);
double windf1(double M, double L, double R, double Mc, int K);
double windf2(double M, double L, double R, double Mc, int K);
double windf3(double M, double L, double R, double Mc, int K);
int perturb(double M, double Mc, double *L, double *R, double t, int K);
double tnf(double M, int K);

double Lzamsf(double M);                                         /* Main Sequance functions */
double Rzamsf(double M);
double thookf(double M);
double tmsf(double M);
double Lmsf(double M, double t);
double Rmsf(double M, double t);
double Ltmsf(double M);
double Rtmsf(double M);

double thgf(double M);                                           /* Hartzprung Gap functions */
double Mchgf(double M, double t);
double Lhgf(double M, double t);
double Rhgf(double M, double Mtmp, double t);
double Mcehgf(double M);
double Lehgf(double M);
double Rehgf(double M, double Mtmp);

double tbgbf(double M);                                          /* Red Giant Branch functions */
double Mcbgbf1a(double M);                                       /* real core mass on BRGB */
double Mcbgbf1b(double M);                                       /* dummy core mass on BRGB */
double Lbgbf(double M);
double Ahpf(double M);
double Bf(double M);
double Df(double M);
double pf(double M);
double qf(double M);
double tinf1f(double M);
double tinf2f(double M);
double txf(double M);
double Mxf(double M);
double Lxf(double M);
double Mcgbf1a(double M, double t);                              /* real core mass on RGB */
double Mcgbf1b(double M, double t);                              /* dummy core mass on RGB */
double Mcgbf2(double M, double L);
double Lgbf1(double M, double t);
double Lgbf2(double M, double Mc);
double Rgbf(double M, double L);

double Rhe1f(double M);                                          /* Core Helium Burning Phase functions */
double Mche1f(double M);
double Lhe1f(double M);
double the1f(double M);
double thef(double M);
double Mchef(double M, double t);
double Lminhef(double M);
double Lzahbf(double M, double Mc);
double Rzahbf(double M, double Mc);
double Rmhef(double M);
double Rmhefa(double M, double Lzahb);
double fblf(double M);
double taublf(double M);
int CHeBf(double M, double t, double *R, double *L, double *Mc);

double tbagbf(double M);                                         /* Asymptotic Giant Branch functions */
double Mcbagbf(double M);
double Lbagbf(double M);
double Rbagbf(double M);
double Ahef(void);
double tinf1ef(double M);
double tinf2ef(double M);
double txef(double M);
double Mceagbf1(double M, double t);
double Mceagbf2(double M, double L);
double Leagbf1(double M, double t);
double Leagbf2(double M, double Mcco);
double Mcduf(double M);
double Lduf(double M);
double tduf(double M);
double Ahhef(void);
double tinf1tf(double M);
double tinf2tf(double M);
double txtf(double M);
double Mctagbf1a(double M, double t);                            /* real core mass on TP AGB */
double Mctagbf1b(double M, double t);                            /* dummy core mass on TP AGB */
double Mctagbf2(double M, double L);
double Ltagbf1(double M, double t);
double Ltagbf2(double M, double Mc);
double lambdaf(double M);
double Ragbf(double M, double L);

double Lzhsf(double M);                                          /* Naked Helium Stars functions */
double Rzhsf(double M);
double thsmsf(double M);
double Lhsmsf(double M, double t);
double Rhsmsf(double M, double t);
double Lthsf(double M);
double Bhsf(void);
double Dhsf(double M);
double phsf(void);
double qhsf(void);
double tinf1hsf(double M);
double tinf2hsf(double M);
double txhsf(double M);
double Lxhsf(double M);
double Mxhsf(double M);
double thsevolf(double M, double Mc);
double Mchsgbf1(double M, double t);
double Mchsgbf2(double M, double L);
double Lhsgbf(double M, double Mc);
double Rhsgbf(double M, double L, double Lths, double Rzhs, int *K);
double Mcmaxf(double M);

double compactf0(double M, double Mco);
double compactf1(double M, double Mco);
double compactf2(double M, double Mco);
double compactf3(double M, double *Mco, double *Mhe);
double compactf4(double M, double *Mco, double *Mhe);
double compactf5(double M, double *Mco, double *Mhe);
double compactf6(double M, double *Mco, double *Mhe);
double compactf7(double M, double *Mco, double *Mhe);
double compactf8(double M, double *Mco, double *Mhe);
double compactf9(double M, double *Mco, double *Mhe);
double compactf10(double M, double *Mco, double *Mhe);
double compactf12(double M, double *Mco, double *Mhe);
double compactf13(double M, double *Mco, double *Mhe);
double compactf14(double M, double *Mco, double *Mhe);

double Lwdf(double M, double t, int K);
double Rwdf(double M, int K);
double Rnsf(void);
double Rnsf1(double M);
double Rbhf(double M);
double Eiscof(double R, double M, double aspin);
double Riscof(double M, double aspin);
double Mcsnf(double Mcbagb, double Mcco0);
double Mupf(void);
double Mecf(void);

double get_T(double L, double R);                                /* Other functions */
double inter_line(double x1, double y1, double x2, double y2, double x);

                                                                 /* functions for binary.c */
double tnfI(double M, double m0, int K);
double get_t(double Min, double Max);
