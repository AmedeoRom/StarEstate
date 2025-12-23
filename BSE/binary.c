/* software by chris belczynski: last modified on Nov 21, 2022 */
/* program for single and binary evolution */
/* use input from singl.c: single evolution by Hurlay et al. 2000 as implemented by KB 2000/01 */
/* need a header file sinbin.h, where all parameters are set */
/* change Dec02, 2003: error in dMmtf() for CE2=1 both MT and CE systems were assumed mergers -- corrected: */
/* now for CE2=1 only CE systems assumed to be mergers: ultracompact/GC faint/ngc1569 papers may be affected */
/* change Dec15, 2003: Tauris&Takens 1998 speeds for disrupted components at SN (eccentric orbits added) */
/* change Jan04, 2004: beta_wind changed in Xbin(): used to be always 0.125, now can be bigger (upto 7.0) */
/* for most of the stars -> smaller wind acc. rates, smaller Lx for wind systems (upto order of mag.) */
/* change Jan04, 2004: all potential X-ray binaries send to Det* files, no Lxcrit limit imposed */
/* now: make sure when using Det* files to check Lx of each source!!! and impose Lxcrit to your liking */
/* change Jan09, 2004: in dMgainf() K=16 accretor added, before the error messege may have been generated */
/* change Jan09, 2004: AIC full kick option itroduced, see sinbin.h for AICkick parameter */
/* change Jan10, 2004: some clean up, additions on msp counting */
/* change Feb25, 2004: after SN the stellar radius is allowed to exceed lobe radius twice, before */
/* calling it after-SN collision (RF>=Rl -> RF>=2.0*Rl), important change if AIC=1! */
/* change March02, 2004: H accretion onto CO/ONeMg WD changed (lowered slightly) dMgainf() changed */
/* change March03, 2004: velocities in ALL output files in [Rsun/day] */
/* change March08, 2004: output files modified + some code cleaning + all output times (+Tst) */
/* change March09, 2004: no-tides option added, "ultra.dat" removed */
/* change March15, 2004: donor K=16 added in dMtranf() */
/* change March15, 2004: inbinA, inbinB (binary/single flag) added in "sb.dat" */
/* change April06, 2004: CE2,CE3,CE4 properly accounted for now. */
/* change May04, 2004: tidial term added to MT rate calculation */
/* change May04, 2004: action for K=4 in MT and MB calc. changed: M<7.0 -- conv env. and M>=7.0 -- rad env. */
/* change May04, 2004: MT calculation procedure unified for all non-degenerate stars */
/* change May04, 2004: file "mtprob.dat" added to store problematic MT systems */
/* change May12, 2004: in dMgainf() Eddington limit imposed on accretor K=16 */
/* -------------------------------------------------------------------------------------------------------- */
/* change Nov15, 2004: IMFs changed to one unified get_M() function, see new IMF settings in sinbin.h */
/* change Nov16, 2004: LISA "wdwd.dat" output file added (for Ashley) wd-wd binaries only */
/* change Nov17, 2004: BH output added: evol. for each system with BH recorded in "bhsys.dat" file (for Olek) */
/* change Nov17, 2004: New Magnetic Braking added:Ivanova & Taam 2003 (ApJ 599, 516) */
/* Nata supplied the lines of code for tmbf(): used as a STANDARD from now on (MB=3 in sinbin.h) */
/* change Nov18, 2004: AIC marker added to all output files in evolutionary history */
/* both AIC of WD->NS and AIC of NS->BH are marked and recorded */
/* change Dec01, 2004: Bolometric Luminosities for WDs added */
/* change Dec03, 2004: He accumulation on He WD changed, Mout2 added (both here and in singl() */
/* change Dec13, 2004: LISA output extended to all compacts: wd/ns/bh-wd/ns/bh: "lisa.dat" (instead wdwd.dat) */
/* change Dec15, 2004: Spectra output file for Bob "spectra.dat" added */
/* change Jan13, 2005: starsin.c/singlS.c/starsin.h added to the package: population synthesis of single stars */
/* change Jan14, 2005: unexpected R jump (K=3->4) taken care of (all sin*.c files on ~/Working dir. corrected) */
/* it was causing MT to stop when K=3->4 and the radius would sometimes jump by factor 2-4 leading to a merger */
/* correction: change of M->M0 in line: (*L)=Lzahb*pow(Lhe1f(M->M0)/Lzahb,pow(tau,3.0)); in CHeBf() func. */
/* change Jan27, 2005: new way to calculate velocities of SN disrupted stars: from T.Bulik: sinvel() */
/* change Feb01, 2005: K=17 hybrid He/CO WD added */
/* change Feb02, 2005: SN Ia added: "sn1a.dat" */
/* change Feb04, 2005: G.Nelemans CE gamma prescription added */
/* change Feb10, 2005: output for BHs updated/corrected (e.g., no disrupted BH systems in bhsys.dat) */
/* change Feb21, 2005: orbital velocity corrected (increased) in Xbin(): possible decrease of acc_rate/Lx */
/*        correction: v_orbi2=GG*Mb/a -> v_orbi2=GG*(Ma+Mb)/a   where Ma -- donor mass */
/* change Feb24, 2005: AIC on/off switch introduced (see AIC in sinbin.h) */
/* change Mar23, 2005: singl.c in AGBf() -- proper management of while() loop when in R jump added */
/* change Mar23, 2005: binary.c -- systems with bad solution in orb_change() (e.g., spin vel negative): */
/*        skipped and counted and reported in info.dat */
/* change May26, 2005: NS single formation mass added */
/* change May26, 2005: OB stars added */
/* change Jun03, 2005: sub_MCh outburst of CO & NeMg WDs in case of He accre. was overestimated in dMgainf() */
/*                     now proper switch was introduced, when dM goes to higher rates -- reported by Natasha */
/* change Jun06, 2005: std MB corrected: ke3=8.06e-39 -> kw3=619.2 (42 orders of magnitude change!!!) */
/*                     all systems with MB important for evolution are shite .... before this date  */
/* change Jun07, 2005: error for dirupted binaries: SN object was reciving companion velocity and vice versa: */
/*                     in sinvel() in last lines Vsa and Vsb interchanged, and now all is OK (Bulo+Ash) */
/* change Jun18, 2005: in dMgainf() K=11 donor to K=10,17 accretor introduced with eta=100% acc. efficiency */
/*                     possible to form such system with q-reversal (Kdon=11)+stripping of K=8 star (Kacc=17) */
/* change Jun18, 2005: in dMgainf() change of K=10,17 to He star (K=7) allowed only over Macc=0.35 Msun */
/* change Jul28, 2005: small change in sinvel() -- may have been giving problems for systems (e=o, no kick) */
/* change Aug02, 2005: a0,e0 added to lisa.dat file output */
/* change Aug02, 2005: Xray function for CVs added for ash: Xcvf(), output goes to cv.dat */
/* change Sep09, 2005: theoretical alpha=0.1,0.005,0.02 for He,CO,O (Menou Fig.2) was used to get dMtran */
/*                     and it was used wrongly: C=0.1/alpha=1,5,20 => C^0.42=1.0,1.97,3.52 */
/*                     (instead of C=alpha/0.1=1,0.05,0.2): see function dMtranf() for details */
/*                     now we use always alpha=0.1: observational evidence (Menou, S4, first parag.) */
/*                     so the values are C=1.0 in each case, and C^0.42=1.0 => which results in: */
/*                     smaller dMtran by factor: ~1,2,3.5 for He, CO and ONe donors: less transients */
/* change Oct17, 2005: file compact.dat added for compact object binaries: NS/BH-NS/BH */
/* change Nov03, 2005: in dMmtf() tmb_a+tmb_b used instead of only tmb (only for star a) in MT rate calcul */
/*                     same for: tti_a+tti_b  used instead of only tmb (only for star a) in MT rate calcul */
/*                     dMmtf() three more arguments passed: wb,Ib,KTb */
/* change Nov04, 2005: dMgainf() updated for He accumulation on CO WDs (Kato & Hachisu 2004) */
/* change Nov18, 2005: critical Edd. accretion rate changed: epsilon introduced: in dMeddf() */
/*                     crit.acc.Edd rates are now higher by factor 2! changes evol. for some RLOF-BH-XRBs */
/* change Nov20, 2005: Mmaxconvf(), Mminconvf() introduced: max./min MS mass to develope conv. env.- depends on Z */
/* change Nov22, 2005: Mhecon value and usage changed!!! important for all He star evolution -- now: */
/*                     1) convective envelopes only for (K=9 && M<Mhecon=3.0), and it was for all K=8,9 */
/*                     2) Mhecon used for Tides now (but not for MB) and not for CE calcuations (old) */
/* change Dec02, 2005: Hobbs 2005 kicks added: get_Vkick6() function added (now standard kick distr.) */
/* change Dec04, 2005: output for compact.dat changed: sinbin.h: Mminsec added; binary.c output changed, also: */
/*                     dt after SN kept small (was big): sn!=1 added in: else if(Ka>=10 && Kb>=10 && sn!=1) */
/* change Dec05, 2005: lambda dependance (on radius) added for He stars: lambda [0.6-0.004] so smaller Af: */
/*                     more systems will merge during CE now: see lamf() */
/* change Mar20, 2006: sqrt() in eq.(21) with radial damping added: see KTf() -- now more efficient tides for */
/*                     masive stars, Hurley 2002, eq.(42) was missing sqrt() over (MR^2/a^5) and so did mine */
/* change Apr13, 2006: diffrent limits on con/rad envelopes for K=2 and K=4 stars: number of func. changed: */
/*                     Renv_con(),Menv_con(),KTf(),tau_tidef(),tmbf(),dorbdt(),odeint1(),rkqs1(),rkck1(),orb_change() */
/*                     so i can pass luminosity and radius (-> Teff) now needed to get limit of con/rad envs */
/* change Apr21, 2006: Electron capture SNe introduced; in singl.c: AGBF(), HSGBf(), REMNANTf() changed */
/*                     in sinbin.h: ECSlower, ECS added, in binary.c: singl() takes one more paramet., also explode(), */
/*                     copySin() changed, ecssna, ecssnb added */
/* change Apr21, 2006: AIC changed: kick handling (see AIClower in sinbin.h and explode()), mass treshold (MCh->Mecs) */
/*                     in singl(), and effect of mass loss (10% in binding energy) accounted for on binary orbit */
/* change May17, 2006: choice of massive binaries from IMF corrected: it was skewd towards higher q */
/* change May17, 2006: inbinA,inbinB were declared as doubles, now properly integers, so they contained wrong info */
/* change May17, 2006: system param. stored at last RLOF timestep for record in compact.dat */
/* change May25, 2006: corrections/improvements to ECS/AIC handling */
/* change May27, 2006: NS mass spectrum new (FeNi core masses adopted from Timmes et al 1996 in REMNANTf() */
/* change Jun09, 2006: in compactf() 4.823 changed to 4.280 (typo): final CO core mass for Mzams=18 msun */
/* change Jun27, 2006: BBeta changed (Mb removed from the expression) in dMmtf() */
/* change Jun27, 2006: tidal term removed from Mdot calculation: no change since no tides on circ/sync orbit */
/*                     as forced during ongoing MT */
/* change Aug17, 2006: in compact.dat output for initial Ma,Mb,a,e extended to %f.12 */
/* change Oct25, 2006: in singl.c in func compactf() one segment added to eq. calculating Mfeni: */
/*                     just to make it continous, no big change */
/* change Nov08, 2006: there was still CE allowed for K=2 (even if CE2 was set to 1) in case that both stars */
/*                     were over their Roche lobes: double CE. in output files look for string "CE12(X-X; and */
/*                     eliminate the systems that have at least one X=2, the most common are CE12(2-2; CE12(4-2; */
/* change Dec07, 2006: model added for twin binaries: set by MarkRat=4 and qdivide in sinbin.h */
/* change Jan06, 2007: new accretion model added for CE phase: mass accreted from range ADDM1-ADDM2 (unif. dist.) */
/*                     it is new standard prescription: HCE=3  with ADDM1=0.05 and ADDM2=0.1 */
/*                     change in comm_env3d1() and comm_env3b1(): used both in energy and Nelemans CE pres. */
/* change Jan06, 2007: if HCE=3 and CE=2-> Nelemans CE pres. used for all cases (NS/BH comp. included): comm_env3b1() */
/*                     before i was using in case of NS/BH comp. energy balance with B-H accretion: comm_env3d1() */
/* change Feb07, 2007: BH spin evolution added: project with Ron: spin_evol() added */
/* change Mar27, 2007: polar kicks added: explode() */
/* change Apr12, 2007: Simple final masses for NS added: see SimpleNS parameter in sinbin.h */
/* change Apr28, 2007: correction in singl.c for rejuvanation added -- now stars fully rejuvanated: */
/*                     change for con. core MS stars (M<0.35 || M>1.25) and all He MS stars - they live now longer on MS*/
/* change Jul18, 2007: modification of BH spin evolution to include new Ohsuga 2006 paper and new critical Mdot rate */
/* change Oct12, 2007: addition of qratio for TWIN binaries */
/* change Nov16, 2007: BH spin revisions: less mass on (0.5->0.1 Bondi acc) BH; 0.05-0.1Msun -> 0.1 Bondi acc for NS */
/* change Dec13, 2007: updates for H accretion onto CO WD in dMgainf() from Nomoto et al. 2007: ApJ 663, 1269 */
/* change Dec14, 2007: cluster addition: cluster.dat added for project on M67: set CLUSTER to 1 to run it */
/* change Jan21, 2008: output to fp210 (BH spin data) corrected: now a,e,iidd_old and evol history correct */
/*                     before, these data was from before the second compact object formation and iidd_old completely off! */
/* change Jan28, 2008: mt1H[][],mt2H[][] size changed from 17-17 to 18-18 (to include output for Ka=Kb=17!!! */
/*                     ealier, output for MT1 and MT2 with Ka=17 or Kb=17 amy not have shown at all, or showd in error */
/* change Jan31, 2008: output for Livio&Riess SN Ia added: to sn1a.dat */
/* change Feb14, 2008: output for long GRBs added: longgrb.dat (fp400) */
/* change Feb16, 2008: in orb_change() eccentricity over 0.999+ cut down to 0.999: otherwise integration in odeint1() */
/*                     can take very very long (hrs for one system), this increases Rlobe size: but such systems very rare */
/* change Feb27, 2008: ecssna,b and idum_run passed to Xbin() and added to output for X-ray binaries: DetMT.dat and DetWIND.dat */
/* change Mar13, 2008: in get_Vkick6() (Hobbs kicks) zero kicks added explicitely: if Sigma3=0 do Vkick=0 */
/* change Apr26, 2008: in lost_env() change added: if(Mac<0.35) then for Ka==2 || Ka==3 star A becomes Ka=10 (He WD) */
/*                     as such a He core is too low mass to strat burning helium; otherwise what was happening is that afetr */
/*                     CE e.g., CE(11-3;11-7) K=3 would become K=7 with mass say 0.34 and then in the next step K=7 would */
/*                     become K=17 (Hybrid WD) that was incorrect; so some previous K=17 are really K=10 */
/* change Apr27, 2008: effective mass (M0a, M0b) added to file output "sn1a.dat" */
/* change Jul25, 2008: output for AM CVn added: file output "amcvn.dat" */
/* change Feb28, 2009: now binary.c can run single stars: see BINARY in sinbin.h */
/* change Mar11, 2009: in interWD() dMcr was set wrong (condition for low masses Mb<0.6 Msun was ignored): corrected */
/* change Mar18, 2009: wind2f() added in singl.c -- new mass loss rates, set by WIND=2 in sinbin.h */
/* change Mar19, 2009: in Xbin() timestep is either 1 Myr or equals to intrinsic timestep: see Xdtstep in sinbin.h */
/* change May05, 2009: changes in output, better output parameter management (nothing significant) */
/* change Jun17, 2009: error in comm_env3c1 (Nelemans CE)  -- Mbi got wrong value assigned: corrected */
/* change Aug12, 2009: change in singl.c in func. perturb(): for K=6,8,9 if Mc>MCh then the radius drops down to value */
/*                     for CO WD with M=1 msun (R=0.008058 Rsun = 5500 km). Before the radius would drop to a value of */
/*                     0.000014 Rsun=10 km (NS). Mdum=1.0 Msun instead of MCh is now used in Rc and Lc calculation */
/* change Dec01, 2009: error calculating dMwinda, dMwindb: they were calculated in [Msun] and used as rates [Msun/Myr] */
/*                     changed to rates: dMwinda=(Maold-Ma-dMmta*dt)/dt and dMwindb=(Mbold-Mb-dMmtb*dt)/dt */
/*                     affected: orbit calculation orb_change(), and output: symbiotic() and Xbin() functions*/
/* change Dec03, 2009: update: tidal constant E2: func. E2f() added: tides for massive stars VERY WEAK in any case */
/*                     old way use EE2=1 (sinbin.h): Hurley et al. 2002 formula; more effective radiative damping */
/*                     new way use EE2=2 (sinbin.h): Claret 2007 models used: less effective radiative damping */
/*                     new way more realistic, but models only for stars with convective cores and radiative env. */
/* change Jul28, 2010: Xcvf() improved -- better calculation of X-ray luminosities added. needs improvement */
/* change Jul30, 2010: lamf() added: calculation of lambda for CE: Chinese physical lambdas used: see chlambda in sinbin.h */
/*                     qeestionable use of Mzams in lamf(): if star gains mass in RLOF on MS than higher mass than Mzams */
/*                     should be employed, how high depends on the level of rejuvenation, but entire issue is messed up by */
/*                     the wind mass loss on MS. M0 can not be used as it gives current mass of a star on MS */
/*                     considering uncertainties of lambda calculation and ones inherent to Chinese paper: Mzams is OK */
/* change Oct27, 2010: optional new remnant mass calculation added: compactf2() from Chris F. (SASI SN mechanism) */
/*                     turned on by REMNANT=1 in sinbin.h */
/* change Oct29, 2010: syncronization/circularization bug correction: even if tidal forces were very weak if an expanding */
/*                     or shrinking star got rot. velocity = orbital velocity, it would get frozen to that velocity */
/*                     based on assumption that once the system becomes synchronized it stays synchronized: but it is wrong */
/*                     now i allow only such freezeing during ongoing RLOF, otherwise full numerical solution to tidal eqs. */
/*                     is used in orb_change(): it will slow down calculations and alter evolution of systems that were */
/*                     affected by such freezing */
/* change Dec01, 2010: optional new remnant mass calculation added: compactf3() from Chris F. (QUICK SN mechanism) */
/*                     turned on by REMNANT=2 in sinbin.h */
/* change Dec03, 2010: various options on kicks added: see parameters kickback in sinbin.h */
/* change May20, 2011: donor radius and roche lobe added in Xbin output: DetMT.dat and DetWIND.dat */
/* change May20, 2011: fb=0 added for AIC events in singl.c l.215, otherwise error was reported in get_Vkick6() */
/* change Jun20, 2011: in compactf2() and compactf3() in singl.c calculation of Menv corrected */
/* change Jun20, 2011: sinvel() exchanged -- it was giving wrong velocities of disrupted components */
/* change Jul15, 2011: singl(,SNR) argument added to provide proper SN counters recorded in info.dat */
/*                     uptonow the counters were increased too many times: SN callibrations were off */
/* change Jul26, 2011: output added for WD-WD binaries in file wdwd.dat: set WDWD to 1 in sinbin.h to activate */
/* change Jul27, 2011: change of Oct29, 2010 redone: full numerical approach was giving a lot of trouble: odeint() was not */
/*                     able to solve many cases circular and almost synchronized cases...., orb_change() corrected:  */
/*                     Sia, Sib introduced to control circ/sync and switch from full numerical to analytical solution */
/* change Jul28, 2011: in singl.c: CHeBf(), HEf() corrected: when K=3 was with almost no Menv, and it was passed as K=4 to HEf() */
/*                     the newly calculated core mass would exceed star mass and results in "nan" for star properties */
/* change Aug01, 2011: in singl.c: AGBf() corrected: when almost no Menv, sometimes Mcco may have exceeded the mass of star */
/*                     and would result in "nan" for star properties */
/* change Dec08, 2011: output added for paulina for RR Lyrae stars: file RRwd.dat */
/* change Mar13, 2012: fudge factor added to physical lambdas based on Ivanova & Chaichenets 2010, ApJ 731, L36 */
/*                     check golambda in sinbin.h; factor applied to all lambdas and all stars, set golambda=1 if not desired */
/* change Mar14, 2012: output added for Mergers of stars with NS/BH (MOSTLY WITHIN CE): transient project with Chris Fryer */
/*                     set Merger=1 in sinbin.h and output is written into Merger.dat */
/* change Sep17, 2012: in windf2() new winds added for WR stars from Zdziarski et al. 2012: see WRWINDS in sinbin.h */
/* change Dec06, 2012: CE=3 option added: mix of alpha + gamma CE prescription introduced (should be standard for WD simulations) */
/* change Dec14, 2012: correction of Sia, Sib counters: they were not reset to initial value between systems: wrong tidal calculations */
/* change Mar06, 2013: output altered in Xbin(): infiles DetMT.dat DetWIND.dat timestep added and precision of time increased */
/* change Apr24, 2014: output added for CE mergers (see TranCE in sinbin.h), transient project with Chris Kochanek */
/* change Apr24, 2014: famin added to allow for shift in minimum initial separation at ZAMS (check sinbin.h) */
/* change May28, 2014: new Sana et al. 2012 initial parameter distributions added: look at SS in sinbin.h */
/* change may30, 2014: amin in a0=get_a(amin,amax) corrected. for large e=0.999999 amin was larger than amax: */
/*                     now in such case  amin is set to 0.1*amax: and value of amax is given in sinbin.h */
/*                     not very important as for 10^6 systems it is only 10-20 cases like this; they become mergers anyway */
/* change Jan28, 2015: changes with damian: limiting mass of CO WD and ONeMg WD: singl.c in REMNANTf() */
/*                     fitting formulae for eta in binary.c in dMgainf() for He accretion onto CO and ONeMg WDs */
/*                     mass dependent size of He shell for double detonation SN1a from Iben and Tutukov 1991, Apj 370, 615 (eq.13) */
/*                     mass dependent size of He shell for double detonation SN1a from Shen and Bildsten 2009, ApJ 699, 1365 (fig.15) */
/* chanhe Mar12, 2015: in dMgainf we added accretion regimes with various eta: dMcrit6a -- dMcrit6d: from Ruiter et al. 2014 */
/* change Mar12, 2015: output for Transient.dat (either on or off) was HALTING evolution of EACH binary at first CE! */
/* change Mar18, 2015: warning in dMmtf() revised, Kb=0 and Kb=1 allowed: all physical -> leads to merger right away  */
/*                     warning: in dMmtf() non-remnant and not He MS star (K: %d) accretor to WD donor! */
/* change Mar18, 2015: wind mass accretion added (it is not on during ongoing RLOF) */
/* change Mar30, 2015: added functions for atmospheric RLOF from krystian: see dMamtf() */
/* change Mar30, 2015: added functions for wind RLOF from abate 2013 et al. Gbetaf() */
/* change May15, 2015: the warning listed in "change Mar18, 2015" removed entierly; some small stars may be accretors to WDs */
/* change May16, 2015: in HSGBf() time step was negative in some cases dtrg<0 for a star depleted by RLOF to bare CO core */
/*                     in such a case i send such star to REMNANTf() and make it CO WD, instead of having error messege */
/* change May17, 2015: proper resetting of dMallb and dMalla values to zero introduced */
/* change May19, 2015: change of max3() to min3() to select the most important non-RLOF mass transfer: all negative so min needed not max */
/* change Sep03, 2015: correction to a decision tree on wind accretion and atmospheric RLOF, also 2 bugs corrected on the same issues */
/* change Nov02, 2015: idum_run and iidd_old added in spin_evol() error messeges */
/* change Nov13, 2015: wind accretion fixed */
/* change Nov22, 2015: compact.dat only for binaries and with reduced output. for single stars output blocked with: */
/*                     if(0 && COMPACT==1 && ...)  but still in the code; can be unblocked for GAIA projects */
/* change Apr26, 2016: wind accretion done wrong during CE step: e.g., accretion of He instead of H, so acc. stopped alltogether */
/*                     in CE step (just before "Fbetaa=...": else if(inbin==1 && ce!=1) insted of else if(inbin==1)  */
/* change Apr27, 2016: addtion of pair instability pulsations (PIPs) and pair instability supernovas (PISs) */
/*                     to turn it on set: REMNANT=3 in sinbin.h, the change is done in singl.c in compactf4() */
/*                     K=18: new type added -- massless remnant from pair instability supernova */
/* change May01, 2016: file fp333 with evolutionary data replaced by to functions: evroute_add() and evroute_new() to */
/*                     manage the evolutionary history output (file fp333 output was causing problems on supercomputer ATLAS */
/* change Jun07, 2016: POP III initial distribution  added: set SS=2 in sinbin.h */
/* change May21, 2017: in all CE functions, e.g. comm_env3d1(), (*Om)=0.0; line removed -- position of line of nodes */
/*                     does not change in CE, (*om)=0.0; (*tau)=0.0; -- also removed: these are not defined for e=0 as in CE */
/*                     (*e)=(*e); (*i)=(*i); -- also removed as they do nothing */
/* change May22, 2017: Mhe,Mco are now not zeroed in REMNANTf() in singl.c; they keep their values while components are remnants */
/*                     this line was removed: (*Mhe)=(*Mco)=0.0; in REMNANTf() */
/*                     this change is used to initialize BH spins with bhspininit(); see sinbin.h and BHSPIN=2 (new: added) */
/* change May23, 2017: extra (6th) line added in compact.dat: with spins and tilts for BHs; ignore this line for NSs */
/* change May23, 2017: corrected expression for F in dMamtf(); from Krystian Ilkiewicz -- ffects only atmospheric RLOF calculations */
/* change May23, 2017: option added for NS/BH natal kicks without fallback decrease: set kick=7 and desired Sigma3 in sinbin.h */
/* change Jun03, 2017: introduction of new BH spin model: bhspininit1() */
/* change Aug06, 2017: introduction of CO/BH mass increase due to RLOF accretion spin up of progenitor on MS */
/*                     set Spin1 to 1 in sinbin.h: new code in compactf4() in singl.c, and in binary.c */
/* change Aug07, 2017: introduction of new BH spin model: bhspininit2() */
/* change Aug07, 2017: output added to compact.dat for WR formation (from star A and B) -- 10 numbers (1 line) */
/* change Aug08, 2017: output on spin orientations added for Davide Gerosa and Richard simulations */
/*                     one long extra line in compact.dat with 33 numbers */
/* change Aug08, 2017: new natal kick model: kicks decreasing with BH mass (2.5/Mbh): set kick=7 in sinbin.h */
/* change Aug08, 2017: new natal kick model: kicks decreasing with NS/BH mass (Mejecta/Mns or Mejecta/Mbh): set kick=8 in sinbin.h */
/* change Oct08, 2017: correction to the new BH spin model: bhspininit2(): new values of metallicity limits */
/* change Nov27, 2017: correction in sinvel(): center of mass velocity added from 1st SNa (it was not taken into account) */
/* change Dec07, 2017: correction of units in dMamtf():  1.0e-06/Msun -> 1.0e+06/Msun */
/* change Jan09, 2018: BH spin added into DetMT.dat and DetWIND.dat in Xbin() */
/* change Jan09, 2018: ecssna, ecssnb were in wrong order in second call to Xbin(): corrected */
/* change Jan09, 2018: Xout -> changed to Xout1 (output for RLOF XRBs) and Xout2 (output for WIND XRBs) */
/* change Feb28, 2018: ouput to DetMT.dat and DetWIND.dat corrected not to include systems during CE step */
/* change Mar30, 2018: in explode() added: snmerger=0 for exploding second star in disrupted binary */
/*                     otherwise, it was possible that disrupted binary was named SN merger at 2nd SNa */
/* change Aug06, 2018: new BH/NS accretion/X-ray luminosity introduced with Samaresh and JP */
/* change Oct30, 2018: in get_Vkick8() "if(x0<0.0) x0=0.0;" added to protect 3D kick be negative */
/*                     used with new Bray&Eldridge2018 alpha1=100km/s and beta1=-170km/s (sinbin.h) */
/* change Nov01, 2018: in Xbin1() condition added that output to Det*.dat files only if Xout1==1, Xout2==1 */
/* change Nov27, 2018: in singl.c compactf5() added: to do rapid with PPSN/PSN, set by REMNANT=4 in sinbin.h */
/* change Nov30, 2018: accumulation onto NS/BH corrected: corresponds to Samaresh new prescriptions now */
/* change Dec12, 2018: function Fbetaf() modified to include lower beta=0.125 in more continuos way */
/* change Jan04, 2019: in dMmtf(): added dMmtt to protect artificial thermal timescale step at nuclear RLOF end */
/* change Jan07, 2019: evolutionary times (MS, HG, AGB, ...) added in compact.dat output (fp200): one extra line */
/* change Jan16, 2019: introduction of new BH spin model: bhspininit3(), set by BHSPIN==4 in sinbin.h */
/* change Jan17, 2019: introduction of new PPSN/PSN model and smaller neutrino emission for BH mass estimate */
/*                     compactf6(): set by REMNANT==5 in sinbin.h (PPSN BHs more massive, 10%->1% neutrino mass loss */
/* change Feb16, 2019: introduction of new BH spin model: bhspininit10(), set by BHSPIN==10 in sinbin.h */
/*                     based on Jim Fuller calculations all BHs get a_natal=0.01 */
/* change Feb26, 2019: modification of new BH spin model: bhspininit3(), set by BHSPIN==4 in sinbin.h */
/* change Mar11, 2019: PPSN and PSN can be ON or OFF: in sinbin.h check PPSN and PSN parameters */
/* change Mar14, 2019: singl.c AGBf() function corrected to give good BH masses for WIND1=WIND2<1.0 */
/*                     for Mzams>30.0 && t>8.0, evol. stopped and BH formed; otherwise evol. continues: small BH mass */
/* change Apr21, 2019: introduction of new PPSN/PSN model and smaller neutrino emission for BH mass estimate */
/*                     compactf8(): set by REMNANT==7 in sinbin.h (PPSN BHs more massive, 10%->1% neutrino mass loss */
/* change May28, 2019: aspinb added in dJf() to properly calculate dervatives for decide() on CE */
/* change May28, 2019: in dMamtf() used for symbiotic binaries (Krystian Ilkiewicz) some logical (if) statements corrected */
/* change Jun13, 2019: compactfN() functions reorganized and added: N=3,4,5,6,7,8,9: set by REMNANT in sinbin.h */
/* change Jun13, 2019: extra mass loss during massive BH formation allowed now: set by GRB parameter in sinbin.h */
/* change Oct29, 2019: in Xbin1() wind accretion corrected in case of WindRLOF==1 */
/* change Oct30, 2019: in Xbin1() calculates either RLOF ***OR*** WIND luminosity (not the two togehter) */
/*                     same is done with mass accretion: it is either RLOF or WIND taken into acount (RLOF takes precedence) */
/* change Oct30, 2019: Leddf1() was giving factor of 2 smaller value of Ledd: corrected (mass accretion was done correctly) */
/* change Oct30, 2019: accumulation from wind (dMwindacca, dMwindaccb, dMwindrlofa, dMwindrlofb)) set now only by dMgainf() */
/* change Oct30, 2019: in dMgainf() we limit accretion onto normal stars (K<10) to Eddington limit (we still use Fa) */
/* change Nov12, 2019: mt1H[] and mt2H[] -- changed to mt1H and mt2H: better management of evolutionary output */
/* change Nov12, 2019: in several places potential division by zero avoided by added if statements */
/* change Dec09, 2019: in Xbin1() Ledd1 was calculated for the wrong type of donor: Ledd1(Mb,Kb) -> now: Ledd1(Mb,Ka) */
/*                     Leddf1() changed to protect against having NS or BH being donors in this function */
/* change Mar03, 2020: qmin added for initial mass ratio Sana distribution -- see SanaOriginal parameter in sinbin.h */
/* change Sep17, 2020: two new REMNANT mass fucntions added: no PPSN + high PSN, BHs upto 90Msun; set by REMNANT=9,10 */
/* change Feb09, 2021: upper limit on stellar radius added for PopI/II stars; see RMAX in sinbin.h. change in perturb() */
/* change Feb10, 2021: correction in Gbetaf(): x-had wrong definition; affected only wind RLOF calculations */
/* change Mar06, 2021: compactf12() added: extended partial fallback from Mco=11Msun to McoExt set in sinbin.h */
/*                     corresponds to modified compactf8(); set with REMNANT=12 -- delayed SN engine weak PPSN */
/* change Mar08, 2021: fit for Rmax in perturb() revised: it is now for 10% Zsun */
/* change Mar31, 2021: new (shallower) Z-dependence of wind mass loss for hot H-rich stars added from Vink&Sander2021 */
/*                     set by Z33 in windf2(), the dependance for cooler stars can be changed now with Z22 (sinbin.h) */
/* change Apr06, 2021: extended fallback added in function compactf10(), see McoExt and REMNANT=8 in sinbin.h */
/* change Apr17, 2021: correction to introduction of maximum radius (RMAX==1) in perturb() */
/* change Oct05, 2022: combined mass (fallback-decreased) and neutrino NS/BH kicks added: new get_Vkick9() */
/*                     controled by Sigma4 and Sigma5 set in sinbin.h with kick=10 */
/* change Oct12, 2022: added new (Pavlovskii2017/Olka2021) new CE/stableRLOF criteria: CE=4 in sinbin.h */
/* change Oct13, 2022: new remnant mass functions added: compactf13() and compact14(f) set by REMNANT=13/14 in sinbin.h */
/*                     includes smooth transition from rapid-to-delayed engine (fmix and McritNSBH) */
/*                     includes various options on PPSN and PSN (XXPSN, PPSN, PSN, Mppsn1, Mppsn2, Mpsn1, Mpsn2) */
/* change Nov21, 2022: new (H-rich and He-rich stars) winds added: windf3() in singl.c, set by WIND=3 in sinbin.h */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>

#define  External
#include "sinbin.h"


#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
#define JMAX 40
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS2 1.2e-7
#define RNMX (1.0-EPS2)
#define EPS3 1.0e-5
#define JMAX3 2000                      /* my change to integrate Tmerge of high e systems */
#define FUNC(x) ((*func)(x))
#define FUNC1(x,M) ((*func)(x,M))

#define CON 1.4
#define CON2 (CON*CON)
#define BIG1 1.0e30
#define NTAB1 10
#define SAFE1 2.0

#define MAXSTP 10000
#define MAXSTP1 1000000
#define TINY 1.0e-30
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define NR_END 1
#define FREE_ARG char*
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))
#define NRANSI



/*----------------------------------- binary evolution ----------------------------------------------*/

int interaction1(double Mzamsa, double *M0a, double *Ma, int *Ka, int *Kpa, double *tvira,
                 double *La, double *Ra, double *Mca, double *Mhea, double *Mcoa, int *flaga, double *dta, double *tstarta,
                 double Mzamsb, double *M0b, double *Mb, int *Kb, int *Kpb, double *tvirb, double *Lb, double *Rb,
                 double *Mcb, double *Mheb, double *Mcob, int *flagb, double *dtb, double *tstartb,
                 double dMmta, double *dMmtb, double *a, double *e, double *i, double *Om, double *om,
                 double *tau, double t, double dt, double *ta_end, double *tb_end, int *stop1,
                 int *stop2, int *cmt1, double dtmt, double *dMce);

int interaction2(double Mzamsa, double *M0a, double *Ma, int *Ka, int *Kpa, double *tvira, double *La, double *Ra,
                 double *Mca, double *Mhea, double *Mcoa, int *flaga, double *dta, double *tstarta,
                 double Mzamsb, double *M0b, double *Mb, int *Kb, int *Kpb, double *tvirb, double *Lb, double *Rb,
                 double *Mcb, double *Mheb, double *Mcob, int *flagb, double *dtb, double *tstartb, double *a,
                 double *e, double *i, double *Om, double *om, double *tau, double t, double *ta_end, double *tb_end,
                 int *stop1, int *stop2);

void get_SNcounters();

void undo_SNcounters();

void lost_env(double M0a, double Ma, double Mca, double *Mhea, double Mcoa, double *tvira,
              double *tstarta, double *dta, int *Ka, int *Kpa);

void comm_env3b0(double Mzamsa, double *Ma, double Mca, double Ra, int Ka, double *Mb,  double *a, double *e,
                 double *i, double *Om, double *om, double *tau);

void comm_env3b1(double Mzamsa, double *Ma, double Mca, double Ra, int Ka, double *Mb, int *Kb, double *a, double *e,
                 double *i, double *Om, double *om, double *tau, double *dMce);

void comm_env3c0(double Mzamsa, double *Ma, double Mca, double Ra, int Ka, double Mzamsb, double *Mb, double Mcb,
                 double Rb, int Kb, double *a, double *e, double *i, double *Om, double *om, double *tau);

void comm_env3c1(double *Ma, double Mca, double Ra, int Ka, double *Mb, double Mcb, double *a, double *e, double *i,
                 double *Om, double *om, double *tau);

void comm_env3d0(double Mzamsa, double *Ma, double Mca, double Ra, double *Mb, double *a, double *e,
                 double *i, double *Om, double *om, double *tau, double t, int Ka, int *Kb, double *dMce);

void comm_env3d1(double Mzamsa, double *Ma, double Mca, double Ra, double *Mb, double *a, double *e,
                 double *i, double *Om, double *om, double *tau, double t, int Ka, int *Kb, double *dMce);

double bondi(double Mzamsa, double Ma, double Mca, double Ra, double Mb, double a, int Ka, int Kb);

double lamf(double Mzamsa, double Ma, double Mca, double Ra, int Ka);

double roche(double q, double a);

double Aroche(double q, double Rl);

int explode(double Ma, double Maf, int Ka, double fraca, double Mb, double Rb, double *Vsm,
            double *Vsa, double *Vsb, double *Vextr, double *a, double *e, double *i, double *Om, double *om,
            double *tau, double *Vkick, int *cmt3, double *texp, int inbin, int ecssna, int fba,
            double *jxi, double *jyi, double *jzi, double *jxf, double *jyf, double *jzf);

int sinvel(double Ma, double Mapre, double Mb, double Rb, double *Xa, double *Va, double *wkick, double *Vsm, double *Vsa, double *Vsb);

void rotate2for(double *Vi, double *Vf, double c1, double s1, double c2, double s2);

void rotate2back(double *Vi, double *Vf, double c1, double s1, double c2, double s2);

double tmerge(double Ma, double Mb, double a, double e);

void symbiotic(double Ma, double dMwinda, double Mca, double La, double Ra, double Rb,
               double t, double Mb, double a, double e, int Ka, int Kb, double dMmta, double dMmtb,
               double dMtran, int mttype, double Mzamsa, double Mzamsb, double a0, double e0);

void Xbin(double *Lxmt, double dMtran, double Ma, double dMwinda, double Mca, double dMmta, double La, double Ra,
          double t, double dt, double Mb, double Rb, double a, double e, int Ka, int Kb, double *Vsm, double tsn,
          int *mark10, int *mark11, int mttype, int flagbra, int flagbrb, int mtflag1a, int mtflag1b, int mtflag1ab,
          double Maio, double Mbio, double dMcea, double dMceb, double Mzamsa, double Mzamsb, double a0, double e0,
          int ecssna, int ecssnb, double aspinb);

void Xbin1(double t, double dt, double Ma, double Mb, int Ka, int Kb, double Ra, double dMmta, double dMwinda, double aspinb, double a, double e, double La);

double Fbetaf1(double dMwinda, double Ma, double Ra, double Rb, double Mb, double a, double e, int Ka, int Kb, double aspinb);

double Fbetaf(double Ma, double Ra, double Mb, double a, double e, int Ka);

double Gbetaf(double Ma,double Mb,double Ta,double Ra,double Rla);

void Xwdwd(double *Lxwd, double dMmta, double Mb, double Rb, int Ka, int Kb);

double Xcvf(int Ka, int Kb, double Mb, double Rb, double dMmta, double *dMacc, int *cvtype);

void spin_evol(int Ka, double *Ma, double *aspin, double *Mdisa, double Mrest);

double bhspininit1(double Mco);

double bhspininit2(double Mco);

double bhspininit3(double Mco);

double bhspininit10(void);

/*--------------------------------- detailed MT computation ----------------------------------------*/


double dMmtf(double Ma, double Mb, double Ra, double Rb, double La, double Lb, double wa, double Ia, double KTa,
             double wb, double Ib, double KTb, int Ka, int Kb, double a, double e, double t, double dt,
             double *input1, int *input2, double dMmta_old, double dMmtb_old, int *stop2, int *mttype, int mttypelast,
             int *doce, int *merger, double *dMth, double *Tth, int *dec, double aspinb, double Mzamsa, double Mzamsb);

double dMgainf(double dMmta, double dMtran, double Mb, int Ka, int Kb, double Rb, double *Mout1, double *Mout2,
               int *mark_out, int *doce, int *merger, int *mark777, double *Mb0, double aspinb);

double Mshellf(double Mb, int Kb, double dMmta);

double interWD(double Mb, double dMmta);

double dMeddf(int Ka, double Ra, int Kb);

double dMeddf1(double Ma, int Kb);

double Leddf1(double Ma, int Kb);

double dMcrif(int Ka, double Ma);

double dMtranf(double Ma, double Mb, double a, double e, int Ka, int Kb);

int decide(double Ma, double Mb, double a, double e, int Ka, int Kb, double Ra, double Rb, double Zstar, double tth, double aspinb);

double dJf(double Mai, double Mbi, double ai, double Maf, double tth, double Rb, int Ka, int Kb, double *x, double aspinb);

double tthf(double Ma, double Ra, double La, int Ka);

double tmbf(double Ma, double Ra, double wa, double Mb, double a, double La, int Ka);

double Mminconvf(void);

double Mmaxconvf(void);

double dlnRl(int Ka, int Kb, double Ma, double Mb, double a, double dMmta_old, double dMmtb_old);

double dlnRdlnM(double M, double *input1, int * input2, int *stop2);

double funclnRM(double lnM, double *input1, int *input2);

double dlnRdt(double t, double dt, double *input1, int * input2, int *stop2);

double funclnRt(double t, double *input1, int *input2);

double dfridr1(double (*func)(double, double *, int *), double x, double h, double *err, double *input1, int *input2, int *stop2);


/*----------------------------------- orbit computation --------------------------------------------*/

void orbit1(double t0, double m1, double m2, double a, double e, double i, double Om, double om,
            double tau, double *X, double *V);

int orbit2(double texp, double m1, double m2, double *X, double *V, double *a, double *e, double *i,
           double *Om, double *om, double *tau);

double kepler(double M, double e);

double get_ang(double t1, double t2);

double dist(double a1, double a2, double a3, double b1, double b2, double b3);


/*---------------------------------------- tools ---------------------------------------------------*/

void copyV(double *V1, double *V2);

void copySN(double Maold, int Kaold, double Mbold, int Kbold, double aold, double eold, double Ma,
            double *MpgaA, int *KpgaA, double *MpgbA, int *KpgbA, double *apgA, double *epgA, double *MendaA);

void copySin(double Mzamsa, double M0a, double Ma, int Ka, double tbeg, double tvira, double tenda, double La,
             double Ra, double Mca, double Mhea, double Mcoa, int flaga, double dta, double Mprea, int Kpa,
             double tstarta, double fraca, double dMmta, int ecssna, double *sa1, double *sa2, double *sa3, int *sa4,
             double *sa5, double *sa6, double *sa7, double *sa8, double *sa9, double *sa10, double *sa11, double *sa12,
             int *sa13, double *sa14, double *sa15, int *sa16, double *sa17, double *sa18, double *sa19, int *sa20);

void rec_dat1a(FILE *fp1, double Ma, double Mb, double *V1, double *V2, double a, double e,
               double ta_end, double tb_end, double Tmr, int Ka, int Kb,
               double MpgaA, int KpgaA, double MpgbA, int KpgbA, double apgA, double epgA, double tpgA, double MendaA,
               double MpgaB, int KpgaB, double MpgbB, int KpgbB, double apgB, double epgB, double tpgB, double MendaB,
               double Mzamsa, double Mzamsb, double a0, double e0);

void rec_dat1b(FILE *fp1, double Ma, double Mb, double *V1, double *V2, double *Vsa, double *Vsb,
               double *Vextr, double a, double e, double ta_end, double tb_end, double Tmr, int Ka, int Kb,
               double MpgaA, int KpgaA, double MpgbA, int KpgbA, double apgA, double epgA, double tpgA, double MendaA,
               double MpgaB, int KpgaB, double MpgbB, int KpgbB, double apgB, double epgB, double tpgB, double MendaB,
               double Mzamsa, double Mzamsb, double a0, double e0, double dMcea, double dMceb, int inbinA, int inbinB);

void rec_dat2(FILE *fp2, int j, int aicns, int aicbh, int badorb);

double min3(double x1, double x2, double x3);


/*------------------------------ initail binary parameters -----------------------------------------*/

double get_q1(double Min, double Max);
double get_q2(double Min, double Max);
double get_q3(double qmin);
double get_SS(double Min, double Max, double Ex);
double get_M(double Min, double Max);
double get_M_linear(double Min, double Max);
double get_a(double Min, double Max);
double get_e(double Min, double Max);
double get_i(double Min, double Max);
double get_Om(double Min, double Max);
double get_om(double Min, double Max);
double get_texp(double Min, double Max);
void get_Vkick1(double *Vkick);
void get_Vkick2(double *Vkick, double fraca, int fba);
void get_Vkick3(double *Vkick, int Ka);
void get_Vkick4(double *Vkick);
void get_Vkick5(double *Vkick, int Ka, double fraca, int *cmt3, int fba);
void get_Vkick6(double *Vkick, double fraca, int fba);
void get_Vkick7(double *Vkick, double Mf);
void get_Vkick8(double *Vkick, double Msn, double Mf);
void get_Vkick9(double *Vkick, double fraca, int fba);
double func_max1(double v);
double func_max1a(double v);
double func_max1b(double v);
double func_max2(double mtop, double P0);
double func_max2a(double mtop, double P0);
double func_max2b(double mtop, double P0);
double get_flat(double Min, double Max);
double interp(double p, double xPacz[], double yPacz[], int NPacz);


double rtbis1(double (*func)(double, double), double x1, double x2, double xacc, double P0);
double func1(double x, double P0);
int zbrac2(double (*func)(double, double, double, double, double, int, int, double), double *x1, double *x2,
          double Ma, double Mb, double Ia, double Ib, int Ka, int Kb, double Ac);
double rtbis2(double (*func)(double, double, double, double, double, int, int, double), double x1, double x2,
              double xacc, double Ma, double Mb, double Ia, double Ib, int Ka, int Kb, double Ac);
double func2(double x, double Ma, double Mb, double Ia, double Ib, int Ka, int Kb, double Ac);
double erffNR(double x);
double gammp(double a, double x);
void gcf(double *gammcf, double a, double x, double *gln);
void gser(double *gamser, double a, double x, double *gln);
double gammln(double xx);
void nrerror(char error_text[]);
float ran2(long *idum);
float ran3(long *idum);
double func_tmr(double e);
double qtrap(double (*func)(double), double a, double b);
double trapzd(double (*func)(double), double a, double b, int n);


double get_gauss(double xave, double sig, double min1, double max1);
double func_x(double x, double xave, double sig, double P0);
double int_gauss(double x, double xave, double sig);
double rtbis(double (*func)(double, double, double, double),
             double x1, double x2, double xacc, double xave, double sig, double P0);



/*---------------------------------- hyper-critical accretion --------------------------------------*/

void derivs(double, double [], double [], double, double);
void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
     double hmin, int *nok, int *nbad, void (*derivs)(double, double [], double [], double, double),
     void (*rkqs)(double [], double [], int, double *, double, double, double [], double *, double *,
     void (*)(double, double [], double [], double, double), double, double), double Mca, double lambda);
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
     double yscal[], double *hdid, double *hnext,
     void (*derivs)(double, double [], double [], double, double), double Mca, double lambda);
void rkck(double y[], double dydx[], int n, double x, double h, double yout[], double yerr[],
     void (*derivs)(double, double [], double [], double, double), double Mca, double lambda);
double *dvector(int nl, int nh);
void free_dvector(double *v, int nl, int nh);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);


/*------------------------------- MB + GR + Tidial breaking --------------------------------------*/

double get_w1(double M, double R);

double get_w2(double M, double R);

int orb_change(double t1, double t2, double *a, double *e, double *wa, double *wb, double tvira, double tvirb,
               double Ma, double Mb, double M0a, double M0b, double Mzamsa, double Mzamsb, double dMwinda,
               double dMwindb, double Mca, double Mcb, double Ra, double Rb, double Raold, double Rbold,
               double La, double Lb, int Ka, int Kb, int Kaold, int Kbold, int mt, int ce, int mttype,
               double *Iaold, double *Ibold, double *KTa, double *KTb, double dMmta, double dMmtb,
               double Mdisa, double Mdisb, int *darwin);

void dorbdt(double t, double y[], double dydx[], double KTa, double KTb,
            double Ma, double Mb, double dMa, double dMb, double Ia, double Ib,
            double Ra, double Rb, double Rca, double Rcb, double wcrit_a, double wcrit_b,
            double La, double Lb, int Ka, int Kb, int magb_a, int magb_b);

double E2f(int Ka, double Ma, double Mzamsa);

double KTf(double tvira, double Ma, double M0a, double Mzamsa, double Mca, double Ra,
           double Rca, double La, double wa, int Ka, double Mb, double a, double e);

double tau_tidef(double Ma, double Mb, double a, double e, double wa, double Ia, double KTa, double Ra,
                 double La, int Ka);

double Inerf(double Ma, double Mca, double Ra, double Rca, int Ka);

double Renv_con(double tvira, double Ma, double M0a, double Ra, double Rca, double La, int Ka);

double Menv_con(double tvira, double Ma, double M0a, double Mca, double Ra, double La, int Ka);

double Rcf(double tvira, double Ma, double M0a, double Mca, int Ka);

double wcritf(double Ma, int Ka);


int odeint1(double ystart[], int nvar, double x1, double x2, double eps, double h1, double hmin,
             int *nok, int *nbad,
 	     void (*derivs1)(double, double [], double [], double, double, double, double, double,
 	                     double, double, double, double, double, double, double, double, double,
 	                     double, double, int, int, int, int),
 	     int (*rkqs1)(double [], double [], int, double *, double, double, double [], double *, double *,
 	                  void (*)(double, double [], double [], double, double, double, double, double,
 	                           double, double, double, double, double, double, double, double, double,
 	                           double, double, int, int, int, int),
 	                  double, double, double, double, double, double, double, double, double, double,
                          double, double, double, double, double, double, int, int, int, int),
             double KTa, double KTb, double Ma, double Mb, double dMa, double dMb, double Ia, double Ib,
             double Ra, double Rb, double Rca, double Rcb, double wcrit_a, double wcrit_b,
             double La, double Lb, int Ka, int Kb, int magb_a, int magb_b);

int rkqs1(double y[], double dydx[], int n, double *x, double htry, double eps,
          double yscal[], double *hdid, double *hnext,
          void (*derivs1)(double, double [], double [], double, double, double, double, double,
	                  double, double, double, double, double, double, double, double, double,
 	                  double, double, int, int, int, int),
          double KTa, double KTb, double Ma, double Mb, double dMa, double dMb, double Ia, double Ib,
          double Ra, double Rb, double Rca, double Rcb, double wcrit_a, double wcrit_b,
          double La, double Lb, int Ka, int Kb, int magb_a, int magb_b);

void rkck1(double y[], double dydx[], int n, double x, double h, double yout[], double yerr[],
           void (*derivs1)(double, double [], double [], double, double, double, double, double,
 	                   double, double, double, double, double, double, double, double, double,
 	                   double, double, int, int, int, int),
           double KTa, double KTb, double Ma, double Mb, double dMa, double dMb, double Ia, double Ib,
           double Ra, double Rb, double Rca, double Rcb, double wcrit_a, double wcrit_b,
           double La, double Lb, int Ka, int Kb, int magb_a, int magb_b);

/*--------------------------------- atmospheric RLOF: from krystian -------------------------------*/

double dMamtf(double Ma, double Mb, double La, double Ra, double Rla, double a, int Ka);
double ro_ph(double M, double R, double Teff);
double ro_ph_z00000013415(double logg, double logT);
double ro_ph_z000001341446(double logg, double logT);
double ro_ph_z0000134117(double logg, double logT);
double ro_ph_z0000423904(double logg, double logT);
double ro_ph_z0001338401(double logg, double logT);
double ro_ph_z000421151(double logg, double logT);
double ro_ph_z0013113(double logg, double logT);
double ro_ph_z0023094(double logg, double logT);
double ro_ph_z004050(double logg, double logT);
double ro_ph_z0070598(double logg, double logT);
double ro_ph_z0122(double logg, double logT);
double ro_ph_z020267(double logg, double logT);
double ro_ph_z0322626(double logg, double logT);
double ro_ph_z048359(double logg, double logT);
double ro_ph_z067218(double logg, double logT);

/*---------------------------------- evolutionary output: from grzegorz ---------------------------*/
void evroute_add(char *format, ...);
void evroute_new();

/*------------------------------------ global variables -------------------------------------------*/

int iidd=0;                                        /* iidd - counts how many times ran3() was used */
int idum_run;                                                  /* stores value of idum for the run */
int Sia=-1;                             /* star A: counts consecutive steps with syncA+circ binary */
int Sib=-1;                             /* star B: counts consecutive steps with syncB+circ binary */


/*-------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
 FILE *fp1,*fp2,*fp3,*fp141,*fp142,*fp143,*fp144,*fp145,*fp210,*fp220;
 double Ra0,Rb0,Ra,Rb,Rla,Rlb,RFa,RFb,Raold,Rbold;                                                   /* [Rsun] */
 double La,Lb;                                                                                       /* [Lsun] */
 double Ta,Tb;                                                          /* effective temperatures of stars [K] */
 double Mca,Mcb,Mhea,Mcoa,Mheb,Mcob;                                                                 /* [Msun] */
 double q0,q,qmin,numberq;                                /* q0=Mb0/Ma0 - in the whole program (initial value) */
 double Mzamsa,Mzamsb,M0a,M0b,Ma,Mb,Maold,Mbold;                                                     /* [Msun] */
 double Mprea,Mpreb;                                                               /* presupernova mass [Msun] */
 double amin,dper;                                                                                   /* [Rsun] */
 double a0,e0,i0,Om0,om0,tau0,a,e,i,Om,om,tau;
 double wa,wb,waold,wbold;                                         /* component angular spin velocity [Myr^-1] */
 double wbreaka,wbreakb;                                  /* component break up angular spin velocity [Myr^-1] */
 double worb;                                                              /* angular binary velocity [Myr^-1] */
 double Iaold,Ibold,Iatmp,Ibtmp;                                             /* moment of inetria of componets */
 double KTa,KTb;                                                                           /* tidal constatnts */
 double t,tbeg,tend,tenda,tendb,dt,dt1,tvira,tvirb,dta,dtb,tstarta,tstartb,ta_end,tb_end,dtorb; /* times [Myr] */
 double tdis;                                                     /* time of binary disruption in SN explosion */
 double tform;                                                         /* formation time of WD-WD binary [Myr] */
 double Tmr;                                                                              /* merger time [Myr] */
 double Tms,texp;
 double aold,eold,MpgaA,MpgbA,MendaA,MpgaB,MpgbB,MendbB,apgA,apgB,epgA,epgB,tpgA,tpgB; /* values before SN A,B */
 double MccoaA,McheaA,fracaA,MccobB,MchebB,fracbB;                                      /* if no SN then zeros */
 int fba,fbb;                                /* start: -1, 0: no fall back, 1: partial fall back, 2: direct BH */

 double lmin,lmax,lxave;                               /* variables for initial distribution of a0 for POP III */

 int inbinA,inbinB;                                  /* 0- no binary after, 1-binary survived A/B SN explosion */

 double Maio,Mbio;                                                      /* mass of NS at the time of formation */
 double dMcea,dMceb,dMceaold,dMcebold;                              /* mass accreted by NS/BH in all CE events */
 int aicmark,aicns=0,aicbh=0;                                                              /* AIC housekeeping */
 int aicce1,aicce2;

 double dMmsa,dMmsb;                                  /* mass accreted in RLOF onto star A/B during MS in RLOF */
 double Mina,Minb;                                  /* mass of star A/B (accretor) at the beginning of MS RLOF */
 double Fmsa,Fmsb;                                /* fraction of accreted mass relative to entire MS star mass */
 int mark76a,mark76b,mark77a,mark77b;

 int spinout;                                                                       /* flag for BH spin output */
 double Mrest;                                                 /* rest mass of matter that is accreted onto BH */
 double Mdisa,Mdisb;                             /* mass that disappeared in accretion onto BH: rest M- grav M */
 double aspina,aspinb;                  /* spin of BH from star A, spin of BH from star B: -1.0 if A/B not BHs */
 double iA,iB;                                        /* tilt [rad] after SNA and SNB -- initial tilt is in i0 */
 double aspina0,aspina1,aspina2,aspinb0,aspinb1,aspinb2;              /* 0-formation, 1-after CE, 2-after RLOF */
 double Mbha0,Mbha1,Mbha2,Mbhb0,Mbhb1,Mbhb2;                                       /* a -- star A, b -- star B */
 double tbhenda0,tbhenda1,tbhenda2,tbhendb0,tbhendb1,tbhendb2;

 double RmaxA,RmaxB;

 int comdone;                                            /* controls output for compact binaries "compact.dat" */


 int deriv;                           /* sets what type of call: 0-regular, 1-derivative estimation to singl() */

 double dtbr;                                        /* maximum time step after system detaches in detailed MT */
 int maxbr;                                                                /* kept only for next "maxbr" steps */

 int mtflag1ab,mtflag1a,mtflag1b,mtflag2,mtflag3;        /* show uncircularized and unsynchronized cases at MT */
 int darwin;                           /* system not solved by ODE at Darwin Ins. phase --- send to MT by hand */

 int bhdone;                                                                /* controls output for BH binaries */

 int flagbra,flagbrb;                                                 /* flags star spin over breakup velocity */

 int Is;                                      /* after SN: 1: orbit bound, 0: disrupted, -1: components merged */
 double Vsm[3],Vsm0[3],VsmA[3],VsmB[3];       /* store center of mass velocity at different moments [Rsun/day] */
 double Vsa[3],Vsb[3],Vextr[3],Vkick[3];

 int mark1;                                            /* controls quality of orbital solution in orb_change() */
 int badorb=0;                      /* counts number of bad orb sol. and stopped systems: reported in info.dat */

 double input1[16];
 int input2[4];
 double dtmt;                                                   /* max. time step during MT for X-ray binaries */
 double Lxmt;                                                            /* MT X-ray binary luminosity [erg/s] */
 double Lxwd;                                                             /* MT WD-WD X-ray luminosity [erg/s] */
 double Lwd_a,Lwd_b;                                                   /* bolometric luminosity for WDs [Lsun] */
 double dMmta,dMmtb,dMmta_old,dMmtb_old,dMbh;                     /* RLOF rate (transfer/accretion) [Msun/Myr] */
 double dMalla,dMallb;                                                /* Total mass gain/loss: RLOF: all types */
                                   /* include (only) accretion from wind, as wind mass loss treated in singl() */
 double Fbetaa,Fbetab;                                                    /* fraction of wind captured by wind */
 double Gbetaa,Gbetab;
 double dMatmorlofa,dMwindrlofa,dMwindacca;
 double dMatmorlofb,dMwindrlofb,dMwindaccb;

 double dMtran,dMedd;                                                                   /* MT rates [Msun/Myr] */
 double dMwinda,dMwindb;                                                     /* Wind mass loss rate [Msun/Myr] */
 double dMtha,dMthb,Ttha,Tthb,Ratrue,Rbtrue;                                           /* thermal MT variables */
 double Mout1a,Mout1b;                                          /* [Msun] SN Ia outburst mass for accreting WD */
 double Mout2a,Mout2b;                                                /* [Msun] He WD->He ZAMS transition mass */
 int mark_outa,mark_outb;                                                   /* markers to control Mouta, Moutb */
 double Macc0a,Macc0b;
 int mark777a,mark777b;
 int deca,decb;

 double twra,awra,ewra,Mwra,Mcoma;                       /* system parameters when A becomes WR star (K=7,8,9) */
 double twrb,awrb,ewrb,Mwrb,Mcomb;                       /* system parameters when B becomes WR star (K=7,8,9) */
 int mark88a,mark88b;

 double tlast,Mdonlast,Macclast,alast;                          /* store data on last MT event for compact.dat */
 int Kdonlast,Kacclast,mttypelast;              /* mttypelast: 0-noMT, 5-deg. donor, 4-nuc/GR/MB, 6-SCE, 7-DCE */

 int Kawd,Kbwd,Kawdold,Kbwdold;                                 /* tells if A or B is WD for CV classification */
 int cvdon;                                                      /* mark the donor in CV, 1-A donor, 2-B donor */
 int cvtype;                                             /* 0-not a CV, 1-intermediate polar, 2-CV (not an IP) */
 double dMacc,Lxcv;                                               /* accretion rate and X-ray luminosity of CV */

 double Mwda_0,Mwdb_0;                                                        /* formation mass of WDa and WDb */

 double Mheaold,Mcoaold,Mhebold,Mcobold;

 int mark10,mark11;                                                          /* output controls for X-ray bin. */
 int filledX[100];                                              /* controls detailed output for X-ray binaries */
 int mttype,mttype_old;                                      /* tells what MT type took place for X-ray binary */

 int j,jj,ii,ll,kk,max;
 char his[100];



 double atmp1,atmp2;

 int ecssna,ecssnb;                                               /* 1- ECS SN, 0-(either no SN or regular SN) */

 double fraca,fracb;                           /* amount of fall back of A,B in their SNs: (-1)-not applicable */
                                                 /* [0.0:1.0]: 0.0-no fb,1.0 complete fb, inbetween partial fb */
 int inbin;                           /* mark if stars A,B are still in the binary (1) or are single stars (0) */
 int sna,snb;                               /* mark if for a given star SN already took (=1) place or not (=0) */
 int sn;                                                                             /* SN step: 0-no SN, 1-SN */
 int wda,wdb;                                       /* mark if a given star become already WD (=1) or not (=0) */
 int WDout;                                                               /* controls output for WDWD binaries */
 int SPCout;                                                               /* controls output bhnswd.dat file */
 int OBstar;                                                                 /* controls output OBbin.dat file */
 int flaga,flagb;                                                              /* flags needed to run single() */
 int stop1,stop2;                         /* stop1,stop2=1 evol of binary stoped, stop1,stop2=0 evol continues */
 int cmt1;                             /* control counting not proper treatment of rejuvenation during cons.MT */
 int cmt3;                                                     /* control filling Paczynski distribuant arrays */
 int Ka,Kb,Kpa,Kpb,Kaold,Kbold,KpgaA,KpgbA,KpgaB,KpgbB;                                        /* type of star */
 int PP,pp;

 double tmsa,tmsb;                                                        /* MS lifetime of a given star [Myr] */
 double fmsa,fmsb;                                                    /* fractional age on MS for a given star */
 int typcls;                                                         /* type of binary for claster simulations */


 int mt;                                                                    /* detailed MT step: 0-no MT, 1-MT */
 int ce;                                                                             /* CE step: 0-no CE, 1-CE */
 int doce;                                        /* decides if CE should take place: 0-shouldn't CE, 1-should */
 int merger;                                  /* tells if merger encountered in dMmtf(), 1-merger, 0-no merger */

 int rem;                                                          /* controls the output for to data.dat file */

 int typwd;                                                                                   /* type of SN Ia */

 double sa1,sa2,sa3,sa5,sa6,sa7,sa8,sa9,sa10,sa11,sa12,sa14,sa15,sa17,sa18,sa19;    /* storage of single  data */
 int sa4,sa13,sa16,sa20;
 double sb1,sb2,sb3,sb5,sb6,sb7,sb8,sb9,sb10,sb11,sb12,sb14,sb15,sb17,sb18,sb19;
 int sb4,sb13,sb16,sb20;

 double per0,per;                                                                             /* binary period */
 double J0,J1;                                                                /* binary total angular momentum */
 int maxred;                            /* number of timestep reductions in order to keep change of Jtot small */

 int mt1H,mt2H;                                                  /* controls of evolutionary history of system */
 int qtu1;

 int Kace,Kbce,markce,typce,cea,cedonor;                                    /* pre CE parameters of the system */
 double ace,ece,Mace,Mbce,Race,Rbce,Mcace,Mcbce,Mheace,Mhebce,Mcoace,Mcobce;

 int amcvn,Ka11,Kb11,Ka22,Kb22;                                                           /* AM CVn parameters */
 double t11,Ma11,Mb11,a11,e11,Ra11,Rb11,t22,Ma22,Mb22,a22,e22,Ra22,Rb22;

 double jx_0,jy_0,jz_0;                                     /* L vector initially: binary ang. momentum vector */
 double jx_1,jy_1,jz_1;                                                                  /* L vector after SNA */
 double jx_2,jy_2,jz_2;                                                                  /* L vector after SNB */
 double a_0,e_0,i_0,Om_0,om_0,tau_0;                                   /* orbital parameters just prior to SNA */
 double a_0a,e_0a,i_0a,Om_0a,om_0a,tau_0a;                                /* orbital parameters just after SNA */
 double a_1,e_1,i_1,Om_1,om_1,tau_1;                                   /* orbital parameters just prior to SNB */
 double a_2,e_2,i_2,Om_2,om_2,tau_2;                                   /* orbital parameters just after to SNB */
 double jx_i,jy_i,jz_i,jx_f,jy_f,jz_f;                                                    /* temporary storage */

 double ttms1a,tthg1a,ttrgb1a,ttcheb1a,ttagb1a,tthems1a,tthehg1a,tthergb1a;     /* times of K type change [Myr] */
 double ttms1b,tthg1b,ttrgb1b,ttcheb1b,ttagb1b,tthems1b,tthehg1b,tthergb1b;

 int singsyst_head = 0;                                                         /* 0/1 not printed/printed system output header*/


 fp0=fopen("error.dat","w");
 fp2=fopen("info.dat","w");

 if(SINOUT==1) fp3=fopen("sb.dat","w");
 if(SSout==1) fp108=fopen("symb.dat","w");
 if(Xout2==1)
   fp130=fopen("DetWIND.dat","w");
 if(Xout1==1)
   fp131=fopen("DetMT.dat","w");
 if(BHCAT1==1) {
   fp141=fopen("bhcat1b.dat","w");
   fp142=fopen("bhcat2b.dat","w");
   fp143=fopen("bhcat3b.dat","w");
   fp144=fopen("bhcat4b.dat","w");
   fp145=fopen("bhcat5b.dat","w");
 }
 if(WDWDout==1) fp150=fopen("lisa.dat","w");
 if(CVout==1) fp151=fopen("cv.dat","w");
 if(BHCAT2==1) fp160=fopen("bhsys.dat","w");
 if(BHNSWD==1) fp170=fopen("bhnswd.dat","w");
 if(SN1A==1) fp180=fopen("sn1a.dat","w");
 if(AMCVN==1) fp181=fopen("amcvn.dat","w");
 if(OBCAT==1) fp190=fopen("OBbin.dat","w");
 if(COMPACT==1) fp200=fopen("compact.dat","w");
 if(BHSOUT==1) fp210=fopen("bhbh.spin.05.dat","w");
 if(CLUSTER==1) fp300=fopen("cluster.dat","w");
 if(LONGGRB==1) fp400=fopen("longgrb.dat","w");
 if(BINARY==0) fp500=fopen("sstars.dat","w");
 if(BINARY==1 && BINOUT==1) fp510=fopen("bstars.dat","w");
 if(WDWD==1) fp600=fopen("wdwd.dat","w");
 if(RRwd==1) fp700=fopen("RRwd.dat","w");
 if(Merger==1) fp800=fopen("Merger.dat","w");
 if(TranCE==1) fp900=fopen("Transient.dat","w");

 idum=(long int)-fmod((double)time(NULL),1000000.0);        /* inicjalizacja ran(3) liczba sekund z systemu */
 if(idum>=0.0 || idum<(-1.0e+06))
   nrerror("at start: idum must be negative and smaller then -1000000!\n");
 fprintf(fp2,"idum: %ld\n",idum);                           /* UWAGA!!! idum ma byc nie wiecej niz 6 cyfrowe */
 idum_run=idum;                                             /* record the initialization inside the program */
 fflush(fp2);


 max=num_tested;
 PP=1;                                              /* PP=0 no screen output, PP=1 -- detailed screen output */
 if(PP==1) {                                        /* PP=2 preset system evolved */
   idum=-955025;                                    /* tu wstawic wartosc z runu ktory chce sie sprawdzac */
   idum_run=idum;
   for(ii=0;ii<24496;ii++)                           /* iran: number of ran3() runs upto problematic place */
     if(Random==1) ran3(&idum);
     else ran2(&idum);
   max=1;
 }
 if(PP==2) {                                       /* read idum,iidd from a file called "preset.dat" */
   fp220=fopen("preset.bhbh.dat","r");
   max=0;
   while(fscanf(fp220,"%*d%*d")!=EOF)
     max++;
   rewind(fp220);
 }

 jj=1;
 ss0=ss1a=ss1b=ss1c=ss2a=ss2b=ss2c=ss3a=ss3b=ss3c=ss4a=ss4b=ss4c=0;
 ns1a=ns1b=ns1c=ns2a=ns2b=ns2c=ns3a=ns3b=ns3c=ns4a=ns4b=ns4c=0;
 hce1=hce2=hce3=0;
 cmt3=0;

 M_hook=M_hookf();                                          /* Initializing functions: order of calls important */
 M_HeF=M_HeFf();                                                                        /* sets critical masses */
 M_FGB=M_FGBf();
 coef_aa();                                                                     /* sets ZZ dependent parameters */
 coef_bb();
 init_SN();                                                                  /* initializing SN counting tables */
 maxMB=Mmaxconvf();
 minMB=Mminconvf();

 jx_0=jy_0=jz_0=jx_1=jy_1=jz_1=jx_2=jy_2=jz_2=0.0;

/* -------------------------------------- MAIN LOOP ----------------------------------------- */

 for(j=0;j<max;j++) {
//printf("j: %d\n",j);
   if(PP==2) {                                                                 /* evolve only requested systems */
     fscanf(fp220,"%ld%d",&idum,&pp);
     idum_run=idum;
     iidd=0;
     for(ii=0;ii<pp;ii++)
       if(Random==1) ran3(&idum);
       else ran2(&idum);
   }

   iidd_old=iidd;                                                 /* how many runs of run3() were done till now */
   if((j-jj*NOUT)==0) {                                   /* writes informations to info.dat every 1000 systems */
     rec_dat2(fp2,j,aicns,aicbh,badorb);
     fflush(fp2);
     jj++;
   }

   for(ll=0;ll<100;ll++)
     filledX[ll]=0;

   mt1H=mt2H=0;                                                /* initialize MT history tables */


   Tst=get_t(0.0,dtSFR);                                       /* sets SFR: either outbursting or continuos */
   t_hubble=hub_val-Tst;


   if(BINARY==0) {                           /* start two single stars */
     // Mzamsa=get_M(Mmina,Mmaxa);
     Mzamsa=get_M_linear(Mmina,Mmaxa);               // Linear distribution of initial masses
     Mzamsb=get_M(Mmina,Mmaxa);
     Ra0=Rzamsf(Mzamsa);                     /* radii of components on ZAMS [Rsun] */
     Rb0=Rzamsf(Mzamsb);
     if(RotF==1) {                           /* components initial angular spin velocities */
       wa=get_w1(Mzamsa,Ra0);
       wb=get_w1(Mzamsb,Rb0);
     }
     else if(RotF==2) {
       wa=get_w2(Mzamsa,Ra0);
       wb=get_w2(Mzamsb,Rb0);
     }
     a0=1000000.0;
     e0=0.0;
     q0=i0=Om0=om0=tau0=-1.0;
     inbinA=inbinB=inbin=0;
   }

   else if(BINARY==1) {                      /* start one binary */
     qtu1=0;
     if(SS==2) {                                   /* POP III stars: mock/test simulations */
       Mzamsa=get_gauss(101.4,40.4,50.0,200.0);    /* normal/gauss distribution in range 50-200 Msun */
       q0=get_q1(0.02,0.2);                        /* flat mass ratio */
       Mzamsb=q0*Mzamsa;
       e0=get_gauss(0.537,0.277,0.0,0.95);         /* normal/gauss distribution in range 0.0-0.95 */

       lmin=log(0.0001);                           /* minimum separation is 0.0001AU: 0.022Rsun */
       lmax=log(100.0);                            /* maximum separation is 100AU: 21,500Rsun */
       lxave=log(0.135);                           /* mean value of distr. -2.002 (a=0.135AU=29.0Rsun) */
       a0=get_gauss(lxave,1.862,lmin,lmax);        /* lognormal distribution: natural log of a0 */
       a0=pow(2.718282,a0);                        /* a0 [AU] */
       a0/=(0.00465047);                           /* a0 [Rsun] */

       Ra0=Rzamsf(Mzamsa);                         /* radii of components on ZAMS [Rsun] */
       Rb0=Rzamsf(Mzamsb);

       if(RotF==1) {                               /* components initial angular spin velocities */
         wa=get_w1(Mzamsa,Ra0);
         wb=get_w1(Mzamsb,Rb0);
       }
       else if(RotF==2) {
         wa=get_w2(Mzamsa,Ra0);
         wb=get_w2(Mzamsb,Rb0);
       }                              /* no preMS circ/sync, as eccentricities come from dynamical suimulations */

       i0=get_i(0.0,Pi);
       Om0=get_Om(0.0,2*Pi);
       om0=get_om(0.0,2*Pi);
       tau0=0.0;                                   /* czas przejscia przez perygeum = 0.0 jako czas odniesienia */

       inbinA=inbinB=inbin=1;

     }
     else if(IMF==1 && SS==0) {                /* Mzamsa from Kroupa broken power law, Mzamsb from q-distr. */
       while(1) {
         Mzamsa=get_M(Mmina,Mmaxa);            /* Mzamsa[Mmin-Mmax] and Mzamsb[0.08-Mmax] */
         qmin=0.08/Mzamsa;                     /* qmin allows secondary minimum mass of 0.08 - nuclear burning limit */
         if(MarkRat==1)
           q0=get_q1(qmin,1.0);                /* Rat=0.0:  flat q distribution */
         else if(MarkRat==2)
           q0=get_q2(qmin,1.0);                /* Rat!=0.0: q^(Rat) distribution */
         else if(MarkRat==3)
           q0=get_q3(qmin);                    /* Rat!=0.0: q=const + q^(Rat) distribution */
         else if(MarkRat==4) {
           if(qtu1==0) numberq=get_q1(0.0,1.0);
           qtu1=1;                             /* one numberq per one binary, otherwise more TWINS than predicted by qratio */
           if(numberq<=qratio)                 /* model for TWIN binaries: */
             q0=get_q1(qdivide,1.0);           /* qratio systems with flat q distr in range  (qdivide-1.0) */
           else                                /* model for REGULAR binaries: */
             q0=get_q1(qmin,qdivide);          /* 1-qratio systems with flat q distr in range (qmin-qdivide) */
         }
         Mzamsb=q0*Mzamsa;
         if(Mzamsb>=Mminb)
           break;
       }
     }
     else if(IMF==2 && SS==0) {                /* Mzamsa AND Mzamsb independently from Kroupa broken power law */
       Mzamsa=get_M(Mmina,Mmaxa);              /* Mzamsa[Mmina-Mmaxa] and Mzamsb[Mminb-Mmaxb] */
       Mzamsb=get_M(Mminb,Mmaxb);
       q0=Mzamsb/Mzamsa;
     }
     else if(IMF==3 && SS==0) {                /* Mzamsa AND Mzamsb from linear distributions, but with Mzamsb still folloing q */
       while(1) {
         Mzamsa=get_M_linear(Mmina,Mmaxa);     /* Mzamsa[Mmina-Mmaxa] and Mzamsb[Mminb-Mzamsa] */
         qmin=0.08/Mzamsa;                     /* qmin allows secondary minimum mass of 0.08 - nuclear burning limit */
         if(SanaOriginal==1)
           q0=get_M_linear(0.1,1.0);
         else
           q0=get_M_linear(qmin,1.0);
         Mzamsb=q0*Mzamsa;
         if(Mzamsb>=Mminb)
           break;
       }
     }
     else if(IMF==1 && SS==1) {                /* Mzamsa from Kroupa broken power law, Mzamsb from Sana et al. q-distr. */
       while(1) {
         Mzamsa=get_M(Mmina,Mmaxa);            /* Mzamsa[Mmin-Mmax] and Mzamsb[0.08-Mmax] */
         qmin=0.08/Mzamsa;                     /* qmin allows secondary minimum mass of 0.08 - nuclear burning limit */
         if(SanaOriginal==1)
           q0=get_SS(0.1,1.0,0.0);
         else
           q0=get_SS(qmin,1.0,0.0);
         Mzamsb=q0*Mzamsa;
         if(Mzamsb>=Mminb)
           break;
       }
     }

     if(SS==0) {
       e0=get_e(0.0,1.0);
       if(fabs(e0-1.0)<1.0e-06)          /* zabezpieczam aby e0 nie bylo nigdy 1.0 (skladniki maja obiegac sie) */
         e0=0.999999;
     }
     else if(SS==1) {
       e0=get_SS(0.0,0.9,-0.42);
       if(fabs(e0-1.0)<1.0e-06)          /* zabezpieczam aby e0 nie bylo nigdy 1.0 (skladniki maja obiegac sie) */
         e0=0.999999;
     }

     Ra0=Rzamsf(Mzamsa);                                                /* radii of components on ZAMS [Rsun] */
     Rb0=Rzamsf(Mzamsb);

     if(SS==0) {
       amin=(famin*Ra0+famin*Rb0)/(1.0-e0);      /* zaczynam losowac a0 od amin takiego aby w peryastronie dla orbity o */
       if(amin>amax) amin=amax/10.0;             /* making sure that for e=0.999999 or so, amin is not larger than amax */
       a0=get_a(amin,amax);                     /* amin i e0 dper bylo >= 2*Ra_ZAMS+2*Rb_ZAMS -moze powstac uklad kont. */
       per0=2.0*Pi*a0*sqrt(a0/(G*(Mzamsa+Mzamsb)));                                       /* initial orbital period [d] */
     }
     else if(SS==1) {
       per0=get_SS(0.15,5.5,-0.55);                                                 /* log10 of orbital period */
       per0=pow(10.0,per0);                                                         /* orbital period [day] */
       a0=pow(per0*per0*G*(Mzamsa+Mzamsb)/(4.0*Pi*Pi),(1.0/3.0));                   /* semi-major axis [Rsun] */
     }


     if(RotF==1) {                     /* components initial angular spin velocities */
       wa=get_w1(Mzamsa,Ra0);
       wb=get_w1(Mzamsb,Rb0);
     }
     else if(RotF==2) {
       wa=get_w2(Mzamsa,Ra0);
       wb=get_w2(Mzamsb,Rb0);
     }

     i0=get_i(0.0,Pi);
     Om0=get_Om(0.0,2*Pi);
     om0=get_om(0.0,2*Pi);
     tau0=0.0;                                   /* czas przejscia przez perygeum = 0.0 jako czas odniesienia */

     if(per0<PMStid && Tides==1 && SS==0) {                                         /* pre-MS circularization */
       e0=0.0;
       worb=sqrt(GGG*(Mzamsa+Mzamsb))*pow(a0,-1.5);                 /* mean orbital angular velocity [Myr^-1] */
       wa=wb=worb;                                                                  /* pre-MS synchronization */
     }
     if(per0<4.3 && Tides==1 && SS==1) {                                         /* pre-MS circularization */
       e0=0.0;
       worb=sqrt(GGG*(Mzamsa+Mzamsb))*pow(a0,-1.5);                 /* mean orbital angular velocity [Myr^-1] */
       wa=wb=worb;                                                                  /* pre-MS synchronization */
     }
     inbinA=inbinB=inbin=1;
   }

   q=q0;
   Ma=M0a=Mzamsa;
   Mb=M0b=Mzamsb;
   a=a0;
   e=e0;
   i=i0;
   Om=Om0;
   om=om0;
   tau=tau0;

   Vsm[0]=0.0; Vsm[1]=0.0; Vsm[2]=0.0;             /* final center of mass velocity */
   Vsm0[0]=0.0; Vsm0[1]=0.0; Vsm0[2]=0.0;          /* center of mass velocity after 1st SN */
   VsmA[0]=0.0; VsmA[1]=0.0; VsmA[2]=0.0;          /* center of mass velocity after A go SN */
   VsmB[0]=0.0; VsmB[1]=0.0; VsmB[2]=0.0;          /* center of mass velocity after B go SN */
   Vkick[0]=0.0; Vkick[1]=0.0; Vkick[2]=0.0;
   Vsa[0]=0.0; Vsa[1]=0.0; Vsa[2]=0.0;             /* component A after disruption spatial velocity */
   Vsb[0]=0.0; Vsb[1]=0.0; Vsb[2]=0.0;             /* component B after disruption spatial velocity */
   Vextr[0]=0.0; Vextr[1]=0.0; Vextr[2]=0.0;       /* final velocity of single star after its SN */
   Tmr=0.0;                                        /* merger, starting  time of a given formed system */


 /*--------------------------------- MAIN EVOLUTIONARY LOOP -----------------------------------------*/

   stop1=stop2=0;
   t=tbeg=tend=0.0;
   dta=delms*tmsf(Ma);
   dtb=delms*tmsf(Mb);

   Ra=Ra0;                                                             /* initializes first evolutionary step */
   Rb=Rb0;
   La=Lb=Ta=Tb=0.0;
   Rla=Rlb=0.0;
   Mca=Mcb=Mhea=Mheb=Mcoa=Mcob=Mprea=Mpreb=0.0;
   tvira=tvirb=tstarta=tstartb=ta_end=tb_end=tdis=0.0;
   Ka=Kb=Kpa=Kpb=1;
   KpgaA=KpgbA=KpgaB=KpgbB=0;
   MpgaA=MpgbA=MpgaB=MpgbB=apgA=apgB=epgA=epgB=tpgA=tpgB=MendaA=MendbB=0.0;
   MccoaA=McheaA=fracaA=MccobB=MchebB=fracbB=0.0;
   flaga=flagb=0;
   wda=wdb=0;
   WDout=SPCout=OBstar=0;
   ecssna=ecssnb=0;
   fraca=fracb=-1.0;                 /* -1 -- not applicable if not set in compactf1/2() functions of singl.c */
   fba=fbb=-1;
   sn=sna=snb=0;
   cmt1=0;
   atmp1=atmp2=0.0;
   mark1=0;
   mark10=mark11=-1;
   Iaold=Ibold=Iatmp=Ibtmp=-2.0;
   Maio=Mbio=dMcea=dMceb=dMceaold=dMcebold=0.0;
   dMwinda=dMwindb=dMmta=dMmtb=dMmta_old=dMmtb_old=dMtran=dMedd=0.0;
   dMalla=dMallb=0.0;
   dMatmorlofa=dMwindrlofa=dMwindacca=0.0;
   dMatmorlofb=dMwindrlofb=dMwindaccb=0.0;
   dMmsa=dMmsb=Fmsa=Fmsb=0.0;
   mark76a=mark76b=mark77a=mark77b=0;
   mark88a=mark88b=0;

   dMtha=dMthb=Ttha=Tthb=0.0;
   Mout1a=Mout1b=Mout2a=Mout2b=100.0;
   mark_outa=mark_outb=0;
   Macc0a=Macc0b=0.0;
   mark777a=mark777b=0;
   deca=decb=0;
   Lxmt=0.0;
   dtorb=0.0;
   maxred=0;
   aicmark=0;
   bhdone=0;
   mtflag1ab=mtflag1a=mtflag1b=mtflag2=mtflag3=flagbra=flagbrb=darwin=0;
   maxbr=0;
   dtbr=dtmt=0.0;
   rem=0;
   mttype=mttype_old=-1;
   mt=ce=0;
   texp=0.0;
   comdone=0;
   tlast=Mdonlast=Macclast=alast=0.0;
   Kdonlast=Kacclast=mttypelast=0;
   spinout=0;
   aspina=aspinb=-1.0;
   Mwda_0=Mwdb_0=-1.0;
   aspina0=aspina1=aspina2=aspinb0=aspinb1=aspinb2=-1.0;
   Mbha0=Mbha1=Mbha2=Mbhb0=Mbhb1=Mbhb2=-1.0;
   tbhenda0=tbhenda1=tbhenda2=tbhendb0=tbhendb1=tbhendb2=-1.0;
   Kace=Kbce=0;
   ace=Mace=Mbce=Race=Rbce=Mcace=Mcbce=Mheace=Mhebce=Mcoace=Mcobce=0.0;
   amcvn=0;
   Sia=Sib=-1;
   twra=awra=ewra=Mwra=Mcoma=twrb=awrb=ewrb=Mwrb=Mcomb=-1.0;

   ttms1a=tthg1a=ttrgb1a=ttcheb1a=ttagb1a=tthems1a=tthehg1a=tthergb1a=-1.0;
   ttms1b=tthg1b=ttrgb1b=ttcheb1b=ttagb1b=tthems1b=tthehg1b=tthergb1b=-1.0;


   dper=a*(1.0-e);                                               /* separation of stars at periastron */
   Rla=roche(Mb/Ma,dper);                                        /* minimal roche lobe sizes -- at peryastron */
   Rlb=roche(Ma/Mb,dper);

   evroute_new();                                                /* evolutionary history: cleans the storage memory */

                                                                 /* MAIN EVOLUTIONARY LOOP */
   while(stop1==0 && stop2==0 && t<t_hubble) {                   /* stop1=1 formation of interesting us binary */

     aicce1=aicce2=0;
     if(ce==1 && Kaold==13 && Ka==14) aicce1=1;                  /* AIC NS-->BH in CE was encountered: record below */
     if(ce==1 && Kbold==13 && Kb==14) aicce2=1;

     if(t==0){
       RmaxA=0;
       RmaxB=0;
     }

     else{
       if(RmaxA<Ra){
         RmaxA=Ra;
       }
       if(RmaxB<Rb){
         RmaxB=Rb;
       }
     }

     Maold=Ma;
     Mbold=Mb;                                                   /* actual mass before time step */
     Raold=Ra;
     Rbold=Rb;
     Kaold=Ka;
     Kbold=Kb;
     aold=a;
     eold=e;
     waold=wa;
     wbold=wb;
     Mheaold=Mhea;
     Mhebold=Mheb;
     Mcoaold=Mcoa;
     Mcobold=Mcob;
     get_SNcounters();
     deriv=0;
     Lxmt=dMtran=KTa=KTb=0.0;
     doce=merger=0;
     Mdisa=Mdisb=0.0;

     if(inbin==1) {                     /* timestep for binaries */

       if(sn==1)                        /* in case of SN in previous step, make dt very small, for proper functioning of */
         dt=1.0e-06;                                                                            /* orb_change() function */
       else {
         dt=min(dta,dtb);
         dt=min(dt,100.0);
       }

       if(dMmta>acc || dMmtb>acc)                                   /* step during MT X-ray phase or degenerate MT phase */
         dt=dtmt;                       /* singl() forced to move both components by dtmt, disragarding change in M or R */
       else if(((Ka==13 || Ka==14) && Kb<10) || ((Kb==13 || Kb==14) && Ka<10))  /* max step at WIND X-ray phase: 100 Myr */
         dt=min(100.0,dt);
       else if(Ka>=10 && Kb>=10 && sn!=1) {                  /* timestep for detached COB is in range (10 yrs -- 100Myr) */
         Tmr=tmerge(Ma,Mb,a,e);                                         /* sn!=1: do not make big 100 Myr step after SN: */
         dt=min(100.0,0.01*Tmr);                                             /* so there is corect output to compact.dat */
         dt=max(0.00001,dt);
         if(Tmr<0.00002 && (Ka==13 || Ka==14) && (Kb==13 || Kb==14)) {          /* if merger time of two compact objects */
           stop1=1;                                                                            /* is smaller than 20 yrs */
           break;                                                                         /* STOP: evolution, COB merger */
         }
         else if(Tmr<0.00002) {                                  /* bring larger component to contact, system very tight */
           dt=1.0e-12;                                                                 /* ODE may have a problem with it */
           e=0.0;                                                          /* larger component may be WD in WD-NS binary */
           atmp1=0.99*Aroche(Mb/Ma,Ra);                            /* then NS becomes donor, and program halts evolution */
           atmp2=0.99*Aroche(Ma/Mb,Rb);                                            /* before calling MT function dMmtf() */
           a=max(atmp1,atmp2);
         }
       }

       if(mt==1) {                                                                    /* this set in last detailed MT step */
         dtbr=dt;
         maxbr=10;
       }
       if(maxbr>0 && mt==0) {                                                           /* if in next step system detaches */
         dt=min(dt,dtbr);                                             /* keep dt="timestep from last MT" for 10 next steps */
         maxbr-=1;                                                     /* or keep it smaller if set above to smaller value */
       }

       if(dtorb>acc)              /* dcreased dt in case of big total binary angular mom. change */
         dt=dtorb;

     }
     else {                       /* timestep for single stars */
       if(sn==1)
         dt=1.0e-06;
       else if(Ka<10 || Kb<10)
         dt=min(dta,dtb);
       else
         dt=100.0;
     }


     if(mt==1 && amcvn==0 && (((Ka==10 || Ka==11 || Ka==12 || Ka==16 || Ka==17) && (Kb==7 || Kb==8 || Kb==9 || Kb==10 || Kb==17))
       || ((Kb==10 || Kb==11 || Kb==12 || Kb==16 || Kb==17) && (Ka==7 || Ka==8 || Ka==9 || Ka==10 || Ka==17)))) {
       amcvn=1;
       t11=Tst+t;
       Ka11=Ka;
       Kb11=Kb;
       Ma11=Ma;
       Mb11=Mb;
       a11=a;
       e11=e;
       Ra11=Ra;
       Rb11=Rb;
     }

     if(mt==1 && amcvn==1 && (((Ka==10 || Ka==11 || Ka==12 || Ka==16 || Ka==17) && (Kb==7 || Kb==8 || Kb==9 || Kb==10 || Kb==17))
       || ((Kb==10 || Kb==11 || Kb==12 || Kb==16 || Kb==17) && (Ka==7 || Ka==8 || Ka==9 || Ka==10 || Ka==17)))) {
       t22=Tst+t;
       Ka22=Ka;
       Kb22=Kb;
       Ma22=Ma;
       Mb22=Mb;
       a22=a;
       e22=e;
       Ra22=Ra;
       Rb22=Rb;
     }


     if(BHCAT1==1) {
       if(t<TIMEOUT1) dt=min(dt,(TIMEOUT1-t+0.00001));
       else if(t<TIMEOUT2) dt=min(dt,(TIMEOUT2-t+0.00001));
       else if(t<TIMEOUT3) dt=min(dt,(TIMEOUT3-t+0.00001));
       else if(t<TIMEOUT4) dt=min(dt,(TIMEOUT4-t+0.00001));
       else if(t<TIMEOUT5) dt=min(dt,(TIMEOUT5-t+0.00001));
     }

     if(dta<0.0 || dtb<0.0 || dt<0.0) {
       fprintf(fp0,"error: wrong sign of timestep in binary.c: dta: %g, dtb: %g, dt: %g, %d\n",dta,dtb,dt,iidd_old);
       fflush(fp0);
     }

     tbeg=t;
     tenda=tendb=t+dt;

     copySin(Mzamsa,M0a,Ma,Ka,tbeg,tvira,tenda,La,Ra,Mca,Mhea,Mcoa,flaga,dta,Mprea,Kpa,tstarta,fraca,dMalla,ecssna,
             &sa1,&sa2,&sa3,&sa4,&sa5,&sa6,&sa7,&sa8,&sa9,&sa10,&sa11,&sa12,&sa13,&sa14,&sa15,&sa16,&sa17,&sa18,&sa19,&sa20);
     copySin(Mzamsb,M0b,Mb,Kb,tbeg,tvirb,tendb,Lb,Rb,Mcb,Mheb,Mcob,flagb,dtb,Mpreb,Kpb,tstartb,fracb,dMallb,ecssnb,
             &sb1,&sb2,&sb3,&sb4,&sb5,&sb6,&sb7,&sb8,&sb9,&sb10,&sb11,&sb12,&sb13,&sb14,&sb15,&sb16,&sb17,&sb18,&sb19,&sb20);

     singl(&Mzamsa,&M0a,&Ma,&Ka,&tbeg,&tvira,&tenda,&La,&Ra,&Mca,&Mhea,&Mcoa,&flaga,&dta,&Mprea,&Kpa,&tstarta,&fraca,dMalla,deriv,Mout1a,Mout2a,&ecssna,&fba,1,Fmsa);
     singl(&Mzamsb,&M0b,&Mb,&Kb,&tbeg,&tvirb,&tendb,&Lb,&Rb,&Mcb,&Mheb,&Mcob,&flagb,&dtb,&Mpreb,&Kpb,&tstartb,&fracb,dMallb,deriv,Mout1b,Mout2b,&ecssnb,&fbb,1,Fmsb);

     if(fabs(t+dt-tenda)>10.0*acc && fabs(t+dt-tendb)>10.0*acc) { /* redo single evol. with smaller time step for A,B */
       dt=min(tenda-t,tendb-t);                                   /* must be here, and not from copies below */
       copySin(sb1,sb2,sb3,sb4,sb5,sb6,sb7,sb8,sb9,sb10,sb11,sb12,sb13,sb14,sb15,sb16,sb17,sb18,sb19,sb20,
               &Mzamsb,&M0b,&Mb,&Kb,&tbeg,&tvirb,&tendb,&Lb,&Rb,&Mcb,&Mheb,&Mcob,&flagb,&dtb,&Mpreb,&Kpb,&tstartb,&fracb,&dMallb,&ecssnb);
       copySin(sa1,sa2,sa3,sa4,sa5,sa6,sa7,sa8,sa9,sa10,sa11,sa12,sa13,sa14,sa15,sa16,sa17,sa18,sa19,sa20,
               &Mzamsa,&M0a,&Ma,&Ka,&tbeg,&tvira,&tenda,&La,&Ra,&Mca,&Mhea,&Mcoa,&flaga,&dta,&Mprea,&Kpa,&tstarta,&fraca,&dMalla,&ecssna);
       tenda=tendb=t+dt;                                    /* new values (must be here), not from the above copies!!! */
       singl(&Mzamsb,&M0b,&Mb,&Kb,&tbeg,&tvirb,&tendb,&Lb,&Rb,&Mcb,&Mheb,&Mcob,&flagb,&dtb,&Mpreb,&Kpb,&tstartb,&fracb,dMallb,deriv,Mout1b,Mout2b,&ecssnb,&fbb,0,Fmsb);
       singl(&Mzamsa,&M0a,&Ma,&Ka,&tbeg,&tvira,&tenda,&La,&Ra,&Mca,&Mhea,&Mcoa,&flaga,&dta,&Mprea,&Kpa,&tstarta,&fraca,dMalla,deriv,Mout1a,Mout2a,&ecssna,&fba,0,Fmsa);
     }
     else if(fabs(t+dt-tenda)>10.0*acc) {  /* redo the single evol. with smaller time step for star B */
       dt=tenda-t;
       copySin(sb1,sb2,sb3,sb4,sb5,sb6,sb7,sb8,sb9,sb10,sb11,sb12,sb13,sb14,sb15,sb16,sb17,sb18,sb19,sb20,
               &Mzamsb,&M0b,&Mb,&Kb,&tbeg,&tvirb,&tendb,&Lb,&Rb,&Mcb,&Mheb,&Mcob,&flagb,&dtb,&Mpreb,&Kpb,&tstartb,&fracb,&dMallb,&ecssnb);
       tendb=t+dt;                     /* new value (must be here), not from the above copy!!! */
       singl(&Mzamsb,&M0b,&Mb,&Kb,&tbeg,&tvirb,&tendb,&Lb,&Rb,&Mcb,&Mheb,&Mcob,&flagb,&dtb,&Mpreb,&Kpb,&tstartb,&fracb,dMallb,deriv,Mout1b,Mout2b,&ecssnb,&fbb,0,Fmsb);
     }
     else if(fabs(t+dt-tendb)>10.0*acc) {  /* redo the single evol. with smaller time step for star A */
       dt=tendb-t;
       copySin(sa1,sa2,sa3,sa4,sa5,sa6,sa7,sa8,sa9,sa10,sa11,sa12,sa13,sa14,sa15,sa16,sa17,sa18,sa19,sa20,
               &Mzamsa,&M0a,&Ma,&Ka,&tbeg,&tvira,&tenda,&La,&Ra,&Mca,&Mhea,&Mcoa,&flaga,&dta,&Mprea,&Kpa,&tstarta,&fraca,&dMalla,&ecssna);
       tenda=t+dt;                     /* new value (must be here), not from the above copy!!! */
       singl(&Mzamsa,&M0a,&Ma,&Ka,&tbeg,&tvira,&tenda,&La,&Ra,&Mca,&Mhea,&Mcoa,&flaga,&dta,&Mprea,&Kpa,&tstarta,&fraca,dMalla,deriv,Mout1a,Mout2a,&ecssna,&fba,0,Fmsa);
     }

     if(Ka!=Kaold || Kb!=Kbold) {                /* change of stellar type management */
       dt1=1.0e-12;                              /* if star has changed its type, move system a bit forward */
       tbeg=t+dt;                                /* from t+dt to t+dt+dt1: to get components proper parameters */
       tenda=tendb=t+dt+dt1;                     /* but do not change old values!!!, also do not apply mass loss/gain here */
       singl(&Mzamsa,&M0a,&Ma,&Ka,&tbeg,&tvira,&tenda,&La,&Ra,&Mca,&Mhea,&Mcoa,&flaga,&dta,&Mprea,&Kpa,&tstarta,&fraca,0.0,deriv,Mout1a,Mout2a,&ecssna,&fba,0,Fmsa);
       singl(&Mzamsb,&M0b,&Mb,&Kb,&tbeg,&tvirb,&tendb,&Lb,&Rb,&Mcb,&Mheb,&Mcob,&flagb,&dtb,&Mpreb,&Kpb,&tstartb,&fracb,0.0,deriv,Mout1b,Mout2b,&ecssnb,&fbb,0,Fmsb);
       dt=dt+dt1;                                /* dt is bigger now by dt1 */
       mt1H=0;                                   /* controls for evolutionary output */
       mt2H=0;
     }
     t+=dt;                                      /* t: time at the end of current evolutionary step: [t-dt,t] */



     if(Ka==0 || Ka==1) ttms1a=t;                /* set times of various evolutionary stages [Myr] */
     else if(Ka==2) tthg1a=t;
     else if(Ka==3) ttrgb1a=t;
     else if(Ka==4) ttcheb1a=t;
     else if(Ka==5 || Ka==6) ttagb1a=t;
     else if(Ka==7) tthems1a=t;
     else if(Ka==8) tthehg1a=t;
     else if(Ka==9) tthergb1a=t;

     if(Kb==0 || Kb==1) ttms1b=t;
     else if(Kb==2) tthg1b=t;
     else if(Kb==3) ttrgb1b=t;
     else if(Kb==4) ttcheb1b=t;
     else if(Kb==5 || Kb==6) ttagb1b=t;
     else if(Kb==7) tthems1b=t;
     else if(Kb==8) tthehg1b=t;
     else if(Kb==9) tthergb1b=t;


     if(mt==1 && Ka==1 && dMalla<-1.0e-15) {      /* RLOF on, star A on MS, star A accrets (negative sign) */
       dMmsa+=(-dMalla*dt);                       /* accreted mass [Msun] */
       if(mark76a==0) {
         Mina=Ma;                                 /* accretor mass at the onset of RLOF */
         mark76a=1;
       }
     }
     if(mt==1 && Kb==1 && dMallb<-1.0e-15) {      /* RLOF on, star B on MS, star B accrets (negative sign) */
       dMmsb+=(-dMallb*dt);                       /* accreted mass [Msun] */
       if(mark76b==0) {
         Minb=Mb;                                 /* accretor mass at the onset of RLOF */
         mark76b=1;
       }
     }
     if(Ka==2 && mark77a==0) {                    /* end of MS */
       if(fabs(Mina)<acc) Fmsa=0.0;
       else Fmsa=dMmsa/Mina;                      /* fraction of star mass that was accreted in RLOF during MS */
       mark77a=1;
     }
     if(Kb==2 && mark77b==0) {                    /* end of MS */
       if(fabs(Minb)<acc) Fmsb=0.0;
       else Fmsb=dMmsb/Minb;                      /* fraction of star mass that was accreted in RLOF during MS */
       mark77b=1;
     }


     if(Ka>=7 && Ka<=9 && Kaold>=7 && Kaold<=9 && mark88a==0) {   /* record formation of helium star A */
       twra=t;
       awra=a;
       ewra=e;
       Mwra=Ma;
       Mcoma=Mb;
       mark88a=1;
     }
     if(Kb>=7 && Kb<=9 && Kbold>=7 && Kbold<=9 && mark88b==0) {   /* record formation of helium star B */
       twrb=t;
       awrb=a;
       ewrb=e;
       Mwrb=Mb;
       Mcomb=Ma;
       mark88b=1;
     }


     if(Ka==18 || Kb==18) {                   /* either A or B exploded as pair instability SN: no remnant, no binary */
       inbin=0;
       a=1.0e+50;
       e=0.0;
       if(Ka==18) {
        ta_end=t;
        Ra=RFa=Ma=La=wa=0.0;
       }
       else {
        tb_end=t;
        Rb=RFb=Mb=Lb=wb=0.0;
       }
       if(PP==1) //printf("pair instability SN: stop evolution (Ka: %d, Kb: %d)\n",Ka,Kb);
       break;                                 /* break out of while loop; go to next binary */
     }

     if(Kaold==3 && Ka==4 && mt==1 && Ra>1.2*Raold) {
       fprintf(fp0,"warning: R jump star A (%d %d  %.3f %.3f  %.1f %.1f %.1f  %f %f  %d %d)\n",Ka,Kb,Ma,Mb,a,Tst,t,Raold,Ra,idum_run,iidd_old);
       fflush(fp0);
     }
     if(Kbold==3 && Kb==4 && mt==1 && Rb>1.2*Rbold) {
       fprintf(fp0,"warning: R jump star B (%d %d  %.3f %.3f  %.1f %.1f %.1f  %f %f  %d %d)\n",Ka,Kb,Ma,Mb,a,Tst,t,Rbold,Rb,idum_run,iidd_old);
       fflush(fp0);
     }


     if(Ka==14 && aspina<-0.5) {                                 /* initialize spin/mass of component A BH */
       if(BHSPIN==1)
         aspina0=aspina=aspin_init;                              /* set by hand in sinbin.h */
       else if(BHSPIN==2)
         aspina0=aspina=bhspininit1(Mcoa);                       /* set in singl.c by bhspininit1() */
       else if(BHSPIN==3)
         aspina0=aspina=bhspininit2(Mcoa);                       /* set in singl.c by bhspininit2() */
       else if(BHSPIN==4)
         aspina0=aspina=bhspininit3(Mcoa);                       /* set in singl.c by bhspininit3() */
       else if(BHSPIN==10)
         aspina0=aspina=bhspininit10();                          /* set in singl.c by bhspininit10() */
       Mbha0=Ma;
       tbhenda0=t;
     }
     if(Kb==14 && aspinb<-0.5) {
       if(BHSPIN==1)
         aspinb0=aspinb=aspin_init;
       else if(BHSPIN==2)
         aspinb0=aspinb=bhspininit1(Mcob);
       else if(BHSPIN==3)
         aspinb0=aspinb=bhspininit2(Mcob);
       else if(BHSPIN==4)
         aspinb0=aspinb=bhspininit3(Mcob);
       else if(BHSPIN==10)
         aspinb0=aspinb=bhspininit10();
       Mbhb0=Mb;
       tbhendb0=t;
     }
     if(BHSPIN>0 && ((Kaold==14 && Ka==14) || (Kbold==14 && Kb==14)) && (mt==1 || ce==1 || (Ka==14 && Kb>=0 && Kb<=9) || (Kb==14 && Ka>=0 && Ka<=9))) {
       if(Ka==14) {                                                /* spin up (potentially accreting) BH */
         if(mt==1) {
           Mrest=-dMmta*dt;                                        /* RLOF: rest mass (positive) of accreted material [Msun] */
           Ma=Maold;                                               /* use pre-accretion values of M: Mold */
         }
         else if(ce==1) {
           Mrest=dMcea-dMceaold;                                   /* CE: rest mass (positive) of accreted material [Msun] */
           Ma-=Mrest;                                              /* use pre-accretion values of M: Mold not good here */
         }                                                         /* it was already changed to new: post-CE mass */
         else if(Kb>=0 && Kb<=9) {                                 /* wind accretion */
           Mrest=Ma-Maold;
           Ma=Maold;
         }
         else
           fprintf(fp0,"error before spin_evol() star A: %d %d\n",idum_run,iidd_old);
         spin_evol(Ka,&Ma,&aspina,&Mdisa,Mrest);                   /* changes BH spin and mass: */
         if(mt==1) {
           aspina2=aspina;
           Mbha2=Ma;
           tbhenda2=t;
         }
         else {
           Mdisa=0.0;                                              /* Mdis -- is not needed in case of CE */
           aspina1=aspina;                                         /* M is going to be smaller than one calculated in singl() */
           Mbha1=Ma;                                               /* since it is grav. mass now not rest mass */
           tbhenda1=t;
         }
       }
       else {
         if(mt==1) {
           Mrest=-dMmtb*dt;
           Mb=Mbold;
         }
         else if(ce==1) {
           Mrest=dMceb-dMcebold;
           Mb-=Mrest;
         }
         else if(Ka>=0 && Ka<=9) {                                 /* wind accretion */
           Mrest=Mb-Mbold;
           Mb=Mbold;
         }
         else
           fprintf(fp0,"error before spin_evol() star B: %d %d\n",idum_run,iidd_old);
         spin_evol(Kb,&Mb,&aspinb,&Mdisb,Mrest);
         if(mt==1) {
           aspinb2=aspinb;
           Mbhb2=Mb;
           tbhendb2=t;
         }
         else {
           Mdisb=0.0;
           aspinb1=aspinb;
           Mbhb1=Mb;
           tbhendb1=t;
         }
       }
     }                                                                              /* write out at the formation of DCO */
     if(BHSOUT==1 && spinout==0 && (((Ka==14 && Kaold==14) && Kbold>9) || ((Kb==14 && Kbold==14) && Kaold>9))) {
       spinout=1;
       fprintf(fp210,"%f  %f %f %f %f %f %f  %d %d  %f %f %f  %f %f %f %f %f %f  %f %f %f %f %f %f  %f %f %f %f %f %f  %f %f  %d %d  %d %s\n",
               t,aspina,aspinb,Ma,Mb,a,e,Ka,Kb,i0,iA,iB,aspina0,aspina1,aspina2,aspinb0,aspinb1,aspinb2,
               Mbha0,Mbha1,Mbha2,Mbhb0,Mbhb1,Mbhb2,tbhenda0,tbhenda1,tbhenda2,tbhendb0,tbhendb1,tbhendb2,
               ta_end,tb_end,idum_run,iidd_old,nevroute,evroute);
       fflush(fp210);
     }


     if((Kaold==10 || Kaold==17) && Ka==7) {          /* accretion induced transition of He or hybrid WD -> ZAMS He star */
       ta_end=0.0;
       wda=0;
     }
     if((Kbold==10 || Kbold==17) && Kb==7) {
       tb_end=0.0;
       wdb=0;
     }

     if(ce==1)                                   /* finish recording CE */
       evroute_add("%d-%d) ",Ka,Kb);

     aicmark=0;
     if(Kaold==12 && Ka==13) {                   /* AIC: WD->NS due to MT, AIC taken care of in singl() */
       aicns++;
       aicmark=1;
       Maio=Ma;                                  /* records the formation mass of NS */
       sna=0;                                    /* send system to explode(): mass loss/potential kick */
       evroute_add("AICNS1 ");                   /* record in evolutionary description */
       if(SN1A==1) {
         typwd=1;
         fprintf(fp180,"%d  %f %f  %d %d  %f %f %f %f  %f %f  -1 -1 -1 -1 -1 -1   %f %f %f %f %f %f  %f %f %f %f    -1 -1 -1 -1 -1 -1 -1  %d %d  %d %s\n",
                 typwd,t+Tst,Tst,Kaold,Kb,Maold,Mb,a,e,Raold,Rb,
                 M0a,M0b,Mzamsa,Mzamsb,a0,e0,Mwda_0,Mwdb_0,ta_end,tb_end,
                 idum_run,iidd_old,nevroute,evroute);
         fflush(fp180);
       }
     }
     if(Kbold==12 && Kb==13) {
       aicns++;
       aicmark=1;
       Mbio=Mb;
       snb=0;                                    /* send system to explode(): mass loss/potential kick */
       evroute_add("AICNS2 ");                   /* record in evolutionary description */
       if(SN1A==1) {
         typwd=1;
         fprintf(fp180,"%d  %f %f  %d %d  %f %f %f %f  %f %f  -1 -1 -1 -1 -1 -1   %f %f %f %f %f %f  %f %f %f %f    -1 -1 -1 -1 -1 -1 -1  %d %d  %d %s\n",
                 typwd,t+Tst,Tst,Kbold,Ka,Mbold,Ma,a,e,Rbold,Ra,
                 M0b,M0a,Mzamsb,Mzamsa,a0,e0,Mwdb_0,Mwda_0,tb_end,ta_end,
                 idum_run,iidd_old,nevroute,evroute);
         fflush(fp180);
       }
     }

     if((Kaold==13 && Ka==14) || aicce1==1) {    /* AIC: NS->BH due to MT: M,K set in single() */
       aicbh++;                                  /* AIC: NS->BH due to CE, M,K set in comm_env() functions */
       evroute_add("AICBH1 ");                   /* record in evolutionary description */
     }
     if((Kbold==13 && Kb==14) || aicce2==1) {
       aicbh++;
       evroute_add("AICBH2 ");
     }

     if(Kaold!=15 && Ka==15) {                   /* either WD blow up or SN Ia from regular star: same effect */
       if(SN1A==1) {
         typwd=2;
         fprintf(fp180,"%d  %f %f  %d %d  %f %f %f %f  %f %f  -1 -1 -1 -1 -1 -1   %f %f %f %f %f %f  %f %f %f %f    -1 -1 -1 -1 -1 -1 -1  %d %d  %d %s\n",
                 typwd,t+Tst,Tst,Kaold,Kb,Ma,Mb,a,e,Raold,Rb,
                 M0a,M0b,Mzamsa,Mzamsb,a0,e0,Mwda_0,Mwdb_0,ta_end,tb_end,
                 idum_run,iidd_old,nevroute,evroute);
         fflush(fp180);
       }
       sna=1;
       inbin=0;
       a=1.0e+50;
       e=0.0;
       dta=1.0e+50;
       if(fabs(ta_end)<acc) ta_end=t;     /* if it was WD (e.g. K=11) exploded in SN Ia (K=15); keep the end time of WD */
       Ra=RFa=Ma=La=wa=0.0;
     }
     if(Kbold!=15 && Kb==15) {
       if(SN1A==1) {
         typwd=2;
         fprintf(fp180,"%d  %f %f  %d %d  %f %f %f %f  %f %f  -1 -1 -1 -1 -1 -1   %f %f %f %f %f %f  %f %f %f %f    -1 -1 -1 -1 -1 -1 -1  %d %d  %d %s\n",
                 typwd,t+Tst,Tst,Kbold,Ka,Mb,Ma,a,e,Rbold,Ra,
                 M0b,M0a,Mzamsb,Mzamsa,a0,e0,Mwdb_0,Mwda_0,tb_end,ta_end,
                 idum_run,iidd_old,nevroute,evroute);
         fflush(fp180);
       }
       snb=1;
       inbin=0;
       a=1.0e+50;
       e=0.0;
       dtb=1.0e+50;
       if(fabs(tb_end)<acc) tb_end=t;
       Rb=RFb=Mb=Lb=wb=0.0;
     }
     if(inbin==0) {                                          /* clean markers/variables needed only for binary evolution */
       mt=ce=doce=merger=mttype=0;
       dMmta=dMmtb=0.0;
     }
     if(Ka==15 && Kb==15) {                                                                           /* no binary left! */
       stop2=1;
       break;
     }


     if(ce==0 && Ka<10) {                                            /* wind mass loss rate [Msun/Myr] (positive number) */
       dMwinda=(Maold-Ma-dMmta*dt)/dt;                                            /* not calculated after SN or CE event */
       if(dMwinda<acc) dMwinda=0.0;                                                          /* to avoid numerical noise */
     }
     else
       dMwinda=0.0;
     if(ce==0 && Kb<10) {
       dMwindb=(Mbold-Mb-dMmtb*dt)/dt;
       if(dMwindb<acc) dMwindb=0.0;
     }
     else
       dMwindb=0.0;


     if(Kaold==-1 && (Ka==13 || Ka==14)) RFa=0.01;                  /* CO core which will go SN, revealed by MT */
     else RFa=Ra;                                                /* radii approximated by WD radii of 0.01 Rsun */
     if(Kbold==-1 && (Kb==13 || Kb==14)) RFb=0.01;         /* for any other star type its actual radius is used */
     else RFb=Rb;

                                                    /* compact objects (WD/NS/BH): no spin evolution allowed!!! */
     if(Ka>=10) wa=0.0;                                  /* keep spins of WD/NS/BH equal to zero (neglect them) */
     if(Kb>=10) wb=0.0;


           /* check if system merged after SN: if any star went over its Roche lobe or NS/BH cuts through comp. */
     if(sn==1 && (((RFa+RFb)>=a*(1.0-e) || RFa>=2.0*Rla || RFb>=2.0*Rlb) || Is==-1)) {             /* after SN: */
       if(fabs(ta_end)<acc) ta_end=t;                          /* set evolutionary end time of A if not yet set */
       if(fabs(tb_end)<acc) tb_end=t;                          /* set evolutionary end time of B if not yet set */
       if((Ka>=7 && Ka<=9) || (Kb>=7 && Kb<=9))                                             /* SN helium merger */
         stop1=1;
       else if(((Ka>=10 && Ka<=12) || Ka==16 || Ka==17) || ((Kb>=10 && Kb<=12) || Kb==16 || Kb==17))
       	 stop1=1;                                                                         /* SN WD+NS/BH merger */
       else if((Ka==13 || Ka==14) && (Kb==13 || Kb==14))                            /* SN compact object merger */
         stop1=1;
       else if(Ka==0 || Ka==1 || Kb==0 || Kb==1)                                           /* merger of MS+NS/BH */
         stop2=1;
       else                                                                   /* merger which do not interest us */
         stop2=1;
       typce=6;                                         /* components parameters passed to record in Merger file */
       break;
     }
                                    /* if one or two of the star(s) overfill Roche Lobe after CE merger assumed: */
     else if(inbin==1 && ce==1 && (RFa>=Rla || RFb>=Rlb)) {                         /* CE inspiral didn't finish */
       if(fabs(ta_end)<acc) ta_end=t;                           /* set evolutionary end time of A if not yet set */
       if(fabs(tb_end)<acc) tb_end=t;                           /* set evolutionary end time of B if not yet set */

       if(SN1A==1) {
         if((Kace>9 && (Kbce==2 || Kbce==3 || Kbce==5 || Kbce==6 || Kbce==8 || Kbce==9) && (Mace+Mcbce)>1.38) ||
            (Kbce>9 && (Kace==2 || Kace==3 || Kace==5 || Kace==6 || Kace==8 || Kace==9) && (Mbce+Mcace)>1.38) ||
            ((Kace==2 || Kace==3 || Kace==5 || Kace==6 || Kace==8 || Kace==9) &&
             (Kbce==2 || Kbce==3 || Kbce==5 || Kbce==6 || Kbce==8 || Kbce==9) && (Mcace+Mcbce)>1.38)) {
           typwd=6;                                  /* Livio&Riess2003 scenario: WD + degenerate core merge in CE */
           fprintf(fp180,"%d  %f %f  %d %d  %f %f %f %f  %f %f  %f %f %f %f %f %f   %f %f %f %f %f %f  %f %f %f %f    %d %d %f %f %f %f %f  %d %d  %d %s\n",
                   typwd,t+Tst,Tst,Kace,Kbce,Mace,Mbce,ace,ece,Race,Rbce,Mheace,Mcoace,Mhebce,Mcobce,Mcace,Mcbce,
                   M0a,M0b,Mzamsa,Mzamsb,a0,e0,Mwda_0,Mwdb_0,ta_end,tb_end,
                   Ka,Kb,Ma,Mb,Ra,Rb,a,idum_run,iidd_old,nevroute,evroute);
           fflush(fp180);
         }
       }
       if(Kaold==-1) {                                         /* K=-1 denote CO core, which came from the giant */
         Ka=7;                                    /* to count these as Helium mergers I assume core type to be 7 */
         Ma=Maold;            /* to get mass of CO core and not mass of NS/BH which maight have formed out of it */
         Ra=RFa;                                                 /* radius of CO core is assumed to be WD radius */
         undo_SNcounters();                                                                  /* merger before SN */
       }
       if(Kbold==-1) {                                         /* K=-1 denote CO core, which came from the giant */
         Kb=7;                                    /* to count these as Helium mergers I assume core type to be 7 */
         Mb=Mbold;            /* to get mass of CO core and not mass of NS/BH which maight have formed out of it */
         Rb=RFb;                                                 /* radius of CO core is assumed to be WD radius */
         undo_SNcounters();                                                                  /* merger before SN */
       }
       if((Ka>=7 && Ka<=9 && (Kb==13 || Kb==14)) || (Kb>=7 && Kb<=9  && (Ka==13 || Ka==14)))    /* helium merger */
         stop1=1;
       else if(((Ka==13 || Ka==14) && ((Kb>=10 && Kb<=14) || Kb==16 || Kb==17)) || ((Kb==13 || Kb==14) && ((Ka>=10 && Ka<=14) || Ka==16 || Ka==17))) {
         stop1=1;                                                  /* A:NS,BH--B:WD,NS,BH or B:NS,BH--A:WD,NS,BH */
         fprintf(fp0,"error: should not be here 10: %d\n",iidd_old);              /* STOP: compact object merger */
       }
       else                                                                   /* merger which do not interest us */
         stop2=1;
       if(TranCE==1) {
         typce=1;                                       /* components parameters passed to record in Merger file */
         fprintf(fp900,"  0\n");
         break;
       }
     }
     else if(TranCE==1 && inbin==1 && ce==1 && RFa<Rla && RFb<Rlb) {           /* CE inspiral finish succesfully */
       fprintf(fp900,"  1\n");
     }
     else if(inbin==1 && mttype!=4 && (RFa>=2.0*Rla || RFb>=2.0*Rlb)) {    /* star too big; do nothing at the MT */
       if(((Ka>=13 && Ka<=14 && Kb>=2 && Kb<=6) || (Kb>=13 && Kb<=14 && Ka>=2 && Ka<=6))) {         /* He merger */
         if(Kb==2 || Kb==3 || Kb==4) Kb=7;                    /* tells recording function how to classify merger */
         else if(Kb==5 || Kb==6) Kb=8;
         if(Ka==2 || Ka==3 || Ka==4) Ka=7;
         else if(Ka==5 || Ka==6) Ka=8;
         stop1=1;
       }
       else if((Ka>=13 && Ka<=14 && Kb>=7 && Kb<=9) || (Kb>=13 && Kb<=14 && Ka>=7 && Ka<=9))        /* He merger */
         stop1=1;
       else if(Ka>=10 && Kb>=10) {                                                               /* Remnant merger */
         if(SN1A==1 && (Ka==10 || Ka==11 || Ka==12 || Ka==16 || Ka==17) && (Kb==10 || Kb==11 || Kb==12 || Kb==16 || Kb==17)) {
           typwd=3;
           fprintf(fp180,"%d  %f %f  %d %d  %f %f %f %f  %f %f  -1 -1 -1 -1 -1 -1   %f %f %f %f %f %f  %f %f %f %f    -1 -1 -1 -1 -1 -1 -1  %d %d  %d %s\n",
                  typwd,t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,Ra,Rb,
                  M0a,M0b,Mzamsa,Mzamsb,a0,e0,Mwda_0,Mwdb_0,ta_end,tb_end,
                  idum_run,iidd_old,nevroute,evroute);
           fflush(fp180);
         }
         stop1=1;
       }
       else
         stop2=1;                                                                           /* all other mergers */
       if(TranCE==1) {
         typce=3;
         fprintf(fp900,"%d  %f %f  %d %d  %f %f %f %f  %f %f  %f %f %f %f %f %f   %f %f %f %f %f %f  %d %d  %d %s ",
                 4,t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,Ra,Rb,Mhea,Mcoa,Mheb,Mcob,Mca,Mcb,
                 M0a,M0b,Mzamsa,Mzamsb,a0,e0,idum_run,iidd_old,nevroute,evroute);
         fprintf(fp900," 0 ");
         fprintf(fp900,"\n");
         fflush(fp900);
         break;
       }
     }

     if(((Ka>=10 && Ka<=12) || Ka==16 || Ka==17) && wda==0)                                      /* A becomes WD */
       {wda=1; ta_end=t; Mwda_0=Ma;}
     if(((Kb>=10 && Kb<=12) || Kb==16 || Kb==17) && wdb==0)                                      /* B becomes WD */
       {wdb=1; tb_end=t; Mwdb_0=Mb;}



     Kace=Kbce=0;
     ace=Mace=Mbce=Race=Rbce=Mcace=Mcbce=Mheace=Mhebce=Mcoace=Mcobce=0.0;

     if(inbin==1) {                                                                 /* perform only for binaries */

       J0=Ma*Mb*sqrt((GGG*a)/(Ma+Mb))+Iaold*wa+Ibold*wb;
       Iatmp=Iaold;
       Ibtmp=Ibold;
       // if(PP==1 && dMmta>acc)
       //   printf("RR MT1: dM_trans=%g [Msun/yr], dM_acc=%g [Msun/yr], mt type: %d\n",dMmta*1.0e-06,dMmtb*1.0e-06,mttype);
       // if(PP==1 && dMmtb>acc)
       //   printf("RR MT2: dM_trans=%g [Msun/yr], dM_acc=%g [Msun/yr], mt type: %d\n",dMmtb*1.0e-06,dMmta*1.0e-06,mttype);
       mark1=orb_change(t-dt,t,&a,&e,&wa,&wb,tvira,tvirb,Ma,Mb,M0a,M0b,Mzamsa,Mzamsb,dMwinda,dMwindb,Mca,Mcb,Ra,Rb,Raold,
                        Rbold,La,Lb,Ka,Kb,Kaold,Kbold,mt,ce,mttype,&Iaold,&Ibold,&KTa,&KTb,dMmta,dMmtb,Mdisa,Mdisb,&darwin);


       if(mark1==0 || wa<(-acc) || wb<(-acc) || e>1.0 || e<-acc) {         /* system not calculated properly: skip it */
         stop2=1;
         badorb++;
         typce=10;                                                      /* to protect writing this system as a merger */
         break;
       }

       dtorb=0.0;
       J1=Ma*Mb*sqrt((GGG*a)/(Ma+Mb))+Iaold*wa+Ibold*wb;
       if(mark1!=2 && maxred<4 && (t-dt)>acc && ce==0 && mt==0 && sn==0 && Ka==Kaold && Kb==Kbold && fabs(J1-J0)>0.01*fabs(J0)) {
         dtorb=dt/5.0;                                                    /* decrease timestep: no more than 5^4=625 times */
         t-=dt;                                                           /* return to starting parameter values */
         Iaold=Iatmp;
         Ibold=Ibtmp;
         a=aold;
         e=eold;
         wa=waold;
         wb=wbold;
         copySin(sa1,sa2,sa3,sa4,sa5,sa6,sa7,sa8,sa9,sa10,sa11,sa12,sa13,sa14,sa15,sa16,sa17,sa18,sa19,sa20,
                 &Mzamsa,&M0a,&Ma,&Ka,&tbeg,&tvira,&tenda,&La,&Ra,&Mca,&Mhea,&Mcoa,&flaga,&dta,&Mprea,&Kpa,&tstarta,&fraca,&dMmta,&ecssna);
         copySin(sb1,sb2,sb3,sb4,sb5,sb6,sb7,sb8,sb9,sb10,sb11,sb12,sb13,sb14,sb15,sb16,sb17,sb18,sb19,sb20,
                 &Mzamsb,&M0b,&Mb,&Kb,&tbeg,&tvirb,&tendb,&Lb,&Rb,&Mcb,&Mheb,&Mcob,&flagb,&dtb,&Mpreb,&Kpb,&tstartb,&fracb,&dMmtb,&ecssnb);
         undo_SNcounters();
         maxred++;
         dper=a*(1.0-e);                                                      /* separation of stars at peryastron */
         Rla=roche(Mb/Ma,dper);                                       /* minimal roche lobe sizes -- at peryastron */
         Rlb=roche(Ma/Mb,dper);
         continue;
       }
       if(maxred==5) {                                                          /* problem with keeping del(Jtot) small */
         fprintf(fp0,"error: problem with keeping del(Jtot) small: %g, %d\n",fabs(J1-J0)/fabs(J0),iidd_old);
         fflush(fp0);
       }
       maxred=0;
       dper=a*(1.0-e);                                                      /* separation of stars at peryastron */
       Rla=roche(Mb/Ma,dper);                                       /* minimal roche lobe sizes -- at peryastron */
       Rlb=roche(Ma/Mb,dper);

       if(Ra>=Rla || Rb>=Rlb) {                         /* if any contact: circularization and synchronization assumed */
         if(e>0.01 && Ra>=Rla && Rb>=Rlb) mtflag1ab+=1;
         else if(e>0.01 && Ra>=Rla) mtflag1a+=1;
         else if(e>0.01 && Rb>=Rlb) mtflag1b+=1;
         a=a*(1.0-e);                                                    /* circularization at periastron assumed */
         e=0.0;                                                          /* periastron assumed: so at any MT func. e=0 */
         worb=sqrt(GGG*(Ma+Mb))*pow(a,-1.5);
         Rla=roche(Mb/Ma,a);
         Rlb=roche(Ma/Mb,a);
         if(Ka<10 && fabs(wa-worb)>0.05*worb) mtflag2+=1;
         if(Kb<10 && fabs(wb-worb)>0.05*worb) mtflag3+=1;
         if(Ka<10) wa=worb;
         if(Kb<10) wb=worb;                                  /* synchronized to mean orbital angular velocity [Myr^-1] */
       }


       if(mark1==3 || ((Ra>Rla || Rb>Rlb) && Ma*Mb*sqrt((GGG*a)/(Ma+Mb))<=3.0*(Iaold*wa+Ibold*wb))) { /* MT sys. Darwin inst. */
         if((Ra>Rla && Ka>=2 && Ka<=9 && Ka!=7) || (Rb>Rlb && Kb>=2 && Kb<=9 && Kb!=7))               /* CE evolution */
           doce=1;
         else {
           if(SN1A==1 && (Ka==10 || Ka==11 || Ka==12 || Ka==16 || Ka==17) && (Kb==10 || Kb==11 || Kb==12 || Kb==16 || Kb==17)) {
             typwd=4;
             fprintf(fp180,"%d  %f %f  %d %d  %f %f %f %f  %f %f  -1 -1 -1 -1 -1 -1   %f %f %f %f %f %f  %f %f %f %f    -1 -1 -1 -1 -1 -1 -1  %d %d  %d %s\n",
                    typwd,t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,Ra,Rb,
                    M0a,M0b,Mzamsa,Mzamsb,a0,e0,Mwda_0,Mwdb_0,ta_end,tb_end,
                    idum_run,iidd_old,nevroute,evroute);
             fflush(fp180);
           }                                                                                           /* merger assumed */
           stop2=1;
           if(TranCE==1) {
             typce=2;                                           /* components parameters passed to record in Merger file */
             fprintf(fp900,"%d  %f %f  %d %d  %f %f %f %f  %f %f  %f %f %f %f %f %f   %f %f %f %f %f %f  %d %d  %d %s\n",
                     5,t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,Ra,Rb,Mhea,Mcoa,Mheb,Mcob,Mca,Mcb,
                     M0a,M0b,Mzamsa,Mzamsb,a0,e0,idum_run,iidd_old,nevroute,evroute);
             fflush(fp900);
             break;
           }
         }
       }

       if(dMmta>acc && mttype==4) Ratrue=1.01*Rla;    /* at thermal MT, true radius of donor is equal to roche lobe radius */
       else Ratrue=Ra;
       if(dMmtb>acc && mttype==4) Rbtrue=1.01*Rlb;
       else Rbtrue=Rb;
       wbreaka=sqrt(GGG*Ma/pow(Ratrue,3.0));
       wbreakb=sqrt(GGG*Mb/pow(Rbtrue,3.0));
       if(wa>wbreaka) {                          /* one component spins over breakup velocity, flag */
         flagbra+=1;                             /* slow it down to berakup velocity: star would shed little of */
         wa=wbreaka;                             /* mass with high J loss and would slow down */
       }
       if(wb>wbreakb) {
       	 flagbrb+=1;
       	 wb=wbreakb;
       }

       dMmta_old=dMmta;
       dMmtb_old=dMmtb;
       dMmta=dMmtb=0.0;                                       /* if needed they are filled below, otherwise zeros */
       dMatmorlofa=dMatmorlofb=0.0;

       mt=0;
       if(Ra<Rla || mttype!=4) {
         dMtha=Ttha=0.0;
         deca=0;
       }
       if(Rb<Rlb || mttype!=4) {
       	 dMthb=Tthb=0.0;
       	 decb=0;
       }
       mttype_old=mttype;

       mttype=0;
       if((Ka==13 || Ka==14) && Kb!=13 && Kb!=14)
         dMtran=dMtranf(Mb,Ma,a,e,Kb,Ka);
       else if((Kb==13 || Kb==14) && Ka!=13 && Ka!=14)
         dMtran=dMtranf(Ma,Mb,a,e,Ka,Kb);
       else
         dMtran=0.0;

       if(Ra>=Rla && doce==0 && (!((Kb==13 || Kb==14) && Kbold<=9))) {               /* skip: CE and SN steps */
         if(Ka==13 || Ka==14) {                    /* STOP: donor NS or BH, companion must be remnant as well */
           stop1=1;                                   /* NS in tight WD-NS binary, brought by hand to contact */
           typce=10;
           break;                     /* NS lower mass, same radius as massive WD, gets to contact first (GR) */
         }
         copySin(Mzamsa,M0a,Ma,Ka,tbeg,tvira,tenda,La,Ra,Mca,Mhea,Mcoa,flaga,dta,Mprea,Kpa,tstarta,fraca,dMmta,ecssna,
                 &input1[0],&input1[1],&input1[2],&input2[0],&input1[3],&input1[4],&input1[5],&input1[6],&input1[7],
                 &input1[8],&input1[9],&input1[10],&input2[1],&input1[11],&input1[12],&input2[2],&input1[13],
                 &input1[14],&input1[15],&input2[3]);
         dMmta=dMmtf(Ma,Mb,Ra,Rb,La,Lb,wa,Iaold,KTa,wb,Ibold,KTb,Ka,Kb,a,e,t,dt,input1,input2,dMmta_old,dMmtb_old,
                     &stop2,&mttype,mttypelast,&doce,&merger,&dMtha,&Ttha,&deca,aspinb,Mzamsa,Mzamsb);
         dMmtb=dMgainf(dMmta,dMtran,Mb,Ka,Kb,Rb,&Mout1b,&Mout2b,&mark_outb,&doce,&merger,&mark777b,&Macc0b,aspinb);  /* acc rate onto B */
         tlast=t; Kdonlast=Ka; Kacclast=Kb; mttypelast=mttype; Mdonlast=Ma; Macclast=Mb; alast=a;
         if(dMmta<1.0e-20)                      /* if dMmt is numerical zero, we protect from division by zero */
           dtmt=1.0e+10;
         else
           dtmt=0.01*Ma/dMmta;		      /* next time step such that max. 1% of donor mass is transffered */
         dtmt=min(dtmt,dta);                       /* also the stars natural timescales are taken into account */
         dtmt=min(dtmt,dtb);
         dtmt=min(dtmt,100.0);                                                    /* step smaller than 100 Myr */
         if(mttype_old==4 && mttype_old!=mttype && Ra>1.99*Rla) {      /* aftre thermal MT, star still too big */
           if(Ka==0 || Ka==1 || Ka==7 || Ka>9)                       /* do merger for MS, HeMS, remannt donors */
             merger=1;
           else                                                                 /* do CE for giant-like donors */
             doce=1;
         }
         if(merger==1) {                                                   /* dyn. istab., merger of some type */
           if(((Ka>=10 && Ka<=12) || Ka==16 || Ka==17) && ((Kb>=10 && Kb<=14) || Kb==16 || Kb==17)) {
             if(SN1A==1 && (Ka==10 || Ka==11 || Ka==12 || Ka==16 || Ka==17) && (Kb==10 || Kb==11 || Kb==12 || Kb==16 || Kb==17)) {
               typwd=5;
               fprintf(fp180,"%d  %f %f  %d %d  %f %f %f %f  %f %f  -1 -1 -1 -1 -1 -1   %f %f %f %f %f %f  %f %f %f %f    -1 -1 -1 -1 -1 -1 -1  %d %d  %d %s\n",
                      typwd,t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,Ra,Rb,
                      M0a,M0b,Mzamsa,Mzamsb,a0,e0,Mwda_0,Mwdb_0,ta_end,tb_end,
                      idum_run,iidd_old,nevroute,evroute);
               fflush(fp180);
             }
             stop1=1;                                       /* STOP: WD+WD/NS/BH merger, compact object merger */
           }
           else                                                                                /* other merger */
             stop2=1;
           if(TranCE==1) {
             fprintf(fp900,"%d  %f %f  %d %d  %f %f %f %f  %f %f  %f %f %f %f %f %f   %f %f %f %f %f %f  %d %d  %d %s ",
                     6,t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,Ra,Rb,Mhea,Mcoa,Mheb,Mcob,Mca,Mcb,
                     M0a,M0b,Mzamsa,Mzamsb,a0,e0,idum_run,iidd_old,nevroute,evroute);
             fprintf(fp900," 0 ");
             fprintf(fp900,"\n");
             fflush(fp900);
             typce=3;                                 /* components parameters passed to record in Merger file */
             break;
           }
         }
         if(dMmta>acc && doce==0) {                                  /* detailed MT calculation in coming step */
           mt=1;
           if(mt1H==0) {
             evroute_add("MT1(%d-%d) ",Ka,Kb);
             mt1H=1;
           }
         }
       }
       else if(Rb>=Rlb && doce==0 && (!((Ka==13 || Ka==14) && Kaold<=9))) {          /* skip: CE and SN steps */
         if(Kb==13 || Kb==14) {                    /* STOP: donor NS or BH, companion must be remnant as well */
           stop1=1;
           typce=10;
           break;
         }
         copySin(Mzamsb,M0b,Mb,Kb,tbeg,tvirb,tendb,Lb,Rb,Mcb,Mheb,Mcob,flagb,dtb,Mpreb,Kpb,tstartb,fracb,dMmtb,ecssnb,
                 &input1[0],&input1[1],&input1[2],&input2[0],&input1[3],&input1[4],&input1[5],&input1[6],&input1[7],
       	         &input1[8],&input1[9],&input1[10],&input2[1],&input1[11],&input1[12],&input2[2],&input1[13],
                 &input1[14],&input1[15],&input2[3]);
         dMmtb=dMmtf(Mb,Ma,Rb,Ra,Lb,La,wb,Ibold,KTb,wa,Iaold,KTa,Kb,Ka,a,e,t,dt,input1,input2,dMmtb_old,dMmta_old,
                     &stop2,&mttype,mttypelast,&doce,&merger,&dMthb,&Tthb,&decb,aspina,Mzamsb,Mzamsa);
         dMmta=dMgainf(dMmtb,dMtran,Ma,Kb,Ka,Ra,&Mout1a,&Mout2a,&mark_outa,&doce,&merger,&mark777a,&Macc0a,aspina);
         tlast=t; Kdonlast=Kb; Kacclast=Ka; mttypelast=mttype; Mdonlast=Mb; Macclast=Ma; alast=a;
         if(dMmtb<1.0e-20)                      /* if dMmt is numerical zero, we protect from division by zero */
           dtmt=1.0e+10;
         else
           dtmt=0.01*Mb/dMmtb;		      /* next time step such that only 1% of donor mass is transffered */
         dtmt=min(dtmt,dta);                       /* also the stars natural timescales are taken into account */
         dtmt=min(dtmt,dtb);
         dtmt=min(dtmt,100.0);                                                    /* step smaller than 100 Myr */
         if(mttype_old==4 && mttype_old!=mttype && Rb>1.99*Rlb) {      /* aftre thermal MT, star still too big */
           if(Kb==0 || Kb==1 || Kb==7 || Kb>9)                       /* do merger for MS, HeMS, remannt donors */
             merger=1;
           else                                                                 /* do CE for giant-like donors */
             doce=1;
         }
	 if(merger==1) {                                                   /* dyn. istab., merger of some type */
           if(((Kb>=10 && Kb<=12) || Kb==16 || Kb==17) && ((Ka>=10 && Ka<=14) || Ka==16 || Ka==17)) {
             if(SN1A==1 && (Ka==10 || Ka==11 || Ka==12 || Ka==16 || Ka==17) && (Kb==10 || Kb==11 || Kb==12 || Kb==16 || Kb==17)) {
               typwd=5;
               fprintf(fp180,"%d  %f %f  %d %d  %f %f %f %f  %f %f  -1 -1 -1 -1 -1 -1   %f %f %f %f %f %f  %f %f %f %f    -1 -1 -1 -1 -1 -1 -1  %d %d  %d %s\n",
                      typwd,t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,Ra,Rb,
                      M0a,M0b,Mzamsa,Mzamsb,a0,e0,Mwda_0,Mwdb_0,ta_end,tb_end,
                      idum_run,iidd_old,nevroute,evroute);
               fprintf(fp180,"\n");
               fflush(fp180);
             }
             stop1=1;                                       /* STOP: WD+WD/NS/BH merger, compact object merger */
           }
           else if(Kb==7 && (Ka==13 || Ka==14))                                              /* He star merger */
             stop1=1;
           else                                                                             /* MS+BH/NS merger */
             stop2=1;
           if(TranCE==1) {
             typce=3;                                   /* components parameters passed to record in Merger file */
             fprintf(fp900,"%d  %f %f  %d %d  %f %f %f %f  %f %f  %f %f %f %f %f %f   %f %f %f %f %f %f  %d %d  %d %s ",
                     7,t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,Ra,Rb,Mhea,Mcoa,Mheb,Mcob,Mca,Mcb,
                     M0a,M0b,Mzamsa,Mzamsb,a0,e0,idum_run,iidd_old,nevroute,evroute);
             fprintf(fp900," 0 ");
             fprintf(fp900,"\n");
             fflush(fp900);
             break;
           }
         }
         if(dMmtb>acc && doce==0) {                                    /* detailed MT calculation in coming step */
           mt=1;
           if(mt2H==0) {
             evroute_add("MT2(%d-%d) ",Ka,Kb);
             mt2H=1;
           }
         }
       }
       else if(AtmoRLOF==1 && Ka<=9 && Ra>=0.8*Rla && (Ra/Rla)>(Rb/Rlb) && (!((Kb==13 || Kb==14) && Kbold<=9))) {   /* atmospheric RLOF, skip: SN steps */
         dMatmorlofa=dMamtf(Ma,Mb,La,Ra,Rla,a,Ka);
         if(Kb==10 || Kb==11 || Kb==12 || Kb==17)       /* change accretion into accumulations for WD accretors */
           dMatmorlofb=dMgainf(dMatmorlofa,dMtran,Mb,Ka,Kb,Rb,&Mout1b,&Mout2b,&mark_outb,&doce,&merger,&mark777b,&Macc0b,aspinb);  /* acc rate onto B */
         else                                           /* for all other stars accrete ALL, and not only Fa*ALL as in case of regular RLOF */
           dMatmorlofb=-dMatmorlofa;
         mt=1;
         mttype=6;
         tlast=t; Kdonlast=Ka; Kacclast=Kb; mttypelast=mttype; Mdonlast=Ma; Macclast=Mb; alast=a;
         dtmt=0.01*Ma/dMatmorlofa;		      /* next time step such that max. 1% of donor mass is transffered */
         dtmt=min(dtmt,dta);                          /* also the stars natural timescales are taken into account */
         dtmt=min(dtmt,dtb);
         dtmt=min(dtmt,100.0);                                                    /* step smaller than 100 Myr */
       }
       else if(AtmoRLOF==1 && Kb<=9 && Rb>=0.8*Rlb && (!((Ka==13 || Ka==14) && Kaold<=9))) {     /* atmospheric RLOF, skip: SN steps */
         dMatmorlofb=dMamtf(Mb,Ma,Lb,Rb,Rlb,a,Kb);
         if(Ka==10 || Ka==11 || Ka==12 || Ka==17)
           dMatmorlofa=dMgainf(dMatmorlofb,dMtran,Ma,Kb,Ka,Ra,&Mout1a,&Mout2a,&mark_outa,&doce,&merger,&mark777a,&Macc0a,aspina);
         else                                           /* for all other stars accrete ALL, and not only Fa*ALL as in case of regular RLOF */
           dMatmorlofa=-dMatmorlofb;
         mt=1;
         mttype=6;
         tlast=t; Kdonlast=Kb; Kacclast=Ka; mttypelast=mttype; Mdonlast=Mb; Macclast=Ma; alast=a;
         dtmt=0.01*Mb/dMatmorlofb;		        /* next time step such that only 1% of donor mass is transffered */
         dtmt=min(dtmt,dta);                            /* also the stars natural timescales are taken into account */
         dtmt=min(dtmt,dtb);
         dtmt=min(dtmt,100.0);                                                    /* step smaller than 100 Myr */
       }
                                                         /* output for compact object binaries: NS/BH-NS/BH */
       // if(comdone==0 && COMPACT==1 && (Ka==13 || Ka==14) && (Kb==13 || Kb==14)) {
       //   comdone=1;
       //
       //   fprintf(fp200,"%f  %d %d  %f %f  %f %f %f %f  %f %f  %d\n",t,Ka,Kb,Ma,Mb,RmaxA,RmaxB,a,e,ta_end,tb_end,inbin);
       //   fprintf(fp200,"%f %f %f %f  %f\n",a0,e0,Mzamsa,Mzamsb,ZZ);
       //   fprintf(fp200,"\n");
       //   fflush(fp200);
       // }

                                                                            /* output for BH-WD and BH-NS systems at formation */
       if(inbin==1 && BHNSWD==1 && SPCout==0 && (((Kaold==13 || Kaold==14) && (Kb==10 || Kb==11 || Kb==12 || Kb==16 || Kb==17)) || ((Kbold==13 || Kbold==14) && (Ka==10 || Ka==11 || Ka==12 || Ka==16 || Ka==17)))) {
         SPCout=1;
         fprintf(fp170,"%d %d %.3f %.3f   %f %f   %.3f %.3f %.3f   %d %d   ",Ka,Kb,Ma,Mb,a,e,ta_end,tb_end,t,idum_run,iidd_old);
         fprintf(fp170,"%d %s\n",nevroute,evroute);
         fflush(fp170);
       }

       Kawd=Kbwd=0;
       if(Ka==10 || Ka==11 || Ka==12 || Ka==16 || Ka==17) Kawd=1;                    /* CV output/classification */
       if(Kb==10 || Kb==11 || Kb==12 || Kb==16 || Kb==17) Kbwd=1;
       Kawdold=Kbwdold=0;
       if(Kaold==10 || Kaold==11 || Kaold==12 || Kaold==16 || Kaold==17) Kawdold=1;
       if(Kbold==10 || Kbold==11 || Kbold==12 || Kbold==16 || Kbold==17) Kbwdold=1;

       if(CVout==1 && Tcv>(t-dt+Tst) && Tcv<=(t+Tst) && mt==1 && ((Kawd==1 && dMmtb>0.0) || (Kbwd==1 && dMmta>0.0))) {
         if(Ra>=Rla && Rb<Rlb) cvdon=1;         /* A donor, B accretor */
         else if (Rb>=Rlb && Ra<Rla) cvdon=2;   /* B donor, A accretor */
         else fprintf(fp0,"error: unknown type of donor in CV classification (%d %d)\n",idum_run,iidd_old);

         if(cvdon==1)                                         /* A donor, B accretor; calculate Lxcv */
           Lxcv=Xcvf(Ka,Kb,Mb,Rb,dMmta,&dMacc,&cvtype);
         else                                                 /* B donor, A accretor; calculate Lxcv */
           Lxcv=Xcvf(Kb,Ka,Ma,Ra,dMmtb,&dMacc,&cvtype);
         per=2.0*Pi*a*sqrt(a/(G*(Ma+Mb)));                                               /* binary period [days] */

         fprintf(fp151,"%.1f %.1f  %d %d  %d %d  %.3f %.3f %.3f %.3f %f  %.1f %.1f  %d %g  %f %.3f %f %.3f  %.2f %.3f %.2f %.2f  %g %g  %d %d  %d %s\n",
                        t+Tst,Tst,cvdon,cvtype,Ka,Kb,Ma,Mb,a,e,per,ta_end+Tst,tb_end+Tst,mt,Lxcv,Ra,Rla,Rb,Rlb,a0,e0,Mzamsa,Mzamsb,dMmta,dMmtb,
                        idum_run,iidd_old,nevroute,evroute);
         fflush(fp151);
       }


       if(Kawdold==1 && Kawd==1 && Kbwdold==1 && Kbwd==1 && WDWD==1 && Merger==0) {   /* output for WD-WD at the formation */
         per=2.0*Pi*a*sqrt(a/(G*(Ma+Mb)));                                               /* binary period [days] */
         Tmr=tmerge(Ma,Mb,a,e);

         if(ta_end<=tb_end)                                                     /* star A first turned to remnant */
           tend=tb_end;
         else                                                                   /* star B first turned to remnant */
           tend=ta_end;

         fprintf(fp600,"%f %f %f  %d %d %f %f %f %f %f %f %f  %f %f %g  %f %f %f %f %d %d  %d %s\n",
                 tend,Ma,Mb,Ka,Kb,a,e,per,ta_end,tb_end,Mwda_0,Mwdb_0,Ra,Rb,Tmr,
                 a0,e0,Mzamsa,Mzamsb,idum_run,iidd_old,nevroute,evroute);
         fflush(fp600);
         break;
       }
                        /* if doing WDWD output: do not continue evolution once NS or BH or single SN Ia is formed */
                        /* do continue evolution if Merger output required, typce not required here before break */
       if(WDWD==1 && Merger==0 && ((Kaold==13 && Ka==13) || (Kaold==14 && Ka==14) || (Kaold==15 && Ka==15)))
         break;
       if(WDWD==1 && Merger==0 && ((Kbold==13 && Kb==13) || (Kbold==14 && Kb==14) || (Kbold==15 && Kb==15)))
         break;



       if(XEdd==2 && (Kb==13 || Kb==14) && (Kbold==13 || Kbold==14) && doce==0 && ce==0)       /* B: NS/BH potential accretor */
         Xbin1(t,dt,Ma,Mb,Ka,Kb,Ra,dMmta,dMwinda,aspinb,a,e,La);
       else if(XEdd==2 && (Ka==13 || Ka==14) && (Kaold==13 || Kaold==14) && doce==0 && ce==0)  /* A: NS/BH potential accretor */
         Xbin1(t,dt,Mb,Ma,Kb,Ka,Rb,dMmtb,dMwindb,aspina,a,e,Lb);


                                                                                /* do for any binary with NS or BH!!! */
       if(XEdd==1 && (Kb==13 || Kb==14) && (Kbold==13 || Kbold==14) && doce==0 && ce==0) {          /* B: NS/BH potential accretor */
         Xbin(&Lxmt,dMtran,Ma,dMwinda,Mca,dMmta,La,Ra,t,dt,Mb,Rb,a,e,Ka,Kb,Vsm,tb_end,&mark10,&mark11,mttype,
              flagbra,flagbrb,mtflag1a,mtflag1b,mtflag1ab,Maio,Mbio,dMcea,dMceb,Mzamsa,Mzamsb,a0,e0,ecssna,ecssnb,aspinb);
       }
       else if(XEdd==1 && (Ka==13 || Ka==14) && (Kaold==13 || Kaold==14) && doce==0 && ce==0) {     /* A: NS/BH potential accretor */
         Xbin(&Lxmt,dMtran,Mb,dMwindb,Mcb,dMmtb,Lb,Rb,t,dt,Ma,Ra,a,e,Kb,Ka,Vsm,ta_end,&mark10,&mark11,mttype,
              flagbrb,flagbra,mtflag1a,mtflag1b,mtflag1ab,Mbio,Maio,dMceb,dMcea,Mzamsb,Mzamsa,a0,e0,ecssnb,ecssna,aspina);
       }

       if(SSout==1 && (t+Tst)>=(TIMEOUT-DTOUT) && (t+Tst)<=(TIMEOUT+DTOUT) && doce==0) {       /* output symbiotics */
         symbiotic(Ma,dMwinda,Mca,La,Ra,Rb,t,Mb,a,e,Ka,Kb,dMmta,dMmtb,dMtran,mttype,Mzamsa,Mzamsb,a0,e0);
         symbiotic(Mb,dMwindb,Mcb,Lb,Rb,Ra,t,Ma,a,e,Kb,Ka,dMmtb,dMmta,dMtran,mttype,Mzamsb,Mzamsa,a0,e0);
       }


       if(BHCAT2==1 && ((Ka==14 && Kaold>=10) || (Kb==14 && Kbold>=10)) && doce==0) {  /* BH sys. detailed output */
         per=2.0*Pi*a*sqrt(a/(G*(Ma+Mb)));                                               /* binary period [days] */
         if(mt==1) dMbh=max(dMmta,dMmtb);                                           /* mass transfer is positive */
         else dMbh=0.0;
         fprintf(fp160,"%.1f %.1f %d %d  %.3f %.3f %.3f %.3f %f  %d %g %g  %f %.3f %f %.3f  %.2f %.2f %.2f %.2f  %.1f %.1f %.1f  %d %d\n",
                 t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,per,mt,Lxmt,dMbh,Ra,Rla,Rb,Rlb,Mzamsa,Mzamsb,a0,e0,Vsm[0],Vsm[1],Vsm[2],idum_run,iidd_old);
         fflush(fp160);
       }


       tform=max(Tst+ta_end,Tst+tb_end);                                /* output for WD/NS/BH-WD/NS/BH binaries */
       if(WDWDout==1 && WDout==0 && tform<Twdwd && (t+Tst)>=Twdwd && Ka>=10 && Ka<=17 && Ka!=15 && Kb>=10 && Kb<=17 && Kb!=15) {
         WDout=1;
         per=2.0*Pi*a*sqrt(a/(G*(Ma+Mb)));                                               /* binary period [days] */
         if(dMmta>acc && doce==0)                                                        /* A donor, B accretor */
           if(Kb==13 || Kb==14)
             Lxwd=Lxmt;
           else
             Xwdwd(&Lxwd,dMmta,Mb,Rb,Ka,Kb);
         else if(dMmtb>acc && doce==0)
           if(Ka==13 || Ka==14)
             Lxwd=Lxmt;
           else
             Xwdwd(&Lxwd,dMmtb,Ma,Ra,Kb,Ka);
         else
           Lxwd=0.0;
         if(Ka==13 || Ka==14)
           Lwd_a=0.0;
         else
           Lwd_a=Lwdf(Ma,t-ta_end,Ka);                                                /* Luminosity of WD A [Lsun] */
         if(Kb==13 || Kb==14)
           Lwd_b=0.0;
         else
           Lwd_b=Lwdf(Mb,t-tb_end,Kb);
         fprintf(fp150,"%.1f %.1f %d %d  %.3f %.3f %.3f %.3f %f  %.1f %.1f  %d %g  %f %.3f %f %.3f  %.2f %.3f %.2f %.2f  %g %g  %d %d  %d %s\n",
                t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,per,ta_end+Tst,tb_end+Tst,mt,Lxwd,Ra,Rla,Rb,Rlb,a0,e0,Mzamsa,Mzamsb,Lwd_a,Lwd_b,idum_run,iidd_old,nevroute,evroute);
         fflush(fp150);
       }
     }



     if(BHCAT1==1 && TIMEOUT1>(t-dt+Tst) && TIMEOUT1<=(t+Tst) && ((Ka==14 && Kaold>=10) || (Kb==14 && Kbold>=10))) {
       if(inbin==1) {
         if(Ka==14) dMedd=dMeddf(Ka,Ra,Kb);
         else dMedd=dMeddf(Kb,Rb,Ka);
         per=2.0*Pi*a*sqrt(a/(G*(Ma+Mb)));                                                /* binary period [days] */
         per*=(24.0);                                                                     /* [days->hrs] */
         fprintf(fp141,"%.3f %.1f %g %d %d %.3f %.3f %.1e %.3f %.1f %.1f %.1f %.3f %d %d %d %.3f %.1e %.1e %.1e %.1e %.1e %.1f %.1f %.1f %f %.1f %f %.1f %.3f %.3f\n",
                        t+Tst,Tst,ZZ,Ka,Kb,Ma,Mb,a,e,Mzamsa,Mzamsb,a0,e0,inbin,iidd_old,
                        mttype,per,dMmta,dMmtb,dMedd,dMtran,Lxmt,Vsm[0],Vsm[1],Vsm[2],Ra,Rla,Rb,Rlb,ta_end+Tst,tb_end+Tst);
       }
       else
         fprintf(fp141,"%.3f %.1f %g %d %d %.3f %.3f %.1e %.3f %.1f %.1f %.1f %.3f %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.3f %.3f %.3f\n",
                        t+Tst,Tst,ZZ,Ka,Kb,Ma,Mb,a,e,Mzamsa,Mzamsb,a0,e0,inbin,iidd_old,
                        Vsa[0],Vsa[1],Vsa[2],Vsb[0],Vsb[1],Vsb[2],ta_end+Tst,tb_end+Tst,tdis+Tst);
       fflush(fp141);
     }

     if(BHCAT1==1 && TIMEOUT2>(t-dt+Tst) && TIMEOUT2<=(t+Tst) && ((Ka==14 && Kaold>=10) || (Kb==14 && Kbold>=10))) {
       if(inbin==1) {
         if(Ka==14) dMedd=dMeddf(Ka,Ra,Kb);
         else dMedd=dMeddf(Kb,Rb,Ka);
         per=2.0*Pi*a*sqrt(a/(G*(Ma+Mb)));                                                /* binary period [days] */
         per*=(24.0);                                                                     /* [days->hrs] */
         fprintf(fp142,"%.3f %.1f %g %d %d %.3f %.3f %.1e %.3f %.1f %.1f %.1f %.3f %d %d %d %.3f %.1e %.1e %.1e %.1e %.1e %.1f %.1f %.1f %f %.1f %f %.1f %.3f %.3f\n",
                        t+Tst,Tst,ZZ,Ka,Kb,Ma,Mb,a,e,Mzamsa,Mzamsb,a0,e0,inbin,iidd_old,
                        mttype,per,dMmta,dMmtb,dMedd,dMtran,Lxmt,Vsm[0],Vsm[1],Vsm[2],Ra,Rla,Rb,Rlb,ta_end+Tst,tb_end+Tst);
       }
       else
         fprintf(fp142,"%.3f %.1f %g %d %d %.3f %.3f %.1e %.3f %.1f %.1f %.1f %.3f %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.3f %.3f %.3f\n",
                        t+Tst,Tst,ZZ,Ka,Kb,Ma,Mb,a,e,Mzamsa,Mzamsb,a0,e0,inbin,iidd_old,
                        Vsa[0],Vsa[1],Vsa[2],Vsb[0],Vsb[1],Vsb[2],ta_end+Tst,tb_end+Tst,tdis+Tst);
       fflush(fp142);
     }

     if(BHCAT1==1 && TIMEOUT3>(t-dt+Tst) && TIMEOUT3<=(t+Tst) && ((Ka==14 && Kaold>=10) || (Kb==14 && Kbold>=10))) {
       if(inbin==1) {
         if(Ka==14) dMedd=dMeddf(Ka,Ra,Kb);
         else dMedd=dMeddf(Kb,Rb,Ka);
         per=2.0*Pi*a*sqrt(a/(G*(Ma+Mb)));                                                /* binary period [days] */
         per*=(24.0);                                                                     /* [days->hrs] */
         fprintf(fp143,"%.3f %.1f %g %d %d %.3f %.3f %.1e %.3f %.1f %.1f %.1f %.3f %d %d %d %.3f %.1e %.1e %.1e %.1e %.1e %.1f %.1f %.1f %f %.1f %f %.1f %.3f %.3f\n",
                        t+Tst,Tst,ZZ,Ka,Kb,Ma,Mb,a,e,Mzamsa,Mzamsb,a0,e0,inbin,iidd_old,
                        mttype,per,dMmta,dMmtb,dMedd,dMtran,Lxmt,Vsm[0],Vsm[1],Vsm[2],Ra,Rla,Rb,Rlb,ta_end+Tst,tb_end+Tst);
       }
       else
         fprintf(fp143,"%.3f %.1f %g %d %d %.3f %.3f %.1e %.3f %.1f %.1f %.1f %.3f %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.3f %.3f %.3f\n",
                        t+Tst,Tst,ZZ,Ka,Kb,Ma,Mb,a,e,Mzamsa,Mzamsb,a0,e0,inbin,iidd_old,
                        Vsa[0],Vsa[1],Vsa[2],Vsb[0],Vsb[1],Vsb[2],ta_end+Tst,tb_end+Tst,tdis+Tst);
       fflush(fp143);
     }

     if(BHCAT1==1 && TIMEOUT4>(t-dt+Tst) && TIMEOUT4<=(t+Tst) && ((Ka==14 && Kaold>=10) || (Kb==14 && Kbold>=10))) {
       if(inbin==1) {
         if(Ka==14) dMedd=dMeddf(Ka,Ra,Kb);
         else dMedd=dMeddf(Kb,Rb,Ka);
         per=2.0*Pi*a*sqrt(a/(G*(Ma+Mb)));                                                /* binary period [days] */
         per*=(24.0);                                                                     /* [days->hrs] */
         fprintf(fp144,"%.3f %.1f %g %d %d %.3f %.3f %.1e %.3f %.1f %.1f %.1f %.3f %d %d %d %.3f %.1e %.1e %.1e %.1e %.1e %.1f %.1f %.1f %f %.1f %f %.1f %.3f %.3f\n",
                        t+Tst,Tst,ZZ,Ka,Kb,Ma,Mb,a,e,Mzamsa,Mzamsb,a0,e0,inbin,iidd_old,
                        mttype,per,dMmta,dMmtb,dMedd,dMtran,Lxmt,Vsm[0],Vsm[1],Vsm[2],Ra,Rla,Rb,Rlb,ta_end+Tst,tb_end+Tst);
       }
       else
         fprintf(fp144,"%.3f %.1f %g %d %d %.3f %.3f %.1e %.3f %.1f %.1f %.1f %.3f %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.3f %.3f %.3f\n",
                        t+Tst,Tst,ZZ,Ka,Kb,Ma,Mb,a,e,Mzamsa,Mzamsb,a0,e0,inbin,iidd_old,
                        Vsa[0],Vsa[1],Vsa[2],Vsb[0],Vsb[1],Vsb[2],ta_end+Tst,tb_end+Tst,tdis+Tst);
       fflush(fp144);
     }

     if(BHCAT1==1 && TIMEOUT5>(t-dt+Tst) && TIMEOUT5<=(t+Tst) && ((Ka==14 && Kaold>=10) || (Kb==14 && Kbold>=10))) {
       if(inbin==1) {
         if(Ka==14) dMedd=dMeddf(Ka,Ra,Kb);
         else dMedd=dMeddf(Kb,Rb,Ka);
         per=2.0*Pi*a*sqrt(a/(G*(Ma+Mb)));                                                /* binary period [days] */
         per*=(24.0);                                                                     /* [days->hrs] */
         fprintf(fp145,"%.3f %.1f %g %d %d %.3f %.3f %.1e %.3f %.1f %.1f %.1f %.3f %d %d %d %.3f %.1e %.1e %.1e %.1e %.1e %.1f %.1f %.1f %f %.1f %f %.1f %.3f %.3f\n",
                        t+Tst,Tst,ZZ,Ka,Kb,Ma,Mb,a,e,Mzamsa,Mzamsb,a0,e0,inbin,iidd_old,
                        mttype,per,dMmta,dMmtb,dMedd,dMtran,Lxmt,Vsm[0],Vsm[1],Vsm[2],Ra,Rla,Rb,Rlb,ta_end+Tst,tb_end+Tst);
       }
       else
         fprintf(fp145,"%.3f %.1f %g %d %d %.3f %.3f %.1e %.3f %.1f %.1f %.1f %.3f %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.3f %.3f %.3f\n",
                        t+Tst,Tst,ZZ,Ka,Kb,Ma,Mb,a,e,Mzamsa,Mzamsb,a0,e0,inbin,iidd_old,
                        Vsa[0],Vsa[1],Vsa[2],Vsb[0],Vsb[1],Vsb[2],ta_end+Tst,tb_end+Tst,tdis+Tst);
       fflush(fp145);
     }

     if(OBstar==0 && OBCAT==1 && Tst<TIMEOB && (t+Tst)>TIMEOB && ((Ma>=8.0 && Ka<=9) || (Mb>=8.0 && Kb<=9))) { /* OB star output */
       OBstar=1;
       per=2.0*Pi*a*sqrt(a/(G*(Ma+Mb)));                                                /* binary period [days] */
       fprintf(fp190,"%.1f %.1f %d %d  %.3f %.3f %.3f %.3f %f  %.1f %.1f  %d %d  %f %.3f %f %.3f  %.2f %.2f  %g %g  %d %d  %.3f %.3f %.3f  %.3f %.3f %.3f  %.3f %.3f %.3f     %d %s\n",
                t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,per,ta_end+Tst,tb_end+Tst,mt,inbin,Ra,Rla,Rb,Rlb,Mzamsa,Mzamsb,La,Lb,idum_run,iidd_old,
                Vsm[0],Vsm[1],Vsm[2],Vsa[0],Vsa[1],Vsa[2],Vsb[0],Vsb[1],Vsb[2],nevroute,evroute);
       fflush(fp190);
     }



                  /* it is assumed that in one step only one process takes place: either SN or MT in that order */
                     /* this order is taken as SN takes place much faster than any MT, so even if MT begun just */
               /* prior to SN, only small amount of mass would have been lost or transfered before SN explosion */

     sn=ce=0;                                                /* in one step there may be 2 SNs, but only one MT */
     if((Ka>=13 && Ka<=15 && sna==0) || (Kb>=13 && Kb<=15 && snb==0)) {
       sn=1;
       mt=0;
       dMmta=dMmtb=0.0;                /* in SN step do not do MT component updates (in following singl() calls */

       if((Ka==15 && Kaold!=15) || (Kb==15 && Kbold!=15))         /* A/B goes SN and it leaves massless remnant */
         ;                                    /* evolve other component: no system after SN explosion of A or B */

       if((Ka==13 || Ka==14) && sna==0) {                                                /* component A goes SN */
         if(aicmark==0) {
           evroute_add("SN1 ");
           ta_end=t;                                        /* set only end time in case of SN, but not in AIC! */
         }
         tau=0.0;

         if(Ka==14 && aicmark==0 && LONGGRB==1) {                                  /* report potential long GRB */
           fprintf(fp400,"1  %.3f %.3f   %d %d %d %d  %f %f %f %f %g %f  %f %d  %f %f %f %f  %f %f %f %f  %d %d  %d %s",
                   t+Tst,Tst,Ka,Kaold,Kb,Kbold,Ma,Maold,Mb,Mbold,a,e,fraca,inbin,Mzamsa,Mzamsb,a0,e0,
                   Mheaold,Mcoaold,Mhebold,Mcobold,idum_run,iidd_old,nevroute,evroute);
           fflush(fp400);
         }

         a_0=a; e_0=e; i_0=i; Om_0=Om; om_0=om; tau_0=tau;                                 /* spin misalignemnt */

         Is=explode(Maold,Ma,Ka,fraca,Mb,Rb,Vsm,Vsa,Vsb,Vextr,&a,&e,&i,&Om,&om,&tau,Vkick,&cmt3,&texp,inbin,ecssna,fba,&jx_i,&jy_i,&jz_i,&jx_f,&jy_f,&jz_f);

         jx_0=jx_i; jy_0=jy_i; jz_0=jz_i; jx_1=jx_f; jy_1=jy_f; jz_1=jz_f;
         a_0a=a; e_0a=e; i_0a=i; Om_0a=Om; om_0a=om; tau_0a=tau;

         iA=i;                                                                      /* tilt of orbit after SN A */
         if(Ka==13) Maio=Ma;                                                /* records the formation mass of NS */

         if(BINARY==0) {                                                               /* single star evolution */
           Vsa[0]=Vkick[0];
           Vsa[1]=Vkick[1];
           Vsa[2]=Vkick[2];
         }
         else if(Is==0 && inbin==1) {         /* STOP binary evolution: system torn out or merged after SN of A */
           inbin=0;                            /* disruption: START single evolution: both stars are single now */
           a=1.0e+50;                                       /* big 'a' so program will never go to MT functions */
           e=0.0;                     /* for merger in the coming step: program stopped, after calls to sigle() */
           tdis=t;                                                                   /* records disruption time */
         }
         else if(inbin==1) {
           copyV(Vsm,VsmA);
         }
         inbinA=inbin;
         sna=1;
         tpgA=t;                                                              /* time of SNA, and preSNA data */
         copySN(Maold,Kaold,Mbold,Kbold,aold,eold,Ma,&MpgaA,&KpgaA,&MpgbA,&KpgbA,&apgA,&epgA,&MendaA);
         McheaA=Mheaold;
         MccoaA=Mcoaold;
         fracaA=fraca;
       }

       if((Kb==13 || Kb==14) && snb==0) {                                                /* component B goes SN */
         if(aicmark==0) {
           evroute_add("SN2 ");
           tb_end=t;                                        /* set only end time in case of SN, but not in AIC! */
         }
         tau=0.0;

         if(Kb==14 && aicmark==0 && LONGGRB==1) {                                  /* report potential long GRB */
           fprintf(fp400,"2  %.3f %.3f   %d %d %d %d  %f %f %f %f %g %f  %f %d  %f %f %f %f  %f %f %f %f  %d %d  %d %s\n",
                   t+Tst,Tst,Ka,Kaold,Kb,Kbold,Ma,Maold,Mb,Mbold,a,e,fracb,inbin,Mzamsa,Mzamsb,a0,e0,
                   Mheaold,Mcoaold,Mhebold,Mcobold,idum_run,iidd_old,nevroute,evroute);
           fflush(fp400);
         }

         a_1=a; e_1=e; i_1=i; Om_1=Om; om_1=om; tau_1=tau;                                 /* spin misalignemnt */

         Is=explode(Mbold,Mb,Kb,fracb,Ma,Ra,Vsm,Vsb,Vsa,Vextr,&a,&e,&i,&Om,&om,&tau,Vkick,&cmt3,&texp,inbin,ecssnb,fbb,&jx_i,&jy_i,&jz_i,&jx_f,&jy_f,&jz_f);

         jx_2=jx_f; jy_2=jy_f; jz_2=jz_f; a_2=a; e_2=e; i_2=i; Om_2=Om; om_2=om; tau_2=tau;

         iB=i;                                                                      /* tilt of orbit after SN B */
         if(Kb==13) Mbio=Mb;

         if(BINARY==0) {                                                               /* single star evolution */
           Vsb[0]=Vkick[0];
           Vsb[1]=Vkick[1];
           Vsb[2]=Vkick[2];
         }
         else if(Is==0 && inbin==1) {         /* STOP binary evolution: system torn out or merged after SN of B */
           inbin=0;                            /* disruption: START single evolution: both stars are single now */
           a=1.0e+50;                                       /* big 'a' so program will never go to MT functions */
           e=0.0;                     /* for merger in the coming step: program stopped, after calls to sigle() */
           tdis=t;                                                                   /* records disruption time */
         }
         else if(inbin==1) {
           copyV(Vsm,VsmB);
         }
         inbinB=inbin;
         snb=1;
         tpgB=t;                                                                /* time of SNB, and preSNB data */
         copySN(Maold,Kaold,Mbold,Kbold,aold,eold,Mb,&MpgaB,&KpgaB,&MpgbB,&KpgbB,&apgB,&epgB,&MendbB);
         MchebB=Mhebold;
         MccobB=Mcobold;
         fracbB=fracb;
       }
     }
     else if(inbin==1) {                                                           /* do only for binary */

       typce=1;                  /* components begin interaction, stellar and orbital parameters changed */
       Kace=Ka;
       Kbce=Kb;
       ace=a;
       ece=e;
       Mace=Ma;
       Mbce=Mb;
       Mheace=Mhea;
       Mhebce=Mheb;
       Mcoace=Mcoa;
       Mcobce=Mcob;
       Race=Ra;
       Rbce=Rb;
       Mcace=Mca;
       Mcbce=Mcb;
       cedonor=0;

       markce=0;
       if((doce==1 && ((Ka==2 && CE2==1) || (Ka==3 && CE3==1) || (Ka==4 && CE4==1))) ||
          (doce==1 && ((Kb==2 && CE2==1) || (Kb==3 && CE3==1) || (Kb==4 && CE4==1)))) {
         markce=1;
         if(SN1A==1) {
           if((Kace>9 && (Kbce==2 || Kbce==3 || Kbce==5 || Kbce==6 || Kbce==8 || Kbce==9) && (Mace+Mcbce)>1.38) ||
              (Kbce>9 && (Kace==2 || Kace==3 || Kace==5 || Kace==6 || Kace==8 || Kace==9) && (Mbce+Mcace)>1.38) ||
              ((Kace==2 || Kace==3 || Kace==5 || Kace==6 || Kace==8 || Kace==9) &&
               (Kbce==2 || Kbce==3 || Kbce==5 || Kbce==6 || Kbce==8 || Kbce==9) && (Mcace+Mcbce)>1.38)) {
             typwd=7;                                  /* Livio&Riess2003 scenario: WD + degenerate core merge in CE */
             fprintf(fp180,"%d  %f %f  %d %d  %f %f %f %f  %f %f  %f %f %f %f %f %f   %f %f %f %f %f %f  %f %f %f %f    %d %d %f %f %f %f %f  %d %d  %d %s\n",
                     typwd,t+Tst,Tst,Kace,Kbce,Mace,Mbce,ace,ece,Race,Rbce,Mheace,Mcoace,Mhebce,Mcobce,Mcace,Mcbce,
                     M0a,M0b,Mzamsa,Mzamsb,a0,e0,Mwda_0,Mwdb_0,ta_end,tb_end,
                     Ka,Kb,Ma,Mb,Ra,Rb,a,idum_run,iidd_old,nevroute,evroute);
             fflush(fp180);
           }
         }
       }


       if(markce==1) {                      /* CE aborted */
         stop2=1;                           /* STOP: non-compact object merger, non-He merger */
         typce=4;                           /* components parameters passed to record in Merger file */
         break;
       }


       dMceaold=dMcea;
       dMcebold=dMceb;
       if(Ra>=Rla && Rb<Rlb && doce==1) {                                        /* star A goes over its Roche lobe */
         evroute_add("CE1(%d-%d;",Ka,Kb);
         ce=1;
         cedonor=1;
         interaction1(Mzamsa,&M0a,&Ma,&Ka,&Kpa,&tvira,&La,&Ra,&Mca,&Mhea,&Mcoa,&flaga,&dta,&tstarta,
                      Mzamsb,&M0b,&Mb,&Kb,&Kpb,&tvirb,&Lb,&Rb,&Mcb,&Mheb,&Mcob,&flagb,&dtb,&tstartb,
                      dMmta,&dMmtb,&a,&e,&i,&Om,&om,&tau,t,dt,&ta_end,&tb_end,&stop1,&stop2,&cmt1,
                      dtmt,&dMceb);
         tlast=t; Kdonlast=Ka; Kacclast=Kb; mttypelast=6; Mdonlast=Ma; Macclast=Mb; alast=a;
         // if(PP==1) printf("RR CE1\n");
       }
       else if(Rb>=Rlb && Ra<Rla && doce==1) {                                     /* star B goes over its Roche lobe */
         evroute_add("CE2(%d-%d;",Ka,Kb);
         ce=1;
         cedonor=2;
         interaction1(Mzamsb,&M0b,&Mb,&Kb,&Kpb,&tvirb,&Lb,&Rb,&Mcb,&Mheb,&Mcob,&flagb,&dtb,&tstartb,
                      Mzamsa,&M0a,&Ma,&Ka,&Kpa,&tvira,&La,&Ra,&Mca,&Mhea,&Mcoa,&flaga,&dta,&tstarta,
                      dMmtb,&dMmta,&a,&e,&i,&Om,&om,&tau,t,dt,&tb_end,&ta_end,&stop1,&stop2,&cmt1,
                      dtmt,&dMcea);
         tlast=t; Kdonlast=Kb; Kacclast=Ka; mttypelast=6; Mdonlast=Mb; Macclast=Ma; alast=a;
         // if(PP==1) printf("CE2\n");
       }
       else if(Ra>=Rla && Rb>=Rlb) {                                /* star A and B go together over their Roche lobe */
         evroute_add("CE12(%d-%d;",Ka,Kb);
         ce=1;
         cedonor=3;
         if(((Ka==2 || Kb==2) && CE2==1) || ((Ka==3 || Kb==3) && CE3==1) || ((Ka==4 || Kb==4) && CE4==1)) {   /* CE aborted */
           stop2=1;                                                 /* STOP: non-compact object merger, non-He merger */
           typce=4;                                                 /* components parameters passed to record in Merger file */
           break;
         }
         interaction2(Mzamsa,&M0a,&Ma,&Ka,&Kpa,&tvira,&La,&Ra,&Mca,&Mhea,&Mcoa,&flaga,&dta,&tstarta,
                      Mzamsb,&M0b,&Mb,&Kb,&Kpb,&tvirb,&Lb,&Rb,&Mcb,&Mheb,&Mcob,&flagb,&dtb,&tstartb,
                      &a,&e,&i,&Om,&om,&tau,t,&ta_end,&tb_end,&stop1,&stop2);
         tlast=t; Kdonlast=Ka; Kacclast=Kb; mttypelast=7; Mdonlast=Ma; Macclast=Mb; alast=a;
         // if(PP==1) printf("RR CE12\n");
       }

       if(TranCE==1 && ce==1) {
         fprintf(fp900,"%d  %f %f  %d %d  %f %f %f %f  %f %f  %f %f %f %f %f %f   %f %f %f %f %f %f  %d %d  %d %s\n",
                 cedonor,t+Tst,Tst,Kace,Kbce,Mace,Mbce,ace,ece,Race,Rbce,Mheace,Mcoace,Mhebce,Mcobce,Mcace,Mcbce,
                 M0a,M0b,Mzamsa,Mzamsb,a0,e0,idum_run,iidd_old,nevroute,evroute);
         fflush(fp900);
       }

     }
     Rla=roche(Mb/Ma,a*(1.0-e));
     Rlb=roche(Ma/Mb,a*(1.0-e));

  //    if(PP==1) {
	// printf("END: a: %g, e: %f, wa: %g, wb: %g, Ma: %g, Mb: %f, Ka: %d, Kb: %d\n\n",a,e,wa,wb,Ma,Mb,Ka,Kb);
  //      fflush(NULL);
  //    }



     if(Ka>9 && Kb>9 && rem==0 && stop1==0 && stop2==0 && inbin==1) {     /* remnant binary: WD/NS/BH+WD/NS/BH */
                                                                                 /* sets merger time [10^6yrs] */
       if(((Ka==13 || Ka==14) && Kb>=7 && Kb<=9) || ((Kb==13 || Kb==14) && Ka>=7 && Ka<=9))   /* Helium merger */
         Tmr=0.0;
       else
         Tmr=tmerge(Ma,Mb,a,e);

       rem=1;
       if(ta_end<=tb_end)                                                     /* star A first turned to remnant */
         copyV(VsmA,Vsm0);                            /* fills Vsm0[] with center of mass velocity after 1st SN */
       else                                                                   /* star B first turned to remnant */
         copyV(VsmB,Vsm0);
     }


     if(Ma<0.08 && (Ka==0 || Ka==1)) {                           /* A on MS and M<0.08 Msun: becomes degenerate */
       if(Mzamsa<0.7)                                           /* donor becomes H-rich WD, with that ZAMS mass */
         Ka=16;                                                        /* no Helium was generated in the center */
       else                                                   /* donor becomes He-rich WD, as for higher masses */
         Ka=10;                                         /* possibly there is a significant amount of He in core */
       continue;
     }
     if(Mb<0.08 && (Kb==0 || Kb==1)) {                     /* B on MS and M<0.08 Msun, donor becomes He or H WD */
       if(Mzamsb<0.7)
         Kb=16;
       else
         Kb=10;
       continue;
     }
     if(Ma<0.35 && Ka==7) {                       /* A on He MS and M<0.35 Msun, donor becomes hybrid He/CO WD */
       Ka=17;
       continue;
     }
     if(Mb<0.35 && Kb==7) {                       /* A on He MS and M<0.35 Msun, donor becomes hybrid He/CO WD */
       Kb=17;
       continue;
     }
     if((Ma<0.01 && ((Ka>=10 && Ka<=12) || Ka==16 || Ka==17)) || (Mb<0.01 && ((Kb>=10 && Kb<=12) || Kb==16 || Kb==17))) {
       stop2=1;              /* STOP: A or B is WD and M<0.01Msun, becomes undegenerate?, merger with accretor */
       typce=5;                                       /* components parameters passed to record in Merger file */
       break;
     }


     if(inbin==1 && mttype!=6 && (fabs(dMmta)>acc || fabs(dMmtb)>acc)) {    /* mass loss/gain rate from RLOF */
       dMalla=dMmta;
       dMallb=dMmtb;
     }
     else if(inbin==1 && ce!=1) {                /* mass gain from WIND (regular and gravitationaly lensed and atmo. RLOF accretion */

       if(XEdd==2) Fbetaa=Fbetaf1(dMwinda,Ma,Ra,Rb,Mb,a,e,Ka,Kb,aspinb);           /* fraction of wind from Ma, accreted onto Mb */
       else Fbetaa=Fbetaf(Ma,Ra,Mb,a,e,Ka);
       dMwindaccb=Fbetaa*dMwinda;                /* POSITIVE number (mass loss); changed to NEGATIVE (accumulation) by dMgainf() below */
       if(XEdd==2) Fbetab=Fbetaf1(dMwindb,Mb,Rb,Ra,Ma,a,e,Kb,Ka,aspina);           /* fraction of wind from Mb, accreted onto Ma */
       else Fbetab=Fbetaf(Mb,Rb,Ma,a,e,Kb);
       dMwindacca=Fbetab*dMwindb;

       if(WindRLOF==1 && Ma<8.0 && Ka>=3 && Ka<=6) {   /* valid only for low to intermediate--mass giants: Abate et al 2013 */
         Ta=get_T(La,Ra);
         Gbetaa=Gbetaf(Ma,Mb,Ta,Ra,Rla);               /* fraction of wind from Ma, accreted onto Mb */
         dMwindrlofb=Gbetaa*dMwinda;                   /* POSITIVE number here->changed to NEGATIVE (accumulation) by dMgainf() below */
       }
       else dMwindrlofb=0.0;

       if(WindRLOF==1 && Mb<8.0 && Kb>=3 && Kb<=6) {   /* valid only for low to intermediate--mass giants: Abate et al 2013 */
         Tb=get_T(Lb,Rb);
         Gbetab=Gbetaf(Mb,Ma,Tb,Rb,Rlb);               /* fraction of wind from Mb, accreted onto Ma */
         dMwindrlofa=Gbetab*dMwindb;
       }
       else dMwindrlofa=0.0;

                                                       /* if atmospheric RLOF is present it is added/subtracted to more important wind */
       if(fabs(dMwindacca)>fabs(dMwindrlofa)) {        /* chose the more important wind accretion for star A */
         dMwindacca=dMgainf(dMwindacca,dMtran,Ma,Kb,Ka,Ra,&Mout1a,&Mout2a,&mark_outa,&doce,&merger,&mark777a,&Macc0a,aspina);
         dMalla=dMatmorlofa+dMwindacca;                /* atmospheric RLOF is already accumulation in case of accretor, and it is mass loss rate for donor */
       }
       else {
         dMwindrlofa=dMgainf(dMwindrlofa,dMtran,Ma,Kb,Ka,Ra,&Mout1a,&Mout2a,&mark_outa,&doce,&merger,&mark777a,&Macc0a,aspina);
         dMalla=dMatmorlofa+dMwindrlofa;
       }


       if(fabs(dMwindaccb)>fabs(dMwindrlofb)) {              /* chose the more important wind accretion for star B */
         dMwindaccb=dMgainf(dMwindaccb,dMtran,Mb,Ka,Kb,Rb,&Mout1b,&Mout2b,&mark_outb,&doce,&merger,&mark777b,&Macc0b,aspinb);  /* acc rate onto B */
         dMallb=dMatmorlofb+dMwindaccb;
       }
       else {
         dMwindrlofb=dMgainf(dMwindrlofb,dMtran,Mb,Ka,Kb,Rb,&Mout1b,&Mout2b,&mark_outb,&doce,&merger,&mark777b,&Macc0b,aspinb);  /* acc rate onto B */
         dMallb=dMatmorlofb+dMwindrlofb;
       }
                                                          /* atmospheric RLOF (through L1): dMatmorlof */
                                                          /* wind RLOF -- gravitational lensing: dMwindrlofa */
                                                          /* wind accreted in regular way: Bondi-Hoyle: dMwindacca */
     }
     else                                 /* disrupted binary */
       dMalla=dMallb=0.0;

   } /* end of WHILE loop */


   if(BINARY==0) {                                                                   /* output for single star evolution */
     fprintf(fp500,"%f %f  %f %d %f  %f %f  %f %d %f %f %f %d  %f %f %f  %d %d\n",
             t+Tst,Tst,Ma,Ka,Ra,ta_end,Mzamsa,MpgaA,KpgaA,McheaA,MccoaA,fracaA,ecssna,
             Vsa[0],Vsa[1],Vsa[2],idum_run,iidd_old);

     fprintf(fp500,"%f %f  %f %d %f  %f %f  %f %d %f %f %f %d  %f %f %f  %d %d\n",
             t+Tst,Tst,Mb,Kb,Rb,tb_end,Mzamsb,MpgbB,KpgbB,MchebB,MccobB,fracbB,ecssnb,
             Vsb[0],Vsb[1],Vsb[2],idum_run,iidd_old);
     fflush(fp500);
   }

   if(BINARY==1 && BINOUT==1) {                                                      /* output for binary star evolution */
     fprintf(fp510,"%f %f  %f %d %f  %f %f  %f %d %f %f %f %d  %f %f  %f %f %f  %f %f %f  %f %f %f  %f %f %f %f  %d %d  %d %d\n",
             t+Tst,Tst,Ma,Ka,Ra,ta_end,Mzamsa,MpgaA,KpgaA,McheaA,MccoaA,fracaA,ecssna,
             Tst+tpgA,MendaA,Vsm[0],Vsm[1],Vsm[2],VsmA[0],VsmA[1],VsmA[2],Vextr[0],Vextr[1],Vextr[2],
             a,e,a0,e0,inbinA,inbin,idum_run,iidd_old);

     fprintf(fp510,"%f %f  %f %d %f  %f %f  %f %d %f %f %f %d  %f %f  %f %f %f  %f %f %f  %f %f %f  %f %f %f %f  %d %d  %d %d\n",
             t+Tst,Tst,Mb,Kb,Rb,tb_end,Mzamsb,MpgbB,KpgbB,MchebB,MccobB,fracbB,ecssnb,
             Tst+tpgB,MendbB,Vsm[0],Vsm[1],Vsm[2],VsmB[0],VsmB[1],VsmB[2],Vextr[0],Vextr[1],Vextr[2],
             a,e,a0,e0,inbinB,inbin,idum_run,iidd_old);
     fflush(fp510);
   }



   if(AMCVN==1 && amcvn==1) {        /* output for AM CVn systems */
     fprintf(fp181,"%f %f %d %d %f %f %g %f %f %f    %f %f %f %f %f %f   %f %f %f %f   %f %d %d %f %f %f %f %f %f   %f %d %d %f %f %f %f %f %f   %d %d  %d %s\n",
                 t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,Ra,Rb,
                 M0a,M0b,Mzamsa,Mzamsb,a0,e0,
                 Mwda_0,Mwdb_0,ta_end,tb_end,
                 t11,Ka11,Kb11,Ma11,Mb11,a11,e11,Ra11,Rb11,
                 t22,Ka22,Kb22,Ma22,Mb22,a22,e22,Ra22,Rb22,
                 idum_run,iidd_old,nevroute,evroute);
         fflush(fp181);
   }


   if(CLUSTER==1) {           /* report population of stars from binaries for a given cluster */

     if((stop1==1 || stop2==1) && t<t_hubble && Ka!=15 && Kb!=15)         /* merger encountered */
       typcls=2;
     else if(inbin==0)                                                    /* disrupted binary */
       typcls=3;
     else                                                                 /* marks surviving binary */
       typcls=1;

     tmsa=tmsb=fmsa=fmsb=-1.0;
     if(Ka==0 || Ka==1) {
       tmsa=tmsf(Ma);                              /* MS lifetime for star A [Myr] */
       fmsa=tvira/tmsa;                            /* fractional MS lifetime for star A */
     }
     if(Kb==0 || Kb==1) {
       tmsb=tmsf(Mb);
       fmsb=tvirb/tmsb;
     }

     fprintf(fp300,"%d %f  %d %d %f %f  %f %f %f %f  %f %f %d  %f %f %f %f  %d %d  %d %d     %d %s\n",
             typcls,t, Ka,Kb,Ma,Mb, Mzamsa,Mzamsb,a0,e0, ta_end,tb_end,inbin, tmsa,fmsa,tmsb,fmsb,  stop1,stop2,
             idum_run,iidd_old,nevroute,evroute);
     fflush(fp300);
   }

                                                               /* merger encountered: record for Merger catalog */
   if(typce>0 && typce<10 && (stop1==1 || stop2==1) && t<t_hubble && Ka!=15 && Kb!=15 && Merger==1) {
     if((Kace==13 || Kace==14) && Kbce>=0 && Kbce<=9) {
       fprintf(fp800,"%d  %f %f  %d %d  %f %f %f %f  %f %f  %f %f %f %f %f %f   %f %f %f %f %f %f  %d %d  %d %s\n",
               typce,t+Tst,Tst,Kace,Kbce,Mace,Mbce,ace,ece,Race,Rbce,Mheace,Mcoace,Mhebce,Mcobce,Mcace,Mcbce,
               M0a,M0b,Mzamsa,Mzamsb,a0,e0,idum_run,iidd_old,nevroute,evroute);
       fflush(fp800);
     }
     else if((Ka==13 || Ka==14) && Kb>=0 && Kb<=9) {
       fprintf(fp800,"%d  %f %f  %d %d  %f %f %f %f  %f %f  %f %f %f %f %f %f   %f %f %f %f %f %f  %d %d  %d %s\n",
               typce,t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,Ra,Rb,Mhea,Mcoa,Mheb,Mcob,Mca,Mcb,
               M0a,M0b,Mzamsa,Mzamsb,a0,e0,idum_run,iidd_old,nevroute,evroute);
       fflush(fp800);
     }
     else if((Ka==13 || Ka==14) && Kb==-1)
       fprintf(fp800,"-1 %d  %d %d\n",typce,idum_run,iidd_old);


     if((Kbce==13 || Kbce==14) && Kace>=0 && Kace<=9) {
       fprintf(fp800,"%d  %f %f  %d %d  %f %f %f %f  %f %f  %f %f %f %f %f %f   %f %f %f %f %f %f  %d %d  %d %s\n",
               typce,t+Tst,Tst,Kbce,Kace,Mbce,Mace,ace,ece,Rbce,Race,Mhebce,Mcobce,Mheace,Mcoace,Mcbce,Mcace,
               M0b,M0a,Mzamsb,Mzamsa,a0,e0,idum_run,iidd_old,nevroute,evroute);
       fflush(fp800);
     }
     else if((Kb==13 || Kb==14) && Ka>=0 && Ka<=9) {
       fprintf(fp800,"%d  %f %f  %d %d  %f %f %f %f  %f %f  %f %f %f %f %f %f   %f %f %f %f %f %f  %d %d  %d %s\n",
               typce,t+Tst,Tst,Kb,Ka,Mb,Ma,a,e,Rb,Ra,Mheb,Mcob,Mhea,Mcoa,Mcb,Mca,
               M0b,M0a,Mzamsb,Mzamsa,a0,e0,idum_run,iidd_old,nevroute,evroute);
       fflush(fp800);
     }
     else if((Kb==13 || Kb==14) && Ka==-1)
       fprintf(fp800,"-1 %d  %d %d\n",typce,idum_run,iidd_old);
   }


   if((stop1==1 || stop2==1) && t<t_hubble && Ka!=15 && Kb!=15 && BHCAT1==1) {            /* merger encountered */
     inbin=2;                                                                             /* record for BH catalog */
     a=e=0.0;
     if(Ka<10) ta_end=t;
     if(Kb<10) tb_end=t;
     if(TIMEOUT1>(t-dt+Tst) )                                    /* merger before TIMEOUT1, merger of ANY type */
       fprintf(fp141,"%.3f %.1f %g %d %d %.3f %.3f %.1e %.3f %.1f %.1f %.1f %.3f %d %d %.1f %.1f %.1f %f %f %.3f %.3f\n",
                      t+Tst,Tst,ZZ,Ka,Kb,Ma,Mb,a,e,Mzamsa,Mzamsb,a0,e0,inbin,iidd_old,Vsm[0],Vsm[1],Vsm[2],Ra,Rb,ta_end+Tst,tb_end+Tst);
     if(TIMEOUT2>(t-dt+Tst) )
       fprintf(fp142,"%.3f %.1f %g %d %d %.3f %.3f %.1e %.3f %.1f %.1f %.1f %.3f %d %d %.1f %.1f %.1f %f %f %.3f %.3f\n",
                      t+Tst,Tst,ZZ,Ka,Kb,Ma,Mb,a,e,Mzamsa,Mzamsb,a0,e0,inbin,iidd_old,Vsm[0],Vsm[1],Vsm[2],Ra,Rb,ta_end+Tst,tb_end+Tst);
     if(TIMEOUT3>(t-dt+Tst) )
       fprintf(fp143,"%.3f %.1f %g %d %d %.3f %.3f %.1e %.3f %.1f %.1f %.1f %.3f %d %d %.1f %.1f %.1f %f %f %.3f %.3f\n",
                      t+Tst,Tst,ZZ,Ka,Kb,Ma,Mb,a,e,Mzamsa,Mzamsb,a0,e0,inbin,iidd_old,Vsm[0],Vsm[1],Vsm[2],Ra,Rb,ta_end+Tst,tb_end+Tst);
     if(TIMEOUT4>(t-dt+Tst) )
       fprintf(fp144,"%.3f %.1f %g %d %d %.3f %.3f %.1e %.3f %.1f %.1f %.1f %.3f %d %d %.1f %.1f %.1f %f %f %.3f %.3f\n",
                      t+Tst,Tst,ZZ,Ka,Kb,Ma,Mb,a,e,Mzamsa,Mzamsb,a0,e0,inbin,iidd_old,Vsm[0],Vsm[1],Vsm[2],Ra,Rb,ta_end+Tst,tb_end+Tst);
     if(TIMEOUT5>(t-dt+Tst) )
       fprintf(fp145,"%.3f %.1f %g %d %d %.3f %.3f %.1e %.3f %.1f %.1f %.1f %.3f %d %d %.1f %.1f %.1f %f %f %.3f %.3f\n",
                      t+Tst,Tst,ZZ,Ka,Kb,Ma,Mb,a,e,Mzamsa,Mzamsb,a0,e0,inbin,iidd_old,Vsm[0],Vsm[1],Vsm[2],Ra,Rb,ta_end+Tst,tb_end+Tst);
     fflush(fp141);
     fflush(fp142);
     fflush(fp143);
     fflush(fp144);
     fflush(fp145);
   }


   if(((((Ka>=10 && Ka<=14) || Ka==16 || Ka==17) && ((Kb>=10 && Kb<=14) || Kb==16 || Kb==17)) || stop1==1) && stop2==0 && t<t_hubble) {

     if(inbin==1 && rem==0) {                                                                       /* binaries */
       if(((Ka==13 || Ka==14) && Kb>=7 && Kb<=9) || ((Kb==13 || Kb==14) && Ka>=7 && Ka<=9))   /* Helium merger */
         Tmr=0.0;
       else
         Tmr=tmerge(Ma,Mb,a,e);

       if(ta_end<=tb_end)                                                     /* star A first turned to remnant */
         copyV(VsmA,Vsm0);                            /* fills Vsm0[] with center of mass velocity after 1st SN */
       else                                                                   /* star B first turned to remnant */
         copyV(VsmB,Vsm0);
     }
   }

   if(SINOUT==1 && inbin==0) {                                         /* single stars originating in binaries */
     Tmr=-1.0;
     if(ta_end<=tb_end)                                                     /* star A first turned to remnant */
       copyV(VsmA,Vsm0);                            /* fills Vsm0[] with center of mass velocity after 1st SN */
     else                                                                   /* star B first turned to remnant */
       copyV(VsmB,Vsm0);
     rec_dat1b(fp3,Ma,Mb,Vsm0,Vsm,Vsa,Vsb,Vextr,a,e,ta_end,tb_end,Tmr,Ka,Kb,MpgaA,KpgaA,MpgbA,KpgbA,apgA,
               epgA,tpgA,MendaA,MpgaB,KpgaB,MpgbB,KpgbB,apgB,epgB,tpgB,MendbB,Mzamsa,Mzamsb,a0,e0,dMcea,dMceb,
               inbinA,inbinB);
     fflush(fp3);
   }

                                                                    /* disrupted binaries with two compact obects */
   // if(0 && COMPACT==1 && (Ka==13 || Ka==14) && (Kb==13 || Kb==14) && inbin==0) {              /* output to compact.dat */
   //   per=-1.0;                                                                           /* binary period [days] */
   //   Tms=fabs(tb_end-ta_end);                                   /* time between formation of two compact objects */
   //   Tmr=-1.0;
   //
   //   if(ta_end<=tb_end)                                                     /* star A first turned to remnant */
   //     copyV(VsmA,Vsm0);                            /* fills Vsm0[] with center of mass velocity after 1st SN */
   //   else                                                                   /* star B first turned to remnant */
   //     copyV(VsmB,Vsm0);
   //
   //   fprintf(fp200,"%f %.3f  %d %d  %f %f %g %f %f  %f %f  %d %d %d  %d\n",
   //           t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,per,ta_end,tb_end,mtflag1a,mtflag1b,mtflag1ab,inbin);
   //   fprintf(fp200,"%.12f %.12f %.12f %.12f  %d %d\n",
   //           a0,e0,Mzamsa,Mzamsb,idum_run,iidd_old);
   //   fprintf(fp200,"%e %e %f %f %f %f %f %f\n",
   //           Tms,Tmr,Vsm0[0],Vsm0[1],Vsm0[2],Vsm[0],Vsm[1],Vsm[2]);
   //   fprintf(fp200,"%f %f %d %d %g %f %f   %f %f  %d\n",
   //           MpgaA,MpgbA,KpgaA,KpgbA,apgA,epgA,Tst+tpgA,MendaA,dMcea,ecssna);
   //   fprintf(fp200,"%f %f %d %d %g %f %f   %f %f  %d\n",
   //           MpgaB,MpgbB,KpgaB,KpgbB,apgB,epgB,Tst+tpgB,MendbB,dMceb,ecssnb);
   //   fprintf(fp200,"%f %d %d %d %f %f %f\n",
   //           Tst+tlast,Kdonlast,Kacclast,mttypelast,Mdonlast,Macclast,alast);
   //
   //                                                                         /* extra line of output for disrupted binaries */
   //   fprintf(fp200,"%.3f %.3f %.3f  %.3f %.3f %.3f  %.3f %.3f %.3f  %d %d\n",
   //           Vsm0[0]+Vsa[0],Vsm0[1]+Vsa[1],Vsm0[2]+Vsa[2],                       /* component A velocity after disruption */
   //           Vsm0[0]+Vsb[0],Vsm0[1]+Vsb[1],Vsm0[2]+Vsb[2],                       /* component B velocity after disruption */
   //           Vextr[0],Vextr[1],Vextr[2],inbinA,inbinB);                     /* final velocity of single star after 2nd SN */
   //
   //   fprintf(fp200,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
   //           ttms1a,tthg1a,ttrgb1a,ttcheb1a,ttagb1a,tthems1a,tthehg1a,tthergb1a,ttms1b,tthg1b,ttrgb1b,ttcheb1b,ttagb1b,tthems1b,tthehg1b,tthergb1b);
   //
   //   fprintf(fp200,"%d %s\n",nevroute,evroute);
   //   fprintf(fp200,"\n");
   //   fflush(fp200);
   // }


 } /* end of FOR loop */



 // if(PP==1) {
 //   printf("\n t: %.1f, Ma: %.3f, Mb: %.3f, Ra: %f, Rb: %f, a: %g, e: %f, Ka: %d, Kb: %d\n",t,Ma,Mb,Ra,Rb,a,e,Ka,Kb);
 //   printf(" mtflag1a: %d, mtflag1b: %d, mtflag1ab: %d, flagbra: %d, flagbrb: %d\n",mtflag1a,mtflag1b,mtflag1ab,flagbra,flagbrb);
 //   printf(" darwin: %d\n",darwin);
 //   for(j=0;j<20;j++)
 //     printf("%d ",sIbc[j]);
 //   printf("\n");
 //   for(j=0;j<20;j++)
 //     printf("%d ",sII[j]);
 //   for(j=0;j<20;j++)
 //      printf("%d ",sFa[j]);
 //   printf("\n");
 //   fflush(NULL);
 // }
 rec_dat2(fp2,j,aicns,aicbh,badorb);
 fflush(fp2);


 fclose(fp0);
 fclose(fp2);
 if(SINOUT==1) fclose(fp3);
 if(SSout==1) fclose(fp108);
 if(Xout2==1)
   fclose(fp130);
 if(Xout1==1)
   fclose(fp131);
 if(BHCAT1==1) {
   fclose(fp141);
   fclose(fp142);
   fclose(fp143);
   fclose(fp144);
   fclose(fp145);
 }
 if(WDWDout==1) fclose(fp150);
 if(CVout==1) fclose(fp151);
 if(BHCAT2==1) fclose(fp160);
 if(BHNSWD==1) fclose(fp170);
 if(SN1A==1) fclose(fp180);
 if(AMCVN==1) fclose(fp181);
 if(OBCAT==1) fclose(fp190);
 if(COMPACT==1) fclose(fp200);
 if(BHSOUT==1) fclose(fp210);
 if(CLUSTER==1) fclose(fp300);
 if(LONGGRB==1) fclose(fp400);
 if(BINARY==0) fclose(fp500);
 if(BINARY==1 && BINOUT==1) fclose(fp510);
 if(WDWD==1) fclose(fp600);
 if(RRwd==1) fclose(fp700);
 if(Merger==1) fclose(fp800);
 if(TranCE==1) fclose(fp900);
 if(PP==2) fclose(fp220);
 remove("tmp5467.dat");
 return 0;
}



/*------------------------------ BINARY EVOLUTION -----------------------------------*/

int interaction1(double Mzamsa, double *M0a, double *Ma, int *Ka, int *Kpa, double *tvira,
                 double *La, double *Ra, double *Mca, double *Mhea, double *Mcoa, int *flaga, double *dta, double *tstarta,
                 double Mzamsb, double *M0b, double *Mb, int *Kb, int *Kpb, double *tvirb, double *Lb,
                 double *Rb, double *Mcb, double *Mheb, double *Mcob, int *flagb, double *dtb, double *tstartb,
                 double dMmta, double *dMmtb, double *a, double *e, double *i, double *Om, double *om,
                 double *tau, double t, double dt, double *ta_end, double *tb_end, int *stop1, int *stop2,
                 int *cmt1, double dtmt, double *dMce)
{ /* does only CE evolution -- in one time step */
  /* star A goes over its roche lobe -- donor, star B is inside its Roche lobe: at start of MT */
  /* Ka,Kb shall not be 15 */
 int Kaold,Kbold;

 Kaold=(*Ka);   /* values at the beginig of this function */
 Kbold=(*Kb);


 if((*Ka)==15 || (*Kb)==15)
   fprintf(fp0,"error: in function interaction1() Ka or Kb equals 15!\n");

 else if((*Ka)==0 || (*Ka)==1 || (*Ka)==7) {                                                        /* donor:MS/HMS */
   (*stop2)=1;
 }

 else if((*Ka)>=2 && (*Ka)<=9 && (*Ka)!=7) {                                        /* donor: giant -- CE evolution */
   lost_env(*M0a,*Ma,*Mca,Mhea,*Mcoa,tvira,tstarta,dta,Ka,Kpa);

   if((*Kb)>=2 && (*Kb)<=9 && (*Kb)!=7) {                                                        /* double comm env */
     lost_env(*M0b,*Mb,*Mcb,Mheb,*Mcob,tvirb,tstartb,dtb,Kb,Kpb);                                      /* B changes */
     if(CE==1)
       comm_env3c0(Mzamsa,Ma,*Mca,*Ra,Kaold,Mzamsb,Mb,*Mcb,*Rb,Kbold,a,e,i,Om,om,tau);
     else if(CE==2)
       comm_env3c1(Ma,*Mca,*Ra,Kaold,Mb,*Mcb,a,e,i,Om,om,tau);
     else if(CE==3)
       comm_env3c1(Ma,*Mca,*Ra,Kaold,Mb,*Mcb,a,e,i,Om,om,tau);
     else if(CE==4)
       comm_env3c0(Mzamsa,Ma,*Mca,*Ra,Kaold,Mzamsb,Mb,*Mcb,*Rb,Kbold,a,e,i,Om,om,tau);
   }
   else if((*Kb)==0 || (*Kb)==1 || (*Kb)==7 || (*Kb)==10 || (*Kb)==11 || (*Kb)==12 || (*Kb)==16 || (*Kb)==17) {
     if(CE==1)
       comm_env3b0(Mzamsa,Ma,*Mca,*Ra,Kaold,Mb,a,e,i,Om,om,tau);                      /* single CE: B do not change */
     else if(CE==2)
       comm_env3b1(Mzamsa,Ma,*Mca,*Ra,Kaold,Mb,Kb,a,e,i,Om,om,tau,dMce);
     else if(CE==3) {
       if(*Kb<10)
         comm_env3b1(Mzamsa,Ma,*Mca,*Ra,Kaold,Mb,Kb,a,e,i,Om,om,tau,dMce);
       else
         comm_env3b0(Mzamsa,Ma,*Mca,*Ra,Kaold,Mb,a,e,i,Om,om,tau);
     }
     else if(CE==4)
       comm_env3b0(Mzamsa,Ma,*Mca,*Ra,Kaold,Mb,a,e,i,Om,om,tau);
   }
   else if((*Kb)==13 || (*Kb)==14) {                                                                   /* B changes */
     if(CE==1)
       comm_env3d1(Mzamsa,Ma,*Mca,*Ra,Mb,a,e,i,Om,om,tau,t,Kaold,Kb,dMce);     /* super acc. common env., A giant, B NS/BH */
     else if(CE==2)
       comm_env3b1(Mzamsa,Ma,*Mca,*Ra,Kaold,Mb,Kb,a,e,i,Om,om,tau,dMce);
     else if(CE==3)
       comm_env3d1(Mzamsa,Ma,*Mca,*Ra,Mb,a,e,i,Om,om,tau,t,Kaold,Kb,dMce);
     else if(CE==4)
       comm_env3d1(Mzamsa,Ma,*Mca,*Ra,Mb,a,e,i,Om,om,tau,t,Kaold,Kb,dMce);
   }
   else
     fprintf(fp0,"error: in interaction1() unexpected Kb: %d accretor type\n",(*Kb));
 }

 else if(((*Ka)>=10 && (*Ka)<=14) || (*Ka)==16 || (*Ka)==17) {    /* for WD I shouldn,t be here, NS/BH shouldn't be donor */
   fprintf(fp0,"error: in interaction1() compact donor, Ka: %d, Kb: %d, %d\n",(*Ka),(*Kb),iidd_old);
   fflush(fp0);
   (*stop1)=1;
 }
 else {
   fprintf(fp0,"error: in interaction1() unexpected types, Ka: %d, Kb: %d, %d\n",(*Ka),(*Kb),iidd_old);
   fflush(fp0);
   (*stop1)=1;
 }


 return 0;
}


int interaction2(double Mzamsa, double *M0a, double *Ma, int *Ka, int *Kpa, double *tvira, double *La, double *Ra,
                 double *Mca, double *Mhea, double *Mcoa, int *flaga, double *dta, double *tstarta,
                 double Mzamsb, double *M0b, double *Mb, int *Kb, int *Kpb, double *tvirb, double *Lb, double *Rb,
                 double *Mcb, double *Mheb, double *Mcob, int *flagb, double *dtb, double *tstartb,
                 double *a, double *e, double *i, double *Om, double *om, double *tau,
                 double t, double *ta_end, double *tb_end, int *stop1, int *stop2)
{ /* both star A and star B go over their roche lobe -- donors: at start of MT  */
  /* Ka,Kb shall not be 15 */
 int Kaold,Kbold;


 Kaold=(*Ka);
 Kbold=(*Kb);

 if((*Ka)==15 || (*Kb)==15)
   fprintf(fp0,"error: in function interaction2() Ka or Kb equals 15!\n");

 else if((*Ka)>=2 && (*Ka)<=9 && (*Ka)!=7 && (*Kb)>=2 && (*Kb)<=9 && (*Kb)!=7) {           /* two giants: double CE */
   lost_env(*M0a,*Ma,*Mca,Mhea,*Mcoa,tvira,tstarta,dta,Ka,Kpa);
   lost_env(*M0b,*Mb,*Mcb,Mheb,*Mcob,tvirb,tstartb,dtb,Kb,Kpb);
   if(CE==1)
     comm_env3c0(Mzamsa,Ma,*Mca,*Ra,Kaold,Mzamsb,Mb,*Mcb,*Rb,Kbold,a,e,i,Om,om,tau);
   else if(CE==2)
     comm_env3c1(Ma,*Mca,*Ra,Kaold,Mb,*Mcb,a,e,i,Om,om,tau);
   else if(CE==3)
     comm_env3c1(Ma,*Mca,*Ra,Kaold,Mb,*Mcb,a,e,i,Om,om,tau);
   else if(CE==4)
     comm_env3c0(Mzamsa,Ma,*Mca,*Ra,Kaold,Mzamsb,Mb,*Mcb,*Rb,Kbold,a,e,i,Om,om,tau);
 }
 else {                                                                            /* merger/contact system assumed */
   if((*ta_end)<acc)
     (*ta_end)=t;
   if((*tb_end)<acc)
     (*tb_end)=t;
   (*stop2)=1;
 }

 return 0;
}


void get_SNcounters()
{
 st0=ss0;
 st1a=ss1a; nt1a=ns1a;
 st1b=ss1b; nt1b=ns1b;
 st1c=ss1c; nt1c=ns1c;
 st2a=ss2a; nt2a=ns2a;
 st2b=ss2b; nt2b=ns2b;
 st2c=ss2c; nt2c=ns2c;
 st3a=ss3a; nt3a=ns3a;
 st3b=ss3b; nt3b=ns3b;
 st3c=ss3c; nt3c=ns3c;
 st4a=ss4a; nt4a=ns4a;
 st4b=ss4b; nt4b=ns4b;
 st4c=ss4c; nt4c=ns4c;
}

void undo_SNcounters()
{
 ss0=st0;
 ss1a=st1a; ns1a=nt1a;
 ss1b=st1b; ns1b=nt1b;
 ss1c=st1c; ns1c=nt1c;
 ss2a=st2a; ns2a=nt2a;
 ss2b=st2b; ns2b=nt2b;
 ss2c=st2c; ns2c=nt2c;
 ss3a=st3a; ns3a=nt3a;
 ss3b=st3b; ns3b=nt3b;
 ss3c=st3c; ns3c=nt3c;
 ss4a=st4a; ns4a=nt4a;
 ss4b=st4b; ns4b=nt4b;
 ss4c=st4c; ns4c=nt4c;
}


void lost_env(double M0a, double Ma, double Mca, double *Mhea, double Mcoa, double *tvira,
              double *tstarta, double *dta, int *Ka, int *Kpa)
{ /* change the parameters of star A when it loses its envelope during giant stages K=2,3,4,5,6,8,9 */
  /* parameters changed: Ka,Kpa,dta,tvira,tstarta and Mhea changed down to Mcoa for K=8,9 */
  /* important: mass Ma of losing envelope star is not changed!!! */
  /* as it is to be changed in comm_env() functions */
 double thsms,tx,tinf1,tinf2;

 if((*Ka)==2 || (*Ka)==3) {                                                                  /* donor:   HG or RG */
   if(Mca<0.35) {                                                                              /* becomes a He WD */
     (*Ka)=10;                                         /* otehrwise A becomes He MS star (7) and then Hyb WD (17) */
     (*dta)=1.0e-06;                              /* and that would be incorrect as he never burns in such a case */
   }
   else if(M0a>M_HeF) {                                                                       /* becomes: ZAMS HS */
     (*Ka)=7;
     (*tvira)=(*tstarta)=0.0;
     (*dta)=delms*thsmsf(Mca);
   }
   else {                                                                               /* becomes remnant: He WD */
     (*Kpa)=(*Ka);
     (*Ka)=-1;
     (*dta)=1.0e-06;                                            /* causes to call single() to create remnant of A */
   }
 }
 else if((*Ka)==4) {                                                                               /* donor: CHeB */
   (*Ka)=7;                                                                             /* becomes: evolved MS HS */
   thsms=thsmsf(Mca);
   (*tstarta)=thsms*((*tvira)-the1f(M0a))/thef(M0a);                      /* tvira here is from last step for K=4 */
   (*tvira)=0.0;                                                          /* tvira here is for first step for K=7 */
   (*dta)=delms*thsms;
 }
 else if((*Ka)==5) {                                                                               /* donor: EAGB */
   (*Ka)=8;                                                                             /* becomes: evolved GB HS */
   (*tstarta)=thsevolf(*Mhea,Mcoa);
   (*tvira)=0.0;
   tx=txhsf(Mca);
   tinf1=tinf1hsf(Mca);
   tinf2=tinf2hsf(Mca);
   if((*tstarta)<=tx)
     (*dta)=delrg*(tinf1-(*tstarta));
   else
     (*dta)=delrg*(tinf2-(*tstarta));
 }
 else if((*Ka)==6) {                                         /* donor: TPAGB becomes: remnant: CO/ONe WD,NS,BH,15 */
   (*Kpa)=(*Ka);       /* doesn't become HS giant as during TPAGB He and CO cores grow together and once envelope */
   (*Ka)=-1;                     /* is lost, we have bare CO core, so no He layer to form He envelope of HS giant */
   (*dta)=1.0e-06;                                              /* causes to call single() to create remnant of A */
 }
 else if((*Ka)==8 || (*Ka)==9) {                                                                  /* donor: GB HS */
   (*Kpa)=(*Ka);                                                          /* becomes: remnant: CO/ONe WD,NS,BH,15 */
   (*Ka)=-1;
   (*dta)=1.0e-06;                                              /* causes to call single() to create remnant of A */
   (*Mhea)=Mcoa;                  /* He env lost, and now for consistency with HSGBf() Mhea gets value of C0 core */
 }
 else
   fprintf(fp0,"error: in lost_env() encountered nonexistent K object's type\n");
}


void comm_env3b0(double Mzamsa, double *Ma, double Mca, double Ra, int Ka, double *Mb, double *a, double *e,
                 double *i, double *Om, double *om, double *tau)
{
/* SINGLE COMMON ENVELOPE: standard alpha/lambda prescription */
/* donor:giant/supergiant -> He core,  companion:MS,He star -> MS,He star */
/* my new compilation */
/* utrata masy z Ma, z Ma zostaje core Mca, a Mb pozostaje niezmieniony, nic nie akreuje */
 double MaEnv,R1,tmp,lambda;

 lambda=lamf(Mzamsa,*Ma,Mca,Ra,Ka);
 MaEnv=(*Ma)-Mca;                            /* envelope mass around core of star A (donor) */
 R1=roche((*Mb)/(*Ma),(*a));                 /* initial roche lobe radius of star A (donor) */
 tmp=2.0/(lambda*Alfa);

 (*a)=(Mca*(*Mb)*(*a))/(tmp*(*Ma)*MaEnv*(*a)/R1+(*Ma)*(*Mb));

 (*Ma)=Mca;                                  /* final mass of donor */
}


void comm_env3b1(double Mzamsa, double *Ma, double Mca, double Ra, int Ka, double *Mb, int *Kb, double *a, double *e, double *i,
                 double *Om, double *om, double *tau, double *dMce)
{ /* SINGLE COMMON ENVELOPE: Nelemans prescription */
  /* no hyper critical accretion, but some accretion allowed on NS/BH if set by HCE=3 */
  /* donor(A): giant -> core; companion(B): MS,HeMS,WD,NS,BH -> MS,HeMS,WD,NS,BH */
 double Mai,Mbi,Maf,Mbf,MaEnv,ai,af,dMrc;

 Mai=(*Ma);
 Maf=Mca;
 MaEnv=Mai-Maf;                                     /* available envelope mass */
 Mbi=(*Mb);
 if(HCE==3 && ((*Kb)==13 || (*Kb)==14)) {
   if((*Kb)==14)
     dMrc=BHBOND*bondi(Mzamsa,*Ma,Mca,Ra,*Mb,*a,Ka,*Kb);   /* BH accrete fraction (BHBOND) of B-H rate */
   else {
     if(NSACC==1)
       dMrc=NSBOND*bondi(Mzamsa,*Ma,Mca,Ra,*Mb,*a,Ka,*Kb); /* NS accrete fraction (NSBOND) of B-H rate */
     else
       dMrc=get_t(ADDM1,ADDM2);                     /* for NS: accrete in range ADDM1-ADDM2: set in sinbin.h */
   }
   dMrc=min(dMrc,0.5*MaEnv);                        /* do not allow accretion of more than 50% of env. mass */
   MaEnv-=(dMrc);                                   /* decrease envelope mass due to mass accreted onto NS/BH  */
 }
 else dMrc=0.0;
 Mbf=(*Mb)+dMrc;
 ai=(*a);

 af=(ai*pow(1.0-(Gamma*MaEnv)/(Mai+Mbi),2.0)*(Maf+Mbf)*Mai*Mai*Mbi*Mbi)/((Mai+Mbi)*Maf*Maf*Mbf*Mbf);

 (*dMce)+=(Mbf-Mbi);                         /* in case of more than one CE, I add up the accreted mass "+=" */
 if((*Kb)==13 && Mbf>Mmaxns)                 /* if NS accreated enough mass it becomes BH */
   (*Kb)=14;

 (*a)=af;
 (*Ma)=Maf;                                  /* final mass of donor */
 (*Mb)=Mbf;                                  /* final mass of companion */
}


void comm_env3c0(double Mzamsa, double *Ma, double Mca, double Ra, int Ka, double Mzamsb, double *Mb, double Mcb,
                 double Rb, int Kb, double *a, double *e, double *i, double *Om, double *om, double *tau)
{ /* DOUBLE COMMON ENVELOPE: standard alpha/lambda prescription */
  /* donors: giant/supergiant + giant/supergiant -> He star + He star  */
  /* my new compilation */
  /* utrata masy z Ma i z Mb, z Ma zostaje core Mca, z Mb zostaje core Mcb */
 double MaEnv,MbEnv,R1,R2,tmp1,tmp2,lambda1,lambda2;

 MaEnv=(*Ma)-Mca;                                     /* envelope mass around helium core of star A (donor) */
 MbEnv=(*Mb)-Mcb;

 lambda1=lamf(Mzamsa,*Ma,Mca,Ra,Ka);                  /* lambda for star A */
 lambda2=lamf(Mzamsb,*Mb,Mcb,Rb,Kb);                  /* lambda for star B */
 R1=roche((*Mb)/(*Ma),*a);                            /* initial roche lobe radius of star A [as a fraction of a0] */
 R2=roche((*Ma)/(*Mb),*a);
 tmp1=2.0/(lambda1*Alfa);
 tmp2=2.0/(lambda2*Alfa);

 (*a)=Mca*Mcb/(tmp1*(*Ma)*MaEnv/R1+tmp2*(*Mb)*MbEnv/R2+(*Ma)*(*Mb)/(*a));

 (*Ma)=Mca;                                           /* final mass of first donor */
 (*Mb)=Mcb;                                           /* final mass of second donor */
}


void comm_env3c1(double *Ma, double Mca, double Ra, int Ka, double *Mb, double Mcb, double *a, double *e, double *i,
                 double *Om, double *om, double *tau)
{ /* DOUBLE COMMON ENVELOPE: Nelemans prescription (no hyper critical accretion)  */
  /* donor: giant -> core,  companion: giant -> core */
 double Mai,Mbi,Maf,Mbf,MaEnv,ai,af;

 Mai=(*Ma);
 Maf=Mca;
 Mbi=(*Mb);
 Mbf=Mcb;
 MaEnv=(Mai-Maf)+(Mbi-Mbf);     /* double envelope */
 ai=(*a);

 af=(ai*pow(1.0-(Gamma*MaEnv)/(Mai+Mbi),2.0)*(Maf+Mbf)*Mai*Mai*Mbi*Mbi)/((Mai+Mbi)*Maf*Maf*Mbf*Mbf);

 (*a)=af;
 (*Ma)=Maf;                                  /* final mass of donor */
 (*Mb)=Mbf;                                  /* final mass of companion */
}


void comm_env3d0(double Mzamsa, double *Ma, double Mca, double Ra, double *Mb, double *a, double *e,
                 double *i, double *Om, double *om, double *tau, double t, int Ka, int *Kb, double *dMce)
{ /* SINGLE COMMON ENVELOPE WITH SUPER EDINGTON ACCRETION: standard alpha/lambda prescription */
  /* donor: giant/supergiant   accretor: BH or NS star  */
  /* my compilation, Bondi-Hoyle accretion (Bethe,Brown 1998, ApJ 506, 780 */
  /* transfer z Ma --> Mb, z Ma zostaje core Mca, a Mb zwieksza mase wynik jest wpisywany tez w Mb */
  /* WARNING: at March 1 2001 there was a change in tmp= (instead of Mca now is MaEnv=Ma-Mca */
  /* all results were till than WRONG!!! */
 double R1,tmp,MaEnv,Maold,Mbold,lambda;


 if((*Kb)!=13 && (*Kb)!=14)
   fprintf(fp0,"error: wrong Kb: %d in comm_env3d0()\n",(*Kb));

 Maold=(*Ma);
 Mbold=(*Mb);
 MaEnv=(*Ma)-Mca;                                     /* envelope mass around helium core of star A (donor) */

 lambda=lamf(Mzamsa,*Ma,Mca,Ra,Ka);
 R1=roche((*Mb)/(*Ma),1.0);          /* initial (in CE phase) roche lobe radius of star A [as a fraction of a0] */
 tmp=2.0*((*Ma)-Mca)/(lambda*Alfa*(*Mb)*R1)+1.0;


 if(HCE==1 && (*Kb)==13 && Ka>=0 && Ka<=6) {                    /* NS do not accreate anything from H-rich env. */
   ;
 }
 else if(HCE==2) {                                            /* do nothing to NS or BH companion--no accretion */
   ;
 }
 else {                                        /* new mass of B as calculated from formulae of Bethe,Brown 1998 */
   (*Mb)=pow(tmp,1.0/6.0)*(*Mb);
 }


 if((*Kb)==13 && (*Mb)>Mmaxns)                /* if NS accreated enough mass it becomes BH */
   (*Kb)=14;


 (*a)=(pow(tmp,-5.0/6.0)*Mca/(*Ma))*(*a);

 (*Ma)=Mca;                                /* final mass of first donor */

 if(Mbold<Mmaxns && (*Mb)<Mmaxns) {          /* NS-->NS */
   hce1++;
 }
 else if(Mbold<Mmaxns && (*Mb)>=Mmaxns) {    /* NS-->BH */
   hce2++;
 }
 else if(Mbold>=Mmaxns && (*Mb)>=Mmaxns) {   /* BH-->BH */
   hce3++;
 }

 (*dMce)+=((*Mb)-Mbold);     /* in case of more than one CE, I add up the accreted mass "+=" */
}


void comm_env3d1(double Mzamsa, double *Ma, double Mca, double Ra, double *Mb, double *a, double *e,
                 double *i, double *Om, double *om, double *tau, double t, int Ka, int *Kb, double *dMce)
{ /* SINGLE COMMON ENVELOPE WITH SUPER EDINGTON ACCRETION:  EXACT SOLUTION */
  /* 1 March 2001: standard alpha/lambda prescription */
  /* donor: giant/supergiant   accretor: BH or NS star  */
  /* transfer from Ma->Mb, Ma goes down to its core Mca */
  /* Mb increases its mass with new mass put into same Mb */
 double *ystart;
 double h1,hmin;
 double eps=1.0e-06;
 double Rce,anew,dMrc,Mbnew;
 double Maold,Mbold,MaEnv,tmp,lambda;
 int nvar,nok,nbad;


 if((*Kb)!=13 && (*Kb)!=14)
   fprintf(fp0,"error: wrong Kb: %d in comm_env3d1()\n",(*Kb));

 Maold=(*Ma);
 Mbold=(*Mb);

 lambda=lamf(Mzamsa,*Ma,Mca,Ra,Ka);
 if(HCE==3) {
   if((*Kb)==14)
     dMrc=BHBOND*bondi(Mzamsa,*Ma,Mca,Ra,*Mb,*a,Ka,*Kb);   /* for BH, accrete fraction (BHBOND) of B-H rate */
   else {
     if(NSACC==1)
       dMrc=NSBOND*bondi(Mzamsa,*Ma,Mca,Ra,*Mb,*a,Ka,*Kb); /* NS accrete fraction (NSBOND) of B-H rate */
     else
       dMrc=get_t(ADDM1,ADDM2);                     /* for NS: accrete in range ADDM1-ADDM2: set in sinbin.h */
   }
   MaEnv=(*Ma)-Mca;                                 /* available envelope mass */
   dMrc=min(dMrc,0.5*MaEnv);                        /* do not allow accretion of more than 50% of env. mass */
   MaEnv-=(dMrc);                                   /* decrease envelope mass due to mass accreted onto NS/BH  */

   Rce=roche((*Mb)/(*Ma),(*a));                     /* initial roche lobe radius of star A (donor) */
   tmp=2.0/(lambda*Alfa);
   Mbnew=(*Mb)+dMrc;                                /* new (after CE) mass of NS/BH */
   (*a)=(Mca*Mbnew*(*a))/(tmp*(*Ma)*MaEnv*(*a)/Rce+(*Ma)*(*Mb));
   (*Ma)=Mca;                                       /* Ma change */
   (*Mb)=Mbnew;                                     /* Mb change */
 }
 else if(HCE==1 && (*Kb)==13 && Ka>=0 && Ka<=6) {    /* HCE=1 mod: NS do not accreate anything from H-rich env. */
   MaEnv=(*Ma)-Mca;                             /* do single CE as in comm_env3b0() */
   Rce=roche((*Mb)/(*Ma),(*a));                 /* initial roche lobe radius of star A (donor) */
   tmp=2.0/(lambda*Alfa);
   (*a)=(Mca*(*Mb)*(*a))/(tmp*(*Ma)*MaEnv*(*a)/Rce+(*Ma)*(*Mb));
   (*Ma)=Mca;                                   /* Ma change and Mb do not change !!! */
 }
 else if(HCE==2) {                                 /* HCE=2 mod: NS or BH do not accreate anything at all */
   MaEnv=(*Ma)-Mca;                                /* do single CE as in comm_env3b0() */
   Rce=roche((*Mb)/(*Ma),(*a));                    /* initial roche lobe radius of star A (donor) */
   tmp=2.0/(lambda*Alfa);
   (*a)=(Mca*(*Mb)*(*a))/(tmp*(*Ma)*MaEnv*(*a)/Rce+(*Ma)*(*Mb));
   (*Ma)=Mca;                                 /* Ma change and Mb do not change !!! */
 }
 else {
  h1=((*Ma)-Mca)/1000.0;    /* first step */
  hmin=0.0;                 /* minimum step */
  kmax=0;                   /* no intermidiate data storage */
  nvar=2;
  ystart=dvector(1,nvar);
  Rce=roche((*Mb)/(*Ma),(*a));             /* initial (in CE phase) roche lobe radius of star A */
  anew=(2.0*((*Ma)-Mca)+lambda*Alfa*(*Mb))/(2.0*((*Ma)-Mca)/Rce+lambda*Alfa*(*Mb)/(*a));
  ystart[1]=anew;
  ystart[2]=(*Mb);
  odeint(ystart,nvar,*Ma,Mca,eps,h1,hmin,&nok,&nbad,derivs,rkqs,Mca,lambda);
  (*a)=ystart[1];
  (*Mb)=ystart[2];                            /* Mb change */
  (*Ma)=Mca;                                  /* Ma change */
 }

 if((*Kb)==13 && (*Mb)>Mmaxns)               /* if NS accreated enough mass it becomes BH */
   (*Kb)=14;

 if(Mbold<Mmaxns && (*Mb)<Mmaxns) {         /* NS-->NS */
   hce1++;
 }
 else if(Mbold<Mmaxns && (*Mb)>=Mmaxns) {   /* NS-->BH */
   hce2++;
 }
 else if(Mbold>=Mmaxns && (*Mb)>=Mmaxns) {  /* BH-->BH */
   hce3++;
 }

 (*dMce)+=((*Mb)-Mbold);     /* in case of more than one CE, I add up the accreted mass "+=" */
}


double bondi(double Mzamsa, double Ma, double Mca, double Ra, double Mb, double a, int Ka, int Kb)
{ /* returns estimate (exact solution) of mass accreted in CE [Msun]: Bondi-Hoyle rate used  */
  /* A: donor: giant/supergiant,  B: accretor: BH or NS star  */
  /* no parameters are changed here */
 double *ystart,h1,hmin;
 double eps=1.0e-06,Rce,anew;
 double dMce,Mbold,Menv,lambda;
 int nvar,nok,nbad;


 if(Kb!=13 && Kb!=14)
   fprintf(fp0,"error: wrong Kb: %d in bondi()\n",Kb);

 Mbold=Mb;
 Menv=Ma-Mca;
 lambda=lamf(Mzamsa,Ma,Mca,Ra,Ka);

 h1=(Ma-Mca)/1000.0;       /* first step */
 hmin=0.0;                 /* minimum step */
 kmax=0;                   /* no intermidiate data storage */
 nvar=2;
 ystart=dvector(1,nvar);
 Rce=roche(Mb/Ma,a);             /* initial (in CE phase) roche lobe radius of star A */
 anew=(2.0*(Ma-Mca)+lambda*Alfa*Mb)/(2.0*(Ma-Mca)/Rce+lambda*Alfa*Mb/a);
 ystart[1]=anew;
 ystart[2]=Mb;
 odeint(ystart,nvar,Ma,Mca,eps,h1,hmin,&nok,&nbad,derivs,rkqs,Mca,lambda);
 a=ystart[1];
 Mb=ystart[2];                            /* Mb change */
 Ma=Mca;                                  /* Ma change */
 dMce=(Mb-Mbold);     /* estimate of mass accreted onto NS/BH */

 if(dMce>Menv)
   fprintf(fp0,"error: bondi() accreted mass (%f) larger than envelope mass (%f) (%d)\n",dMce,Menv,iidd_old);

 return dMce;
}


double lamf(double Mzamsa, double Ma, double Mca, double Ra, int Ka)
{ /* returns lambda for CE energy prescription for star A */
  /* H stars: either uses set lambda (see Lambda value in sinbin.h) or calculates lambda from Chinese: */
  /* X.-J. Xu and X.-D. Li arXiv:1004.4957 (v1, 28Apr2010) -- set by chlambda in sinbin.h */
  /* He: stars: from Natasha calculations */
  /* all above lambdas are multiplied by fudge factor golambda (see sinbin.h) */
 double lambda,lamb,lamg,maxb,maxg;
 double a0,a1,a2,a3,a4,a5;
 double b0,b1,b2,b3,b4,b5;
 double x,y1,y2,Zlimit,Menva;
 double Rmin,Rmax;
 int manual;


 if(Ka==8 || Ka==9) {                                                          /* Helium stars: always Natasha fit */
   Rmin=0.25;                                                                  /* minimum considered radius: Natasha */
   Rmax=120.0;                                                                 /* maximum considered radius: Natasha */
   if(Ra<Rmin)
     lambda=0.3*pow(Rmin,-0.8);
   else if(Ra>Rmax)
     lambda=0.3*pow(Rmax,-0.8);
   else
     lambda=0.3*pow(Ra,-0.8);
 }
 else if(chlambda>=1 && chlambda<=3 && Ka>=2 && Ka<=6) {                       /* H-rich stars NEW APPROACH: Chinese */
   Zlimit=0.0105;                                                              /* limiting value of Z for Chinese models */
   Menva=Ma-Mca;                                                               /* current He core mass */
   manual=0;                                                                   /* then approximated by hand: manual=1 */

   if(ZZ>Zlimit)  {                                                            /* Z>0.5 Zsun: popI */
     if(Mzamsa<1.5) {
       maxb=2.5;
       maxg=1.5;

       if(Ra>200.0) {
         lamb=0.05;
         lamg=0.05;
         manual=1;
       }
       else if((Ka==2 || Ka==3) && Ra > 2.7) {
         lamb=-Ra*9.18e-03+2.33;
         lamg=-Ra*4.59e-03+1.12;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=8.35897;  a1=-18.89048; a2=10.47651; a3=0.99352;  a4=0.0; a5=0.0;  /* coeffs: with internal */
         b0=17.58328; b1=-34.84355; b2=10.70536; b3=8.49042;  b4=0.0; b5=0.0;  /* coeffs: no internal */
       }
       else if(Ka==4) {
         a0=46.00978; a1=-298.64993; a2=727.40936; a3=-607.66797; a4=0.0; a5=0.0;
         b0=63.61259; b1=-399.89494; b2=959.62055; b3=-795.20699; b4=0.0; b5=0.0;
       }
       else {
         lamb=-Ra*3.57e-04+0.1;
         lamg=-Ra*3.57e-04+0.1;
         manual=1;
       }
     }
     else if(Mzamsa<2.5) {
       maxb=4.0;
       maxg=2.0;

       if(Ra>340.0) {
         lamb=3.589970;
         lamg=0.514132;
         manual=1;
       }
       else if(Ka==4 && Ra > 8.5 && Ra < 60.0) {
         lamb=3.0;
         lamg=1.2;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=2.05363; a1=-0.00685; a2=-3.42739e-04; a3=3.93987e-06; a4=-1.18237e-08; a5=0.0;
         b0=1.07658; b1=-0.01041; b2=-4.90553e-05; b3=1.13528e-06; b4=-3.91609e-09; b5=0.0;
       }
       else if(Ka==4) {
         a0=34.41826; a1=-6.65259; a2=0.43823; a3=-0.00953; a4=0.0; a5=0.0;
         b0=13.66058; b1=-2.48031; b2=0.15275; b3=-0.00303; b4=0.0; b5=0.0;
       }
       else {
         a0=0.88954; a1=0.0098;  a2=-3.1411e-05;  a3=7.66979e-08; a4=0.0;         a5=0.0;
         b0=0.48271; b1=0.00584; b2=-6.22051e-05; b3=2.41531e-07; b4=-3.1872e-10; b5=0.0;
       }
     }
     else if(Mzamsa<3.5) {
       maxb=500.0;
       maxg=10.0;

       if(Ra>400.0) {
         lamb=116.935557;
         lamg=0.848808;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=2.40831; a1=-0.42459; a2=0.03431; a3=-9.26879e-04; a4=8.24522e-06; a5=0.0;
         b0=1.30705; b1=-0.22924; b2=0.01847; b3=-5.06216e-04; b4=4.57098e-06; b5=0.0;
         maxb=2.5; maxg=1.5;
       }
       else if(Ka==4) {
         a0=-42.98513; a1=7.90134; a2=-0.54646; a3=0.01863; a4=-3.13101e-04; a5=2.07468e-06;
         b0=-6.73842;  b1=1.06656; b2=-0.05344; b3=0.00116; b4=-9.34446e-06; b5=0.0;
         maxb=2.5; maxg=1.5;
       }
       else {
         a0=-0.04669; a1=0.00764; a2=-4.32726e-05; a3=9.31942e-08; a4=0.0;         a5=0.0;
         b0=0.44889;  b1=0.01102; b2=-6.46629e-05; b3=5.66857e-09; b4=7.21818e-10; b5=-1.2201e-12;
       }
     }
     else if(Mzamsa<4.5) {
       maxb=1000.0;
       maxg=8.0;

       if(Ra>410.0) {
         lamb=52.980056;
         lamg=1.109736;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.8186;  a1=-0.17464; a2=0.00828; a3=-1.31727e-04; a4=7.08329e-07; a5=0.0;
         b0=1.02183; b1=-0.1024;  b2=0.00493; b3=-8.16343e-05; b4=4.55426e-07; b5=0.0;
         maxb=2.5; maxg=1.5;
       }
       else if(Ka==4) {
         a0=-7.3098;  a1=0.56647; a2=-0.01176; a3=7.90112e-05; a4=0.0; a5=0.0;
         b0=-3.80455; b1=0.29308; b2=-0.00603; b3=4.00471e-05; b4=0.0; b5=0.0;
         maxb=2.5; maxg=1.5;
       }
       else {
         a0=-0.37322; a1=0.00943; a2=-3.26033e-05; a3=5.37823e-08; a4=0.0; a5=0.0;
         b0=0.13153;  b1=0.00984; b2=-2.89832e-05; b3=2.63519e-08; b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<5.5) {
       maxb=1000.0;
       maxg=8.0;

       if(Ra>430.0) {
         lamb=109.593522;
         lamg=1.324248;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.52581; a1=-0.08125; a2=0.00219; a3=-2.0527e-05;  a4=6.79169e-08; a5=0.0;
         b0=0.85723; b1=-0.04922; b2=0.00137; b3=-1.36163e-05; b4=4.68683e-08; b5=0.0;
       }
       else if(Ka==4) {
         a0=-9.93647; a1=0.42831; a2=-0.00544; a3=2.25848e-05; a4=0.0; a5=0.0;
         b0=-5.33279; b1=0.22728; b2=-0.00285; b3=1.16408e-05; b4=0.0; b5=0.0;
       }
       else {
         a0=-0.80011; a1=0.00992; a2=-3.03247e-05; a3=5.26235e-08;  a4=0.0; a5=0.0;
         b0=-0.00456; b1=0.00426; b2=4.71117e-06;  b3=-1.72858e-08; b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<6.5) {
       maxb=25.5;
       maxg=5.0;

       if(Ra>440.0) {
         lamb=16.279603;
         lamg=1.352166;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.41601; a1=-0.04965; a2=8.51527e-04; a3=-5.54384e-06; a4=1.32336e-08; a5=0.0;
         b0=0.78428; b1=-0.02959; b2=5.2013e-04;  b3=-3.45172e-06; b4=8.17248e-09; b5=0.0;
       }
       else if(Ka==4) {
         a0=13.91465; a1=-0.55579; a2=0.00809; a3=-4.94872e-05; a4=1.08899e-07; a5=0.0;
         b0=7.68768;  b1=-0.30723; b2=0.00445; b3=-2.70449e-05; b4=5.89712e-08; b5=0.0;
       }
       else {
         a0=-2.7714; a1=0.06467;  a2=-4.01537e-04; a3=7.98466e-07;  a4=0.0; a5=0.0;
         b0=0.23083; b1=-0.00266; b2=2.21788e-05;  b3=-2.35696e-08; b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<7.5) {
       maxb=9.0;
       maxg=3.0;

       if(Ra>420.0) {
         lamb=5.133959;
         lamg=1.004036;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.38344; a1=-0.04093; a2=5.78952e-04; a3=-3.19227e-06; a4=6.40902e-09; a5=0.0;
         b0=0.76009; b1=-0.02412; b2=3.47104e-04; b3=-1.92347e-06; b4=3.79609e-09; b5=0.0;
       }
       else if(Ka==4) {
         a0=4.12387; a1=-0.12979; a2=0.00153;     a3=-7.43227e-06; a4=1.29418e-08; a5=0.0;
         b0=2.18952; b1=-0.06892; b2=8.00936e-04; b3=-3.78092e-06; b4=6.3482e-09;  b5=0.0;
       }
       else {
         a0=-0.63266; a1=0.02054;  a2=-1.3646e-04; a3=2.8661e-07;   a4=0.0; a5=0.0;
         b0=0.26294;  b1=-0.00253; b2=1.32272e-05; b3=-7.12205e-09; b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<8.5) {
       maxb=7.0;
       maxg=3.0;

       if(Ra>490.0) {
         lamb=4.342985;
         lamg=0.934659;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.35516; a1=-0.03414; a2=4.02065e-04; a3=-1.85931e-06; a4=3.08832e-09; a5=0.0;
         b0=0.73826; b1=-0.01995; b2=2.37842e-04; b3=-1.09803e-06; b4=1.79044e-09; b5=0.0;
       }
       else if(Ka==4) {
         a0=-3.89189; a1=0.19378; a2=-0.0032;  a3=2.39504e-05; a4=-8.28959e-08; a5=1.07843e-10;
         b0=-2.24354; b1=0.10918; b2=-0.00179; b3=1.33244e-05; b4=-4.57829e-08; b5=5.90313e-11;
         maxb=1.0; maxg=0.5;
       }
       else {
         a0=-0.1288; a1=0.0099;   a2=-6.71455e-05;  a3=1.33568e-07;   a4=0.0; a5=0.0;
         b0=0.26956; b1=-0.00219; b2=7.97743e-06;   b3=-1.53296e-09;  b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<9.5) {
       maxb=4.0;
       maxg=2.0;

       if(Ra>530.0) {
         lamb=2.441672;
         lamg=0.702310;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.32549; a1=-0.02845; a2=2.79097e-04; a3=-1.07254e-06; a4=1.46801e-09; a5=0.0;
         b0=0.71571; b1=-0.01657; b2=1.64607e-04; b3=-6.31935e-07; b4=8.52082e-10; b5=0.0;
       }
       else if(Ka==4) {
         a0=0.86369; a1=-0.00995; a2=4.80837e-05;  a3=-6.10454e-08; a4=-2.79504e-12; a5=0.0;
         b0=-0.7299; b1=0.0391;   b2=-5.78132e-04; b3=3.7072e-06;   b4=-1.07036e-08; b5=1.14833e-11;
       }
       else {
         a0=1.19804; a1=-0.01961; a2=1.28222e-04; a3=-3.41278e-07; a4=3.35614e-10; a5=0.0;
         b0=0.40587; b1=-0.0051;  b2=2.73866e-05; b3=-5.74476e-08; b4=4.90218e-11; b5=0.0;
       }
     }
     else if(Mzamsa<11.0) {
       maxb=3.0;
       maxg=1.5;

       if(Ra>600.0) {
         lamb=1.842314;
         lamg=0.593854;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.29312; a1=-0.02371; a2=1.93764e-04; a3=-6.19576e-07; a4=7.04227e-10; a5=0.0;
         b0=0.69245; b1=-0.01398; b2=1.17256e-04; b3=-3.81487e-07; b4=4.35818e-10; b5=0.0;
         maxb=1.0; maxg=0.6;
       }
       else if(Ka==4) {
         a0=0.74233; a1=-0.00623; a2=2.04197e-05; a3=-1.30388e-08; a4=0.0; a5=0.0;
         b0=0.36742; b1=-0.00344; b2=1.27838e-05; b3=-1.0722e-08;  b4=0.0; b5=0.0;
       }
       else {
         a0=0.3707;  a1=2.67221e-04; a2=-9.86464e-06; a3=2.26185e-08; a4=0.0; a5=0.0;
         b0=0.25549; b1=-0.00152;    b2=3.35239e-06;  b3=2.24224e-10; b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<13.0) {
       maxb=1.5;
       maxg=1.0;

       if(Ra>850.0) {
         lamb=0.392470;
         lamg=0.176660;
         manual=1;
       }
       else if(Ra>0.0 && Ra<=350.0) {
         a0=1.28593; a1=-0.02209; a2=1.79764e-04; a3=-6.21556e-07; a4=7.59444e-10; a5=0.0;
         b0=0.68544; b1=-0.01394; b2=1.20845e-04; b3=-4.29071e-07; b4=5.29169e-10; b5=0.0;
       }
       else if(Ra>350.0 && Ra<=600.0) {
         a0=-11.99537; a1=0.0992;  a2=-2.8981e-04; a3=3.62751e-07;  a4=-1.65585e-10; a5=0.0;
         b0=0.46156;   b1=-0.0066; b2=3.9625e-05;  b3=-9.98667e-08; b4=-8.84134e-11; b5=0.0;
       }
       else {
         a0=-58.03732; a1=0.23633; a2=-3.20535e-04; a3=1.45129e-07; a4=0.0; a5=0.0;
         b0=-15.11672; b1=0.06331; b2=-8.81542e-05; b3=4.0982e-08;  b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<15.0) {
       maxb=1.5;
       maxg=1.0;

       if(Ra>1000.0) {
         lamb=0.414200;
         lamg=0.189008;
         manual=1;
       }
       else if(Ka==4 && Ra > 69.0 && Ra < 126.0) {
         lamb=Ra*(-8.77e-04)+0.5;
         lamg=0.18;
         manual=1;
       }
       else if((Ka==2 || Ka==3) && Ra>190.0 && Ra<600.0) {
         lamb=0.15;
         lamg=0.15;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.39332;  a1=-0.0318; a2=3.95917e-04; a3=-2.23132e-06; a4=4.50831e-09; a5=0.0;
         b0=0.78215; b1=-0.02326; b2=3.25984e-04; b3=-1.94991e-06; b4=4.08044e-09; b5=0.0;
       }
       else if(Ka==4) {
         a0=1.12889; a1=-0.00901; a2=3.04077e-05; a3=-4.31964e-08; a4=2.14545e-11; a5=0.0;
         b0=0.568;   b1=-0.0047;  b2=1.57818e-05; b3=-2.21207e-08; b4=1.08472e-11; b5=0.0;
       }
       else {
         a0=-106.90553; a1=0.36469; a2=-4.1472e-04;  a3=1.57349e-07; a4=0.0; a5=0.0;
         b0=-39.93089;  b1=0.13667; b2=-1.55958e-04; b3=5.94076e-08; b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<18.0) {
       maxb=1.5;
       maxg=1.0;

       if(Ra>1050.0) {
         lamb=0.2;
         lamg=0.1;
         manual=1;
       }
       else if((Ka==2 || Ka==3) && Ra>120.0 && Ra<170.0) {
         lamb=0.2;
         lamg=0.2;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.43177; a1=-0.03533; a2=5.11128e-04; a3=-3.57633e-06; a4=9.36778e-09; a5=0.0;
         b0=0.85384; b1=-0.03086; b2=5.50878e-04; b3=-4.37671e-06; b4=1.25075e-08; b5=0.0;
       }
       else if(Ka==4) {
         a0=0.84143; a1=-0.00576; a2=1.68854e-05; a3=-2.0827e-08;  a4=8.97813e-12; a5=0.0;
         b0=0.36014; b1=-0.00254; b2=7.49639e-06; b3=-9.20103e-09; b4=3.93828e-12; b5=0.0;
       }
       else {
         a0=-154.70559; a1=0.46718; a2=-4.70169e-04; a3=1.57773e-07; a4=0.0; a5=0.0;
         b0=-65.39602;  b1=0.19763; b2=-1.99078e-04; b3=6.68766e-08; b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<35.0) {
       maxb=1.5;
       maxg=1.0;

       if(Ra>1200.0) {
         lamb=0.05;
         lamg=0.05;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         lamb=1.2*exp(-Ra/90.0);
         lamg=0.55*exp(-Ra/160.0);
         manual=1;
       }
       else if(Ka==4) {
         a0=0.48724; a1=-0.00177;    a2=2.60254e-06; a3=-1.25824e-09; a4=0.0; a5=0.0;
         b0=0.22693; b1=-8.7678e-04; b2=1.28852e-06; b3=-6.12912e-10; b4=0.0; b5=0.0;
       }
       else {
         a0=-260484.85724; a1=4.26759e+06; a2=-2.33016e+07; a3=4.24102e+07; a4=0.0; a5=0.0;
         b0=-480055.67991; b1=7.87484e+6;  b2=-4.30546e+07; b3=7.84699e+07; b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<75.0) {
       maxb=1.0;
       maxg=0.5;

       a0=0.31321; a1=-7.50384e-04; a2=5.38545e-07; a3=-1.16946e-10; a4=0.0; a5=0.0;
       b0=0.159;   b1=-3.94451e-04; b2=2.88452e-07; b3=-6.35132e-11; b4=0.0; b5=0.0;
     }
     else {
       maxb=1.0;
       maxg=0.5;

       a0=0.376;  a1=-0.0018;  a2=2.81083e-06; a3=-1.67386e-09; a4=3.35056e-13; a5=0.0;
       b0=0.2466; b1=-0.00121; b2=1.89029e-06; b3=-1.12066e-09; b4=2.2258e-13;  b5=0.0;
     }

   }
   else {                                                                       /* Z<=0.5 Zsun: popI and popII */
     if(Mzamsa<1.5) {
       maxb=2.0;
       maxg=1.5;

       if(Ra>160.0) {
         lamb=0.05;
         lamg=0.05;
         manual=1;
       }
       else if((Ka==2 || Ka==3) && Ra>12.0) {
         lamb=1.8*exp(-Ra/80.0);
         lamg=exp(-Ra/45.0);
         manual=1;
       }
       else if(Ka==4) {
         a0=0.37294; a1=-0.05825; a2=0.00375; a3=-7.59191e-05; a4=0.0; a5=0.0;
         b0=0.24816; b1=-0.04102; b2=0.0028;  b3=-6.20419e-05; b4=0.0; b5=0.0;
       }
       else {
         a0=0.24012; a1=-0.01907; a2=6.09529e-04; a3=-8.17819e-06; a4=4.83789e-08; a5=-1.04568e-10;
         b0=0.15504; b1=-0.01238; b2=3.96633e-04; b3=-5.3329e-06;  b4=3.16052e-08; b5=-6.84288e-11;
       }
     }
     else if(Mzamsa<2.5) {
       maxb=4.0;
       maxg=2.0;

       if(Ra>350.0) {
         lamb=2.868539;
         lamg=0.389991;
         manual=1;
       }
       else if((Ka==2 || Ka==3) && Ra>22.0 && Ra<87.0) {
         lamb=1.95;
         lamg=0.85;
         manual=1;
       }
       else if(Ka==4 && Ra>6.0 && Ra<50.0) {
         lamb=0.8;
         lamg=0.35;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=2.56108; a1=-0.75562; a2=0.1027;  a3=-0.00495; a4=8.05436e-05; a5=0.0;
         b0=1.41896; b1=-0.4266;  b2=0.05792; b3=-0.00281; b4=4.61e-05;    b5=0.0;
       }
       else if(Ka==4) {
         a0=-103.92538; a1=25.37325; a2=-2.03273; a3=0.0543;  a4=0.0; a5=0.0;
         b0=-56.03478;  b1=13.6749;  b2=-1.09533; b3=0.02925; b4=0.0; b5=0.0;
       }
       else {
         a0=0.5452;  a1=0.00212;      a2=6.42941e-05; a3=-1.46783e-07; a4=0.0;        a5=0.0;
         b0=0.30594; b1=-9.58858e-04; b2=1.12174e-04; b3=-1.04079e-06; b4=3.4564e-09; b5=-3.91536e-12;
       }
     }
     else if(Mzamsa<3.5) {
       maxb=600.0;
       maxg=2.0;

       if(Ra>400.0) {
         lamb=398.126442;
         lamg=0.648560;
         manual=1;
       }
       else if(Ka==4 && Ra>36.0 && Ra<53.0) {
         lamb=1.0;
         lamg=1.0;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.7814;  a1=-0.17138; a2=0.00754; a3=-9.02652e-05; a4=0.0; a5=0.0;
         b0=0.99218; b1=-0.10082; b2=0.00451; b3=-5.53632e-05; b4=0.0; b5=0.0;
       }
       else if(Ka==4) {
         a0=-12.40832; a1=1.59021; a2=-0.06494; a3=8.69587e-04; a4=0.0; a5=0.0;
         b0=-6.47476;  b1=0.8328;  b2=-0.03412; b3=4.58399e-04; b4=0.0; b5=0.0;
       }
       else {
         a0=-0.475;  a1=-0.00328; a2=1.31101e-04; a3=-6.03669e-07; a4=8.49549e-10; a5=0.0;
         b0=0.05434; b1=0.0039;   b2=9.44609e-06; b3=-3.87278e-08; b4=0.0;         b5=0.0;
       }
     }
     else if(Mzamsa<4.5) {
       maxb=600.0;
       maxg=2.0;

       if(Ra>410.0) {
         lamb=91.579093;
         lamg=1.032432;
         manual=1;
       }
       else if(Ka==4 && Ra>19.0 && Ra<85.0) {
         lamb=0.255;
         lamg=0.115;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.65914; a1=-0.10398; a2=0.0029;  a3=-2.24862e-05; a4=0.0; a5=0.0;
         b0=0.92172; b1=-0.06187; b2=0.00177; b3=-1.42677e-05; b4=0.0; b5=0.0;
       }
       else if(Ka==4) {
         a0=-5.89253; a1=0.54296; a2=-0.01527; a3=1.38354e-04; a4=0.0; a5=0.0;
         b0=-3.21299; b1=0.29583; b2=-0.00833; b3=7.55646e-05; b4=0.0; b5=0.0;
       }
       else {
         a0=-0.2106; a1=-0.01574; a2=2.01107e-04; a3=-6.90334e-07; a4=7.92713e-10; a5=0.0;
         b0=0.36779; b1=-0.00991; b2=1.19411e-04; b3=-3.59574e-07; b4=3.33957e-10; b5=0.0;
       }
     }
     else if(Mzamsa<5.5) {
       maxb=10.0;
       maxg=3.0;

       if(Ra>320) {
         lamb=7.618019;
         lamg=1.257919;
         manual=1;
       }
       else if(Ka==4 && Ra>85.0 && Ra<120.0) {
         lamb=0.4;
         lamg=0.1;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.58701; a1=-0.06897; a2=0.00129;     a3=-6.99399e-06; a4=0.0; a5=0.0;
         b0=0.87647; b1=-0.04103; b2=7.91444e-04; b3=-4.41644e-06; b4=0.0; b5=0.0;
       }
       else if(Ka==4) {
         a0=-0.67176; a1=0.07708; a2=-0.00175;    a3=1.1991e-05;  a4=0.0; a5=0.0;
         b0=-0.38561; b1=0.0427;  b2=-9.6948e-04; b3=6.64455e-06; b4=0.0; b5=0.0;
       }
       else {
         a0=-0.12027; a1=0.01981;  a2=-2.27908e-04; a3=7.55556e-07;  a4=0.0; a5=0.0;
         b0=0.31252;  b1=-0.00527; b2=3.60348e-05;  b3=-3.22445e-08; b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<6.5) {
       maxb=4.0;
       maxg=1.5;

       if(Ra>330.0) {
         lamb=2.390575;
         lamg=0.772091;
         manual=1;
       }
       else if(Ka==4 && Ra>115.0 && Ra<165.0) {
         lamb=0.2;
         lamg=0.1;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.527;   a1=-0.04738; a2=6.1373e-04;  a3=-2.36835e-06; a4=0.0; a5=0.0;
         b0=0.83636; b1=-0.02806; b2=3.73346e-04; b3=-1.47016e-06; b4=0.0; b5=0.0;
       }
       else if(Ka==4) {
         a0=0.30941; a1=0.00965; a2=-2.31975e-04; a3=1.26273e-06; a4=0.0; a5=0.0;
         b0=0.14576; b1=0.00562; b2=-1.30273e-04; b3=7.06459e-07; b4=0.0; b5=0.0;
       }
       else {
         a0=0.26578; a1=0.00494;  a2=-7.02203e-05; a3=2.25289e-07; a4=0.0; a5=0.0;
         b0=0.26802; b1=-0.00248; b2=6.45229e-06;  b3=1.69609e-08; b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<7.5) {
       maxb=2.5;
       maxg=1.0;

       if(Ra>360.0) {
         lamb=1.878174;
         lamg=0.646353;
         manual=1;
       }
       else if(Ka==4 && Ra>150.0 && Ra<210.0) {
         lamb=0.2;
         lamg=0.1;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.49995; a1=-0.03921; a2=4.2327e-04; a3=-1.37646e-06; a4=0.0; a5=0.0;
         b0=0.81688; b1=-0.02324; b2=2.5804e-04; b3=-8.54696e-07; b4=0.0; b5=0.0;
       }
       else if(Ka==4) {
         a0=0.44862; a1=0.00234; a2=-9.23152e-05; a3=4.67797e-07; a4=0.0; a5=0.0;
         b0=0.21873; b1=0.00154; b2=-5.18806e-05; b3=2.60283e-07; b4=0.0; b5=0.0;
       }
       else {
         a0=0.8158;  a1=-0.01633; a2=1.46552e-04; a3=-5.75308e-07; a4=8.77711e-10; a5=0.0;
         b0=0.26883; b1=-0.00219; b2=4.12941e-06; b3=1.33138e-08;  b4=0.0;         b5=0.0;
       }
     }
     else if(Mzamsa<8.5) {
       maxb=2.0;
       maxg=1.0;

       if(Ra>400.0) {
         lamb=1.517662;
         lamg=0.553169;
         manual=1;
       }
       else if(Ka==4 && Ra>190.0 && Ra<260.0) {
         lamb=0.2;
         lamg=0.1;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.46826; a1=-0.03184; a2=2.85622e-04; a3=-7.91228e-07; a4=0.0; a5=0.0;
         b0=0.79396; b1=-0.01903; b2=1.77574e-04; b3=-5.04262e-07; b4=0.0; b5=0.0;
       }
       else if(Ka==4) {
         a0=0.50221; a1=-3.19021e-04; a2=-3.81717e-05; a3=1.80726e-07; a4=0.0; a5=0.0;
         b0=0.24748; b1=-9.9338e-05;  b2=-1.99272e-05; b3=9.47504e-08; b4=0.0; b5=0.0;
       }
       else {
         a0=0.74924; a1=-0.01233; a2=9.55715e-05; a3=-3.37117e-07; a4=4.67367e-10; a5=0.0;
         b0=0.25249; b1=-0.00161; b2=8.35478e-07; b3=1.25999e-08;  b4=0.0;         b5=0.0;
       }
     }
     else if(Mzamsa<9.5) {
       maxb=1.6;
       maxg=1.0;

       if(Ra>440.0) {
         lamb=1.136394;
         lamg=0.478963;
         manual=1;
       }
       else if(Ka==4 && Ra>180.0 && Ra<300.0) {
         lamb=0.15;
         lamg=0.08;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.49196; a1=-0.03247; a2=3.08066e-04; a3=-9.53247e-07; a4=0.0; a5=0.0;
         b0=0.805;   b1=-0.02;    b2=2.01872e-04; b3=-6.4295e-07;  b4=0.0; b5=0.0;
       }
       else if(Ka==4) {
         a0=0.39342; a1=0.00259;     a2=-4.97778e-05; a3=1.69533e-07; a4=0.0; a5=0.0;
         b0=0.20796; b1=6.62921e-04; b2=-1.84663e-05; b3=6.58983e-08; b4=0.0; b5=0.0;
       }
       else {
         a0=0.73147; a1=-0.01076; a2=7.54308e-05; a3=-2.4114e-07;  a4=2.95543e-10; a5=0.0;
         b0=0.31951; b1=-0.00392; b2=2.31815e-05; b3=-6.59418e-08; b4=7.99575e-11; b5=0.0;
       }
     }
     else if(Mzamsa<11.0) {
       maxb=1.6;
       maxg=1.0;

       if(Ra>500.0) {
         lamb=1.068300;
         lamg=0.424706;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         lamb=1.75*exp(-Ra/35.);
         lamg=0.9*exp(-Ra/35.);
         manual=1;
       }
       else if(Ka==4) {
         a0=0.75746; a1=-0.00852; a2=3.51646e-05; a3=-4.57725e-08; a4=0.0; a5=0.0;
         b0=0.35355; b1=-0.00388; b2=1.56573e-05; b3=-1.98173e-08; b4=0.0; b5=0.0;
       }
       else {
         a0=-9.26519; a1=0.08064;  a2=-2.30952e-04; a3=2.21986e-07; a4=0.0; a5=0.0;
         b0=0.81491;  b1=-0.00161; b2=-8.13352e-06; b3=1.95775e-08; b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<13.0) {
       maxb=1.6;
       maxg=1.0;

       if(Ra>600.0) {
         lamb=0.537155;
         lamg=0.211105;
         manual=1;
       }
       else if((Ka==4 && Ra>200.0 && Ra<410.0) || (Ka==5 && Ra>390.0 && Ra<460.0)) {
         lamb=0.08;
         lamg=0.05;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.63634; a1=-0.04646; a2=7.49351e-04; a3=-5.23622e-06; a4=0.0; a5=0.0;
         b0=1.17934; b1=-0.08481; b2=0.00329;     b3=-4.69096e-05; b4=0.0; b5=0.0;
       }
       else if(Ka==4) {
         a0=0.85249; a1=-0.00861; a2=2.99246e-05; a3=-3.21416e-08; a4=0.0; a5=0.0;
         b0=0.37188; b1=-0.00365; b2=1.24944e-05; b3=-1.32388e-08; b4=0.0; b5=0.0;
       }
       else {
         a0=-51.15252; a1=0.30238; a2=-5.95397e-04; a3=3.91798e-07; a4=0.0; a5=0.0;
         b0=-13.44;    b1=0.08141; b2=-1.641e-04;   b3=1.106e-07;   b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<15.0) {
       maxb=1.6;
       maxg=1.0;

       if(Ra>650.0) {
         lamb=0.3;
         lamg=0.160696;
         manual=1;
       }
       else if((Ka==4 && Ra>250.0 && Ra<490.0) || (Ka==5 && Ra>480.0 && Ra<540.0)) {
         lamb=0.06;
         lamg=0.05;
         manual=1;
       }
       else if(Ka==5 && Ra>540.0 && Ra<650.0) {
         lamb=Ra*1.8e-03-0.88;
         lamg=Ra*9.1e-04-0.43;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.45573; a1=-0.00937; a2=-0.00131; a3=3.07004e-05;  a4=0.0; a5=0.0;
         b0=1.19526; b1=-0.08503; b2=0.00324;  b3=-4.58919e-05; b4=0.0; b5=0.0;
       }
       else if(Ka==4) {
         a0=0.85271; a1=-0.00793; a2=2.5174e-05;  a3=-2.4456e-08;  a4=0.0; a5=0.0;
         b0=0.36163; b1=-0.00328; b2=1.03119e-05; b3=-9.92712e-09; b4=0.0; b5=0.0;
       }
       else {
         a0=-140.0;   a1=0.7126;  a2=-0.00121;     a3=6.846e-07;   a4=0.0; a5=0.0;
         b0=-44.1964; b1=0.22592; b2=-3.85124e-04; b3=2.19324e-07; b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<18.0) {
       maxb=1.5;
       maxg=1.0;

       if(Ra>750.0) {
         lamb=0.5;
         lamg=0.204092;
         manual=1;
       }
       else if((Ka==4 && Ra>200.0 && Ra<570.0) || (Ka==5 && Ra>560.0 && Ra<650.0)) {
         lamb=0.1;
         lamg=0.05;
         manual=1;
       }
       else if(Ka==5 && Ra>650.0 && Ra<750.0) {
         lamb=Ra*4.0e-03-2.5;
         lamg=Ra*1.5e-03-0.93;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.33378; a1=0.01274;  a2=-0.00234; a3=4.6036e-05;   a4=0.0; a5=0.0;
         b0=1.17731; b1=-0.07834; b2=0.00275;  b3=-3.58108e-05; b4=0.0; b5=0.0;
       }
       else if(Ka==4) {
         a0=0.83254; a1=-0.00696; a2=1.9597e-05;  a3=-1.67985e-08; a4=0.0; a5=0.0;
         b0=0.34196; b1=-0.0028;  b2=7.82865e-06; b3=-6.66684e-09; b4=0.0; b5=0.0;
       }
       else {
         a0=-358.4;     a1=1.599;   a2=-0.00238;    a3=1.178e-06;   a4=0.0; a5=0.0;
         b0=-118.13757; b1=0.52737; b2=-7.8479e-04; b3=3.89585e-07; b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<35.0) {
       maxb=1.5;
       maxg=1.0;

       if(Ra>900.0) {
         lamb=0.2;
         lamg=0.107914;
         manual=1;
       }
       else if((Ka==4 && Ra>230.0 && Ra<755.0) || (Ka==5 && Ra>725.0 && Ra<850.0)) {
         lamb=0.1;
         lamg=0.05;
         manual=1;
       }
       else if(Ka==5 && Ra>850.0 && Ra<900.0) {
         lamb=Ra*2.0e-03-1.6;
         lamg=Ra*1.0e-03-0.8;
         manual=1;
       }
       else if(Ka==2 || Ka==3) {
         a0=1.27138; a1=0.00538;  a2=-0.0012; a3=1.80776e-05;  a4=0.0; a5=0.0;
         b0=1.07496; b1=-0.05737; b2=0.00153; b3=-1.49005e-05; b4=0.0; b5=0.0;
       }
       else if(Ka==4) {
         a0=0.69746; a1=-0.0043;  a2=8.97312e-06; a3=-5.83402e-09; a4=0.0; a5=0.0;
         b0=0.26691; b1=-0.00161; b2=3.3378e-06;  b3=-2.1555e-09;  b4=0.0; b5=0.0;
       }
       else {
         a0=-436.00777; a1=1.41375; a2=-0.00153;     a3=5.47573e-07; a4=0.0; a5=0.0;
         b0=-144.53456; b1=0.46579; b2=-4.99197e-04; b3=1.78027e-07; b4=0.0; b5=0.0;
       }
     }
     else if(Mzamsa<75.0) {
       maxb=20.0;
       maxg=3.0;

       a0=0.821;   a1=-0.00669; a2=1.57665e-05; a3=-1.3427e-08;  a4=3.74204e-12; a5=0.0;
       b0=0.49287; b1=-0.00439; b2=1.06766e-05; b3=-9.22015e-09; b4=2.58926e-12; b5=0.0;
     }
     else {
       maxb=4.0;
       maxg=2.0;

       a0=1.25332; a1=-0.02065; a2=1.3107e-04;  a3=-3.67006e-07; a4=4.58792e-10; a5=-2.09069e-13;
       b0=0.81716; b1=-0.01436; b2=9.31143e-05; b3=-2.6539e-07;  b4=3.30773e-10; b5=-1.51207e-13;
     }
   }

   if(manual==0) {
     if(ZZ>Zlimit && (Mzamsa<1.5 || (Mzamsa>=18.0 && Mzamsa<25.0 && Ka>=5 && Ka<=6))) {
       x=Menva/Ma;
       y1=a0+a1*x+a2*x*x+a3*x*x*x+a4*x*x*x*x+a5*x*x*x*x*x;
       y2=b0+b1*x+b2*x*x+b3*x*x*x+b4*x*x*x*x+b5*x*x*x*x*x;
       lamb=1.0/y1;
       lamg=1.0/y2;
     }
     else if((ZZ>Zlimit && Mzamsa>=2.5 && Mzamsa<5.5 && Ka>=5 && Ka<=6) || (ZZ<=Zlimit && Mzamsa>=2.5 && Mzamsa<4.5 && Ka>=5 && Ka<=6)) {
       x=Ra;
       y1=a0+a1*x+a2*x*x+a3*x*x*x+a4*x*x*x*x+a5*x*x*x*x*x;
       y2=b0+b1*x+b2*x*x+b3*x*x*x+b4*x*x*x*x+b5*x*x*x*x*x;
       lamb=pow(10.0,y1);
       lamg=y2;
     }
     else {
       x=Ra;
       y1=a0+a1*x+a2*x*x+a3*x*x*x+a4*x*x*x*x+a5*x*x*x*x*x;
       y2=b0+b1*x+b2*x*x+b3*x*x*x+b4*x*x*x*x+b5*x*x*x*x*x;
       lamb=y1;
       lamg=y2;
     }
   }

   if(lamb<0.05)
     lamb=0.05;
   if(lamg<0.05)
     lamg=0.05;
   if(lamb>maxb)
     lamb=maxb;
   if(lamg>maxg)
     lamg=maxg;


   if(chlambda==1)                                                              /* the average value */
     lambda=(lamb+lamg)/2.0;
   else if(chlambda==2)                                                         /* lambda with full internal energy */
     lambda=lamb;
   else                                                                         /* lambda with no internal energy */
     lambda=lamg;
 }
 else if(chlambda==0 && Ka>=2 && Ka<=6)                                         /* H-eich stars: OLD APPROACH */
     lambda=Lambda;                                                             /* single value of lambda from sinbin.h */
 else
   fprintf(fp0,"error: unknown choice to calcuate CE lambda in lamf()\n");


 lambda*=golambda;                                      /* increase of lambda by fudge factor due to enthalpy arguments */

 return lambda;
}


double roche(double q, double a)
{ /* zwraca promien strefy Rocha dla M1, gdzie q=M2/M1 -Eggleton formula */
 double tmp=1/q;

 return (0.49*a)/(0.6+pow(tmp,-2.0/3.0)*log(1+pow(tmp,1.0/3.0)));
}


double Aroche(double q, double Rl)
{ /* returns separation for system with q=M2/M1 in which star M1 has roche lobe radius Rl */
 double tmp=1/q;

 return Rl*(0.6+pow(tmp,-2.0/3.0)*log(1+pow(tmp,1.0/3.0)))/0.49;
}


int explode(double Ma, double Maf, int Ka, double fraca, double Mb, double Rb, double *Vsm,
            double *Vsa, double *Vsb, double *Vextr, double *a, double *e, double *i, double *Om, double *om,
            double *tau, double *Vkick, int *cmt3, double *texp, int inbin, int ecssna, int fba,
            double *jxi, double *jyi, double *jzi, double *jxf, double *jyf, double *jzf)
{ /* wybucha gwiazda o masie Ma i pozostaje po niej obiekt zwarty o Maf i typie Ka, towarzysz byl i pozostaje Mb */
  /* orbita a,e,i,Om,om,tau ulega zmianie, nowe wartosci wpisane w te same zmienne */
  /* dodatkowa predkosc srodka masy ukladu dodawana do poprzedniej wartosci i zapisana w Vsm */
  /* if binary disrupted, velocities are recorded in Vsa, Vsb [Rsun/day] */
  /* all velocities in [Rsun/day] here */

 double X[3],V[3],v[3],tmp1,tmp2,P;
 double Xpre[3],Vpre[3],Vsmpre[3];
 double jx,jy,jz,J,Vtot;
 int is,snmerger,j;

 copyV(Vsm,Vsmpre);                    /* stores preSN velocity, no velocity gain assumed during SN merger */


 P=2.0*Pi*(*a)*sqrt((*a)/(G*(Ma+Mb)));                                               /* okres orbit. [dni] */
 (*texp)=get_texp(0.0,P);                                        /* czas wybuchu [dni], liczony od tau=0.0 */

 if(kick==1)
   get_Vkick1(Vkick);                                              /* kicks: -- composition of 2 Maxwlians */
 else if(kick==2)
   get_Vkick2(Vkick,fraca,fba);                          /* kicks: BH kicks decreased + no direct BH kicks */
 else if(kick==3)
   get_Vkick3(Vkick,Ka);                                              /* kicks: prescription2, no BH kicks */
 else if(kick==4)
   get_Vkick4(Vkick);                                                                    /* no NS,BH kicks */
 else if(kick==5)                                                   /* Paczynski kicks, BH kicks decreased */
   get_Vkick5(Vkick,Ka,fraca,cmt3,fba);
 else if(kick==6 || kick==7)            /* kick=6: Hobbs kicks: NS/BH kicks decreased + no direct BH kicks */
   get_Vkick6(Vkick,fraca,fba);                             /* kick=7: Hobbs kicks: full NS/BH natal kicks */
 else if(kick==8)                         /* kick=8: Hobbs kicks: full NS, and BH kicks lowered by BH mass */
   get_Vkick7(Vkick,Maf);                               /* fallback information ignored in this kick model */
 else if(kick==9)                                                    /* kick=9: Bray & Eldridge 2016 kicks */
   get_Vkick8(Vkick,Ma,Maf);
 else if(kick==10)                        /* kick=10: combined mass (fallback-decreased) and neutrino kick */
   get_Vkick9(Vkick,fraca,fba);


 if(ecssna==1)                     /* lower value of kicks by ECSlower for ECS SN */
   for(j=0;j<3;j++)
     Vkick[j]*=(ECSlower);
 if(ecssna==2)                     /* lower value of kicks by AIClower for ECS AIC */
   for(j=0;j<3;j++)
     Vkick[j]*=(AIClower);

 if(BINARY==0)                 /* for single star evolution just kick is needed */
   return 0;                   /* terminate here */


 if(inbin==1) {           /* 1st or 2nd SN, bound orbit at start */
   orbit1(*texp,Ma,Mb,*a,*e,*i,*Om,*om,*tau,X,v);
   copyV(X,Xpre);         /* position of exploding component before SN */
   copyV(v,Vpre);         /* velocity of exploding component before SN */

   *jxi=X[1]*v[2]-X[2]*v[1];
   *jyi=X[2]*v[0]-X[0]*v[2];
   *jzi=X[0]*v[1]-X[1]*v[0];
   J=sqrt((*jxi)*(*jxi)+(*jyi)*(*jyi)+(*jzi)*(*jzi));
   *jxi=(*jxi/J);
   *jyi=(*jyi/J);
   *jzi=(*jzi/J);

   if(polar==1) {
     jx=X[1]*v[2]-X[2]*v[1];
     jy=X[2]*v[0]-X[0]*v[2];
     jz=X[0]*v[1]-X[1]*v[0];
     J=sqrt(jx*jx+jy*jy+jz*jz);
     jx/=J;
     jy/=J;
     jz/=J;
     Vtot=sqrt(Vkick[0]*Vkick[0]+Vkick[1]*Vkick[1]+Vkick[2]*Vkick[2]);
     Vkick[0]=jx*Vtot;
     Vkick[1]=jy*Vtot;
     Vkick[2]=jz*Vtot;
   }
   V[0]=v[0]+Vkick[0];
   V[1]=v[1]+Vkick[1];
   V[2]=v[2]+Vkick[2];
   tmp1=Maf/(Maf+Mb);
   tmp2=(Mb*(Maf-Ma))/((Mb+Maf)*(Mb+Ma));

   is=orbit2(*texp,Maf,Mb,X,V,a,e,i,Om,om,tau);

   snmerger=0;
   if(is==1) {                          /* orbit bound */
     *jxf=X[1]*V[2]-X[2]*V[1];
     *jyf=X[2]*V[0]-X[0]*V[2];
     *jzf=X[0]*V[1]-X[1]*V[0];
     J=sqrt((*jxf)*(*jxf)+(*jyf)*(*jyf)+(*jzf)*(*jzf));
     *jxf=(*jxf/J);
     *jyf=(*jyf/J);
     *jzf=(*jzf/J);

     Vsm[0]=Vsm[0]+tmp1*Vkick[0]+tmp2*v[0];
     Vsm[1]=Vsm[1]+tmp1*Vkick[1]+tmp2*v[1];
     Vsm[2]=Vsm[2]+tmp1*Vkick[2]+tmp2*v[2];
     if(Rb>(*a)*(1.0-(*e)))
       snmerger=1;
   }
   else {                                        /* orbit disrupted (is=0) */
     snmerger=sinvel(Maf,Ma,Mb,Rb,Xpre,Vpre,Vkick,Vsm,Vsa,Vsb);
     Vsm[0]=Vsm[1]=Vsm[2]=0.0;                   /* no binary --> no center of mass velocity */

     *jxf=0.0;
     *jyf=0.0;
     *jzf=0.0;
   }
 }
 else {                                 /* exploding single star, may happen only at 2nd SN */
   Vextr[0]=Vsa[0]+Vkick[0];
   Vextr[1]=Vsa[1]+Vkick[1];
   Vextr[2]=Vsa[2]+Vkick[2];
   snmerger=0;                          /* mereger not possible: already disrupted binary */
 }

 if(snmerger==1) {             /* merger during SN. NS/BH shot into companion */
   copyV(Vsmpre,Vsa);
   copyV(Vsmpre,Vsb);
   is=-1;
 }

 return is;
}


int sinvel(double Ma, double Mapre, double Mb, double Rb, double *Xa, double *Va, double *wkick, double *Vsm, double *Vsa, double *Vsb)
{ /* calculates velocities of SN disrupted components: T.Bulik 2011, June 20 */
  /* and returns: 1- if a SN remnant (NS/BH) was kicked into the companion, 0- if not */
  /* A- SN star, B-companion, Ma=NS/BH mass, Mb=Mb (not changed) -- masses after SN explosion */
  /* Xa[3], Va[3] -- position and velocity of exploding star on relative orbit (prior to explosion) */
  /* wkick[3]: kick velocity imparted to star A, Vsm[3] -- center of mass velocity of pre-SN binary */
  /* fills in Vsa[3] -- velocity of single NS/BH, Vsb[3] -- velocity of single companion [Rsun/day] */
  /* units: velocities [Rsun/day], mass [Msun], distance: [Rsun]  */
 double V[3],Xpl[3],Vfpl[3],Vf[3],Vt1[3],Vt2[3],Vt3[3];
 double mi,s1,s2,c1,c2,sf,cf,E,alfa,eps,pe,jv[3],Jv2,Jvp,Rmin;
 double vfin,rzero,cosfizero,sinfizero,tmp,cosfifin,sinfifin,fizero;
 double Vapre[3],Vbpre[3],Vcm1[3],Vt[3];
 int i,snmerger;

 mi=Ma*Mb/(Mb+Ma);
 alfa=G*Ma*Mb;

 for(i=0;i<3;i++) {
   Vapre[i]=Mb*Va[i]/(Mapre+Mb);                           /* center of mass velocity of A: pre SN */
   Vbpre[i]=Vapre[i]-Va[i];                                /* center of mass velocity of B: pre SN */
   V[i]=Va[i]+wkick[i];                                    /* relative speed after SN explosion */
   Vcm1[i]=((Vapre[i]+wkick[i])*Ma+Vbpre[i]*Mb)/(Ma+Mb);   /* center of mass velocity after SN explosion */
 }

 jv[0]=mi*(Xa[1]*V[2]-Xa[2]*V[1]);                      /* angular momentum */
 jv[1]=mi*(Xa[2]*V[0]-Xa[0]*V[2]);
 jv[2]=mi*(Xa[0]*V[1]-Xa[1]*V[0]);

 Jv2=sqrt(jv[0]*jv[0]+jv[1]*jv[1]+jv[2]*jv[2]);
 Jvp=sqrt(jv[0]*jv[0]+jv[1]*jv[1]);
 if(Jvp>0.0000001) {                                    /* so no division by zero, and Jvp should be positive */
   s1=-jv[1]/Jvp;
   c1=jv[0]/Jvp;
 }
 else {
   s1=0.0;
   c1=1.0;
 }
 s2=Jvp/Jv2;
 c2=jv[2]/Jv2;

 rotate2for(Xa,Xpl,c1,s1,c2,s2);       /* rotate around z by c1, s1 and around y by c2 s2 to coordinate system where */
                                       /* vel. and dist. are confined to xy plane, and ang momentum is along z */
 rotate2for(V,Vt,c1,s1,c2,s2);
 rzero=sqrt(Xa[0]*Xa[0]+Xa[1]*Xa[1]+Xa[2]*Xa[2]);       /* inital distance (initial orbit) */
 E=0.5*mi*(V[0]*V[0]+V[1]*V[1]+V[2]*V[2])-alfa/rzero;   /* E - energy - must be positive for disrupted binary */
 vfin=sqrt(2.0*E/mi);                                   /* value of the final energy */
 eps=sqrt(1.0+2.0*E*Jv2*Jv2/(alfa*alfa*mi));            /* eps, pe - hyperbola parameters */
 pe=Jv2*Jv2/(alfa*mi);

 cosfizero=(pe/rzero-1.0)/eps;                          /* fizero -inital postion on the orbital hyperbola */

 if(cosfizero>0.99999) {                                /* to protect sinfizero become negative: */
   cosfizero=1.0;
   sinfizero=0.0;
 }
 else
   sinfizero=sqrt(1.0-cosfizero*cosfizero);             /* here if cosfizero is 1.0000000001! */
 tmp=Xa[0]*V[0]+Xa[1]*V[1]+Xa[2]*V[2];
 if(tmp<0.0)
   sinfizero=-sinfizero;
 snmerger=0;
                                                        /* after SN: NS approches on hiperbolic orbit its companion */
 Rmin=pe/(1.0+eps);                                     /* minimum distance between NS and companion */
 if(tmp<0.0) {
   if(Rmin<Rb)
     snmerger=1;                                        /* merger assumed */
 }

 cosfifin=-1.0/eps;                                     /* fifin -final poistion on the hyperbola, when r-> infty */
 sinfifin=sqrt(1.0-cosfifin*cosfifin);
 cf=cosfifin*cosfizero+sinfifin*sinfizero;              /* cf, sf - sine and cosine of the angle fifin-fizero */
 sf=sinfifin*cosfizero-cosfifin*sinfizero;
 for(i=0;i<3;i++)                                       /* normalized vector to the initial direction */
   Xpl[i]=Xpl[i]/rzero;

 Vfpl[0]=vfin*(Xpl[0]*cf-Xpl[1]*sf);                    /* final velocity - value and direction */
 Vfpl[1]=vfin*(+Xpl[0]*sf+Xpl[1]*cf);
 Vfpl[2]=0.0;

 rotate2back(Vfpl,Vf,c1,s1,c2,s2);                      /* Rotate back */

 for(i=0;i<3;i++) {                                     /* transform back to the initial coordinate system */
   Vsa[i]=Vsm[i]+Vcm1[i]+mi*Vf[i]/Ma;                          /* and here are your final velocities v1f and v2f */
   Vsb[i]=Vsm[i]+Vcm1[i]-mi*Vf[i]/Mb;
 }

 return snmerger;
}


void rotate2for(double *Vi, double *Vf, double c1, double s1, double c2, double s2)
{ /* rotate 3-component vector Vi from relative orbit coordinates to: */
  /* Tauris&Takens cordinates and fills in values to Vf */
 double m[3][3];
 int i;

 m[0][0]=c2*c1;
 m[1][0]=-c2*s1;
 m[2][0]=-s2;

 m[0][1]=s1;
 m[1][1]=c1;
 m[2][1]=0.0;

 m[0][2]=s2*c1;
 m[1][2]=-s2*s1;
 m[2][2]=c2;

 for(i=0;i<3;i++)
   Vf[i]=m[0][i]*Vi[0]+m[1][i]*Vi[1]+m[2][i]*Vi[2];
}


void rotate2back(double *Vi, double *Vf, double c1, double s1, double c2, double s2)
{ /* rotate 3-component vector Vi from Tauris&Takens cordinates back to: */
  /* relative orbit coordinates (fills in Vf): rotation by 3 angles: zeta1, zeta2, zeta3 */
 double m[3][3];
 int i;

 m[0][0]=c1*c2;
 m[1][0]=s1;
 m[2][0]=c1*s2;

 m[0][1]=-c2*s1;
 m[1][1]=c1;
 m[2][1]=-s1*s2;

 m[0][2]=-s2;
 m[1][2]=0.0;
 m[2][2]=c2;

 for(i=0;i<3;i++)
   Vf[i]=m[0][i]*Vi[0]+m[1][i]*Vi[1]+m[2][i]*Vi[2];
}


double tmerge(double Ma, double Mb, double a, double e)
{ /* calculates exact merging time of binary due to emission of grawitational waves */
  /* If e>0.9999 qtrap() may be unable to integrate given system, so I set Tmerge=-1.0 */
  /* as for idum=-194953, iidd_old=1972851 where e=0.999993933701 => [3] floating point exception */
  /* also if qtrap() returns -1, which means that it is unable to integrate system in JMAX3 steps */
  /* I set Tmerge=-1.0 */
  /* using equations of Peters P.C. 1964, Phys.Rev. 136, B1224 */
  /* Ma,Mb [M_sun], a [R_sun] --> Tmr [10^6 yrs] */
  /* needs functions func_tmr(), qtrap(), trapzd() */
 double tm,tmp,c0,beta0;
 double c=2.99792458e+10;                                         /* speed of light [cm/s] */
                                                       /* tr4 below change it to [Rsun/day] */

 if(e>0.9999)                 /* qtrap() may be unable to integrate this system with higher e */
   e=0.9999;                 /* so calculate the lower limit on T_merger */

 beta0=64.0*G*G*G*Ma*Mb*(Ma+Mb)/(5.0*pow(tr4*c,5.0));                        /* after eq. 5.9 */

 if(fabs(e)<acc) {                                                             /* circular orbit, eq. 5.10 */
   tm=a*a*a*a/(4.0*beta0);
   tm/=(365.25*(1.0e+06));                                                   /* now in [10^6yrs] */
 }
 else {                                                                      /* eccentric orbit */
   c0=a*(1.0-e*e)/(pow(e,12.0/19.0)*pow(1.0+121.0*e*e/304.0,870.0/2299.0));  /* from eq. 5.11 */
   tmp=qtrap(func_tmr,0.0,e);                                              /* integral from eq. 5.11 */
   if(fabs(tmp+1.0)<acc) {
     tm=-1.0;                                              /* qtrap() unable to integrate this system */
   }
   else {
     tm=12.0*c0*c0*c0*c0*tmp/(19.0*beta0);                                     /* whole eq. 5.11 */
     tm/=(365.25*(1.0e+06));                                                   /* now in [10^6yrs] */
   }
 }

 return tm;
}


void symbiotic(double Ma, double dMwinda, double Mca, double La, double Ra, double Rb,
               double t, double Mb, double a, double e, int Ka, int Kb, double dMmta, double dMmtb,
               double dMtran, int mttype, double Mzamsa, double Mzamsb, double a0, double e0)
{ /* symbiotic classification + output for symbiotics */
  /* A: giant donor, B: compact accretor */
 double Tsymb,dMsymb,Psymb,dMwind,dMedd,Rla;
 double beta_wind,alfa_wind,v_wind2,v_orbi2,v_2,dMa,Ta,p;
 char his[100];
 int kk;

 Tsymb=5000.0;                                      /* corresponds to types later than mid-G for stars III-I */
 dMsymb=1.0e-09;                                /* [Msun/yr]  critical lower acc. limit to produce symbiotic */
 Psymb=10.0;                                                        /* [yr]  maximum symbiotic binary period */
 p=2.0*Pi*a*sqrt(a/(G*(Ma+Mb)));                                                      /* binary period [day] */
 p/=365.25;                                                                                     /* day -> yr */
 Ta=get_T(La,Ra);                                                              /* effective temperature of A */
 Rla=roche(Mb/Ma,a*(1.0-e));

 if((Ka==3 || Ka==5 || Ka==6) && ((Kb>=10 && Kb<=14) || Kb==16 || Kb==17) && p<Psymb && Ta<=Tsymb) {

   dMedd=dMeddf(Kb,Rb,Ka);
   dMa=dMwinda/1.0e+06;                                                               /* change to [Msun/yr] */
   beta_wind=0.125;                                                      /* parametrization of wind velocity */
   alfa_wind=1.5;                                                         /* Bondi-Hoyle spherical accretion */
   v_wind2=2.0*beta_wind*GG*Ma/Ra;                                             /* 2nd power of wind velocity */
   v_orbi2=GG*(Ma+Mb)/a;                                                       /* 2nd power of B orbital velocity */
   v_2=v_orbi2/v_wind2;
   dMwind=(pow(GG*Mb/v_wind2,2.0)*alfa_wind*dMa)/(sqrt(1.0-e*e)*2.0*a*a*pow(1.0+v_2,1.5));
   dMwind=min(0.5*dMa,SymdM*dMwind);                                          /* wind accre. on Mb, enhanced */
                                                     /* hard to accrete more than 50% what was lost by comp. */
   if(dMmta>acc) {                                                                         /* RLOF symbiotic */
     fprintf(fp108,"1 %.3f %.1f %d %d %.3f %.3f %.3f %.3f %.3f %.3f %.1e %.1e %.0f %.3f %.1e %.1e %.1e %.1e %d %.1f %.1f %.1f %.2f %d         %d %s\n",
                      t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,p,Rla,Ra,La,Ta,Mca,dMmta*1.0e-06,dMmtb*1.0e-06,
                      dMedd*1.0e-06,dMtran*1.0e-06,mttype,Mzamsa,Mzamsb,a0,e0,iidd_old,nevroute,evroute);
   }
   else if(dMwind>dMsymb) {                                                                /* WIND symbiotic */
     fprintf(fp108,"2 %.3f %.1f %d %d %.3f %.3f %.3f %.3f %.3f %.3f %.1e %.1e %.0f %.3f %.1e %.1e %.1e %.1e %d %.1f %.1f %.1f %.2f %d         %d %s\n",
                      t+Tst,Tst,Ka,Kb,Ma,Mb,a,e,p,Rla,Ra,La,Ta,Mca,dMa,dMwind,
                      dMedd*1.0e-06,dMtran*1.0e-06,mttype,Mzamsa,Mzamsb,a0,e0,iidd_old,nevroute,evroute);
   }
   fflush(fp108);   /* output: period [yr], dM [Msun/yr] */
 }
}


void Xbin(double *Lxmt, double dMtran, double Ma, double dMwinda, double Mca, double dMmta, double La, double Ra,
          double t, double dt, double Mb, double Rb, double a, double e, int Ka, int Kb, double *Vsm, double tsn,
          int *mark10, int *mark11, int mttype, int flagbra, int flagbrb, int mtflag1a, int mtflag1b, int mtflag1ab,
          double Maio, double Mbio, double dMcea, double dMceb, double Mzamsa, double Mzamsb, double a0, double e0,
          int ecssna, int ecssnb, double aspinb)
{ /* Xray_binary: class+output (Jan 9, 2018) */
  /* A: donor, B: compact accretor */
  /* accretion units here are [Msun/Myr], only output changed to [Msun/yr] !!! */
 double p,dMa,dMwind,dMedd,Lx,Lxcrit,Ledd,dMdonor,Ys,etaA;
 double beta_wind,alfa_wind,v_wind2,v_orbi2,v_2,Rla,Rlb;
 double Lxdonor,dMcapt,Lxcapt,Racc,Tmr;
 double a1,b1,u1,u2,w1,w2;
 int tran,i,kk;
 char his[100];

 Ys=31557600;                                                                             /* year in seconds */
 dMdonor=dMmta;                                                 /* donor mass loss rate which feeds the disk */
                                                                     /* to be compared to dMtran! [Msun/Myr] */
 Rla=roche(Mb/Ma,a*(1.0-e));
 Rlb=roche(Ma/Mb,a*(1.0-e));

 dMedd=dMeddf(Kb,Rb,Ka);
 Tmr=tmerge(Ma,Mb,a,e);                /* formal merger time [Myr] for two point masses Ma,Mb at orbit (a,e) */

 if(Kb==13) {
   etaA=1.0;                                                                    /* surface accretion onto NS */
   Racc=Rb;
 }
 else if(Kb==14) {
   etaA=0.5;                                                         /* disk accretion onto BH */
   Racc=3.0*Rb;                                       /* disk ends at 3 Swartzhild radii of BH */
 }
 Ledd=etaA*GG*Mb*(dMedd/1.0e+06)/Racc;
 Ledd*=(Rsun*Rsun*Msun/pow(Ys,3.0));                           /* Eddington Luminosity for the NS/BH [erg/s] */

 Lxcrit=-1.0e+50;                       /* e.g.: 1.0e+27 [erg/s] critical Lx, over which code record sources */
                          /* with large negative value it records every potential source, both WIND and RLOF */
 p=2.0*Pi*a*sqrt(a/(G*(Ma+Mb)));                                                      /* binary period [day] */

 (*Lxmt)=0.0;
 if(dMmta>acc) {                                                                          /* MT X-ray system */

   Lxdonor=etaA*GG*Mb*(dMmta/1.0e+06)/Racc;                       /* not Eddington limitted X-ray luminosity */
   Lxdonor*=(Rsun*Rsun*Msun/pow(Ys,3.0));

   dMmta=min(dMedd,dMmta);                                                 /* Eddington limitted MT accretion */
   (*Lxmt)=etaA*GG*Mb*(dMmta/1.0e+06)/Racc;                 /* [Rsun^2 Msun/yr^3] Hurley et al. 2002, eq. 105 */
   (*Lxmt)*=(Rsun*Rsun*Msun/pow(Ys,3.0));                         /* [erg/s], W=J/s=m^2 kg /s^3, 1J=10^7 ergs */

   tran=0;
   if(dMdonor<=dMtran && PT==1) {                                                         /* transient source */
     tran=1;
     (*Lxmt)=Ledd;                 /* for transient source we assume that it emmits with Eddington Luminosity */
   }                                  /* decision if the source is on or off, must be done in post-processing */

   if((*Lxmt)>=Lxcrit) {                                                              /* record Lx MT evolution */
     i=(int)floor(t+Tst);                                                         /* time coordinate t, dt=1Myr */

     if((Xdtstep==1 && i!=(*mark11) && Xout1==1) || (Xdtstep==0 && Xout1==1)) {
       (*mark11)=i;                                                                /* block the time t=i+1Myr   */
       fprintf(fp131,"%f %f %.3f %.3f %d %d %.3f %.3f %f %f %.1e %.1e %.1e %.1e %.1e %.1e %d %d %d %d %d %d %d %.1f %.1f %.1f %.1f %.1f %.2f %.2f %.2f %.2f %.2f %.2f %.1f %.2f %.1e %d %d %.4f %d %d  %d %s\n",
                      t+Tst,dt,Ma,Mb,Ka,Kb,a,e,Ra,Rla,
                      dMedd,dMtran,dMmta,*Lxmt,dMdonor,Lxdonor,mttype,tran,
                      flagbra,flagbrb,mtflag1a,mtflag1b,mtflag1ab,Tst,tsn+Tst,Vsm[0],Vsm[1],Vsm[2],
                      Maio,Mbio,dMcea,dMceb,Mzamsa,Mzamsb,a0,e0,Tmr,ecssna,ecssnb,aspinb,idum_run,iidd_old,nevroute,evroute);
       fflush(fp131);
     }
   }

 }
 else {                                                                   /* check for WIND-fed X-ray system */

   if(Ka>=2 && Ka<=6 && Ra>=900.0)                    /* parametrization of wind velocity, after Jarrod binary paper */
     beta_wind=0.125;                                 /* small speeds for big (cool) supergiants */
   else if(Ka==0 || Ka==1) {                          /* for MS stars drops from 7->0.5 from O to F spc. type */
     if(Ma>=120.0) beta_wind=7.0;
     else if(Ma<=1.4) beta_wind=0.5;
     else {
       u1=120.0; w1=7.0;
       u2=1.4;   w2=0.5;
       a1=(w1-w2)/(u1-u2);
       b1=w1-a1*u1;
       beta_wind=a1*Ma+b1;
     }
   }
   else {                                                 /* for other stars drops from 7->0.5 from O to F spc. type */
     if(Ma>=140.0) beta_wind=7.0;                         /* masses for O,F spc. types in interpolation for supergiants */
     else if(Ma<=10.0) beta_wind=0.5;
     else {
       u1=140.0; w1=7.0;
       u2=10.0;   w2=0.5;
       a1=(w1-w2)/(u1-u2);
       b1=w1-a1*u1;
       beta_wind=a1*Ma+b1;
     }
   }

   dMa=dMwinda;                                             /* wind dM/dt for A [Msun/Myr] (positive number) */
   alfa_wind=1.5;                                                           /* Bondi-Hoyle spherical accretion */
   v_wind2=2.0*beta_wind*GG*Ma/Ra;                                               /* 2nd power of wind velocity */
   v_orbi2=GG*(Ma+Mb)/a;                                                         /* 2nd power of B orbital velocity */
   v_2=v_orbi2/v_wind2;
   dMwind=(pow(GG*Mb/v_wind2,2.0)*alfa_wind*dMa)/(sqrt(1.0-e*e)*2.0*a*a*pow(1.0+v_2,1.5));
   dMwind=min(0.8*dMa,dMwind);                              /* orbit avaraged wind accretion onto B [Msun/Myr] */
   dMcapt=dMwind;
   dMwind=min(dMedd,dMwind);                  /* to protect accretion luminosity going over the Eddington rate */

   Lxcapt=etaA*GG*Mb*(dMcapt/1.0e+06)/Racc;                         /* not Eddington limitted X-ray luminosity */
   Lxcapt*=(Rsun*Rsun*Msun/pow(Ys,3.0));

   Lx=etaA*GG*Mb*(dMwind/1.0e+06)/Racc;                      /* [Rsun^2 Msun/yr^3] Hurley et al. 2002, eq. 105 */
   Lx*=(Rsun*Rsun*Msun/pow(Ys,3.0));                               /* [erg/s], W=J/s=m^2 kg /s^3, 1J=10^7 ergs */

   if(Lx>=Lxcrit) {                      /* record Lx wind-fed evolution */
     i=(int)floor(t+Tst);                /* time coordinate t, dt=1Myr */

     if((Xdtstep==1 && i!=(*mark10) && Xout2==1) || (Xdtstep==0 && Xout2==1)) {
       (*mark10)=i;                      /* block the time t=i+1Myr   */
       fprintf(fp130,"%f %f %.3f %.3f %d %d %.3f %.3f %f %f %.1e %.1e %.1e %.1e %.1e %.1e -10 -10 %d %d %d %d %d %.1f %.1f %.1f %.1f %.1f %.2f %.2f %.2f %.2f %.2f %.2f %.1f %.2f %.1e %d %d %.4f %d %d  %d %s\n",
                      t+Tst,dt,Ma,Mb,Ka,Kb,a,e,Ra,Rla,
                      dMedd,dMtran,dMwind,Lx,dMcapt,Lxcapt,flagbra,flagbrb,
                      mtflag1a,mtflag1b,mtflag1ab,Tst,tsn+Tst,Vsm[0],Vsm[1],Vsm[2],Maio,Mbio,dMcea,dMceb,
                      Mzamsa,Mzamsb,a0,e0,Tmr,ecssna,ecssnb,aspinb,idum_run,iidd_old,nevroute,evroute);
       fflush(fp130);
     }
   }
 }
}


void Xbin1(double t, double dt, double Ma, double Mb, int Ka, int Kb, double Ra, double dMmta, double dMwinda, double aspinb, double a, double e, double La)
{ /* Xray_binary: class+output (Oct 29, 2019): calculates either RLOF ***OR*** WIND luminosity (not the two togehter) */
  /* A: donor, B: compact accretor (only for NS and BH) */
  /* all accretion/mass loss units are [Msun/Myr] on input and output of this function */
  /* dMmta [Msun/Myr]: RLOF mass loss (positive number) */
  /* dMwinda [Msun/Myr]: wind mass loss (positive number) */
 double dMedd1,Ledd1;
 double f1,dM0rlof,Lxiso;
 double etaA,Rns,Risco;
 double b,dm0,Lxbeam;
 double fwind,dMcap,dM0wind;
 double Ta,Rla,Gbetaa,dper;


 if(dMmta>acc) {                  /* RLOF X-ray system */

   dMedd1=dMeddf1(Mb,Ka);         /* [Msun/Myr]: Eddington accretion critical rate; positive */
   Ledd1=Leddf1(Mb,Ka);           /* [erg/s]: Eddington critical luminosity */

   f1=1.0;
   dM0rlof=f1*dMmta;              /* accretion rate at spherizaion radius (Rsph) [Msun/Myr]; positive */

   if(dM0rlof>dMedd1)
     Lxiso=Ledd1*(1.0+log(dM0rlof/dMedd1));    /* isotropic X-ray luminosity [erg/s] */
   else {
     if(Kb==13) {
       Rns=Rnsf1(Mb);                         /* NS stars radius [Rsun]: assumed to be NS accretion radius */
       etaA=(GGG*Mb)/(CCC*CCC*Rns);           /* efficiency of gravitational energy release for NS [unitless] */
     }
     else if(Kb==14) {
       Risco=Riscof(Mb,aspinb);               /* ISCO radius [Rsun]: assumed to be BH accretion radius */
       etaA=1.0-Eiscof(Risco,Mb,aspinb);      /* efficiency of gravitational energy release for BH [unitless] */
     }
     Lxiso=etaA*dM0rlof*CCC*CCC;              /* isotropic X-ray luminosity [Rsun^2 Msun Myr^-3] */
     Lxiso=tr1*Lxiso;                         /* [erg s^-1] */
   }

   dm0=dM0rlof/dMedd1;                        /* mass transfer rate at Rsph in Eddington transer rate units [unitless] */
   if(dm0>=150.0)
     b=3.2e-03;
   else if(dm0>=8.5)
     b=73.0/(dm0*dm0);
   else
     b=1.0;
   b=min(b,1.0);
   Lxbeam=Lxiso/b;                            /* beamed X-ray luminosity [erg s^-1] */

   if(Xout1==1 && Lxbeam>Lxcrit1) {
     fprintf(fp131,"%f %f %.3f %.3f %d %d %.3f %.3f %g %g %g %d %d  %d %s\n",
             t,dt,Ma,Mb,Ka,Kb,a,e,dMmta,Lxiso,Lxbeam,idum_run,iidd_old,nevroute,evroute);
     fflush(fp131);
   }
 }


 else if(dMwinda>acc) {             /* WIND */

   fwind=Fbetaf(Ma,Ra,Mb,a,e,Ka);   /* capture fraction of A star wind into disk around compact object B */
   if(WindRLOF==1 && Ma<8.0 && Ka>=3 && Ka<=6) {   /* valid only for low to intermediate--mass giants: Abate et al 2013 */
      dper=a*(1.0-e);                                               /* separation of stars at periastron */
      Rla=roche(Mb/Ma,dper);
      Ta=get_T(La,Ra);
      Gbetaa=Gbetaf(Ma,Mb,Ta,Ra,Rla);               /* fraction of wind from Ma, accreted onto Mb */
      fwind=max(fwind,Gbetaa);
   }

   dMcap=fwind*dMwinda;             /* orbit avaraged wind accretion into disk around B [Msun/Myr] */


   dMedd1=dMeddf1(Mb,Ka);         /* [Msun/Myr]: Eddington accretion critical rate; positive */
   Ledd1=Leddf1(Mb,Ka);           /* [erg/s]: Eddington critical luminosity */

   f1=1.0;
   dM0wind=f1*dMcap;              /* accretion rate at spherizaion radius (Rsph) [Msun/Myr]; positive */

   if(dM0wind>dMedd1)
     Lxiso=Ledd1*(1.0+log(dM0wind/dMedd1));    /* isotropic X-ray luminosity [erg/s] */
   else {
     if(Kb==13) {
       Rns=Rnsf1(Mb);                         /* NS stars radius [Rsun]: assumed to be NS accretion radius */
       etaA=(GGG*Mb)/(CCC*CCC*Rns);           /* efficiency of gravitational energy release for NS [unitless] */
     }
     else if(Kb==14) {
       Risco=Riscof(Mb,aspinb);               /* ISCO radius [Rsun]: assumed to be BH accretion radius */
       etaA=1.0-Eiscof(Risco,Mb,aspinb);      /* efficiency of gravitational energy release for BH [unitless] */
     }
     Lxiso=etaA*dM0wind*CCC*CCC;              /* isotropic X-ray luminosity [[Rsun^2 Msun Myr^-3]] */
     Lxiso=tr1*Lxiso;                         /* [erg s^-1] */
   }

   dm0=dM0wind/dMedd1;                        /* mass transfer rate at Rsph in Eddington transer rate units [unitless] */
   if(dm0>=150.0)
     b=3.2e-03;
   else if(dm0>=8.5)
     b=73.0/(dm0*dm0);
   else
     b=1.0;
   b=min(b,1.0);
   Lxbeam=Lxiso/b;                            /* beamed X-ray luminosity [erg/s] */


   if(Xout2==1 && Lxbeam>Lxcrit2) {
     fprintf(fp130,"%f %f %.3f %.3f %d %d %.3f %.3f %g %g %g %g %d %d  %d %s\n",
             t,dt,Ma,Mb,Ka,Kb,a,e,dMwinda,dMcap,Lxiso,Lxbeam,idum_run,iidd_old,nevroute,evroute);
     fflush(fp130);
   }
 }

}


double Fbetaf1(double dMwinda, double Ma, double Ra, double Rb, double Mb, double a, double e, int Ka, int Kb, double aspinb)
{ /* calculates fraction of wind mass loss rate that is captured on the companion */
  /* A: donor, B: wind accretor,  accretion units here are [Msun/Myr] */
  /* July 31, 2018; new accretion scheme SAMARESH+JP */
 double beta_wind,alfa_wind,v_wind2,v_orbi2,v_2;
 double a1,b1,u1,u2,w1,w2;
 double Fbeta,dMcap,dM0wind,dMedd,dMedd1,fwind;
 double f1,f2,Rs,Rsph,Racc,dMgain;


 fwind=Fbetaf(Ma,Ra,Mb,a,e,Ka);   /* capture fraction of A star wind into disk around compact object B */
 dMcap=fwind*dMwinda;             /* orbit avaraged wind accretion into disk around B [Msun/Myr] */

 if(Kb==13 || Kb==14) {
   dMedd1=dMeddf1(Mb,Ka);                             /* Eddington limit [Msun/Myr]; positive */
   f1=1.0;                                            /* (1-f1): wind mass loss from outer disk */
   dM0wind=f1*dMcap;                                  /* accretion at Rsph [Msun/Myr]; positive */
   Rs=(2.0*GGG*Mb)/(CCC*CCC);                         /* Schwarzschild radius for NS/BH [Rsun] */
   Rsph=(27.0*dM0wind*Rs)/(4.0*dMedd1);               /* spherization radius for NS/BH disk [Rsun */
   if(Kb==13)                                         /* accretion radius for NS or BH */
     Racc=Rnsf1(Mb);                                  /* NS radius: surface accretion */
   else
     Racc=Riscof(Mb,aspinb);                          /* ISCO radius: truncated disk for BH */
   if(dM0wind<=dMedd1)                                /* (1-f2): wind mass loss from inner disk */
     f2=1.0;
   else
     f2=Racc/Rsph;
   dMgain=-f1*f2*dM0wind;                               /* accumulation rate on NS/BH [Msun/Myr]; negative sign */
   Fbeta=fabs(dMgain/dMwinda);
 }
 else if(Kb==10 || Kb==11 || Kb==12 || Kb==16 || Kb==17) {  /* Edd. limit imposed for WD accret. */
   dMedd=dMeddf(Kb,Rb,Ka);                                  /* Eddington limit [Msun/Myr]; positive */
   if(fabs(dMcap)<=dMedd) {
     dMgain=dMcap;
     Fbeta=fabs(dMgain/dMwinda);
   }
   else {
     dMgain=-dMedd;                                         /* WD accretes maximum at Eddington limit */
     Fbeta=fabs(dMedd/dMwinda);
   }
 }
 else {                                    /* any other non-compact accretor; catches all it can from wind */
   Fbeta=fabs(dMcap/dMwinda);
 }

 return Fbeta;
}


double Fbetaf(double Ma, double Ra, double Mb, double a, double e, int Ka)
{ /* calculates fraction of wind mass loss rate that is captured on the companion */
  /* A: donor, B: wind accretor,  accretion units here are [Msun/Myr] */
  /* corrected on Dec 12 2018: beta=0.125 for some stars was producing sharp jump in accretion from winds */
  /* it was not physical: simple correction (to be checked and discussed with wind experts) */
 double beta_wind,alfa_wind,v_wind2,v_orbi2,v_2;
 double a1,b1,u1,u2,w1,w2;
 double Fbeta,dMcap;


 if(Ka==0 || Ka==1) {                               /* for MS stars drops from 7->0.5 from O to F spc. type */
   if(Ma>=120.0) beta_wind=7.0;
   else if(Ma<=1.4) beta_wind=0.5;
   else {
     u1=120.0; w1=7.0;
     u2=1.4;   w2=0.5;
     a1=(w1-w2)/(u1-u2);
     b1=w1-a1*u1;
     beta_wind=a1*Ma+b1;
   }
 }
 else {                                                 /* for other stars drops from 7->0.5 from O to F spc. type */
   if(Ma>=140.0) beta_wind=7.0;                         /* masses for O,F spc. types in interpolation for supergiants */
   else if(Ma<=10.0) beta_wind=0.125;
   else {
     u1=140.0; w1=7.0;
     u2=10.0;   w2=0.125;
     a1=(w1-w2)/(u1-u2);
     b1=w1-a1*u1;
     beta_wind=a1*Ma+b1;
   }
 }

 alfa_wind=1.5;                                                           /* Bondi-Hoyle spherical accretion */
 v_wind2=2.0*beta_wind*GG*Ma/Ra;                                               /* 2nd power of wind velocity */
 v_orbi2=GG*(Ma+Mb)/a;                                                         /* 2nd power of B orbital velocity */
 v_2=v_orbi2/v_wind2;
 Fbeta=(pow(GG*Mb/v_wind2,2.0)*alfa_wind)/(sqrt(1.0-e*e)*2.0*a*a*pow(1.0+v_2,1.5));
 if(Fbeta>0.8) Fbeta=0.8;                       /* to protect accreting over 80% for highly eccentric cases */

 return Fbeta;
}


double Gbetaf(double Ma,double Mb,double Ta,double Ra,double Rla)
{ /* Wind Roche Lobe Overflow: gravitational lensing: returns accretion efficiency (beta) */
  /* Ma-donor mass [Msun]; Mb-accretor mass [Msun]; Ta-donor effective temperature [K]; */
  /* Rd-donor radius [Rsun]; Rld - donor Roche lobe radius [Rsun]*/
  /* Based on Abate et al., 2013, 552, 26 */
 double x,beta;

 x=(0.5*Ra*pow(Ta/1500.,2.5))/Rla;  /* (Dust formation radius/Roche Lobe Radius) calculated Eq.1 in Abate */
                                    /* Höfner (2007); assumed dust condensation temperature = 1500 K  */

 beta=(25.0/9.0)*pow(Mb/Ma,2.0)*(-0.284*x*x+0.918*x-0.234);                   /* Eq. 9 in Abate et al.*/
 if(beta>0.5)                                                                 /*Eq. 9 continued*/
   beta=0.5;
 else if(beta<0.)
   beta=0.0;

 return beta;
}


void Xwdwd(double *Lxwd, double dMmta, double Mb, double Rb, int Ka, int Kb)
{ /* WD Xray binary (Nov 16, 2004): A: WD donor, B: WD accretor */
  /* fills in Lxwd -- X-ray luminosity [erg/s] */
 double dMedd,etaA,Racc,Ys;

 Ys=31557600;                                                                             /* year in seconds */
 dMedd=dMeddf(Kb,Rb,Ka);
 etaA=1.0;                                /* surface accretion onto WD */
 Racc=Rb;

 (*Lxwd)=0.0;
 if(dMmta>acc) {                                                                          /* MT X-ray system */
   dMmta=min(dMedd,dMmta);                                                 /* Eddington limitted MT accretion */
   (*Lxwd)=etaA*GG*Mb*(dMmta/1.0e+06)/Racc;                 /* [Rsun^2 Msun/yr^3] Hurley et al. 2002, eq. 105 */
   (*Lxwd)*=(Rsun*Rsun*Msun/pow(Ys,3.0));                         /* [erg/s], W=J/s=m^2 kg /s^3, 1J=10^7 ergs */
 }
}


double Xcvf(int Ka, int Kb, double Mb, double Rb, double dMmta, double *dMacc, int *cvtype)
{ /* returns X-ray Luminosity (Lx) of a CV/IP system [ergs/sec] and fills in cvtype */
  /* it needs to be expanded to work better for CVs */
  /* cvtype=0 -- not a CV, cvtype=Ys=31557600;1 -- IP, cvtype=2 -- other type of CV */
  /* A is donor, B is accretor */
 double dMedd,Lacc,Lx,tmp1;


 if((Kb==10 || Kb==11 || Kb==12 || Kb==16 || Kb==17) && dMmta>acc)  {                  /* CV/IP: WD accretor */
   tmp1=get_flat(0.0,1.0);
   if(tmp1<=IPfrac)
     (*cvtype)=1;                                                         /* system is an IP: Lx calculation */
   else
     (*cvtype)=2;                                               /* system is a CV but not IP: Lx calculation */

   if((*cvtype)==1) {                                                           /* IP: accreting magnetic WD */
     dMedd=dMeddf(Kb,Rb,Ka);                             /* critical Eddington mass transfer rate [Msun/Myr] */
     (*dMacc)=min(dMmta,dMedd);                /* Eddington limited mass accretion rate onto a WD [Msun/Myr] */
     (*dMacc)=etaIP*(*dMacc);                         /* IP limited mass accretion rate onto a WD [Msun/Myr] */
     Lacc=GGG*Mb*(*dMacc)/Rb;                         /* bolometric accretion luminosity [Rsun^2 Msun/Myr^3] */
     Lacc*=tr1;                                                                                   /* [erg/s] */
     Lx=etaBE*etaGEO*etaXIP*Lacc;                                     /* IP X-ray Chandra luminosity [erg/s] */
   }
   else if((*cvtype)==2) {                                           /* regular CV (not an IP): accreting WD */
     dMedd=dMeddf(Kb,Rb,Ka);                             /* critical Eddington mass transfer rate [Msun/Myr] */
     (*dMacc)=min(dMmta,dMedd);                /* Eddington limited mass accretion rate onto a WD [Msun/Myr] */
     Lacc=GGG*Mb*(*dMacc)/Rb;                         /* bolometric accretion luminosity [Rsun^2 Msun/Myr^3] */
     Lacc*=tr1;                                                                                   /* [erg/s] */
     Lx=etaXCV*Lacc;                                                          /* CV X-ray luminosity [erg/s] */
   }
 }
 else {
   (*cvtype)=0;
   Lx=0.0;
   fprintf(fp0,"error: in Xcvf() unknown accretor (Kb: %d) and/or no ongoing RLOF (dMmta: %g)\n",Kb,dMmta);
 }

 return Lx;                                                /* output goes to Lxcv in cvcut.dat and allcv.dat */
}


void spin_evol(int Ka, double *Ma, double *aspin, double *Mdisa, double Mrest)
{ /* A-accreting BH, changes BH spin and mass due to accretion (Mrest=-dMmta*dt) */
  /* Mrest: rest mass (positive) of transferred material [Msun] */
  /* formalism from: Brown et al. 2000, New Astronomy, 5, 191 */
  /* aspin=Jc/(M^2G) [unitless]: BH spin parameter */
  /* dMmta [Msun Myr^(-1)] has negative sign! dt -- timestep [Myr], Ma [Msun] BH mass */
 double Rbh,J,Aspin,faspin;
 double E,zz1,zz2,Rms,Rmb,Rlso;
 double lacc,delJ,J_new;
 double Eacc,delE,delM;
 double delMrest;
 int i,n;

 if(Ka!=14)
   fprintf(fp0,"error: wrong accretor type in spin_evol()!\n");

 if(Mrest<=0.1) n=1;                  /* accrete onto BH no more than 0.1 Msun at one step */
 else if(Mrest>0.1 && Mrest<=1.0) n=10;
 else if(Mrest>1.0 && Mrest<=10.0) n=100;
 else if(Mrest>10.0 && Mrest<=100.0) n=1000;
 else if(Mrest>100.0) {
  n=10000;
  fprintf(fp0,"warning: in spin_evol() too much accretion, Mrest: %f!, idum: %d, iidd_old: %d\n",Mrest,idum_run,iidd_old);
 }

 (*Mdisa)=0.0;
 delMrest=Mrest/(double)n;
 for(i=0;i<n;i++) {

   Rbh=Rbhf(*Ma);                     /* Schwarzschild radius of BH [Rsun] */
   J=(*aspin)*(*Ma)*(*Ma)*GGG/CCC;    /* BH spin ang. momentum [Msun Rsun^2 Myr^(-1)] */
   Aspin=J/((*Ma)*CCC);               /* other spin parameter: Aspin=J/Mc=aspin(GM/c^2) [Rsun] */
   faspin=1.0-sqrt(0.5*(1.0+sqrt(1.0-(*aspin)*(*aspin))));   /* BH energy spin function, eq.12 [unitless] */
   E=faspin*(*Ma)*CCC*CCC;            /* BH spin energy [Msun Rsun^2 Myr^(-2)] = [erg] */

                                      /* zz1, zz2 -- unitless parameters */
   zz1=1.0+pow(1.0-(4.0*Aspin*Aspin/(Rbh*Rbh)),0.333333)*(pow(1.0+(2.0*Aspin/Rbh),0.333333)+pow(1.0-(2.0*Aspin/Rbh),0.333333));
   zz2=sqrt(12.0*Aspin*Aspin/(Rbh*Rbh)+zz1*zz1);
   Rms=0.5*Rbh*(3.0+zz2-sqrt((3.0-zz1)*(3.0+zz1+2.0*zz2)));  /* radius of marginaly stable orbit [Rsun]: LARGER */
   Rmb=Rbh-Aspin+sqrt(Rbh*(Rbh-2.0*Aspin));                  /* radius of marginaly bound orbit [Rsun]: SMALLER */
   Rlso=Rms;
                                      /* specific ang. momentum of accreted test particle [Rsun^2 Myr^(-1)] */
   lacc=CCC*sqrt(0.5*Rbh*Rlso)*((Rlso*Rlso-Aspin*sqrt(2.0*Rbh*Rlso)+Aspin*Aspin)/(Rlso*sqrt(Rlso*Rlso-1.5*Rbh*Rlso+Aspin*sqrt(2.0*Rbh*Rlso))));
   delJ=lacc*delMrest;                /* change of BH spin angular momentum [Msun Rsun^2 Myr^(-1)] */
   J_new=J+delJ;
                                      /* specific energy  of accreted test particle [Rsun^2 Myr^(-2)] */
   Eacc=CCC*CCC*((Rlso*Rlso-Rbh*Rlso+Aspin*sqrt(0.5*Rbh*Rlso))/(Rlso*sqrt(Rlso*Rlso-1.5*Rbh*Rlso+Aspin*sqrt(2.0*Rbh*Rlso))));
   delE=Eacc*delMrest;                /* change of BH spin energy [Msun Rsun^2 Myr^(-2)] = [erg] */
   delM=delE/(CCC*CCC);               /* change of BH mass (grav mass difference) [Msun] */

   (*Ma)=(*Ma)+delM;
   (*aspin)=J_new*CCC/((*Ma)*(*Ma)*GGG);
   (*Mdisa)+=(delMrest-delM);         /* gravitational mass (positive) that disappeared from the system [Msun] */

   if(delM<-0.0000001 || delM>delMrest || (*aspin)>1.000001 || (*aspin)<-0.0000001 || (*Mdisa)<0.0)
     fprintf(fp0,"error: wrong calculation in spin_evol(). idum: %d, iidd_old: %d\n",idum_run,iidd_old);
 }

}


double bhspininit1(double Mco)
{ /* returns initial a_spin for BH formed out of a star with CO core mass of Mco */
  /* June 3, 2017: Georges Meynet and Chris Fryer simulations */
 double aspin,f1,f2;

 if(Mco<1.0 || Mco>500.0) {
   fprintf(fp0,"error: in bhspininit1(); Mco mass wrong: %f\n",Mco);
   fflush(fp0);
 }

 if(Mco<=13.0)
   aspin=0.9;
 else if(Mco>=27.0)
   aspin=0.0;
 else
   aspin=-0.064*Mco+1.736;

 return aspin;
}


double bhspininit2(double Mco)
{ /* returns initial a_spin for BH formed out of a star with CO core mass of Mco */
  /* Oct 8, 2017: Georges Meynet and Chris Fryer simulations, and dependence on Z included */
  /* Genva modles with meridional currents and with initial rotation=40% breakup velocity */
 double aspin,m11,m22,a11,b11,alow;

 if(Mco<1.0 || Mco>500.0) {
   fprintf(fp0,"error: in bhspininit2(); Mco mass wrong: %f\n",Mco);
   fflush(fp0);
 }

 a11=-0.088;

 if(ZZ>=0.00916) {                     /* Z=0.014 */
   b11=2.258;
   m11=16.0;
   m22=24.205;
   alow=0.128;
 }
 else if(ZZ<0.00916 && ZZ>=0.00346) {  /* Z=0.006 */
   b11=3.578;
   m11=31.0;
   m22=37.818;
   alow=0.25;
 }
 else if(ZZ<0.00346 && ZZ>=0.00089) {  /* Z=0.002 */
   b11=2.434;
   m11=18.0;
   m22=27.659;
   alow=0.0;
 }
 else if(ZZ<0.00089) {                 /* Z=0.0004 */
   b11=3.666;
   m11=32.0;
   m22=38.818;
   alow=0.25;
 }
 else
   fprintf(fp0,"error: in bhspininit2 -- unrecognized metallicity: Jakub challenge\n");

 if(Mco<=m11)
   aspin=0.85;
 else if(Mco>=m22)
   aspin=alow;
 else
   aspin=a11*Mco+b11;

 return aspin;
}


double bhspininit3(double Mco)
{ /* returns initial a_spin for BH formed out of a star with CO core mass of Mco */
  /* Feb 26, 2019: Carl Fields and Chris Fryer simulations, and dependence on Z included */
  /* MESA modles with Magnetic dynamo and with initial rotation=40% breakup velocity */
 double aspin,m11,a11,b11,a22,b22;

 if(Mco<1.0 || Mco>500.0) {
   fprintf(fp0,"error: in bhspininit2(); Mco mass wrong: %f\n",Mco);
   fflush(fp0);
 }


 if(ZZ>=0.00916) {                     /* Z=0.014 */
   a11=-0.0016;
   b11=0.115;
   aspin=a11*Mco+b11;
 }
 else if(ZZ<0.00916 && ZZ>=0.00346) {  /* Z=0.006 */
   a11=-0.0006;
   b11=0.105;
   aspin=a11*Mco+b11;
 }
 else if(ZZ<0.00346 && ZZ>=0.00089) {  /* Z=0.002 */
   m11=12.09;
   a11=0.0076;
   b11=0.050;
   a22=-0.0019;
   b22=0.165;
   if(Mco<=m11) aspin=a11*Mco+b11;
   else aspin=a22*Mco+b22;
 }
 else if(ZZ<0.00089) {                 /* Z=0.0004 */
   a11=-0.0010;
   b11=0.125;
   aspin=a11*Mco+b11;
 }
 else
   fprintf(fp0,"error: in bhspininit3 -- unrecognized metallicity: Jakub challenge\n");

 return aspin;
}


double bhspininit10(void)
{ /* returns initial a_spin for BH formed out of a star */
  /* Jim Fuller model that does not depend on star properties */
 double aspin;

 aspin=0.01;

 return aspin;
}


/*--------------------------------- detailed MT computation ----------------------------------------*/

double dMmtf(double Ma, double Mb, double Ra, double Rb, double La, double Lb, double wa, double Ia, double KTa,
             double wb, double Ib, double KTb, int Ka, int Kb, double a, double e, double t, double dt,
             double *input1, int *input2, double dMmta_old, double dMmtb_old, int *stop2, int *mttype, int mttypelast,
             int *doce, int *merger, double *dMth, double *Tth, int *dec, double aspinb, double Mzamsa, double Mzamsb)
{ /* calculates mass loss rate from donor: dMmt [Msun/Myr] POSITIVE NUMBER */
  /* star A donor, star B accretor, needs previous step accretion rate (or if none then give 0.0)  */
  /* here it is set to BBeta=specific ang.momentum of compact accretor */
  /* FFa -- fraction of mass lost by donor attached to companion */
  /* here BBeta and FFa are not free parameters!!! */
 double dMmt,dMmtt,dMedd,Rla,Rlb,DDIcrit;
 double Zstar,Zlobe,dlnR,tmb_a,tmb_b,tth,tgr,tmt;
 double q,per,BBeta,FFa,Jgr,D;
 double c=2.99792458e+10;                       /* speed of light [cm/s], tr4 below change it to [Rsun/day] */

 dMedd=dMeddf(Kb,Rb,Ka);
 Rla=roche(Mb/Ma,a);
 Rlb=roche(Ma/Mb,a);

 per=2.0*Pi*a*sqrt(a/(G*(Ma+Mb)));                                                                   /* [days] */
 Jgr=(-32.0*pow(2.0*Pi,8.0/3.0)*pow(G,5.0/3.0)*Ma*Mb*pow(Ma+Mb,-1.0/3.0)*pow(per,-8.0/3.0))/(5.0*pow(tr4*c,5.0));
 Jgr*=(365.25*1.0e+06);                                                           /* Jgr -- dJ_gr/J in [1/Myr] */

 if(CE==1 || CE==2 || CE==3) {
   if((Ka>=10 && Ka<=12) || Ka==16 || Ka==17) {                                /* calculates MT loss from WD donor */
     q=Mb/Ma;
     BBeta=Ma*Ma/pow(Ma+Mb,2.0);                  /* specific ang.momentum of accretor: WD/NS/BH and lost material */
                                                  /* BBeta above assumed for remnant accretor, so check here */

     if(dMmta_old<acc)                  /* FFa=1 cons MT, FFa=0 completely non-cons MT */
       FFa=1.0;                         /* assumes conservative evolution for the first time step */
     else
       FFa=min(1.0,fabs(dMmtb_old)/dMmta_old);    /* min(1.0,...) just assures that I am not slightly above 1.0 */
                                                  /* due to numerical inaccuracies, because dMmtb_old<=dMmta_old! */
                                                  /* Eddington limit was imposed already in dMgainf()! */
     Zstar=dlnRdlnM(Ma,input1,input2,stop2);
     D=(5.0/6.0+Zstar/2.0)-((1.0-FFa)/(3.0*(1.0+q)))-(((1.0-FFa)*BBeta*(1.0+q)+FFa)/q);
     dMmt=-Ma*Jgr/D;                                                                                 /* [Msun/Myr] */
     (*mttype)=5;

     if(D<0.0 || (dMmt>=2.0*dMedd*Rlb/Rb && TRwd==1)) {        /* dynamicaly unstable MT: merger */
       dMmt=0.0;                                               /* D<0.0 -- regular dynamical instability */
       (*merger)=1;                               /* dMmt>... -- trapping radius exceeds the roche lobe radius (CE) */
     }
   }

   else {                                         /* non-degenerate donor */
     dlnR=dlnRdt(t,dt,input1,input2,stop2);
     Zstar=dlnRdlnM(Ma,input1,input2,stop2);
     Zlobe=dlnRl(Ka,Kb,Ma,Mb,a,dMmta_old,dMmtb_old);

     if(fabs(Zstar-Zlobe)<acc) {                       /* stop MT functions */
       dMmt=0.0;                                       /* 2) if Zstar-Zlobe=0, do not calculate dMmt because it is a denominator */
       return dMmt;
     }

     tgr=-(1.0/Jgr);
     if((*Tth)<acc) {                                    /* first step of thermal MT */
       (*Tth)=tth=tthf(Ma,Ra,La,Ka);
       (*dMth)=Ma/tth;
       tmb_a=tmbf(Ma,Ra,wa,Mb,a,La,Ka);
       tmb_b=tmbf(Mb,Rb,wb,Ma,a,Lb,Kb);
     }
     else {                                              /* one of subsequent steps of thermal MT, dMth is set */
       tth=(*Tth);
       tmb_a=tmbf(Ma,1.01*Rla,wa,Mb,a,La,Ka);                 /* during thermal MT donor has radius equal to roche lobe radius */
       tmb_b=tmbf(Mb,1.01*Rlb,wb,Ma,a,Lb,Kb);
     }
     dMmtt=(Ma*(dlnR+2.0/(tmb_a+tmb_b)+2.0/tgr))/(Zstar-Zlobe);   /* nuc/MB/GR/Tide MT rate [Msun/Myr] */
     tmt=Ma/dMmtt;                                                /* and timescale [Myr] */

     if(Ra<=1.1*Rla)                                     /* Rdonor almost = Rlobe (within 50%), Thermal MT is allowed */
       (*dec)=0;                                         /* to switch to diffrent MT, donor ragained its equilibrum */

     if(Ka>=0 && Ka<=6) DDIcrit=DDIcrit1;
     else if(Ka==7) DDIcrit=DDIcrit2;
     else if(Ka==8 || Ka==9) DDIcrit=DDIcrit3;
     else fprintf(fp0,"unknown type Ka: %d for DDI, %d\n",Ka,iidd_old);

                                                         /* for all non-degen. donors, Kdonor=0,1,2,3,4,5,6,7,8,9 */
     if(Ma>=DDIcrit*Mb) {                                /* deleyed dynamical or dynamical instability */
       dMmt=0.0;
       if(Ka>=2 && Ka<=9 && Ka!=7)                       /* donor Ka=2,3,4,5,6,8,9 giant -- CE evolution */
         (*doce)=1;
       else
         (*merger)=1;                                    /* donor K=0,1,7: merger assumed (no defined core) */
     }
     else if(dMmtt>0.0 && tmt>tth && (*dec)==0) {         /* nuclear/MB/GR/Tide timescale MT */
       dMmt=dMmtt;
       (*mttype)=1;
     }
     else if(dMmtt>0.0 && tmt<=tth && (*dec)==0) {	       /* thermal timescale MT */
       dMmt=(*dMth);                                     /* used also if the GR or MB or Tide dominates over Thermal */
       (*dec)++;
       (*mttype)=4;
     }
     else {                                               /* dMmt has wrong sign, checking what to do: */
       (*doce)=decide(Ma,Mb,a,e,Ka,Kb,Ra,Rb,Zstar,*dMth,aspinb);
       (*dec)++;
       if((*doce)==1 && (Ka>=2 && Ka<=9 && Ka!=7))        /* dynamically unstable MT  -> CE (giant-like donor) */
         dMmt=0.0;
       else if((*doce)==1) {                              /* merger (nongiant donor) */
         dMmt=0.0;
         (*doce)=0;
         (*merger)=1;
       }
       else {                                             /* thermal timescale MT */
         dMmt=(*dMth);
         (*mttype)=4;
       }
     }

     if(dMmt>=2.0*dMedd*Rlb/Rb) {            /* trapping radius exceeded roche lobe radius of accretor */
       dMmt=0.0;
       (*mttype)=0;
       if(Ka==2 || Ka==3 || Ka==4 || Ka==5 || Ka==6 || Ka==8 || Ka==9)        /* donor giant: do CE */
         (*doce)=1;
       else                                                                   /* donor non-giant: merger */
         (*merger)=1;
     }
   }

   if(dMmt<-acc) {                           /* stop evolution of system: wrong sign of dMmt!!! */
     (*stop2)=1;
     fprintf(fp0,"error: in dMmtf() dMmt has wrong sign!!!, Ka: %d, Kb: %d, %d\n",Ka,Kb,iidd_old);
     fflush(fp0);
   }
                                               /* to protect turning on thermal MT, right at the end */
				               /* of nuclear MT when dlnR/dt becomes negative: star shrinks */
				               /* then dMmtt is negative and instead of stopping MT, we used */
				               /* to get one step of thermal MT */
				               /* only sick case: if mtype=1 changes to 4 naturally, then we */
				               /* treat the first thermal time step as mtype=4, but with dMmt=nuclear */
   if(mttypelast==1 && (*mttype)==4) dMmt=fabs(dMmtt);
 }
 else if(CE==4) {
   if((Ka>=10 && Ka<=12) || Ka==16 || Ka==17) {                                /* calculates MT loss from WD donor */
     q=Mb/Ma;
     BBeta=Ma*Ma/pow(Ma+Mb,2.0);                  /* specific ang.momentum of accretor: WD/NS/BH and lost material */
							/* BBeta above assumed for remnant accretor, so check here */

     if(dMmta_old<acc)                                              /* FFa=1 cons MT, FFa=0 completely non-cons MT */
       FFa=1.0;                                          /* assumes conservative evolution for the first time step */
     else
       FFa=min(1.0,fabs(dMmtb_old)/dMmta_old);    /* min(1.0,...) just assures that I am not slightly above 1.0 */
						  /* due to numerical inaccuracies, because dMmtb_old<=dMmta_old! */
						  /* Eddington limit was imposed already in dMgainf()! */
     Zstar=dlnRdlnM(Ma,input1,input2,stop2);
     D=(5.0/6.0+Zstar/2.0)-((1.0-FFa)/(3.0*(1.0+q)))-(((1.0-FFa)*BBeta*(1.0+q)+FFa)/q);
     dMmt=-Ma*Jgr/D;                                                                                 /* [Msun/Myr] */
     (*mttype)=5;

     if(D<0.0 || (dMmt>=2.0*dMedd*Rlb/Rb && TRwd==1)) {        /* dynamicaly unstable MT: merger */
       dMmt=0.0;                                               /* D<0.0 -- regular dynamical instability */
       (*merger)=1;                               /* dMmt>... -- trapping radius exceeds the roche lobe radius (CE) */
     }
   }

   else {                                           /* non-degenerate donor */
     dlnR=dlnRdt(t,dt,input1,input2,stop2);
     Zstar=dlnRdlnM(Ma,input1,input2,stop2);
     Zlobe=dlnRl(Ka,Kb,Ma,Mb,a,dMmta_old,dMmtb_old);

     if(fabs(Zstar-Zlobe)<acc) {                    /* stop MT functions */
       dMmt=0.0;                                    /* 2) if Zstar-Zlobe=0, do not calculate dMmt because it is a denominator */
       return dMmt;
     }

     tgr=-(1.0/Jgr);
     if((*Tth)<acc) {                               /* first step of thermal MT */
       (*Tth)=tth=tthf(Ma,Ra,La,Ka);
       (*dMth)=Ma/tth;
       tmb_a=tmbf(Ma,Ra,wa,Mb,a,La,Ka);
       tmb_b=tmbf(Mb,Rb,wb,Ma,a,Lb,Kb);
     }
     else {                                                   /* one of subsequent steps of thermal MT, dMth is set */
       tth=(*Tth);
       tmb_a=tmbf(Ma,1.01*Rla,wa,Mb,a,La,Ka);                 /* during thermal MT donor has radius equal to roche lobe radius */
       tmb_b=tmbf(Mb,1.01*Rlb,wb,Ma,a,Lb,Kb);
     }
     dMmtt=(Ma*(dlnR+2.0/(tmb_a+tmb_b)+2.0/tgr))/(Zstar-Zlobe);   /* nuc/MB/GR/Tide MT rate [Msun/Myr] */
     tmt=Ma/dMmtt;                                                /* and timescale [Myr] */

     if(Ra<=1.1*Rla)                                        /* Rdonor almost = Rlobe (within 50%), Thermal MT is allowed */
       (*dec)=0;                                            /* to switch to diffrent MT, donor ragained its equilibrum */

     if(Ka>=0 && Ka<=6) DDIcrit=DDIcrit1;
     else if(Ka==7) DDIcrit=DDIcrit2;
     else if(Ka==8 || Ka==9) DDIcrit=DDIcrit3;
     else fprintf(fp0,"unknown type Ka: %d for DDI, %d\n",Ka,iidd_old);

     if(Ka>1 && Ka<7 && Mzamsa>18.0) {                      /* stable/unstable MT, stability diagram for giant type stars */
       if(Mzamsa>18.0 && Mzamsa<60.0) {
         if(ZZ>=Zrlof)
	   DDIcrit=19.6/7.0;
         if(ZZ<Zrlof)
           DDIcrit=19.6/7.0;
       }
       if(Mzamsa>=60.0 && Mzamsa<80.0) {
         if(ZZ>=Zrlof)
	   DDIcrit=41.0/12.0;
         if(ZZ<Zrlof)
	   DDIcrit=56.8/12.0;
       }
       if(Mzamsa>=80.0)
         DDIcrit=74.4/14.0;

       if(((Ma>=DDIcrit*Mb) && ((ZZ>=Zrlof && Ra>(62.32*Ma-515.25)) || (ZZ<Zrlof && Ra<((-0.293*Ma*Ma)+30.289*Ma-498.218)) || (ZZ<Zrlof && Ra>(26.28*Ma+262.27)))) || (Ma/Mb>8.0))
         (*doce)=1;
       else if(dMmt>=2.0*dMedd*Rlb/Rb) {	            /* trapping radius exceeded roche lobe radius of accretor */
         (*doce)=1;
         evroute_add("TRAP ");
       }
       else if(dMmtt>0.0 && tmt>tth && (*dec)==0) {         /* nuclear/MB/GR/Tide timescale MT */
         dMmt=dMmtt;
         (*mttype)=1;
       }
       else if(dMmtt>0.0 && tmt<=tth && (*dec)==0) {	    /* thermal timescale MT */
         dMmt=(*dMth);                                      /* used also if the GR or MB or Tide dominates over Thermal */
         (*dec)++;
         (*mttype)=4;
       }
       else {
         (*dec)++;
         dMmt=(*dMth);                                      /* used also if the GR or MB or Tide dominates over Thermal */
         (*mttype)=4;
       }
     }
     else if((Ka>1 && Ka<7 && Mzamsa<=18.0) || (Ka>7 && Ka<10)){
						            /* for all non-degen. donors, Kdonor=0,1,2,3,4,5,6,7,8,9 */
       if(Ma>=DDIcrit*Mb) {                                 /* deleyed dynamical or dynamical instability */
         dMmt=0.0;
         (*doce)=1;
       }
       else if(dMmtt>0.0 && tmt>tth && (*dec)==0) {         /* nuclear/MB/GR/Tide timescale MT */
         dMmt=dMmtt;
         (*mttype)=1;
       }
       else if(dMmtt>0.0 && tmt<=tth && (*dec)==0) {	    /* thermal timescale MT */
         dMmt=(*dMth);                                      /* used also if the GR or MB or Tide dominates over Thermal */
         (*dec)++;
         (*mttype)=4;
       }
       else {                                               /* dMmt has wrong sign, checking what to do: */
         (*doce)=decide(Ma,Mb,a,e,Ka,Kb,Ra,Rb,Zstar,*dMth,aspinb);
         (*dec)++;
         if((*doce)==1)                                     /* dynamically unstable MT  -> CE (giant-like donor) */
           dMmt=0.0;
         else {                                             /* thermal timescale MT */
           dMmt=(*dMth);
           (*mttype)=4;
         }
       }

       if(dMmt>=2.0*dMedd*Rlb/Rb) {                      /* trapping radius exceeded roche lobe radius of accretor */
         dMmt=0.0;
         (*mttype)=0;
         (*doce)=1;
       }
     }
     else if(Ka<2 || Ka==7) {
                                                         /* for all non-degen. donors, Kdonor=0,1,7 */
       if(Ma>=DDIcrit*Mb) {                              /* deleyed dynamical or dynamical instability */
         dMmt=0.0;
         (*merger)=1;                                    /* donor K=0,1,7: merger assumed (no defined core) */
       }
       else if (dMmt>=2.0*dMedd*Rlb/Rb) {                /* trapping radius exceeded roche lobe radius of accretor */
         dMmt=0.0;
         (*merger)=1;
       }
       else if(dMmtt>0.0 && tmt>tth && (*dec)==0) {      /* nuclear/MB/GR/Tide timescale MT */
         dMmt=dMmtt;
         (*mttype)=1;
       }
       else if(dMmtt>0.0 && tmt<=tth && (*dec)==0) {	 /* thermal timescale MT */
         dMmt=(*dMth);                                   /* used also if the GR or MB or Tide dominates over Thermal */
         (*dec)++;
         (*mttype)=4;
       }
       else {                                       /* dMmt has wrong sign, checking what to do (stable MT or merger): */
         (*doce)=decide(Ma,Mb,a,e,Ka,Kb,Ra,Rb,Zstar,*dMth,aspinb);
         (*dec)++;
         if((*doce)==1) {                           /* merger  */
           dMmt=0.0;
           (*doce)=0;
           (*merger)=1;
         }
         else {                                     /* thermal timescale MT */
           dMmt=(*dMth);
           (*mttype)=4;
         }
       }
     }
   }
   if(dMmt<-acc) {                                  /* stop evolution of system: wrong sign of dMmt!!! */
     (*stop2)=1;
     fprintf(fp0,"error: in dMmtf() dMmt has wrong sign!!!, Ka: %d, Kb: %d, %d\n",Ka,Kb,iidd_old);
     fflush(fp0);
   }
                                                    /* to protect turning on thermal MT, right at the end */
					            /* of nuclear MT when dlnR/dt becomes negative: star shrinks */
						    /* then dMmtt is negative and instead of stopping MT, we used */
						    /* to get one step of thermal MT */
						    /* only sick case: if mtype=1 changes to 4 naturally, then we */
						    /* treat the first thermal time step as mtype=4, but with dMmt=nuclear */
   if(mttypelast==1 && (*mttype)==4) dMmt=fabs(dMmtt);
 }


 return dMmt;    /* [Msun/Myr] */
}


double dMgainf(double dMmta, double dMtran, double Mb, int Ka, int Kb, double Rb, double *Mout1, double *Mout2,
               int *mark_out, int *doce, int *merger, int *mark777, double *Mb0, double aspinb)
{ /* returns mass accretion (gain) [Msun/Myr] onto accretor B from donor A loosing mass at rate dMmta */
  /* accretion rate has negative sign!!! */
  /* version used after Mar 02, 2004 -- correction for H accretion onto Kb=11,12 */
 double dMgain,dMedd,dMcri,tmp1,eta,Q,dMcr,lam,l0,Ys,X,a,b;
 double dMcrit4,dMcrit5,dMcrit4a,dMcrit5a,dMcrit4b,dMcrit5b,aa1,aa2,aa3,bb1,bb2;
 double ldMmta,Mshell;
 double dMcrit6a,dMcrit6b,dMcrit6c,dMcrit6d;
 int mark1;
 double dMedd1,dM0rlof,Racc,f1,f2,Rs,Rsph;

 eta=0.0;
 if(dMmta>acc && Kb==13 && dMmta<dMtran && TRA==0)   /* transient source, NS accretor only (BH accrets all material) */
   dMgain=0.0;                                       /* NS accretor not allowed to accumulate any material */
 else if(dMmta>acc) {                          /* Ma donor, both Ma,Mb will be changed due to MT in coming singl() */
   if(Kb==13 || Kb==14)                        /* compact object accretor Kb=13,14 */
     eta=1.0;
   else if(Kb==10 || Kb==17) {                       /* He WD/Hybrid WD accretor Kb=10 */
     (*Mout1)=(*Mout2)=100.0;
     if((Ka>=0 && Ka<=6) || Ka==16) {

       Ys=31557600.0;                                /* year in seconds */
       Q=6.4e+18;                                    /* energy yield of H burning [erg/g]: Nomoto etal.2007:ApJ663,1269 */
       Q=(Q*Msun)/(Ys*Lsun);                         /* [Lsun yr Msun^(-1)] */
       if(ZZ>0.01) {                                 /* Pop I stars */
         lam=8.0;
         l0=1995262.3;                               /* [Lsun] */
         X=0.7;
       }
       else {                                        /* Pop II stars */
         lam=5.0;
         l0=31622.8;
         X=0.8;
       }
       dMcr=l0*pow(Mb,lam)/(X*Q);                     /* [Msun/yr] based on Ritter 99, MN 309, 360: eq.10,12 + table 2 */
       dMcr*=1.0e+06;                                 /* [Msun/Myr] */
       if(dMmta<dMcr)                                 /* nova explosions, material lost with ang. mom. specific to accretor */
         eta=0.0;
       else {                                         /* CE -- accumulation of material on the WD */
         if((Ka>=0 && Ka<=2) || Ka==16)               /* leads to merger (R.Taam) */
           (*merger)=1;
         else                                         /* do CE with entire donor envelope */
           (*doce)=1;
       }
     }
     else if((Ka>=7 && Ka<=10) || Ka==17) {           /* He accretion: Hashimoto et al. 1986, ApJ 307, 687 */
       if(dMmta<=0.02) {                              /* if dMmta<=0.02 SN Ia for sub-Chandrasekhar mass: Mout1 */
         a=-400.0;
         b=1.34;
         (*Mout1)=a*dMmta+b;                          /* if extrapolated close to dM<=0.02 results in Mout1<0.0 then: */
         (*Mout1)=max(Mb,*Mout1);                     /* for high dM, but smaller than 0.02, blow up right away */
         eta=1.0;                                     /* Mass at outburst can't be smaller than current accretor mass */
       }
       else {                                         /* if dMmta>0.02 [Msun/Myr] then acculmulation till He flash */
         a=-0.0778;                                   /* which breaks degeneracy and changes accretor WD -> He MS star */
         b=0.128;
         (*Mout2)=Mb+(a*dMmta+b);
         (*Mout2)=max(Mb,*Mout2);                     /* if Mout2 becomes negative (for dMmta>=1.64 Msun/Myr) then Mout=Mb */
         (*Mout2)=max(0.35,*Mout2);                   /* He MS star is not supposed to be burning below 0.35 Msun, so even */
         eta=1;                                       /* if He flash breaks degeneracy below that mass it will cool off to */
       }                                              /* become WD again: i am skipping the part of He MS below 0.35 Msun */
     }
     else if(Ka==11) {                                /* rare cases: q reversal, formation of low mass (0.2Msun) C0 WD first */
       eta=1.0;                                       /* then companion becomes He MS star, goes <0.35 (in MT) becomes Hyb WD */
     }                                                /* then GR orbit decay: CO WD transfers to Hyb WD - assumed full accum. */
     else {
       fprintf(fp0,"warning: dMgainf() unexpected type of donor1, Ka: %d!, %d\n",Ka,iidd_old);
       fflush(fp0);
     }
   }
   else if(Kb==11 || Kb==12) {                                 /* CO or ONe WD accretor Kb=11,12 */

     if(newAcc1==0) {

       if((Ka>=0 && Ka<=6) || Ka==16)                            /* H accretion */
         eta=interWD(Mb,dMmta);

       else if((Ka>=7 && Ka<=10) || Ka==17) {                    /* He accretion */
         mark1=0;
         if(Mb<=0.7) {
           dMcrit4=dMcrit5=1.0e+06*pow(10.0,-7.6);               /* [Msun/Myr] */
           aa1=aa2=aa3=0.0;                                      /* should not be used */
         }
         else if(Mb>0.7 && Mb<=0.85) {
           dMcrit4=1.0e+06*pow(10.0,-6.34);                      /* [Msun/Myr] */
           dMcrit5=1.0e+06*pow(10.0,-6.5);
           aa1=-0.35; aa2=6.1; aa3=1.02;
         }
         else if(Mb>0.85 && Mb<=0.95) {
           dMcrit4=1.0e+06*pow(10.0,-6.05);                      /* [Msun/Myr] */
           dMcrit5=1.0e+06*pow(10.0,-6.88);
           aa1=-0.35; aa2=5.6; aa3=1.07;
         }
         else if(Mb>0.95 && Mb<=1.05) {
           dMcrit4=1.0e+06*pow(10.0,-5.93);                      /* [Msun/Myr] */
           dMcrit5=1.0e+06*pow(10.0,-6.92);
           aa1=-0.35; aa2=5.6; aa3=1.01;
         }
         else if(Mb>1.05 && Mb<=1.25) {
           mark1=1;
           dMcrit4=1.0e+06*pow(10.0,-5.76);                      /* [Msun/Myr] */
           dMcrit5=1.0e+06*pow(10.0,-7.06);
           aa1=-0.54; aa2=5.6; aa3=1.01;
           dMcrit4a=1.0e+06*pow(10.0,-5.76);                      /* [Msun/Myr] */
           dMcrit5a=1.0e+06*pow(10.0,-5.95);
           bb1=0.54; bb2=4.16;
           dMcrit4b=1.0e+06*pow(10.0,-5.95);
           dMcrit5b=1.0e+06*pow(10.0,-7.06);
         }
         else if(Mb>1.25 && Mb<=1.325) {
           dMcrit4=1.0e+06*pow(10.0,-5.83);                      /* [Msun/Myr] */
           dMcrit5=1.0e+06*pow(10.0,-7.35);
           aa1=-0.175; aa2=5.35; aa3=1.03;
         }
         else if(Mb>1.325) {
           dMcrit4=1.0e+06*pow(10.0,-6.05);                      /* [Msun/Myr] */
           dMcrit5=1.0e+06*pow(10.0,-7.4);
           aa1=-0.115; aa2=5.7; aa3=1.01;
         }
         else {
           fprintf(fp0,"something is wrong with He acc onto CO WD in dMgainf(), %d\n",iidd_old);
           fflush(fp0);
         }
         if(dMmta>=dMcrit4) {                                      /* accumulation of material until MCh */
           eta=1.0;                                              /* --> then SN Ia for Kb=11, AIC for Kb=12 */
           (*Mout1)=100.0;                                       /* to protect WD from SNIa in case dMmta (system) came from (88) */
           (*mark_out)=0;
         }
         else if(dMmta<dMcrit4 && dMmta>dMcrit5) {                 /* same as above but with smaller efficiency */
           if(mark1==0)
             eta=aa1*pow(log10(dMmta*1.0e-06)+aa2,2.0)+aa3;
           else if(dMmta<=dMcrit4a && dMmta>dMcrit5a && mark1==1)
             eta=aa1*pow(log10(dMmta*1.0e-06)+aa2,2.0)+aa3;
           else if(dMmta<=dMcrit4b && dMmta>=dMcrit5b && mark1==1)
             eta=bb1*log10(dMmta*1.0e-06)+bb2;
           else {
             fprintf(fp0,"warning: dMgainf() unexpected dMmta range, %d\n",iidd_old);
             fflush(fp0);
           }
           (*Mout1)=100.0;                                       /* to protect WD from SNIa in case dMmta (system) came from (88) */
           (*mark_out)=0;
         }
         else {                                                  /* (88) accumulation of 0.1Msun -then sub-Ch mass SN Ia */
           eta=1.0;                                              /* is MCh is reached ealier then: */
           if((*mark777)==0) {
             (*Mb0)=Mb;
             (*mark777)=1;
           }
           if((*mark_out)==0 && (newOut==0 || newOut==1 || newOut==3)) {              /* SN Ia for Kb=11 and  AIC for Kb=12 */
             Mshell=Mshellf(Mb,Kb,dMmta);
             (*Mout1)=Mb+Mshell;                                 /* only if entire explosion mass (Mout1) is accumulated in */
             (*mark_out)=1;                                      /* this dMmta regime WD goes sub_MCh SN Ia in singl.c */
           }
           else if(newOut==2) {
             Mshell=Mshellf(*Mb0,Kb,dMmta);
             (*Mout1)=(*Mb0)+Mshell;
           }
         }                                                       /* the accumulation to Mout1 is followed without interuption */
       }
       else if(Ka==11)                                  /* CO accretion: accumulation of material until MCh */
         eta=1;                                         /* --> then SN Ia for Kb=11, AIC for Kb=12 */
       else {
         fprintf(fp0,"warning: dMgainf() unexpected type of donor2, %d\n",iidd_old);
         fflush(fp0);
       }
     }
     else if(newAcc1==1) {
       ldMmta=log10(dMmta*1.0e-06);                  /* Msun/Myr -> msun/yr and log needed in the fit below */

       dMcrit6a=-6.84+1.349*Mb;                      /* log(dM/dt) above which we have Red Giant Configuration */
       dMcrit6a=pow(10.0,dMcrit6a)*1.0e+06;          /* [Msun/Myr] */
       dMcrit6b=-8.115+2.290*Mb;                     /* log(dM/dt) above which we have Steady Accretion */
       dMcrit6b=pow(10.0,dMcrit6b)*1.0e+06;          /* [Msun/Myr] */
       dMcrit6c=-8.233+2.022*Mb;                     /* log(dM/dt) above which we have Mild Flashes */
       dMcrit6c=pow(10.0,dMcrit6c)*1.0e+06;          /* [Msun/Myr] */
       dMcrit6d=-8.313+1.018*Mb;                     /* log(dM/dt) below which we have double detonations: Piersanti eta l. 2014, Tab A1 */
       dMcrit6d=pow(10.0,dMcrit6d)*1.0e+06;          /* [Msun/Myr] */

       if(dMmta>dMcrit6a) {                 /* Red Giant Configuration: WD allowed to reach MCh, no double detonations */
         eta=min(1.0,fabs(dMcrit6a)/fabs(dMmta));
         (*Mout1)=100.0;                    /* to protect WD from SNIa in case dMmta (system) came from (88) */
         (*mark_out)=0;
       }
       else if(dMmta>dMcrit6b) {            /* Steady Accretion: WD allowed to reach MCh, no double detonations */
         eta=1.0;
         (*Mout1)=100.0;                    /* to protect WD from SNIa in case dMmta (system) came from (88) */
         (*mark_out)=0;
       }
       else if(dMmta>dMcrit6c) {            /* Mild Flashes: WD allowed to reach MCh, no double detonations */
         eta=1.0;
         (*Mout1)=100.0;                    /* to protect WD from SNIa in case dMmta (system) came from (88) */
         (*mark_out)=0;
       }
       else if(dMmta>dMcrit6d) {                            /* WD allowed to reach MCh, no double detonations */
         eta=61.8651+12.5532*ldMmta-31.0676*Mb-3.15082*ldMmta*Mb+0.631685*ldMmta*ldMmta+4.12687*Mb*Mb;
         if(eta<0.0) eta=0.0;                             /* Damian interpolation used only in this regime */
         if(eta>1.0) eta=1.0;
         (*Mout1)=100.0;                       /* to protect WD from SNIa in case dMmta (system) came from (88) */
         (*mark_out)=0;
       }
       else {                                  /* double detonation regime: dMmta<=dMcrit6d */
         eta=1.0;
         if((*mark777)==0) {
           (*Mb0)=Mb;
           (*mark777)=1;
         }
         if((*mark_out)==0 && (newOut==0 || newOut==1 || newOut==3)) {              /* SN Ia for Kb=11 and  AIC for Kb=12 */
           Mshell=Mshellf(Mb,Kb,dMmta);
           (*Mout1)=Mb+Mshell;                                 /* only if entire explosion mass (Mout1) is accumulated in */
           (*mark_out)=1;                                      /* this dMmta regime WD goes sub_MCh SN Ia in singl.c */
         }
         else if(newOut==2) {
           Mshell=Mshellf(*Mb0,Kb,dMmta);
           (*Mout1)=(*Mb0)+Mshell;
         }
       }
     }
     else {
       fprintf(fp0,"warning: dMgainf() unknown mode of He accretion onto CO or ONeMg WD, %d\n",iidd_old);
       fflush(fp0);
     }
   }
   else if(Kb>=0 && Kb<=9)                           /* regular accretor Kb=0,1,2,3,4,5,6,7,8,9 */
     eta=Fa;
   else if(Kb==16)                                   /* H WD (Kb=16) added with no basis for it -- just my guess (KB) */
     eta=Fa;                                         /* I just let accrete any type of material ...... */
   else {
     fprintf(fp0,"error: dMgainf() unexpected type of accretor, Kb: %d!, %d\n",Kb,iidd_old);
     fflush(fp0);
   }
   dMgain=-eta*dMmta;                          /* accretion rate, NEGATIVE */

                                       /* BY HAND WE LIMIT ACCRETION (FROM WIND AND RLOF) ONTO NORMAL STARS TO EDDINGTON LIMIT */
   if(Kb>=0 && Kb<=9) {                          /* regular accretor Kb=0,1,2,3,4,5,6,7,8,9 */
     dMedd=-dMeddf(Kb,Rb,Ka);                    /* our definition of Eddington rate, POSITIVE */
     dMgain=max(dMgain,dMedd);
   }


   if(XEdd==1) {
     if(Kb==14) {                              /* Ohsuga 2006 critical limit imposed on BH accretor */
       dMcri=dMcrif(Kb,Mb);                        /* critical (Oshuga) acc. rate on BH [Msun Myr^(-1)], POSITIVE */
       if(fabs(dMgain)<=dMcri)
         dMgain=dMgain;
       else if(fabs(dMgain)<=10.0*dMcri) {
         tmp1=0.544*log10(fabs(dMgain)/dMcri);
         dMgain=-dMcri*pow(10.0,tmp1);
       }
       else {
         tmp1=0.934*log10(fabs(dMgain)/dMcri)-0.380;
         dMgain=-dMcri*pow(10.0,tmp1);
       }
     }
     else if((Kb>=10 && Kb<=13) || Kb==16 || Kb==17) {  /* Edd. limit imposed for WD/NS accret. */
       dMedd=dMeddf(Kb,Rb,Ka);                          /* our definition of Eddington rate, POSITIVE */
       if(fabs(dMgain)<=dMedd)
         dMgain=dMgain;
       else
         dMgain=-dMedd;                                 /* WD/NS accretes maximum at Eddington limit */
     }
   }
   else if(XEdd==2) {                                   /* new 2018 Samaresh/JP prescriptions on NS/BH accretion */
     if(Kb==13 || Kb==14) {
       dMedd1=dMeddf1(Mb,Ka);                             /* Eddington limit [Msun/Myr]; positive */
       f1=1.0;                                            /* (1-f1): wind mass loss from outer disk */
       dM0rlof=f1*dMmta;                                  /* accretion at Rsph [Msun/Myr]; positive */
       Rs=(2.0*GGG*Mb)/(CCC*CCC);                         /* Schwarzschild radius for NS/BH [Rsun] */
       Rsph=(27.0*dM0rlof*Rs)/(4.0*dMedd1);               /* spherization radius for NS/BH disk [Rsun */
       if(Kb==13)                                         /* accretion radius for NS or BH */
         Racc=Rnsf1(Mb);                                  /* NS radius: surface accretion */
       else
         Racc=Riscof(Mb,aspinb);                          /* ISCO radius: truncated disk for BH */
       if(dM0rlof<=dMedd1)                                /* (1-f2): wind mass loss from inner disk */
         f2=1.0;
       else
         f2=Racc/Rsph;
       dMgain=-f1*f2*dM0rlof;                               /* accumulation rate on NS/BH [Msun/Myr]; negative sign */
     }
     else if(Kb==10 || Kb==11 || Kb==12 || Kb==16 || Kb==17) {  /* Edd. limit imposed for WD accret. */
       dMedd=dMeddf(Kb,Rb,Ka);                                  /* Eddington limit [Msun/Myr]; positive */
       if(fabs(dMgain)<=dMedd)
         dMgain=dMgain;
       else
         dMgain=-dMedd;                                         /* WD accretes maximum at Eddington limit */
     }
   }
 }
 else
   dMgain=0.0;

 return dMgain;
}


double Mshellf(double Mb, int Kb, double dMmta)
{ /* calculates He shell mass required for double detonation SN1a on K=11/12 */
  /* used in dMgainf() -- from Ruiter et al. 2014, and from Damian fits to Shen&Bildsten 2009 */
 double Rb,Mshell;

 if(newOut==0)
   Mshell=Msh;
 else if(newOut==1 || newOut==2) {     /* system parameters taken at the onset of RLOF to determine Mshell */
   Rb=Rwdf(Mb,Kb);
   Mshell=pow(10.0,6.65)*pow(Rb,3.75)*pow(Mb,-0.3)*pow(dMmta*100.0,-0.57);      /* dMmta*100: 10^-8 Msun/yr */
 }                                      /* Ash refers to this as CONSTANT Mshell determination in Ruiter et al. 2014 */
 else if(newOut==3) {                   /* Damian fir to Shen&Bildsten 2009 shell mass model depends only on the mass of WD */
   if(Mb<=1.13)
      Mshell=0.620404*exp(-1.79348*Mb-1.51894*pow(Mb,2.0));
   else if(Mb>1.13)
      Mshell=7.17406*pow(10.0,-8.0)*exp(27.3697*Mb-14.8498*pow(Mb,2.0));
 }

 return Mshell;
}


double interWD(double Mb, double dMmta)
{ /* returns H accumulation efficiency for CO/ONeMg WDs */
  /* accreting WD of mass Mb [Msun], mass transfer rate is dMmta [Msun/Myr] */
  /* bilinear interpolation in table 1 Prialnik & Kovetz 1995, ApJ 445, 789 */
  /* updates from Nomoto et al. 2007: ApJ 663, 1269 */
 double dMcr,dMcr6,eta,x1,x2,y1,y2,y3,y4,t,u;
 double ya[4][6]={                                    /* accumulation efficiencies */
   {0.42,0.44,0.45,0.52,1.00,1.00},
   {0.40,0.43,0.46,0.50,0.62,1.00},
   {0.39,0.43,0.47,0.48,0.52,1.00},
   {0.35,0.42,0.47,0.50,0.59,1.00}
 };
 double x1a[4]={0.65,1.0,1.25,1.40};                   /* WD masses [Msun] */
 double x2a[6]={0.00001,0.0001,0.001,0.01,0.1,1.0};    /* mass transfer rates [Msun/Myr] */
 int i,ii,n,nn,j,k;

 if(Mb<0.6) dMcr=0.69e-01;                       /* Nomoto et al 2007 eq (6) shitty for low Mb masses */
 else dMcr=6.682e-01*(Mb-0.4453);                /* [Msun/Myr] */

 if(Mb<0.6) dMcr6=0.165e-01;                     /* Nomoto et al 2007 eq (5) breaks down below Mb=0.5357 Msun */
 else dMcr6=3.066e-01*(Mb-0.5357);               /* [Msun/Myr] */

 if(dMmta<=0.00001)                              /* strong nova explosions (outside of Table 1 limit) */
   eta=0.0;
 else if(dMmta>=dMcr6 || dMmta>=dMcr) {          /* superwind incorporated to cut into Table 1 */
   if(dMmta>=dMcr)
     eta=dMcr/dMmta;
   else
     eta=1.0;
 }
 else {                                          /* unstable H burning, some material ejections, Table 1 */
   i=4;
   n=6;
   if(Mb<=x1a[0]) x1=x1a[0]+acc;                 /* putting the accretor mass on the grid of Table 1 */
   else if(Mb>=x1a[3]) x1=x1a[3]-acc;
   else x1=Mb;
   x2=dMmta;
   j=-1;
   for(ii=0;ii<i-1;ii++)
     if(x1>=x1a[ii] && x1<x1a[ii+1]) {
       j=ii;
       break;
     }
   if(j==-1) fprintf(fp0,"error01: something is wrong in interWD(), %d\n",iidd_old);
   k=-1;
   for(nn=0;nn<n-1;nn++)
     if(x2>=x2a[nn] && x2<x2a[nn+1]) {
       k=nn;
       break;
     }
   if(k==-1) fprintf(fp0,"error02: something is wrong in interWD(), %d\n",iidd_old);
   y1=ya[j][k];
   y2=ya[j+1][k];
   y3=ya[j+1][k+1];
   y4=ya[j][k+1];
   t=(x1-x1a[j])/(x1a[j+1]-x1a[j]);
   u=(x2-x2a[k])/(x2a[k+1]-x2a[k]);
   eta=(1.0-t)*(1.0-u)*y1+t*(1.0-u)*y2+t*u*y3+(1.0-t)*u*y4;
 }

 return eta;
}


double dMeddf(int Ka, double Ra, int Kb) {
/* calculates Eddington limit on mass transfer rate [Msun/Myr] onto A */
/* A -- accretor, B -- donor, returns 0.0 if accretor type is unknown */
/* POSITIVE sign */
double X,Redd,dMedd,etaA;

 if((Kb>=0 && Kb<=6) || Kb==16)                                                     /* content of H in donor */
   X=0.7;                                                                                   /* H-rich donors */
 else
   X=0.0;                                                                              /* H-deficient donors */

 if(Ka<10) {Redd=Ra; etaA=1.0;}
 else if(Ka==10 || Ka==11 || Ka==12 || Ka==16 || Ka==17) {Redd=Ra; etaA=1.0;}
 else if(Ka==13) {Redd=Ra; etaA=1.0;}
 else if(Ka==14) {Redd=3.0*Ra; etaA=0.5;}
 else Redd=etaA=0.0;                                                         /* unknown accretor, return 0.0 */

 dMedd=(2.08797e-03*Redd)/(etaA*(1.0+X));        /* [Msun/yr] Edd. crit. acc. rate, eq. 13.7 Vanbeveren book */
 dMedd*=1.0e+06;                                 /* [Msun/Myr], (3.0*Rb) -- accretion to last stable orbit */

 return dMedd;
}


double dMeddf1(double Ma, int Kb) {
/* calculates Eddington limit on mass transfer rate [Msun/Myr] onto A (POSITIVE sign) */
/* A -- accretor, B -- donor */
double X,Redd,dMedd,etaA;

 if((Kb>=0 && Kb<=6) || Kb==16)                                                     /* content of H in donor */
   X=0.7;                                                                                   /* H-rich donors */
 else
   X=0.0;                                                                              /* H-deficient donors */

 dMedd=(4.44e-08*Ma)/(1.0+X);                                              /* [Msun/yr] Edd. crit. acc. rate */
 dMedd*=1.0e+06;                                                                               /* [Msun/Myr] */

 return dMedd;
}


double Leddf1(double Ma, int Kb) {
/* calculates Eddington limit on luminosity [erg/sec] onto A */
/* A -- accretor, B -- donor */
double Ledd,X;

 if((Kb>=0 && Kb<=6) || Kb==16)                                                     /* content of H in donor */
   X=0.7;                                                                                   /* H-rich donors */
 else if(Kb!=13 && Kb!=14)                                                             /* H-deficient donors */
   X=0.0;
 else
   fprintf(fp0,"error: wrong donor type in Leddf1(): %d\n",Kb);

 Ledd=(2.5141e+38*Ma)/(1.0+X);         /* [erg/sec] Edd. crit. luminosity: formula after eq.7 in Samaresh paper */

 return Ledd;
}


double dMcrif(int Ka, double Ma)
{ /* calculates critical (Eddington) accretion (Ohsuga 2006) rate onto a BH [Msun/Myr] */
  /* POSITIVE sign, A -- accreting BH */
 double dMcrit;

 if(Ka!=14)
   printf("error: not a BH in dMcritf()\n");

 dMcrit=2.6e-02*(Ma/10.0);        /* [Msun/Myr] */

 return dMcrit;
}


double dMtranf(double Ma, double Mb, double a, double e, int Ka, int Kb)
{ /* A: donor, B: compact accretor, calculates dMtran, source transient/persistent */
  /* if dM (mass transfer rate) > dMtran source is persistent, otherwise transient */
  /* correction Sep 09, 2005: theoretical alpha=0.1,0.005,0.02 for He,CO,O -- see Menou Fig.2 */
  /* but we use always alpha=0.1 -- observational evidence -- see Menou, S4, first parag. */
 double Rdisk,Rlb,C,Ys,dMtran;

 if(Kb==13 || Kb==14) {

   Ys=31557600;                                                                           /* year in seconds */
   Rlb=roche(Ma/Mb,a*(1.0-e));
                                                   /* if dM>dMtran source is persistent, otherwise transient */
   Rdisk=0.666667*Rlb;                                                                             /* [Rsun] */
   Rdisk*=6.9599;                                                                              /* [10^10 cm] */

   if((Ka>=0 && Ka<=6) || Ka==16) {                                                          /* H-rich donor */
     C=1.0;                                                                    /* expressed in 5.0e-04 units */
     dMtran=1.5e+15*pow(Mb,-0.4)*pow(Rdisk,2.1)*pow(C,-0.5);      /* [g/sec], Dubus et al. 1999, MN 303, 139 */
     dMtran*=(Ys*1.0e+06/Msun);                                                                /* [Msun/Myr] */
   }
   else if((Ka>=7 && Ka<=10) || Ka==17) {                                                   /* He-rich donor */
     C=1.0;                                                         /* alpha_0.1=alpha/0.1, for He alpha=0.1 */
     dMtran=5.9e+16*pow(Mb,-0.87)*pow(Rdisk,2.62)*pow(C,0.41);   /* [g/sec], Menou et al. 2002, ApJ 564, L81 */
     dMtran*=(Ys*1.0e+06/Msun);                                                                /* [Msun/Myr] */
   }
   else if(Ka==11) {                                                                             /* CO donor */
     C=1.0;                                                         /* alpha_0.1=alpha/0.1, for CO alpha=0.1 */
     dMtran=1.2e+16*pow(Mb,-0.74)*pow(Rdisk,2.21)*pow(C,0.42);   /* [g/sec], Menou et al. 2002, ApJ 564, L81 */
     dMtran*=(Ys*1.0e+06/Msun);                                                                /* [Msun/Myr] */
   }
   else if(Ka==12) {              	                                                        /* ONe donor */
     C=1.0;                                                          /* alpha_0.1=alpha/0.1, for O alpha=0.1 */
     dMtran=5.0e+16*pow(Mb,-0.68)*pow(Rdisk,2.05)*pow(C,0.45);   /* [g/sec], Menou et al. 2002, ApJ 564, L81 */
     dMtran*=(Ys*1.0e+06/Msun);                                                                /* [Msun/Myr] */
   }
   else {
     fprintf(fp0,"error: wrong type of donor in dMtranf(), Ka: %d, %d\n",Ka,iidd_old);
     fflush(fp0);
   }
 }
 else                                                                                      /* not applicable */
   dMtran=0.0;

 return dMtran;
}


int decide(double Ma, double Mb, double a, double e, int Ka, int Kb, double Ra, double Rb, double Zstar, double dMth, double aspinb)
{ /* A -- donor, B-accretor: decides either thermal MT or CE is encountered */
  /* Zstar is dlnRa/dlnMa */
 double Ma1,Mb1,q1,a1,Rl1,Ra1;
 double Ma2,q2,Rl2,Ra2;
 double Mai,Mbi,Maf,Mbf,ai,af,q,Rl,delM;
 double dJ,Jorbi,x,cr;
 double y1,y2,y3,y4;
 int i,mark,ce;

 if(e>acc) {
   fprintf(fp0,"error: in decide() e is not zero, e: %g!!!, %d\n",e,iidd_old);
   fflush(fp0);
 }

 if(Ka>=3 && Ka<=6)
   cr=CR2;
 else
   cr=CR1;


 Ma1=Ma;                              /* current position on Ma-Ra plane */
 Mb1=Mb;
 q1=Ma1/Mb1;
 a1=a;
 Rl1=roche(Mb1/Ma1 ,a1);
 Ra1=1.0*Rl1;      /* I don't use equlibrum donor radius Ra, as it may exceed by large the real  donor radius */
                   /* real donor radius should stay within its Roche lobe */
 mark=0;
 delM=Ma/1000.0;
 Maf=Ma1;
 Mbf=Mb1;
 ai=a1;
 Rl2=Rl1;
 q2=q1;
 for(i=0;i<998;i++) {
   Mai=Maf;
   Mbi=Mbf;
   Maf-=delM;
   dJ=dJf(Mai,Mbi,ai,Maf,dMth,Rb,Ka,Kb,&x,aspinb);
   Mbf=Mbf+x*delM;
   q=Maf/Mbf;
   Jorbi=Mai*Mbi*sqrt(ai*GGG*(Mai+Mbi))/(Mai+Mbi);
   af=((Maf+Mbf)*pow(Jorbi+dJ,2.0))/(GGG*Maf*Maf*Mbf*Mbf);                  /* non-conservative MT assumption */
   Rl=roche(Mbf/Maf,af);
   if(i==0) y1=Rl;
   else if(i==1) y2=Rl;
   else if(i==950) y3=Rl;
   else if(i==951) y4=Rl;
   if(Rl<Rl2) {
     Rl2=Rl;
     q2=q;
     Ma2=Maf;
     Ra2=exp(Zstar*(log(Ma2)-log(Ma1))+log(Ra1));
   }
   ai=af;
 }

 if(!(y2<y1 && y3<y4))             /* do thermal MT: no minimum, no q2 -- orbit expands */
   ce=0;
 else if(q1>=cr*q2)                /* send system to CE */
   ce=1;
 else                              /* no obvious CE, do thermal MT */
   ce=0;

 return ce;
}


double dJf(double Mai, double Mbi, double ai, double Maf, double dMth, double Rb, int Ka, int Kb, double *x, double aspinb)
{ /* A-donor, B-accretor, calculates the amount of angular momentum loss when Mdonor: Mai->Maf */
  /* fills in x: which is how much [0-1] donor material is attached to companion */
 double Jorb,Jmt,dMgain,dMtran,FFa,worb,Rcom;
 double dumb1,dumb2,dumb7;
 int dumb3,dumb4,dumb5,dumb6;

 if((Kb==13 || Kb==14) && Ka!=13 && Ka!=14) dMtran=dMtranf(Mai,Mbi,ai,0.0,Ka,Kb);
 else dMtran=0.0;
 dumb1=100.0;
 dumb2=dumb3=dumb4=dumb6=dumb7=0;
 dMgain=dMgainf(dMth,dMtran,Mbi,Ka,Kb,Rb,&dumb1,&dumb2,&dumb3,&dumb4,&dumb5,&dumb6,&dumb7,aspinb);     /* accretion rate onto B */
 FFa=min(1.0,fabs(dMgain)/dMth);
 (*x)=FFa;

 if((Kb>=10 && Kb<=14) || Kb==16 || Kb==17) {           /* compact object accretor Kb=10,11,12,13,14,16,17 */
   worb=sqrt(GGG*(Mai+Mbi))*pow(ai,-1.5);               /* mean orbital angular velocity [Myr^-1] */
   Rcom=Mai*ai/(Mai+Mbi);                               /* orbital angular momentum loss due to MT: */
   Jmt=(1.0-FFa)*(Maf-Mai)*Rcom*Rcom*worb;              /* loss with specific J of compact accretor */
 }
 else {                                                 /* regular accretor Kb=0,1,2,3,4,5,6,7,8,9 */
   Jorb=Mai*Mbi*sqrt((GGG*ai)/(Mai+Mbi));		/* orbital angular momentum loss due to MT: */
   Jmt=Beta*((1.0-FFa)*(Maf-Mai)/(Mai+Mbi))*Jorb;	/* loss with specific ang. momentum of orbital J */
 }

 return Jmt;                      /* negative as this is ang. mom. loss, conetcted to mass loss  Maf<Mai */
}


double tthf(double Ma, double Ra, double La, int Ka)
{ /* estimates thermal timescale for star A [Myrs] */
 double tth;

 tth=30.0*Ma*Ma/(Ra*La);  /* [Myrs], see eq.(2) in ApJ 458, 301: Kalogera et al. 1996 */

 return tth;
}


double tmbf(double Ma, double Ra, double wa, double Mb, double a, double La, int Ka)
{ /* estimates magnetic breaking timescale for star A [Myrs] */
  /* if MB is not important, return large number */
  /* applied for stars with convective envelopes */
 double wcrit_a,kw1,kw2,dJmagb_a,Jbin,tmb,kw3,wsun,wx,Ta;
 int magb_a;

 Ta=get_T(La,Ra);

 if(((Ka==0 || Ka==1) && Ma>=minMB && Ma<=maxMB) || Ka==3 || ((Ka==2 || Ka==4) && Ta<5370.0) || Ka==5 || Ka==6) {
   magb_a=1;
   wcrit_a=wcritf(Ma,Ka);
 }
 else                              /* no MB */
   return 10000000000000000000000000.0;


 kw1=8.8801236797060555969e-22;    /* 2.7e+47: constant: [g cm^2 s] -> [Msun Rsun^2 Myr] (New MB) */
 kw2=5.8329147475726926586e-22;    /* 3.8e-30: constant: [s cm^-2] -> [Myr Rsun^-2] (Old MB) */

 kw3=619.2;                        /* Ivanova&Taam: this is Kj=6.e+30 [dyn cm] -> units changed to [Msun Rsun^2 Myr^-2] */
 wsun=9.46e+07;                    /* sun rotational velocity [Myr^-1] */
 wx=9.46e+08;                      /* critcal saturation velocity= 10*w_sun, w_sun=3*10^-6 sec^-1=9.46*10^7 Myr^-1 */
                                   /* I use their eq.(4) with Td/Tdsun=1.0 -- from Ron */

 if(magb_a==1 && wa<=wcrit_a && MB==1)          /* change of angular momentum due to MB [(Msun Rsun^2/Myr)/Myr] */
   dJmagb_a=-kw1*sqrt(Ra/Ma)*wa*wa*wa;
 else if(magb_a==1 && wa>wcrit_a && MB==1)
   dJmagb_a=-kw1*sqrt(Ra/Ma)*wa*wcrit_a*wcrit_a;
 else if(magb_a==1 && MB==2)
   dJmagb_a=-kw2*Ma*pow(Ra,gamMB)*wa*wa*wa;
 else if(magb_a==1 && wa<=wx && MB==3)
   dJmagb_a=-kw3*pow(Ra,4.0)*pow(wa/wsun,3.0);
 else if(magb_a==1 && wa>wx && MB==3)
   dJmagb_a=-kw3*pow(Ra,4.0)*pow(wa/wsun,1.3)*pow(wx/wsun,1.7);
 else {
   fprintf(fp0,"error: unknown case in tmbf(), iidd_old: %d\n",iidd_old);
   fflush(fp0);
 }

 Jbin=Ma*Mb*sqrt(GGG*a/(Ma+Mb));	     /* bin. ang. momentum [Msun Rsun^2/Myr], rot. J of stars neglected */
 tmb=-Jbin/dJmagb_a;                         /* [Myr] */

 return tmb;
}


double Mminconvf(void)
{ /* returns value minMB: min. mass of MS star over which conv. env. developes */
  /* MS star: below minMB fully convective, over maxMB no conv env. */

 return 0.35;
}


double Mmaxconvf(void)
{ /*returns value maxMB: max. mass of MS star below which  conv. env. developes (depends on ZZ) */

 return 0.747883+55.730610*ZZ-1532.116926*ZZ*ZZ;
}


double dlnRl(int Ka, int Kb, double Ma, double Mb, double a, double dMmta_old, double dMmtb_old)
{ /* returns partial derivative (Roche Lobe Radius) dlnRl/dlnM at mass M for constant star composition */
  /* Ma -- donor, Mb -- accretor */
  /* dMedd -eddington limit on accretion rate, positive, dMrlof -RLOF mass loss from donor, positive */
  /* at first time step, when dMrlof=0 take beta=1.0 -- conservative MT, all dM acreted onto companion */
  /* later: suuply dMrlof from previous time step -- TO BE CORRECTED!!!!: iteretion needed */
 double dlnR,F,q,beta;
 double Ma1,Mb1,a1,Jorb1,Rl1,Ma2,Mb2,a2,Jorb2,Rl2,dJ,delM;

 if(Kb>=10) {                            /* compact accretor, material lost with J of compact accretor */
   q=Ma/Mb;
   F=(0.4*q+(0.333333*pow(q,0.6666667))/(1.0+pow(q,0.333333)))/(0.6*q+pow(q,0.333333)*log(1.0+pow(q,0.333333)));

   if(Ka<10)                             /* non-compact donor */
     beta=Fa;
   else if(dMmta_old<acc)
     beta=1.0;
   else
     beta=min(1.0,fabs(dMmtb_old)/dMmta_old);

   dlnR=(2.0*q*q-(1.0-beta)*q-2.0)/(1.0+q)+(1.0+beta*q)*(0.6666667-F);
 }
 else {                                  /* regular star accretor, material lost with J of binary */

   Ma1=Ma;                               /* starting point */
   Mb1=Mb;
   a1=a;
   Jorb1=Ma1*Mb1*sqrt(a1*GGG*(Ma1+Mb1))/(Ma1+Mb1);
   Rl1=roche(Mb1/Ma1,a1);

   delM=0.01*Ma;                         /* 1% mass decrease of donor */
   Ma2=Ma-delM;
   Mb2=Mb+Fa*delM;

   dJ=Beta*((1.0-Fa)*(Ma2-Ma1)/(Ma1+Mb1))*Jorb1;
   a2=((Ma2+Mb2)*pow(Jorb1+dJ,2.0))/(GGG*Ma2*Ma2*Mb2*Mb2);               /* non-conservative MT assumption */
   Jorb2=Ma2*Mb2*sqrt(a2*GGG*(Ma2+Mb2))/(Ma2+Mb2);
   Rl2=roche(Mb2/Ma2,a2);
   dlnR=(log(Rl2)-log(Rl1))/(log(Ma2)-log(Ma1));
 }

 return dlnR;
}


double dlnRdlnM(double M, double *input1, int * input2, int *stop2)
{ /* returns partial derivative dlnR/dlnM at mass M for constant star composition */
  /* change of star mass by 1% as the first mass step: h */
  /* input1[] --array of doubles, input2[] --array of ints for single() */
 double dlnR,delM,error,h,rad,err_min,rad_min;
 int i;

 for(i=1;i<=7;i++) {
   delM=1.0/(pow(10.0,(double)i));                        /* first mass step is 0.1 M */
   h=log(1.0+delM);
   h=max(1.0e-10,h);                                    /* h must be bigger than zero */
   rad=dfridr1(funclnRM,log(M-1.2*delM*M),h,&error,input1,input2,stop2);
   if(error<err_min || i==1) {
     err_min=error;
     rad_min=rad;
   }
 }

 dlnR=rad_min;
 error=err_min;
 if(error>0.25 && input2[0]==2) {                        /* allows for a bigger error for HG stars: big radius jump */
   fprintf(fp0,"warning: derivative from dlnRdlnM() not accurate, error: %g, K: %d, %d\n\n",error,input2[0],iidd_old);
   fflush(fp0);                                                 /* error>0.1 when K=2 and is about to change to K=3 */
 }
 else if(error>0.01 && input2[0]!=2) {
   fprintf(fp0,"warning: derivative from dlnRdlnM() not accurate, error: %g, K: %d, %d\n\n",error,input2[0],iidd_old);
   fflush(fp0);
 }

 return dlnR;
}


double funclnRM(double lnM, double *input1, int *input2)
{ /* returns ln of star radius at mass M for constant star composition, star given by input1 and input2 */
  /* function takes ln(M) and not M!!! */

 double Mzams,M0,M,tbeg,tvir,tend,L,R,Mc,Mhe,Mco,dt,Mpre,tstart,frac,dMmt;
 double Mtmp,Mout1,Mout2;
 int K,flag,Kp,deriv,dummy1,dummy2;

 Mzams=input1[0];  M0=input1[1];       M=input1[2];
 K=input2[0];      tbeg=input1[3];     tvir=input1[4];
 tend=input1[5];   L=input1[6];        R=input1[7];
 Mc=input1[8];     Mhe=input1[9];      Mco=input1[10];
 flag=input2[1];   dt=input1[11];      Mpre=input1[12];
 Kp=input2[2];     tstart=input1[13];  frac=input1[14];
 Mout1=Mout2=100.0;


 Mtmp=exp(lnM);                                          /* stellar mass for derivative estimation  */
 if(Mtmp<0.0 || Mtmp>=M) {                               /* singl() cannot evolve star with negative mass */
   fprintf(fp0,"error: infunclnR() wrong Mass (lnM) supplied, %d\n",iidd_old);
   fflush(fp0);
 }

 tbeg=tend;            /* point to which star is evolved: starting (time,mass) for derivative estimation */
 tend=tbeg+1.0e-06;    /* tend-time changed only very slightly [+1yr], just to get ln(R) for given mass */

 dMmt=(M-Mtmp)/(tend-tbeg);           /* proper dM to change star mass M to required Mtmp (positive) [Msun/Myr] */
 deriv=1;                             /* tells singl() that this is only derivative estimation run */

 singl(&Mzams,&M0,&M,&K,&tbeg,&tvir,&tend,&L,&R,&Mc,&Mhe,&Mco,&flag,&dt,&Mpre,&Kp,&tstart,&frac,dMmt,deriv,Mout1,Mout2,&dummy1,&dummy2,0,0.0);

 return log(R);
}


double dlnRdt(double t, double dt, double *input1, int * input2, int *stop2)
{ /* returns partial derivative dlnR/dt at time t for constant mass of star */
  /* dt -- time step for current step in while loop */
  /* input1[] --array of doubles, input2[] --array of ints for single() */
 double dlnR,error,h,rad,err_min,rad_min;
 int i;

 for(i=1;i<=5;i++) {
   h=(1.0/(pow(10.0,(double)i)))*dt;                       /* first time step is 0.1 dt */
   h=max(1.0e-10,h);
   rad=dfridr1(funclnRt,t+2.0*h,h,&error,input1,input2,stop2);
   if(error<err_min || i==1) {
     err_min=error;
     rad_min=rad;
   }
 }

 dlnR=rad_min;
 error=err_min;
                                      /* do nothing for CHeB stars, dlnR has proper value even if error is large */
 if(error>0.1 && input2[0]==4)        /* happans for big massive stars when they almost exhausted their envelope */
   ;
 if(error>0.1 && input2[0]==2) {      /* allows for a bigger error for HG stars: rapid radiu increase */
   fprintf(fp0,"warning: derivative from dlnRdt() not accurate, error: %g, K: %d, %d\n\n",error,input2[0],iidd_old);
   fflush(fp0);
 }
 else if(error>0.01 && input2[0]!=2 && input2[0]!=4) {
   fprintf(fp0,"warning: derivative from dlnRdt() not accurate, error: %g, K: %d, %d\n\n",error,input2[0],iidd_old);
   fflush(fp0);
 }

 return dlnR;
}


double funclnRt(double t, double *input1, int *input2)
{ /* returns ln of star radius at time t for constant star mass, star given by input1 and input2 */

 double Mzams,M0,M,tbeg,tvir,tend,L,R,Mc,Mhe,Mco,dt,Mpre,tstart,frac,dMmt,Mout1,Mout2;
 int K,flag,Kp,deriv,dummy1,dummy2;

 Mzams=input1[0];  M0=input1[1];       M=input1[2];
 K=input2[0];      tbeg=input1[3];     tvir=input1[4];
 tend=input1[5];   L=input1[6];        R=input1[7];
 Mc=input1[8];     Mhe=input1[9];      Mco=input1[10];
 flag=input2[1];   dt=input1[11];      Mpre=input1[12];
 Kp=input2[2];     tstart=input1[13];  frac=input1[14];
 dMmt=0.0;
 Mout1=Mout2=100.0;

 if(t<tend) {                                            /* singl() cannot evolve star backward */
   fprintf(fp0,"error: infunclnR() wrong t supplied, %d\n",iidd_old);
   fflush(fp0);
 }

 tbeg=tend;          /* point to which star is evolved: starting time for derivative estimation */
 tend=t;               /* t- point at which dfridr1() want to evaluate log(R) to get derivative */
 deriv=1;                          /* tells singl() that this is only derivative estimation run */

 singl(&Mzams,&M0,&M,&K,&tbeg,&tvir,&tend,&L,&R,&Mc,&Mhe,&Mco,&flag,&dt,&Mpre,&Kp,&tstart,&frac,dMmt,deriv,Mout1,Mout2,&dummy1,&dummy2,0,0.0);

 return log(R);
}


double dfridr1(double (*func)(double, double *, int *), double x, double h, double *err, double *input1, int *input2, int *stop2)
{ /* returns derivative of function func in point x, h-init.step, err-error of estimate */

	int i,j;
	double errt,fac,hh,**a,ans;

	if (h == 0.0) {
	  fprintf(fp0,"h must be nonzero in dfridr1 (evol. of system halted), iidd_old: %d\n",iidd_old);
          fflush(fp0);
          (*stop2)=1;                                                                    /* stop evolution of system */
          return 0.0;
        }
	a=dmatrix(1,NTAB1,1,NTAB1);
	hh=h;
	a[1][1]=((*func)(x+hh,input1,input2)-(*func)(x-hh,input1,input2))/(2.0*hh);
	*err=BIG1;
	for (i=2;i<=NTAB1;i++) {
		hh /= CON;
		a[1][i]=((*func)(x+hh,input1,input2)-(*func)(x-hh,input1,input2))/(2.0*hh);
		fac=CON2;
		for (j=2;j<=i;j++) {
			a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
			fac=CON2*fac;
			errt=DMAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
			if (errt <= *err) {
				*err=errt;
				ans=a[j][i];
			}
		}
		if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE1*(*err)) {
			free_dmatrix(a,1,NTAB1,1,NTAB1);
			return ans;
		}
	}
	free_dmatrix(a,1,NTAB1,1,NTAB1);
	return ans;
}


/*------------------------------------ ORBIT COMPUTATION -------------------------------------------*/

void orbit1(double t0, double m1, double m2, double a, double e, double i,
            double Om, double om, double tau, double *X, double *V)
{ /* taking elements of relative orbit (a,e,i,Om,om,tau) of m2 aroun m1,
     orbit1 calculates for time t0 (time elapsed from periastron passage tau)
     position X[x1,x2,x3] and velocity V[v1,v2,v3] of m1 on its orbit. */
 double E,dE,M,P,n,b,dd;
 double l1,k1,n1,l2,k2,n2;


 P=2.0*Pi*a*sqrt(a/(G*(m1+m2)));
 n=(2.0*Pi)/P;
 M=n*(t0-tau);
 E=kepler(M,e);

 l1=cos(Om)*cos(om)-sin(Om)*sin(om)*cos(i);
 k1=sin(Om)*cos(om)+cos(Om)*sin(om)*cos(i);
 n1=sin(om)*sin(i);
 l2=-cos(Om)*sin(om)-sin(Om)*cos(om)*cos(i);
 k2=-sin(Om)*sin(om)+cos(Om)*cos(om)*cos(i);
 n2=cos(om)*sin(i);

 dd=1.0-e*e;
 if(dd<0.0 && fabs(dd)<acc)                  /* ensures that argument of sqrt() when 0.0 is positive zero */
   dd=acc;                                   /* and not some number like: -0.0000000000000000000000002344 */
 b=a*sqrt(dd);                               /* semiminor axis of orbit */
 X[0]=a*l1*cos(E)+b*l2*sin(E)-a*e*l1;
 X[1]=a*k1*cos(E)+b*k2*sin(E)-a*e*k1;
 X[2]=a*n1*cos(E)+b*n2*sin(E)-a*e*n1;

 dE=n/(1.0-e*cos(E));                    /* derivative of eccentric anomaly */
 V[0]=-a*l1*sin(E)*dE+b*l2*cos(E)*dE;
 V[1]=-a*k1*sin(E)*dE+b*k2*cos(E)*dE;
 V[2]=-a*n1*sin(E)*dE+b*n2*cos(E)*dE;
}


int orbit2(double texp, double m1, double m2, double *X, double *V,
           double *a, double *e, double *i, double *Om, double *om,
           double *tau)
{ /* counts from position (at time texp)  X[x1,x2,x3] and velocity
     V[v1,v2,v3] of m2, its relative orbit around m1, and fills 6
     orbit parameters: (a,e,i,Om,om,tau) where tau is time of
     perihelion passage. returns 0 if orbit not bound, and 1 if it
     is bound */
 double Un,Tn,En,J,J2,j1,j2,j3,alpha,mi,r,dr;
 double A,B,C,tmp1,tmp2,tmp10;
 double p,b,P,Teta,teta,E,M,n,dd,ff;
 int is;


 mi=(m1*m2)/(m1+m2);
 r=dist(0.0,0.0,0.0,X[0],X[1],X[2]);
 r=max(acc,r);
 alpha=G*m1*m2;

 A=X[1]*V[2]-X[2]*V[1];
 B=X[2]*V[0]-X[0]*V[2];
 C=X[0]*V[1]-X[1]*V[0];
 j1=mi*A;
 j2=mi*B;
 j3=mi*C;
 J2=j1*j1+j2*j2+j3*j3;
 if(J2<0 && fabs(J2)<acc)                    /* ensures that argument of sqrt() when 0.0 is positive zero */
   J2=acc;                                   /* and not some number like: -0.0000000000000000000000002344 */
 J=sqrt(J2);

 Un=(-alpha/r);
 Tn=0.5*mi*(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
 En=Tn+Un;



 if(En<0.0) {
   is=1;
   p=J2/(alpha*mi);

   dd=1.0+(2.0*En*J2)/(alpha*alpha*mi);
   if(dd<0 && fabs(dd)<acc)                    /* ensures that argument of sqrt() when 0.0 is positive zero */
     dd=acc;                                   /* and not some number like: -0.0000000000000000000000002344 */
   (*e)=sqrt(dd);                              /* small positive number so sqrt() would return a real value */
   if(fabs((*e)-1.0)<acc)
     (*e)=1.0-acc;                          /* ensures that e in never equal exactly 1.0, but a bit smaller */
   (*e)=max(acc,*e);                     /* ensures that e is never equal 0.0, to protect division by e=0.0 */
                                         /* as could have been encountered below */
   (*a)=p/(1.0-(*e)*(*e));
   b=p/sqrt(1.0-(*e)*(*e));
   P=2.0*Pi*(*a)*sqrt((*a)/(G*(m1+m2)));

   ff=A*A+B*B+C*C;
   if(ff<0 && fabs(ff)<acc)                    /* ensures that argument of sqrt() when 0.0 is positive zero */
     ff=acc;                                    /* and not some number like: -0.0000000000000000000000002344 */
   (*i)=acos(C/sqrt(ff));
   tmp10=sin(*i);                            /* ensures sin(i) is never equal 0.0,protect division by sin(i)=0.0 */
   tmp10=max(acc,tmp10);
   (*Om)=get_ang(A/(sqrt(ff)*tmp10),-B/(sqrt(ff)*tmp10));

   tmp2=(p-r)/(r*(*e));
   dr=(X[0]*V[0]+X[1]*V[1]+X[2]*V[2])/r;
   tmp1=(mi*p*dr)/(J*(*e));
   teta=get_ang(tmp1,tmp2);
   tmp2=(X[0]*cos((*Om))+X[1]*sin((*Om)))/r;
   tmp1=X[2]/(r*tmp10);
   Teta=get_ang(tmp1,tmp2);
   (*om)=Teta-teta;
   if((*om)<0.0)
     (*om)+=(2.0*Pi);

   E=2.0*atan(tan(teta/2.0)*sqrt((1.0-(*e))/(1.0+(*e))));
   if(E<0.0) E=2.0*Pi+E;
   M=E-(*e)*sin(E);
   n=(2.0*Pi)/P;
   (*tau)=(texp-M/n);
   if((*tau)<0.0) (*tau)=P+(*tau);
                   /* time of perihelion passage as reffered to given texp
                      (texp=100,peryhelion 20 dni po texp, to tau=120)
                      could be only + (positive) because calculated from point
                      of texp along the orbit (in direction of m2 moving)
                      - means that m2 will pass perihelion in tau days  */
 }
 else {
   is=0;
 }

 return is;
}


double kepler(double M, double e)
{
 double accuracy=0.000001;     /* defines accuracy with wich E is found */
 double E,E0,M0,tmp1,tmp2,delE,delM;
 int i;

 tmp1=1000000.0;
 for(i=0;i<360;i++) {
   E0=(i*Pi)/180.0;
   M0=E0-e*sin(E0);
   delM=M-M0;
   if(fabs(delM)<tmp1) {
     tmp1=delM;
     tmp2=E0;
   }
 }
 E0=E=tmp2;
 M0=E0-e*sin(E0);
 delM=M-M0;

 while(fabs(delM)>accuracy) {
   delE=delM/(1.0-e*cos(E0));
   E=E0+delE;
   M0=E-e*sin(E);
   delM=M-M0;
   E0=E;
 }

 return E;
}


double get_ang(double t1,double t2)
{ /* calculate angle: ang [rad], from the t1=sin(ang) and t2=cos(ang) */
 double ang;


 if(t1>=0.0 && t2>=0.0) {
   ang=asin(t1);
 }
 else if(t1>=0.0 && t2<=0.0) {
   ang=Pi-asin(t1);
 }
 else if(t1<=0.0 && t2<=0.0) {
   ang=Pi-asin(t1);
 }
 else if(t1<=0.0 && t2>=0.0) {
   ang=2*Pi+asin(t1);
 }
 else
   fprintf(fp0,"something is wrong in get_ang()\n");

 return ang;
}


double dist(double a1,double a2,double a3,double b1,double b2,double b3)
{
 return sqrt(pow(a1-b1,2.0)+pow(a2-b2,2.0)+pow(a3-b3,2.0));
}



/*-------------------------------------- TOOLS ---------------------------------------------*/

void copyV(double *V1, double *V2)
{ /* copy array V1[3] to array V2[3] */

 V2[0]=V1[0];
 V2[1]=V1[1];
 V2[2]=V1[2];
}


void copySN(double Maold, int Kaold, double Mbold, int Kbold, double aold, double eold, double Ma,
            double *MpgaA, int *KpgaA, double *MpgbA, int *KpgbA, double *apgA, double *epgA, double *MendaA)
{ /* copy values of first 6 parameters to 6 next (used to store values before last SN) */

 (*MpgaA)=Maold;
 (*KpgaA)=Kaold;
 (*MpgbA)=Mbold;
 (*KpgbA)=Kbold;
 (*apgA)=aold;
 (*epgA)=eold;
 (*MendaA)=Ma;
}


void copySin(double Mzamsa, double M0a, double Ma, int Ka, double tbeg, double tvira, double tenda, double La,
             double Ra, double Mca, double Mhea, double Mcoa, int flaga, double dta, double Mprea, int Kpa,
             double tstarta, double fraca, double dMmta, int ecssna, double *sa1, double *sa2, double *sa3, int *sa4,
             double *sa5, double *sa6, double *sa7, double *sa8, double *sa9, double *sa10, double *sa11, double *sa12,
             int *sa13, double *sa14, double *sa15, int *sa16, double *sa17, double *sa18, double *sa19, int *sa20)
{ /* copy variables for singl() */

 (*sa1)=Mzamsa;
 (*sa2)=M0a;
 (*sa3)=Ma;
 (*sa4)=Ka;
 (*sa5)=tbeg;
 (*sa6)=tvira;
 (*sa7)=tenda;
 (*sa8)=La;
 (*sa9)=Ra;
 (*sa10)=Mca;
 (*sa11)=Mhea;
 (*sa12)=Mcoa;
 (*sa13)=flaga;
 (*sa14)=dta;
 (*sa15)=Mprea;
 (*sa16)=Kpa;
 (*sa17)=tstarta;
 (*sa18)=fraca;
 (*sa19)=dMmta;
 (*sa20)=ecssna;
}


void rec_dat1a(FILE *fp1, double Ma, double Mb, double *V1, double *V2, double a, double e,
              double ta_end, double tb_end, double Tmr, int Ka, int Kb,
              double MpgaA, int KpgaA, double MpgbA, int KpgbA, double apgA, double epgA, double tpgA, double MendaA,
              double MpgaB, int KpgaB, double MpgbB, int KpgbB, double apgB, double epgB, double tpgB, double MendbB,
              double Mzamsa, double Mzamsb, double a0, double e0)
{ /* prints binary data to file "data.dat" */
  /* systype() must be called prior to rec_dat1a() so Tmr is calculated in systype() and passed here */

 double Tms;                                         /* czas pomiedzy koncem ewol. skl. A i skl. B [10^6yrs] */
 double tend;                                                        /* czas konca ewolucji, utworzenie UPOZ */
 int kk;
 char his[100];

 Tms=fabs(tb_end-ta_end);                                                                       /* [10^6yrs] */
 tend=max(ta_end,tb_end);

 fprintf(fp1,"%.2f %.2f %.2f %.2f %.3f %d\n",Mzamsa,Mzamsb,a0,e0,Tst,iidd_old);                   /* starting data */
 fprintf(fp1,"%.2f %.2f %d %d %.2f %.2f %.4f %.2f\n",MpgaA,MpgbA,KpgaA,KpgbA,apgA,epgA,Tst+tpgA,MendaA);   /* SN A */
 fprintf(fp1,"%.2f %.2f %d %d %.2f %.2f %.4f %.2f\n",MpgaB,MpgbB,KpgaB,KpgbB,apgB,epgB,Tst+tpgB,MendbB);   /* SN B */
                                                                                      /* if no SN, then here zeros */
 fprintf(fp1,"%e %.3f %.3f %.3f\n",Tms,V1[0],V1[1],V1[2]);                      /* systemic velocity between 2 SNs */
 fprintf(fp1,"%e %.3f %.3f %.3f\n",Tmr,V2[0],V2[1],V2[2]);                              /* final systemic velocity */
                                                                                             /* final remnant data */
 fprintf(fp1,"%d %d %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",
         Ka,Kb,Ma,Mb,a,e,Tst+tend,Tst+ta_end,Tst+tb_end);

 fprintf(fp1,"%d %s\n",nevroute,evroute);
 fprintf(fp1,"\n");
 fflush(fp1);

 /* time [Myr], velocity [Rsun/day], Mass [Msun], distance 'a' [Rsun] */
}


void rec_dat1b(FILE *fp1, double Ma, double Mb, double *V1, double *V2, double *Vsa, double *Vsb,
               double *Vextr, double a, double e, double ta_end, double tb_end, double Tmr, int Ka, int Kb,
               double MpgaA, int KpgaA, double MpgbA, int KpgbA, double apgA, double epgA, double tpgA, double MendaA,
               double MpgaB, int KpgaB, double MpgbB, int KpgbB, double apgB, double epgB, double tpgB, double MendbB,
               double Mzamsa, double Mzamsb, double a0, double e0, double dMcea, double dMceb, int inbinA, int inbinB)
{ /* prints single data to file "sb.dat" */
  /* systype() must be called prior to rec_dat1b() so Tmr is calculated in systype() and passed here */

 double Tms;                                         /* czas pomiedzy koncem ewol. skl. A i skl. B [10^6yrs] */
 double tend;                                                        /* czas konca ewolucji, utworzenie UPOZ */
 int kk;
 char his[100];

 Tms=fabs(tb_end-ta_end);                                                                       /* [10^6yrs] */
 tend=max(ta_end,tb_end);

 fprintf(fp1,"%.2f %.2f %.2f %.2f %.3f %d\n",Mzamsa,Mzamsb,a0,e0,Tst,iidd_old);                   /* starting data */
 fprintf(fp1,"%.2f %.2f %d %d %.2e %.2f %.4f %.2f %d\n",                                                   /* SN A */
         MpgaA,MpgbA,KpgaA,KpgbA,apgA,epgA,Tst+tpgA,MendaA,inbinA);
 fprintf(fp1,"%.2f %.2f %d %d %.2e %.2f %.4f %.2f %d\n",                                                   /* SN B */
         MpgaB,MpgbB,KpgaB,KpgbB,apgB,epgB,Tst+tpgB,MendbB,inbinB);
                                                                                      /* if no SN, then here zeros */
 fprintf(fp1,"%e %.3f %.3f %.3f\n",Tms,V1[0],V1[1],V1[2]);                      /* systemic velocity between 2 SNs */
 fprintf(fp1,"%e %.3f %.3f %.3f\n",Tmr,V2[0],V2[1],V2[2]);                              /* final systemic velocity */

 fprintf(fp1,"%.3f %.3f %.3f\n",V1[0]+Vsa[0],V1[1]+Vsa[1],V1[2]+Vsa[2]);  /* component A velocity after disruption */
 fprintf(fp1,"%.3f %.3f %.3f\n",V1[0]+Vsb[0],V1[1]+Vsb[1],V1[2]+Vsb[2]);  /* component B velocity after disruption */
 fprintf(fp1,"%.3f %.3f %.3f\n",Vextr[0],Vextr[1],Vextr[2]);         /* final velocity of single star after 2nd SN */

                                                                                             /* final remnant data */
 fprintf(fp1,"%d %d %.3f %.3f %.3e %.3f %.3f %.3f %.3f %.3f %.3f\n",
         Ka,Kb,Ma,Mb,a,e,Tst+tend,Tst+ta_end,Tst+tb_end,dMcea,dMceb);

 fprintf(fp1,"%d %s\n",nevroute,evroute);
 fprintf(fp1,"\n");
 fflush(fp1);

 /* time [Myr], velocity [Rsun/day], Mass [Msun], distance 'a' [Rsun] */
}


void rec_dat2(FILE *fp2, int j, int aicns, int aicbh, int badorb)
{ /* dump sn counters to info.dat */
 int i;

 fprintf(fp2,"aicns: %d, aicbh: %d\n",aicns,aicbh);
 fprintf(fp2,"badorb: %d (bad orbital solution in orb_change(): systems skipped\n",badorb);
 for(i=0;i<20;i++)
   fprintf(fp2,"sIbc: %d ",sIbc[i]);
 fprintf(fp2,"\n");
 for(i=0;i<20;i++)
   fprintf(fp2,"sII:  %d ",sII[i]);
 fprintf(fp2,"\n");
 for(i=0;i<20;i++)
   fprintf(fp2,"sFa:  %d ",sFa[i]);
 fprintf(fp2,"\n\n");

 fprintf(fp2,"  j: %d (Total number of tested binaries)\n\n\n",j);
}


double min3(double x1, double x2, double x3)
{ /* returns minimum of the three argumens */

 if(x1<x2 && x1<x3) return x1;
 else if(x2<x3) return x2;
 else return x3;
}


/*-------------------------------------- GETTING INITIAL VALUES ----------------------------*/

double get_t(double Min, double Max)
{ /* obtain hubble time=evolutionary time from flat distribution */
  /* to simulate constant star formation rate */

 if(Random==1)
   return fabs(Max-Min)*(double)ran3(&idum)+Min;
 else
   return fabs(Max-Min)*(double)ran2(&idum)+Min;
}


double get_q1(double Min, double Max)
{ /* zadajemy plaski losowy rozklad q, pomiedzy [Min,Max] */
 double q;

 if(Random==1)
   q=fabs(Max-Min)*(double)ran3(&idum)+Min;
 else
   q=fabs(Max-Min)*(double)ran2(&idum)+Min;

 return q;
}


double get_q2(double Min, double Max)
{ /* q drawn from distribution P(q)=A*q^(Rat) in range q: [Min,Max] */
  /* Global Variable Rat -exponent of q-distribution, Rat>0.0 works fine */
  /* for Rat<0.0, let say Rat=-2.7 q drops down so fast that Mb almost always <0.5 Msun */
  /* use then get_q3() function */
  /* It never shall be Rat=-1.0 as then different integrations are needed */
 double A,q0,P0,tmp;

 if(fabs(Min)<acc) Min=0.000001;   /* if M==0 and Rat<0 then pow(Min,tmp) would be NaN */

 if(fabs(Rat-1.0)<acc)
   fprintf(fp0,"error: in get_q2() wrong integration for exponent Rat=-1.0!!!");

 tmp=Rat+1.0;
 A=tmp/(pow(Max,tmp)-pow(Min,tmp));
 if(Random==1)
   P0=(double)ran3(&idum);
 else
   P0=(double)ran2(&idum);
 q0=pow(tmp*P0/A+pow(Min,tmp),1.0/tmp);

 return q0;
}


double get_q3(double qmin)
{ /* distribution is ~c2*q^{Rat} for q>0.2, for qmin<=q<=0.2 distribution is flat at the level of c1*/
  /* to match ~c2*q^{Rat} distribution at q=0.2, prepared for Rat<0.0 */
 double P0,tmp,limit,c1,c2;

 if(fabs(Rat-1.0)<acc) {
   fprintf(fp0,"error: in get_q2() wrong integration for exponent Rat=-1.0!!!");
   fflush(fp0);
 }

 tmp=Rat+1.0;
 c2=1.0/(pow(0.2,Rat)*(0.2-qmin)+(1.0-pow(0.2,tmp))/tmp);
 c1=c2*pow(0.2,Rat);
 limit=c1*(0.2-qmin);

 if(Random==1)
   P0=(double)ran3(&idum);          /* P0: [0:1] */
 else
   P0=(double)ran2(&idum);

 if(P0<=limit)                    /* flat distribution normalized to c1: f(q)=c1 */
   return (P0/c1)+qmin;
 else
   return pow(tmp*(P0-c1*(0.2-qmin))/c2+pow(0.2,tmp),1.0/tmp);
}


double get_SS(double Min, double Max, double Ex)
{ /* q drawn from distribution P(q)=A*q^(Ex) in range q: [Min,Max] */
  /* It never shall be Ex=-1.0 as then different integrations are needed */
 double A,q0,P0,tmp;

 if(fabs(Min)<acc) Min=0.000001;   /* if M==0 and Ex<0 then pow(Min,tmp) would be NaN */

 if(fabs(Ex-1.0)<acc)
   fprintf(fp0,"error: in get_SS() wrong integration for exponent Ex=-1.0!!!");

 tmp=Ex+1.0;
 A=tmp/(pow(Max,tmp)-pow(Min,tmp));
 if(Random==1)
   P0=(double)ran3(&idum);
 else
   P0=(double)ran2(&idum);
 q0=pow(tmp*P0/A+pow(Min,tmp),1.0/tmp);

 return q0;
}

double get_M_linear(double Min, double Max) {
    /* M drawn from a linear (uniform) distribution in the range [Min:Max] */
    /* Assumes Min < Max */

    // Generate a random double between 0.0 and 1.0 (exclusive of 1.0)
    // Using rand() for simplicity. For higher quality randomness, consider other methods.
    double u = (double)rand() / RAND_MAX;

    // Scale and shift to the desired range [Min, Max)
    double M = Min + u * (Max - Min);

    return M;
}

double get_M(double Min, double Max)
{ /* M drawn from merger of three distributions (Kroupa, Taut Gilmore 1993, MNRAS 262, 545): */
  /* M drawn from the range [Min:Max], Min must be >=0.08, Max must be >1.0 */
  /*  P(M)=A1*M^(-1.3) in range M: [0.08,0.5] */
  /*  P(M)=A2*M^(-2.2) in range M: [0.5,1.0] */
  /*  P(M)=A3*M^(Sal) in range M: [1.0,Max],  e.g., Max=100 Msun */
  /* Sal=-2.7 or even higher (-3.2) for field stars, Sal=-2.35 for cluster stars */
  /* A1,A2,A3 choosen so the distributions meet at 0.5 and 1.0, and their integral is = 1.0 */
  /* Min value achived through draw repetition until M is higher than Min */
 double m,n,a,b,c,A1,A2,A3,tot;
 double P0,M,tmp;

 if(Min<(0.08-acc) || Max<(1.0+acc) || Min>=Max) {
   fprintf(fp0,"error: in get_M() wrong choice for either Min or/and Max value\n");
   fflush(fp0);
 }
 if(fabs(Sal-1.0)<acc) {
   fprintf(fp0,"error: in get_M() wrong integration for IMF exponent Sal=-1.0!!!");
   fflush(fp0);
 }

 tmp=Sal+1.0;
 m=pow(0.5,-1.3)/pow(0.5,-2.2);
 n=pow(1.0,-2.2)/pow(1.0,Sal);
 a=(1.0/(-0.3))*(pow(0.5,-0.3)-pow(0.08,-0.3));
 b=(1.0/(-1.2))*(pow(1.0,-1.2)-pow(0.5,-1.2));
 c=(1.0/(tmp))*(pow(Max,tmp)-pow(1.0,tmp));
 A1=1.0/(a+m*b+n*m*c);
 A2=m*A1;
 A3=m*n*A1;
 tot=A1*a+A2*b+A3*c;

 M=-1.0;
 while(M<Min) {
   if(Random==1)
     P0=(double)ran3(&idum);
   else
     P0=(double)ran2(&idum);

   if(P0<=(A1*a))
     M=pow(-0.3*P0/A1+pow(0.08,-0.3),-1.0/0.3);
   else if(P0<=(A1*a+A2*b))
     M=pow((P0-A1*a)*(-1.2)/A2+pow(0.5,-1.2),1.0/-1.2);
   else
     M=pow((P0-A1*a-A2*b)*(tmp)/A3+pow(1.0,tmp),1.0/tmp);
 }

 return M;
}


double get_a(double Min, double Max)
{ /* losujemy log10(a) z rozkladu plaskiego */
 double tmp1,tmp2,tmp3,a0;

 tmp1=log10(Min);
 tmp2=log10(Max);

 if(Random==1)
   tmp3=fabs(tmp2-tmp1)*(double)ran3(&idum)+tmp1;
 else
   tmp3=fabs(tmp2-tmp1)*(double)ran2(&idum)+tmp1;
 a0=pow(10.0,tmp3);

 return a0;
}


double get_e(double Min, double Max)
{ /* rozklad f(e)=2e */
 double P0,e0;

 if(Random==1)
   P0=(double)ran3(&idum);
 else
   P0=(double)ran2(&idum);
 e0=sqrt(P0);

 return e0;
}


double get_i(double Min, double Max)
{ /* rozklad plaski, zwraca 'i' w radianach */

 if(Random==1)
   return fabs(Max-Min)*(double)ran3(&idum)+Min;
 else
   return fabs(Max-Min)*(double)ran2(&idum)+Min;
}


double get_Om(double Min, double Max)
{ /* rozklad plaski */

 if(Random==1)
   return fabs(Max-Min)*(double)ran3(&idum)+Min;
 else
   return fabs(Max-Min)*(double)ran2(&idum)+Min;
}


double get_om(double Min, double Max)
{ /* rozklad plaski */

 if(Random==1)
   return fabs(Max-Min)*(double)ran3(&idum)+Min;
 else
   return fabs(Max-Min)*(double)ran2(&idum)+Min;
}


double get_texp(double Min, double Max)
{ /* rozklad plaski  */

 if(Random==1)
   return fabs(Max-Min)*(double)ran3(&idum)+Min;
 else
   return fabs(Max-Min)*(double)ran2(&idum)+Min;
}


void get_Vkick1(double *Vkick)
{ /* losowanie ze zlozenia 2 gausjanow; 1 z Sigma1(=175 km/s) dla (Divide)*100% kickow,
     2 z Sigma2(=700 km/s) dla (1-Divide)100% kickow, returns in units of [Rsun/day] */
 double sig,xacc,tmp1,P0;
 double x1,x2,f1,f2,x0;
 int i;


 if(Random==1)
   tmp1=(double)ran3(&idum);           /* tmp1: [0:1] */
 else
   tmp1=(double)ran2(&idum);
 if(tmp1<=Divide) sig=Sigma1;
 else sig=Sigma2;

 for(i=0;i<3;i++) {
   if(Random==1)
     P0=(double)ran3(&idum);             /* P0: [0:1] */
   else
     P0=(double)ran2(&idum);
   x1=-30.0;
   x2=30.0;
                                         /* granice [x1:x2] w jakich szukamy roz. dla */
                                  /* naszego gausjanu granica[sigmy]=-/+x2*(1/sqrt(2))*/
   f1=-0.5*erffNR(-x1)+0.5;
   f2=0.5*erffNR(x2)+0.5;
   if(P0<=f1) {x0=f1;}
   else if(P0>=f2) {x0=f2;}
   else {
     xacc=(1.0e-06)*(fabs(x1)+fabs(x2))/2.0;
     x0=rtbis1(func1,x1,x2,xacc,P0);
   }
   (*(Vkick+i))=sqrt(2.0)*sig*x0;
   (*(Vkick+i))*=tr2;
 }

}


void get_Vkick2(double *Vkick, double fraca, int fba)
{ /* losowanie ze zlozenia 2 gausjanow; 1 z Sigma1(=175 km/s) dla (Divide)*100% kickow, */
  /* 2 z Sigma2(=700 km/s) dla (1-Divide)100% kickow, returns in units of [Rsun/day] */
  /* No Fall Back -- full kick, Complete FB -- no kick */
  /* partial FB kick decreased propotionaly to amount of falling back material fraca */
  /* fraca passed form singl() */
 double sig,xacc,tmp1,P0;
 double x1,x2,f1,f2,x0,lower;
 int i;

 if(Random==1)
   tmp1=(double)ran3(&idum);           /* tmp1: [0:1] */
 else
   tmp1=(double)ran2(&idum);
 if(tmp1<=Divide) sig=Sigma1;
 else sig=Sigma2;

 for(i=0;i<3;i++) {
   if(Random==1)
     P0=(double)ran3(&idum);             /* P0: [0:1] */
   else
     P0=(double)ran2(&idum);
   x1=-30.0;
   x2=30.0;
                                           /* granice [x1:x2] w jakich szukamy roz. dla */
                                    /* naszego gausjanu granica[sigmy]=-/+x2*(1/sqrt(2))*/
   f1=-0.5*erffNR(-x1)+0.5;
   f2=0.5*erffNR(x2)+0.5;
   if(P0<=f1) {x0=f1;}
   else if(P0>=f2) {x0=f2;}
   else {
     xacc=(1.0e-06)*(fabs(x1)+fabs(x2))/2.0;
     x0=rtbis1(func1,x1,x2,xacc,P0);
   }
   (*(Vkick+i))=sqrt(2.0)*sig*x0;
   (*(Vkick+i))*=tr2;

   if(fba==0) lower=0.0;                     /* full kick */
   else if(fba==1 || fba==2) lower=fraca;    /* decreased kick */
   else fprintf(fp0,"error: in get_Vkick2() wrong fba: %d\n",fba);

   (*(Vkick+i))*=(1.0-lower);       /* No FB full kick, Complete FB no kick, partial FB kick decreased */
 }                                            /* propotionaly to amount of falling back material fraca */

}


void get_Vkick3(double *Vkick, int Ka)
{ /* losowanie ze zlozenia 2 gausjanow; 1 z Sigma1(=175 km/s) dla (Divide)*100% kickow, */
  /* 2 z Sigma2(=700 km/s) dla (1-Divide)100% kickow, returns in units of [Rsun/day] */
  /* if oubursting star is black hole Ka=14 then kick is assumed to be 0.0 */
 double sig,xacc,tmp1,P0;
 double x1,x2,f1,f2,x0;
 int i;

 if(Random==1)
   tmp1=(double)ran3(&idum);           /* tmp1: [0:1] */
 else
   tmp1=(double)ran2(&idum);
 if(tmp1<=Divide) sig=Sigma1;
 else sig=Sigma2;

 for(i=0;i<3;i++) {
   if(Ka==13) {
     if(Random==1) {
       P0=(double)ran3(&idum);             /* P0: [0:1] */
     }
     else {
       P0=(double)ran2(&idum);
     }
     x1=-30.0;
     x2=30.0;
                                                 /* granice [x1:x2] w jakich szukamy roz. dla */
                                          /* naszego gausjanu granica[sigmy]=-/+x2*(1/sqrt(2))*/
     f1=-0.5*erffNR(-x1)+0.5;
     f2=0.5*erffNR(x2)+0.5;
     if(P0<=f1) {x0=f1;}
     else if(P0>=f2) {x0=f2;}
     else {
       xacc=(1.0e-06)*(fabs(x1)+fabs(x2))/2.0;
       x0=rtbis1(func1,x1,x2,xacc,P0);
     }
     (*(Vkick+i))=sqrt(2.0)*sig*x0;
     (*(Vkick+i))*=tr2;
   }
   else if(Ka==14) {                      /* if compact object is BH, then its kick is ZERO */
     (*(Vkick+i))=0.0;
   }
 }

}


void get_Vkick4(double *Vkick)
{ /* No kicks */
 int i;

 for(i=0;i<3;i++)
   (*(Vkick+i))=0.0;
}


void get_Vkick5(double *Vkick, int Ka, double fraca, int *cmt3, int fba)
{ /* Paczynski 1990 (ApJ 348, 485 kicks with sigma=(V^*)=600 km/s (instead of original 270 km/s) */
  /* sigma =600 km/s suggested by Hartmann 97 (Vicky comment) */
  /* I draw value of kick u from eq.(3), and then draw random kick direction: */
  /* geting from uniform distribution fi,cos(teta) to get uniform distribution over sphere, */
  /* and thus I obtain u_x,u_y,u_z form Paczynski distribution */
  /* needs global xPacz[],yPacz defined in sinbin.h */
  /* if oubursting star is black hole Ka=14 then kick is decreased */
 double teta,fi,tmp1;
 double u,v,sigma,pu,P0,del,lower;
 int i,j,NPacz;

 if((*cmt3)==0) {                               /* first time in get_Vkick5 so fills xPacz,yPacz */
   (*cmt3)=1;
   j=0;
   del=0.01;
   u=0.0;                                       /* make sure that sizes of xPacz,yPacz are at least Npacz */
   while(u<=101.0) {                            /* this 101 ensures that P0 from [0,Max=100] will be braceted */
     if(u>5.0) {del=0.1;}                       /* in interp() below */
     pu=2.0*(u/(1.0+u*u)+atan(u))/(Pi);
     xPacz[j]=u;
     yPacz[j]=pu;
     u+=del;
     j++;
   }
 }
 NPacz=1460;                                    /* how many points in tables xPacz,yPacz -- set by del */

 sigma=600.0;                                   /* (V^*) from paczynski eq.(3) */

 if(Random==1)
   fi=2.0*Pi*(double)ran3(&idum);                 /* fi: [0:2*Pi] */
 else
   fi=2.0*Pi*(double)ran2(&idum);
 if(Random==1)
   tmp1=2.0*(double)ran3(&idum)-1.0;              /* cos(teta): [-1:1] */
 else
   tmp1=2.0*(double)ran2(&idum)-1.0;
 teta=acos(tmp1);


 if(Random==1)
   P0=(double)ran3(&idum);                        /* random uniform number from [0,1] */
 else
   P0=(double)ran2(&idum);
 u=interp(P0,xPacz,yPacz,NPacz);                /* xPacz,yPacz store dystribuant eq.(7) of Pacz. dystribution */
 v=u*sigma;

 (*(Vkick+0))=v*sin(teta)*cos(fi);
 (*(Vkick+1))=v*sin(teta)*sin(fi);
 (*(Vkick+2))=v*cos(teta);

 for(i=0;i<3;i++) {
   (*(Vkick+i))*=tr2;

   if(fba==0) lower=0.0;                     /* full kick */
   else if(fba==1 || fba==2) lower=fraca;    /* decreased kick */
   else fprintf(fp0,"error: in get_Vkick5() wrong fba: %d\n",fba);

   (*(Vkick+i))*=(1.0-lower);       /* No FB full kick, Complete FB no kick, partial FB kick decreased */
 }                                            /* propotionaly to amount of falling back material fraca */

}


void get_Vkick6(double *Vkick, double fraca, int fba)
{ /* single Maxwellian Hobbs et al. distribution with Sigma3=265 km/s (3d speed): func_max2() */
  /* returns three components (Vx,Vy,Vz) in units of [Rsun/day] */
  /* kick==6: No Fall Back -- full kick, Complete FB -- no kick */
  /* partial FB kick decreased propotionaly to amount of falling back material fraca */
  /* kick==7: no natal NS/BH kick */
  /* fraca passed form singl() */
 double xacc,P0,x1,x2,f1,f2,x0;
 double cteta,teta,fi,lower;
 int i;


 if(fabs(Sigma3)<0.001) {                        /* zero kick model */
   for(i=0;i<3;i++)
      (*(Vkick+i))=0.0;
 }
 else {
   if(Random==1)
     P0=(double)ran3(&idum);                     /* P0: [0:1] flat ditr. */
   else
     P0=(double)ran2(&idum);

   x1=0.0000001;
   x2=1000000.0;
   f1=func_max2(x1,0.0);
   f2=func_max2(x2,0.0);
   xacc=0.0001;

   if(P0<=f1+acc) x0=0.0;
   else if(P0>=f2-acc) x0=1000000.0;
   else x0=rtbis1(func_max2,x1,x2,xacc,P0);      /* magnitude of kick set */

   cteta=get_flat(-1.0,1.0);                     /* cos(teta) from flat distr. */
   teta=acos(cteta);
   fi=get_flat(-Pi,Pi);                          /* fi from flat distr, and direction of kick set */

   (*(Vkick+0))=x0*sin(teta)*cos(fi);            /* components of kick velocity set */
   (*(Vkick+1))=x0*sin(teta)*sin(fi);
   (*(Vkick+2))=x0*cos(teta);

   for(i=0;i<3;i++) {
     (*(Vkick+i))*=tr2;                          /* [km/s] -> [Rsun/day] */

     if(fba==0) lower=0.0;                       /* full kick */
     else if(fba==1 || fba==2) lower=fraca;      /* decreased kick */
     else fprintf(fp0,"error: in get_Vkick6() wrong fba: %d (%d %d)\n",fba,idum_run,iidd_old);

     if(kick==6)                                 /* No FB: full kick; full FB: no kick; partial FB: kick decreased */
       (*(Vkick+i))*=(1.0-lower);                /* propotionaly to amount of falling back material: fraca for NS/BH */
     else
       ;                                         /* kick==7: no change of kick by fall back: full natal NS/BH kick from Sigma3 */
   }
 }
}


void get_Vkick7(double *Vkick, double Mf)
{ /* single Maxwellian Hobbs et al. distribution with Sigma3=265 km/s (3d speed): func_max2() */
  /* returns three components (Vx,Vy,Vz) in units of [Rsun/day] */
  /* NS kicks directly from Sigma3, BH kicks: lowered proportionaly to BH mass */
  /* Mf is final NS or BH mass. Adopted from Rodriguez et al. 2016, Apj 832, L2: eq.2 */
 double xacc,P0,x1,x2,f1,f2,x0;
 double cteta,teta,fi,lower;
 int i;


 if(fabs(Sigma3)<0.001) {                        /* zero kick model */
   for(i=0;i<3;i++)
      (*(Vkick+i))=0.0;
 }
 else {
   if(Random==1)
     P0=(double)ran3(&idum);                     /* P0: [0:1] flat ditr. */
   else
     P0=(double)ran2(&idum);

   x1=0.0000001;
   x2=1000000.0;
   f1=func_max2(x1,0.0);
   f2=func_max2(x2,0.0);
   xacc=0.0001;

   if(P0<=f1+acc) x0=0.0;
   else if(P0>=f2-acc) x0=1000000.0;
   else x0=rtbis1(func_max2,x1,x2,xacc,P0);      /* magnitude of kick set */

   cteta=get_flat(-1.0,1.0);                     /* cos(teta) from flat distr. */
   teta=acos(cteta);
   fi=get_flat(-Pi,Pi);                          /* fi from flat distr, and direction of kick set */

   (*(Vkick+0))=x0*sin(teta)*cos(fi);            /* components of kick velocity set */
   (*(Vkick+1))=x0*sin(teta)*sin(fi);
   (*(Vkick+2))=x0*cos(teta);

   for(i=0;i<3;i++)
     (*(Vkick+i))*=tr2;                          /* [km/s] -> [Rsun/day] */

   if(Mf>Mmaxns)                                  /* BH */
     for(i=0;i<3;i++)
       (*(Vkick+i))=(2.5/Mf)*(*(Vkick+i));       /* [use of 2.5 maximizies BH kick; use of 1.4 would lower BH kick */
 }
}


void get_Vkick8(double *Vkick, double Msn, double Mf)
{ /* returns three components (Vx,Vy,Vz) in units of [Rsun/day]; does not use Sigma/Maxwelian */
  /* NS and BH kicks: lowered proportionaly to NS/BH mass */
  /* Msn: mass of star at the time of core-collpse/SN, Mf is final NS or BH mass */
  /* adopted from Bray & Eldridge 2016, MNRAS 461, 3747: eq.1 -- done based on NS/pulsar velocities */
 double xacc,P0,x1,x2,f1,f2,x0;
 double cteta,teta,fi,lower;
 double Mejecta;
 int i;


 if(Random==1)
   P0=(double)ran3(&idum);                     /* P0: [0:1] flat ditr. */
 else
   P0=(double)ran2(&idum);

 Mejecta=Msn-Mf;                               /* SN ejected mass [Msun] */
 x0=alpha1*(Mejecta/Mf)+beta1;                 /* magnitude of kick [km/s] */
 if(x0<0.0) x0=0.0;

 cteta=get_flat(-1.0,1.0);                     /* cos(teta) from flat distr. */
 teta=acos(cteta);
 fi=get_flat(-Pi,Pi);                          /* fi from flat distr, and direction of kick set */

 (*(Vkick+0))=x0*sin(teta)*cos(fi);            /* components of kick velocity set */
 (*(Vkick+1))=x0*sin(teta)*sin(fi);
 (*(Vkick+2))=x0*cos(teta);

 for(i=0;i<3;i++)
   (*(Vkick+i))*=tr2;                          /* [km/s] -> [Rsun/day] */
}


void get_Vkick9(double *Vkick, double fraca, int fba)
{ /* kick==10: fall-back (mass) decreased kick with Sigma4 and full (neutrino) kick with Sigma5 */
  /* 2 single Maxwellians: one with Sigma4 and one with Sigma5: set in sinbin.h */
  /* employs: func_max2a() and func_max2b() */
  /* returns three components (Vx,Vy,Vz) of the combined kick in units of [Rsun/day] */
  /* amount of falling back material fraca passed form singl() */
 double xacc,P0,x1,x2,f1,f2,x0;
 double cteta,teta,fi,lower;
 double V1[3],V2[3];
 int i;


 if(fabs(Sigma4)<0.001) {                        /* pick asymmetric mass ejecttion (fall-back decreased) kick */
   for(i=0;i<3;i++)                              /* first random orientation */
      V1[i]=0.0;
 }
 else {

   if(Random==1)
     P0=(double)ran3(&idum);                     /* P0: [0:1] flat ditr. */
   else
     P0=(double)ran2(&idum);

   x1=0.0000001;
   x2=1000000.0;
   f1=func_max2a(x1,0.0);
   f2=func_max2a(x2,0.0);
   xacc=0.0001;

   if(P0<=f1+acc) x0=0.0;
   else if(P0>=f2-acc) x0=1000000.0;
   else x0=rtbis1(func_max2a,x1,x2,xacc,P0);      /* magnitude of kick set */

   cteta=get_flat(-1.0,1.0);                     /* cos(teta) from flat distr. */
   teta=acos(cteta);
   fi=get_flat(-Pi,Pi);                          /* fi from flat distr, and direction of kick set */

   V1[0]=x0*sin(teta)*cos(fi);                   /* components of kick velocity set */
   V1[1]=x0*sin(teta)*sin(fi);
   V1[2]=x0*cos(teta);

   for(i=0;i<3;i++) {
     if(fba==0) lower=0.0;                       /* full kick */
     else if(fba==1 || fba==2) lower=fraca;      /* decreased kick */
     else fprintf(fp0,"error: in get_Vkick9() wrong fba: %d (%d %d)\n",fba,idum_run,iidd_old);
     V1[i]*=tr2;                                 /* [km/s] -> [Rsun/day] */
     V1[i]*=(1.0-lower);                         /* decrease kick magnitude: fraca for NS/BH */
   }
 }


 if(fabs(Sigma5)<0.001) {                        /* pick asymmetric neutrino emission (mass-independent) kick */
   for(i=0;i<3;i++)                              /* second random orientation */
      V2[i]=0.0;
 }
 else {

   if(Random==1)
     P0=(double)ran3(&idum);                     /* P0: [0:1] flat ditr. */
   else
     P0=(double)ran2(&idum);

   x1=0.0000001;
   x2=1000000.0;
   f1=func_max2b(x1,0.0);
   f2=func_max2b(x2,0.0);
   xacc=0.0001;

   if(P0<=f1+acc) x0=0.0;
   else if(P0>=f2-acc) x0=1000000.0;
   else x0=rtbis1(func_max2b,x1,x2,xacc,P0);     /* magnitude of kick set */

   cteta=get_flat(-1.0,1.0);                     /* cos(teta) from flat distr. */
   teta=acos(cteta);
   fi=get_flat(-Pi,Pi);                          /* fi from flat distr, and direction of kick set */

   V2[0]=x0*sin(teta)*cos(fi);                   /* components of kick velocity set */
   V2[1]=x0*sin(teta)*sin(fi);
   V2[2]=x0*cos(teta);

   for(i=0;i<3;i++)                              /* full kick */
     V2[i]*=tr2;                                 /* [km/s] -> [Rsun/day] */
 }

 (*(Vkick+0))=V1[0]+V2[0];                       /* components of combined (mass+neutrino) kick velocity */
 (*(Vkick+1))=V1[1]+V2[1];
 (*(Vkick+2))=V1[2]+V2[2];
}


double func_max1(double v)
{ /* returns value of Maxwellian */

 if(fabs(v)>10.0*Sigma3)     /* for huge kicks return 0.0 by hand, otherwise */
   return 0.0;               /* the return value is larger than allocated memory for double  */
 else
   return sqrt(2.0/Pi)*v*v*exp(-v*v/(2.0*Sigma3*Sigma3))/pow(Sigma3,3.0);
}

double func_max1a(double v)
{ /* returns value of Maxwellian */

 if(fabs(v)>10.0*Sigma4)     /* for huge kicks return 0.0 by hand, otherwise */
   return 0.0;               /* the return value is larger than allocated memory for double  */
 else
   return sqrt(2.0/Pi)*v*v*exp(-v*v/(2.0*Sigma4*Sigma4))/pow(Sigma4,3.0);
}

double func_max1b(double v)
{ /* returns value of Maxwellian */

 if(fabs(v)>10.0*Sigma5)     /* for huge kicks return 0.0 by hand, otherwise */
   return 0.0;               /* the return value is larger than allocated memory for double  */
 else
   return sqrt(2.0/Pi)*v*v*exp(-v*v/(2.0*Sigma5*Sigma5))/pow(Sigma5,3.0);
}

double func_max2(double mtop, double P0)
{ /* returns value of area of Maxwelian at mtop minus P0 */
  double a1;

  a1=qtrap(func_max1,0.000000001,mtop);

  return a1-P0;
}

double func_max2a(double mtop, double P0)
{ /* returns value of area of Maxwelian at mtop minus P0 */
  double a1;

  a1=qtrap(func_max1a,0.000000001,mtop);

  return a1-P0;
}

double func_max2b(double mtop, double P0)
{ /* returns value of area of Maxwelian at mtop minus P0 */
  double a1;

  a1=qtrap(func_max1b,0.000000001,mtop);

  return a1-P0;
}


double get_flat(double Min, double Max)
{ /* returns random number in range: [Min,Max] from flat probablility distr. */
 double flat;

 if(Random==1)
   flat=fabs(Max-Min)*(double)ran3(&idum)+Min;
 else
   flat=fabs(Max-Min)*(double)ran2(&idum)+Min;

 return flat;
}


double interp(double p, double xPacz[], double yPacz[], int NPacz)
{ /* znajduje argument dystrybunaty zapisany w xPacz[],yPacz[] dla danej wartosci p */
  /* NPacz -- dokad zapisane sa tablice xPacz[],yPacz[] */
 double x1,y1,a,b,up,tmp1,tmp2;
 int mark1=0;
 int mark2=0;
 int i;

 for(i=0;i<NPacz;i++) {
   if(yPacz[i]>=p && mark1==0) {
     fprintf(fp0,"error: in interp() p is below available distribuant\n");
     fflush(fp0);
     exit(-1);
   }
   if(yPacz[i]>=p) {
     mark2=1;
     break;
   }
   tmp1=xPacz[i];
   tmp2=yPacz[i];
   mark1=1;
 }

 if(mark2==0) {
   fprintf(fp0,"error: in interp() p is over available distribuant\n");
   fflush(fp0);
   exit(-1);
 }

 x1=tmp1;                                    /* linear interpolation */
 y1=tmp2;
 a=(yPacz[i]-y1)/(xPacz[i]-x1);
 b=yPacz[i]-a*xPacz[i];
 up=(p-b)/a;

 return up;
}



/*-------------------------- FROM NUMERICAL RECIPES -------------------------*/

double rtbis1(double (*func)(double, double), double x1, double x2,
              double xacc, double P0)
{ /* przepisane rtbis.c z metod Num. zmiany float na double, i postac func */
 int j;
 double dx,f,fmid,xmid,rtb;

 f=(*func)(x1,P0);
 fmid=(*func)(x2,P0);
 if(f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis1");
 rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
 for(j=1;j<=JMAX;j++) {
   fmid=(*func)(xmid=rtb+(dx *= 0.5),P0);
   if(fmid <= 0.0) rtb=xmid;
   if(fabs(dx) < xacc || fmid == 0.0) return rtb;
 }

 nrerror("Too many bisections in rtbis1");
 return 0.0;
}


double func1(double x, double P0)
{ /* zwraca wartosc calki gausjanu[o sig=1/sqrt(2.0)] od [-nies:x] minus P0 */

 if(x<0.0)
   return 0.5-0.5*erffNR(-x)-P0;
 else
   return 0.5+0.5*erffNR(x)-P0;
}


int zbrac2(double (*func)(double, double, double, double, double, int, int, double), double *x1, double *x2,
          double Ma, double Mb, double Ia, double Ib, int Ka, int Kb, double Ac)
{ /* tries to bracket the root of "a" in Jtot(final)-Jtot(initial)=0 equation */
 int ntry=50;
 double factor=1.6;
 int j;
 double f1,f2;

 if (*x1 == *x2) nrerror("Bad initial range in zbrac");
 if(*x1<acc || *x2<acc) return 0;	/* no root found, and "a" become negative! */
 f1=(*func)(*x1,Ma,Mb,Ia,Ib,Ka,Kb,Ac);
 f2=(*func)(*x2,Ma,Mb,Ia,Ib,Ka,Kb,Ac);
 for (j=1;j<=ntry;j++) {
   if(*x1<acc || *x2<acc) return 0;    /* no root found, and "a" become negative! */
   if (f1*f2 < 0.0) return 1;
   if (fabs(f1) < fabs(f2))
     f1=(*func)(*x1 += factor*(*x1-*x2),Ma,Mb,Ia,Ib,Ka,Kb,Ac);
   else
     f2=(*func)(*x2 += factor*(*x2-*x1),Ma,Mb,Ia,Ib,Ka,Kb,Ac);
 }

 return 0;
}


double rtbis2(double (*func)(double, double, double, double, double, int, int, double), double x1, double x2,
              double xacc, double Ma, double Mb, double Ia, double Ib, int Ka, int Kb, double Ac)
{ /* przepisane rtbis.c z metod Num. zmiany float na double, i postac func */
 int j;
 double dx,f,fmid,xmid,rtb;

 f=(*func)(x1,Ma,Mb,Ia,Ib,Ka,Kb,Ac);
 fmid=(*func)(x2,Ma,Mb,Ia,Ib,Ka,Kb,Ac);
 if(f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis2");
 rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
 for(j=1;j<=JMAX;j++) {
   fmid=(*func)(xmid=rtb+(dx *= 0.5),Ma,Mb,Ia,Ib,Ka,Kb,Ac);
   if(fmid <= 0.0) rtb=xmid;
   if(fabs(dx) < xacc || fmid == 0.0) return rtb;
 }

 nrerror("Too many bisections in rtbis2");
 return 0.0;
}


double func2(double x, double Ma, double Mb, double Ia, double Ib, int Ka, int Kb, double Ac)
{ /* angular momentum conservation for binary: returns Jtot(initial)-Jtot(final) */
  /* for given semi-major axix "a" (which here is called "x") */

  if(Ka>=10) Ia=0.0;   /* Ia is not zero, but it mimics the wa=0.0, which is true for WD/NS/BH */
  if(Kb>=10) Ib=0.0;   /* we have assumed zero spins for remnants, and can not have them synchronized */
                       /* in the equation below: Jorb+Ja+Jb-Ac=Jorb+Ia*wa+Ib*wb-Ac, where for all non */
                       /* remnants donors it is assumed wa=wb=worb, but now we get solution for */
                       /* w=0.0 for K>=10 with means of Ia */

  return Ma*Mb*sqrt(GGG*x/(Ma+Mb))+sqrt(GGG*(Ma+Mb))*(Ia+Ib)*pow(x,-1.5)-Ac;
}


double erffNR(double x)
{
	return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}


double gammp(double a, double x)
{
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammp");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}


double gammq(double a, double x)
{
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}


void gcf(double *gammcf, double a, double x, double *gln)
{
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}


void gser(double *gamser, double a, double x, double *gln)
{
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}


double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(fp0,"Not only Numerical Recipes run-time error...\n");
  fprintf(fp0,"%s\n",error_text);
  fprintf(fp0,"...now exiting to system...\n");
  exit(1);
}


float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

        iidd++;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}


float ran3(long *idum)
{ /* z metod numerycznych */
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

        iidd++;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
                *idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}


double func_tmr(double e)
{ /* returns value of function which needs to be integrated in tmerge() */

 return pow(e,29.0/19.0)*pow(1.0+121.0*e*e/304.0,1181.0/2299.0)/pow(1-e*e,1.5);
}


double qtrap(double (*func)(double), double a, double b)
{ /* 'Numerical Recipes in C': needed for integration in tmerge() */
  /* WARNING: if qtrap() is unable to integrate in JMAX3 steps then it returns -1!!! */

 int j;
 double s,olds;

 olds = -1.0e30;
 for (j=1;j<=JMAX3;j++) {
   s=trapzd(func,a,b,j);
   if (fabs(s-olds) < EPS3*fabs(olds)) return s;
   olds=s;
 }
 fprintf(fp0,"iidd_old: %d, Too many steps in routine qtrap (tmerge=-1.0)",iidd_old);
 fflush(fp0);
 return -1.0;                                           /* my change!!!!!! */
}


double trapzd(double (*func)(double), double a, double b, int n)
{ /* 'Numerical Recipes in C': needed for integration in tmerge() */
 double x,tnm,sum,del;
 static double s;
 int it,j;

 if(n==1) {
   return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
 }
 else {
   for (it=1,j=1;j<n-1;j++) it <<= 1;
   tnm=it;
   del=(b-a)/tnm;
   x=a+0.5*del;
   for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
   s=0.5*(s+(b-a)*sum/tnm);
   return s;
 }
}


double get_gauss(double xave, double sig, double min1, double max1)
{ /* returns x drawn from normal/gauss distribution with: (xave,sig) */
 double x,P0,x1,x2,f1,f2;
 double xacc;
 int i;

 do {

   P0=get_t(0.0,1.0);                         /* get a number from a uniform distr: 0.0-1.0 */

   x1=xave-1000.0;
   x2=xave+1000.0;
   f1=int_gauss(x1,xave,sig);
   f2=int_gauss(x2,xave,sig);

   if(P0<f1) P0=int_gauss(x1+1.0,xave,sig);
   if(P0>f2) P0=int_gauss(x1-1.0,xave,sig);

   xacc=1.0e-06;
   x=rtbis(func_x,x1,x2,xacc,xave,sig,P0);
 } while(x<min1 || x>max1);

 return x;
}

double func_x(double x, double xave, double sig, double P0)
{ /* returns P0 - integrated gauss from -inf to x gauss: (xave, sig) */
 double tmp1,tmp2,tmp3;

 tmp1=(x-xave)/(sig*sqrt(2.0));
 if(x<xave)
  tmp2=0.5-0.5*erffNR(-tmp1);     /* integrated gauss */
 else
   tmp2=0.5+0.5*erffNR(tmp1);     /* integrated gauss */

 tmp3=P0-tmp2;

 return tmp3;
}

double int_gauss(double x, double xave, double sig)
{ /* returns integrated gauss from -inf to x gauss: (xave, sig) */
 double tmp1,tmp2;

 tmp1=(x-xave)/(sig*sqrt(2.0));
 if(x<xave)
  tmp2=0.5-0.5*erffNR(-tmp1);
 else
   tmp2=0.5+0.5*erffNR(tmp1);

 return tmp2;
}

double rtbis(double (*func)(double, double, double, double),
             double x1, double x2, double xacc, double xave, double sig, double P0)
{
 int j;
 double dx,f,fmid,xmid,rtb;

 f=(*func)(x1,xave,sig,P0);
 fmid=(*func)(x2,xave,sig,P0);
 if(f*fmid >= 0.0) {
   printf("Root must be bracketed for bisection in rtbis\n");
   exit(-1);
 }
 rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
 for(j=1;j<=JMAX;j++) {
   fmid=(*func)(xmid=rtb+(dx *= 0.5),xave,sig,P0);
   if(fmid <= 0.0) rtb=xmid;
   if(fabs(dx) < xacc || fmid == 0.0) return rtb;
 }

 printf("Too many bisections in rtbis1\n");
 exit(-1);
 return 0.0;
}


void derivs(double x, double y[], double dydx[], double Mca, double lambda)
{ /* return values of derivatives */
 double tmp1,Ma,Mb,a;

 a=y[1];
 Ma=y[2];
 Mb=x;

 tmp1=Cd*(Mb+Ma)-Mb;
 dydx[1]=(Ma*Mb/tmp1+2.0*(2.0*Mb-Mca)/(lambda*Alfa)+Ma)/
         (Ma*Mb*Mb/(tmp1*a)+(Mb/a)*(2.0*(Mb-Mca)/(lambda*Alfa)+Ma));
 dydx[2]=(Ma-Ma*Mb*dydx[1]/a)/tmp1;
}


void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad, void (*derivs)(double, double [], double [], double, double),
	void (*rkqs)(double [], double [], int, double *, double, double, double [], double *, double *,
	void (*)(double, double [], double [], double, double), double, double), double Mca, double lambda)
{
	int nstp,i;
	double xsav,x,hnext,hdid,h;
	double *yscal,*y,*dydx;

	yscal=dvector(1,nvar);
	y=dvector(1,nvar);
	dydx=dvector(1,nvar);
	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = (*nbad) = kount = 0;
	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=1;nstp<=MAXSTP;nstp++) {
		(*derivs)(x,y,dydx,Mca,lambda);
		for (i=1;i<=nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
			xp[++kount]=x;
			for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			xsav=x;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		(*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs,Mca,lambda);
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=1;i<=nvar;i++) ystart[i]=y[i];
			if (kmax) {
				xp[++kount]=x;
				for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			}
			free_dvector(dydx,1,nvar);
			free_dvector(y,1,nvar);
			free_dvector(yscal,1,nvar);
			return;
		}
		if (fabs(hnext) <= hmin)
		  fprintf(fp0,"error: Step size too small in odeint,iidd_old: %d\n",iidd_old);
		h=hnext;
	}
	fprintf(fp0,"error: Too many steps in routine odeint,iidd_old: %d\n",iidd_old);
}


void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double [], double, double), double Mca, double lambda)
{
	int i;
	double errmax,h,xnew,*yerr,*ytemp;

	yerr=dvector(1,n);
	ytemp=dvector(1,n);
	h=htry;
	for (;;) {
		rkck(y,dydx,n,*x,h,ytemp,yerr,derivs,Mca,lambda);
		errmax=0.0;
		for (i=1;i<=n;i++) errmax=DMAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax > 1.0) {
			h=SAFETY*h*pow(errmax,PSHRNK);
			if (h < 0.1*h) h *= 0.1;
			xnew=(*x)+h;
			if (xnew == *x)
			  fprintf(fp0,"error: stepsize underflow in rkqs,iidd_old: %d\n",iidd_old);
			continue;
		} else {
			if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
			else *hnext=5.0*h;
			*x += (*hdid=h);
			for (i=1;i<=n;i++) y[i]=ytemp[i];
			break;
		}
	}
	free_dvector(ytemp,1,n);
	free_dvector(yerr,1,n);
}


void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
	double yerr[], void (*derivs)(double, double [], double [], double, double), double Mca, double lambda)
{
	int i;
	static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
		b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
		dc5 = -277.0/14336.0;
	double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0,dc6=c6-0.25;
	double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
	double tmp1;

	ak2=dvector(1,n);
	ak3=dvector(1,n);
	ak4=dvector(1,n);
	ak5=dvector(1,n);
	ak6=dvector(1,n);
	ytemp=dvector(1,n);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+b21*h*dydx[i];
	(*derivs)(x+a2*h,ytemp,ak2,Mca,lambda);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	(*derivs)(x+a3*h,ytemp,ak3,Mca,lambda);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	(*derivs)(x+a4*h,ytemp,ak4,Mca,lambda);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	(*derivs)(x+a5*h,ytemp,ak5,Mca,lambda);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);

	(*derivs)(x+a6*h,ytemp,ak6,Mca,lambda);
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=1;i<=n;i++) {
		if(dydx[i]<1.0e-100)           /* fortran compilator issue */
		   tmp1=0.0;
		else
		   tmp1=dydx[i];
		yerr[i]=h*(dc1*tmp1+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);

	}
	free_dvector(ytemp,1,n);
	free_dvector(ak6,1,n);
	free_dvector(ak5,1,n);
	free_dvector(ak4,1,n);
	free_dvector(ak3,1,n);
	free_dvector(ak2,1,n);
}


double *dvector(int nl, int nh)
{ /* allocate a double vector with subscript range v[nl..nh] */
 double *v;

 v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
 if (!v)
   fprintf(fp0,"error: allocation failure in dvector()\n");
 return v-nl+NR_END;
}


void free_dvector(double *v, int nl, int nh)
{ /* free a double vector allocated with dvector() */
 free((FREE_ARG) (v+nl-NR_END));
}


double **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}


double get_w1(double M, double R)
{ /* returns angular rot. velocity [Myr^-1] of ZAMS star of mass M [Msun] and radius R [Rsun] */
  /* Andreas Zezas 2002, compilation of current observational data (assumed <sin(i)>=0.707) */
 double vrot;

 if(M<=6.34)
   vrot=(10.0*pow(M,0.0354))/(0.0389+pow(M,-7.95));   /* [km/s] */
 else
   vrot=(13.4*pow(M,-0.12))/(0.0389+pow(M,-7.95));

 return 1.0e+06*45.35*vrot/R;                         /* [Myr^-1] */
}


double get_w2(double M, double R)
{ /* returns angular rot. velocity [Myr^-1] of ZAMS star of mass M [Msun] and radius R [Rsun] */
  /* Hurley et al. 2000, eq.107+108 */
 double vrot;

 vrot=330.0*pow(M,3.3)/(15.0+pow(M,3.45));            /* [km/s] */
 return 1.0e+06*45.35*vrot/R;                         /* [Myr^-1] */
}


int orb_change(double t1, double t2, double *a, double *e, double *wa, double *wb, double tvira, double tvirb,
               double Ma, double Mb, double M0a, double M0b, double Mzamsa, double Mzamsb, double dMwinda,
               double dMwindb, double Mca, double Mcb, double Ra, double Rb, double Raold, double Rbold,
               double La, double Lb, int Ka, int Kb, int Kaold, int Kbold, int mt, int ce, int mttype,
               double *Iaold, double *Ibold, double *KTa, double *KTb, double dMmta, double dMmtb,
               double Mdisa, double Mdisb, int *darwin)
{ /* changes (a,e,wa,wb) due to MT,MB,GR and tidial interactions in binary */
  /* returns 1 if calculation is OK, 2-Darwin instability ODE time step, no MT (violate 1% Jtot max change rule in this case) */
  /* 3-Darwin instability during MT (Exact step), 0 if not sucessful (stop of calculation for this system) */
  /* magbreak1.c: old binary.c with all old magnetic breaking functions (including ODEs,dR,R(t)): Feb 14, 2002 */
 double dt,Rca,Rcb,Ia,Ib,Jawind,Jbwind,wcrit_a,wcrit_b;
 double *ystart,h1,hmin,eps,worb;
 double dJmagb_a,dJmagb_b,Jmb,Jgr,Jorb,J1,J2,c,Maold,Mbold,kw1,kw2,Ac;
 double xacc,x1,x2,del1,del2,atmp1,atmp2;
 double Jmt,Rcom,Rla,Rlb,Ramod,Rbmod;
 double kw3,wsun,wx;
 double Ta,Tb;
 int flag2;
 int mia,mib;
 int magb_a,magb_b,nvar,nok,nbad,mark,i;
 int mark000;
 int ode,dar;

 Ta=get_T(La,Ra);
 Tb=get_T(Lb,Rb);

 if(Ka==-1 || Kb==-1 || Ka==15 || Kb==15) {
   fprintf(fp0,"error: wrong star type (Ka,Kb: %d,%d) in orb_change(), iidd_old: %d\n",Ka,Kb,iidd_old);
   fflush(fp0);
 }

 dt=t2-t1;                                     /* current time step [Myr] */

 if((*e)>0.999) (*e)=0.999;              /* for extremely eccentric systems the integration of odeint1() */
                                         /* goes very slow, so the eccentricity is cut down to e=0.999 */

 Rla=roche(Mb/Ma,(*a)*(1-(*e)));
 Rlb=roche(Ma/Mb,(*a)*(1-(*e)));
 Ramod=Ra;                               /* Hurley equilibrium model star radius: wrong if thermal MT is on */
 Rbmod=Rb;
 if(dMmta>acc && mttype==4)              /* during thermal MT use Roche lobe radius as donor radius!!! */
   Ra=1.01*Rla;
 if(dMmtb>acc && mttype==4)
   Rb=1.01*Rlb;


 /* mass transfer treated first */
 if(dMmta>acc) {         /* Ma donor, both Ma,Mb are already changed due to MT in singl() */
   if((Kb>=10 && Kb<=14) || Kb==16 || Kb==17) {
     worb=sqrt(GGG*(Ma+Mb))*pow(*a,-1.5);            /* mean orbital angular velocity [Myr^-1] */
     Rcom=Ma*(*a)/(Ma+Mb);                           /* orbital angular momentum loss due to MT: */
     Jmt=-(dMmta+dMmtb)*dt*Rcom*Rcom*worb;           /* loss with specific J of compact accretor */
   }
   else {                                            /* regular accretor Kb=0,1,2,3,4,5,6,7,8,9 */
     Jorb=Ma*Mb*sqrt((GGG*(*a))/(Ma+Mb));            /* orbital angular momentum loss due to MT: */
     Jmt=-Beta*((dMmta+dMmtb)*dt/(Ma+Mb))*Jorb;      /* loss with specific ang. momentum of orbital J */
   }
 }
 else if(dMmtb>acc) {         /* Mb donor, both Ma,Mb are already changed due to MT in singl() */
   if((Ka>=10 && Ka<=14) || Ka==16 || Ka==17) {      /* compact object accretor Ka=10,11,12,13,14,16,17 */
     worb=sqrt(GGG*(Ma+Mb))*pow(*a,-1.5);            /* mean orbital angular velocity [Myr^-1] */
     Rcom=Mb*(*a)/(Ma+Mb);                           /* orbital angular momentum loss due to MT: */
     Jmt=-(dMmtb+dMmta)*dt*Rcom*Rcom*worb;           /* loss with specific J of compact accretor */
   }
   else {                                            /* regular accretor Ka=0,1,2,3,4,5,6,7,8,9 */
     Jorb=Ma*Mb*sqrt((GGG*(*a))/(Ma+Mb));            /* orbital angular momentum loss due to MT: */
     Jmt=-Beta*((dMmtb+dMmta)*dt/(Ma+Mb))*Jorb;      /* loss with specific ang. momentum of orbital J */
   }
 }
 else                                                /* no MT associated ang. momentum loss */
   Jmt=0.0;

 mark000=0;
 if(Kaold==-1 || ce==1) {                      /* Kold==-1: don't do anything, let the star become remnant and update its I*/
   Rca=Rcf(tvira,Ma,M0a,Mca,Ka);               /* ce==1: for CE, do not change w of star, stripped of its envelope */
   (*Iaold)=Inerf(Ma,Mca,Ra,Rca,Ka);           /* core radii and inertia calculated for time t2 */
   mark000=1;                                  /* Mca, Mcb correspond to Mche for K=2,3,4,5,6, to */
 }                                             /* Mcco for K=8,9 and are 0 for the rest: good for windf() */
 if(Kbold==-1 || ce==1) {
   Rcb=Rcf(tvirb,Mb,M0b,Mcb,Kb);
   (*Ibold)=Inerf(Mb,Mcb,Rb,Rcb,Kb);
   mark000=1;
 }
 if(mark000==1) {
   (*KTa)=(*KTb)=0.0;
   return 1;
 }

 if(((Ka==0 || Ka==1) && Ma>=minMB && Ma<=maxMB) || Ka==3 || ((Ka==4 || Ka==2) && Ta<5370.0) || Ka==5 || Ka==6) {
   magb_a=1;
   wcrit_a=wcritf(Ma,Ka);
 }
 else {magb_a=0; wcrit_a=0.0;}
 if(((Kb==0 || Kb==1) && Mb>=minMB && Mb<=maxMB) || Kb==3 || ((Kb==4 || Kb==2) && Tb<5370.0) || Kb==5 || Kb==6) {
   magb_b=1;
   wcrit_b=wcritf(Mb,Kb);
 }
 else {magb_b=0; wcrit_b=0.0;}

 Rca=Rcf(tvira,Ma,M0a,Mca,Ka);            /* core radii and inertia calculated for time t2 */
 Rcb=Rcf(tvirb,Mb,M0b,Mcb,Kb);            /* Mca, Mcb correspond to Mche for K=2,3,4,5,6, to */
 Ia=Inerf(Ma,Mca,Ra,Rca,Ka);              /* Mcco for K=8,9 and are 0 for the rest: good for windf() */
 Ib=Inerf(Mb,Mcb,Rb,Rcb,Kb);
 (*KTa)=KTf(tvira,Ma,M0a,Mzamsa,Mca,Ra,Rca,La,*wa,Ka,Mb,*a,*e);
 (*KTb)=KTf(tvirb,Mb,M0b,Mzamsb,Mcb,Rb,Rcb,Lb,*wb,Kb,Ma,*a,*e);


 if((*Iaold)<-1.0) (*Iaold)=Ia;           /* first time step, fill in Iold variables */
 if((*Ibold)<-1.0) (*Ibold)=Ib;

 mia=1;
 if(Ka>=2 && Ka<=9 && Ka!=7 && fabs(Ma-Mca)<0.1*Ma && Raold>Ramod)         /* for small env. mass of giants */
   mia=0;                                                                  /* do not update spin of star */
 if(mt==0 && ce==0 && Ka==7 && (Kaold==2 || Kaold==3 || Kaold==4))         /* giant lost its envelope in wind */
   mia=0;                                                                  /* and becomes Helium MS star */
 if(mt==0 && ce==0 && Ka==8 && Kaold==5)                                   /* giant lost its envelope in wind */
   mia=0;                                                                  /* and become Helium giant star */
                        /* loss of envelope in wind  from K=6,8,9 always give remnant which is excluded below */

 mib=1;
 if(Kb>=2 && Kb<=9 && Kb!=7 && fabs(Mb-Mcb)<0.1*Mb && Rbold>Rbmod)        /* for small env. mass of giants */
   mib=0;                                                                  /* than do not update spin of star */
 if(mt==0 && ce==0 && Kb==7 && (Kbold==2 || Kbold==3 || Kbold==4))         /* giant lost its envelope in wind */
   mib=0;                                                                  /* and become Helium MS star */
 if(mt==0 && ce==0 && Kb==8 && Kbold==5)                                   /* giant lost its envelope in wind */
   mib=0;                                                                  /* and become Helium giant star */
                        /* loss of envelope in wind  from K=6,8,9 always give remnant which is excluded below */


 if((*e)<0.0001) (*e)=0.0;           /* to make sure that the system becomes circular when e drops to almost zero */
 worb=sqrt(GGG*(Ma+Mb))*pow(*a,-1.5);       /* mean orbital angular velocity [Myr^-1] */


 if((*e)<acc && (fabs((*wa)-worb)<0.05*worb || Ka>=10)) {      /* syncA+circ case */
   if(Sia==-1) Sia=1;
   else Sia++;
   if(Ka>=10) Sia=1000;
   if(Ra>0.5*Rla) Sia=1000;                                    /* that may be too strong of an assumption? */
 }
 else
   Sia=-1;

 if((*e)<acc && (fabs((*wb)-worb)<0.05*worb || Kb>=10)) {      /* syncB+circ case */
   if(Sib==-1) Sib=1;
   else Sib++;
   if(Kb>=10) Sib=1000;
   if(Rb>0.5*Rlb) Sib=1000;                                    /* that may be too strong of an assumption? */
 }
 else
   Sib=-1;

 ode=1;
 dar=1;
 if((Sia>10 && Sib>10) || dMmta>acc || dMmtb>acc) {        /* sync+circ case or ongoing RLOF (otherwise use full numerical solution) */
   kw1=8.8801236797060555969e-22;    /* 2.7e+47: constant: [g cm^2 s] -> [Msun Rsun^2 Myr] (New MB) */
   kw2=5.8329147475726926586e-22;    /* 3.8e-30: constant: [s cm^-2] -> [Myr Rsun^-2] (Old MB) */
   kw3=619.2;                        /* Ivanova&Taam: this is Kj=6.e+30 [dyn cm] -> units changed to [Msun Rsun^2 Myr^-2] */
   wsun=9.46e+07;                    /* sun rotational velocity [Myr^-1] */
   wx=9.46e+08;                      /* critcal saturation velocity= 10*w_sun, w_sun=3*10^-6 sec^-1=9.46*10^7 Myr^-1 */
                                     /* I use their eq.(4) with Td/Tdsun=1.0 -- from Ron */

   if(magb_a==1 && (*wa)<=wcrit_a && MB==1)
     dJmagb_a=-kw1*sqrt(Ra/Ma)*(*wa)*(*wa)*(*wa);
   else if(magb_a==1 && (*wa)>wcrit_a && MB==1)
     dJmagb_a=-kw1*sqrt(Ra/Ma)*(*wa)*wcrit_a*wcrit_a;
   else if(magb_a==1 && MB==2)
     dJmagb_a=-kw2*Ma*pow(Ra,gamMB)*(*wa)*(*wa)*(*wa);
   else if(magb_a==1 && (*wa)<=wx && MB==3)
     dJmagb_a=-kw3*pow(Ra,4.0)*pow((*wa)/wsun,3.0);
   else if(magb_a==1 && (*wa)>wx && MB==3)
     dJmagb_a=-kw3*pow(Ra,4.0)*pow((*wa)/wsun,1.3)*pow(wx/wsun,1.7);
   else
     dJmagb_a=0.0;
   if(magb_b==1 && (*wb)<=wcrit_b && MB==1)
     dJmagb_b=-kw1*sqrt(Rb/Mb)*(*wb)*(*wb)*(*wb);
   else if(magb_b==1 && (*wb)>wcrit_b && MB==1)
     dJmagb_b=-kw1*sqrt(Rb/Mb)*(*wb)*wcrit_b*wcrit_b;
   else if(magb_b==1 && MB==2)
     dJmagb_b=-kw2*Mb*pow(Rb,gamMB)*(*wb)*(*wb)*(*wb);
   else if(magb_b==1 && (*wb)<=wx && MB==3)
     dJmagb_b=-kw3*pow(Rb,4.0)*pow((*wb)/wsun,3.0);
   else if(magb_b==1 && (*wb)>wx && MB==3)
     dJmagb_b=-kw3*pow(Rb,4.0)*pow((*wb)/wsun,1.3)*pow(wx/wsun,1.7);
   else
     dJmagb_b=0.0;


   Jmb=dt*(dJmagb_a+dJmagb_b);
   c=13593198857139.902344;                                                /* speed of light [cm/s] -> [Rsun/Myr] */
   Jgr=-(dt*32.0*pow(GGG,3.5)*Ma*Ma*Mb*Mb*sqrt(Ma+Mb))/(5.0*pow(c,5.0)*pow(*a,3.5));
   Maold=Ma+dMwinda*dt+dMmta*dt+Mdisa;                        /* corrected for MT mass loss/gain + WIND mass loss */
   Mbold=Mb+dMwindb*dt+dMmtb*dt+Mdisb;
   Jorb=Maold*Mbold*sqrt((GGG*(*a))/(Maold+Mbold));
   J1=(*Iaold)*(*wa);
   J2=(*Ibold)*(*wb);
   Ac=Jorb+J1+J2+Jmb+Jgr+Jmt;                /* Jtide not accounted for here since here we are in sync+circ case: */
                                             /* and tidial forces do not operate here */

   flag2=0;
   del1=(0.5*(*a))/100.0;
   del2=(*a)/100.0;
   x1=x2=(*a);
   i=0;
   while(flag2==0 && i<100) {                      /* expansive bracketing with initial brackets within: [0.5a:2a] */
     x1-=del1;                                     /* do maksimum 100 loops */
     x2+=del2;
     flag2=zbrac2(func2,&x1,&x2,Ma,Mb,Ia,Ib,Ka,Kb,Ac);
     i++;
   }

   if(flag2==0 && (dMmta>acc || dMmtb>acc)) {      /* MT is going on, Darwin instability encountered */
     ode=0;
     dar=3;
   }
   else if(flag2==0) {                             /* no root Jtot_f-Jtot_i=0: Darwin instability, very probable merger */
     ode=1;                                        /* unless ensuing later MT will stabalize the instability */
     dar=2;
   }                                               /* pass the system to ODEs */
   else {
     ode=0;                                        /* doing exact solution */
     dar=1;
     xacc=1.0e-06*(*a);
     (*a)=rtbis2(func2,x1,x2,xacc,Ma,Mb,Ia,Ib,Ka,Kb,Ac);
     worb=sqrt(GGG*(Ma+Mb))*pow(*a,-1.5);
     if(Ka<10)                               /* keep synchronized only for non-remnants, remnant velocity not changed */
       (*wa)=worb;
     if(Kb<10)
       (*wb)=worb;
   }
 }

 if(ode==1) {
   if((ce==0 && Ka<10 && (*Iaold)>-1.0) || (mt==0 && (*Iaold)>-1.0 && Ka<10 && mia==1)) {             /* don't update for: */
     Jawind=0.66666667*pow((Raold+Ra)/2.0,2.0)*dMwinda*dt*(*wa);             /* first step, after SN,CE, remnats, CO cores */
     (*wa)=((*wa)*(*Iaold)-Jawind)/Ia;
     (*wa)=max(0.0,*wa);                                                          /* star should not spin down below zero! */
   }
   if((ce==0 && Kb<10 && (*Ibold)>-1.0) || (mt==0 && (*Ibold)>-1.0 && Kb<10 && mib==1)) {
     Jbwind=0.66666667*pow((Rbold+Rb)/2.0,2.0)*dMwindb*dt*(*wb);                    /* Iaold, Ibold from previous timestep */
     (*wb)=((*wb)*(*Ibold)-Jbwind)/Ib;                                           /* update during detailed MT calculations */
     (*wb)=max(0.0,*wb);                                                  	        /* see Hurley at al. 2000, eq. 110 */
   }
   nvar=4;
   ystart=dvector(1,nvar);
   eps=1.0e-03;
   h1=0.01*(t2-t1);                /* first time step */
   hmin=0.0;                       /* minimum step */
   kmax=0;                         /* no intermidiate data storage */
   i=0;
   mark=1;

   while(mark==1) {
     ystart[1]=(*a); ystart[2]=(*e); ystart[3]=(*wa); ystart[4]=(*wb);
     mark=odeint1(ystart,nvar,t1,t2,eps,h1,hmin,&nok,&nbad,dorbdt,rkqs1,*KTa,*KTb,Ma,Mb,dMwinda,dMwindb,Ia,Ib,Ra,Rb,Rca,Rcb,
                  wcrit_a,wcrit_b,La,Lb,Ka,Kb,magb_a,magb_b);

     if(dar==2 && i==2) {            /* Darwin Instability System, ODE can not solve it for some reason */
       (*e)=0.0;                     /* bring it to contact, by hand */
       atmp1=0.99*Aroche(Mb/Ma,Ra);
       atmp2=0.99*Aroche(Ma/Mb,Rb);
       (*a)=max(atmp1,atmp2);
       worb=sqrt(GGG*(Ma+Mb))*pow((*a),-1.5);
       if(Ka<10) (*wa)=worb;
       if(Kb<10) (*wb)=worb;
       (*Iaold)=Ia;
       (*Ibold)=Ib;
       (*darwin)+=1;                 /* i==2 --- do not attampt 5 calc. tries, they take a long time */
       return dar;                   /* exit function here, and do MT */
     }

     if(i>=4 || mark==-1 || mark==-2) {
       fprintf(fp0,"error: odeint1 can't solve the case (i: %d, mark: %d), idum_run: %d, iidd_old: %d\n\n",i,mark,idum_run,iidd_old);
       fflush(fp0);
       return 0;
     }
     if(isnan(ystart[1])!=0 || isnan(ystart[2])!=0 || isnan(ystart[3])!=0 || isnan(ystart[4])!=0)
       mark=1;
     if(isinf(ystart[1])!=0 || isinf(ystart[2])!=0 || isinf(ystart[3])!=0 || isinf(ystart[4])!=0)
       mark=1;
     i++;
     h1=h1/pow(10.0,(double)i);
   }

   (*a)=ystart[1];
   (*e)=ystart[2];
   (*wa)=ystart[3];
   (*wb)=ystart[4];
   free_dvector(ystart,1,nvar);
 }

 (*Iaold)=Ia;
 (*Ibold)=Ib;
  return dar;
}


void dorbdt(double t, double y[], double dydx[], double KTa, double KTb,
            double Ma, double Mb, double dMa, double dMb, double Ia, double Ib,
            double Ra, double Rb, double Rca, double Rcb, double wcrit_a, double wcrit_b,
            double La, double Lb, int Ka, int Kb, int magb_a, int magb_b)
{ /* calculates right hand sides of equations for orbit and spin evolution for time t: equations don't depand */
  /*  directly on t, and are calculated for star parameters corresponding to t2, and supplied (a,e,wa,wb) */
  /* (in a given interval: t1--t2): all stellar paramaters are assumed constant within this interval */
  /* only change of stellar spins: (wa,wb) and binary: (a,e)  is calculated */
  /* R,M,dM,Mc,tvir,K,M0,I for each star are constant during t1--t2 */

 double a,e,wa,wb;
 double dJtide_a,dJtide_b,dJmagb_a,dJmagb_b;
 double datide_a,datide_b,detide_a,detide_b,dagrav,degrav;
 double dawind_ab,dewind_ab;
 double rga,rgb,qa,qb,worb,SRa,SRb;
 double f1,f2,f3,f4,f5,e1,e2,kw1,kw2,c;
 double con1,con2,con3,con4,con5,con6,con7,con8,con9,con10,con11;
 double kw3,wsun,wx;
 double Ta,Tb;

 a=y[1];
 e=y[2];
 wa=y[3];
 wb=y[4];

 Ta=get_T(La,Ra);
 Tb=get_T(Lb,Rb);

 if(e<acc) {
   con1=con2=con3=con4=0.0;
   f1=f2=f3=f4=f5=1.0;
 }
 else {
   con1=e*e;             /* e^2 */
   con2=con1*con1;       /* e^4 */
   con3=con1*con2;       /* e^6 */
   con4=con1*con3;       /* e^8 */
   f1=1.0+15.5*con1+31.875*con2+11.5625*con3+0.390625*con4;
   f2=1.0+7.5*con1+5.625*con2+0.3125*con3;
   f3=1.0+3.75*con1+1.875*con2+0.078125*con3;
   f4=1.0+1.5*con1+0.125*con2;
   f5=1.0+3.0*con1+0.375*con2;
 }

 con5=pow(Ra/a,6.0);                      /* (Ra/a)^6 */
 con6=con5*(Ra/a)*(Ra/a);                 /* (Ra/a)^8 */
 con7=pow(Rb/a,6.0);                      /* (Rb/a)^6 */
 con8=con7*(Rb/a)*(Rb/a);                 /* (Rb/a)^8 */

 kw1=8.8801236797060555969e-22;    /* 2.7e+47: constant: [g cm^2 s] -> [Msun Rsun^2 Myr] (New MB) */
 kw2=5.8329147475726926586e-22;    /* 3.8e-30: constant: [s cm^-2] -> [Myr Rsun^-2] (Old MB) */
 c=13593198857139.902344;          /* speed of light [cm/s] -> [Rsun/Myr] */
 kw3=619.2;                        /* Ivanova&Taam: this is Kj=6.e+30 [dyn cm] -> units changed to [Msun Rsun^2 Myr^-2] */
 wsun=9.46e+07;                    /* sun rotational velocity [Myr^-1] */
 wx=9.46e+08;                      /* critcal saturation velocity= 10*w_sun, w_sun=3*10^-6 sec^-1=9.46*10^7 Myr^-1 */
                                               /* I use their eq.(4) with Td/Tdsun=1.0 -- from Ron */

 con9=pow(GGG,3.0)*(Ma+Mb)*(Ma*Mb)/(pow(c,5.0)*a*a*a);    /* for GR decay */
 e1=1.0+3.0416667*con1+0.38541667*con2;
 e2=1.0+0.39802632*con1;


 if(((Ka==0 || Ka==1) && Ma>maxMB) || ((Ka==2 || Ka==4) && Ta>=5370.0) || Ka==7 || Ka==8 || (Ka==9 && Ma>=Mhecon))  /* tides with radiative damping */
   rga=0.31622777;                                                /* rga=k2'=sqrt(0.1) */
 else                                                             /* convective tidial damping */
   rga=sqrt(Ia/Ma)/Ra;                                            /* rga: gyration radius of A */
 if(((Kb==0 || Kb==1) && Mb>maxMB) || ((Kb==2 || Kb==4) && Tb>=5370.0) || Kb==7 || Kb==8 || (Kb==9 && Mb>=Mhecon))
   rgb=0.31622777;
 else
   rgb=sqrt(Ib/Mb)/Rb;                                            /* rgb: gyration radius of B */


 worb=sqrt(GGG*(Ma+Mb))*pow(a,-1.5);       /* mean orbital angular velocity [Myr^-1] */
 SRa=wa/worb;
 SRb=wb/worb;
 qa=Mb/Ma;
 qb=Ma/Mb;

 if(Ka>=0 && Ka<=9)
   dJtide_a=3.0*Ia*KTa*qa*qa*con5*worb*pow(1.0-con1,-6.0)*(f2-pow(1.0-con1,1.5)*f5*SRa)/(rga*rga);
 else dJtide_a=0.0;
 if(Kb>=0 && Kb<=9)
   dJtide_b=3.0*Ib*KTb*qb*qb*con7*worb*pow(1.0-con1,-6.0)*(f2-pow(1.0-con1,1.5)*f5*SRb)/(rgb*rgb);
 else dJtide_b=0.0;

 if(magb_a==1 && wa<=wcrit_a && MB==1)
   dJmagb_a=-kw1*sqrt(Ra/Ma)*wa*wa*wa;
 else if(magb_a==1 && wa>wcrit_a && MB==1)
   dJmagb_a=-kw1*sqrt(Ra/Ma)*wa*wcrit_a*wcrit_a;
 else if(magb_a==1 && MB==2)
   dJmagb_a=-kw2*Ma*pow(Ra,gamMB)*wa*wa*wa;
 else if(magb_a==1 && wa<=wx && MB==3)
   dJmagb_a=-kw3*pow(Ra,4.0)*pow(wa/wsun,3.0);
 else if(magb_a==1 && wa>wx && MB==3)
   dJmagb_a=-kw3*pow(Ra,4.0)*pow(wa/wsun,1.3)*pow(wx/wsun,1.7);
 else
   dJmagb_a=0.0;

 if(magb_b==1 && wb<=wcrit_b && MB==1)
   dJmagb_b=-kw1*sqrt(Rb/Mb)*wb*wb*wb;
 else if(magb_b==1 && wb>wcrit_b && MB==1)
   dJmagb_b=-kw1*sqrt(Rb/Mb)*wb*wcrit_b*wcrit_b;
 else if(magb_b==1 && MB==2)
   dJmagb_b=-kw2*Mb*pow(Rb,gamMB)*wb*wb*wb;
 else if(magb_b==1 && wb<=wx && MB==3)
   dJmagb_b=-kw3*pow(Rb,4.0)*pow(wb/wsun,3.0);
 else if(magb_b==1 && wb>wx && MB==3)
   dJmagb_b=-kw3*pow(Rb,4.0)*pow(wb/wsun,1.3)*pow(wx/wsun,1.7);
 else
   dJmagb_b=0.0;



 if(Ka>=0 && Ka<=9) {
   con10=KTa*qa*(1.0+qa)*con6;
   con11=pow(1.0-con1,1.5)*SRa;
   datide_a=-6.0*con10*a*pow(1.0-con1,-7.5)*(f1-con11*f2);
   detide_a=-27.0*con10*e*pow(1.0-con1,-6.5)*(f3-0.61111111*con11*f4);
 }
 else
   datide_a=detide_a=0.0;

 if(Kb>=0 && Kb<=9) {
   con10=KTb*qb*(1.0+qb)*con8;
   con11=pow(1.0-con1,1.5)*SRb;
   datide_b=-6.0*con10*a*pow(1.0-con1,-7.5)*(f1-con11*f2);
   detide_b=-27.0*con10*e*pow(1.0-con1,-6.5)*(f3-0.61111111*con11*f4);
 }
 else
   datide_b=detide_b=0.0;


 dagrav=-12.8*con9*e1/(pow(1.0-con1,3.5));
 degrav=-20.266667*con9*e*e2/(a*pow(1.0-con1,2.5));

 dawind_ab=a*(dMa+dMb)/(Ma+Mb);   /* assumed circular case: good assumption even for e!=0 */
 dewind_ab=0.0;

 /* 1-a, 2-e, 3-Omspin_a, 4-Omspin_b  */

 dydx[1]=datide_a+datide_b+dagrav+dawind_ab;

 dydx[2]=detide_a+detide_b+degrav+dewind_ab;

 dydx[3]=(dJtide_a+dJmagb_a)/Ia;

 dydx[4]=(dJtide_b+dJmagb_b)/Ib;

}


double E2f(int Ka, double Ma, double Mzamsa)
{ /* calculates tidial coefficient E2 for radiative damping */
  /* if EE2=1 uses  old Zhan/Hurley 2002 formula: eq.43 */
  /* if EE2=2 uses new (my own) fit to Claret 2007 data (A&A 467, 1389): */
  /* by using Mzams in EE2=2 we ignore effects of de/rejuvanation in RLOF on stellar structure */
  /* no actual fit was attempted: the "fit" is by eye only */
 double E2;

 if(Ka>9) {
   E2=0.0;
   fprintf(fp0,"error: func. E2f() should not be here\n");
 }

 if(EE2==1)
   E2=1.592e-09*pow(Ma,2.84);
 else {
   if(Mzamsa<2.3) {
     if(Ka==1) E2=1.0e-08;
     else if(Ka==2) E2=1.0e-10;
     else if(Ka==3 || Ka==4) E2=1.0e-12;
     else if(Ka==5 || Ka==6) E2=1.0e-14;
     else if(Ka>=7 && Ka<=9) E2=1.0e-15;   /* if helium star is encountered here it must be massive */
   }                                       /* with radiative envelope: as set in KTf() */
   else if(Mzamsa>=2.3 && Mzamsa<7.1) {
     if(Ka==1) E2=1.0e-07;
     else if(Ka==2) E2=1.0e-09;
     else if(Ka==3 || Ka==4) E2=1.0e-11;
     else if(Ka==5 || Ka==6) E2=1.0e-13;
     else if(Ka>=7 && Ka<=9) E2=1.0e-15;
   }
   else if(Mzamsa>=7.1 && Mzamsa<35.5) {
     if(Ka==1) E2=1.0e-06;
     else if(Ka==2) E2=1.0e-09;
     else if(Ka==3 || Ka==4) E2=1.0e-11;
     else if(Ka==5 || Ka==6) E2=1.0e-13;
     else if(Ka>=7 && Ka<=9) E2=1.0e-15;
   }
   else if(Mzamsa>=35.5) {
     if(Ka==1) E2=1.0e-06;
     else if(Ka==2) E2=1.0e-09;
     else if(Ka>=3 && Ka<=6) E2=1.0e-12;    /* usually at these masses only K=4 is encountered */
     else if(Ka>=7 && Ka<=9) E2=1.0e-15;
   }
 }

 return E2;
}


double KTf(double tvira, double Ma, double M0a, double Mzamsa, double Mca, double Ra,
           double Rca, double La, double wa, int Ka, double Mb, double a, double e)
{ /* calculates coefficent (k/T)_c [Myr^-1] in tidal equations (25,26) as given in Hurley et al. 2002 */
  /* calculates (k/T) both for convective tides as for radiative damping tides */
  /* tvira is needed only for MS and HG stars */
 double kta,fconv,tconv,Menv,Renv;
 double Porb,Pspin,Ptid;
 double qa,E2,cal,Ta;

 Ta=get_T(La,Ra);

 if(Tides==0)                                                                /* calculation with no tides */
   kta=0.0;
 else if(Ka<0 || Ka>=10)             /* remnants and CO cores: no tidal or radiative damping interactions */
   kta=0.0;
 else if(((Ka==0 || Ka==1) && Ma>maxMB) || ((Ka==2 || Ka==4) && Ta>=5370.0) || Ka==7 || Ka==8 || (Ka==9 && Ma>=Mhecon)) {
   qa=Mb/Ma;                                                        /* sync/circ with radiative dampining */
   cal=Calrad;
   E2=E2f(Ka,Ma,Mzamsa);
   kta=cal*1.0e+06*19782.0*sqrt(Ma*Ra*Ra/(a*a*a*a*a))*pow(1.0+qa,5.0/6.0)*E2;      /* [yr^-1] -> [Myr^-1] */
 }
 else {                                                  /*  sync/circ with tides in convective envelopes */
   Renv=Renv_con(tvira,Ma,M0a,Ra,Rca,La,Ka);
   Menv=Menv_con(tvira,Ma,M0a,Mca,Ra,La,Ka);
   if(Renv<0.0001*Ra || Menv<0.000001*Ma)                                /* almost no convective envelope */
     if(Ka==0 || Ka==1 || Ka==7) {
       qa=Mb/Ma;
       cal=Calrad;
       E2=E2f(Ka,Ma,Mzamsa);
       kta=cal*1.0e+06*19782.0*sqrt(Ma*Ra*Ra/(a*a*a*a*a))*pow(1.0+qa,5.0/6.0)*E2;  /* [yr^-1] -> [Myr^-1] */
     }
     else
       kta=0.0;
   else {
     tconv=0.4311*pow(Menv*Renv*(Ra-0.5*Renv)/(3.0*La),1.0/3.0);      /* [yr], input M,R,L in solar units */
     Porb=2.0*Pi*a*sqrt(a/(G*(Ma+Mb)))/365.25;                                           /* [day] -> [yr] */
     if(wa<acc)                                                     /* to protect division by 0 if wa=0.0 */
       Ptid=Porb;
     else {
       Pspin=1.0e+06*2.0*Pi/wa;                                                          /* [Myr] -> [yr] */
       if(fabs(Porb-Pspin)<1.0e-12)                                          /* to avoid division by zero */
         Ptid=1.0e+12;
       else
         Ptid=1.0/fabs(1.0/Porb-1.0/Pspin);                                                       /* [yr] */
     }
     fconv=min(1.0,pow(Ptid/(2.0*tconv),2.0));                                              /* [unitless] */
     cal=Caltid;
     kta=cal*1.0e+06*0.095238095*fconv*Menv/(tconv*Ma);                            /* [yr^-1] -> [Myr^-1] */
   }
 }

 return kta;
}


double tau_tidef(double Ma, double Mb, double a, double e, double wa, double Ia, double KTa, double Ra, double La, int Ka)
{ /* calculates timescale for tidal interaction of star A with orbit */
  /* timescale: positive if orbit shrinks, negative if orbit expands: for direct use in dMmtf() */
  /* WATCH OUT: make fabs() if you want real timecale */
 double dJtide_a,Jorb,rga,qa,worb,SRa,tau_tide;
 double f2,f5,con1,con2,con3,con5,Ta;

 Ta=get_T(La,Ra);

 if(e<acc) {
   con1=con2=con3=0.0;
   f2=f5=1.0;
 }
 else {
   con1=e*e;             /* e^2 */
   con2=con1*con1;       /* e^4 */
   con3=con1*con2;       /* e^6 */
   f2=1.0+7.5*con1+5.625*con2+0.3125*con3;
   f5=1.0+3.0*con1+0.375*con2;
 }
 con5=pow(Ra/a,6.0);                      /* (Ra/a)^6 */

 if(((Ka==0 || Ka==1) && Ma>maxMB) || ((Ka==2 || Ka==4) && Ta>=5370.0) || Ka==7  || Ka==8 || (Ka==9 && Ma>=Mhecon))  /* tides with radiative damping */
   rga=0.31622777;                                                /* rga=k2'=sqrt(0.1) */
 else                                                             /* convective tidial damping */
   rga=sqrt(Ia/Ma)/Ra;                                            /* rga: gyration radius of A */

 Jorb=Ma*Mb*sqrt((GGG*a)/(Ma+Mb))*sqrt(1.0-e*e);                  /* orbital angular momentum */
 worb=sqrt(GGG*(Ma+Mb))*pow(a,-1.5);                              /* mean orbital angular velocity [Myr^-1] */
 SRa=wa/worb;
 qa=Mb/Ma;

 if(Ka>=0 && Ka<=9)
   dJtide_a=3.0*Ia*KTa*qa*qa*con5*worb*pow(1.0-con1,-6.0)*(f2-pow(1.0-con1,1.5)*f5*SRa)/(rga*rga);
 else
   dJtide_a=0.0;

 if(fabs(dJtide_a)<acc)
   tau_tide=10000000000000000000000000.0;                         /* no tide */
 else
   tau_tide=(-Jorb/dJtide_a);

 return tau_tide;                                                 /* tidal interaction timescale [Myr] */
}


double Inerf(double Ma, double Mca, double Ra, double Rca, int Ka)
{ /* calculates moment of inertia I [Msun Rsun^2] for star A, see eq. 109 of Hurley et al. 2000 */
  /* also Lai et al. 1993, ApJS 88, 205, table 1 */
  /* core and envelope treated speparately I=Ienv+Icore */
 double k2,k3,k4;

 if(Ma<=1.0) k3=0.204596;                     /* polytrop with n=1.5 (low mass MS stars), Lai */
 else k3=0.075356;                            /* polytrop with n=3.0 (high mass MS stars), Lai */
 k2=0.1;                                      /* eq. 109 of Hurley et al. 2000*/
 k4=0.204596;                                 /* k3,k2=2/5*kn: as needed for I=2/5*kn*M*R^2 */

 if(Ka==0 || Ka==1 || Ka==7 || Ka==-1 || Ka>=10)   /* MS,HeMS,remnants, CO core: do not have cores */
   return k3*Ma*Ra*Ra;
 else
   return k2*(Ma-Mca)*Ra*Ra+k4*Mca*Rca*Rca;        /* giants with cores and envelopes */
}


double Renv_con(double tvira, double Ma, double M0a, double Ra, double Rca, double La, int Ka)
{ /* returns radial extent (depth) of convective envelope [Msun] for a star of Ma,Mca,Ka */
  /* tvira is needed only for MS and HG stars */
 double Rconv0,Rconv,R1,tms,tms1,tau,tt,Ta;

 Ta=get_T(La,Ra);   /* effective temperature of A [K] */

 if(Ka==-1 || Ka>=10 || Ka==7 || Ka==8 || (Ka==9 && Ma>=Mhecon) || ((Ka==2 || Ka==4) && Ta>=5370.0) || ((Ka==0 || Ka==1) && Ma>maxMB))
   Rconv=0.0;                                                                   /* no conv. envelope */
 else if(Ka==3 || Ka==4 || Ka==5 || Ka==6 || Ka==9)
   Rconv=Ra-Rca;

 else if((Ka==0 || Ka==1) && Ma<minMB) {
   tms=tmsf(Ma);
   Rconv0=Ra;
   tau=tvira/tms;                                   /* tvira gives proper evolutionary time for star */
   Rconv=Rconv0*pow(1.0-tau,0.25);
 }
 else if((Ka==0 || Ka==1) && Ma<=maxMB) {
   tms=tmsf(Ma);
   tms1=tmsf(minMB);
   tau=tvira/tms;
   tt=tau*tms1;
   R1=Rmsf(minMB,tt);
   Rconv0=R1*pow((maxMB-Ma)/0.9,0.5);
   Rconv=Rconv0*pow(1.0-tau,0.25);
 }
 else if(Ka==2) {
   tms=tmsf(M0a);
   tau=(tvira-tms)/(tbgbf(M0a)-tms);
   Rconv=pow(tau,0.5)*(Ra-Rca);
 }
 else {
   fprintf(fp0,"error: in Renv_con() unknown Ka type: %d, iidd_old: %d",Ka,iidd_old);
   fflush(fp0);
 }

 return Rconv;
}


double Menv_con(double tvira, double Ma, double M0a, double Mca, double Ra, double La, int Ka)
{ /* returns mass of convective envelope [Msun] for a star of Ma,Mca,Ka */
  /* tvira is needed only for MS and HG stars */
 double Mconv0,Mconv,tms,tau,Ta;

 Ta=get_T(La,Ra);   /* effective temperature of A [K] */

 if(Ka==-1 || Ka>=10 || Ka==7 || Ka==8 || (Ka==9 && Ma>=Mhecon) || ((Ka==2 || Ka==4) && Ta>=5370.0) || ((Ka==0 || Ka==1) && Ma>maxMB))
   Mconv=0.0;                                                                        /* no conv. envelope */
 else if(Ka==3 || Ka==4 || Ka==5 || Ka==6 || Ka==9)
   Mconv=Ma-Mca;
 else if((Ka==0 || Ka==1) && Ma<minMB) {
   Mconv0=Ma;
   tau=tvira/tmsf(Ma);                                   /* tvira gives proper evolutionary time for star */
   Mconv=Mconv0*pow(1.0-tau,0.25);
 }
 else if((Ka==0 || Ka==1) && Ma<=maxMB) {
   Mconv0=minMB*pow((maxMB-Ma)/0.9,2.0);
   tau=tvira/tmsf(Ma);                                   /* tvira gives proper evolutionary time for star */
   Mconv=Mconv0*pow(1.0-tau,0.25);
 }
 else if(Ka==2) {
   tms=tmsf(M0a);
   tau=(tvira-tms)/(tbgbf(M0a)-tms);
   Mconv=tau*(Ma-Mca);
 }
 else {
   fprintf(fp0,"error: in Menv_con() unknown Ka type: %d, iidd_old: %d",Ka,iidd_old);
   fflush(fp0);
 }

 return Mconv;
}


double Rcf(double tvira, double Ma, double M0a, double Mca, int Ka)
{ /* core radius estimate for tvira -- corresponding to the end of a given general binary timestep */
  /* Feb 13, 2002: calculates radius of the core (similar to perturb func.) */
  /* uses hrdiag.f l.626,628 which defines what shall happen for K=5 if remnant is K=8 */
  /* m0 should be filled with M0a in binary.c */
 double tbagb,the1,the,tau,L,Lc,Lths,Rc,Rzhs,Mche,Mcco,beta;
 int Kdummy;

 if(Ka==0 || Ka==1 || Ka==7 || Ka==-1 || Ka>=10)   /* MS,HeMS,remnants, CO core: do not have cores */
   Rc=0.0;
 else if(Ka==2 || Ka==3) {
   if(M0a>M_HeF)
     Rc=Rzhsf(Mca);
   else
     Rc=Rwdf(Mca,10);
 }
 else if(Ka==4) {
   the1=the1f(M0a);
   the=thef(M0a);
   tau=(tvira-the1)/the;
   beta=max(0.0,0.4-0.22*log10(Mca));
   Rzhs=Rzhsf(Mca);
   Rc=Rzhs*(1.0+beta*tau-beta*pow(tau,6.0));
 }
 else if(Ka==5) {
   tbagb=tbagbf(M0a);
   tau=3.0*(tvira-tbagb)/(tnfI(M0a,M0a,Ka)-tbagb);
   Mche=Mcbagbf(M0a);
   Mcco=Mceagbf1(M0a,tvira);
   L=Leagbf2(M0a,Mcco);
   Lc=Lhsgbf(Mche,Mcco);                       /* from paper0, luminosity of Helium Giant */
   if(tau<1.0) {                               /* from hrdiag.f, looks like Luminosity of Helium Hartprung Gap */
     Lths=Lthsf(Mche);
     Lc=Lths*pow(Lc/Lths,tau);
   }
   Lths=Lthsf(Mche);                           /* Mche is M0 for the hipotetical Helium Star here */
   Rzhs=Rzhsf(Mche);
   Rc=Rhsgbf(Mche,Lc,Lths,Rzhs,&Kdummy);       /* K dummy variable here,changed to 8 or 9, */
 }                                             /* but only in the scope of this func */
 else {                                        /* they assume that for K=6,8,9, remnant will be CO WD: K=11 */
   if(Mca>=MCh)                                /* WATCH!!!: perturbation as applied is only my approximation! */
     Rc=Rwdf(MCh-acc,11);                      /* Jarrod do not define pert. for star of K=6,8,9 if its core */
   else                                        /* mass is greater then MCh -- so I leave perturbation on level */
     Rc=Rwdf(Mca,11);                          /* of MCh-acc and as for WD, although it shall be for massive */
 }                                             /* C0 star -- models of Pools to come */

 return Rc;
}

double wcritf(double Ma, int Ka)
{ /* returns critical mag break. staruration angular rot. velocity of star of mass Ma [Msun] in [Myr^-1] */
 double wcrit;
 double Ds=86400.0;		/* day in [s] */


 if(Ka!=0 && Ka!=1)                     /* for star not on MS I do not have values of wcrit!!! */
   wcrit=1000000.0;                     /* huge value assumed, not to allow star enter mb saturation phase */
 else if(Ma>=0.1 && Ma<=0.5)
   wcrit=-0.0000005640+0.0000144400*Ma;
 else if(Ma>0.5 && Ma<=0.6)
   wcrit=-0.0000174400+0.0000484000*Ma;
 else if(Ma>0.6 && Ma<=0.7)
   wcrit=-0.0000058000+0.0000290000*Ma;
 else if(Ma>0.7 && Ma<=0.8)
   wcrit=-0.0000114000+0.0000370000*Ma;
 else if(Ma>0.8 && Ma<=0.9)
   wcrit=-0.0000178000+0.0000450000*Ma;
 else if(Ma>0.9 && Ma<=1.0)
   wcrit=-0.0000430000+0.0000730000*Ma;
 else if(Ma>1.0 && Ma<=1.1)
   wcrit=-0.0001090000+0.0001390000*Ma;
 else if(Ma>1.1 && Ma<=1.2)
   wcrit=-0.0004258000+0.0004270000*Ma;
 else if(Ma>1.2 && Ma<=2.0)
   wcrit=1000000.0;                 /* huge value assumed, not to allow star enter mb saturation phase */
 else {
   fprintf(fp0,"error: donor mass Ma out of range in mag_break0()\n");
   fflush(fp0);
 }
 wcrit=wcrit*Ds*1000000.0*365.25;   /* [s^-1] -> [Myr^-1] */

 return wcrit;
}


int odeint1(double ystart[], int nvar, double x1, double x2, double eps, double h1, double hmin,
             int *nok, int *nbad,
 	     void (*derivs1)(double, double [], double [], double, double, double, double, double,
 	                     double, double, double, double, double, double, double, double, double,
 	                     double, double, int, int, int, int),
 	     int (*rkqs1)(double [], double [], int, double *, double, double, double [], double *, double *,
 	                  void (*)(double, double [], double [], double, double, double, double, double,
 	                           double, double, double, double, double, double, double, double, double,
 	                           double, double, int, int, int, int),
 	                  double, double, double, double, double, double, double, double, double, double,
 	                  double, double, double, double, double, double, int, int, int, int),
             double KTa, double KTb, double Ma, double Mb, double dMa, double dMb, double Ia, double Ib,
             double Ra, double Rb, double Rca, double Rcb, double wcrit_a, double wcrit_b,
             double La, double Lb, int Ka, int Kb, int magb_a, int magb_b)
{ /* NR ver. 2, the odeint() may be something else */

	int nstp,i,mark;
	double xsav,x,hnext,hdid,h;
	double *yscal,*y,*dydx;
        double Rla,Rlb,dper;

	yscal=dvector(1,nvar);
	y=dvector(1,nvar);
	dydx=dvector(1,nvar);
	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = (*nbad) = kount = 0;
	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=1;nstp<=MAXSTP1;nstp++) {

	        dper=y[1]*(1.0-y[2]);
	        Rla=roche(Mb/Ma,dper);
	        Rlb=roche(Ma/Mb,dper);
	        if(Ra>=1.01*Rla || Rb>1.01*Rlb) {      /* my stop: if A or B overfill roche lobe by more than 1% */
	      	  for(i=1;i<=nvar;i++) ystart[i]=y[i];
		  if(kmax) {
		    xp[++kount]=x;
		    for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
		  }
		  free_dvector(dydx,1,nvar);
		  free_dvector(y,1,nvar);
		  free_dvector(yscal,1,nvar);
	          return 0;
		}

	 	(*derivs1)(x,y,dydx,KTa,KTb,Ma,Mb,dMa,dMb,Ia,Ib,Ra,Rb,Rca,Rcb,wcrit_a,wcrit_b,
	 	           La,Lb,Ka,Kb,magb_a,magb_b);
		for (i=1;i<=nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
			xp[++kount]=x;
			for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			xsav=x;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		mark=(*rkqs1)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs1,KTa,KTb,Ma,Mb,dMa,dMb,Ia,Ib,Ra,Rb,Rca,Rcb,
		              wcrit_a,wcrit_b,La,Lb,Ka,Kb,magb_a,magb_b);
	        if(mark==-2) return -2;
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=1;i<=nvar;i++) ystart[i]=y[i];
			if (kmax) {
				xp[++kount]=x;
				for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			}
			free_dvector(dydx,1,nvar);
			free_dvector(y,1,nvar);
			free_dvector(yscal,1,nvar);
			return 0;
		}
		if (fabs(hnext) <= hmin) /* error: step size too small in routine odeint1 */
	          return -1;
		h=hnext;
	}
        /* warning: too many steps in routine odeint1 */
	return 1;
}

int rkqs1(double y[], double dydx[], int n, double *x, double htry, double eps,
           double yscal[], double *hdid, double *hnext,
           void (*derivs1)(double, double [], double [], double, double, double, double, double,
 	                   double, double, double, double, double, double, double, double, double,
 	                   double, double, int, int, int, int),
           double KTa, double KTb, double Ma, double Mb, double dMa, double dMb, double Ia, double Ib,
           double Ra, double Rb, double Rca, double Rcb, double wcrit_a, double wcrit_b,
           double La, double Lb, int Ka, int Kb, int magb_a, int magb_b)
{
	int i;
	double errmax,h,htemp,xnew,*yerr,*ytemp;

	yerr=dvector(1,n);
	ytemp=dvector(1,n);
	h=htry;
	for (;;) {
		rkck1(y,dydx,n,*x,h,ytemp,yerr,derivs1,KTa,KTb,Ma,Mb,dMa,dMb,Ia,Ib,Ra,Rb,Rca,Rcb,
		      wcrit_a,wcrit_b,La,Lb,Ka,Kb,magb_a,magb_b);
		if(ytemp[1]<0.0 || ytemp[2]<0.0 || isnan(ytemp[1])!=0 || isnan(ytemp[2])!=0 || isinf(ytemp[1])!=0 || isinf(ytemp[2])!=0) {
		  h/=2.0;
		  continue;
	        }
	        if(y[1]<0.0 || y[2]<0.0) {
	          h/=2.0;
	          continue;
	        }
		errmax=0.0;
		for (i=1;i<=n;i++) errmax=DMAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=(h >= 0.0 ? DMAX(htemp,0.1*h) : DMIN(htemp,0.1*h));
		xnew=(*x)+h;
		if (xnew == *x)  /* error: stepsize underflow in rkqs1 */
		  return -2;
	}
	if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
	else *hnext=5.0*h;
	*x += (*hdid=h);
	for (i=1;i<=n;i++) y[i]=ytemp[i];
	free_dvector(ytemp,1,n);
	free_dvector(yerr,1,n);
        return 0;
}

void rkck1(double y[], double dydx[], int n, double x, double h, double yout[], double yerr[],
           void (*derivs1)(double, double [], double [], double, double, double, double, double,
 	                   double, double, double, double, double, double, double, double, double,
 	                   double, double, int, int, int, int),
           double KTa, double KTb, double Ma, double Mb, double dMa, double dMb, double Ia, double Ib,
           double Ra, double Rb, double Rca, double Rcb, double wcrit_a, double wcrit_b,
           double La, double Lb, int Ka, int Kb, int magb_a, int magb_b)
{
	int i;
	static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
		b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
		dc5 = -277.00/14336.0;
	double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0,dc6=c6-0.25;
	double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;

	ak2=dvector(1,n);
	ak3=dvector(1,n);
	ak4=dvector(1,n);
	ak5=dvector(1,n);
	ak6=dvector(1,n);
	ytemp=dvector(1,n);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+b21*h*dydx[i];
	(*derivs1)(x+a2*h,ytemp,ak2,KTa,KTb,Ma,Mb,dMa,dMb,Ia,Ib,Ra,Rb,Rca,Rcb,wcrit_a,wcrit_b,La,Lb,Ka,Kb,magb_a,magb_b);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	(*derivs1)(x+a3*h,ytemp,ak3,KTa,KTb,Ma,Mb,dMa,dMb,Ia,Ib,Ra,Rb,Rca,Rcb,wcrit_a,wcrit_b,La,Lb,Ka,Kb,magb_a,magb_b);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	(*derivs1)(x+a4*h,ytemp,ak4,KTa,KTb,Ma,Mb,dMa,dMb,Ia,Ib,Ra,Rb,Rca,Rcb,wcrit_a,wcrit_b,La,Lb,Ka,Kb,magb_a,magb_b);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	(*derivs1)(x+a5*h,ytemp,ak5,KTa,KTb,Ma,Mb,dMa,dMb,Ia,Ib,Ra,Rb,Rca,Rcb,wcrit_a,wcrit_b,La,Lb,Ka,Kb,magb_a,magb_b);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	(*derivs1)(x+a6*h,ytemp,ak6,KTa,KTb,Ma,Mb,dMa,dMb,Ia,Ib,Ra,Rb,Rca,Rcb,wcrit_a,wcrit_b,La,Lb,Ka,Kb,magb_a,magb_b);
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=1;i<=n;i++)
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
	free_dvector(ytemp,1,n);
	free_dvector(ak6,1,n);
	free_dvector(ak5,1,n);
	free_dvector(ak4,1,n);
	free_dvector(ak3,1,n);
	free_dvector(ak2,1,n);
}


double dMamtf(double Ma, double Mb, double La, double Ra, double Rla, double a, int Ka)
{       /* atmospheric Roche Lobe Overflow: return s dM/dt [POSITIVE SIGN] */
	/* Returns mass transfer rate from donor to acceptor (roles must be previously derived from Roche lobe filling factor) [Msun/Myr]*/
	/* Input: Ma-donor Mbss [Msun]; Mb-acceptor Mbss [Msun]; La-donor effective luminosity [K]; */
	/* Ra-donor radius [Rsun]; Rla - donor Roche lobe radius [Rsun]; a-orbital separation [R_sun]; Ka-type of donor ; z-metallicity*/
	/* Should be accepted with coution for binaries with 10.<q<0.5 */
	/* Based on Ritter H., 1988, ApJ, 202, 93 */

	double f2,g,F,gamma,u;

	double Ta=get_T(La,Ra);

	if (Ka<=6){
		u=0.612; /*mean molecular weight*/
	}
	else if (Ka>=7 && Ka<=9){
		u=1.34; /*mean molecular weight*/
	}
	else fprintf(fp0,"error: wrong donor type in function dMamtf %d",Ka);


	/* Calculating the F function from Ritter (1988) Appendix */
	if(0.5<(Ma/Mb) && (Ma/Mb)<10.0){
		F=1.23+0.5*log10(Ma/Mb); /* aproxiMbtion - Eq. A9 in Ritter (1988) */
	}
	else {
		f2=Rla/a; /* f2 as defined in eq. (A6) in Ritter 1988 */
		g=(Ma/Mb)/pow(Rla,3.)+1./pow(1.-Rla,3.); /* Eq. A5 in Ritter (1988); X_L1 aproxiMbted with R_L1 */
		F=pow((g-(1.+Ma/Mb))*g,-0.5)*pow(f2,-1.5); /* Eq. A8 in Ritter (1988); X_L1 aproxiMbted with R_L1 */

	}


	double Rgas=8.314462175*1.0e+07;  /* Gas constant	[erg K^-1 mol^-1]*/
	double Gc=6.6725985*1.0e-08;   /* Gravitational constant	 [cm^3 g^-1 s^-2]*/

	double M0=(2.*M_PI/sqrt(M_E))*pow(Rgas*Ta/(u),1.5)*pow(Ra*Rsun,3.)*ro_ph(Ma,Ra,Ta)*F/(Gc*(Ma*Msun)); /* Eq. A11 in Ritter (1988) */

	double Hp=Rgas*Ta*pow(Ra*Rsun,2)/((u)*Gc*Ma*Msun);  /* Eq. 8b in Ritter (1988) */


	if (0.04<(Ma/Mb) && (Ma/Mb)<1.0){
		 gamma=0.954+0.025*log10(Ma/Mb)-0.038*log10(Ma/Mb)*log10(Ma/Mb); /* Eq. 7 in Ritter (1988) */
	}
	else if (1.0<(Ma/Mb) && (Ma/Mb)<20.0){
		 gamma=0.954+0.039*log10(Ma/Mb)+0.114*log10(Ma/Mb)*log10(Ma/Mb); /* Eq. 7 in Ritter (1988) */
	}
	else{
		 gamma=1.; /* since gamMb is always close to unity */
	}

	Hp/=gamma; /* Eq. 8c in Ritter (1988) */
	return M0*exp(-(Rla-Ra)*Rsun/Hp)*(31556926.*1.0e+06/Msun); /* Eq. 11 in Ritter (1988) scaled to Msun/Myr */
}


double ro_ph(double M, double R, double Teff)
{ /* Calculates density in the photosphere [g/cm^2]*/
  /* If avalaible it takes a interpolated density at kappa_ross=2/3 from MARCS model atmospheres */
  /* Input:  M-Mass of the star in sun units, Mcore-Mbss of the star core, R-radius of the star in sun units, Rcore-radius of the star core, Teff-effective temperature*/

	double logT=log10(Teff);
	double logg= log10(27542.29*M/(R*R));   /* log(g) of the donor - scaled with Sun's surface gravity*/
	double ro;

	if (ZZ>0.048359){
		ro=ro_ph_z048359(logg,logT)+(ro_ph_z067218(logg,logT)-ro_ph_z048359(logg,logT))*(ZZ-0.048359)/0.018859; /*interpolating density for the given metallicity*/
	}
	else if (ZZ>0.0322626){
		ro=ro_ph_z0322626(logg,logT)+(ro_ph_z048359(logg,logT)-ro_ph_z0322626(logg,logT))*(ZZ-0.0322626)/0.0160964; /*interpolating density for the given metallicity*/
	}
	else if (ZZ>0.020267){
		ro=ro_ph_z020267(logg,logT)+(ro_ph_z0322626(logg,logT)-ro_ph_z020267(logg,logT))*(ZZ-0.020267)/0.0119956; /*interpolating density for the given metallicity*/
	}
	else if (ZZ>0.0122){
		ro=ro_ph_z0122(logg,logT)+(ro_ph_z020267(logg,logT)-ro_ph_z0122(logg,logT))*(ZZ-0.0122)/0.008067; /*interpolating density for the given metallicity*/
	}
	else if (ZZ>0.0070598){
		ro=ro_ph_z0070598(logg,logT)+(ro_ph_z0122(logg,logT)-ro_ph_z0070598(logg,logT))*(ZZ-0.0070598)/0.0051402; /*interpolating density for the given metallicity*/
	}
	else if (ZZ>0.004050){
		ro=ro_ph_z004050(logg,logT)+(ro_ph_z0070598(logg,logT)-ro_ph_z004050(logg,logT))*(ZZ-0.004050)/0.0030098; /*interpolating density for the given metallicity*/
	}
	else if (ZZ>0.0023094){
		ro=ro_ph_z0023094(logg,logT)+(ro_ph_z004050(logg,logT)-ro_ph_z0023094(logg,logT))*(ZZ-0.0023094)/0.0017406; /*interpolating density for the given metallicity*/
	}
	else if (ZZ>0.0013113){
		ro=ro_ph_z0013113(logg,logT)+(ro_ph_z0023094(logg,logT)-ro_ph_z0013113(logg,logT))*(ZZ-0.0013113)/0.0009981; /*interpolating density for the given metallicity*/
	}
	else if (ZZ>0.000421151){
		ro=ro_ph_z000421151(logg,logT)+(ro_ph_z0013113(logg,logT)-ro_ph_z000421151(logg,logT))*(ZZ-0.000421151)/0.000890149; /*interpolating density for the given metallicity*/
	}
	else if (ZZ>0.0001338401){
		ro=ro_ph_z0001338401(logg,logT)+(ro_ph_z000421151(logg,logT)-ro_ph_z0001338401(logg,logT))*(ZZ-0.0001338401)/0.000287311; /*interpolating density for the given metallicity*/
	}
	else if (ZZ>0.0000423904){
		ro=ro_ph_z0000423904(logg,logT)+(ro_ph_z0001338401(logg,logT)-ro_ph_z0000423904(logg,logT))*(ZZ-0.0000423904)/0.00009145; /*interpolating density for the given metallicity*/
	}
	else if (ZZ>0.0000134117){
		ro=ro_ph_z0000134117(logg,logT)+(ro_ph_z0000423904(logg,logT)-ro_ph_z0000134117(logg,logT))*(ZZ-0.0000134117)/0.000028979; /*interpolating density for the given metallicity*/
	}
	else if (ZZ>0.000001341446){
		ro=ro_ph_z000001341446(logg,logT)+(ro_ph_z0000134117(logg,logT)-ro_ph_z000001341446(logg,logT))*(ZZ-0.000001341446)/0.00001207; /*interpolating density for the given metallicity*/
	}
	else {
		ro=ro_ph_z00000013415(logg,logT)+(ro_ph_z000001341446(logg,logT)-ro_ph_z00000013415(logg,logT))*(ZZ-0.00000013415)/0.000001207; /*interpolating density for the given metallicity*/
	}


	if (ro<0.){ /* if extrapolation was too far*/
		ro=0.;
	}

	return ro;
}


double ro_ph_z00000013415(double logg, double logT)
{ /* Calculate density in the photosphere for the z=0.00000013415 [M/H]=-5*/
        /*density at kappa_ross=2/3 from MbRCS model atmospheres*/
		double log_ro1,log_ro2,log_ro;

		if (logg>=4.5 && 3.49<logT && logT <3.89){
				log_ro1=2.26496527471945e+08+logT*(-4.00169957517793e+08)+pow(logT,2)*(2.86967929977283e+08)+pow(logT,3)*(-1.00248810012839e+08)+pow(logT,4)*(1.24947497383337e+07)+pow(logT,5)*(2.97751676794779e+06)+pow(logT,6)*(-1.45346615288982e+06)+pow(logT,7)*(2.47958474490435e+05)+pow(logT,8)*(-2.03246239591655e+04)+pow(logT,9)*(6.63098991105944e+02);
				log_ro2=6.69883963619169e+07+logT*(-9.61972663361916e+07)+pow(logT,2)*(4.93928274250201e+07)+pow(logT,3)*(-7.79836645166545e+06)+pow(logT,4)*(-1.81566525133977e+06)+pow(logT,5)*(6.26610766369138e+05)+pow(logT,6)*(7.93050903546093e+04)+pow(logT,7)*(-5.63568349732258e+04)+pow(logT,8)*(8.66160274805517e+03)+pow(logT,9)*(-4.55315819993596e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 && 3.477<logT && logT <3.89){
				log_ro1=1.47056003132479e+09+logT*(-2.79485951795742e+09)+pow(logT,2)*(2.18063596991136e+09)+pow(logT,3)*(-8.48504842120225e+08)+pow(logT,4)*(1.32156159430640e+08)+pow(logT,5)*(2.10186544618431e+07)+pow(logT,6)*(-1.39064488129673e+07)+pow(logT,7)*(2.76972246524054e+06)+pow(logT,8)*(-2.64480903487885e+05)+pow(logT,9)*(1.02334418637311e+04);
				log_ro2=2.26496527471945e+08+logT*(-4.00169957517793e+08)+pow(logT,2)*(2.86967929977283e+08)+pow(logT,3)*(-1.00248810012839e+08)+pow(logT,4)*(1.24947497383337e+07)+pow(logT,5)*(2.97751676794779e+06)+pow(logT,6)*(-1.45346615288982e+06)+pow(logT,7)*(2.47958474490435e+05)+pow(logT,8)*(-2.03246239591655e+04)+pow(logT,9)*(6.63098991105944e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 && 3.447<logT && logT <3.90){
				log_ro1=1.81151505279466e+09+logT*(-3.46617388946230e+09)+pow(logT,2)*(2.72364637563480e+09)+pow(logT,3)*(-1.06827958184809e+09)+pow(logT,4)*(1.68548898720548e+08)+pow(logT,5)*(2.62209713283877e+07)+pow(logT,6)*(-1.76941272831115e+07)+pow(logT,7)*(3.55461245806609e+06)+pow(logT,8)*(-3.41963192370421e+05)+pow(logT,9)*(1.33255524681905e+04);
				log_ro2=1.47056003132479e+09+logT*(-2.79485951795742e+09)+pow(logT,2)*(2.18063596991136e+09)+pow(logT,3)*(-8.48504842120225e+08)+pow(logT,4)*(1.32156159430640e+08)+pow(logT,5)*(2.10186544618431e+07)+pow(logT,6)*(-1.39064488129673e+07)+pow(logT,7)*(2.76972246524054e+06)+pow(logT,8)*(-2.64480903487885e+05)+pow(logT,9)*(1.02334418637311e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 && 3.447<logT && logT <3.875){
				log_ro1=5.13801300966598e+07+logT*(-7.53700700463663e+07)+pow(logT,2)*(3.96796109234423e+07)+pow(logT,3)*(-6.57287025671420e+06)+pow(logT,4)*(-1.42754201312264e+06)+pow(logT,5)*(5.28227495299842e+05)+pow(logT,6)*(6.12860293870630e+04)+pow(logT,7)*(-4.72814729141697e+04)+pow(logT,8)*(7.46395776803570e+03)+pow(logT,9)*(-4.00886552386347e+02);
				log_ro2=1.81151505279466e+09+logT*(-3.46617388946230e+09)+pow(logT,2)*(2.72364637563480e+09)+pow(logT,3)*(-1.06827958184809e+09)+pow(logT,4)*(1.68548898720548e+08)+pow(logT,5)*(2.62209713283877e+07)+pow(logT,6)*(-1.76941272831115e+07)+pow(logT,7)*(3.55461245806609e+06)+pow(logT,8)*(-3.41963192370421e+05)+pow(logT,9)*(1.33255524681905e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 && 3.447<logT && logT <3.875){
				log_ro1=2.10251272794756e+09+logT*(-4.06212244501223e+09)+pow(logT,2)*(3.22453472490813e+09)+pow(logT,3)*(-1.27938489251914e+09)+pow(logT,4)*(2.05766605493639e+08)+pow(logT,5)*(3.08429181092962e+07)+pow(logT,6)*(-2.14644875763487e+07)+pow(logT,7)*(4.36568266030419e+06)+pow(logT,8)*(-4.24372881296649e+05)+pow(logT,9)*(1.66990393612679e+04);
				log_ro2=5.13801300966598e+07+logT*(-7.53700700463663e+07)+pow(logT,2)*(3.96796109234423e+07)+pow(logT,3)*(-6.57287025671420e+06)+pow(logT,4)*(-1.42754201312264e+06)+pow(logT,5)*(5.28227495299842e+05)+pow(logT,6)*(6.12860293870630e+04)+pow(logT,7)*(-4.72814729141697e+04)+pow(logT,8)*(7.46395776803570e+03)+pow(logT,9)*(-4.00886552386347e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 && 3.447<logT && logT <3.875){
				log_ro1=-5.67316906924641e+07+logT*(8.21028672673356e+07)+pow(logT,2)*(-4.25552254445725e+07)+pow(logT,3)*(6.85074236533140e+06)+pow(logT,4)*(1.54685916300144e+06)+pow(logT,5)*(-5.49236534272155e+05)+pow(logT,6)*(-6.69382341199912e+04)+pow(logT,7)*(4.91797290426369e+04)+pow(logT,8)*(-7.63699083355760e+03)+pow(logT,9)*(4.04719281681329e+02);
				log_ro2=2.10251272794756e+09+logT*(-4.06212244501223e+09)+pow(logT,2)*(3.22453472490813e+09)+pow(logT,3)*(-1.27938489251914e+09)+pow(logT,4)*(2.05766605493639e+08)+pow(logT,5)*(3.08429181092962e+07)+pow(logT,6)*(-2.14644875763487e+07)+pow(logT,7)*(4.36568266030419e+06)+pow(logT,8)*(-4.24372881296649e+05)+pow(logT,9)*(1.66990393612679e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 && 3.49<logT && logT <3.778){
				log_ro1=-3.90050144534519e+08+logT*(5.68019395479577e+08)+pow(logT,2)*(-2.94477058217596e+08)+pow(logT,3)*(4.56192281514329e+07)+pow(logT,4)*(1.19721754950763e+07)+pow(logT,5)*(-3.96920142757739e+06)+pow(logT,6)*(-5.74463337338530e+05)+pow(logT,7)*(3.88838287479923e+05)+pow(logT,8)*(-6.02953062388368e+04)+pow(logT,9)*(3.21712047684013e+03);
				log_ro2=-5.67316906924641e+07+logT*(8.21028672673356e+07)+pow(logT,2)*(-4.25552254445725e+07)+pow(logT,3)*(6.85074236533140e+06)+pow(logT,4)*(1.54685916300144e+06)+pow(logT,5)*(-5.49236534272155e+05)+pow(logT,6)*(-6.69382341199912e+04)+pow(logT,7)*(4.91797290426369e+04)+pow(logT,8)*(-7.63699083355760e+03)+pow(logT,9)*(4.04719281681329e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (3.53<logT && logT <3.72){
				log_ro1=6.08329341250592e+07+logT*(-6.25764836109800e+07)+pow(logT,2)*(1.66098828573007e+07)+pow(logT,3)*(2.18748794816387e+06)+pow(logT,4)*(-9.83630563216045e+05)+pow(logT,5)*(-2.22987201258889e+05)+pow(logT,6)*(5.57953484571625e+04)+pow(logT,7)*(2.08049176841391e+04)+pow(logT,8)*(-6.71883567133970e+03)+pow(logT,9)*(5.14162180028259e+02);
				log_ro2=-3.90050144534519e+08+logT*(5.68019395479577e+08)+pow(logT,2)*(-2.94477058217596e+08)+pow(logT,3)*(4.56192281514329e+07)+pow(logT,4)*(1.19721754950763e+07)+pow(logT,5)*(-3.96920142757739e+06)+pow(logT,6)*(-5.74463337338530e+05)+pow(logT,7)*(3.88838287479923e+05)+pow(logT,8)*(-6.02953062388368e+04)+pow(logT,9)*(3.21712047684013e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 ){ //same as above interpolations, but takes the lowest avalaible density for T outside of models range
				log_ro1=6.45541e-08;
				log_ro2=1.49844e-07;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 ){
				log_ro1=1.96611e-08;
				log_ro2=6.45541e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 ){
				log_ro1=8.18373e-09;
				log_ro2=1.96611e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 ){
				log_ro1=6.97149e-09;
				log_ro2=8.18373e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 ){
				log_ro1=1.99814e-09;
				log_ro2=6.97149e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 && 3.447<logT && logT <3.875){
				log_ro1=1.22844e-09;
				log_ro2=1.99814e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 ){
				log_ro1=5.22954e-09;
				log_ro2=1.22844e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else {
				log_ro1=8.16445e-09;
				log_ro2=5.22954e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
}


double ro_ph_z000001341446(double logg, double logT)
{ /* Calculate density in the photosphere for the z=0.000001341446 [M/H]=-4*/
        /*density at kappa_ross=2/3 from MbRCS model atmospheres*/
		double log_ro1,log_ro2,log_ro;

		if (logg>=5.	 && 3.50<logT && logT <3.58){
				log_ro1=-2.84041424132356e+08+logT*(5.63065430833740e+08)+pow(logT,2)*(-4.58502910884758e+08)+pow(logT,3)*(1.86902433183293e+08)+pow(logT,4)*(-3.13067304331913e+07)+pow(logT,5)*(-4.36045655166286e+06)+pow(logT,6)*(3.23858481065936e+06)+pow(logT,7)*(-6.75479670823008e+05)+pow(logT,8)*(6.70090261878473e+04)+pow(logT,9)*(-2.68549711678471e+03);
				log_ro2=1.51649996421058e+02+logT*(-1.89708473282296e+01)+pow(logT,2)*(-1.09688467703929e+01)+pow(logT,3)*(-2.46559849274244e+00)+pow(logT,4)*(-1.83052993030065e-01)+pow(logT,5)*(1.12714751911437e-01)+pow(logT,6)*(6.44707112297275e-02)+pow(logT,7)*(1.85904839317130e-02)+pow(logT,8)*(1.60352404149306e-03)+pow(logT,9)*(-1.95745158853105e-03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 && 3.447<logT && logT <3.89){
				log_ro1=-6.18952530003041e+08+logT*(1.21649575985636e+09)+pow(logT,2)*(-9.82633667822500e+08)+pow(logT,3)*(3.97360739920269e+08)+pow(logT,4)*(-6.58405573680519e+07)+pow(logT,5)*(-9.33069755234477e+06)+pow(logT,6)*(6.81664693914618e+06)+pow(logT,7)*(-1.41296254736884e+06)+pow(logT,8)*(1.39505085749942e+05)+pow(logT,9)*(-5.56859969216931e+03);
				log_ro2=-2.84041424132356e+08+logT*(5.63065430833740e+08)+pow(logT,2)*(-4.58502910884758e+08)+pow(logT,3)*(1.86902433183293e+08)+pow(logT,4)*(-3.13067304331913e+07)+pow(logT,5)*(-4.36045655166286e+06)+pow(logT,6)*(3.23858481065936e+06)+pow(logT,7)*(-6.75479670823008e+05)+pow(logT,8)*(6.70090261878473e+04)+pow(logT,9)*(-2.68549711678471e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 && 3.43<logT && logT <3.89){
				log_ro1=4.75900327793851e+08+logT*(-8.95177158113150e+08)+pow(logT,2)*(6.90748927678393e+08)+pow(logT,3)*(-2.65442711384560e+08)+pow(logT,4)*(4.05821977297623e+07)+pow(logT,5)*(6.58801198572858e+06)+pow(logT,6)*(-4.24402141347515e+06)+pow(logT,7)*(8.32772905696067e+05)+pow(logT,8)*(-7.84117149357229e+04)+pow(logT,9)*(2.99126686419748e+03);
				log_ro2=-6.18952530003041e+08+logT*(1.21649575985636e+09)+pow(logT,2)*(-9.82633667822500e+08)+pow(logT,3)*(3.97360739920269e+08)+pow(logT,4)*(-6.58405573680519e+07)+pow(logT,5)*(-9.33069755234477e+06)+pow(logT,6)*(6.81664693914618e+06)+pow(logT,7)*(-1.41296254736884e+06)+pow(logT,8)*(1.39505085749942e+05)+pow(logT,9)*(-5.56859969216931e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 && 3.43<logT && logT <3.90){
				log_ro1=1.63335445531054e+09+logT*(-3.13042690943983e+09)+pow(logT,2)*(2.46466744830715e+09)+pow(logT,3)*(-9.69496223194759e+08)+pow(logT,4)*(1.54214693812367e+08)+pow(logT,5)*(2.32743766946433e+07)+pow(logT,6)*(-1.59586547628222e+07)+pow(logT,7)*(3.21702520603091e+06)+pow(logT,8)*(-3.10133681809411e+05)+pow(logT,9)*(1.21053269904700e+04);
				log_ro2=4.75900327793851e+08+logT*(-8.95177158113150e+08)+pow(logT,2)*(6.90748927678393e+08)+pow(logT,3)*(-2.65442711384560e+08)+pow(logT,4)*(4.05821977297623e+07)+pow(logT,5)*(6.58801198572858e+06)+pow(logT,6)*(-4.24402141347515e+06)+pow(logT,7)*(8.32772905696067e+05)+pow(logT,8)*(-7.84117149357229e+04)+pow(logT,9)*(2.99126686419748e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 && 3.43<logT && logT <3.875){
				log_ro1=3.99441892084331e+07+logT*(-5.88357194600175e+07)+pow(logT,2)*(3.11323803906666e+07)+pow(logT,3)*(-5.21354509427317e+06)+pow(logT,4)*(-1.10981653936927e+06)+pow(logT,5)*(4.17746297710305e+05)+pow(logT,6)*(4.72577007223569e+04)+pow(logT,7)*(-3.72372138931161e+04)+pow(logT,8)*(5.91100023926400e+03)+pow(logT,9)*(-3.18786786811371e+02);
				log_ro2=1.63335445531054e+09+logT*(-3.13042690943983e+09)+pow(logT,2)*(2.46466744830715e+09)+pow(logT,3)*(-9.69496223194759e+08)+pow(logT,4)*(1.54214693812367e+08)+pow(logT,5)*(2.32743766946433e+07)+pow(logT,6)*(-1.59586547628222e+07)+pow(logT,7)*(3.21702520603091e+06)+pow(logT,8)*(-3.10133681809411e+05)+pow(logT,9)*(1.21053269904700e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 && 3.43<logT && logT <3.875){
				log_ro1=1.13161962110323e+09+logT*(-2.19708877565031e+09)+pow(logT,2)*(1.75326185190402e+09)+pow(logT,3)*(-7.00131736928229e+08)+pow(logT,4)*(1.14142421968394e+08)+pow(logT,5)*(1.64145014674356e+07)+pow(logT,6)*(-1.17200150166782e+07)+pow(logT,7)*(2.40029230642886e+06)+pow(logT,8)*(-2.34463630390798e+05)+pow(logT,9)*(9.26463440344015e+03);
				log_ro2=3.99441892084331e+07+logT*(-5.88357194600175e+07)+pow(logT,2)*(3.11323803906666e+07)+pow(logT,3)*(-5.21354509427317e+06)+pow(logT,4)*(-1.10981653936927e+06)+pow(logT,5)*(4.17746297710305e+05)+pow(logT,6)*(4.72577007223569e+04)+pow(logT,7)*(-3.72372138931161e+04)+pow(logT,8)*(5.91100023926400e+03)+pow(logT,9)*(-3.18786786811371e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 && 3.40<logT && logT <3.86){
				log_ro1=3.15382393218734e+09+logT*(-6.15153174631909e+09)+pow(logT,2)*(4.93067262535992e+09)+pow(logT,3)*(-1.97638774480148e+09)+pow(logT,4)*(3.22067815457358e+08)+pow(logT,5)*(4.78470617501024e+07)+pow(logT,6)*(-3.39004898479258e+07)+pow(logT,7)*(6.96798741930292e+06)+pow(logT,8)*(-6.83960012245062e+05)+pow(logT,9)*(2.71702142351027e+04);
				log_ro2=1.13161962110323e+09+logT*(-2.19708877565031e+09)+pow(logT,2)*(1.75326185190402e+09)+pow(logT,3)*(-7.00131736928229e+08)+pow(logT,4)*(1.14142421968394e+08)+pow(logT,5)*(1.64145014674356e+07)+pow(logT,6)*(-1.17200150166782e+07)+pow(logT,7)*(2.40029230642886e+06)+pow(logT,8)*(-2.34463630390798e+05)+pow(logT,9)*(9.26463440344015e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 && 3.46<logT && logT <3.829){
				log_ro1=-1.82384259987877e+08+logT*(2.66338099925826e+08)+pow(logT,2)*(-1.39047040681826e+08)+pow(logT,3)*(2.22948790344844e+07)+pow(logT,4)*(5.30567872230205e+06)+pow(logT,5)*(-1.85795450024108e+06)+pow(logT,6)*(-2.40406713866218e+05)+pow(logT,7)*(1.73317432416182e+05)+pow(logT,8)*(-2.70929486269498e+04)+pow(logT,9)*(1.44897649633218e+03);
				log_ro2=3.15382393218734e+09+logT*(-6.15153174631909e+09)+pow(logT,2)*(4.93067262535992e+09)+pow(logT,3)*(-1.97638774480148e+09)+pow(logT,4)*(3.22067815457358e+08)+pow(logT,5)*(4.78470617501024e+07)+pow(logT,6)*(-3.39004898479258e+07)+pow(logT,7)*(6.96798741930292e+06)+pow(logT,8)*(-6.83960012245062e+05)+pow(logT,9)*(2.71702142351027e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 && 3.477<logT && logT <3.74){
				log_ro1=1.37442023920782e+08+logT*(-1.99716639890121e+08)+pow(logT,2)*(1.03100823127357e+08)+pow(logT,3)*(-1.56566779540175e+07)+pow(logT,4)*(-4.33863950128744e+06)+pow(logT,5)*(1.39782668228052e+06)+pow(logT,6)*(2.14610477771278e+05)+pow(logT,7)*(-1.40974269151719e+05)+pow(logT,8)*(2.18063876902107e+04)+pow(logT,9)*(-1.16451135125656e+03);
				log_ro2=-1.82384259987877e+08+logT*(2.66338099925826e+08)+pow(logT,2)*(-1.39047040681826e+08)+pow(logT,3)*(2.22948790344844e+07)+pow(logT,4)*(5.30567872230205e+06)+pow(logT,5)*(-1.85795450024108e+06)+pow(logT,6)*(-2.40406713866218e+05)+pow(logT,7)*(1.73317432416182e+05)+pow(logT,8)*(-2.70929486269498e+04)+pow(logT,9)*(1.44897649633218e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if ( 3.676<logT && logT <3.72){
				log_ro1=6.69258095990584e+00+logT*(4.08075042146622e-01)+pow(logT,2)*(-1.66485639157397e-01)+pow(logT,3)*(-9.31751552667997e-02)+pow(logT,4)*(-3.12422781323950e-02)+pow(logT,5)*(-8.26204718559186e-03)+pow(logT,6)*(-1.70840712215930e-03)+pow(logT,7)*(-1.95770271571331e-04)+pow(logT,8)*(5.13949559033553e-05)+pow(logT,9)*(5.05441058109099e-05);
				log_ro2=1.37442023920782e+08+logT*(-1.99716639890121e+08)+pow(logT,2)*(1.03100823127357e+08)+pow(logT,3)*(-1.56566779540175e+07)+pow(logT,4)*(-4.33863950128744e+06)+pow(logT,5)*(1.39782668228052e+06)+pow(logT,6)*(2.14610477771278e+05)+pow(logT,7)*(-1.40974269151719e+05)+pow(logT,8)*(2.18063876902107e+04)+pow(logT,9)*(-1.16451135125656e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}


		else if (logg>=5.){
				log_ro1=1.48954e-07;
				log_ro2=0.000184299;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 ){ //same as above interpolations, but takes the lowest avalaible density for T outside of models range
				log_ro1=6.41378e-08;
				log_ro2=1.48954e-07;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0){
				log_ro1=1.94777e-08;
				log_ro2=6.41378e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5){
				log_ro1=8.05169e-09;
				log_ro2=1.94777e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0){
				log_ro1=6.91006e-09;
				log_ro2=8.05169e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5){
				log_ro1=1.28912e-09;
				log_ro2=6.91006e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0){
				log_ro1=1.81527e-09;
				log_ro2=1.28912e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 ){
				log_ro1=1.6886e-09;
				log_ro2=1.81527e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 ){
				log_ro1=5.25196e-09;
				log_ro2=1.6886e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else {
				log_ro1=3.65432e-09;
				log_ro2=5.25196e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}


}


double ro_ph_z0000134117(double logg, double logT)
{ /* Calculate density in the photosphere for the z=0.0000134117 [M/H]=-3*/
        /*density at kappa_ross=2/3 from MbRCS model atmospheres*/
		double log_ro1,log_ro2,log_ro;

		if (logg>=5. && 3.544<logT && logT <3.59){
				log_ro1=-2.64412012920675e+07+logT*(3.72920688870408e+07)+pow(logT,2)*(-1.87471226547541e+07)+pow(logT,3)*(2.83809310927833e+06)+pow(logT,4)*(7.02008199409796e+05)+pow(logT,5)*(-2.28657751309118e+05)+pow(logT,6)*(-3.10881205673838e+04)+pow(logT,7)*(2.07140936223670e+04)+pow(logT,8)*(-3.11679686968858e+03)+pow(logT,9)*(1.61132779430279e+02);
				log_ro2=3.77444238684283e+04+logT*(-1.32240911166001e+04)+pow(logT,2)*(-2.74348292617446e+03)+pow(logT,3)*(1.69588946884629e+02)+pow(logT,4)*(2.35344638264921e+02)+pow(logT,5)*(6.33380360659607e+01)+pow(logT,6)*(2.32325640407499e+00)+pow(logT,7)*(-4.84353865696523e+00)+pow(logT,8)*(-1.68478244009079e+00)+pow(logT,9)*(3.94867870662643e-01);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 && 3.518<logT && logT <3.889){
				log_ro1=2.32646129659295e+09+logT*(-4.44361672177726e+09)+pow(logT,2)*(3.48542348608965e+09)+pow(logT,3)*(-1.36438403525153e+09)+pow(logT,4)*(2.14600473947871e+08)+pow(logT,5)*(3.35664185896256e+07)+pow(logT,6)*(-2.25449233886313e+07)+pow(logT,7)*(4.52036101037484e+06)+pow(logT,8)*(-4.34176726706390e+05)+pow(logT,9)*(1.68941172159142e+04);
				log_ro2=-2.64412012920675e+07+logT*(3.72920688870408e+07)+pow(logT,2)*(-1.87471226547541e+07)+pow(logT,3)*(2.83809310927833e+06)+pow(logT,4)*(7.02008199409796e+05)+pow(logT,5)*(-2.28657751309118e+05)+pow(logT,6)*(-3.10881205673838e+04)+pow(logT,7)*(2.07140936223670e+04)+pow(logT,8)*(-3.11679686968858e+03)+pow(logT,9)*(1.61132779430279e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 && 3.518<logT && logT <3.90){
				log_ro1=2.71136444572400e+09+logT*(-5.19671345626987e+09)+pow(logT,2)*(4.09127867791688e+09)+pow(logT,3)*(-1.60872815496765e+09)+pow(logT,4)*(2.55302208498302e+08)+pow(logT,5)*(3.90022432934040e+07)+pow(logT,6)*(-2.66021571746608e+07)+pow(logT,7)*(5.35986784683926e+06)+pow(logT,8)*(-5.16722227568852e+05)+pow(logT,9)*(2.01730301940591e+04);
				log_ro2=2.32646129659295e+09+logT*(-4.44361672177726e+09)+pow(logT,2)*(3.48542348608965e+09)+pow(logT,3)*(-1.36438403525153e+09)+pow(logT,4)*(2.14600473947871e+08)+pow(logT,5)*(3.35664185896256e+07)+pow(logT,6)*(-2.25449233886313e+07)+pow(logT,7)*(4.52036101037484e+06)+pow(logT,8)*(-4.34176726706390e+05)+pow(logT,9)*(1.68941172159142e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 && 3.40<logT && logT <3.90){
				log_ro1=1.02777800467622e+09+logT*(-1.97645373358802e+09)+pow(logT,2)*(1.56182467914930e+09)+pow(logT,3)*(-6.17258637414476e+08)+pow(logT,4)*(9.92935538473600e+07)+pow(logT,5)*(1.44169063079539e+07)+pow(logT,6)*(-1.01029945587868e+07)+pow(logT,7)*(2.04679210760763e+06)+pow(logT,8)*(-1.97931117867726e+05)+pow(logT,9)*(7.74448610861221e+03);
				log_ro2=2.71136444572400e+09+logT*(-5.19671345626987e+09)+pow(logT,2)*(4.09127867791688e+09)+pow(logT,3)*(-1.60872815496765e+09)+pow(logT,4)*(2.55302208498302e+08)+pow(logT,5)*(3.90022432934040e+07)+pow(logT,6)*(-2.66021571746608e+07)+pow(logT,7)*(5.35986784683926e+06)+pow(logT,8)*(-5.16722227568852e+05)+pow(logT,9)*(2.01730301940591e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 && 3.40<logT && logT <3.875){
				log_ro1=1.96830136786164e+09+logT*(-3.80872514888106e+09)+pow(logT,2)*(3.02860834254621e+09)+pow(logT,3)*(-1.20444910598624e+09)+pow(logT,4)*(1.94880662220742e+08)+pow(logT,5)*(2.85702836603124e+07)+pow(logT,6)*(-2.01240253665703e+07)+pow(logT,7)*(4.10344503723085e+06)+pow(logT,8)*(-3.99474771342413e+05)+pow(logT,9)*(1.57368909994513e+04);
				log_ro2=1.02777800467622e+09+logT*(-1.97645373358802e+09)+pow(logT,2)*(1.56182467914930e+09)+pow(logT,3)*(-6.17258637414476e+08)+pow(logT,4)*(9.92935538473600e+07)+pow(logT,5)*(1.44169063079539e+07)+pow(logT,6)*(-1.01029945587868e+07)+pow(logT,7)*(2.04679210760763e+06)+pow(logT,8)*(-1.97931117867726e+05)+pow(logT,9)*(7.74448610861221e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 && 3.40<logT && logT <3.875){
				log_ro1=1.89902555509771e+09+logT*(-3.68657917785092e+09)+pow(logT,2)*(2.94228346915685e+09)+pow(logT,3)*(-1.17596441651393e+09)+pow(logT,4)*(1.92646717437463e+08)+pow(logT,5)*(2.69931273368581e+07)+pow(logT,6)*(-1.95090701128040e+07)+pow(logT,7)*(4.00126090880617e+06)+pow(logT,8)*(-3.90985433170492e+05)+pow(logT,9)*(1.54497911232653e+04);
				log_ro2=1.96830136786164e+09+logT*(-3.80872514888106e+09)+pow(logT,2)*(3.02860834254621e+09)+pow(logT,3)*(-1.20444910598624e+09)+pow(logT,4)*(1.94880662220742e+08)+pow(logT,5)*(2.85702836603124e+07)+pow(logT,6)*(-2.01240253665703e+07)+pow(logT,7)*(4.10344503723085e+06)+pow(logT,8)*(-3.99474771342413e+05)+pow(logT,9)*(1.57368909994513e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 && 3.40<logT && logT <3.90){
				log_ro1=1.61824289790862e+09+logT*(-3.15741310036298e+09)+pow(logT,2)*(2.53353745333078e+09)+pow(logT,3)*(-1.01905832285274e+09)+pow(logT,4)*(1.68953377161032e+08)+pow(logT,5)*(2.28994588317158e+07)+pow(logT,6)*(-1.69375261196335e+07)+pow(logT,7)*(3.49819136849617e+06)+pow(logT,8)*(-3.43653661359311e+05)+pow(logT,9)*(1.36446322088121e+04);
				log_ro2=1.89902555509771e+09+logT*(-3.68657917785092e+09)+pow(logT,2)*(2.94228346915685e+09)+pow(logT,3)*(-1.17596441651393e+09)+pow(logT,4)*(1.92646717437463e+08)+pow(logT,5)*(2.69931273368581e+07)+pow(logT,6)*(-1.95090701128040e+07)+pow(logT,7)*(4.00126090880617e+06)+pow(logT,8)*(-3.90985433170492e+05)+pow(logT,9)*(1.54497911232653e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 && 3.40<logT && logT <3.829){
				log_ro1=3.69604199294176e+09+logT*(-7.24757169614758e+09)+pow(logT,2)*(5.83807569643517e+09)+pow(logT,3)*(-2.34941526190363e+09)+pow(logT,4)*(3.82253850094554e+08)+pow(logT,5)*(5.90922537405433e+07)+pow(logT,6)*(-4.14303368632767e+07)+pow(logT,7)*(8.54365610321481e+06)+pow(logT,8)*(-8.42597192760214e+05)+pow(logT,9)*(3.36456482766061e+04);
				log_ro2=1.61824289790862e+09+logT*(-3.15741310036298e+09)+pow(logT,2)*(2.53353745333078e+09)+pow(logT,3)*(-1.01905832285274e+09)+pow(logT,4)*(1.68953377161032e+08)+pow(logT,5)*(2.28994588317158e+07)+pow(logT,6)*(-1.69375261196335e+07)+pow(logT,7)*(3.49819136849617e+06)+pow(logT,8)*(-3.43653661359311e+05)+pow(logT,9)*(1.36446322088121e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 && 3.40<logT && logT <3.77){
				log_ro1=-4.73439415142100e+07+logT*(6.82616206744696e+07)+pow(logT,2)*(-3.50259108503443e+07)+pow(logT,3)*(5.37629415944336e+06)+pow(logT,4)*(1.38583673142828e+06)+pow(logT,5)*(-4.54901603281020e+05)+pow(logT,6)*(-6.48002349701028e+04)+pow(logT,7)*(4.33336236406854e+04)+pow(logT,8)*(-6.62469294442309e+03)+pow(logT,9)*(3.48181815695679e+02);
				log_ro2=3.69604199294176e+09+logT*(-7.24757169614758e+09)+pow(logT,2)*(5.83807569643517e+09)+pow(logT,3)*(-2.34941526190363e+09)+pow(logT,4)*(3.82253850094554e+08)+pow(logT,5)*(5.90922537405433e+07)+pow(logT,6)*(-4.14303368632767e+07)+pow(logT,7)*(8.54365610321481e+06)+pow(logT,8)*(-8.42597192760214e+05)+pow(logT,9)*(3.36456482766061e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if ( 3.49<logT && logT <3.74){
				log_ro1=2.00995741069017e+08+logT*(-2.96787982601042e+08)+pow(logT,2)*(1.55885755299594e+08)+pow(logT,3)*(-2.43265319917939e+07)+pow(logT,4)*(-6.60684926170953e+06)+pow(logT,5)*(2.19731463946853e+06)+pow(logT,6)*(3.31244747898281e+05)+pow(logT,7)*(-2.24435725485319e+05)+pow(logT,8)*(3.52629046421978e+04)+pow(logT,9)*(-1.90898169209299e+03);
				log_ro2=-4.73439415142100e+07+logT*(6.82616206744696e+07)+pow(logT,2)*(-3.50259108503443e+07)+pow(logT,3)*(5.37629415944336e+06)+pow(logT,4)*(1.38583673142828e+06)+pow(logT,5)*(-4.54901603281020e+05)+pow(logT,6)*(-6.48002349701028e+04)+pow(logT,7)*(4.33336236406854e+04)+pow(logT,8)*(-6.62469294442309e+03)+pow(logT,9)*(3.48181815695679e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=5. ){
				log_ro1=1.47117e-07;
				log_ro2=9.4015e-05;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5){
				log_ro1=4.64114e-08;
				log_ro2=1.47117e-07;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0){
				log_ro1=1.92753e-08;
				log_ro2=4.64114e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5){
				log_ro1=7.91219e-09;
				log_ro2=1.92753e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 ){
				log_ro1=6.82139e-09;
				log_ro2=7.91219e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 ){
				log_ro1=1.26895e-09;
				log_ro2=6.82139e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 ){
				log_ro1= 4.71465e-10;
				log_ro2=1.26895e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 ){
				log_ro1=1.66401e-09;
				log_ro2= 4.71465e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0){
				log_ro1=2.371e-09;
				log_ro2=1.66401e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else {
				log_ro1=2.34184e-09;
				log_ro2=2.371e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}


}


double ro_ph_z0000423904(double logg, double logT)
{ /* Calculate density in the photosphere for the z=0.0000423904 [M/H]=-2.5*/
        /*density at kappa_ross=2/3 from MbRCS model atmospheres*/
		double log_ro1,log_ro2,log_ro;

		if (logg>=5. 	 && 3.447<logT && logT <3.59){
				log_ro1=-1.34153837479747e+09+logT*(2.61319987024629e+09)+pow(logT,2)*(-2.09233932607284e+09)+pow(logT,3)*(8.38453378371340e+08)+pow(logT,4)*(-1.37223996361958e+08)+pow(logT,5)*(-1.97613092503981e+07)+pow(logT,6)*(1.41745360508961e+07)+pow(logT,7)*(-2.91394388541560e+06)+pow(logT,8)*(2.85702664121162e+05)+pow(logT,9)*(-1.13318527984974e+04);
				log_ro2=1.53462052738403e+07+logT*(-1.62652955306008e+07)+pow(logT,2)*(4.41961245337978e+06)+pow(logT,3)*(6.26501325958978e+05)+pow(logT,4)*(-2.80897852332271e+05)+pow(logT,5)*(-6.80783054049728e+04)+pow(logT,6)*(1.70426968182492e+04)+pow(logT,7)*(6.79876708626968e+03)+pow(logT,8)*(-2.24365126670075e+03)+pow(logT,9)*(1.76887517269456e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 && 3.431<logT && logT <3.89){
				log_ro1=6.60787983416220e+08+logT*(-1.24887304690417e+09)+pow(logT,2)*(9.68518755251412e+08)+pow(logT,3)*(-3.74288586689134e+08)+pow(logT,4)*(5.77329643768065e+07)+pow(logT,5)*(9.25145223522582e+06)+pow(logT,6)*(-6.03759657072898e+06)+pow(logT,7)*(1.19218902528661e+06)+pow(logT,8)*(-1.12888238642684e+05)+pow(logT,9)*(4.33029884440939e+03);
				log_ro2=-1.34153837479747e+09+logT*(2.61319987024629e+09)+pow(logT,2)*(-2.09233932607284e+09)+pow(logT,3)*(8.38453378371340e+08)+pow(logT,4)*(-1.37223996361958e+08)+pow(logT,5)*(-1.97613092503981e+07)+pow(logT,6)*(1.41745360508961e+07)+pow(logT,7)*(-2.91394388541560e+06)+pow(logT,8)*(2.85702664121162e+05)+pow(logT,9)*(-1.13318527984974e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 && 3.431<logT && logT <3.90){
				log_ro1=2.08494725075866e+08+logT*(-3.82795052051248e+08)+pow(logT,2)*(2.87134736691503e+08)+pow(logT,3)*(-1.06458870009644e+08)+pow(logT,4)*(1.52105360422853e+07)+pow(logT,5)*(2.77817713197874e+06)+pow(logT,6)*(-1.61466629881662e+06)+pow(logT,7)*(3.00500061397849e+05)+pow(logT,8)*(-2.68640320580824e+04)+pow(logT,9)*(9.68887764278117e+02);
				log_ro2=6.60787983416220e+08+logT*(-1.24887304690417e+09)+pow(logT,2)*(9.68518755251412e+08)+pow(logT,3)*(-3.74288586689134e+08)+pow(logT,4)*(5.77329643768065e+07)+pow(logT,5)*(9.25145223522582e+06)+pow(logT,6)*(-6.03759657072898e+06)+pow(logT,7)*(1.19218902528661e+06)+pow(logT,8)*(-1.12888238642684e+05)+pow(logT,9)*(4.33029884440939e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 && 3.40<logT && logT <3.875){
				log_ro1=9.59831121384848e+08+logT*(-1.84233211326044e+09)+pow(logT,2)*(1.45281219852882e+09)+pow(logT,3)*(-5.72720941826118e+08)+pow(logT,4)*(9.16923733823480e+07)+pow(logT,5)*(1.34659495303048e+07)+pow(logT,6)*(-9.35658957572846e+06)+pow(logT,7)*(1.88970860861918e+06)+pow(logT,8)*(-1.82265368855739e+05)+pow(logT,9)*(7.11368861745044e+03);
				log_ro2=2.08494725075866e+08+logT*(-3.82795052051248e+08)+pow(logT,2)*(2.87134736691503e+08)+pow(logT,3)*(-1.06458870009644e+08)+pow(logT,4)*(1.52105360422853e+07)+pow(logT,5)*(2.77817713197874e+06)+pow(logT,6)*(-1.61466629881662e+06)+pow(logT,7)*(3.00500061397849e+05)+pow(logT,8)*(-2.68640320580824e+04)+pow(logT,9)*(9.68887764278117e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 && 3.40<logT && logT <3.875){
				log_ro1=1.04712227769785e+09+logT*(-2.01228940959711e+09)+pow(logT,2)*(1.58822813959535e+09)+pow(logT,3)*(-6.26134732175276e+08)+pow(logT,4)*(9.98020857178538e+07)+pow(logT,5)*(1.50807837443239e+07)+pow(logT,6)*(-1.03576892729453e+07)+pow(logT,7)*(2.09047296549717e+06)+pow(logT,8)*(-2.01714264221231e+05)+pow(logT,9)*(7.87843205077166e+03);
				log_ro2=9.59831121384848e+08+logT*(-1.84233211326044e+09)+pow(logT,2)*(1.45281219852882e+09)+pow(logT,3)*(-5.72720941826118e+08)+pow(logT,4)*(9.16923733823480e+07)+pow(logT,5)*(1.34659495303048e+07)+pow(logT,6)*(-9.35658957572846e+06)+pow(logT,7)*(1.88970860861918e+06)+pow(logT,8)*(-1.82265368855739e+05)+pow(logT,9)*(7.11368861745044e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 && 3.40<logT && logT <3.875){
				log_ro1=1.53580908088827e+09+logT*(-2.97813797408398e+09)+pow(logT,2)*(2.37399456581823e+09)+pow(logT,3)*(-9.47480880865258e+08)+pow(logT,4)*(1.54831608146973e+08)+pow(logT,5)*(2.18157623831985e+07)+pow(logT,6)*(-1.56966502036681e+07)+pow(logT,7)*(3.21394635909716e+06)+pow(logT,8)*(-3.13607449648730e+05)+pow(logT,9)*(1.23753878432468e+04);
				log_ro2=1.04712227769785e+09+logT*(-2.01228940959711e+09)+pow(logT,2)*(1.58822813959535e+09)+pow(logT,3)*(-6.26134732175276e+08)+pow(logT,4)*(9.98020857178538e+07)+pow(logT,5)*(1.50807837443239e+07)+pow(logT,6)*(-1.03576892729453e+07)+pow(logT,7)*(2.09047296549717e+06)+pow(logT,8)*(-2.01714264221231e+05)+pow(logT,9)*(7.87843205077166e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 && 3.40<logT && logT <3.875){
				log_ro1=1.07772252870003e+09+logT*(-2.09916714896423e+09)+pow(logT,2)*(1.68026975368127e+09)+pow(logT,3)*(-6.72858292535950e+08)+pow(logT,4)*(1.09862968043288e+08)+pow(logT,5)*(1.59695610364872e+07)+pow(logT,6)*(-1.13937672234275e+07)+pow(logT,7)*(2.33901888938685e+06)+pow(logT,8)*(-2.29081160039768e+05)+pow(logT,9)*(9.07615039724846e+03);
				log_ro2=1.53580908088827e+09+logT*(-2.97813797408398e+09)+pow(logT,2)*(2.37399456581823e+09)+pow(logT,3)*(-9.47480880865258e+08)+pow(logT,4)*(1.54831608146973e+08)+pow(logT,5)*(2.18157623831985e+07)+pow(logT,6)*(-1.56966502036681e+07)+pow(logT,7)*(3.21394635909716e+06)+pow(logT,8)*(-3.13607449648730e+05)+pow(logT,9)*(1.23753878432468e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 && 3.40<logT && logT <3.845){
				log_ro1=-8.45684672526155e+07+logT*(1.62754238216092e+08)+pow(logT,2)*(-1.28970141943845e+08)+pow(logT,3)*(5.11698412947106e+07)+pow(logT,4)*(-8.21296717035675e+06)+pow(logT,5)*(-1.27488658525130e+06)+pow(logT,6)*(8.83820123725904e+05)+pow(logT,7)*(-1.81603510687343e+05)+pow(logT,8)*(1.79039548316565e+04)+pow(logT,9)*(-7.16711076583342e+02);
				log_ro2=1.07772252870003e+09+logT*(-2.09916714896423e+09)+pow(logT,2)*(1.68026975368127e+09)+pow(logT,3)*(-6.72858292535950e+08)+pow(logT,4)*(1.09862968043288e+08)+pow(logT,5)*(1.59695610364872e+07)+pow(logT,6)*(-1.13937672234275e+07)+pow(logT,7)*(2.33901888938685e+06)+pow(logT,8)*(-2.29081160039768e+05)+pow(logT,9)*(9.07615039724846e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 && 3.40<logT && logT <3.796){
				log_ro1=-2.29098195636935e+07+logT*(3.34507404151832e+07)+pow(logT,2)*(-1.74484357206619e+07)+pow(logT,3)*(2.78473711769804e+06)+pow(logT,4)*(6.71490531158088e+05)+pow(logT,5)*(-2.33279197237067e+05)+pow(logT,6)*(-3.06669644505166e+04)+pow(logT,7)*(2.18931703434023e+04)+pow(logT,8)*(-3.41653602011542e+03)+pow(logT,9)*(1.82553328623131e+02);
				log_ro2=-8.45684672526155e+07+logT*(1.62754238216092e+08)+pow(logT,2)*(-1.28970141943845e+08)+pow(logT,3)*(5.11698412947106e+07)+pow(logT,4)*(-8.21296717035675e+06)+pow(logT,5)*(-1.27488658525130e+06)+pow(logT,6)*(8.83820123725904e+05)+pow(logT,7)*(-1.81603510687343e+05)+pow(logT,8)*(1.79039548316565e+04)+pow(logT,9)*(-7.16711076583342e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5 && 3.40<logT && logT <3.74){
				log_ro1=-1.24653527290972e+08+logT*(1.85585191586531e+08)+pow(logT,2)*(-9.86938573910418e+07)+pow(logT,3)*(1.60234787633453e+07)+pow(logT,4)*(3.97910463531006e+06)+pow(logT,5)*(-1.40433812265045e+06)+pow(logT,6)*(-1.90255700882485e+05)+pow(logT,7)*(1.37965277747599e+05)+pow(logT,8)*(-2.19761121847792e+04)+pow(logT,9)*(1.19931852134230e+03);
				log_ro2=-2.29098195636935e+07+logT*(3.34507404151832e+07)+pow(logT,2)*(-1.74484357206619e+07)+pow(logT,3)*(2.78473711769804e+06)+pow(logT,4)*(6.71490531158088e+05)+pow(logT,5)*(-2.33279197237067e+05)+pow(logT,6)*(-3.06669644505166e+04)+pow(logT,7)*(2.18931703434023e+04)+pow(logT,8)*(-3.41653602011542e+03)+pow(logT,9)*(1.82553328623131e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if ( 3.40<logT && logT <3.53){
				log_ro1=-2.57906931065356e+01+logT*(-2.67997421573416e+00)+pow(logT,2)*(2.13667982097364e-01)+pow(logT,3)*(2.42728836883762e-01)+pow(logT,4)*(9.45444039037490e-02)+pow(logT,5)*(2.69435219557533e-02)+pow(logT,6)*(5.70299977750937e-03)+pow(logT,7)*(5.22858842616709e-04)+pow(logT,8)*(-3.12376739535123e-04)+pow(logT,9)*(-2.60661797349995e-04);
				log_ro2=-1.24653527290972e+08+logT*(1.85585191586531e+08)+pow(logT,2)*(-9.86938573910418e+07)+pow(logT,3)*(1.60234787633453e+07)+pow(logT,4)*(3.97910463531006e+06)+pow(logT,5)*(-1.40433812265045e+06)+pow(logT,6)*(-1.90255700882485e+05)+pow(logT,7)*(1.37965277747599e+05)+pow(logT,8)*(-2.19761121847792e+04)+pow(logT,9)*(1.19931852134230e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=5.){
				log_ro1=1.45148e-07;
				log_ro2=6.19146e-05;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5){
				log_ro1= 4.59313e-08;
				log_ro2=1.45148e-07;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0){
				log_ro1=1.90742e-08;
				log_ro2=4.59313e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 ){
				log_ro1=7.77521e-09;
				log_ro2=1.90742e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 ){
				log_ro1=6.73275e-09;
				log_ro2=7.77521e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 ){
				log_ro1=1.24856e-09;
				log_ro2=6.73275e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 ){
				log_ro1=1.16234e-09;
				log_ro2=1.24856e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 ){
				log_ro1=1.12224e-09;
				log_ro2=1.16234e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 ){
				log_ro1=1.55187e-09;
				log_ro2=1.12224e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5){
				log_ro1=2.30214e-09;
				log_ro2=1.55187e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else {
				log_ro1=1.01267e-08;
				log_ro2=2.30214e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}

}


double ro_ph_z0001338401(double logg, double logT)
{ /* Calculate density in the photosphere for the z=0.0001338401 [M/H]=-2*/
        /*density at kappa_ross=2/3 from MbRCS model atmospheres*/
		double log_ro1,log_ro2,log_ro;

		if (logg>=5.	 && 3.40<logT && logT <3.59){
				log_ro1=5.30366093813980e+09+logT*(-7.98514736118541e+09)+pow(logT,2)*(4.25938130073495e+09)+pow(logT,3)*(-6.57180781608924e+08)+pow(logT,4)*(-1.99018845255812e+08)+pow(logT,5)*(6.44791280945082e+07)+pow(logT,6)*(1.09310785681452e+07)+pow(logT,7)*(-7.19431354751120e+06)+pow(logT,8)*(1.14684247141821e+06)+pow(logT,9)*(-6.32971275959849e+04);
				log_ro2=-1.29153010654010e+09+logT*(2.52618215239678e+09)+pow(logT,2)*(-2.03141601370292e+09)+pow(logT,3)*(8.18116790075141e+08)+pow(logT,4)*(-1.35119831274356e+08)+pow(logT,5)*(-1.90041444370151e+07)+pow(logT,6)*(1.38613386735035e+07)+pow(logT,7)*(-2.86459936539110e+06)+pow(logT,8)*(2.81999376043980e+05)+pow(logT,9)*(-1.12253359476931e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 && 3.40<logT && logT <3.89){
				log_ro1=-1.53838290252781e+09+logT*(3.00726627122255e+09)+pow(logT,2)*(-2.41704644587486e+09)+pow(logT,3)*(9.72995935802510e+08)+pow(logT,4)*(-1.60640212964903e+08)+pow(logT,5)*(-2.25866054129733e+07)+pow(logT,6)*(1.64704025638955e+07)+pow(logT,7)*(-3.40318861739694e+06)+pow(logT,8)*(3.34980398751799e+05)+pow(logT,9)*(-1.33335208996748e+04);
				log_ro2=-1.29153010654010e+09+logT*(2.52618215239678e+09)+pow(logT,2)*(-2.03141601370292e+09)+pow(logT,3)*(8.18116790075141e+08)+pow(logT,4)*(-1.35119831274356e+08)+pow(logT,5)*(-1.90041444370151e+07)+pow(logT,6)*(1.38613386735035e+07)+pow(logT,7)*(-2.86459936539110e+06)+pow(logT,8)*(2.81999376043980e+05)+pow(logT,9)*(-1.12253359476931e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 && 3.40<logT && logT <3.89){
				log_ro1=-2.43319363357935e+08+logT*(4.94943297366238e+08)+pow(logT,2)*(-4.14012682765776e+08)+pow(logT,3)*(1.74094893959945e+08)+pow(logT,4)*(-3.08472205494937e+07)+pow(logT,5)*(-3.64541790523608e+06)+pow(logT,6)*(3.04644643242255e+06)+pow(logT,7)*(-6.57068745752758e+05)+pow(logT,8)*(6.68463038097033e+04)+pow(logT,9)*(-2.73970549044038e+03);
				log_ro2=-1.53838290252781e+09+logT*(3.00726627122255e+09)+pow(logT,2)*(-2.41704644587486e+09)+pow(logT,3)*(9.72995935802510e+08)+pow(logT,4)*(-1.60640212964903e+08)+pow(logT,5)*(-2.25866054129733e+07)+pow(logT,6)*(1.64704025638955e+07)+pow(logT,7)*(-3.40318861739694e+06)+pow(logT,8)*(3.34980398751799e+05)+pow(logT,9)*(-1.33335208996748e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 && 3.40<logT && logT <3.90){
				log_ro1=4.28973981187798e+08+logT*(-8.13647468051304e+08)+pow(logT,2)*(6.33290174409065e+08)+pow(logT,3)*(-2.45804555993348e+08)+pow(logT,4)*(3.82947894887328e+07)+pow(logT,5)*(5.93978135372895e+06)+pow(logT,6)*(-3.94367912591377e+06)+pow(logT,7)*(7.81467947170254e+05)+pow(logT,8)*(-7.41193275841470e+04)+pow(logT,9)*(2.84532761804913e+03);
				log_ro2=-2.43319363357935e+08+logT*(4.94943297366238e+08)+pow(logT,2)*(-4.14012682765776e+08)+pow(logT,3)*(1.74094893959945e+08)+pow(logT,4)*(-3.08472205494937e+07)+pow(logT,5)*(-3.64541790523608e+06)+pow(logT,6)*(3.04644643242255e+06)+pow(logT,7)*(-6.57068745752758e+05)+pow(logT,8)*(6.68463038097033e+04)+pow(logT,9)*(-2.73970549044038e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 && 3.40<logT && logT <3.875){
				log_ro1=-9.81443746776757e+07+logT*(2.09694525666690e+08)+pow(logT,2)*(-1.83625430263160e+08)+pow(logT,3)*(8.07836525324496e+07)+pow(logT,4)*(-1.51721406101994e+07)+pow(logT,5)*(-1.61793049608174e+06)+pow(logT,6)*(1.50586599056738e+06)+pow(logT,7)*(-3.37500976534987e+05)+pow(logT,8)*(3.54116033726517e+04)+pow(logT,9)*(-1.49191945940075e+03);
				log_ro2=4.28973981187798e+08+logT*(-8.13647468051304e+08)+pow(logT,2)*(6.33290174409065e+08)+pow(logT,3)*(-2.45804555993348e+08)+pow(logT,4)*(3.82947894887328e+07)+pow(logT,5)*(5.93978135372895e+06)+pow(logT,6)*(-3.94367912591377e+06)+pow(logT,7)*(7.81467947170254e+05)+pow(logT,8)*(-7.41193275841470e+04)+pow(logT,9)*(2.84532761804913e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 && 3.40<logT && logT <3.875){
				log_ro1=-9.31551958423463e+07+logT*(1.85524830912723e+08)+pow(logT,2)*(-1.52244641716738e+08)+pow(logT,3)*(6.28566648484347e+07)+pow(logT,4)*(-1.08700030023711e+07)+pow(logT,5)*(-1.35181905130131e+06)+pow(logT,6)*(1.08214610136006e+06)+pow(logT,7)*(-2.30592855153915e+05)+pow(logT,8)*(2.32750411963378e+04)+pow(logT,9)*(-9.48574596479680e+02);
				log_ro2=-9.81443746776757e+07+logT*(2.09694525666690e+08)+pow(logT,2)*(-1.83625430263160e+08)+pow(logT,3)*(8.07836525324496e+07)+pow(logT,4)*(-1.51721406101994e+07)+pow(logT,5)*(-1.61793049608174e+06)+pow(logT,6)*(1.50586599056738e+06)+pow(logT,7)*(-3.37500976534987e+05)+pow(logT,8)*(3.54116033726517e+04)+pow(logT,9)*(-1.49191945940075e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 && 3.40<logT && logT <3.875){
				log_ro1=-1.94592514191920e+09+logT*(3.79705581537419e+09)+pow(logT,2)*(-3.04586936946702e+09)+pow(logT,3)*(1.22299174906365e+09)+pow(logT,4)*(-2.00608602244115e+08)+pow(logT,5)*(-2.89269275678663e+07)+pow(logT,6)*(2.08054749421456e+07)+pow(logT,7)*(-4.28769593460184e+06)+pow(logT,8)*(4.21468052097039e+05)+pow(logT,9)*(-1.67611297814167e+04);
				log_ro2=-9.31551958423463e+07+logT*(1.85524830912723e+08)+pow(logT,2)*(-1.52244641716738e+08)+pow(logT,3)*(6.28566648484347e+07)+pow(logT,4)*(-1.08700030023711e+07)+pow(logT,5)*(-1.35181905130131e+06)+pow(logT,6)*(1.08214610136006e+06)+pow(logT,7)*(-2.30592855153915e+05)+pow(logT,8)*(2.32750411963378e+04)+pow(logT,9)*(-9.48574596479680e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 && 3.40<logT && logT <3.845){
				log_ro1=-2.16294210144690e+09+logT*(4.24523321960924e+09)+pow(logT,2)*(-3.42422941125624e+09)+pow(logT,3)*(1.38133889385702e+09)+pow(logT,4)*(-2.26578655584832e+08)+pow(logT,5)*(-3.38557974088208e+07)+pow(logT,6)*(2.41521582430286e+07)+pow(logT,7)*(-4.99724987736610e+06)+pow(logT,8)*(4.93786117394663e+05)+pow(logT,9)*(-1.97471704000374e+04);
				log_ro2=-1.94592514191920e+09+logT*(3.79705581537419e+09)+pow(logT,2)*(-3.04586936946702e+09)+pow(logT,3)*(1.22299174906365e+09)+pow(logT,4)*(-2.00608602244115e+08)+pow(logT,5)*(-2.89269275678663e+07)+pow(logT,6)*(2.08054749421456e+07)+pow(logT,7)*(-4.28769593460184e+06)+pow(logT,8)*(4.21468052097039e+05)+pow(logT,9)*(-1.67611297814167e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 && 3.462<logT && logT <3.778){
				log_ro1=3.93398589324670e+08+logT*(-5.75965019741409e+08)+pow(logT,2)*(3.00698596229962e+08)+pow(logT,3)*(-4.74156399923983e+07)+pow(logT,4)*(-1.20443905164906e+07)+pow(logT,5)*(4.09530135770282e+06)+pow(logT,6)*(5.70351502346786e+05)+pow(logT,7)*(-3.97135472782886e+05)+pow(logT,8)*(6.20573059530690e+04)+pow(logT,9)*(-3.32954262001412e+03);
				log_ro2=-2.16294210144690e+09+logT*(4.24523321960924e+09)+pow(logT,2)*(-3.42422941125624e+09)+pow(logT,3)*(1.38133889385702e+09)+pow(logT,4)*(-2.26578655584832e+08)+pow(logT,5)*(-3.38557974088208e+07)+pow(logT,6)*(2.41521582430286e+07)+pow(logT,7)*(-4.99724987736610e+06)+pow(logT,8)*(4.93786117394663e+05)+pow(logT,9)*(-1.97471704000374e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5 && 3.462<logT && logT <3.72){
				log_ro1=-1.92455895363827e+08+logT*(2.89529404553129e+08)+pow(logT,2)*(-1.55548874541239e+08)+pow(logT,3)*(2.54966965763610e+07)+pow(logT,4)*(6.40675399775166e+06)+pow(logT,5)*(-2.28102112469152e+06)+pow(logT,6)*(-3.12957824914168e+05)+pow(logT,7)*(2.28769832893375e+05)+pow(logT,8)*(-3.67798368794954e+04)+pow(logT,9)*(2.02591835446732e+03);
				log_ro2=3.93398589324670e+08+logT*(-5.75965019741409e+08)+pow(logT,2)*(3.00698596229962e+08)+pow(logT,3)*(-4.74156399923983e+07)+pow(logT,4)*(-1.20443905164906e+07)+pow(logT,5)*(4.09530135770282e+06)+pow(logT,6)*(5.70351502346786e+05)+pow(logT,7)*(-3.97135472782886e+05)+pow(logT,8)*(6.20573059530690e+04)+pow(logT,9)*(-3.32954262001412e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if ( 3.415<logT && logT <3.60){
				log_ro1=1.20739091115304e+09+logT*(-1.82824579445255e+09)+pow(logT,2)*(9.82376185496623e+08)+pow(logT,3)*(-1.54301301888691e+08)+pow(logT,4)*(-4.54529346404357e+07)+pow(logT,5)*(1.51068641356843e+07)+pow(logT,6)*(2.45072988500779e+06)+pow(logT,7)*(-1.66199388731518e+06)+pow(logT,8)*(2.67195687150555e+05)+pow(logT,9)*(-1.48378230594028e+04);
				log_ro2=-1.92455895363827e+08+logT*(2.89529404553129e+08)+pow(logT,2)*(-1.55548874541239e+08)+pow(logT,3)*(2.54966965763610e+07)+pow(logT,4)*(6.40675399775166e+06)+pow(logT,5)*(-2.28102112469152e+06)+pow(logT,6)*(-3.12957824914168e+05)+pow(logT,7)*(2.28769832893375e+05)+pow(logT,8)*(-3.67798368794954e+04)+pow(logT,9)*(2.02591835446732e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=5.){
				log_ro1=1.41915e-07;
				log_ro2=3.67914e-05;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5){
				log_ro1= 6.12001e-08;
				log_ro2=1.41915e-07;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0){
				log_ro1=1.87436e-08;
				log_ro2= 6.12001e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 ){
				log_ro1=7.69745e-09;
				log_ro2=1.87436e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 ){
				log_ro1=6.60063e-09;
				log_ro2=7.69745e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 ){
				log_ro1=1.21276e-09;
				log_ro2=6.60063e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 ){
				log_ro1=1.15056e-09;
				log_ro2=1.21276e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 ){
				log_ro1=1.09961e-09;
				log_ro2=1.15056e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 ){
				log_ro1=2.1252e-09;
				log_ro2=1.09961e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5){
				log_ro1=3.40883e-09;
				log_ro2=2.1252e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else {
				log_ro1=1.01267e-08;
				log_ro2=3.40883e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
}


double ro_ph_z000421151(double logg, double logT)
{ /* Calculate density in the photosphere for the z=0.000421151 [M/H]=-1.5*/
        /*density at kappa_ross=2/3 from MbRCS model atmospheres*/
		double log_ro1,log_ro2,log_ro;

		if (logg>=5.  && 3.41<logT && logT <3.59){
				log_ro1=-1.60143139830878e+09+logT*(3.12129731907923e+09)+pow(logT,2)*(-2.50089314645101e+09)+pow(logT,3)*(1.00318224456784e+09)+pow(logT,4)*(-1.64655139833208e+08)+pow(logT,5)*(-2.34344572106874e+07)+pow(logT,6)*(1.69146441333668e+07)+pow(logT,7)*(-3.48128386633615e+06)+pow(logT,8)*(3.41538971715379e+05)+pow(logT,9)*(-1.35523100879237e+04);
				log_ro2=-4.00455284365350e+09+logT*(6.00048448428426e+09)+pow(logT,2)*(-3.18093504681206e+09)+pow(logT,3)*(4.82749801717103e+08)+pow(logT,4)*(1.50588861352572e+08)+pow(logT,5)*(-4.79344106722217e+07)+pow(logT,6)*(-8.22393866364294e+06)+pow(logT,7)*(5.34288494024399e+06)+pow(logT,8)*(-8.47179960554860e+05)+pow(logT,9)*(4.65541964692112e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 && 3.40<logT && logT <3.89){
				log_ro1=-1.11994082503943e+09+logT*(2.19107328289649e+09)+pow(logT,2)*(-1.76246228257718e+09)+pow(logT,3)*(7.10093562025287e+08)+pow(logT,4)*(-1.17387276249307e+08)+pow(logT,5)*(-1.64625120124356e+07)+pow(logT,6)*(1.20304154332110e+07)+pow(logT,7)*(-2.48764112922311e+06)+pow(logT,8)*(2.45001425577568e+05)+pow(logT,9)*(-9.75674258817085e+03);
				log_ro2=-1.60143139830878e+09+logT*(3.12129731907923e+09)+pow(logT,2)*(-2.50089314645101e+09)+pow(logT,3)*(1.00318224456784e+09)+pow(logT,4)*(-1.64655139833208e+08)+pow(logT,5)*(-2.34344572106874e+07)+pow(logT,6)*(1.69146441333668e+07)+pow(logT,7)*(-3.48128386633615e+06)+pow(logT,8)*(3.41538971715379e+05)+pow(logT,9)*(-1.35523100879237e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 && 3.40<logT && logT <3.89){
				log_ro1=-7.88663238298543e+08+logT*(1.54332016919209e+09)+pow(logT,2)*(-1.24217400344902e+09)+pow(logT,3)*(5.01256738038886e+08)+pow(logT,4)*(-8.34208179592856e+07)+pow(logT,5)*(-1.13073182015308e+07)+pow(logT,6)*(8.40468072592893e+06)+pow(logT,7)*(-1.74214308520999e+06)+pow(logT,8)*(1.71751028684749e+05)+pow(logT,9)*(-6.84363081444503e+03);
				log_ro2=-1.11994082503943e+09+logT*(2.19107328289649e+09)+pow(logT,2)*(-1.76246228257718e+09)+pow(logT,3)*(7.10093562025287e+08)+pow(logT,4)*(-1.17387276249307e+08)+pow(logT,5)*(-1.64625120124356e+07)+pow(logT,6)*(1.20304154332110e+07)+pow(logT,7)*(-2.48764112922311e+06)+pow(logT,8)*(2.45001425577568e+05)+pow(logT,9)*(-9.75674258817085e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 && 3.40<logT && logT <3.90){
				log_ro1=-8.61561158744175e+08+logT*(1.68207160844439e+09)+pow(logT,2)*(-1.35066980070375e+09)+pow(logT,3)*(5.43641062953791e+08)+pow(logT,4)*(-9.01106509750280e+07)+pow(logT,5)*(-1.23170907405023e+07)+pow(logT,6)*(9.09092355224066e+06)+pow(logT,7)*(-1.87963935425103e+06)+pow(logT,8)*(1.84933466996771e+05)+pow(logT,9)*(-7.35557679035834e+03);
				log_ro2=-7.88663238298543e+08+logT*(1.54332016919209e+09)+pow(logT,2)*(-1.24217400344902e+09)+pow(logT,3)*(5.01256738038886e+08)+pow(logT,4)*(-8.34208179592856e+07)+pow(logT,5)*(-1.13073182015308e+07)+pow(logT,6)*(8.40468072592893e+06)+pow(logT,7)*(-1.74214308520999e+06)+pow(logT,8)*(1.71751028684749e+05)+pow(logT,9)*(-6.84363081444503e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 && 3.40<logT && logT <3.875){
				log_ro1=-2.22014972503264e+09+logT*(4.32111065763447e+09)+pow(logT,2)*(-3.45692396543057e+09)+pow(logT,3)*(1.38383158285749e+09)+pow(logT,4)*(-2.25904759941718e+08)+pow(logT,5)*(-3.28521207148870e+07)+pow(logT,6)*(2.34416572740276e+07)+pow(logT,7)*(-4.81460035482276e+06)+pow(logT,8)*(4.71864064642752e+05)+pow(logT,9)*(-1.87120054224842e+04);
				log_ro2=-8.61561158744175e+08+logT*(1.68207160844439e+09)+pow(logT,2)*(-1.35066980070375e+09)+pow(logT,3)*(5.43641062953791e+08)+pow(logT,4)*(-9.01106509750280e+07)+pow(logT,5)*(-1.23170907405023e+07)+pow(logT,6)*(9.09092355224066e+06)+pow(logT,7)*(-1.87963935425103e+06)+pow(logT,8)*(1.84933466996771e+05)+pow(logT,9)*(-7.35557679035834e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 && 3.40<logT && logT <3.875){
				log_ro1=-1.80546969549446e+09+logT*(3.51307903508489e+09)+pow(logT,2)*(-2.81087697965818e+09)+pow(logT,3)*(1.12671839851194e+09)+pow(logT,4)*(-1.85438778846778e+08)+pow(logT,5)*(-2.57707754587150e+07)+pow(logT,6)*(1.87751406659506e+07)+pow(logT,7)*(-3.86431891888820e+06)+pow(logT,8)*(3.78797346090972e+05)+pow(logT,9)*(-1.50147579711720e+04);
				log_ro2=-2.22014972503264e+09+logT*(4.32111065763447e+09)+pow(logT,2)*(-3.45692396543057e+09)+pow(logT,3)*(1.38383158285749e+09)+pow(logT,4)*(-2.25904759941718e+08)+pow(logT,5)*(-3.28521207148870e+07)+pow(logT,6)*(2.34416572740276e+07)+pow(logT,7)*(-4.81460035482276e+06)+pow(logT,8)*(4.71864064642752e+05)+pow(logT,9)*(-1.87120054224842e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 && 3.40<logT && logT <3.89){
				log_ro1=-2.97715360479537e+09+logT*(5.81453818213418e+09)+pow(logT,2)*(-4.66922639009506e+09)+pow(logT,3)*(1.87799021333330e+09)+pow(logT,4)*(-3.09761348862776e+08)+pow(logT,5)*(-4.35514771780760e+07)+pow(logT,6)*(3.17241793458955e+07)+pow(logT,7)*(-6.55021861677196e+06)+pow(logT,8)*(6.44325281901522e+05)+pow(logT,9)*(-2.56312983970531e+04);
				log_ro2=-1.80546969549446e+09+logT*(3.51307903508489e+09)+pow(logT,2)*(-2.81087697965818e+09)+pow(logT,3)*(1.12671839851194e+09)+pow(logT,4)*(-1.85438778846778e+08)+pow(logT,5)*(-2.57707754587150e+07)+pow(logT,6)*(1.87751406659506e+07)+pow(logT,7)*(-3.86431891888820e+06)+pow(logT,8)*(3.78797346090972e+05)+pow(logT,9)*(-1.50147579711720e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 && 3.45<logT && logT <3.845){
				log_ro1=1.98131866003905e+08+logT*(-2.89027144151047e+08)+pow(logT,2)*(1.50912398898792e+08)+pow(logT,3)*(-2.43777300880704e+07)+pow(logT,4)*(-5.63621612071392e+06)+pow(logT,5)*(2.00110208080790e+06)+pow(logT,6)*(2.50368364691496e+05)+pow(logT,7)*(-1.83613523254479e+05)+pow(logT,8)*(2.87247880174912e+04)+pow(logT,9)*(-1.53502177931919e+03);
				log_ro2=-2.97715360479537e+09+logT*(5.81453818213418e+09)+pow(logT,2)*(-4.66922639009506e+09)+pow(logT,3)*(1.87799021333330e+09)+pow(logT,4)*(-3.09761348862776e+08)+pow(logT,5)*(-4.35514771780760e+07)+pow(logT,6)*(3.17241793458955e+07)+pow(logT,7)*(-6.55021861677196e+06)+pow(logT,8)*(6.44325281901522e+05)+pow(logT,9)*(-2.56312983970531e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 && 3.45<logT && logT <3.81){
				log_ro1=2.06514055490731e+07+logT*(-2.85982627381769e+07)+pow(logT,2)*(1.39679138566065e+07)+pow(logT,3)*(-1.95027822428679e+06)+pow(logT,4)*(-5.50611092630284e+05)+pow(logT,5)*(1.58756871801144e+05)+pow(logT,6)*(2.55084096390758e+04)+pow(logT,7)*(-1.46877714148177e+04)+pow(logT,8)*(2.07727399107991e+03)+pow(logT,9)*(-1.01135179040267e+02);
				log_ro2=1.98131866003905e+08+logT*(-2.89027144151047e+08)+pow(logT,2)*(1.50912398898792e+08)+pow(logT,3)*(-2.43777300880704e+07)+pow(logT,4)*(-5.63621612071392e+06)+pow(logT,5)*(2.00110208080790e+06)+pow(logT,6)*(2.50368364691496e+05)+pow(logT,7)*(-1.83613523254479e+05)+pow(logT,8)*(2.87247880174912e+04)+pow(logT,9)*(-1.53502177931919e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5 && 3.40<logT && logT <3.74){
				log_ro1=-1.94132614638330e+08+logT*(2.91668661883867e+08)+pow(logT,2)*(-1.56723932399327e+08)+pow(logT,3)*(2.59310928389571e+07)+pow(logT,4)*(6.28550370277051e+06)+pow(logT,5)*(-2.27619763701574e+06)+pow(logT,6)*(-2.99535137869096e+05)+pow(logT,7)*(2.23540455613914e+05)+pow(logT,8)*(-3.59665767984159e+04)+pow(logT,9)*(1.97869614683105e+03);
				log_ro2=2.06514055490731e+07+logT*(-2.85982627381769e+07)+pow(logT,2)*(1.39679138566065e+07)+pow(logT,3)*(-1.95027822428679e+06)+pow(logT,4)*(-5.50611092630284e+05)+pow(logT,5)*(1.58756871801144e+05)+pow(logT,6)*(2.55084096390758e+04)+pow(logT,7)*(-1.46877714148177e+04)+pow(logT,8)*(2.07727399107991e+03)+pow(logT,9)*(-1.01135179040267e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if ( 3.40<logT && logT <3.72){
				log_ro1=-3.66370410642678e+08+logT*(5.49854041170391e+08)+pow(logT,2)*(-2.94751558017287e+08)+pow(logT,3)*(4.82312565506473e+07)+pow(logT,4)*(1.20758755505501e+07)+pow(logT,5)*(-4.29528916000656e+06)+pow(logT,6)*(-5.86069265088283e+05)+pow(logT,7)*(4.28391377600239e+05)+pow(logT,8)*(-6.87589551090302e+04)+pow(logT,9)*(3.78072199965761e+03);
				log_ro2=-1.94132614638330e+08+logT*(2.91668661883867e+08)+pow(logT,2)*(-1.56723932399327e+08)+pow(logT,3)*(2.59310928389571e+07)+pow(logT,4)*(6.28550370277051e+06)+pow(logT,5)*(-2.27619763701574e+06)+pow(logT,6)*(-2.99535137869096e+05)+pow(logT,7)*(2.23540455613914e+05)+pow(logT,8)*(-3.59665767984159e+04)+pow(logT,9)*(1.97869614683105e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}

		else if (logg>=5.){
				log_ro1=1.36792e-07;
				log_ro2=2.04871e-05;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 ){ //same as above interpolations, but takes the lowest avalaible density for T outside of models range
				log_ro1=5.91634e-08;
				log_ro2=1.36792e-07;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0){
				log_ro1=1.82371e-08;
				log_ro2=5.91634e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 ){
				log_ro1=7.21528e-09;
				log_ro2=1.82371e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 ){
				log_ro1=6.40678e-09;
				log_ro2=7.21528e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5){
				log_ro1=1.15809e-09;
				log_ro2=6.40678e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 ){
				log_ro1=7.0538e-10;
				log_ro2=1.15809e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 ){
				log_ro1=1.06465e-09;
				log_ro2=7.0538e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 ){
				log_ro1=1.00043e-09;
				log_ro2=1.06465e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5){
				log_ro1=2.06868e-09;
				log_ro2=1.00043e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else {
				log_ro1=1.4269e-09;
				log_ro2=2.06868e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}

}


double ro_ph_z0013113(double logg, double logT)
{ /* Calculate density in the photosphere for the z=0.0013113 [M/H]=-1*/
        /*density at kappa_ross=2/3 from MbRCS model atmospheres*/
		double log_ro1,log_ro2,log_ro;

		if (logg>=5. && 3.40<logT && logT <3.59){
				log_ro1=-7.22831412724337e+08+logT*(1.40115894737500e+09)+pow(logT,2)*(-1.11618775846086e+09)+pow(logT,3)*(4.44794860314638e+08)+pow(logT,4)*(-7.22100487982564e+07)+pow(logT,5)*(-1.05104072266150e+07)+pow(logT,6)*(7.44497950211812e+06)+pow(logT,7)*(-1.52112436140007e+06)+pow(logT,8)*(1.48315193496994e+05)+pow(logT,9)*(-5.85083329178044e+03);
				log_ro2=1.05203709668485e+10+logT*(-1.58369057099340e+10)+pow(logT,2)*(8.44501288178149e+09)+pow(logT,3)*(-1.30104231028877e+09)+pow(logT,4)*(-3.95611309150949e+08)+pow(logT,5)*(1.28023938597028e+08)+pow(logT,6)*(2.16882723343679e+07)+pow(logT,7)*(-1.42757166133652e+07)+pow(logT,8)*(2.27564175374246e+06)+pow(logT,9)*(-1.25597356338408e+05);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 && 3.40<logT && logT <3.89){
				log_ro1=-9.40036612370576e+08+logT*(1.82092364155756e+09)+pow(logT,2)*(-1.45007396154280e+09)+pow(logT,3)*(5.78154956103051e+08)+pow(logT,4)*(-9.43477101037182e+07)+pow(logT,5)*(-1.33241871880222e+07)+pow(logT,6)*(9.56682448670978e+06)+pow(logT,7)*(-1.95735150345627e+06)+pow(logT,8)*(1.90885134722200e+05)+pow(logT,9)*(-7.52916167066142e+03);
				log_ro2=-7.22831412724337e+08+logT*(1.40115894737500e+09)+pow(logT,2)*(-1.11618775846086e+09)+pow(logT,3)*(4.44794860314638e+08)+pow(logT,4)*(-7.22100487982564e+07)+pow(logT,5)*(-1.05104072266150e+07)+pow(logT,6)*(7.44497950211812e+06)+pow(logT,7)*(-1.52112436140007e+06)+pow(logT,8)*(1.48315193496994e+05)+pow(logT,9)*(-5.85083329178044e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 && 3.40<logT && logT <3.90){
				log_ro1=-1.45110551701117e+09+logT*(2.80789017702065e+09)+pow(logT,2)*(-2.23355561355412e+09)+pow(logT,3)*(8.89436935831511e+08)+pow(logT,4)*(-1.44856871180526e+08)+pow(logT,5)*(-2.05400658251164e+07)+pow(logT,6)*(1.46981727639087e+07)+pow(logT,7)*(-3.00342341863156e+06)+pow(logT,8)*(2.92600191173562e+05)+pow(logT,9)*(-1.15303624384956e+04);
				log_ro2=-9.40036612370576e+08+logT*(1.82092364155756e+09)+pow(logT,2)*(-1.45007396154280e+09)+pow(logT,3)*(5.78154956103051e+08)+pow(logT,4)*(-9.43477101037182e+07)+pow(logT,5)*(-1.33241871880222e+07)+pow(logT,6)*(9.56682448670978e+06)+pow(logT,7)*(-1.95735150345627e+06)+pow(logT,8)*(1.90885134722200e+05)+pow(logT,9)*(-7.52916167066142e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 && 3.40<logT && logT <3.90){
				log_ro1=-2.53020402054224e+09+logT*(4.90270347313902e+09)+pow(logT,2)*(-3.90555086278121e+09)+pow(logT,3)*(1.55781523136025e+09)+pow(logT,4)*(-2.54402320241198e+08)+pow(logT,5)*(-3.58713241310023e+07)+pow(logT,6)*(2.57902288621538e+07)+pow(logT,7)*(-5.27954112699799e+06)+pow(logT,8)*(5.15129984417861e+05)+pow(logT,9)*(-2.03288087940676e+04);
				log_ro2=-1.45110551701117e+09+logT*(2.80789017702065e+09)+pow(logT,2)*(-2.23355561355412e+09)+pow(logT,3)*(8.89436935831511e+08)+pow(logT,4)*(-1.44856871180526e+08)+pow(logT,5)*(-2.05400658251164e+07)+pow(logT,6)*(1.46981727639087e+07)+pow(logT,7)*(-3.00342341863156e+06)+pow(logT,8)*(2.92600191173562e+05)+pow(logT,9)*(-1.15303624384956e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 && 3.40<logT && logT <3.875){
				log_ro1=-3.84850089936235e+09+logT*(7.47356763804564e+09)+pow(logT,2)*(-5.96444354320511e+09)+pow(logT,3)*(2.38086267989003e+09)+pow(logT,4)*(-3.86767232766188e+08)+pow(logT,5)*(-5.68402623370048e+07)+pow(logT,6)*(4.02156902116478e+07)+pow(logT,7)*(-8.23294054221369e+06)+pow(logT,8)*(8.04674076776603e+05)+pow(logT,9)*(-3.18265298440947e+04);
				log_ro2=-2.53020402054224e+09+logT*(4.90270347313902e+09)+pow(logT,2)*(-3.90555086278121e+09)+pow(logT,3)*(1.55781523136025e+09)+pow(logT,4)*(-2.54402320241198e+08)+pow(logT,5)*(-3.58713241310023e+07)+pow(logT,6)*(2.57902288621538e+07)+pow(logT,7)*(-5.27954112699799e+06)+pow(logT,8)*(5.15129984417861e+05)+pow(logT,9)*(-2.03288087940676e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 && 3.40<logT && logT <3.875){
				log_ro1=-2.98157140979820e+09+logT*(5.80126570261354e+09)+pow(logT,2)*(-4.64130212889555e+09)+pow(logT,3)*(1.86024949818506e+09)+pow(logT,4)*(-3.06187132758583e+08)+pow(logT,5)*(-4.24783335855675e+07)+pow(logT,6)*(3.09596235799518e+07)+pow(logT,7)*(-6.37043330125614e+06)+pow(logT,8)*(6.24221206230494e+05)+pow(logT,9)*(-2.47318481981680e+04);
				log_ro2=-3.84850089936235e+09+logT*(7.47356763804564e+09)+pow(logT,2)*(-5.96444354320511e+09)+pow(logT,3)*(2.38086267989003e+09)+pow(logT,4)*(-3.86767232766188e+08)+pow(logT,5)*(-5.68402623370048e+07)+pow(logT,6)*(4.02156902116478e+07)+pow(logT,7)*(-8.23294054221369e+06)+pow(logT,8)*(8.04674076776603e+05)+pow(logT,9)*(-3.18265298440947e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 && 3.40<logT && logT <3.875){
				log_ro1=-4.53944324619369e+09+logT*(8.86248105433636e+09)+pow(logT,2)*(-7.11244328970050e+09)+pow(logT,3)*(2.85707716773577e+09)+pow(logT,4)*(-4.69030882294049e+08)+pow(logT,5)*(-6.74299827277714e+07)+pow(logT,6)*(4.85701181743210e+07)+pow(logT,7)*(-1.00104956271546e+07)+pow(logT,8)*(9.83866306416533e+05)+pow(logT,9)*(-3.91160006632733e+04);
				log_ro2=-2.98157140979820e+09+logT*(5.80126570261354e+09)+pow(logT,2)*(-4.64130212889555e+09)+pow(logT,3)*(1.86024949818506e+09)+pow(logT,4)*(-3.06187132758583e+08)+pow(logT,5)*(-4.24783335855675e+07)+pow(logT,6)*(3.09596235799518e+07)+pow(logT,7)*(-6.37043330125614e+06)+pow(logT,8)*(6.24221206230494e+05)+pow(logT,9)*(-2.47318481981680e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 && 3.40<logT && logT <3.86){
				log_ro1=-3.86093072956542e+09+logT*(7.56871669746887e+09)+pow(logT,2)*(-6.09810117090841e+09)+pow(logT,3)*(2.45831245855784e+09)+pow(logT,4)*(-4.04173973104209e+08)+pow(logT,5)*(-5.91037535661389e+07)+pow(logT,6)*(4.24807056548110e+07)+pow(logT,7)*(-8.78326377739373e+06)+pow(logT,8)*(8.66447908490990e+05)+pow(logT,9)*(-3.45803860623623e+04);
				log_ro2=-4.53944324619369e+09+logT*(8.86248105433636e+09)+pow(logT,2)*(-7.11244328970050e+09)+pow(logT,3)*(2.85707716773577e+09)+pow(logT,4)*(-4.69030882294049e+08)+pow(logT,5)*(-6.74299827277714e+07)+pow(logT,6)*(4.85701181743210e+07)+pow(logT,7)*(-1.00104956271546e+07)+pow(logT,8)*(9.83866306416533e+05)+pow(logT,9)*(-3.91160006632733e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 && 3.40<logT && logT <3.81){
				log_ro1=-2.71236421487407e+07+logT*(4.25923054618106e+07)+pow(logT,2)*(-2.41041135832965e+07)+pow(logT,3)*(4.45149232567611e+06)+pow(logT,4)*(8.67105460908200e+05)+pow(logT,5)*(-3.73425881010508e+05)+pow(logT,6)*(-3.74848981103662e+04)+pow(logT,7)*(3.46016647414107e+04)+pow(logT,8)*(-5.83057266726148e+03)+pow(logT,9)*(3.30622674578749e+02);
				log_ro2=-3.86093072956542e+09+logT*(7.56871669746887e+09)+pow(logT,2)*(-6.09810117090841e+09)+pow(logT,3)*(2.45831245855784e+09)+pow(logT,4)*(-4.04173973104209e+08)+pow(logT,5)*(-5.91037535661389e+07)+pow(logT,6)*(4.24807056548110e+07)+pow(logT,7)*(-8.78326377739373e+06)+pow(logT,8)*(8.66447908490990e+05)+pow(logT,9)*(-3.45803860623623e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5 && 3.40<logT && logT <3.74){
				log_ro1=-2.65475349133105e+08+logT*(3.97974150192618e+08)+pow(logT,2)*(-2.13305480295165e+08)+pow(logT,3)*(3.51296075860888e+07)+pow(logT,4)*(8.56634743627171e+06)+pow(logT,5)*(-3.08237538240698e+06)+pow(logT,6)*(-4.08599194232160e+05)+pow(logT,7)*(3.02730288417383e+05)+pow(logT,8)*(-4.85837785012797e+04)+pow(logT,9)*(2.66729390061113e+03);
				log_ro2=-2.71236421487407e+07+logT*(4.25923054618106e+07)+pow(logT,2)*(-2.41041135832965e+07)+pow(logT,3)*(4.45149232567611e+06)+pow(logT,4)*(8.67105460908200e+05)+pow(logT,5)*(-3.73425881010508e+05)+pow(logT,6)*(-3.74848981103662e+04)+pow(logT,7)*(3.46016647414107e+04)+pow(logT,8)*(-5.83057266726148e+03)+pow(logT,9)*(3.30622674578749e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0 && 3.40<logT && logT <3.70){
				log_ro1=-6.20657891466988e+08+logT*(9.31386867004959e+08)+pow(logT,2)*(-4.98898608710686e+08)+pow(logT,3)*(8.12386144117074e+07)+pow(logT,4)*(2.06519260064808e+07)+pow(logT,5)*(-7.28402120226985e+06)+pow(logT,6)*(-1.01393647837290e+06)+pow(logT,7)*(7.33290214029206e+05)+pow(logT,8)*(-1.17581595289450e+05)+pow(logT,9)*(6.46532833404009e+03);
				log_ro2=-2.65475349133105e+08+logT*(3.97974150192618e+08)+pow(logT,2)*(-2.13305480295165e+08)+pow(logT,3)*(3.51296075860888e+07)+pow(logT,4)*(8.56634743627171e+06)+pow(logT,5)*(-3.08237538240698e+06)+pow(logT,6)*(-4.08599194232160e+05)+pow(logT,7)*(3.02730288417383e+05)+pow(logT,8)*(-4.85837785012797e+04)+pow(logT,9)*(2.66729390061113e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else if ( 3.40<logT && logT <3.53){
				log_ro1=-3.91897630561175e+00+logT*(-8.96223063299385e-01)+pow(logT,2)*(-1.95610077512887e-01)+pow(logT,3)*(-3.95411115542704e-02)+pow(logT,4)*(-6.87845905931033e-03)+pow(logT,5)*(-7.71193229546380e-04)+pow(logT,6)*(1.02189223838343e-04)+pow(logT,7)*(1.16278040448750e-04)+pow(logT,8)*(5.67231758144790e-05)+pow(logT,9)*(2.25458896730292e-05);
				log_ro2=-6.20657891466988e+08+logT*(9.31386867004959e+08)+pow(logT,2)*(-4.98898608710686e+08)+pow(logT,3)*(8.12386144117074e+07)+pow(logT,4)*(2.06519260064808e+07)+pow(logT,5)*(-7.28402120226985e+06)+pow(logT,6)*(-1.01393647837290e+06)+pow(logT,7)*(7.33290214029206e+05)+pow(logT,8)*(-1.17581595289450e+05)+pow(logT,9)*(6.46532833404009e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}

		else if (logg>=5.){
				log_ro1=1.29102e-07;
				log_ro2=1.08195e-05;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 ){ //same as above interpolations, but takes the lowest avalaible density for T outside of models range
				log_ro1=4.14416e-08;
				log_ro2=1.29102e-07;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 ){
				log_ro1=1.75382e-08;
				log_ro2=4.14416e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 ){
				log_ro1=6.80528e-09;
				log_ro2=1.75382e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0){
				log_ro1=6.11551e-09;
				log_ro2=6.80528e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5){
				log_ro1=1.09481e-09;
				log_ro2=6.11551e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0){
				log_ro1=1.06461e-09;
				log_ro2=1.09481e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5){
				log_ro1=6.70297e-10;
				log_ro2=1.06461e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0){
				log_ro1=9.45123e-10;
				log_ro2=6.70297e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5){
				log_ro1=1.94203e-09;
				log_ro2=9.45123e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0){
				log_ro1=1.76923e-09;
				log_ro2=1.94203e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else {
				log_ro1=1.24326e-09;
				log_ro2=1.76923e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}
}


double ro_ph_z0023094(double logg, double logT)
{ /* Calculate density in the photosphere for the z=0.0023094 [M/H]=-0.75*/
        /*density at kappa_ross=2/3 from MbRCS model atmospheres*/
		double log_ro1,log_ro2,log_ro;

		if (logg>=5. && 3.39<logT && logT <3.59){
				log_ro1=-3.03326406497823e+08+logT*(5.80330823831566e+08)+pow(logT,2)*(-4.55830276730867e+08)+pow(logT,3)*(1.78686938739821e+08)+pow(logT,4)*(-2.82033318537911e+07)+pow(logT,5)*(-4.34379921991862e+06)+pow(logT,6)*(2.93546372373704e+06)+pow(logT,7)*(-5.88374685443020e+05)+pow(logT,8)*(5.64255896805625e+04)+pow(logT,9)*(-2.19043427353118e+03);
				log_ro2=1.15615691974722e+10+logT*(-1.74041456003089e+10)+pow(logT,2)*(9.28065772044609e+09)+pow(logT,3)*(-1.42975608522836e+09)+pow(logT,4)*(-4.34759049132503e+08)+pow(logT,5)*(1.40689375870728e+08)+pow(logT,6)*(2.38344706666148e+07)+pow(logT,7)*(-1.56880317977623e+07)+pow(logT,8)*(2.50075056817451e+06)+pow(logT,9)*(-1.38020502633656e+05);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 && 3.39<logT && logT <3.89){
				log_ro1=-8.90293746730805e+08+logT*(1.71758520060073e+09)+pow(logT,2)*(-1.36192110028288e+09)+pow(logT,3)*(5.40352616835075e+08)+pow(logT,4)*(-8.74636078808950e+07)+pow(logT,5)*(-1.25603849751991e+07)+pow(logT,6)*(8.89270714239832e+06)+pow(logT,7)*(-1.80955460706363e+06)+pow(logT,8)*(1.75664225870155e+05)+pow(logT,9)*(-6.89878687022667e+03);
				log_ro2=-3.03326406497823e+08+logT*(5.80330823831566e+08)+pow(logT,2)*(-4.55830276730867e+08)+pow(logT,3)*(1.78686938739821e+08)+pow(logT,4)*(-2.82033318537911e+07)+pow(logT,5)*(-4.34379921991862e+06)+pow(logT,6)*(2.93546372373704e+06)+pow(logT,7)*(-5.88374685443020e+05)+pow(logT,8)*(5.64255896805625e+04)+pow(logT,9)*(-2.19043427353118e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 && 3.39<logT && logT <3.90){
				log_ro1=-1.70284642892472e+09+logT*(3.29022203826585e+09)+pow(logT,2)*(-2.61318852831014e+09)+pow(logT,3)*(1.03877126965855e+09)+pow(logT,4)*(-1.68679214441409e+08)+pow(logT,5)*(-2.40638724111136e+07)+pow(logT,6)*(1.71320411720564e+07)+pow(logT,7)*(-3.49382292867403e+06)+pow(logT,8)*(3.39806279380297e+05)+pow(logT,9)*(-1.33693405072193e+04);
				log_ro2=-8.90293746730805e+08+logT*(1.71758520060073e+09)+pow(logT,2)*(-1.36192110028288e+09)+pow(logT,3)*(5.40352616835075e+08)+pow(logT,4)*(-8.74636078808950e+07)+pow(logT,5)*(-1.25603849751991e+07)+pow(logT,6)*(8.89270714239832e+06)+pow(logT,7)*(-1.80955460706363e+06)+pow(logT,8)*(1.75664225870155e+05)+pow(logT,9)*(-6.89878687022667e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 && 3.39<logT && logT <3.90){
				log_ro1=-3.83519896821103e+09+logT*(7.41963895888342e+09)+pow(logT,2)*(-5.90013496870418e+09)+pow(logT,3)*(2.34796108365813e+09)+pow(logT,4)*(-3.81379177244428e+08)+pow(logT,5)*(-5.47818916905151e+07)+pow(logT,6)*(3.89556150589090e+07)+pow(logT,7)*(-7.95322998483950e+06)+pow(logT,8)*(7.74585314381339e+05)+pow(logT,9)*(-3.05201822221351e+04);
				log_ro2=-1.70284642892472e+09+logT*(3.29022203826585e+09)+pow(logT,2)*(-2.61318852831014e+09)+pow(logT,3)*(1.03877126965855e+09)+pow(logT,4)*(-1.68679214441409e+08)+pow(logT,5)*(-2.40638724111136e+07)+pow(logT,6)*(1.71320411720564e+07)+pow(logT,7)*(-3.49382292867403e+06)+pow(logT,8)*(3.39806279380297e+05)+pow(logT,9)*(-1.33693405072193e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 && 3.39<logT && logT <3.875){
				log_ro1=-3.81322983155095e+09+logT*(7.40003834154415e+09)+pow(logT,2)*(-5.90142372741789e+09)+pow(logT,3)*(2.35368909568797e+09)+pow(logT,4)*(-3.81794815014553e+08)+pow(logT,5)*(-5.62767668913336e+07)+pow(logT,6)*(3.97179428473836e+07)+pow(logT,7)*(-8.12293510163050e+06)+pow(logT,8)*(7.93242556192374e+05)+pow(logT,9)*(-3.13485468472434e+04);
				log_ro2=-3.83519896821103e+09+logT*(7.41963895888342e+09)+pow(logT,2)*(-5.90013496870418e+09)+pow(logT,3)*(2.34796108365813e+09)+pow(logT,4)*(-3.81379177244428e+08)+pow(logT,5)*(-5.47818916905151e+07)+pow(logT,6)*(3.89556150589090e+07)+pow(logT,7)*(-7.95322998483950e+06)+pow(logT,8)*(7.74585314381339e+05)+pow(logT,9)*(-3.05201822221351e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 && 3.39<logT && logT <3.875){
				log_ro1=-3.44951416035320e+09+logT*(6.71074500410413e+09)+pow(logT,2)*(-5.36793539443262e+09)+pow(logT,3)*(2.15089703352343e+09)+pow(logT,4)*(-3.53768162088854e+08)+pow(logT,5)*(-4.92176667896804e+07)+pow(logT,6)*(3.58135353656206e+07)+pow(logT,7)*(-7.36650132764398e+06)+pow(logT,8)*(7.21644528691532e+05)+pow(logT,9)*(-2.85856184569697e+04);
				log_ro2=-3.81322983155095e+09+logT*(7.40003834154415e+09)+pow(logT,2)*(-5.90142372741789e+09)+pow(logT,3)*(2.35368909568797e+09)+pow(logT,4)*(-3.81794815014553e+08)+pow(logT,5)*(-5.62767668913336e+07)+pow(logT,6)*(3.97179428473836e+07)+pow(logT,7)*(-8.12293510163050e+06)+pow(logT,8)*(7.93242556192374e+05)+pow(logT,9)*(-3.13485468472434e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 && 3.39<logT && logT <3.875){
				log_ro1=-4.34089082433413e+09+logT*(8.47682997990928e+09)+pow(logT,2)*(-6.80488465660159e+09)+pow(logT,3)*(2.73477017158491e+09)+pow(logT,4)*(-4.49606015072714e+08)+pow(logT,5)*(-6.42215348268630e+07)+pow(logT,6)*(4.64113997261114e+07)+pow(logT,7)*(-9.57067469346154e+06)+pow(logT,8)*(9.40864433810968e+05)+pow(logT,9)*(-3.74115933450074e+04);
				log_ro2=-3.44951416035320e+09+logT*(6.71074500410413e+09)+pow(logT,2)*(-5.36793539443262e+09)+pow(logT,3)*(2.15089703352343e+09)+pow(logT,4)*(-3.53768162088854e+08)+pow(logT,5)*(-4.92176667896804e+07)+pow(logT,6)*(3.58135353656206e+07)+pow(logT,7)*(-7.36650132764398e+06)+pow(logT,8)*(7.21644528691532e+05)+pow(logT,9)*(-2.85856184569697e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 && 3.39<logT && logT <3.845){
				log_ro1=-4.29702629776543e+09+logT*(8.41950054405086e+09)+pow(logT,2)*(-6.77834319974799e+09)+pow(logT,3)*(2.72820835865869e+09)+pow(logT,4)*(-4.45783156038439e+08)+pow(logT,5)*(-6.70949253115187e+07)+pow(logT,6)*(4.75532379646051e+07)+pow(logT,7)*(-9.81116453394117e+06)+pow(logT,8)*(9.66994597570047e+05)+pow(logT,9)*(-3.85743299729791e+04);
				log_ro2=-4.34089082433413e+09+logT*(8.47682997990928e+09)+pow(logT,2)*(-6.80488465660159e+09)+pow(logT,3)*(2.73477017158491e+09)+pow(logT,4)*(-4.49606015072714e+08)+pow(logT,5)*(-6.42215348268630e+07)+pow(logT,6)*(4.64113997261114e+07)+pow(logT,7)*(-9.57067469346154e+06)+pow(logT,8)*(9.40864433810968e+05)+pow(logT,9)*(-3.74115933450074e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 && 3.39<logT && logT <3.813){
				log_ro1=-3.80791745238996e+07+logT*(5.88979890501125e+07)+pow(logT,2)*(-3.28128226392863e+07)+pow(logT,3)*(5.91366860675555e+06)+pow(logT,4)*(1.19051071227975e+06)+pow(logT,5)*(-4.94601300355215e+05)+pow(logT,6)*(-5.18150376866676e+04)+pow(logT,7)*(4.57884246843095e+04)+pow(logT,8)*(-7.62173905440892e+03)+pow(logT,9)*(4.28214100181299e+02);
				log_ro2=-4.29702629776543e+09+logT*(8.41950054405086e+09)+pow(logT,2)*(-6.77834319974799e+09)+pow(logT,3)*(2.72820835865869e+09)+pow(logT,4)*(-4.45783156038439e+08)+pow(logT,5)*(-6.70949253115187e+07)+pow(logT,6)*(4.75532379646051e+07)+pow(logT,7)*(-9.81116453394117e+06)+pow(logT,8)*(9.66994597570047e+05)+pow(logT,9)*(-3.85743299729791e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5 && 3.39<logT && logT <3.72){
				log_ro1=-3.25292091988802e+08+logT*(4.87296291657233e+08)+pow(logT,2)*(-2.60624982340208e+08)+pow(logT,3)*(4.24520307612982e+07)+pow(logT,4)*(1.06960296350542e+07)+pow(logT,5)*(-3.77903713082907e+06)+pow(logT,6)*(-5.20406175272028e+05)+pow(logT,7)*(3.77284823065130e+05)+pow(logT,8)*(-6.03876172160358e+04)+pow(logT,9)*(3.31283415000969e+03);
				log_ro2=-3.80791745238996e+07+logT*(5.88979890501125e+07)+pow(logT,2)*(-3.28128226392863e+07)+pow(logT,3)*(5.91366860675555e+06)+pow(logT,4)*(1.19051071227975e+06)+pow(logT,5)*(-4.94601300355215e+05)+pow(logT,6)*(-5.18150376866676e+04)+pow(logT,7)*(4.57884246843095e+04)+pow(logT,8)*(-7.62173905440892e+03)+pow(logT,9)*(4.28214100181299e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0 && 3.39<logT && logT <3.70){
				log_ro1=-3.86004804011404e+08+logT*(5.78955755470543e+08)+pow(logT,2)*(-3.09693679219801e+08)+pow(logT,3)*(5.01019183702376e+07)+pow(logT,4)*(1.29679975770378e+07)+pow(logT,5)*(-4.52750490189598e+06)+pow(logT,6)*(-6.42206067201940e+05)+pow(logT,7)*(4.59210505372117e+05)+pow(logT,8)*(-7.34921888003251e+04)+pow(logT,9)*(4.03698496831106e+03);
				log_ro2=-3.25292091988802e+08+logT*(4.87296291657233e+08)+pow(logT,2)*(-2.60624982340208e+08)+pow(logT,3)*(4.24520307612982e+07)+pow(logT,4)*(1.06960296350542e+07)+pow(logT,5)*(-3.77903713082907e+06)+pow(logT,6)*(-5.20406175272028e+05)+pow(logT,7)*(3.77284823065130e+05)+pow(logT,8)*(-6.03876172160358e+04)+pow(logT,9)*(3.31283415000969e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else if (3.39<logT && logT <3.53){
				log_ro1=-3.96183650190111e+00+logT*(-9.05898401311133e-01)+pow(logT,2)*(-1.97679219826871e-01)+pow(logT,3)*(-3.99443141443818e-02)+pow(logT,4)*(-6.94285145194954e-03)+pow(logT,5)*(-7.75863610719533e-04)+pow(logT,6)*(1.04562219123693e-04)+pow(logT,7)*(1.17956869571073e-04)+pow(logT,8)*(5.74728983506185e-05)+pow(logT,9)*(2.28329313611662e-05);
				log_ro2=-3.86004804011404e+08+logT*(5.78955755470543e+08)+pow(logT,2)*(-3.09693679219801e+08)+pow(logT,3)*(5.01019183702376e+07)+pow(logT,4)*(1.29679975770378e+07)+pow(logT,5)*(-4.52750490189598e+06)+pow(logT,6)*(-6.42206067201940e+05)+pow(logT,7)*(4.59210505372117e+05)+pow(logT,8)*(-7.34921888003251e+04)+pow(logT,9)*(4.03698496831106e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=5.){
				log_ro1=1.24517e-07;
				log_ro2= 8.68569e-06;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 ){ //same as above interpolations, but takes the lowest avalaible density for T outside of models range
				log_ro1=4.01877e-08;
				log_ro2=1.24517e-07;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 ){
				log_ro1=1.71032e-08;
				log_ro2=4.01877e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 ){
				log_ro1=6.57562e-09;
				log_ro2=1.71032e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0){
				log_ro1=5.93883e-09;
				log_ro2=6.57562e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5){
				log_ro1=1.06434e-09;
				log_ro2=5.93883e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0){
				log_ro1=1.03395e-09;
				log_ro2=1.06434e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5){
				log_ro1=9.81293e-10;
				log_ro2=1.03395e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0){
				log_ro1=9.13576e-10;
				log_ro2=9.81293e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5){
				log_ro1=2.57425e-09;
				log_ro2=9.13576e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0){
				log_ro1=1.61692e-09;
				log_ro2=2.57425e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else {
				log_ro1=1.02685e-09;
				log_ro2=1.61692e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}

}


double ro_ph_z004050(double logg, double logT)
{ /* Calculate density in the photosphere for the z=0.004050 [M/H]=-0.5*/
        /*density at kappa_ross=2/3 from MbRCS model atmospheres*/
		double log_ro1,log_ro2,log_ro;

		if (logg>=5. && 3.39<logT && logT <3.59){
				log_ro1=1.18929264932437e+08+logT*(-2.45713119697450e+08)+pow(logT,2)*(2.08580940856639e+08)+pow(logT,3)*(-8.89886349968122e+07)+pow(logT,4)*(1.60444835007922e+07)+pow(logT,5)*(1.86193520094000e+06)+pow(logT,6)*(-1.59936432866088e+06)+pow(logT,7)*(3.49347588693230e+05)+pow(logT,8)*(-3.59318845649456e+04)+pow(logT,9)*(1.48776534452295e+03);
				log_ro2=2.54119141561598e+08+logT*(-3.85178440559439e+08)+pow(logT,2)*(2.07078688514593e+08)+pow(logT,3)*(-3.24335447468611e+07)+pow(logT,4)*(-9.66845697697444e+06)+pow(logT,5)*(3.19467745733725e+06)+pow(logT,6)*(5.29221409131947e+05)+pow(logT,7)*(-3.55959498870293e+05)+pow(logT,8)*(5.72451841674433e+04)+pow(logT,9)*(-3.18304005463089e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 && 3.39<logT && logT <3.89){
				log_ro1=-8.31305188135024e+08+logT*(1.59649609569858e+09)+pow(logT,2)*(-1.25976522350178e+09)+pow(logT,3)*(4.97026873520763e+08)+pow(logT,4)*(-7.96941545438536e+07)+pow(logT,5)*(-1.16668971580259e+07)+pow(logT,6)*(8.12827895372840e+06)+pow(logT,7)*(-1.64349986727521e+06)+pow(logT,8)*(1.58680688975928e+05)+pow(logT,9)*(-6.19955322505916e+03);
				log_ro2=1.18929264932437e+08+logT*(-2.45713119697450e+08)+pow(logT,2)*(2.08580940856639e+08)+pow(logT,3)*(-8.89886349968122e+07)+pow(logT,4)*(1.60444835007922e+07)+pow(logT,5)*(1.86193520094000e+06)+pow(logT,6)*(-1.59936432866088e+06)+pow(logT,7)*(3.49347588693230e+05)+pow(logT,8)*(-3.59318845649456e+04)+pow(logT,9)*(1.48776534452295e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 && 3.39<logT && logT <3.90){
				log_ro1=-1.85898831979880e+09+logT*(3.58721504250685e+09)+pow(logT,2)*(-2.84508753805849e+09)+pow(logT,3)*(1.12913410047753e+09)+pow(logT,4)*(-1.82858151060000e+08)+pow(logT,5)*(-2.62318422958657e+07)+pow(logT,6)*(1.85887218313920e+07)+pow(logT,7)*(-3.78398187059921e+06)+pow(logT,8)*(3.67457382228166e+05)+pow(logT,9)*(-1.44359327545677e+04);
				log_ro2=-8.31305188135024e+08+logT*(1.59649609569858e+09)+pow(logT,2)*(-1.25976522350178e+09)+pow(logT,3)*(4.97026873520763e+08)+pow(logT,4)*(-7.96941545438536e+07)+pow(logT,5)*(-1.16668971580259e+07)+pow(logT,6)*(8.12827895372840e+06)+pow(logT,7)*(-1.64349986727521e+06)+pow(logT,8)*(1.58680688975928e+05)+pow(logT,9)*(-6.19955322505916e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 && 3.39<logT && logT <3.90){
				log_ro1=-4.19949863659828e+09+logT*(8.12310691238750e+09)+pow(logT,2)*(-6.45836480012668e+09)+pow(logT,3)*(2.56954572618072e+09)+pow(logT,4)*(-4.17210646824113e+08)+pow(logT,5)*(-5.99771231198673e+07)+pow(logT,6)*(4.26209155364340e+07)+pow(logT,7)*(-8.69912892536417e+06)+pow(logT,8)*(8.47023579829742e+05)+pow(logT,9)*(-3.33663867793518e+04);
				log_ro2=-1.85898831979880e+09+logT*(3.58721504250685e+09)+pow(logT,2)*(-2.84508753805849e+09)+pow(logT,3)*(1.12913410047753e+09)+pow(logT,4)*(-1.82858151060000e+08)+pow(logT,5)*(-2.62318422958657e+07)+pow(logT,6)*(1.85887218313920e+07)+pow(logT,7)*(-3.78398187059921e+06)+pow(logT,8)*(3.67457382228166e+05)+pow(logT,9)*(-1.44359327545677e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 && 3.39<logT && logT <3.875){
				log_ro1=-3.96582588933917e+09+logT*(7.69574148108211e+09)+pow(logT,2)*(-6.13681881670762e+09)+pow(logT,3)*(2.44735094870772e+09)+pow(logT,4)*(-3.96920033534321e+08)+pow(logT,5)*(-5.85274762028468e+07)+pow(logT,6)*(4.12936672664154e+07)+pow(logT,7)*(-8.44407209451226e+06)+pow(logT,8)*(8.24502671136206e+05)+pow(logT,9)*(-3.25799044640849e+04);
				log_ro2=-4.19949863659828e+09+logT*(8.12310691238750e+09)+pow(logT,2)*(-6.45836480012668e+09)+pow(logT,3)*(2.56954572618072e+09)+pow(logT,4)*(-4.17210646824113e+08)+pow(logT,5)*(-5.99771231198673e+07)+pow(logT,6)*(4.26209155364340e+07)+pow(logT,7)*(-8.69912892536417e+06)+pow(logT,8)*(8.47023579829742e+05)+pow(logT,9)*(-3.33663867793518e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 && 3.39<logT && logT <3.875){
				log_ro1=-3.56183843010269e+09+logT*(6.92884562741695e+09)+pow(logT,2)*(-5.54204420762321e+09)+pow(logT,3)*(2.22054245519742e+09)+pow(logT,4)*(-3.65239403117845e+08)+pow(logT,5)*(-5.07718726601483e+07)+pow(logT,6)*(3.69531664639346e+07)+pow(logT,7)*(-7.60032757403412e+06)+pow(logT,8)*(7.44463908931482e+05)+pow(logT,9)*(-2.94855450087247e+04);
				log_ro2=-3.96582588933917e+09+logT*(7.69574148108211e+09)+pow(logT,2)*(-6.13681881670762e+09)+pow(logT,3)*(2.44735094870772e+09)+pow(logT,4)*(-3.96920033534321e+08)+pow(logT,5)*(-5.85274762028468e+07)+pow(logT,6)*(4.12936672664154e+07)+pow(logT,7)*(-8.44407209451226e+06)+pow(logT,8)*(8.24502671136206e+05)+pow(logT,9)*(-3.25799044640849e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 && 3.39<logT && logT <3.875){
				log_ro1=-4.07722907324678e+09+logT*(7.96353430956592e+09)+pow(logT,2)*(-6.39402539496387e+09)+pow(logT,3)*(2.57013356568292e+09)+pow(logT,4)*(-4.22652085317985e+08)+pow(logT,5)*(-6.03408339343433e+07)+pow(logT,6)*(4.36249339909009e+07)+pow(logT,7)*(-8.99726538844196e+06)+pow(logT,8)*(8.84574130484556e+05)+pow(logT,9)*(-3.51756516102802e+04);
				log_ro2=-3.56183843010269e+09+logT*(6.92884562741695e+09)+pow(logT,2)*(-5.54204420762321e+09)+pow(logT,3)*(2.22054245519742e+09)+pow(logT,4)*(-3.65239403117845e+08)+pow(logT,5)*(-5.07718726601483e+07)+pow(logT,6)*(3.69531664639346e+07)+pow(logT,7)*(-7.60032757403412e+06)+pow(logT,8)*(7.44463908931482e+05)+pow(logT,9)*(-2.94855450087247e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 && 3.39<logT && logT <3.845){
				log_ro1=-3.93323226599362e+09+logT*(7.70478435051065e+09)+pow(logT,2)*(-6.20057281652383e+09)+pow(logT,3)*(2.49383044928311e+09)+pow(logT,4)*(-4.06407441715303e+08)+pow(logT,5)*(-6.18808001045684e+07)+pow(logT,6)*(4.36056434676535e+07)+pow(logT,7)*(-8.98754250771501e+06)+pow(logT,8)*(8.85348312297001e+05)+pow(logT,9)*(-3.53036906618351e+04);
				log_ro2=-4.07722907324678e+09+logT*(7.96353430956592e+09)+pow(logT,2)*(-6.39402539496387e+09)+pow(logT,3)*(2.57013356568292e+09)+pow(logT,4)*(-4.22652085317985e+08)+pow(logT,5)*(-6.03408339343433e+07)+pow(logT,6)*(4.36249339909009e+07)+pow(logT,7)*(-8.99726538844196e+06)+pow(logT,8)*(8.84574130484556e+05)+pow(logT,9)*(-3.51756516102802e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 && 3.39<logT && logT <3.813){
				log_ro1=-3.60400791971980e+07+logT*(5.59673357346569e+07)+pow(logT,2)*(-3.13076931720821e+07)+pow(logT,3)*(5.67621225184161e+06)+pow(logT,4)*(1.13518366158384e+06)+pow(logT,5)*(-4.75580021666528e+05)+pow(logT,6)*(-4.93773014462077e+04)+pow(logT,7)*(4.40718461947579e+04)+pow(logT,8)*(-7.35868944889098e+03)+pow(logT,9)*(4.14423296803469e+02);
				log_ro2=-3.93323226599362e+09+logT*(7.70478435051065e+09)+pow(logT,2)*(-6.20057281652383e+09)+pow(logT,3)*(2.49383044928311e+09)+pow(logT,4)*(-4.06407441715303e+08)+pow(logT,5)*(-6.18808001045684e+07)+pow(logT,6)*(4.36056434676535e+07)+pow(logT,7)*(-8.98754250771501e+06)+pow(logT,8)*(8.85348312297001e+05)+pow(logT,9)*(-3.53036906618351e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5 && 3.39<logT && logT <3.72){
				log_ro1=-2.34496583086725e+08+logT*(3.52024368591689e+08)+pow(logT,2)*(-1.88720760107833e+08)+pow(logT,3)*(3.08706460645591e+07)+pow(logT,4)*(7.73677347495752e+06)+pow(logT,5)*(-2.74888584340017e+06)+pow(logT,6)*(-3.76200355038895e+05)+pow(logT,7)*(2.74403114589544e+05)+pow(logT,8)*(-4.40133641757468e+04)+pow(logT,9)*(2.41856081441414e+03);
				log_ro2=-3.60400791971980e+07+logT*(5.59673357346569e+07)+pow(logT,2)*(-3.13076931720821e+07)+pow(logT,3)*(5.67621225184161e+06)+pow(logT,4)*(1.13518366158384e+06)+pow(logT,5)*(-4.75580021666528e+05)+pow(logT,6)*(-4.93773014462077e+04)+pow(logT,7)*(4.40718461947579e+04)+pow(logT,8)*(-7.35868944889098e+03)+pow(logT,9)*(4.14423296803469e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0 && 3.39<logT && logT <3.72){
				log_ro1=-2.18607247341939e+08+logT*(3.29372106720625e+08)+pow(logT,2)*(-1.77329336539126e+08)+pow(logT,3)*(2.92427134681881e+07)+pow(logT,4)*(7.25075686638147e+06)+pow(logT,5)*(-2.60564046049661e+06)+pow(logT,6)*(-3.51501752448444e+05)+pow(logT,7)*(2.59814201562458e+05)+pow(logT,8)*(-4.18607407249876e+04)+pow(logT,9)*(2.30857006730773e+03);
				log_ro2=-2.34496583086725e+08+logT*(3.52024368591689e+08)+pow(logT,2)*(-1.88720760107833e+08)+pow(logT,3)*(3.08706460645591e+07)+pow(logT,4)*(7.73677347495752e+06)+pow(logT,5)*(-2.74888584340017e+06)+pow(logT,6)*(-3.76200355038895e+05)+pow(logT,7)*(2.74403114589544e+05)+pow(logT,8)*(-4.40133641757468e+04)+pow(logT,9)*(2.41856081441414e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else if (3.39<logT && logT <3.556){
				log_ro1=-1.73887026140743e+05+logT*(6.28965587130427e+04)+pow(logT,2)*(1.27501433893869e+04)+pow(logT,3)*(-9.38544909299892e+02)+pow(logT,4)*(-1.15376918510545e+03)+pow(logT,5)*(-3.04181117495871e+02)+pow(logT,6)*(-8.61481141863517e+00)+pow(logT,7)*(2.45504004068473e+01)+pow(logT,8)*(8.36069993067116e+00)+pow(logT,9)*(-2.02626322305965e+00);
				log_ro2=-2.18607247341939e+08+logT*(3.29372106720625e+08)+pow(logT,2)*(-1.77329336539126e+08)+pow(logT,3)*(2.92427134681881e+07)+pow(logT,4)*(7.25075686638147e+06)+pow(logT,5)*(-2.60564046049661e+06)+pow(logT,6)*(-3.51501752448444e+05)+pow(logT,7)*(2.59814201562458e+05)+pow(logT,8)*(-4.18607407249876e+04)+pow(logT,9)*(2.30857006730773e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=5.){
				log_ro1=1.19231e-07;
				log_ro2=6.96534e-06;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 ){ //same as above interpolations, but takes the lowest avalaible density for T outside of models range
				log_ro1=3.87098e-08;
				log_ro2=1.19231e-07;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 ){
				log_ro1=1.6574e-08;
				log_ro2=3.87098e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 ){
				log_ro1=6.30751e-09;
				log_ro2=1.6574e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0){
				log_ro1=5.73135e-09;
				log_ro2=6.30751e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5){
				log_ro1=1.02891e-09;
				log_ro2=5.73135e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0){
				log_ro1=9.97546e-10;
				log_ro2=1.02891e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5){
				log_ro1=9.44497e-10;
				log_ro2=9.97546e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0){
				log_ro1=8.77202e-10;
				log_ro2=9.44497e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5){
				log_ro1=2.34221e-09;
				log_ro2=8.77202e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0){
				log_ro1=1.12094e-09;
				log_ro2=2.34221e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else {
				log_ro1=8.62083e-10;
				log_ro2=1.12094e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}

}


double ro_ph_z0070598(double logg, double logT)
{ /* Calculate density in the photosphere for the z=0.0070598 [M/H]=-0.25*/
        /*density at kappa_ross=2/3 from MbRCS model atmospheres*/
		double log_ro1,log_ro2,log_ro;

		if (logg>=5.	 && 3.39<logT && logT <3.59){
				log_ro1=5.12327999813832e+08+logT*(-1.01505440812654e+09)+pow(logT,2)*(8.27170971759619e+08)+pow(logT,3)*(-3.38105658649033e+08)+pow(logT,4)*(5.71975963297892e+07)+pow(logT,5)*(7.64153606342723e+06)+pow(logT,6)*(-5.81792652386035e+06)+pow(logT,7)*(1.22128553946297e+06)+pow(logT,8)*(-1.21778047161494e+05)+pow(logT,9)*(4.90543178391408e+03);
				log_ro2=-2.74364460900072e+09+logT*(4.15172876184446e+09)+pow(logT,2)*(-2.22737720783382e+09)+pow(logT,3)*(3.47308052288035e+08)+pow(logT,4)*(1.04105219395761e+08)+pow(logT,5)*(-3.41965549888098e+07)+pow(logT,6)*(-5.70188756222936e+06)+pow(logT,7)*(3.81073033993717e+06)+pow(logT,8)*(-6.11177474148552e+05)+pow(logT,9)*(3.39029460196043e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 && 3.39<logT && logT <3.89){
				log_ro1=-7.60997085699224e+08+logT*(1.45405341163392e+09)+pow(logT,2)*(-1.14108434367161e+09)+pow(logT,3)*(4.47329206852131e+08)+pow(logT,4)*(-7.09430992654316e+07)+pow(logT,5)*(-1.06178422153952e+07)+pow(logT,6)*(7.26236671954630e+06)+pow(logT,7)*(-1.45748217964055e+06)+pow(logT,8)*(1.39815644499047e+05)+pow(logT,9)*(-5.42853010924804e+03);
				log_ro2=5.12327999813832e+08+logT*(-1.01505440812654e+09)+pow(logT,2)*(8.27170971759619e+08)+pow(logT,3)*(-3.38105658649033e+08)+pow(logT,4)*(5.71975963297892e+07)+pow(logT,5)*(7.64153606342723e+06)+pow(logT,6)*(-5.81792652386035e+06)+pow(logT,7)*(1.22128553946297e+06)+pow(logT,8)*(-1.21778047161494e+05)+pow(logT,9)*(4.90543178391408e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 && 3.39<logT && logT <3.90){
				log_ro1=-1.86206037881216e+09+logT*(3.58747588250676e+09)+pow(logT,2)*(-2.84048576969271e+09)+pow(logT,3)*(1.12510469197357e+09)+pow(logT,4)*(-1.81605223707766e+08)+pow(logT,5)*(-2.62289861241107e+07)+pow(logT,6)*(1.84815366418488e+07)+pow(logT,7)*(-3.75374063372581e+06)+pow(logT,8)*(3.63824222047111e+05)+pow(logT,9)*(-1.42670586800344e+04);
				log_ro2=-7.60997085699224e+08+logT*(1.45405341163392e+09)+pow(logT,2)*(-1.14108434367161e+09)+pow(logT,3)*(4.47329206852131e+08)+pow(logT,4)*(-7.09430992654316e+07)+pow(logT,5)*(-1.06178422153952e+07)+pow(logT,6)*(7.26236671954630e+06)+pow(logT,7)*(-1.45748217964055e+06)+pow(logT,8)*(1.39815644499047e+05)+pow(logT,9)*(-5.42853010924804e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 && 3.39<logT && logT <3.90){
				log_ro1=-4.02299531275594e+09+logT*(7.77912494987396e+09)+pow(logT,2)*(-6.18251305962967e+09)+pow(logT,3)*(2.45855285507876e+09)+pow(logT,4)*(-3.98738882890828e+08)+pow(logT,5)*(-5.75258311698963e+07)+pow(logT,6)*(4.07877150603910e+07)+pow(logT,7)*(-8.31968848050072e+06)+pow(logT,8)*(8.09689497862818e+05)+pow(logT,9)*(-3.18816714959335e+04);
				log_ro2=-1.86206037881216e+09+logT*(3.58747588250676e+09)+pow(logT,2)*(-2.84048576969271e+09)+pow(logT,3)*(1.12510469197357e+09)+pow(logT,4)*(-1.81605223707766e+08)+pow(logT,5)*(-2.62289861241107e+07)+pow(logT,6)*(1.84815366418488e+07)+pow(logT,7)*(-3.75374063372581e+06)+pow(logT,8)*(3.63824222047111e+05)+pow(logT,9)*(-1.42670586800344e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 && 3.39<logT && logT <3.875){
				log_ro1=-3.71251623453310e+09+logT*(7.20082380301755e+09)+pow(logT,2)*(-5.73974400775105e+09)+pow(logT,3)*(2.28842271330620e+09)+pow(logT,4)*(-3.71435826074440e+08)+pow(logT,5)*(-5.43737482122841e+07)+pow(logT,6)*(3.84607193641520e+07)+pow(logT,7)*(-7.86321245329721e+06)+pow(logT,8)*(7.67401871105443e+05)+pow(logT,9)*(-3.03051818393353e+04);
				log_ro2=-4.02299531275594e+09+logT*(7.77912494987396e+09)+pow(logT,2)*(-6.18251305962967e+09)+pow(logT,3)*(2.45855285507876e+09)+pow(logT,4)*(-3.98738882890828e+08)+pow(logT,5)*(-5.75258311698963e+07)+pow(logT,6)*(4.07877150603910e+07)+pow(logT,7)*(-8.31968848050072e+06)+pow(logT,8)*(8.09689497862818e+05)+pow(logT,9)*(-3.18816714959335e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 && 3.39<logT && logT <3.875){
				log_ro1=-3.33717684442625e+09+logT*(6.49185117524158e+09)+pow(logT,2)*(-5.19249781935120e+09)+pow(logT,3)*(2.08048487930235e+09)+pow(logT,4)*(-3.42226488588645e+08)+pow(logT,5)*(-4.75428624923832e+07)+pow(logT,6)*(3.46103028423251e+07)+pow(logT,7)*(-7.11813373012673e+06)+pow(logT,8)*(6.97173957799295e+05)+pow(logT,9)*(-2.76096827694539e+04);
				log_ro2=-3.71251623453310e+09+logT*(7.20082380301755e+09)+pow(logT,2)*(-5.73974400775105e+09)+pow(logT,3)*(2.28842271330620e+09)+pow(logT,4)*(-3.71435826074440e+08)+pow(logT,5)*(-5.43737482122841e+07)+pow(logT,6)*(3.84607193641520e+07)+pow(logT,7)*(-7.86321245329721e+06)+pow(logT,8)*(7.67401871105443e+05)+pow(logT,9)*(-3.03051818393353e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 && 3.39<logT && logT <3.875){
				log_ro1=-3.57738315673669e+09+logT*(6.98840323862895e+09)+pow(logT,2)*(-5.61106120061376e+09)+pow(logT,3)*(2.25444216605134e+09)+pow(logT,4)*(-3.69738572761457e+08)+pow(logT,5)*(-5.35665332125390e+07)+pow(logT,6)*(3.84639872383982e+07)+pow(logT,7)*(-7.92635238234224e+06)+pow(logT,8)*(7.79109847475579e+05)+pow(logT,9)*(-3.09800028744490e+04);
				log_ro2=-3.33717684442625e+09+logT*(6.49185117524158e+09)+pow(logT,2)*(-5.19249781935120e+09)+pow(logT,3)*(2.08048487930235e+09)+pow(logT,4)*(-3.42226488588645e+08)+pow(logT,5)*(-4.75428624923832e+07)+pow(logT,6)*(3.46103028423251e+07)+pow(logT,7)*(-7.11813373012673e+06)+pow(logT,8)*(6.97173957799295e+05)+pow(logT,9)*(-2.76096827694539e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 && 3.39<logT && logT <3.875){
				log_ro1=-3.41560857675066e+09+logT*(6.69264063472343e+09)+pow(logT,2)*(-5.39075940119090e+09)+pow(logT,3)*(2.17385642038511e+09)+pow(logT,4)*(-3.58770650158321e+08)+pow(logT,5)*(-5.12415751396536e+07)+pow(logT,6)*(3.72046984848083e+07)+pow(logT,7)*(-7.69721818362850e+06)+pow(logT,8)*(7.59021026191802e+05)+pow(logT,9)*(-3.02710699100750e+04);
				log_ro2=-3.57738315673669e+09+logT*(6.98840323862895e+09)+pow(logT,2)*(-5.61106120061376e+09)+pow(logT,3)*(2.25444216605134e+09)+pow(logT,4)*(-3.69738572761457e+08)+pow(logT,5)*(-5.35665332125390e+07)+pow(logT,6)*(3.84639872383982e+07)+pow(logT,7)*(-7.92635238234224e+06)+pow(logT,8)*(7.79109847475579e+05)+pow(logT,9)*(-3.09800028744490e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 && 3.39<logT && logT <3.813){
				log_ro1=-1.12857741672386e+07+logT*(1.93014822433621e+07)+pow(logT,2)*(-1.18354218538920e+07)+pow(logT,3)*(2.44375269358945e+06)+pow(logT,4)*(4.06400758270174e+05)+pow(logT,5)*(-2.07055929876577e+05)+pow(logT,6)*(-1.68816884071648e+04)+pow(logT,7)*(1.91962554414115e+04)+pow(logT,8)*(-3.39637526856570e+03)+pow(logT,9)*(1.99332914343933e+02);
				log_ro2=-3.41560857675066e+09+logT*(6.69264063472343e+09)+pow(logT,2)*(-5.39075940119090e+09)+pow(logT,3)*(2.17385642038511e+09)+pow(logT,4)*(-3.58770650158321e+08)+pow(logT,5)*(-5.12415751396536e+07)+pow(logT,6)*(3.72046984848083e+07)+pow(logT,7)*(-7.69721818362850e+06)+pow(logT,8)*(7.59021026191802e+05)+pow(logT,9)*(-3.02710699100750e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5 && 3.39<logT && logT <3.76){
				log_ro1=-1.42393409670356e+08+logT*(2.15084605853568e+08)+pow(logT,2)*(-1.16406513474488e+08)+pow(logT,3)*(1.96345964742703e+07)+pow(logT,4)*(4.56029757200867e+06)+pow(logT,5)*(-1.70172012369979e+06)+pow(logT,6)*(-2.12764228972925e+05)+pow(logT,7)*(1.64542206317499e+05)+pow(logT,8)*(-2.66653846988809e+04)+pow(logT,9)*(1.47337291515449e+03);
				log_ro2=-1.12857741672386e+07+logT*(1.93014822433621e+07)+pow(logT,2)*(-1.18354218538920e+07)+pow(logT,3)*(2.44375269358945e+06)+pow(logT,4)*(4.06400758270174e+05)+pow(logT,5)*(-2.07055929876577e+05)+pow(logT,6)*(-1.68816884071648e+04)+pow(logT,7)*(1.91962554414115e+04)+pow(logT,8)*(-3.39637526856570e+03)+pow(logT,9)*(1.99332914343933e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0 && 3.39<logT && logT <3.70){
				log_ro1=-5.56147225730504e+07+logT*(8.49707575446362e+07)+pow(logT,2)*(-4.63777143629895e+07)+pow(logT,3)*(7.77468237965476e+06)+pow(logT,4)*(1.92382730244288e+06)+pow(logT,5)*(-7.03279095509931e+05)+pow(logT,6)*(-9.48904938825268e+04)+pow(logT,7)*(7.12130753313696e+04)+pow(logT,8)*(-1.15825706738056e+04)+pow(logT,9)*(6.44093768957911e+02);
				log_ro2=-1.42393409670356e+08+logT*(2.15084605853568e+08)+pow(logT,2)*(-1.16406513474488e+08)+pow(logT,3)*(1.96345964742703e+07)+pow(logT,4)*(4.56029757200867e+06)+pow(logT,5)*(-1.70172012369979e+06)+pow(logT,6)*(-2.12764228972925e+05)+pow(logT,7)*(1.64542206317499e+05)+pow(logT,8)*(-2.66653846988809e+04)+pow(logT,9)*(1.47337291515449e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else if ( 3.39<logT && logT <3.59){
				log_ro1=-2.79439840938689e+07+logT*(2.91840929945048e+07)+pow(logT,2)*(-7.73278674175929e+06)+pow(logT,3)*(-1.15436458350101e+06)+pow(logT,4)*(4.85054662101620e+05)+pow(logT,5)*(1.22672648439058e+05)+pow(logT,6)*(-2.88851157783858e+04)+pow(logT,7)*(-1.20559982228738e+04)+pow(logT,8)*(3.87036944393164e+03)+pow(logT,9)*(-3.00559106485898e+02);
				log_ro2=-5.56147225730504e+07+logT*(8.49707575446362e+07)+pow(logT,2)*(-4.63777143629895e+07)+pow(logT,3)*(7.77468237965476e+06)+pow(logT,4)*(1.92382730244288e+06)+pow(logT,5)*(-7.03279095509931e+05)+pow(logT,6)*(-9.48904938825268e+04)+pow(logT,7)*(7.12130753313696e+04)+pow(logT,8)*(-1.15825706738056e+04)+pow(logT,9)*(6.44093768957911e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=5.){
				log_ro1=1.13177e-07;
				log_ro2=5.59252e-06;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 ){ //same as above interpolations, but takes the lowest avalaible density for T outside of models range
				log_ro1=3.69984e-08;
				log_ro2=1.13177e-07;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 ){
				log_ro1=1.59402e-08;
				log_ro2=3.69984e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 ){
				log_ro1=6.0003e-09;
				log_ro2=1.59402e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0){
				log_ro1=5.49184e-09;
				log_ro2=6.0003e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5){
				log_ro1=9.87811e-10;
				log_ro2=5.49184e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0){
				log_ro1=9.55222e-10;
				log_ro2=9.87811e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5){
				log_ro1=3.73616e-10;
				log_ro2=9.55222e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0){
				log_ro1=8.35599e-10;
				log_ro2=3.73616e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5){
				log_ro1=1.13367e-09;
				log_ro2=8.35599e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.){
				log_ro1=1.27143e-09;
				log_ro2=1.13367e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else {
				log_ro1=7.417e-10;
				log_ro2=1.27143e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}


}


double ro_ph_z0122(double logg, double logT)
{ /* Calculate density in the photosphere for the z=0.0122 [M/H]=0*/
        /*density at kappa_ross=2/3 from MbRCS model atmospheres*/
		double log_ro1,log_ro2,log_ro;

		if (logg>=5.  && 3.39<logT && logT <3.60){
				log_ro1=3.72562944722486e+08+logT*(-7.42982049129247e+08)+pow(logT,2)*(6.09559646482419e+08)+pow(logT,3)*(-2.51154235687816e+08)+pow(logT,4)*(4.31760917138012e+07)+pow(logT,5)*(5.45843043871341e+06)+pow(logT,6)*(-4.30479040832773e+06)+pow(logT,7)*(9.11111254044130e+05)+pow(logT,8)*(-9.13429703400950e+04)+pow(logT,9)*(3.69566868242510e+03);  /* log density for lower logg*/
				log_ro2= -1.63577117474210e+08+logT*(2.44351967610962e+08)+pow(logT,2)*(-1.29065441173767e+08)+pow(logT,3)*(1.94843710478105e+07)+pow(logT,4)*(6.07200668755127e+06)+pow(logT,5)*(-1.91448286533038e+06)+pow(logT,6)*(-3.33613304682813e+05)+pow(logT,7)*(2.13651766747734e+05)+pow(logT,8)*(-3.36563264596142e+04)+pow(logT,9)*(1.83839433612598e+03); /* log density for higher logg*/
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5; /* interpolation*/
				return pow(10,log_ro);
		}
		else if (logg>=4.5  && 3.39<logT && logT <3.9 ){
				log_ro1=-5.51431378784203e+08+logT*(1.04237761345322e+09)+pow(logT,2)*(-8.08438815315778e+08)+pow(logT,3)*(3.12519821060760e+08)+pow(logT,4)*(-4.83559378544162e+07)+pow(logT,5)*(-7.59936558442846e+06)+pow(logT,6)*(4.99170173610124e+06)+pow(logT,7)*(-9.84781833565036e+05)+pow(logT,8)*(9.30517052545713e+04)+pow(logT,9)*(-3.55925537049821e+03);
				log_ro2=3.72562944722486e+08+logT*(-7.42982049129247e+08)+pow(logT,2)*(6.09559646482419e+08)+pow(logT,3)*(-2.51154235687816e+08)+pow(logT,4)*(4.31760917138012e+07)+pow(logT,5)*(5.45843043871341e+06)+pow(logT,6)*(-4.30479040832773e+06)+pow(logT,7)*(9.11111254044130e+05)+pow(logT,8)*(-9.13429703400950e+04)+pow(logT,9)*(3.69566868242510e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5; /* interpolation*/
				return pow(10,log_ro);
		}
		else if (logg>=4.0  && 3.39<logT && logT <3.9 ){
				log_ro1=-1.88760411157266e+09+logT*(3.63365344143062e+09)+pow(logT,2)*(-2.87446491982868e+09)+pow(logT,3)*(1.13737677447832e+09)+pow(logT,4)*(-1.83261695735710e+08)+pow(logT,5)*(-2.65641431267486e+07)+pow(logT,6)*(1.86610479642314e+07)+pow(logT,7)*(-3.78563619156599e+06)+pow(logT,8)*(3.66537274148040e+05)+pow(logT,9)*(-1.43592214231473e+04);
				log_ro2=-5.51431378784203e+08+logT*(1.04237761345322e+09)+pow(logT,2)*(-8.08438815315778e+08)+pow(logT,3)*(3.12519821060760e+08)+pow(logT,4)*(-4.83559378544162e+07)+pow(logT,5)*(-7.59936558442846e+06)+pow(logT,6)*(4.99170173610124e+06)+pow(logT,7)*(-9.84781833565036e+05)+pow(logT,8)*(9.30517052545713e+04)+pow(logT,9)*(-3.55925537049821e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.0)/0.5; /* interpolation*/
				return pow(10,log_ro);
		}
		else if (logg>=3.5  && 3.39<logT && logT <3.9 ){
				log_ro1=-2.94157766969062e+09+logT*(5.68743774067660e+09)+pow(logT,2)*(-4.51998431432145e+09)+pow(logT,3)*(1.79785428195545e+09)+pow(logT,4)*(-2.92142649194326e+08)+pow(logT,5)*(-4.16716325369616e+07)+pow(logT,6)*(2.96918071422458e+07)+pow(logT,7)*(-6.05828170307005e+06)+pow(logT,8)*(5.89490317997546e+05)+pow(logT,9)*(-2.32025929407420e+04);
				log_ro2=-1.88760411157266e+09+logT*(3.63365344143062e+09)+pow(logT,2)*(-2.87446491982868e+09)+pow(logT,3)*(1.13737677447832e+09)+pow(logT,4)*(-1.83261695735710e+08)+pow(logT,5)*(-2.65641431267486e+07)+pow(logT,6)*(1.86610479642314e+07)+pow(logT,7)*(-3.78563619156599e+06)+pow(logT,8)*(3.66537274148040e+05)+pow(logT,9)*(-1.43592214231473e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5; /* interpolation*/
				return pow(10,log_ro);
		}
		else if (logg>=3.0  && 3.39<logT && logT <3.9 ){
				log_ro1=-3.24254110862828e+09+logT*(6.28755433627256e+09)+pow(logT,2)*(-5.01227260646349e+09)+pow(logT,3)*(2.00073994971836e+09)+pow(logT,4)*(-3.27126482294010e+08)+pow(logT,5)*(-4.59945417155445e+07)+pow(logT,6)*(3.31362093336355e+07)+pow(logT,7)*(-6.78740855652699e+06)+pow(logT,8)*(6.62523300225756e+05)+pow(logT,9)*(-2.61536399583295e+04);
				log_ro2=-2.94157766969062e+09+logT*(5.68743774067660e+09)+pow(logT,2)*(-4.51998431432145e+09)+pow(logT,3)*(1.79785428195545e+09)+pow(logT,4)*(-2.92142649194326e+08)+pow(logT,5)*(-4.16716325369616e+07)+pow(logT,6)*(2.96918071422458e+07)+pow(logT,7)*(-6.05828170307005e+06)+pow(logT,8)*(5.89490317997546e+05)+pow(logT,9)*(-2.32025929407420e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5; /* interpolation*/
				return pow(10,log_ro);
		}
		else if (logg>=2.5  && 3.39<logT && logT <3.9 ){
				log_ro1=-3.26348578900182e+09+logT*(6.34699513661612e+09)+pow(logT,2)*(-5.07520983736141e+09)+pow(logT,3)*(2.03270751261246e+09)+pow(logT,4)*(-3.34060844653509e+08)+pow(logT,5)*(-4.65575925146776e+07)+pow(logT,6)*(3.38270133972195e+07)+pow(logT,7)*(-6.95358133978459e+06)+pow(logT,8)*(6.80812665880863e+05)+pow(logT,9)*(-2.69530345812752e+04);
				log_ro2=-3.24254110862828e+09+logT*(6.28755433627256e+09)+pow(logT,2)*(-5.01227260646349e+09)+pow(logT,3)*(2.00073994971836e+09)+pow(logT,4)*(-3.27126482294010e+08)+pow(logT,5)*(-4.59945417155445e+07)+pow(logT,6)*(3.31362093336355e+07)+pow(logT,7)*(-6.78740855652699e+06)+pow(logT,8)*(6.62523300225756e+05)+pow(logT,9)*(-2.61536399583295e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5; /* interpolation*/
				return pow(10,log_ro);
		}
		else if (logg>=2.0  && 3.39<logT && logT <3.89 ){
				log_ro1=-2.70957771732986e+09+logT*(5.29697973542320e+09)+pow(logT,2)*(-4.25727708607853e+09)+pow(logT,3)*(1.71374875602300e+09)+pow(logT,4)*(-2.83047339657916e+08)+pow(logT,5)*(-3.96518673892605e+07)+pow(logT,6)*(2.89490416362777e+07)+pow(logT,7)*(-5.97987578073138e+06)+pow(logT,8)*(5.88314065115381e+05)+pow(logT,9)*(-2.34029290084411e+04);
				log_ro2=-3.26348578900182e+09+logT*(6.34699513661612e+09)+pow(logT,2)*(-5.07520983736141e+09)+pow(logT,3)*(2.03270751261246e+09)+pow(logT,4)*(-3.34060844653509e+08)+pow(logT,5)*(-4.65575925146776e+07)+pow(logT,6)*(3.38270133972195e+07)+pow(logT,7)*(-6.95358133978459e+06)+pow(logT,8)*(6.80812665880863e+05)+pow(logT,9)*(-2.69530345812752e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5; /* interpolation*/
				return pow(10,log_ro);
		}
		else if (logg>=1.5  && 3.39<logT && logT <3.86){
				log_ro1=-2.56443660633276e+09+logT*(5.02882846822589e+09)+pow(logT,2)*(-4.05268364969600e+09)+pow(logT,3)*(1.63397932136100e+09)+pow(logT,4)*(-2.68649713827052e+08)+pow(logT,5)*(-3.92924736589749e+07)+pow(logT,6)*(2.82354139011039e+07)+pow(logT,7)*(-5.83667978060229e+06)+pow(logT,8)*(5.75608719814695e+05)+pow(logT,9)*(-2.29644110909452e+04);
				log_ro2=-2.70957771732986e+09+logT*(5.29697973542320e+09)+pow(logT,2)*(-4.25727708607853e+09)+pow(logT,3)*(1.71374875602300e+09)+pow(logT,4)*(-2.83047339657916e+08)+pow(logT,5)*(-3.96518673892605e+07)+pow(logT,6)*(2.89490416362777e+07)+pow(logT,7)*(-5.97987578073138e+06)+pow(logT,8)*(5.88314065115381e+05)+pow(logT,9)*(-2.34029290084411e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5; /* interpolation*/
				return pow(10,log_ro);
		}
		else if (logg>=1.0  && 3.39<logT && logT <3.83){
				log_ro1=2.19999478999620e+07+logT*(-3.02733121355389e+07)+pow(logT,2)*(1.46749934539426e+07)+pow(logT,3)*(-2.02451682468426e+06)+pow(logT,4)*(-5.73190836575056e+05)+pow(logT,5)*(1.62531908888755e+05)+pow(logT,6)*(2.63109907991942e+04)+pow(logT,7)*(-1.48405850685612e+04)+pow(logT,8)*(2.06757532348205e+03)+pow(logT,9)*(-9.90128283799830e+01);
				log_ro2=-2.56443660633276e+09+logT*(5.02882846822589e+09)+pow(logT,2)*(-4.05268364969600e+09)+pow(logT,3)*(1.63397932136100e+09)+pow(logT,4)*(-2.68649713827052e+08)+pow(logT,5)*(-3.92924736589749e+07)+pow(logT,6)*(2.82354139011039e+07)+pow(logT,7)*(-5.83667978060229e+06)+pow(logT,8)*(5.75608719814695e+05)+pow(logT,9)*(-2.29644110909452e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5; /* interpolation*/
				return pow(10,log_ro);
		}
		else if (logg>=0.5   && 3.39<logT && logT <3.74){
				log_ro1=-6.68183809531999e+06+logT*(1.24093221707470e+07)+pow(logT,2)*(-8.09240435310487e+06)+pow(logT,3)*(1.75960470677663e+06)+pow(logT,4)*(2.96374847123585e+05)+pow(logT,5)*(-1.57322842472075e+05)+pow(logT,6)*(-1.32960693749467e+04)+pow(logT,7)*(1.53728690703489e+04)+pow(logT,8)*(-2.77138369499400e+03)+pow(logT,9)*(1.65374744173045e+02);
				log_ro2=2.19999478999620e+07+logT*(-3.02733121355389e+07)+pow(logT,2)*(1.46749934539426e+07)+pow(logT,3)*(-2.02451682468426e+06)+pow(logT,4)*(-5.73190836575056e+05)+pow(logT,5)*(1.62531908888755e+05)+pow(logT,6)*(2.63109907991942e+04)+pow(logT,7)*(-1.48405850685612e+04)+pow(logT,8)*(2.06757532348205e+03)+pow(logT,9)*(-9.90128283799830e+01);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5; /* interpolation*/
				return pow(10,log_ro);
		}
		else if (logg>=0.0  && 3.39<logT && logT <3.698){
				log_ro1=2.45790860948633e+08+logT*(-3.65318160687998e+08)+pow(logT,2)*(1.93438568998122e+08)+pow(logT,3)*(-3.07268312939569e+07)+pow(logT,4)*(-8.12707184158528e+06)+pow(logT,5)*(2.77124518751393e+06)+pow(logT,6)*(4.03112491972697e+05)+pow(logT,7)*(-2.81107599378426e+05)+pow(logT,8)*(4.45799154145304e+04)+pow(logT,9)*(-2.43116158308377e+03);
				log_ro2=-6.68183809531999e+06+logT*(1.24093221707470e+07)+pow(logT,2)*(-8.09240435310487e+06)+pow(logT,3)*(1.75960470677663e+06)+pow(logT,4)*(2.96374847123585e+05)+pow(logT,5)*(-1.57322842472075e+05)+pow(logT,6)*(-1.32960693749467e+04)+pow(logT,7)*(1.53728690703489e+04)+pow(logT,8)*(-2.77138369499400e+03)+pow(logT,9)*(1.65374744173045e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5; /* interpolation*/
				return pow(10,log_ro);
		}
		else if ( 3.39<logT && logT <3.628){
				log_ro1=1.02356781929277e+08+logT*(-1.40925231511771e+08)+pow(logT,2)*(6.72145393784369e+07)+pow(logT,3)*(-8.08965988637997e+06)+pow(logT,4)*(-3.17856866292147e+06)+pow(logT,5)*(7.59443916398889e+05)+pow(logT,6)*(1.72088712399419e+05)+pow(logT,7)*(-8.31357868066995e+04)+pow(logT,8)*(1.12031958820879e+04)+pow(logT,9)*(-5.23749016174935e+02);
				log_ro2=2.45790860948633e+08+logT*(-3.65318160687998e+08)+pow(logT,2)*(1.93438568998122e+08)+pow(logT,3)*(-3.07268312939569e+07)+pow(logT,4)*(-8.12707184158528e+06)+pow(logT,5)*(2.77124518751393e+06)+pow(logT,6)*(4.03112491972697e+05)+pow(logT,7)*(-2.81107599378426e+05)+pow(logT,8)*(4.45799154145304e+04)+pow(logT,9)*(-2.43116158308377e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5; /* interpolation*/
				return pow(10,log_ro);

		}
		else if (logg>=5.){
				log_ro1=8.47326e-08;
				log_ro2=4.50867e-06;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 ){ //same as above interpolations, but takes the lowest avalaible density for T outside of models range
				log_ro1=3.6536e-08;
				log_ro2=8.47326e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 ){
				log_ro1=1.5424e-08;
				log_ro2=3.6536e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 ){
				log_ro1=6.3246e-09;
				log_ro2=1.5424e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0){
				log_ro1=2.51057e-09;
				log_ro2=6.3246e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5){
				log_ro1=9.72314e-10;
				log_ro2=2.51057e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0){
				log_ro1=5.91848e-10;
				log_ro2=9.72314e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5){
				log_ro1=5.8506e-10;
				log_ro2=5.91848e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0){
				log_ro1=5.58732e-10;
				log_ro2=5.8506e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5){
				log_ro1=1.46834e-09;
				log_ro2=5.58732e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0){
				log_ro1=1.11871e-09;
				log_ro2=1.46834e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else {
				log_ro1=7.10175e-10;
				log_ro2=1.11871e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}

}


double ro_ph_z020267(double logg, double logT)
{ /* Calculate density in the photosphere for the z=0.020267 [M/H]=0.25*/
        /*density at kappa_ross=2/3 from MbRCS model atmospheres*/
		double log_ro1,log_ro2,log_ro;

		if (logg>=5.  && 3.38<logT && logT <3.59){
				log_ro1=1.30402255202162e+09+logT*(-2.56256427098979e+09)+pow(logT,2)*(2.07078298408106e+09)+pow(logT,3)*(-8.38610819495776e+08)+pow(logT,4)*(1.39788094401400e+08)+pow(logT,5)*(1.92678864698285e+07)+pow(logT,6)*(-1.42871913245132e+07)+pow(logT,7)*(2.97041868754819e+06)+pow(logT,8)*(-2.93868766522628e+05)+pow(logT,9)*(1.17519786703001e+04);
				log_ro2=-6.85457835827644e+08+logT*(1.03113556835141e+09)+pow(logT,2)*(-5.49326263578823e+08)+pow(logT,3)*(8.44421187833202e+07)+pow(logT,4)*(2.57483418622117e+07)+pow(logT,5)*(-8.30718796153958e+06)+pow(logT,6)*(-1.41211600888382e+06)+pow(logT,7)*(9.26341004272580e+05)+pow(logT,8)*(-1.47437795568011e+05)+pow(logT,9)*(8.12594339488941e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 && 3.38<logT && logT <3.89){
				log_ro1=-5.19090551816300e+08+logT*(9.72877215342903e+08)+pow(logT,2)*(-7.47347998969125e+08)+pow(logT,3)*(2.85563818313146e+08)+pow(logT,4)*(-4.32584242308917e+07)+pow(logT,5)*(-7.08279109442964e+06)+pow(logT,6)*(4.49820455202969e+06)+pow(logT,7)*(-8.74207308006581e+05)+pow(logT,8)*(8.14849269378181e+04)+pow(logT,9)*(-3.07399559686238e+03);
				log_ro2=1.30402255202162e+09+logT*(-2.56256427098979e+09)+pow(logT,2)*(2.07078298408106e+09)+pow(logT,3)*(-8.38610819495776e+08)+pow(logT,4)*(1.39788094401400e+08)+pow(logT,5)*(1.92678864698285e+07)+pow(logT,6)*(-1.42871913245132e+07)+pow(logT,7)*(2.97041868754819e+06)+pow(logT,8)*(-2.93868766522628e+05)+pow(logT,9)*(1.17519786703001e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 && 3.38<logT && logT <3.90){
				log_ro1=-1.86580451813831e+09+logT*(3.58679102039791e+09)+pow(logT,2)*(-2.83321761633473e+09)+pow(logT,3)*(1.11913388965081e+09)+pow(logT,4)*(-1.79796095669848e+08)+pow(logT,5)*(-2.62176920530368e+07)+pow(logT,6)*(1.83259204711049e+07)+pow(logT,7)*(-3.71023333206884e+06)+pow(logT,8)*(3.58620783730947e+05)+pow(logT,9)*(-1.40259113879676e+04);
				log_ro2=-5.19090551816300e+08+logT*(9.72877215342903e+08)+pow(logT,2)*(-7.47347998969125e+08)+pow(logT,3)*(2.85563818313146e+08)+pow(logT,4)*(-4.32584242308917e+07)+pow(logT,5)*(-7.08279109442964e+06)+pow(logT,6)*(4.49820455202969e+06)+pow(logT,7)*(-8.74207308006581e+05)+pow(logT,8)*(8.14849269378181e+04)+pow(logT,9)*(-3.07399559686238e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 && 3.38<logT && logT <3.90){
				log_ro1=-2.72199639609406e+09+logT*(5.26319717593510e+09)+pow(logT,2)*(-4.18306524654994e+09)+pow(logT,3)*(1.66398834526459e+09)+pow(logT,4)*(-2.70485131052354e+08)+pow(logT,5)*(-3.85104237758109e+07)+pow(logT,6)*(2.74620944744486e+07)+pow(logT,7)*(-5.60359966112207e+06)+pow(logT,8)*(5.45221586751240e+05)+pow(logT,9)*(-2.14582023883681e+04);
				log_ro2=-1.86580451813831e+09+logT*(3.58679102039791e+09)+pow(logT,2)*(-2.83321761633473e+09)+pow(logT,3)*(1.11913388965081e+09)+pow(logT,4)*(-1.79796095669848e+08)+pow(logT,5)*(-2.62176920530368e+07)+pow(logT,6)*(1.83259204711049e+07)+pow(logT,7)*(-3.71023333206884e+06)+pow(logT,8)*(3.58620783730947e+05)+pow(logT,9)*(-1.40259113879676e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 && 3.38<logT && logT <3.90){
				log_ro1=-2.76227969884688e+09+logT*(5.35765636987763e+09)+pow(logT,2)*(-4.27207672630149e+09)+pow(logT,3)*(1.70582223073538e+09)+pow(logT,4)*(-2.79135907177692e+08)+pow(logT,5)*(-3.91132865534508e+07)+pow(logT,6)*(2.82285877124681e+07)+pow(logT,7)*(-5.78369043441449e+06)+pow(logT,8)*(5.64597111382800e+05)+pow(logT,9)*(-2.22880964643258e+04);
				log_ro2=-2.72199639609406e+09+logT*(5.26319717593510e+09)+pow(logT,2)*(-4.18306524654994e+09)+pow(logT,3)*(1.66398834526459e+09)+pow(logT,4)*(-2.70485131052354e+08)+pow(logT,5)*(-3.85104237758109e+07)+pow(logT,6)*(2.74620944744486e+07)+pow(logT,7)*(-5.60359966112207e+06)+pow(logT,8)*(5.45221586751240e+05)+pow(logT,9)*(-2.14582023883681e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 && 3.38<logT && logT <3.90){
				log_ro1=-2.52314880447664e+09+logT*(4.91293878489423e+09)+pow(logT,2)*(-3.93320653634114e+09)+pow(logT,3)*(1.57742286385573e+09)+pow(logT,4)*(-2.59844309058614e+08)+pow(logT,5)*(-3.59966566048517e+07)+pow(logT,6)*(2.62660877730261e+07)+pow(logT,7)*(-5.40632285092662e+06)+pow(logT,8)*(5.29827471938384e+05)+pow(logT,9)*(-2.09927395122448e+04);
				log_ro2=-2.76227969884688e+09+logT*(5.35765636987763e+09)+pow(logT,2)*(-4.27207672630149e+09)+pow(logT,3)*(1.70582223073538e+09)+pow(logT,4)*(-2.79135907177692e+08)+pow(logT,5)*(-3.91132865534508e+07)+pow(logT,6)*(2.82285877124681e+07)+pow(logT,7)*(-5.78369043441449e+06)+pow(logT,8)*(5.64597111382800e+05)+pow(logT,9)*(-2.22880964643258e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 && 3.38<logT && logT <3.90){
				log_ro1=-1.36600613808314e+09+logT*(2.68186925173523e+09)+pow(logT,2)*(-2.16512590817092e+09)+pow(logT,3)*(8.76242991862758e+08)+pow(logT,4)*(-1.46330029602088e+08)+pow(logT,5)*(-1.97699889756391e+07)+pow(logT,6)*(1.47629831747251e+07)+pow(logT,7)*(-3.06640595733440e+06)+pow(logT,8)*(3.02792339746741e+05)+pow(logT,9)*(-1.20812767137178e+04);
				log_ro2=-2.52314880447664e+09+logT*(4.91293878489423e+09)+pow(logT,2)*(-3.93320653634114e+09)+pow(logT,3)*(1.57742286385573e+09)+pow(logT,4)*(-2.59844309058614e+08)+pow(logT,5)*(-3.59966566048517e+07)+pow(logT,6)*(2.62660877730261e+07)+pow(logT,7)*(-5.40632285092662e+06)+pow(logT,8)*(5.29827471938384e+05)+pow(logT,9)*(-2.09927395122448e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 && 3.38<logT && logT <3.875){
				log_ro1=-1.78006379551631e+09+logT*(3.49493096724952e+09)+pow(logT,2)*(-2.82060918652445e+09)+pow(logT,3)*(1.13975928321418e+09)+pow(logT,4)*(-1.88685006823635e+08)+pow(logT,5)*(-2.67859882045320e+07)+pow(logT,6)*(1.95461071083159e+07)+pow(logT,7)*(-4.05090959710901e+06)+pow(logT,8)*(3.99983319145710e+05)+pow(logT,9)*(-1.59696664885254e+04);
				log_ro2=-1.36600613808314e+09+logT*(2.68186925173523e+09)+pow(logT,2)*(-2.16512590817092e+09)+pow(logT,3)*(8.76242991862758e+08)+pow(logT,4)*(-1.46330029602088e+08)+pow(logT,5)*(-1.97699889756391e+07)+pow(logT,6)*(1.47629831747251e+07)+pow(logT,7)*(-3.06640595733440e+06)+pow(logT,8)*(3.02792339746741e+05)+pow(logT,9)*(-1.20812767137178e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 && 3.38<logT && logT <3.83){
				log_ro1=2.14157576557511e+07+logT*(-2.99297762407491e+07)+pow(logT,2)*(1.48094885530827e+07)+pow(logT,3)*(-2.14202716453107e+06)+pow(logT,4)*(-5.70843681158635e+05)+pow(logT,5)*(1.72983111311192e+05)+pow(logT,6)*(2.59617746776702e+04)+pow(logT,7)*(-1.58029886516734e+04)+pow(logT,8)*(2.28639060018568e+03)+pow(logT,9)*(-1.13760415808922e+02);
				log_ro2=-1.78006379551631e+09+logT*(3.49493096724952e+09)+pow(logT,2)*(-2.82060918652445e+09)+pow(logT,3)*(1.13975928321418e+09)+pow(logT,4)*(-1.88685006823635e+08)+pow(logT,5)*(-2.67859882045320e+07)+pow(logT,6)*(1.95461071083159e+07)+pow(logT,7)*(-4.05090959710901e+06)+pow(logT,8)*(3.99983319145710e+05)+pow(logT,9)*(-1.59696664885254e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5 && 3.38<logT && logT <3.76){
				log_ro1=3.06956216320742e+07+logT*(-4.34977652558241e+07)+pow(logT,2)*(2.18181352881293e+07)+pow(logT,3)*(-3.17136128844836e+06)+pow(logT,4)*(-8.90220935585520e+05)+pow(logT,5)*(2.71079819859289e+05)+pow(logT,6)*(4.26559073931948e+04)+pow(logT,7)*(-2.62676150992574e+04)+pow(logT,8)*(3.89953833756298e+03)+pow(logT,9)*(-1.99992826093036e+02);
				log_ro2=2.14157576557511e+07+logT*(-2.99297762407491e+07)+pow(logT,2)*(1.48094885530827e+07)+pow(logT,3)*(-2.14202716453107e+06)+pow(logT,4)*(-5.70843681158635e+05)+pow(logT,5)*(1.72983111311192e+05)+pow(logT,6)*(2.59617746776702e+04)+pow(logT,7)*(-1.58029886516734e+04)+pow(logT,8)*(2.28639060018568e+03)+pow(logT,9)*(-1.13760415808922e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0 && 3.38<logT && logT <3.70){
				log_ro1=6.92898244373036e+08+logT*(-1.03415990408406e+09)+pow(logT,2)*(5.50155197521192e+08)+pow(logT,3)*(-8.81119333127619e+07)+pow(logT,4)*(-2.30954011626835e+07)+pow(logT,5)*(7.95935810732258e+06)+pow(logT,6)*(1.14495336660363e+06)+pow(logT,7)*(-8.07621195324286e+05)+pow(logT,8)*(1.28635099352682e+05)+pow(logT,9)*(-7.03976625677868e+03);
				log_ro2=3.06956216320742e+07+logT*(-4.34977652558241e+07)+pow(logT,2)*(2.18181352881293e+07)+pow(logT,3)*(-3.17136128844836e+06)+pow(logT,4)*(-8.90220935585520e+05)+pow(logT,5)*(2.71079819859289e+05)+pow(logT,6)*(4.26559073931948e+04)+pow(logT,7)*(-2.62676150992574e+04)+pow(logT,8)*(3.89953833756298e+03)+pow(logT,9)*(-1.99992826093036e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else if (3.38<logT && logT <3.65){
				log_ro1=-5.62459300000717e+07+logT*(5.92280081797911e+07)+pow(logT,2)*(-1.61423355035231e+07)+pow(logT,3)*(-2.13080850558392e+06)+pow(logT,4)*(9.96427081545735e+05)+pow(logT,5)*(2.27232185599153e+05)+pow(logT,6)*(-5.90140658283231e+04)+pow(logT,7)*(-2.21313881436735e+04)+pow(logT,8)*(7.34502885591805e+03)+pow(logT,9)*(-5.75375654290358e+02);
				log_ro2=6.92898244373036e+08+logT*(-1.03415990408406e+09)+pow(logT,2)*(5.50155197521192e+08)+pow(logT,3)*(-8.81119333127619e+07)+pow(logT,4)*(-2.30954011626835e+07)+pow(logT,5)*(7.95935810732258e+06)+pow(logT,6)*(1.14495336660363e+06)+pow(logT,7)*(-8.07621195324286e+05)+pow(logT,8)*(1.28635099352682e+05)+pow(logT,9)*(-7.03976625677868e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=5.){
				log_ro1=9.87142e-08;
				log_ro2=3.30872e-06;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 ){ //same as above interpolations, but takes the lowest avalaible density for T outside of models range
				log_ro1=3.28221e-08;
				log_ro2=9.87142e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 ){
				log_ro1=1.43618e-08;
				log_ro2=3.28221e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 ){
				log_ro1=5.27932e-09;
				log_ro2=1.43618e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0){
				log_ro1=2.29884e-09;
				log_ro2=5.27932e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5){
				log_ro1=8.8786e-10;
				log_ro2=2.29884e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0){
				log_ro1=3.20505e-10;
				log_ro2=8.8786e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5){
				log_ro1=3.29842e-10;
				log_ro2=3.20505e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0){
				log_ro1=4.98129e-10;
				log_ro2=3.29842e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5){
				log_ro1=9.26915e-10;
				log_ro2=4.98129e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0){
				log_ro1=8.76866e-10;
				log_ro2=9.26915e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else {
				log_ro1=4.50841e-10;
				log_ro2=8.76866e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}
}


double ro_ph_z0322626(double logg, double logT)
{ /* Calculate density in the photosphere for the z=0.0322626 [M/H]=0.5*/
        /*density at kappa_ross=2/3 from MbRCS model atmospheres*/
		double log_ro1,log_ro2,log_ro;

		if (logg>=5. && 3.40<logT && logT <3.59){
				log_ro1=1.33035798301649e+09+logT*(-2.61772537269369e+09)+pow(logT,2)*(2.11811762975080e+09)+pow(logT,3)*(-8.58982042726724e+08)+pow(logT,4)*(1.43493839511541e+08)+pow(logT,5)*(1.96910225988867e+07)+pow(logT,6)*(-1.46554153319821e+07)+pow(logT,7)*(3.05099627042869e+06)+pow(logT,8)*(-3.02153422015561e+05)+pow(logT,9)*(1.20943721123011e+04);
				log_ro2=-1.53002214401450e+09+logT*(2.30753845815003e+09)+pow(logT,2)*(-1.23310536930720e+09)+pow(logT,3)*(1.90752150010461e+08)+pow(logT,4)*(5.77254302446854e+07)+pow(logT,5)*(-1.87730315147328e+07)+pow(logT,6)*(-3.16393872161793e+06)+pow(logT,7)*(2.09278829550933e+06)+pow(logT,8)*(-3.34238216247867e+05)+pow(logT,9)*(1.84752779667628e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 && 3.40<logT && logT <3.89){
				log_ro1=-4.99517620690132e+08+logT*(9.31521139902224e+08)+pow(logT,2)*(-7.11530849876706e+08)+pow(logT,3)*(2.69980227947857e+08)+pow(logT,4)*(-4.03656572240896e+07)+pow(logT,5)*(-6.77591343417576e+06)+pow(logT,6)*(4.21659265327343e+06)+pow(logT,7)*(-8.11802004988973e+05)+pow(logT,8)*(7.50099114883699e+04)+pow(logT,9)*(-2.80423487014122e+03);
				log_ro2=1.33035798301649e+09+logT*(-2.61772537269369e+09)+pow(logT,2)*(2.11811762975080e+09)+pow(logT,3)*(-8.58982042726724e+08)+pow(logT,4)*(1.43493839511541e+08)+pow(logT,5)*(1.96910225988867e+07)+pow(logT,6)*(-1.46554153319821e+07)+pow(logT,7)*(3.05099627042869e+06)+pow(logT,8)*(-3.02153422015561e+05)+pow(logT,9)*(1.20943721123011e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 && 3.40<logT && logT <3.90){
				log_ro1=-1.62378336638158e+09+logT*(3.11693177338894e+09)+pow(logT,2)*(-2.45812155463453e+09)+pow(logT,3)*(9.69136566320211e+08)+pow(logT,4)*(-1.55192669720879e+08)+pow(logT,5)*(-2.27806106893682e+07)+pow(logT,6)*(1.58352629940297e+07)+pow(logT,7)*(-3.19878106836422e+06)+pow(logT,8)*(3.08583850053133e+05)+pow(logT,9)*(-1.20461643535228e+04);
				log_ro2=-4.99517620690132e+08+logT*(9.31521139902224e+08)+pow(logT,2)*(-7.11530849876706e+08)+pow(logT,3)*(2.69980227947857e+08)+pow(logT,4)*(-4.03656572240896e+07)+pow(logT,5)*(-6.77591343417576e+06)+pow(logT,6)*(4.21659265327343e+06)+pow(logT,7)*(-8.11802004988973e+05)+pow(logT,8)*(7.50099114883699e+04)+pow(logT,9)*(-2.80423487014122e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 && 3.40<logT && logT <3.90){
				log_ro1=-2.78065317187590e+09+logT*(5.37661899227576e+09)+pow(logT,2)*(-4.27319524872587e+09)+pow(logT,3)*(1.69982652389114e+09)+pow(logT,4)*(-2.76304739306947e+08)+pow(logT,5)*(-3.93409957842498e+07)+pow(logT,6)*(2.80531251025577e+07)+pow(logT,7)*(-5.72407844101364e+06)+pow(logT,8)*(5.56932682294314e+05)+pow(logT,9)*(-2.19186512939806e+04);
				log_ro2=-1.62378336638158e+09+logT*(3.11693177338894e+09)+pow(logT,2)*(-2.45812155463453e+09)+pow(logT,3)*(9.69136566320211e+08)+pow(logT,4)*(-1.55192669720879e+08)+pow(logT,5)*(-2.27806106893682e+07)+pow(logT,6)*(1.58352629940297e+07)+pow(logT,7)*(-3.19878106836422e+06)+pow(logT,8)*(3.08583850053133e+05)+pow(logT,9)*(-1.20461643535228e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 && 3.40<logT && logT <3.90){
				log_ro1=-1.64013216927517e+09+logT*(3.18485030976049e+09)+pow(logT,2)*(-2.54232537133545e+09)+pow(logT,3)*(1.01629550345387e+09)+pow(logT,4)*(-1.66609443643469e+08)+pow(logT,5)*(-2.32381897587812e+07)+pow(logT,6)*(1.68244941185127e+07)+pow(logT,7)*(-3.44998954321679e+06)+pow(logT,8)*(3.36952743103081e+05)+pow(logT,9)*(-1.33059909121005e+04);
				log_ro2=-2.78065317187590e+09+logT*(5.37661899227576e+09)+pow(logT,2)*(-4.27319524872587e+09)+pow(logT,3)*(1.69982652389114e+09)+pow(logT,4)*(-2.76304739306947e+08)+pow(logT,5)*(-3.93409957842498e+07)+pow(logT,6)*(2.80531251025577e+07)+pow(logT,7)*(-5.72407844101364e+06)+pow(logT,8)*(5.56932682294314e+05)+pow(logT,9)*(-2.19186512939806e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 && 3.40<logT && logT <3.90){
				log_ro1=-1.47016255806370e+09+logT*(2.87068398535683e+09)+pow(logT,2)*(-2.30463952484058e+09)+pow(logT,3)*(9.27034323369823e+08)+pow(logT,4)*(-1.53406160518903e+08)+pow(logT,5)*(-2.10556586546872e+07)+pow(logT,6)*(1.54832294976994e+07)+pow(logT,7)*(-3.19562032687481e+06)+pow(logT,8)*(3.13843518180201e+05)+pow(logT,9)*(-1.24583417365978e+04);
				log_ro2=-1.64013216927517e+09+logT*(3.18485030976049e+09)+pow(logT,2)*(-2.54232537133545e+09)+pow(logT,3)*(1.01629550345387e+09)+pow(logT,4)*(-1.66609443643469e+08)+pow(logT,5)*(-2.32381897587812e+07)+pow(logT,6)*(1.68244941185127e+07)+pow(logT,7)*(-3.44998954321679e+06)+pow(logT,8)*(3.36952743103081e+05)+pow(logT,9)*(-1.33059909121005e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 && 3.40<logT && logT <3.89){
				log_ro1=-5.39503869005832e+08+logT*(1.06938531570490e+09)+pow(logT,2)*(-8.71275405112784e+08)+pow(logT,3)*(3.55844360434289e+08)+pow(logT,4)*(-6.01154149473987e+07)+pow(logT,5)*(-8.01439787147283e+06)+pow(logT,6)*(6.08700438460290e+06)+pow(logT,7)*(-1.27394227036272e+06)+pow(logT,8)*(1.26586169967809e+05)+pow(logT,9)*(-5.07900765707611e+03);
				log_ro2=-1.47016255806370e+09+logT*(2.87068398535683e+09)+pow(logT,2)*(-2.30463952484058e+09)+pow(logT,3)*(9.27034323369823e+08)+pow(logT,4)*(-1.53406160518903e+08)+pow(logT,5)*(-2.10556586546872e+07)+pow(logT,6)*(1.54832294976994e+07)+pow(logT,7)*(-3.19562032687481e+06)+pow(logT,8)*(3.13843518180201e+05)+pow(logT,9)*(-1.24583417365978e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 && 3.40<logT && logT <3.86){
				log_ro1=-5.04592626502153e+08+logT*(1.00027699188384e+09)+pow(logT,2)*(-8.14558776888360e+08)+pow(logT,3)*(3.31961805064759e+08)+pow(logT,4)*(-5.54568874973636e+07)+pow(logT,5)*(-7.86307682806966e+06)+pow(logT,6)*(5.79657277878218e+06)+pow(logT,7)*(-1.20881326267473e+06)+pow(logT,8)*(1.19992146554115e+05)+pow(logT,9)*(-4.81318473453511e+03);
				log_ro2=-5.39503869005832e+08+logT*(1.06938531570490e+09)+pow(logT,2)*(-8.71275405112784e+08)+pow(logT,3)*(3.55844360434289e+08)+pow(logT,4)*(-6.01154149473987e+07)+pow(logT,5)*(-8.01439787147283e+06)+pow(logT,6)*(6.08700438460290e+06)+pow(logT,7)*(-1.27394227036272e+06)+pow(logT,8)*(1.26586169967809e+05)+pow(logT,9)*(-5.07900765707611e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 && 3.40<logT && logT <3.83){
				log_ro1=2.61186296400981e+07+logT*(-3.71377760847601e+07)+pow(logT,2)*(1.87826604429447e+07)+pow(logT,3)*(-2.84308968611274e+06)+pow(logT,4)*(-7.17792803799119e+05)+pow(logT,5)*(2.31982181730601e+05)+pow(logT,6)*(3.24652139041828e+04)+pow(logT,7)*(-2.13083589587864e+04)+pow(logT,8)*(3.19373512098751e+03)+pow(logT,9)*(-1.64387287951415e+02);
				log_ro2=-5.04592626502153e+08+logT*(1.00027699188384e+09)+pow(logT,2)*(-8.14558776888360e+08)+pow(logT,3)*(3.31961805064759e+08)+pow(logT,4)*(-5.54568874973636e+07)+pow(logT,5)*(-7.86307682806966e+06)+pow(logT,6)*(5.79657277878218e+06)+pow(logT,7)*(-1.20881326267473e+06)+pow(logT,8)*(1.19992146554115e+05)+pow(logT,9)*(-4.81318473453511e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5 && 3.40<logT && logT <3.78){
				log_ro1=6.77869887832016e+07+logT*(-9.90477597706682e+07)+pow(logT,2)*(5.16476687151910e+07)+pow(logT,3)*(-8.18334737641531e+06)+pow(logT,4)*(-2.02667649113262e+06)+pow(logT,5)*(6.95210787443338e+05)+pow(logT,6)*(9.42454870510987e+04)+pow(logT,7)*(-6.62719814593568e+04)+pow(logT,8)*(1.03361130964135e+04)+pow(logT,9)*(-5.52780118668655e+02);
				log_ro2=2.61186296400981e+07+logT*(-3.71377760847601e+07)+pow(logT,2)*(1.87826604429447e+07)+pow(logT,3)*(-2.84308968611274e+06)+pow(logT,4)*(-7.17792803799119e+05)+pow(logT,5)*(2.31982181730601e+05)+pow(logT,6)*(3.24652139041828e+04)+pow(logT,7)*(-2.13083589587864e+04)+pow(logT,8)*(3.19373512098751e+03)+pow(logT,9)*(-1.64387287951415e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0 && 3.40 <logT && logT <3.70){
				log_ro1=1.05821794990942e+09+logT*(-1.58074477275269e+09)+pow(logT,2)*(8.41724035843423e+08)+pow(logT,3)*(-1.35039983571834e+08)+pow(logT,4)*(-3.53204985391349e+07)+pow(logT,5)*(1.21973423668651e+07)+pow(logT,6)*(1.75218409631451e+06)+pow(logT,7)*(-1.23820619951356e+06)+pow(logT,8)*(1.97370401938842e+05)+pow(logT,9)*(-1.08083330564766e+04);
				log_ro2=6.77869887832016e+07+logT*(-9.90477597706682e+07)+pow(logT,2)*(5.16476687151910e+07)+pow(logT,3)*(-8.18334737641531e+06)+pow(logT,4)*(-2.02667649113262e+06)+pow(logT,5)*(6.95210787443338e+05)+pow(logT,6)*(9.42454870510987e+04)+pow(logT,7)*(-6.62719814593568e+04)+pow(logT,8)*(1.03361130964135e+04)+pow(logT,9)*(-5.52780118668655e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else if ( 3.518<logT && logT <3.65){
				log_ro1=-1.52584172962875e+08+logT*(1.59012693660538e+08)+pow(logT,2)*(-4.24927639578706e+07)+pow(logT,3)*(-5.91476225409080e+06)+pow(logT,4)*(2.61044955176232e+06)+pow(logT,5)*(6.21336023359721e+05)+pow(logT,6)*(-1.53137472918152e+05)+pow(logT,7)*(-5.99677308915515e+04)+pow(logT,8)*(1.94684527135826e+04)+pow(logT,9)*(-1.50942338605147e+03);
				log_ro2=1.05821794990942e+09+logT*(-1.58074477275269e+09)+pow(logT,2)*(8.41724035843423e+08)+pow(logT,3)*(-1.35039983571834e+08)+pow(logT,4)*(-3.53204985391349e+07)+pow(logT,5)*(1.21973423668651e+07)+pow(logT,6)*(1.75218409631451e+06)+pow(logT,7)*(-1.23820619951356e+06)+pow(logT,8)*(1.97370401938842e+05)+pow(logT,9)*(-1.08083330564766e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=5.){
				log_ro1=9.08959e-08;
				log_ro2=2.45826e-06;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 ){ //same as above interpolations, but takes the lowest avalaible density for T outside of models range
				log_ro1=3.05497e-08;
				log_ro2=9.08959e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 ){
				log_ro1=1.34936e-08;
				log_ro2=3.05497e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 ){
				log_ro1=4.90345e-09;
				log_ro2=1.34936e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0){
				log_ro1=2.16249e-09;
				log_ro2=4.90345e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5){
				log_ro1=8.349e-10;
				log_ro2=2.16249e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0){
				log_ro1=5.02332e-10;
				log_ro2=8.349e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5){
				log_ro1=4.86622e-10;
				log_ro2=5.02332e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0){
				log_ro1=4.57292e-10;
				log_ro2=4.86622e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5){
				log_ro1=5.95463e-10;
				log_ro2=4.57292e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0){
				log_ro1=6.81198e-10;
				log_ro2=5.95463e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else {
				log_ro1=3.61416e-10;
				log_ro2=6.81198e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}
}


double ro_ph_z048359(double logg, double logT)
{ /* Calculate density in the photosphere for the z=0.048359 [M/H]=0.75*/
        /*density at kappa_ross=2/3 from MbRCS model atmospheres*/
		double log_ro1,log_ro2,log_ro;

		if (logg>=5.	 && 3.40<logT && logT <3.59){
				log_ro1=1.11892691053898e+09+logT*(-2.20771234341245e+09)+pow(logT,2)*(1.79126555399042e+09)+pow(logT,3)*(-7.28586742349360e+08)+pow(logT,4)*(1.22272954011643e+08)+pow(logT,5)*(1.66195066980619e+07)+pow(logT,6)*(-1.24691272492868e+07)+pow(logT,7)*(2.60327980577652e+06)+pow(logT,8)*(-2.58399946717365e+05)+pow(logT,9)*(1.03640816043229e+04);
				log_ro2=-2.36838144896111e+09+logT*(3.57843923018477e+09)+pow(logT,2)*(-1.91637563945484e+09)+pow(logT,3)*(2.97741284569635e+08)+pow(logT,4)*(8.96339147143956e+07)+pow(logT,5)*(-2.93099217267947e+07)+pow(logT,6)*(-4.91091405396368e+06)+pow(logT,7)*(3.26674076745484e+06)+pow(logT,8)*(-5.22934439920137e+05)+pow(logT,9)*(2.89617055477018e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 && 3.40<logT && logT <3.889){
				log_ro1=-4.90056113900056e+08+logT*(9.13791733567291e+08)+pow(logT,2)*(-6.97910759241884e+08)+pow(logT,3)*(2.64775379773060e+08)+pow(logT,4)*(-3.95773664347985e+07)+pow(logT,5)*(-6.64667874437207e+06)+pow(logT,6)*(4.13461698318072e+06)+pow(logT,7)*(-7.95892815741018e+05)+pow(logT,8)*(7.35305147995215e+04)+pow(logT,9)*(-2.74863139998055e+03);
				log_ro2=1.11892691053898e+09+logT*(-2.20771234341245e+09)+pow(logT,2)*(1.79126555399042e+09)+pow(logT,3)*(-7.28586742349360e+08)+pow(logT,4)*(1.22272954011643e+08)+pow(logT,5)*(1.66195066980619e+07)+pow(logT,6)*(-1.24691272492868e+07)+pow(logT,7)*(2.60327980577652e+06)+pow(logT,8)*(-2.58399946717365e+05)+pow(logT,9)*(1.03640816043229e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 && 3.40<logT && logT <3.90){
				log_ro1=-1.14865636393885e+09+logT*(2.19621418829377e+09)+pow(logT,2)*(-1.72455601545637e+09)+pow(logT,3)*(6.76466979011767e+08)+pow(logT,4)*(-1.07373282495493e+08)+pow(logT,5)*(-1.60457955219010e+07)+pow(logT,6)*(1.09882762279631e+07)+pow(logT,7)*(-2.20610363077460e+06)+pow(logT,8)*(2.11685980150734e+05)+pow(logT,9)*(-8.22054218715309e+03);
				log_ro2=-4.90056113900056e+08+logT*(9.13791733567291e+08)+pow(logT,2)*(-6.97910759241884e+08)+pow(logT,3)*(2.64775379773060e+08)+pow(logT,4)*(-3.95773664347985e+07)+pow(logT,5)*(-6.64667874437207e+06)+pow(logT,6)*(4.13461698318072e+06)+pow(logT,7)*(-7.95892815741018e+05)+pow(logT,8)*(7.35305147995215e+04)+pow(logT,9)*(-2.74863139998055e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 && 3.40<logT && logT <3.90){
				log_ro1=-1.88951544050691e+09+logT*(3.64845717852923e+09)+pow(logT,2)*(-2.89519832702044e+09)+pow(logT,3)*(1.14947065155105e+09)+pow(logT,4)*(-1.86153600834626e+08)+pow(logT,5)*(-2.67629959502864e+07)+pow(logT,6)*(1.89522997223441e+07)+pow(logT,7)*(-3.85807192964051e+06)+pow(logT,8)*(3.74660865882790e+05)+pow(logT,9)*(-1.47185354619849e+04);
				log_ro2=-1.14865636393885e+09+logT*(2.19621418829377e+09)+pow(logT,2)*(-1.72455601545637e+09)+pow(logT,3)*(6.76466979011767e+08)+pow(logT,4)*(-1.07373282495493e+08)+pow(logT,5)*(-1.60457955219010e+07)+pow(logT,6)*(1.09882762279631e+07)+pow(logT,7)*(-2.20610363077460e+06)+pow(logT,8)*(2.11685980150734e+05)+pow(logT,9)*(-8.22054218715309e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 && 3.40<logT && logT <3.90){
				log_ro1=-1.26141142201698e+09+logT*(2.44984641484224e+09)+pow(logT,2)*(-1.95566437066347e+09)+pow(logT,3)*(7.81599442396925e+08)+pow(logT,4)*(-1.27958617764292e+08)+pow(logT,5)*(-1.79738231227986e+07)+pow(logT,6)*(1.29670622819274e+07)+pow(logT,7)*(-2.65736514545156e+06)+pow(logT,8)*(2.59444143990207e+05)+pow(logT,9)*(-1.02419112813295e+04);
				log_ro2=-1.88951544050691e+09+logT*(3.64845717852923e+09)+pow(logT,2)*(-2.89519832702044e+09)+pow(logT,3)*(1.14947065155105e+09)+pow(logT,4)*(-1.86153600834626e+08)+pow(logT,5)*(-2.67629959502864e+07)+pow(logT,6)*(1.89522997223441e+07)+pow(logT,7)*(-3.85807192964051e+06)+pow(logT,8)*(3.74660865882790e+05)+pow(logT,9)*(-1.47185354619849e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 && 3.40<logT && logT <3.90){
				log_ro1=-4.82298469449122e+08+logT*(9.52166114274014e+08)+pow(logT,2)*(-7.72654088498897e+08)+pow(logT,3)*(3.14299000548327e+08)+pow(logT,4)*(-5.28908642170379e+07)+pow(logT,5)*(-7.01444850135997e+06)+pow(logT,6)*(5.30847769037309e+06)+pow(logT,7)*(-1.10646523640039e+06)+pow(logT,8)*(1.09487589753818e+05)+pow(logT,9)*(-4.37453072495502e+03);
				log_ro2=-1.26141142201698e+09+logT*(2.44984641484224e+09)+pow(logT,2)*(-1.95566437066347e+09)+pow(logT,3)*(7.81599442396925e+08)+pow(logT,4)*(-1.27958617764292e+08)+pow(logT,5)*(-1.79738231227986e+07)+pow(logT,6)*(1.29670622819274e+07)+pow(logT,7)*(-2.65736514545156e+06)+pow(logT,8)*(2.59444143990207e+05)+pow(logT,9)*(-1.02419112813295e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 && 3.40<logT && logT <3.90){
				log_ro1=1.76176628909089e+07+logT*(-1.72576525169566e+07)+pow(logT,2)*(1.10137006131467e+05)+pow(logT,3)*(5.88773417212743e+06)+pow(logT,4)*(-2.49195655970986e+06)+pow(logT,5)*(8.47721959045642e+04)+pow(logT,6)*(2.00216126564460e+05)+pow(logT,7)*(-6.05650624482223e+04)+pow(logT,8)*(7.43527704221831e+03)+pow(logT,9)*(-3.47502441334837e+02);
				log_ro2=-4.82298469449122e+08+logT*(9.52166114274014e+08)+pow(logT,2)*(-7.72654088498897e+08)+pow(logT,3)*(3.14299000548327e+08)+pow(logT,4)*(-5.28908642170379e+07)+pow(logT,5)*(-7.01444850135997e+06)+pow(logT,6)*(5.30847769037309e+06)+pow(logT,7)*(-1.10646523640039e+06)+pow(logT,8)*(1.09487589753818e+05)+pow(logT,9)*(-4.37453072495502e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 && 3.40<logT && logT <3.86){
				log_ro1=8.61975403494922e+08+logT*(-1.66956631685583e+09)+pow(logT,2)*(1.32905124872802e+09)+pow(logT,3)*(-5.28874835642582e+08)+pow(logT,4)*(8.52085071754050e+07)+pow(logT,5)*(1.29543153081670e+07)+pow(logT,6)*(-9.01551689208385e+06)+pow(logT,7)*(1.84183097849026e+06)+pow(logT,8)*(-1.79987658782763e+05)+pow(logT,9)*(7.12407034707110e+03);
				log_ro2=1.76176628909089e+07+logT*(-1.72576525169566e+07)+pow(logT,2)*(1.10137006131467e+05)+pow(logT,3)*(5.88773417212743e+06)+pow(logT,4)*(-2.49195655970986e+06)+pow(logT,5)*(8.47721959045642e+04)+pow(logT,6)*(2.00216126564460e+05)+pow(logT,7)*(-6.05650624482223e+04)+pow(logT,8)*(7.43527704221831e+03)+pow(logT,9)*(-3.47502441334837e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 && 3.40<logT && logT <3.83){
				log_ro1=-4.05997960917400e+07+logT*(6.07518679414032e+07)+pow(logT,2)*(-3.26643571132997e+07)+pow(logT,3)*(5.56107197630373e+06)+pow(logT,4)*(1.20117723720358e+06)+pow(logT,5)*(-4.59953541749640e+05)+pow(logT,6)*(-5.27590762488530e+04)+pow(logT,7)*(4.23141941816351e+04)+pow(logT,8)*(-6.82818695171851e+03)+pow(logT,9)*(3.74305621273964e+02);
				log_ro2=8.61975403494922e+08+logT*(-1.66956631685583e+09)+pow(logT,2)*(1.32905124872802e+09)+pow(logT,3)*(-5.28874835642582e+08)+pow(logT,4)*(8.52085071754050e+07)+pow(logT,5)*(1.29543153081670e+07)+pow(logT,6)*(-9.01551689208385e+06)+pow(logT,7)*(1.84183097849026e+06)+pow(logT,8)*(-1.79987658782763e+05)+pow(logT,9)*(7.12407034707110e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5 && 3.40<logT && logT <3.778){
				log_ro1=1.48855492265646e+07+logT*(-2.13707584836559e+07)+pow(logT,2)*(1.09080369553653e+07)+pow(logT,3)*(-1.65596004756983e+06)+pow(logT,4)*(-4.33546592051429e+05)+pow(logT,5)*(1.40251464821950e+05)+pow(logT,6)*(2.03295321754383e+04)+pow(logT,7)*(-1.33896469365286e+04)+pow(logT,8)*(2.03698049155209e+03)+pow(logT,9)*(-1.06691742874476e+02);
				log_ro2=-4.05997960917400e+07+logT*(6.07518679414032e+07)+pow(logT,2)*(-3.26643571132997e+07)+pow(logT,3)*(5.56107197630373e+06)+pow(logT,4)*(1.20117723720358e+06)+pow(logT,5)*(-4.59953541749640e+05)+pow(logT,6)*(-5.27590762488530e+04)+pow(logT,7)*(4.23141941816351e+04)+pow(logT,8)*(-6.82818695171851e+03)+pow(logT,9)*(3.74305621273964e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0 && 3.40<logT && logT <3.70){
				log_ro1=3.59284481384146e+07+logT*(-5.18466906526573e+07)+pow(logT,2)*(2.64605459686339e+07)+pow(logT,3)*(-3.88554575133250e+06)+pow(logT,4)*(-1.13851695640503e+06)+pow(logT,5)*(3.49369730955857e+05)+pow(logT,6)*(5.73078838129681e+04)+pow(logT,7)*(-3.56019667464675e+04)+pow(logT,8)*(5.38472044844455e+03)+pow(logT,9)*(-2.81727628259288e+02);
				log_ro2=1.48855492265646e+07+logT*(-2.13707584836559e+07)+pow(logT,2)*(1.09080369553653e+07)+pow(logT,3)*(-1.65596004756983e+06)+pow(logT,4)*(-4.33546592051429e+05)+pow(logT,5)*(1.40251464821950e+05)+pow(logT,6)*(2.03295321754383e+04)+pow(logT,7)*(-1.33896469365286e+04)+pow(logT,8)*(2.03698049155209e+03)+pow(logT,9)*(-1.06691742874476e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else if ( 3.518<logT && logT <3.60){
				log_ro1=-1.82157829503303e+08+logT*(1.89164776040671e+08)+pow(logT,2)*(-4.96012513267803e+07)+pow(logT,3)*(-7.59139857648139e+06)+pow(logT,4)*(3.10866777935710e+06)+pow(logT,5)*(7.96492959570343e+05)+pow(logT,6)*(-1.83454306246182e+05)+pow(logT,7)*(-7.79549634299875e+04)+pow(logT,8)*(2.47803920657670e+04)+pow(logT,9)*(-1.91432694029750e+03);
				log_ro2=3.59284481384146e+07+logT*(-5.18466906526573e+07)+pow(logT,2)*(2.64605459686339e+07)+pow(logT,3)*(-3.88554575133250e+06)+pow(logT,4)*(-1.13851695640503e+06)+pow(logT,5)*(3.49369730955857e+05)+pow(logT,6)*(5.73078838129681e+04)+pow(logT,7)*(-3.56019667464675e+04)+pow(logT,8)*(5.38472044844455e+03)+pow(logT,9)*(-2.81727628259288e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=5.){
				log_ro1=8.35148e-08;
				log_ro2=1.86795e-06;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 ){ //same as above interpolations, but takes the lowest avalaible density for T outside of models range
				log_ro1=2.84313e-08;
				log_ro2=8.35148e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 ){
				log_ro1=1.26831e-08;
				log_ro2=2.84313e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 ){
				log_ro1=4.55076e-09;
				log_ro2=1.26831e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0){
				log_ro1=2.03554e-09;
				log_ro2=4.55076e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5){
				log_ro1=7.86729e-10;
				log_ro2=2.03554e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0){
				log_ro1=2.83936e-10;
				log_ro2=7.86729e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5){
				log_ro1=4.48868e-10;
				log_ro2=2.83936e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0){
				log_ro1=4.16804e-10;
				log_ro2=4.48868e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5){
				log_ro1=5.18036e-10;
				log_ro2=4.16804e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0){
				log_ro1=4.82163e-10;
				log_ro2=5.18036e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else {
				log_ro1=3.48774e-10;
				log_ro2=4.82163e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}

}


double ro_ph_z067218(double logg, double logT)
{ /* Calculate density in the photosphere for the z=0.067218 [M/H]=1*/
        /*density at kappa_ross=2/3 from MbRCS model atmospheres*/
		double log_ro1,log_ro2,log_ro;

		if (logg>=5.  && 3.38<logT && logT <3.59){
				log_ro1=6.97525868359929e+08+logT*(-1.38679714907556e+09)+pow(logT,2)*(1.13379985794280e+09)+pow(logT,3)*(-4.64947819318958e+08)+pow(logT,4)*(7.90147195135649e+07)+pow(logT,5)*(1.04605064787452e+07)+pow(logT,6)*(-8.02470417017503e+06)+pow(logT,7)*(1.68848093665442e+06)+pow(logT,8)*(-1.68632511866997e+05)+pow(logT,9)*(6.80090148212098e+03);
				log_ro2=-8.02820518407452e+09+logT*(1.20969062668606e+10)+pow(logT,2)*(-6.45793848124384e+09)+pow(logT,3)*(9.97167515194828e+08)+pow(logT,4)*(3.02395560456792e+08)+pow(logT,5)*(-9.81341361300054e+07)+pow(logT,6)*(-1.65749771310160e+07)+pow(logT,7)*(1.09414489004056e+07)+pow(logT,8)*(-1.74617374315698e+06)+pow(logT,9)*(9.64685646609090e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 && 3.38<logT && logT <3.89){
				log_ro1=-4.85003420774362e+08+logT*(9.06837002072860e+08)+pow(logT,2)*(-6.94746764762479e+08)+pow(logT,3)*(2.64589529350976e+08)+pow(logT,4)*(-3.98371545521604e+07)+pow(logT,5)*(-6.59871461804015e+06)+pow(logT,6)*(4.15125430491987e+06)+pow(logT,7)*(-8.03322600500821e+05)+pow(logT,8)*(7.45858728115719e+04)+pow(logT,9)*(-2.80262301904765e+03);
				log_ro2=6.97525868359929e+08+logT*(-1.38679714907556e+09)+pow(logT,2)*(1.13379985794280e+09)+pow(logT,3)*(-4.64947819318958e+08)+pow(logT,4)*(7.90147195135649e+07)+pow(logT,5)*(1.04605064787452e+07)+pow(logT,6)*(-8.02470417017503e+06)+pow(logT,7)*(1.68848093665442e+06)+pow(logT,8)*(-1.68632511866997e+05)+pow(logT,9)*(6.80090148212098e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 && 3.38<logT && logT <3.90){
				log_ro1=-9.03210637005769e+08+logT*(1.72570849375277e+09)+pow(logT,2)*(-1.35401037471350e+09)+pow(logT,3)*(5.30597841126428e+08)+pow(logT,4)*(-8.40720908778061e+07)+pow(logT,5)*(-1.26086813104775e+07)+pow(logT,6)*(8.60870611619805e+06)+pow(logT,7)*(-1.72617777926001e+06)+pow(logT,8)*(1.65448951142799e+05)+pow(logT,9)*(-6.41784757515762e+03);
				log_ro2=-4.85003420774362e+08+logT*(9.06837002072860e+08)+pow(logT,2)*(-6.94746764762479e+08)+pow(logT,3)*(2.64589529350976e+08)+pow(logT,4)*(-3.98371545521604e+07)+pow(logT,5)*(-6.59871461804015e+06)+pow(logT,6)*(4.15125430491987e+06)+pow(logT,7)*(-8.03322600500821e+05)+pow(logT,8)*(7.45858728115719e+04)+pow(logT,9)*(-2.80262301904765e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 && 3.38<logT && logT <3.90){
				log_ro1=-9.94031015235167e+08+logT*(1.91308682394595e+09)+pow(logT,2)*(-1.51259318855571e+09)+pow(logT,3)*(5.97954773234540e+08)+pow(logT,4)*(-9.61446737430369e+07)+pow(logT,5)*(-1.40033820975663e+07)+pow(logT,6)*(9.79628160433393e+06)+pow(logT,7)*(-1.98336447688464e+06)+pow(logT,8)*(1.91658312851268e+05)+pow(logT,9)*(-7.49235560115875e+03);
				log_ro2=-9.03210637005769e+08+logT*(1.72570849375277e+09)+pow(logT,2)*(-1.35401037471350e+09)+pow(logT,3)*(5.30597841126428e+08)+pow(logT,4)*(-8.40720908778061e+07)+pow(logT,5)*(-1.26086813104775e+07)+pow(logT,6)*(8.60870611619805e+06)+pow(logT,7)*(-1.72617777926001e+06)+pow(logT,8)*(1.65448951142799e+05)+pow(logT,9)*(-6.41784757515762e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0 && 3.38<logT && logT <3.90){
				log_ro1=-1.90757089120677e+08+logT*(3.78152796541935e+08)+pow(logT,2)*(-3.07753699386255e+08)+pow(logT,3)*(1.25416827329873e+08)+pow(logT,4)*(-2.11287453370408e+07)+pow(logT,5)*(-2.79456546566399e+06)+pow(logT,6)*(2.11475210305811e+06)+pow(logT,7)*(-4.39907409812815e+05)+pow(logT,8)*(4.34024536637506e+04)+pow(logT,9)*(-1.72761092152689e+03);
				log_ro2=-9.94031015235167e+08+logT*(1.91308682394595e+09)+pow(logT,2)*(-1.51259318855571e+09)+pow(logT,3)*(5.97954773234540e+08)+pow(logT,4)*(-9.61446737430369e+07)+pow(logT,5)*(-1.40033820975663e+07)+pow(logT,6)*(9.79628160433393e+06)+pow(logT,7)*(-1.98336447688464e+06)+pow(logT,8)*(1.91658312851268e+05)+pow(logT,9)*(-7.49235560115875e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5 && 3.38<logT && logT <3.90){
				log_ro1=4.58447261496220e+08+logT*(-8.75216977044250e+08)+pow(logT,2)*(6.86808977975841e+08)+pow(logT,3)*(-2.69514861852303e+08)+pow(logT,4)*(4.28948775015990e+07)+pow(logT,5)*(6.36296640630064e+06)+pow(logT,6)*(-4.38692129998118e+06)+pow(logT,7)*(8.84289431555026e+05)+pow(logT,8)*(-8.52362314734151e+04)+pow(logT,9)*(3.32769459579928e+03);
				log_ro2=-1.90757089120677e+08+logT*(3.78152796541935e+08)+pow(logT,2)*(-3.07753699386255e+08)+pow(logT,3)*(1.25416827329873e+08)+pow(logT,4)*(-2.11287453370408e+07)+pow(logT,5)*(-2.79456546566399e+06)+pow(logT,6)*(2.11475210305811e+06)+pow(logT,7)*(-4.39907409812815e+05)+pow(logT,8)*(4.34024536637506e+04)+pow(logT,9)*(-1.72761092152689e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0 && 3.38<logT && logT <3.89){
				log_ro1=9.63085446845080e+08+logT*(-1.85994301985769e+09)+pow(logT,2)*(1.47665108712451e+09)+pow(logT,3)*(-5.86645306229308e+08)+pow(logT,4)*(9.49797114868869e+07)+pow(logT,5)*(1.37943747012150e+07)+pow(logT,6)*(-9.75122868651731e+06)+pow(logT,7)*(1.98911713148693e+06)+pow(logT,8)*(-1.93707254656269e+05)+pow(logT,9)*(7.63496114464734e+03);
				log_ro2=4.58447261496220e+08+logT*(-8.75216977044250e+08)+pow(logT,2)*(6.86808977975841e+08)+pow(logT,3)*(-2.69514861852303e+08)+pow(logT,4)*(4.28948775015990e+07)+pow(logT,5)*(6.36296640630064e+06)+pow(logT,6)*(-4.38692129998118e+06)+pow(logT,7)*(8.84289431555026e+05)+pow(logT,8)*(-8.52362314734151e+04)+pow(logT,9)*(3.32769459579928e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5 && 3.38<logT && logT <3.845){
				log_ro1=3.95081670363895e+09+logT*(-7.71392970535332e+09)+pow(logT,2)*(6.18817535493346e+09)+pow(logT,3)*(-2.48059019855422e+09)+pow(logT,4)*(4.02215762413971e+08)+pow(logT,5)*(6.18250158810681e+07)+pow(logT,6)*(-4.32332225088074e+07)+pow(logT,7)*(8.88642549485193e+06)+pow(logT,8)*(-8.73605375512091e+05)+pow(logT,9)*(3.47767018387884e+04);
				log_ro2=9.63085446845080e+08+logT*(-1.85994301985769e+09)+pow(logT,2)*(1.47665108712451e+09)+pow(logT,3)*(-5.86645306229308e+08)+pow(logT,4)*(9.49797114868869e+07)+pow(logT,5)*(1.37943747012150e+07)+pow(logT,6)*(-9.75122868651731e+06)+pow(logT,7)*(1.98911713148693e+06)+pow(logT,8)*(-1.93707254656269e+05)+pow(logT,9)*(7.63496114464734e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0 && 3.38<logT && logT <3.83){
				log_ro1=-7.99462958843708e+07+logT*(1.18013089997867e+08)+pow(logT,2)*(-6.24668054175260e+07)+pow(logT,3)*(1.03390753658824e+07)+pow(logT,4)*(2.32037173347018e+06)+pow(logT,5)*(-8.52677684562722e+05)+pow(logT,6)*(-1.02728744042506e+05)+pow(logT,7)*(7.84318661720596e+04)+pow(logT,8)*(-1.24520768292538e+04)+pow(logT,9)*(6.73665383596138e+02);
				log_ro2=3.95081670363895e+09+logT*(-7.71392970535332e+09)+pow(logT,2)*(6.18817535493346e+09)+pow(logT,3)*(-2.48059019855422e+09)+pow(logT,4)*(4.02215762413971e+08)+pow(logT,5)*(6.18250158810681e+07)+pow(logT,6)*(-4.32332225088074e+07)+pow(logT,7)*(8.88642549485193e+06)+pow(logT,8)*(-8.73605375512091e+05)+pow(logT,9)*(3.47767018387884e+04);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5 && 3.38<logT && logT <3.78){
				log_ro1=-1.45767614632323e+07+logT*(2.17322301582162e+07)+pow(logT,2)*(-1.15760266509880e+07)+pow(logT,3)*(1.89151248613568e+06)+pow(logT,4)*(4.59740317723305e+05)+pow(logT,5)*(-1.63663271706492e+05)+pow(logT,6)*(-2.17021299778072e+04)+pow(logT,7)*(1.58529992945062e+04)+pow(logT,8)*(-2.51715185961659e+03)+pow(logT,9)*(1.36663658567866e+02);
				log_ro2=-7.99462958843708e+07+logT*(1.18013089997867e+08)+pow(logT,2)*(-6.24668054175260e+07)+pow(logT,3)*(1.03390753658824e+07)+pow(logT,4)*(2.32037173347018e+06)+pow(logT,5)*(-8.52677684562722e+05)+pow(logT,6)*(-1.02728744042506e+05)+pow(logT,7)*(7.84318661720596e+04)+pow(logT,8)*(-1.24520768292538e+04)+pow(logT,9)*(6.73665383596138e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0 && 3.38<logT && logT <3.70){
				log_ro1=-5.37452077143648e+08+logT*(7.86162296584933e+08)+pow(logT,2)*(-4.07814242259562e+08)+pow(logT,3)*(6.16568603860497e+07)+pow(logT,4)*(1.76794894417029e+07)+pow(logT,5)*(-5.63427573246652e+06)+pow(logT,6)*(-9.00638916422784e+05)+pow(logT,7)*(5.83381466313878e+05)+pow(logT,8)*(-9.03644414916000e+04)+pow(logT,9)*(4.83840940034741e+03);
				log_ro2=-1.45767614632323e+07+logT*(2.17322301582162e+07)+pow(logT,2)*(-1.15760266509880e+07)+pow(logT,3)*(1.89151248613568e+06)+pow(logT,4)*(4.59740317723305e+05)+pow(logT,5)*(-1.63663271706492e+05)+pow(logT,6)*(-2.17021299778072e+04)+pow(logT,7)*(1.58529992945062e+04)+pow(logT,8)*(-2.51715185961659e+03)+pow(logT,9)*(1.36663658567866e+02);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else if ( 3.38<logT && logT <3.60){
				log_ro1=1.59330676303301e+07+logT*(-1.64596752656247e+07)+pow(logT,2)*(4.28367834882511e+06)+pow(logT,3)*(6.60614282008195e+05)+pow(logT,4)*(-2.65105470487017e+05)+pow(logT,5)*(-6.96246445631172e+04)+pow(logT,6)*(1.56359201141847e+04)+pow(logT,7)*(6.78094837961970e+03)+pow(logT,8)*(-2.13581087974401e+03)+pow(logT,9)*(1.64334492989189e+02);
				log_ro2=-5.37452077143648e+08+logT*(7.86162296584933e+08)+pow(logT,2)*(-4.07814242259562e+08)+pow(logT,3)*(6.16568603860497e+07)+pow(logT,4)*(1.76794894417029e+07)+pow(logT,5)*(-5.63427573246652e+06)+pow(logT,6)*(-9.00638916422784e+05)+pow(logT,7)*(5.83381466313878e+05)+pow(logT,8)*(-9.03644414916000e+04)+pow(logT,9)*(4.83840940034741e+03);
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=5.){
				log_ro1=7.71139e-08;
				log_ro2=1.47222e-06;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-5.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.5 ){ //same as above interpolations, but takes the lowest avalaible density for T outside of models range
				log_ro1=2.67054e-08;
				log_ro2=7.71139e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=4.0 ){
				log_ro1=1.20411e-08;
				log_ro2=2.67054e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-4.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.5 ){
				log_ro1=4.26282e-09;
				log_ro2=1.20411e-08;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=3.0){
				log_ro1=1.93682e-09;
				log_ro2=4.26282e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-3.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.5){
				log_ro1=7.50165e-10;
				log_ro2=1.93682e-09;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=2.0){
				log_ro1=3.81273e-10;
				log_ro2=7.50165e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-2.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.5){
				log_ro1=6.25041e-10;
				log_ro2=3.81273e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=1.0){
				log_ro1=3.16787e-10;
				log_ro2=6.25041e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-1.)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.5){
				log_ro1=5.62915e-10;
				log_ro2=3.16787e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.5)/0.5;
				return pow(10,log_ro);
		}
		else if (logg>=0.0){
				log_ro1=3.9606e-10;
				log_ro2=5.62915e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg-0.)/0.5;
				return pow(10,log_ro);
		}
		else {
				log_ro1=2.60952e-10;
				log_ro2=3.9606e-10;
				log_ro=log_ro1+(log_ro2-log_ro1)*(logg+0.5)/0.5;
				return pow(10,log_ro);
		}

}


void evroute_add(char *format, ...)
{ /* adds input to the evolutionary route the same as for printf-like functions */
 char tmp[32];
 va_list ap;

 va_start(ap, format);
 vsprintf(tmp, format, ap);
 strcat(evroute, tmp);
 if (tmp[strlen(tmp)-1]==' ') nevroute++;
}

void evroute_new()
{ /* initializes the evroute variable used for storing the evolutionary route */
  /* should be run at the beginning of the system's evolution */

 evroute[0]='\0';
 nevroute=0;
}



#undef G
#undef GG
#undef GGG
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
#undef JMAX
#undef ITMAX
#undef EPS
#undef FPMIN
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef EPS2
#undef NTAB
#undef NDIV
#undef RNMX
#undef EPS3
#undef JMAX3
#undef FUNC
#undef FUNC1
#undef CON
#undef CON2
#undef BIG1
#undef NTAB1
#undef SAFE1
#undef MAXSTP
#undef MAXSTP1
#undef TINY
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef NR_END
#undef FREE_ARG
#undef NRANSI
