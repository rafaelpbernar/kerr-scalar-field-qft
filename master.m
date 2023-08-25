(* ::Package:: *)

(* ::Title:: *)
(*BWBC1 - Master notebook*)


(* ::Chapter:: *)
(*Initializing...*)


(* ::Input:: *)
(*Quit[]*)


(* ::Text:: *)
(*We load a few convenient stuff that will be used throughout the whole notebook.\[AliasDelimiter]*)


(* ::Input:: *)
(*<< SpinWeightedSpheroidalHarmonics` *)


(* ::Text:: *)
(*Some assumptions about our coordinates and labels. We also do some hacking to make the Derivative operator commute with the Conjugate operator.*)


(* ::Input:: *)
(*$Assumptions={r>=0,ms >= 0,ms^2>a^2, Element[L, Integers], Element[M, Integers], Element[ \[Omega], Reals],Element[{t,r,\[Theta],\[Phi]}, Reals]};*)
(*excluded="ExcludedFunctions"/.("DifferentiationOptions"/.SystemOptions["DifferentiationOptions"]);*)
(*SetSystemOptions["DifferentiationOptions"->"ExcludedFunctions"->Union[excluded,{Conjugate}]];*)
(*Unprotect[Conjugate];*)
(*Conjugate/:D[Conjugate[f_],x__]:=Conjugate[D[f,x]];*)
(*Protect[Conjugate];*)


(* ::Text:: *)
(*We initialize some definitions such as the radial f(r) function, the normalization constants for the field modes and some rules.*)


(* ::Input:: *)
(*aext=ms;*)
(*\[CapitalDelta][r_]:=((r-rp) (r-rm));*)
(*f[r_]:= \[CapitalDelta][r]/(r^2+a^2);*)
(*OmegaH = a/(2 ms rp);*)
(*\[Omega]rule=\[Omega]tilde->\[Omega]-m  OmegaH;*)
(*msrule= ms -> (rp+rm)/2;*)
(*arule=a -> Sqrt[rp]Sqrt[rm];*)
(*rrule = {rp -> ms + Sqrt[ms^2-a^2], rm -> ms - Sqrt[ms^2-a^2]};*)
(*\[Kappa]=Derivative[1][f][rp]/2;*)
(*A["in",\[Omega]_,L_,M_]:=1/Sqrt[4 \[Pi] Abs[\[Omega]]];A["up",\[Omega]_,L_,M_]:=1/Sqrt[4 \[Pi] Abs[\[Omega]-M OmegaH]];*)


(* ::Chapter:: *)
(*Scalar field quantities*)


(* ::Text:: *)
(*In this part we use a set of modes in the frequency domain to compute a few quantities associated to the scalar field observables. Most are contributions from each mode type to the expectation values of several operators in a given state. *)


(* ::Text:: *)
(*We define the prefix for the directory to save us from typing it everywhere.*)


(* ::Input:: *)
(*DirPrefix = NotebookDirectory[]*)


(* ::Text:: *)
(*We load the set of parameters, taken from the auxiliary files in the folders where the modes are stored. We also check if the parameters loaded are consistent with each other.*)


(* ::Input:: *)
(*Parameters1=Import[DirPrefix<>"/a_31_ID-2023-Mar13-173120/parameters_ID-2023-Mar13-173120.txt","Data"];*)
(*Parameters2=Import[DirPrefix<>"/a_31_ID-2023-Mar13-173120/parameters_ID-2023-Mar13-173120.txt","Data"];*)
(*If[Parameters1[[3;;Length[Parameters1]-2]] != Parameters2[[3;;Length[Parameters1]-2]],Print["WARNING! Background parameters are not equal."]];*)
(*If[Parameters1[[1]]==Parameters2[[1]],Print["WARNING! You are loading two sets of the same mode type ("<>Parameters2[[1,2]]<>")."]];*)
(*Print[Parameters1];*)
(*Print[Parameters2];*)


(* ::Text:: *)
(*From the parameters, we initialize the values for each background quantity, together with some details from each set of modes.*)


(* ::Input:: *)
(*rp=Parameters1[[3,2]];*)
(*rm=Parameters1[[4,2]];*)
(*ms=Parameters1[[5,2]];*)
(*a=Rationalize[Parameters1[[6,2]]];*)
(*rinitial=Parameters1[[10,2]];*)
(*Lmin=IntegerPart[Parameters1[[11,2]]];*)
(*Lmax=IntegerPart[Parameters1[[12,2]]];*)
(*name1=Parameters1[[1,2]];*)
(*buildtime1=Parameters1[[13,2]];*)
(*name2=Parameters2[[1,2]];*)
(*buildtime2=Parameters2[[13,2]];*)


(* ::Text:: *)
(*We then add the folders where the modes are stored to the current path to pull them when required.*)


(* ::Input:: *)
(*dir1=DirPrefix<>"a_"<>ToString[IntegerPart[100*a]]<>"_"<>buildtime1;*)
(*dir2=DirPrefix<>"a_"<>ToString[IntegerPart[100*a]]<>"_"<>buildtime2;*)
(*AppendTo[$Path,dir1];*)
(*AppendTo[$Path,dir2];*)


(* ::Text:: *)
(*We initialize the list of frequencies available for the computed modes. There is also a coarser discretization defined for test purposes.*)


(* ::Input:: *)
(*\[Omega]list1=Import["omegalist_"<>buildtime1<>".mx"];*)
(*\[Omega]list2=Import["omegalist_"<>buildtime2<>".mx"];*)
(*If[\[Omega]list1!= \[Omega]list2,Print["WARNING! List of frequencies is not equal."]];*)
(*\[Omega]list=\[Omega]list1;*)
(*\[Omega]listalt=\[Omega]list[[1;;;;100]];*)


(* ::Text:: *)
(*We define some tables where we store the transmission and reflection coefficients. We also need to do some extra word because I forgot to define these tables with the frequency values included.*)


(* ::Input:: *)
(*Do[*)
(*TransmissionTable[L]=Import["in_L_"<>ToString[L]<>"_transmission_coefficient.mx"];*)
(*TRA[L]=Table[{\[Omega]list[[i]],TransmissionTable[L][[i]]},{i,1,Length[\[Omega]list]}];*)
(*ReflectionTable[L]=Import["in_L_"<>ToString[L]<>"_reflection_coefficient.mx"];*)
(*REF[L]=Table[{\[Omega]list[[i]],ReflectionTable[L][[i]]},{i,1,Length[\[Omega]list]}];*)
(*,{L,Lmin,Lmax,1}];*)


(* ::Text:: *)
(*For convenience, we also defined some functions from the above tables, obtained through interpolation.*)


(* ::Input:: *)
(*InterpolatedTransmission[L_]:=Interpolation[TRA[L]]*)
(*InterpolatedReflection[L_]:=Interpolation[REF[L]]*)


(* ::Text:: *)
(*We store the modes in memory. This leads to way faster integration.*)


(* ::Input:: *)
(*AbsoluteTiming[*)
(*Do[*)
(*(*StoredMode["up",\[Omega],L,M]=Import["up"<>"_L_"<>ToString[L]<>"_m_"<>ToString[]<>"_omega_"<>ToString[\[Omega] // N]<>".mx"];*)*)
(*StoredMode["in",\[Omega],L,M]=Import["in"<>"_L_"<>ToString[L]<>"_m_"<>ToString[M]<>"_omega_"<>ToString[\[Omega] // N]<>".mx"];*)
(*,{\[Omega],\[Omega]list},{L,Lmin,Lmax,1},{M,-L,L}*)
(*]*)
(*]*)


(* ::Text:: *)
(*Now the magic happens and the modes, previously computed, are loaded in. We define the modes in terms of the loaded data in memory.*)


(* ::Input:: *)
(*\[Psi][mode_,\[Omega]_,L_,M_][x_]:=( *)
(*Which[mode == "up", *)
(*StoredMode["up",\[Omega],L,M] /. r -> x,*)
(*mode=="in",*)
(*StoredMode["in",\[Omega],L,M] /. r -> x,*)
(*mode=="down",*)
(*Conjugate[StoredMode["up",\[Omega],L,M] /. r -> x],*)
(*mode == "out",*)
(*Conjugate[StoredMode["in",\[Omega],L,M] /. r -> x],*)
(*True,*)
(*Return[$Failed];*)
(*Print["No modes defined with this name"];*)
(*])*)


(* ::Section::Closed:: *)
(*Legacy code*)


(* ::Text:: *)
(*This is a terrible way to load the modes. I am leaving it here so that I never forget.*)


(* ::Input:: *)
(*\[Psi][mode_,\[Omega]_,L_][x_]:=( *)
(*Which[mode == "up", *)
(*Import[mode<>"_L_"<>ToString[L]<>"_omega_"<>ToString[\[Omega] // N]<>".mx"] /. r -> x,*)
(*mode=="in",*)
(*Import[mode<>"_L_"<>ToString[L]<>"_omega_"<>ToString[\[Omega] // N]<>".mx"] /. r -> x,*)
(*mode=="down",*)
(*Conjugate[Import["up"<>"_L_"<>ToString[L]<>"_omega_"<>ToString[\[Omega] // N]<>".mx"] /. r -> x],*)
(*mode == "out",*)
(*Conjugate[Import["in"<>"_L_"<>ToString[L]<>"_omega_"<>ToString[\[Omega] // N]<>".mx"] /. r -> x],*)
(*True,*)
(*Return[$Failed];*)
(*Print["No modes defined with this name"];*)
(*])*)


(* ::Subchapter:: *)
(*Scalar modes contributions*)


(* ::Text:: *)
(*Each scalar field mode contribution to the relevant quantities. Note that we have removed the angular dependency, hence there is no m-dependent part. Even more important is the fact that there is no (2L + 1)/(4 Pi)in front of those quantities, which should come from the m summation.*)


(* ::Input:: *)
(*PhiSquared[mode_,\[Omega]0_,L_,M_]:=Abs[A[mode,\[Omega]0,L,M]]^2 Abs[\[Psi][mode,\[Omega]0,L,M][r]/Sqrt[r^2+a^2]]^2 *)
(*Currenttcomponent[mode_,\[Omega]0_,L_]:=-((q (\[Omega]0-(q Q)/r) Abs[A[mode,\[Omega]0,L]]^2 Abs[\[Psi][mode,\[Omega]0,L][r]/r]^2)/(4 \[Pi] f[r]))*)
(**)
(*Currentrcomponent[mode_,\[Omega]0_,L_]:=-((q f[r])/(4 \[Pi])) Abs[A[mode,\[Omega]0,L]]^2 Im[Conjugate[\[Psi][mode,\[Omega]0,L][r]]/r D[\[Psi][mode,\[Omega]0,L][r]/r,r]]  *)
(**)
(*Stressttcomponent[mode_,\[Omega]0_,L_]:=1/2 ((\[Omega]0-(q Q)/rp)^2+(L (L+1) f[r])/r^2)Abs[A[mode,\[Omega]0,L]]^2 Abs[\[Psi][mode,\[Omega]0,L][r]/r]^2+f[r]^2/2 Abs[A[mode,\[Omega]0,L]]^2 Abs[D[\[Psi][mode,\[Omega]0,L][r]/r,r]]^2 *)
(**)
(**)
(*Stresstrcomponent[mode_,\[Omega]0_,L_]:=-(\[Omega]0-(q Q)/r) Abs[A[mode,\[Omega]0,L]]^2 Im[Conjugate[\[Psi][mode,\[Omega]0,L][r]]/r D[\[Psi][mode,\[Omega]0,L][r]/r,r]]*)
(**)
(*Stressrsupertcomponent[mode_,\[Omega]0_,L_]:=-(\[Omega]0-(q Q)/r)f[r]Abs[A[mode,\[Omega]0,L]]^2 Im[Conjugate[\[Psi][mode,\[Omega]0,L][r]]/r D[\[Psi][mode,\[Omega]0,L][r]/r,r]]*)
(**)
(*Stressrrcomponent[mode_,\[Omega]0_,L_]:=1/2  Abs[A[mode,\[Omega]0,L]]^2 Abs[D[\[Psi][mode,\[Omega]0,L][r]/r,r]]^2+1/(2 f[r]) ((\[Omega]0-(q Q)/rp)^2/f[r]-(L (L+1))/r^2)Abs[A[mode,\[Omega]0,L]]^2 Abs[\[Psi][mode,\[Omega]0,L][r]/r]^2 *)
(*Stressthetathetacomponent[mode_,\[Omega]0_,L_]:=r^2/(2 f[r]) (\[Omega]0-(q Q)/rp)^2 Abs[A[mode,\[Omega]0,L]]^2 Abs[\[Psi][mode,\[Omega]0,L][r]/r]^2-(1/2) (r^2) f[r](Abs[A[mode,\[Omega]0,L]]^2) (Abs[D[\[Psi][mode,\[Omega]0,L][r]/r,r]]^2) *)


(* ::Section:: *)
(*Some examples*)


(* ::Text:: *)
(*The value of |\[CapitalPhi](|^2) for a mode contribution.*)


(* ::Input:: *)
(*Plot[PhiSquared["in",1,0,0] /. r -> x,{x,rinitial,5} ,PlotRange->All]*)


(* ::Text:: *)
(*We may plot the contribution to the current component J^r for each mode. This, multiplied by r^2, should be a conserved quantity.*)


(* ::Input:: *)
(*Plot[r^2 Currentrcomponent["in",1,0] /. r -> x,{x,rinitial,5} ,PlotRange->All]*)


(* ::Input:: *)
(*Plot[r^2 Currentrcomponent["up",1,0] /. r -> x,{x,rinitial,5} ,PlotRange->All]*)


(* ::Subchapter::Closed:: *)
(*Legacy quantities*)


(* ::Text:: *)
(*Quantities that are not needed -- they should be zero. It might be useful for testing purposes though.*)


(* ::Input:: *)
(*Currentthetacomponent[mode_,\[Omega]0_,L_,m_]:=( *)
(*-(q /(4 Pi (r^4) ))Abs[A[mode,\[Omega]0,L]]^2Abs[\[Psi][mode,\[Omega]0,L][r]]^2  Im[Conjugate[SphericalHarmonicY[L,m,\[Theta],\[Phi]]]D[SphericalHarmonicY[L,m,\[Theta],\[Phi]],\[Theta]]] (*// Simplify*)*)
(*)*)
(*Currentphicomponent[mode_,\[Omega]0_,L_(*,m_*)]:=( *)
(*-((m q )/(4 Pi (r^2) (Sin[\[Theta]]^2) ))Abs[A[mode,\[Omega]0,L]]^2 Abs[\[Psi][mode,\[Omega]0,L][r]/r]^2 (*Abs[SphericalHarmonicY[L,m,\[Theta],\[Phi]]]^2*) (*// Simplify*)*)
(*)*)
(*Stresstphicomponent[mode_,\[Omega]0_,L_(*,m_*)]:=( *)
(*- m (\[Omega]-(q Q)/r)Abs[A[mode,\[Omega]0,L]]^2Abs[\[Psi][mode,\[Omega]0,L][r]/r]^2 (*Abs[SphericalHarmonicY[L,m,\[Theta],\[Phi]]]^2*)  (*// Simplify*)*)
(*)*)
(*Stressrthetacomponent[mode_,\[Omega]0_,L_,m_]:= Abs[A[mode,\[Omega]0,L]]^2 Re[\[Psi][mode,\[Omega]0,L][r]/r Conjugate[D[\[Psi][mode,\[Omega]0,L][r]/r,r]]Conjugate[SphericalHarmonicY[L,m,\[Theta],\[Phi]]]D[SphericalHarmonicY[L,m,\[Theta],\[Phi]],\[Theta]]] *)
(**)
(*Stressrphicomponent[mode_,\[Omega]0_,L_(*,m_*)]:=( *)
(*-(m /r) Abs[A[mode,\[Omega]0,L]]^2 Im[\[Psi][mode,\[Omega]0,L][r]/r Conjugate[D[\[Psi][mode,\[Omega]0,L][r]/r,r]]](*Abs[SphericalHarmonicY[L,m,\[Theta],\[Phi]]]^2*)  (*// Simplify*)*)
(*)*)
(*Stressthetaphicomponent[mode_,\[Omega]0_,L_,m_]:=( *)
(*-(m/r^2) Abs[A[mode,\[Omega]0,L]]^2 Abs[\[Psi][mode,\[Omega]0,L][r]]^2 Im[SphericalHarmonicY[L,m,\[Theta],\[Phi]]Conjugate[D[SphericalHarmonicY[L,m,\[Theta],\[Phi]],\[Theta]]]]  (*// Simplify*)*)
(*)*)
(*Stressphiphicomponent[mode_,\[Omega]0_,L_]:=*)
(*(r^2 Sin[\[Theta]]^2)/(2 f[r]) (\[Omega]0-(q Q)/rp)^2 Abs[A[mode,\[Omega]0,L]]^2 Abs[\[Psi][mode,\[Omega]0,L][r]/r]^2-1/2 f[r] r^2 Sin[\[Theta]]^2  Abs[A[mode,\[Omega]0,L]]^2 Abs[D[\[Psi][mode,\[Omega]0,L][r]/r,r]]^2*)
(**)


(* ::Subchapter:: *)
(*Integration*)


(* ::Text:: *)
(*This hacky code defines the relevant quantity to be integrated, according to which states we choose in our subtraction scheme. We note that, besides choosing the correct relevant quantity, we must also define the integration limits to compute the correct values of a given difference of expectation values. *)


(* ::Input:: *)
(*RelevantQuantitygenerator[name0_,quantitystring_]:=( *)
(*Which[quantitystring == "PhiSquared",string="condensate",quantitystring == "Currenttcomponent",string="current_t",quantitystring == "Currentrcomponent",string="current_r",quantitystring == "Stressttcomponent",string="Stress_tt",quantitystring == "Stresstrcomponent",string="Stress_tr",quantitystring == "Stressrsupertcomponent",string="Stress_rsupert",quantitystring == "Stresstphicomponent",string="Stress_tphi",quantitystring == "Stressrrcomponent",string="Stress_rr",quantitystring == "Stressrthetacomponent",string="Stress_rtheta",quantitystring == "Stressrphicomponent",string="Stress_rphi",quantitystring == "Stressthetathetacomponent",string="Stress_thetatheta",quantitystring == "Stressthetaphicomponent",string="Stress_thetaphi",quantitystring == "Stressphiphicomponent",string="Stress_phiphi"];*)
(*name=name0<>"_"<>string;*)
(*quantity=ToExpression[quantitystring];*)
(*Which[name0=="pUnruh_m_pB",*)
(*RelevantQuantity[\[Omega]_,L_]:=quantity["up",\[Omega],L]/(Exp[Abs[(2 \[Pi] (\[Omega]-(q Q)/rp))/\[Kappa]]]-1),*)
(*name0=="pB",*)
(*RelevantQuantity[\[Omega]_,L_]:=1/2 (quantity["in",\[Omega],L]+quantity["up",\[Omega],L]),*)
(*name0=="pB_m_fB",*)
(*RelevantQuantity[\[Omega]_,L_]:=1/2 (quantity["in",\[Omega],L]-quantity["out",\[Omega],L]+quantity["up",\[Omega],L]-quantity["down",\[Omega],L]),*)
(*name0=="pCCH_m_pUnruh",*)
(*RelevantQuantity[\[Omega]_,L_,M_]:=2 quantity["in",\[Omega],L,M]/(Exp[Abs[(2 \[Pi] \[Omega])/\[Kappa]]]-1),*)
(*name0=="fUnruh_m_fB",*)
(*RelevantQuantity[\[Omega]_,L_]:=quantity["down",\[Omega],L]/(Exp[Abs[(2 \[Pi] (\[Omega]-(q Q)/rp))/\[Kappa]]]-1),*)
(*name0=="fCCH_m_fUnruh",*)
(*RelevantQuantity[\[Omega]_,L_]:=quantity["out",\[Omega],L]/(Exp[Abs[(2 \[Pi] \[Omega])/\[Kappa]]]-1),*)
(*name0=="pCCH_m_pB", (*Legacy code*)RelevantQuantity[\[Omega]_,L_]:=quantity["up",\[Omega],L]/(Exp[Abs[2 \[Pi] (\[Omega]-(q Q)/rp)]]-1)+quantity["in",\[Omega],L]/(Exp[Abs[2 \[Pi] \[Omega]]]-1),*)
(*name0 == "B_m_pB",*)
(*RelevantQuantity[\[Omega]_,L_]:=-quantity["up",\[Omega],L],*)
(*name0 == "B_m_fB",*)
(*RelevantQuantity[\[Omega]_,L_]:=-quantity["down",\[Omega],L],*)
(*name0=="FT_m_pUnruh",*)
(*RelevantQuantity[\[Omega]_,L_]:=quantity["in",\[Omega],L]/(Exp[Abs[(2 \[Pi] (\[Omega]-(q Q)/rp))/\[Kappa]]]-1),*)
(*name0=="H_m_FT",*)
(*RelevantQuantity[\[Omega]_,L_]:=-quantity["in",\[Omega],L] Coth[Abs[(\[Pi] (\[Omega]-(q Q)/rp))/\[Kappa]]],*)
(*name0 == "H_m_pUnruh", *)
(*RelevantQuantity[\[Omega]_,L_]:= quantity["in",\[Omega],L]/(Exp[(2 Pi  (\[Omega]-(q Q)/rp))/\[Kappa]]-1)+quantity["in",-\[Omega],L]/(Exp[(2 Pi  (\[Omega]+(q Q)/rp))/\[Kappa]]-1),*)
(*name0 == "HH_m_pB", (*Legacy code*)RelevantQuantity[\[Omega]_,L_]:=quantity["up",\[Omega],L]/(Exp[Abs[2 \[Pi] (\[Omega]-(q Q)/rp)]]-1)+quantity["in",\[Omega],L]/(Exp[Abs[2 \[Pi] (\[Omega]-(q Q)/rp)]]-1),True,*)
(*Print["Undefined vacuum states names"];Return[$Failed]]*)
(*)*)


(* ::Text:: *)
(*We make below a concrete choice of subtraction scheme and quantity that we want to integrate.*)


(* ::Input:: *)
(*RelevantQuantitygenerator["pCCH_m_pUnruh","PhiSquared"];*)
(*Print["You are currently computing quantity defined as "<>name<>" for a="<>ToString[a]<>" rp."];*)


(* ::Text:: *)
(*In case we need to exclude some points.*)


(* ::Input:: *)
(*improved\[Omega]list = Select[\[Omega]list,N[#] != -((q Q)/rp)&];*)


(* ::Text:: *)
(*We start the integration process by creating a list of the relevant quantity as a function of the frequency. The spheroidal harmonics are computed at the equatorial plane and \[Phi]=0.*)


(* ::Input:: *)
(*Print["Generating integrand..."];*)
(*AbsoluteTiming[Do[Integrand[L,M]=Table[{\[Omega], Abs[SpinWeightedSpheroidalHarmonicS[0,L,M, SetPrecision[a \[Omega],20],Pi/2,0]]^2 RelevantQuantity[\[Omega],L,M]},{\[Omega],\[Omega]list}];,{L,Lmin,Lmax,1},{M,-L,L}]]*)


(* ::Text:: *)
(*To better capture the picture near the black hole, we pick more points near the initial radius. For that we create a function that has a higher density of points near the horizon.*)


(* ::Input:: *)
(*radius[n_]:=rinitial+((n-1)/99)^3*(10-rinitial);*)
(*rlist=radius[Range[100]];*)
(*ListPlot[rlist]*)


(* ::Text:: *)
(*We pick the integration limits here. It is hardcoded, currently.*)


(* ::Input:: *)
(*Print["Integrating..."];*)
(*Integral={};*)
(*AbsoluteTiming[*)
(*Do[*)
(*IntegrandSum=Table[{Integrand[0,0][[i,1]],Sum[ Integrand[L,M][[i,2]],{L,Lmin,Lmax,1},{M,-L,L}]},{i,1,Length[Integrand[0,0]]}] /. r-> k;*)
(*IntegrandSuminterpolated=Interpolation[IntegrandSum];*)
(*AppendTo[Integral,{k,NIntegrate[IntegrandSuminterpolated[\[Omega]],{\[Omega],0,0,2},MaxRecursion->25,WorkingPrecision->50,AccuracyGoal->14]}];*)
(*,{k,rlist}];]*)


(* ::Text:: *)
(*We test with the scalar condensate:*)


(* ::Input:: *)
(*ListPlot[Integral,Joined->True,PlotRange->All]*)


(* ::Text:: *)
(*Alternative integration method. There is no need to generate the integrand before this.*)


(* ::Input:: *)
(*Tableo=Table[{\[Omega],Sum[Abs[SpinWeightedSpheroidalHarmonicS[0,L,M, SetPrecision[a \[Omega],20],Pi/2,0]]^2 RelevantQuantity[\[Omega],L,M],{L,Lmin,Lmax},{M,-L,L}]  },{\[Omega],\[Omega]list}];*)


(* ::Input:: *)
(*Print["Integrating..."];*)
(*Integral={};*)
(*AbsoluteTiming[*)
(*Do[*)
(*Integrand[x_]:=Interpolation[Tableo /. r-> k][x];*)
(*solution=NDSolve[{Int'[x]==Integrand[x],Int[\[Omega]list[[1]]]==0},Int[x],{x,\[Omega]list[[1]],4}][[1]];*)
(*AppendTo[Integral,{k,Int[x] /. solution /. x -> 4}];*)
(*,{k,rlist}];]*)


(* ::Input:: *)
(*ListPlot[Integral,Joined->True,PlotRange->All]*)


(* ::Subchapter::Closed:: *)
(*Asymptotic quantities*)


(* ::Input:: *)
(*\[Omega]i = \[Omega]list[[1]];*)
(*\[Omega]f = \[Omega]list[[-1]];*)


(* ::Section:: *)
(*Past Unruh minus past Boulware*)


(* ::Text:: *)
(*Scalar condensate*)


(* ::Input:: *)
(*Sum[(2L+1)NIntegrate[1/(Abs[x- q Q](Exp[Abs[(2 Pi (x -q Q))/\[Kappa]]]-1)) Abs[InterpolatedTransmission[L][x]]^2,{x,\[Omega]i,\[Omega]f}],{L,0,10,1}]*)


(* ::Text:: *)
(*Current operator*)


(* ::Input:: *)
(*Sum[(2L+1)NIntegrate[x/(Abs[x- q Q](Exp[Abs[(2 Pi (x -q Q))/\[Kappa]]]-1)) Abs[InterpolatedTransmission[L][x]]^2,{x,\[Omega]i,\[Omega]f}],{L,0,10,1}]*)


(* ::Text:: *)
(*Stress-energy tensor*)


(* ::Input:: *)
(*Sum[(2L+1)NIntegrate[x^2/(Abs[x- q Q](Exp[Abs[(2 Pi (x -q Q))/\[Kappa]]]-1)) Abs[InterpolatedTransmission[L][x]]^2,{x,\[Omega]i,\[Omega]f}],{L,0,10,1}]*)


(* ::Section:: *)
(*Past Unruh*)


(* ::Text:: *)
(*Flux of charge constant*)


(* ::Input:: *)
(*Sum[(2L+1)NIntegrate[x( Abs[InterpolatedTransmission[L][x]]^2/((x- q Q)(Exp[(2 Pi (x -q Q))/\[Kappa]]-1))-Abs[InterpolatedTransmission[L][-x]]^2/((x+ q Q)(Exp[(2 Pi (x +q Q))/\[Kappa]]-1))),{x,0,\[Omega]f}],{L,0,10,1}]*)


(* ::Text:: *)
(*Flux of energy constant*)


(* ::Input:: *)
(*Sum[(2L+1)NIntegrate[x^2 ( Abs[InterpolatedTransmission[L][x]]^2/((x- q Q)(Exp[(2 Pi (x -q Q))/\[Kappa]]-1))-Abs[InterpolatedTransmission[L][-x]]^2/((x+ q Q)(Exp[(2 Pi (x +q Q))/\[Kappa]]-1))),{x,0,\[Omega]f}],{L,0,10,1}]*)


(* ::Chapter::Closed:: *)
(*Integrated quantities*)


(* ::Text:: *)
(*Note: It is better to use this part of the notebook separately, together with the Initializing section.*)


(* ::Input:: *)
(*DirPrefix = NotebookDirectory[]<>"integration"*)


(* ::Input:: *)
(*SetDirectory[DirPrefix]*)


(* ::Text:: *)
(*All the data that we have, for now, are for Q = 0.5, hence the background parameters are:*)


(* ::Input:: *)
(*rp=1;*)
(*rm=1/4;*)


(* ::Text:: *)
(*We also hardcode the scalar field charge list:*)


(* ::Input:: *)
(*qlist=Range[0,2,1/20];*)


(* ::Text:: *)
(*Functions to load in the data:*)


(* ::Input:: *)
(*Condensate[name_,q_]:=Import["./"<>name<>"/data/"<>name<>"_condensate_Q_50_q_"<>ToString[IntegerPart[100*q]]<>".dat"];*)
(*Currentr[name_,q_]:=Import["./"<>name<>"/data/"<>name<>"_current_r_Q_50_q_"<>ToString[IntegerPart[100*q]]<>".dat"];*)
(*Currentt[name_,q_]:=Import["./"<>name<>"/data/"<>name<>"_current_t_Q_50_q_"<>ToString[IntegerPart[100*q]]<>".dat"];*)
(*Stresstr[name_,q_]:=Import["./"<>name<>"/data/"<>name<>"_Stress_tr_Q_50_q_"<>ToString[IntegerPart[100*q]]<>".dat"];*)
(*Stresstt[name_,q_]:=Import["./"<>name<>"/data/"<>name<>"_Stress_tt_Q_50_q_"<>ToString[IntegerPart[100*q]]<>".dat"];*)
(*Stressthetatheta[name_,q_]:=Import["./"<>name<>"/data/"<>name<>"_Stress_thetatheta_Q_50_q_"<>ToString[IntegerPart[100*q]]<>".dat"];*)
(*Stressrr[name_,q_]:=Import["./"<>name<>"/data/"<>name<>"_Stress_rr_Q_50_q_"<>ToString[IntegerPart[100*q]]<>".dat"];*)


(* ::Text:: *)
(*For better readability, we extract the radial points into a list.*)


(* ::Input:: *)
(*rlist = Condensate["B_m_pB",1/2][[;;,1]];*)


(* ::Text:: *)
(*We define a few derived quantities with that.*)


(* ::Input:: *)
(*rsquaredCurrentr[name_,q_]:=Table[{rlist[[i]],rlist[[i]]^2 Currentr[name,q][[i,2]]},{i,1,Length[rlist]}];*)
(*rsquaredStressrsupert[name_,q_]:=Table[{rlist[[i]],f[rlist[[i]]]rlist[[i]]^2 Stresstr[name,q][[i,2]]},{i,1,Length[rlist]}];*)


(* ::Text:: *)
(*U,V components of the current. Those are actually J^{U/V}/\kappa U/V.*)


(* ::Input:: *)
(*CurrentU[name_,q_]:=Table[{rlist[[i]],-Currentt[name,q][[i,2]]+1/f[rlist[[i]]] Currentr[name,q][[i,2]]},{i,1,Length[rlist]}];*)
(*CurrentV[name_,q_]:=Table[{rlist[[i]],Currentt[name,q][[i,2]]+1/f[rlist[[i]]] Currentr[name,q][[i,2]]},{i,1,Length[rlist]}];*)


(* ::Text:: *)
(*These "modified" quantities have a power of f(r) multiplying the original quantity, in order to cancel the vanishing behavior of U or V.*)


(* ::Input:: *)
(*modifiedCurrentU[name_,q_]:=Module[{temp=CurrentU[name,q]},*)
(*Table[{rlist[[i]],f[rlist[[i]]]^1 temp[[i,2]]},{i,1,Length[rlist]}]]*)
(*modifiedCurrentV[name_,q_]:=Module[{temp=CurrentV[name,q]},*)
(*Table[{rlist[[i]],f[rlist[[i]]]^1 temp[[i,2]]},{i,1,Length[rlist]}]]*)


(* ::Text:: *)
(*U,V components of the stress-energy tensor. These need to be multiplied by the corresponding functions of U and V, which can be seen on the draft. So, for instance, StressUU is actually $\kappa U^2 \langle T_{UU} \rangle$.*)


(* ::Input:: *)
(*StressUU[name_,q_]:=Table[{rlist[[i]],1/4 (Stresstt[name,q][[i,2]]-2 f[rlist[[i]]] Stresstr[name,q][[i,2]]+f[rlist[[i]]]^2 Stressrr[name,q][[i,2]])},{i,1,Length[rlist]}];*)
(*StressUV[name_,q_]:=Table[{rlist[[i]],-(1/4)(Stresstt[name,q][[i,2]]-f[rlist[[i]]]^2 Stressrr[name,q][[i,2]]) },{i,1,Length[rlist]}];*)
(*StressVV[name_,q_]:=Table[{rlist[[i]],1/4 (Stresstt[name,q][[i,2]]+2 f[rlist[[i]]]Stresstr[name,q][[i,2]]+f[rlist[[i]]]^2 Stressrr[name,q][[i,2]] )},{i,1,Length[rlist]}];*)


(* ::Text:: *)
(*Again, "modified" means that we multiply the original quantity by the power of f(r) needed to cancel the vanishing behavior of U or V.*)


(* ::Input:: *)
(*modifiedStressUU[name_,q_]:=Module[{temp=StressUU[name,q]},*)
(*Table[{rlist[[i]],f[rlist[[i]]]^-1 temp[[i,2]]},{i,1,Length[rlist]}]]*)
(*modifiedStressUV[name_,q_]:=Module[{temp=StressUV[name,q]},*)
(*Table[{rlist[[i]],f[rlist[[i]]]^-1 temp[[i,2]]},{i,1,Length[rlist]}]]*)
(*modifiedStressVV[name_,q_]:=Module[{temp=StressVV[name,q]},*)
(*Table[{rlist[[i]],f[rlist[[i]]]^-2 temp[[i,2]]},{i,1,Length[rlist]}]]*)


(* ::Text:: *)
(*We also define ingoing null coordinate system (or ingoing E-F coordinates) components for the SET: *)


(* ::Input:: *)
(*EFStressvv[name_,q_]:=Table[{rlist[[i]],Stresstt[name,q][[i,2]]+2 f[rlist[[i]]] Stresstr[name,q][[i,2]]+ f[rlist[[i]]]^2 Stressrr[name,q][[i,2]]},{i,1,Length[rlist]}];*)
(*EFStressvr[name_,q_]:=Table[{rlist[[i]],(-(1/f[rlist[[i]]])Stresstt[name,q][[i,2]]+ f[rlist[[i]]] Stressrr[name,q][[i,2]])},{i,1,Length[rlist]}];*)
(*EFStressrr[name_,q_]:=Table[{rlist[[i]],Stressrr[name,q][[i,2]]},{i,1,Length[rlist]}];*)


(* ::Text:: *)
(*Use this as a placeholder to compute quantities you want to locally modify (in the notebook scope) before plotting:*)


(* ::Input:: *)
(*LocallymodifiedQuantity[name_,q_]:=Module[{temp=Stressthetatheta[name,q]},*)
(*Table[{rlist[[i]],f[rlist[[i]]]^2 rlist[[i]]^-2 temp[[i,2]]},{i,1,Length[rlist]}]]*)
(*ListPlot[Table[LocallymodifiedQuantity[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Subchapter:: *)
(*Plots*)


(* ::Text:: *)
(*We control the sample of scalar field charges with this variable. Note that these are given in rp units. For a black hole with Q=0.8 M, the numbers below should be multiplied by 1.6 to get the correct q/M number.*)


(* ::Input:: *)
(*qsample={0,1/4,1/2,1};*)


(* ::Text:: *)
(*The sample we are currently using in the paper (with the addition of q=0):*)


(* ::Input:: *)
(*qsample = {0,1/10,2/10,3/10,4/10,5/10}*)


(* ::Text:: *)
(*Generating some legends:*)


(* ::Input:: *)
(*qsamplelegends=LineLegend[Automatic,StringJoin["q=",#]&/@Map[ToString,N[qsample]]];*)


(* ::Section::Closed:: *)
(*Past Unruh minus past Boulware*)


(* ::Input:: *)
(*StateNames="pUnruh_m_pB"*)


(* ::Subsection:: *)
(*Scalar condensate*)


(* ::Text:: *)
(*We note here that the variation in charge has almost no bearing on the changes in the difference in expectation values.*)


(* ::Input:: *)
(*ListPlot[Table[Condensate[StateNames,q],{q,qsample}],Joined->True(*,PlotRange\[Rule]All*)]*)


(* ::Text:: *)
(*A small zoom:*)


(* ::Input:: *)
(*ListPlot[Table[Condensate[StateNames,q],{q,qsample}],Joined->True,PlotRange->{{1.05,1.2},{0,0.005}},PlotLegends->qsamplelegends]*)


(* ::Text:: *)
(*The divergence seems to improve when we multiply by f(r), suggesting that the condensate behaves like f(r)^-1 near the horizon. The sharp peak looks like a numerical error but see below.*)


(* ::Input:: *)
(*LocallymodifiedCondensate[name_,q_]:=Module[{temp=Condensate[name,q]},*)
(*Table[{rlist[[i]],f[rlist[[i]]]^1 temp[[i,2]]},{i,1,Length[rlist]}]]*)
(*ListPlot[Table[LocallymodifiedCondensate[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Text:: *)
(*A little zoom shows the peak to actually be quite smooth.*)


(* ::Input:: *)
(*ListPlot[Table[LocallymodifiedCondensate[StateNames,q],{q,qsample}],Joined->True,PlotRange->{{1,1.5},Automatic},PlotLegends->qsamplelegends]*)


(* ::Subsection:: *)
(*Current*)


(* ::Text:: *)
(*This component helps us check our numerical results (it has to be a straight line). From the conservation equations, expectation values of r^2 J^r have the form of -\[CapitalKappa], which is a measure of the charge flux. In the present case, this shows us that the past Unruh state has a larger charge flux than the past Boulware state.*)


(* ::Input:: *)
(*ListPlot[Table[rsquaredCurrentr[StateNames,q],{q,qsample}],Joined->True(*,PlotRange\[Rule]All*),PlotLegends->qsamplelegends]*)


(* ::Text:: *)
(*The difference in expectation values of the charge density seems to not be well-behaved at the horizon.*)


(* ::Input:: *)
(*ListPlot[Table[Currentt[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Text:: *)
(*When multiplied by one power of f(r), it is better behaved -- the abrupt change near the horizon looks like numerical error. The fact it is negative implies that the charge polarization in the past Unruh state is smaller than the one in the past Boulware state.*)


(* ::Input:: *)
(*LocallymodifiedCurrentt[name_,q_]:=Module[{temp=Currentt[name,q]},*)
(*Table[{rlist[[i]],f[rlist[[i]]]^1 temp[[i,2]]},{i,1,Length[rlist]}]]*)
(*ListPlot[Table[LocallymodifiedCurrentt[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Text:: *)
(*Two powers of f(r) make a prettier plot.*)


(* ::Input:: *)
(*LocallymodifiedCurrentt[name_,q_]:=Module[{temp=Currentt[name,q]},*)
(*Table[{rlist[[i]],f[rlist[[i]]]^2 temp[[i,2]]},{i,1,Length[rlist]}]]*)
(*ListPlot[Table[LocallymodifiedCurrentt[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Text:: *)
(*Both U,V components look divergent at the horizon.*)


(* ::Input:: *)
(*ListPlot[Table[CurrentU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Input:: *)
(*ListPlot[Table[CurrentV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Text:: *)
(*The modified components look better. On the future horizon (U=0), we expect the past Unruh state to be well-behaved, but not the past Boulware state.*)


(* ::Input:: *)
(*ListPlot[Table[modifiedCurrentU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Text:: *)
(*A little zoom. This seems to indicate a switch happens at a certain value of the scalar field charge. Needs further investigation.*)


(* ::Input:: *)
(*ListPlot[Table[modifiedCurrentU[StateNames,q],{q,qsample}],Joined->True,PlotRange -> {{1,1.5},{-10^-6,3*10^-6}}]*)


(* ::Text:: *)
(*In the past horizon, V=0, and we need to multiply our result associated to J^Vby f(r) (recall that CurrentV is actually J^V/V), in order to cancel the vanishing behavior of V. However, we expect both states to be divergent at the past horizon, hence a (possible) divergent behavior here is not unsurprising.*)


(* ::Input:: *)
(*ListPlot[Table[modifiedCurrentV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Subsection:: *)
(*Stress-energy tensor*)


(* ::Text:: *)
(*All associated stress-energy tensor components seem to diverge here, so we have to dig a little further.*)


(* ::Input:: *)
(*ListPlot[Table[StressUU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[StressUV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[StressVV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Text:: *)
(*Even cancelling the vanishing behaviors of U and V at their respective horizons do not improve the situation. This is probably due to the fact that the past Boulware state is divergent at both horizons and there is not much hope for well-behaved results here -- only a fortuitous cancellation of divergences pertaining to both states (such as in, say, the past horizon) would lead to well-behaved quantities.*)


(* ::Input:: *)
(*ListPlot[Table[modifiedStressUU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[modifiedStressUV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[modifiedStressVV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Section::Closed:: *)
(*Past CCH minus past Unruh*)


(* ::Text:: *)
(*The past CCH is one of the states we want to investigate thoroughly, since its properties are not as well-known as, say, the past Boulware one. From the asymptotics we expect its regularity properties to be similar to the properties of the past Unruh state. *)


(* ::Input:: *)
(*StateNames = "pCCH_m_pUnruh"*)


(* ::Subsection:: *)
(*Scalar condensate*)


(* ::Text:: *)
(*There is a large variation when we vary the scalar field charge. The plot suggests that the this difference in expectation values is well-behaved at the horizon. They also tend to constant values, implying the the past CCH is not empty at infinity. We have already seen that this does not happen with the past Unruh state.*)


(* ::Input:: *)
(*ListPlot[Table[Condensate[StateNames,q],{q,qsample}],Joined->True(*,PlotRange\[Rule]All*)]*)


(* ::Subsection:: *)
(*Current*)


(* ::Text:: *)
(*As expected, the quantity associated to the charge flux increases as we increase the charge, which implies the the charge flux of the past CCH state is larger than the one in the past Unruh state.*)


(* ::Input:: *)
(*ListPlot[Table[rsquaredCurrentr[StateNames,q],{q,qsample}],Joined->True(*,PlotRange\[Rule]All*)]*)


(* ::Text:: *)
(*Let us check the degree of divergence in the charge density.*)


(* ::Input:: *)
(*ListPlot[Table[Currentt[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Text:: *)
(*It seems that the charge density difference diverges as (f(r)^-1).*)


(* ::Input:: *)
(*LocallymodifiedCurrentt[name_,q_]:=Module[{temp=Currentt[name,q]},*)
(*Table[{rlist[[i]],f[rlist[[i]]]^1 temp[[i,2]]},{i,1,Length[rlist]}]]*)
(*ListPlot[Table[LocallymodifiedCurrentt[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Text:: *)
(*As par for the course, the CurrentU quantity does not look good at the horizon.*)


(* ::Input:: *)
(*ListPlot[Table[CurrentU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Text:: *)
(*Cancelling the vanishing behavior of U, we seem to get a well-behaved quantity. This implies that this quantity is finite at the future horizon, as expected.*)


(* ::Input:: *)
(*ListPlot[Table[modifiedCurrentU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Text:: *)
(*CurrentV quantity looks fine, apart from the discontinuities near the horizon, which look numerical errors.*)


(* ::Input:: *)
(*ListPlot[Table[CurrentV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Text:: *)
(*A zoom of the discontinuity region might help us:*)


(* ::Input:: *)
(*ListPlot[Table[CurrentV[StateNames,q],{q,qsample}],Joined->True,PlotRange->{{1,1.3},{-10^-6,6*10^-5}}]*)


(* ::Text:: *)
(*Interestingly, the J^V component differences look well-behaved at the past horizon. Write about J^V vanishing at the horizon, which agrees with asymptotics.*)


(* ::Input:: *)
(*ListPlot[Table[modifiedCurrentV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Subsection:: *)
(*Stress-energy tensor*)


(* ::Text:: *)
(*The associated quantities to SET differences all look fine. However, we must cancel the vanishing behavior of U and V at the respective horizons.*)


(* ::Input:: *)
(*ListPlot[Table[StressUU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[StressUV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[StressVV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Text:: *)
(*Cancelling the vanishing behaviors suggests that all the SET components are pathological in, at least, one horizon. The exception could be Subscript[T, UV], but I am not sure about that. *)


(* ::Input:: *)
(*ListPlot[Table[modifiedStressUU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[modifiedStressUV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[modifiedStressVV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Input:: *)
(*LocallymodifiedStressUU[name_,q_]:=Module[{temp=StressUU[name,q]},*)
(*Table[{rlist[[i]],-f[rlist[[i]]]^-1temp[[i,2]]},{i,1,Length[rlist]}]]*)
(*ListPlot[Table[LocallymodifiedStressUU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Section::Closed:: *)
(*Boulware minus past Boulware*)


(* ::Input:: *)
(*StateNames="B_m_pB"*)


(* ::Subsection:: *)
(*Scalar condensate*)


(* ::Input:: *)
(*ListPlot[Table[Condensate[StateNames,q],{q,qsample}],Joined->True(*,PlotRange\[Rule]All*),PlotLegends->qsamplelegends]*)


(* ::Subsection:: *)
(*Current*)


(* ::Input:: *)
(*ListPlot[Table[rsquaredCurrentr[StateNames,q],{q,qsample}],Joined->True(*,PlotRange\[Rule]All*)]*)


(* ::Input:: *)
(*ListPlot[Table[Currentt[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Input:: *)
(*ListPlot[Table[CurrentU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Input:: *)
(*ListPlot[Table[CurrentV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Input:: *)
(*ListPlot[Table[modifiedCurrentU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Input:: *)
(*ListPlot[Table[modifiedCurrentV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Subsection:: *)
(*Stress-energy tensor*)


(* ::Input:: *)
(*ListPlot[Table[StressUU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[StressUV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[StressVV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Input:: *)
(*ListPlot[Table[modifiedStressUU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[modifiedStressUV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[modifiedStressVV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Section::Closed:: *)
(*FT minus past Unruh*)


(* ::Input:: *)
(*StateNames="FT_m_pUnruh"*)


(* ::Subsection:: *)
(*Scalar condensate*)


(* ::Input:: *)
(*ListPlot[Table[Condensate[StateNames,q],{q,qsample}],Joined->True(*,PlotRange\[Rule]All*),PlotLegends->qsamplelegends]*)


(* ::Subsection:: *)
(*Current*)


(* ::Input:: *)
(*ListPlot[Table[rsquaredCurrentr[StateNames,q],{q,qsample}],Joined->True(*,PlotRange\[Rule]All*),PlotLegends->qsamplelegends]*)


(* ::Input:: *)
(*ListPlot[Table[Currentt[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Input:: *)
(*ListPlot[Table[CurrentU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Input:: *)
(*ListPlot[Table[CurrentV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All,PlotLegends->qsamplelegends]*)


(* ::Input:: *)
(*ListPlot[Table[modifiedCurrentU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All,PlotLegends->qsamplelegends]*)


(* ::Input:: *)
(*ListPlot[Table[modifiedCurrentV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All,PlotLegends->qsamplelegends]*)


(* ::Subsection:: *)
(*Stress-energy tensor*)


(* ::Input:: *)
(*ListPlot[Table[StressUU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All,PlotLegends->qsamplelegends]*)
(*ListPlot[Table[StressUV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All,PlotLegends->qsamplelegends]*)
(*ListPlot[Table[StressVV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All,PlotLegends->qsamplelegends]*)


(* ::Input:: *)
(*ListPlot[Table[modifiedStressUU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All,PlotLegends->qsamplelegends]*)
(*ListPlot[Table[modifiedStressUV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All,PlotLegends->qsamplelegends]*)
(*ListPlot[Table[modifiedStressVV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All,PlotLegends->qsamplelegends]*)


(* ::Input:: *)
(*ListPlot[Table[EFStressUU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All,PlotLegends->qsamplelegends]*)
(*ListPlot[Table[EFStressUr[StateNames,q],{q,qsample}],Joined->True,PlotRange->All,PlotLegends->qsamplelegends]*)
(*ListPlot[Table[EFStressrr[StateNames,q],{q,qsample}],Joined->True,PlotRange->All,PlotLegends->qsamplelegends]*)


(* ::Section::Closed:: *)
(*H minus FT*)


(* ::Input:: *)
(*StateNames="H_m_FT"*)


(* ::Subsection:: *)
(*Scalar condensate*)


(* ::Input:: *)
(*ListPlot[Table[Condensate[StateNames,q],{q,qsample}],Joined->True(*,PlotRange\[Rule]All*)]*)


(* ::Subsection:: *)
(*Current*)


(* ::Input:: *)
(*ListPlot[Table[rsquaredCurrentr[StateNames,q],{q,qsample}],Joined->True(*,PlotRange\[Rule]All*)]*)


(* ::Input:: *)
(*ListPlot[Table[Currentt[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Input:: *)
(*ListPlot[Table[CurrentU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Input:: *)
(*ListPlot[Table[CurrentV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Input:: *)
(*ListPlot[Table[modifiedCurrentU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Input:: *)
(*ListPlot[Table[modifiedCurrentV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Subsection:: *)
(*Stress-energy tensor*)


(* ::Input:: *)
(*ListPlot[Table[StressUU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[StressUV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[StressVV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Input:: *)
(*ListPlot[Table[modifiedStressUU[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[modifiedStressUV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)
(*ListPlot[Table[modifiedStressVV[StateNames,q],{q,qsample}],Joined->True,PlotRange->All]*)


(* ::Section:: *)
(*H minus past Unruh*)


(* ::Text:: *)
(*This one is a little different since we can combine results from "H minus FT" and "FT minus past Unruh" together, so no need to import data.*)


(* ::Input:: *)
(*HmpUQuantity[QuantityName_,q_]:=Module[{temp1=QuantityName["H_m_FT",q],temp2=QuantityName["FT_m_pUnruh",q]},*)
(*Table[{rlist[[i]],temp1[[i,2]]+temp2[[i,2]]},{i,1,Length[rlist]}]]*)


(* ::Input:: *)
(*ListPlot[Table[HmpUQuantity[modifiedStressUU,q],{q,qsample}],Joined->True(*,PlotRange\[Rule]All*),PlotLegends->qsamplelegends]*)


(* ::Input:: *)
(*Namegenerator[QuantityName_]:=( *)
(*quantitystring=ToString[QuantityName];*)
(*Which[quantitystring == "Condensate",string="condensate",quantitystring == "Currentt",string="current_t",quantitystring == "Currentr",string="current_r",quantitystring == "Stresstt",string="Stress_tt",quantitystring == "Stresstr",string="Stress_tr",quantitystring == "Stressrsupertcomponent",string="Stress_rsupert",quantitystring == "Stresstphi",string="Stress_tphi",quantitystring == "Stressrr",string="Stress_rr",quantitystring == "Stressrtheta",string="Stress_rtheta",quantitystring == "Stressrphi",string="Stress_rphi",quantitystring == "Stressthetatheta",string="Stress_thetatheta",quantitystring == "Stressthetaphi",string="Stress_thetaphi",quantitystring == "Stressphiphi",string="Stress_phiphi"];*)
(*)*)


(* ::Input:: *)
(*Do[*)
(*name=Stressthetatheta;*)
(*Namegenerator[name];*)
(*Export["H_m_pUnruh_"<>string<>"_Q_50_q_"<>ToString[IntegerPart[100*q]]<>".dat",HmpUQuantity[name,q]];*)
(*,{q,qsample}]*)
