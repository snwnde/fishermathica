(* ::Package:: *)

(* ::Section::Bold:: *)
(*Cosmological Parameters*)


setPScosmoOptions[kdependentGrowth->True, lcdmBool->False, linearBool->True];


(* ::Subsection:: *)
(*Set global cosmological options*)


modelflatness=True;
mdflstr=If[modelflatness==True, "Flat", "NonFlat"];
If[modelflatness!=True,
Print["MGCAMB models do not accept non-flatness"];
Abort[];]


(* ::Subsection:: *)
(*Standard Cosmological Parameters*)


$cosmoParametersProtectedList


nsref=0.96; 
as109ref=2.12605; 
hubbleref=0.67; 
weq0=-1.0; 
weqa=0.;
OmegaMref=0.320; 
Omegabref=0.050;


(* ::Subsection:: *)
(*Transformations into derived parameters*)


omegamref=OmegaMref*hubbleref^2; 
omegabref=Omegabref*hubbleref^2;
OmegaKref=0;
OmegaDEref=Which[
modelflatness==True,
1-OmegaMref;
,
modelflatness==False,
1-OmegaMref-OmegaKref;
];
Asref=as109ref*10^-9;
ln1010as=Log[10^10 ass];
funcAs[as_]=(10^9 as);


(* ::Subsection:: *)
(*Beyond Standard, Modified Gravity Parameters*)


fR0ref=5*10^-5;


logfR0ref=N@Log10[fR0ref]


(*Print["Modified Gravity parameters: ", "E11=", E11ref, ", E22=", E22ref];*)


(* ::Subsection:: *)
(*Neutrino physics parameters*)


Mnutotalref=(*0*)0.06(*0.15*); (*eV*)
Neff = 3.046;


Which[Mnutotalref==0.06,
$massivenus=1; (*1, for Mnu=0.06 *)(*3, for Mnu=0.15 *)(*integer of number of massive neutrinos *)
sigma8ref=0.911031;
,
Mnutotalref==0.,
$massivenus=0; 
sigma8ref=0.831933;
];
Print["Reference sigma8, to be compared to (MG)CAMB results:= ", sigma8ref]


(* ::Subsubsection:: *)
(*Derived Neutrino and Early Universe parameters*)


omgref = (2.469/hubbleref^2)*10^(-5);
omegarref = omgref*(1. + 0.2271*Neff);
omegarh2ref = omegarref*hubbleref^2;
OmegaRDS=omegarref; 

$masslessnus=Neff-$massivenus;
numassfac=94.07; (*eV*)
omeganuh2=Sum[(((Neff/3)^(3/4))*(Mnutotalref/$massivenus))/numassfac, {i,1,$massivenus}];
Omeganuref=omeganuh2/(hubbleref^2);

Omegacref=OmegaMref-Omegabref-Omeganuref;
omegacref=Omegacref*hubbleref^2;


(* ::Subsection:: *)
(*Power Spectrum Parameters*)


$halofitversion=9;


Print["Halofit Version for CAMB: ", $halofitversion];


(* ::Section::Bold:: *)
(*Set CosmologyFunctions rules*)


fiducialModel={OmegaMref,Omegabref,hubbleref,nsref,sigma8ref,logfR0ref};
parnamesModel={Omegam,Omegab , hubble,ns,sigma8,logfR0};
parlabelsModel={"Om","Ob", "h","ns","s8","logfR0"};
parstorun=Range[Length[parnamesModel]];


setParameterOptions[parnamesModel,fiducialModel, parlabelsModel]


Print["Cosmological parameters to be varied: "]
Print[$paramoptions]


setCosmologyFixedValues[internalPkUnits->"h/Mpc", OmegaExtraMatterSpecies->Omeganuref, sigma8reference->sigma8ref, 
Asreference->Asref, internalHubbleUnits->"1/Mpc",
internalDistanceUnits->"Mpc"]


Print["Fixed cosmological values: "]
Print[Options@setCosmologyFixedValues]


(*setORreadTransformedParameters[As, generateInverseFunction->True, explicitTransformFunction->funcAs]
scalarAmplitude[]*)


SetOptions[LCDMCAMBPsPre, activateModifiedGravity->"fofR-HS"];


Options[parametersModifiedGravity]


parametersModifiedGravity[logfR0]


Options[LCDMCAMBPsPre]
