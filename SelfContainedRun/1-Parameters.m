(* ::Package:: *)

(* ::Title:: *)
(*Configuration Parameters*)


(* ::Section:: *)
(*Numerical Parameters*)


interpOrd=3; (* interpolation order*)


interpMeth="Spline";  (* interpolation method for input files*)


(* ::Section::Bold:: *)
(*Cosmological Parameters, choose fiducial values.*)


(* ::Code::Bold:: *)
(**)


(* ::Code::Bold:: *)
(**)


$newrecipelinear=False


If[$newrecipelinear==True,
sigmapnlref=0;
sigmavnlref=0;
,
sigmapnlref=10.728;
sigmavnlref=4.8703;
];



NotebookDirectory[]


NotebookSave[]


(*Intrinsic Alignment*)


$varyIAparams=False;


$includeIAterms=False;


If[$varyIAparams,
$includeIAterms=True; (*just to make sure it's always set when IA params are varied*)
aiafid=$aIAfidu;
etaiafid=$eIAfidu;
betaiafid=$bIAfidu;
$paramfidus = $paramfidus~Join~{aiafid, etaiafid, betaiafid};
$paramlabels = $paramlabels~Join~{"AIA","etaIA","betaIA"};
$paramnames = $paramnames~Join~{$aIA, $eIA, $bIA};
setParameterOptions[$paramnames,$paramfidus, $paramlabels];
];


$paramoptions


nsref=0.966; 
as109ref=2.12605;
hubbleref=0.72; 
weq0=-1.0; weqa=0.;
OmegaMref=0.320; Omegabref=0.050;
omegamref=OmegaMref*hubbleref^2; omegabref=Omegabref*hubbleref^2;
OmegaDEref=1-OmegaMref(*0.68*);
ass=as109ref*10^-9;
$Mnuchoice=0.06;

If[$newrecipelinear==True,
sigmapnlref=0;
sigmavnlref=0;
,
sigmapnlref=10.728;
sigmavnlref=4.8703;
];


Mnutotalref=(*0*)$Mnuchoice(*0.15*); (*eV*)
Which[Mnutotalref==0.06,
massivenus=1; (*1, for Mnu=0.06 *)(*3, for Mnu=0.15 *)(*integer of number of massive neutrinos *)
sigma8ref=0.815583;(* 0.815583 for 1 massive neutrino, 2 massless, M_nu=0.06*)
,
Mnutotalref==0.15,
massivenus=3; 
sigma8ref=0.794958; (* 0.794958 for 3 massive neutrinos, M_nu=0.15*) 
,
Mnutotalref==0.,
massivenus=0; 
sigma8ref=0.830(*0.830 for massless neutrinos*);
];


Neff = 3.046;
omgref = (2.469/hubbleref^2)*10^(-5);
omegarref = omgref*(1. + 0.2271*Neff);
omegarh2ref = omegarref*hubbleref^2;
OmegaRDS=omegarref;

masslessnus=Neff-massivenus;
numassfac=94.07; (*eV*)
omeganuh2=Sum[(((Neff/3)^(3/4))*(Mnutotalref/massivenus))/numassfac, {i,1,massivenus}];
Omeganuref=omeganuh2/(hubbleref^2)

Omegacref=OmegaMref-Omegabref-Omeganuref;
omegacref=Omegacref*hubbleref^2;

Switch[StringMatchQ["Exp*"][$nameOfModel],
True,
parnames={Omegam, Omegab, ns, As109,gamma};
fiducial={OmegaMref, Omegabref, nsref, as109ref,gammaref};
parlabels={"OmegaM", "OmegaB", "ns", "As10E9","gamma"};,
False,
parnames={Omegam, Omegab, ns, As109,hubble};
fiducial={OmegaMref, Omegabref, nsref, as109ref,hubbleref};
parlabels={"OmegaM", "OmegaB", "ns", "As10E9","hubble"};
];

setParameterOptions[parnames,fiducial,parlabels]


funcAs[as_]=(10^9 as);
setORreadTransformedParameters[As, generateInverseFunction->True, explicitTransformFunction->funcAs]


$paramoptions


setPowerSpectrumRelatedSpecs[dzRelativeVelocityError->0.001, APeffect->True, pks8Ratio->False]


Options[setNumericalEpsilon]


$epsilonzstep


setNumericalEpsilon[$epsilonstep, epsilonStepForZdependentParameters->$epsilonzstep];


SetOptions[dObsPowerDZpi, stencilpoints->$stencilpoints];


SetOptions[fGrowthRate,solveDifferentialEquation->True];


setCosmologyFixedValues[sigma8reference->sigma8ref, internalHubbleUnits->"1/Mpc",externalHubbleUnits->"1/Mpc",internalPkUnits->"h/Mpc", externalPkUnits->"h/Mpc", OmegaExtraMatterSpecies->Omeganuref, OmegaRadiation->OmegaRDS,hubbleReference->hubbleref]


setFisherKmaxValues[maxFisherKCut->kmaxLinear, fisherKCutMethod->1];


Print["parameters:  "<>ToString@$paramoptions, "---kmax="<>ToString@$kmaxHard]


$onlyzdeppars=False;


zdepVariablesVector={d2,d1};


FisherBlockBuilder[zdepVariablesVector, fisherBlockDerivativeMethod->"FullNumerical",fbsigma8Variables->True, onlyZdependentParameters->$onlyzdeppars]


Print["z-dependent Vector: "]


Print[$zdependDerivVector];
