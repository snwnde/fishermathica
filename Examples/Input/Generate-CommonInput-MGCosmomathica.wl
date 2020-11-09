(* ::Package:: *)

(* ::Title::Bold:: *)
(*Weak Lensing and Photometric Probes, IST, input files*)


(* ::Subtitle::Bold:: *)
(*Santiago Casas, 2020*)


(* ::Text::Bold:: *)
(**)


(* ::Section:: *)
(*Initialization*)


DirectoryName@$InputFileName
notebookdir=If[aa=DirectoryName@$InputFileName; aa=="", NotebookDirectory[], aa]
SetDirectory[notebookdir];
Get["Init.m"]


initialdate=DateString[];


Print["Initial DateTime: ", initialdate]


(* ::Section:: *)
(*Run type specifications*)


$debugPlots=True;
$debugPrint=False;


$epslist = {0.00625,  0.0125, 0.01875, 0.025, 0.0375, 0.05, 0.10} ;  (*List of epsilon variations for which files will be produced*)


(*$epslist = {0.00625, 0.05, 0.10} ;  (*List of epsilon variations for which files will be produced*)*)


$observableFileType="commonInput"


$fidumodel="-GRlimit"


(* ::Section:: *)
(*Cosmological Parameters, choose fiducial values*)


Which[
$fidumodel=="-GRlimit",
Print[$fidumodel];
Get["GenerationsParameters-GRlimit.wl"]
,
$fidumodel=="-HS5",
Print[$fidumodel];
Get["GenerationsParameters-HS5.wl"]
]


(* ::Section::Bold:: *)
(*Test parameters of cosmology*)


Print["Cosmological species, energy fractions: "]
Print@Thread[Rule[(ToString[#]&/@{OmegaK0Today,OmegaM0Today,OmegaCDM0Today,OmegaBaryon0Today,OmegaDE0Today,1-OmegaDE0Today-OmegaM0Today,OmegaR, Omeganu, hubbleToday}),
{OmegaK0Today[],OmegaM0Today[],OmegaCDM0Today[],OmegaBaryon0Today[],OmegaDE0Today[],1-OmegaDE0Today[]-OmegaM0Today[],$OmegaR,$Omeganu,hubbleToday[]}]]


Print["Options set for modified gravity: "];
Print[Options[parametersModifiedGravity]]


(* ::Section:: *)
(*Choose redshift range*)


minz=0.; maxz=2.5; stepz=100(*302*);


$zrange=Range[minz,maxz,(maxz-minz)/stepz];


Print["Number of redshifts to produce: " , Length@$zrange]


Divisors[Length@$zrange]


(*listof$zrange=Partition[$zrange,101];  (*Partition into three lists, since CAMB only can handle 101 redshifts at a time*)*)


listof$zrange=$zrange


(* ::Section:: *)
(*Choose wavenumber k range*)


$transferKmax=50;
$transferKperlogint = 50;


(* ::Section::Bold:: *)
(*Set accuracy options and test for correct output from CAMB*)


$accboost=2;
$transferHighPrec=True;


SetOptions[LCDMCAMBPsPre, {AccuracyBoost->$accboost, TransferHighPrecision->$transferHighPrec, TransferKperLogInt->$transferKperlogint}]


SetOptions[LCDMCAMBPsPre, {TransferKmax->$transferKmax, HalofitVersion->$halofitversion}]


SetOptions[CAMB, debugIntFloats->False];


SetOptions[LCDMCAMBPsPre,{OmegaNu->$Omeganu, MasslessNeutrinos->$masslessnus, MassiveNeutrinos->$massivenus, NuMassFractions->{1}}]


(*myInterrupt[]*)


(* ::Section::Bold:: *)
(*Compute example run to get sigma8*)


sigma8valofz = LCDMCAMBPsPre[0., 
returnInterpolated->False, returnSigma8->True, checkConsistency->False]


Options[LCDMCAMBPsPre]


psfulltab=LCDMCAMBPsPre[0., 
returnInterpolated->False, returnSigma8->True, checkConsistency->True];


If[$debugPlots,
ListLogLogPlot[{psfulltab[[All,{1,2}]],psfulltab[[All,{1,3}]]}]
]


Print["sigma8 from CAMB, vs sigma8reference defined in parameters: ", " CAMB: ", sigma8valofz, " Reference: ", $sigma8reference]


Options[LCDMCAMBPsPre]


{outZtab,outKtab,HofZtab,
    s8ofZtab,growthTab,fgrowthrateTab,
    psLinTab, psNonLinTab
    }=LCDMCAMBPsPre[listof$zrange, checkConsistency->False];


outZtab//Dimensions


outKtab//Dimensions


HofZtab//MinMax


s8ofZtab//Dimensions


growthTab//Dimensions


outKtab[[1;;-1;;210]]


ListLogLogPlot[Table[Transpose[{outZtab,fgrowthrateTab[[All,ii]]}],{ii,1,Length@outKtab,210}], Joined->True, PlotLegends->outKtab[[1;;-1;;210]]]


fgrowthrateTab//Dimensions


psLinTab//Dimensions


ListLogLogPlot[Table[Transpose[{outKtab,psLinTab[[ii,All]]}],{ii,1,101,20}]]


ListLogLogPlot[Table[Transpose[{outKtab,psNonLinTab[[ii,All]]}],{ii,1,101,20}]]


(* ::Subsection:: *)
(*Produce Files*)


$observableFileType


$filesdir=notebookdir<>"/WP6-fofR-"<>$observableFileType<>"-Cosmomathica-CommonInputFiles-baseline_cosmology-DirectPk-"<>"-"<>replaceStringPoint[Mnutotalref,"Mnu"]<>$fidumodel<>"/";
mkDirectory[$filesdir]


Options[computeEpsilonParameterFiles]


myInterrupt[]


$epslist


$paramoptions


(*{optplus,optminus,step}=numericalDerivativeStep[5,0.1]*)


(*UnsameQ[$paramoptions[[5]][[1]] , sigma8]*)


(*Asoptplus=LCDMCAMBPsPre[0.,  {optplus}~Join~{checkConsistency\[Rule]False, returnRescaledAs->True}]*)


(*Asoptplus/$Asreference*)


parametersModifiedGravity[logfR0]


timetotal={};


Do[
timesec=computeEpsilonParameterFiles[listof$zrange, 
epsii, computeEpsilons->True, computeFiducial->True, 
referenceAmplitude->$Asreference, TransferKmax->$transferKmax, HalofitVersion->$halofitversion, referenceSigma8->$sigma8reference,
exportDirectory->$filesdir, extraNumericCAMBparams->{AccuracyBoost->$accboost, TransferHighPrecision->$transferHighPrec, TransferKperLogInt->$transferKperlogint},
exportDetailedInfo->True,
constantNeutrinoParams->{OmegaNu->$Omeganu, MasslessNeutrinos->$masslessnus, MassiveNeutrinos->$massivenus, NuMassFractions->{1}}];
AppendTo[timetotal, timesec],
{epsii, $epslist}
]


Total[timetotal]/60/60 //N


$epslist


finalDate=DateString[];


Print["Final DateTime: ", finalDate]


timeElapsed[initialdate, finalDate];


Quit[]
