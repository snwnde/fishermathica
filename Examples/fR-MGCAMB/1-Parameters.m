(* ::Package:: *)

(* ::Title:: *)
(*Configuration Parameters*)


(* ::Section:: *)
(*Numerical Parameters*)


$interpOrd=1; (* interpolation order*)


$interpMeth="Spline";  (* interpolation method for input files*)


SetOptions[Interpolation, InterpolationOrder->$interpOrd, Method->$interpMeth]


Options[NIntegrate]


SetOptions[NIntegrate, AccuracyGoal->7, PrecisionGoal->Automatic, MinRecursion->2, MaxRecursion->1000, WorkingPrecision->Floor[$MachinePrecision], Method->{Automatic,"SymbolicProcessing"->0}]


(* ::Section::Bold:: *)
(*Cosmological Parameters, choose fiducial values.*)


(* ::Code::Bold:: *)
(**)


(* ::Code::Bold:: *)
(**)


$newrecipelinear=True


If[$newrecipelinear==True,
sigmapnlref=0;
sigmavnlref=0;
,
sigmapnlref=10.728;
sigmavnlref=4.8703;
];


(*mysigmav = Sqrt[pthetathetaInt[$zmeansurvey]]
mysigmap = Sqrt[pthetathetaInt[$zmeansurvey]]/sigma8ofZ[$zmeansurvey]*)


(*(mysigmav-sigmavnlref)/sigmavnlref*100*)


(*(mysigmap-sigmapnlref)/sigmapnlref*100*)


Get[inputsfolder<>"GenerationsParameters-"<>$fiducialstring<>".wl"]


NotebookDirectory[]


NotebookSave[]


(*Intrinsic Alignment*)


$varyIAparams=False;


$includeIAterms=True;


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
