(* ::Package:: *)

(* ::Subsection:: *)
(*Initialization and Specifications*)


$resultsDirectoryWL=mkDirectory[notebookdir<>"/Results-WL-"<>$dateshort<>"-"<>$stepstring<>"-"<>$surveychosen<>"/"]


$resultsDirectoryWL


specsDir


Options[NIntegrate]


SetOptions[NIntegrate, AccuracyGoal->11, PrecisionGoal->Automatic, MinRecursion->2, MaxRecursion->1000, WorkingPrecision->Floor[$MachinePrecision], Method->{Automatic,"SymbolicProcessing"->0}]


SetOptions[Interpolation, InterpolationOrder->3, Method->"Spline"]


$zminSurvey


$zbinsEquiPopu


DirectoryQ[$resultsDirectoryWL]


$paramoptions


$fisherMatWLResultID="FisherMatrix-WL-"<>$surveychosen<>$stepstring<>"--ellmax-"


(* ::Subsection:: *)
(*Some WL settings*)


$IAswitch


Options[setWLBinsSpecifications]


$kdependentGrowth


SigmaLensing[0.9,1]


$pscosmoopts


$paramoptions


(* ::Input::Initialization:: *)
SigmaTheoryLensing[mu_,eta_]:=(mu/2)*(1+eta);


Options[SigmaLensing]


(*SetOptions[SigmaLensing, theoreticalMGFunction->SigmaFunction]*)


SigmaLensing[0.5]


(* ::Subsection:: *)
(*Set Cij options*)


$stepstring


SetOptions[CijGG, shearshearIntegrand->pijIntegrandGGInterpolFunc]


SetOptions[CijGI, shearIAIntegrand->pijIntegrandGIInterpolFunc]


SetOptions[CijII, IAIAIntegrand->pijIntegrandIIInterpolFunc]


(* ::Subsubsection:: *)
(*Initialize Lensing Kernels*)


(*Table[lensKernelGGInterpol[b1,b2,{}],{b1,1,$numZbinsWL},{b2,1,$numZbinsWL}];//AbsoluteTiming*)


(*Table[lensKernelGIInterpol[b1,b2,{}],{b1,1,$numZbinsWL},{b2,1,$numZbinsWL}];//AbsoluteTiming*)


(*Table[lensKernelIIInterpol[b1,b2,{}],{b1,1,$numZbinsWL},{b2,1,$numZbinsWL}];//AbsoluteTiming*)


(* myInterrupt[] *)


$IAswitch


initializeLensKernels[]


(* ::Subsection:: *)
(*Compute Fisher Matrix*)

derivedParams[fisher_, ellnow_, dampedBool_]:=Module[{newfisherID,hubblefisherID,cplfisherID},

Switch[dampedBool,
True,
newfisherID="NewFisherMatrix-WL-"<>$surveychosen<>"-"<>"-"<>$stepstring<>"--ellmax-"<>ToString@(ellnow)<>"_damped";
hubblefisherID="HubbleFisherMatrix-WL-"<>$surveychosen<>"-"<>"-"<>$stepstring<>"--ellmax-"<>ToString@(ellnow)<>"_damped";
cplfisherID="CPLFisherMatrix-WL-"<>$surveychosen<>"-"<>"-"<>$stepstring<>"--ellmax-"<>ToString@(ellnow)<>"_damped";,
False,
newfisherID="NewFisherMatrix-WL-"<>$surveychosen<>"-"<>"-"<>$stepstring<>"--ellmax-"<>ToString@(ellnow);
hubblefisherID="HubbleFisherMatrix-WL-"<>$surveychosen<>"-"<>"-"<>$stepstring<>"--ellmax-"<>ToString@(ellnow);
cplfisherID="CPLFisherMatrix-WL-"<>$surveychosen<>"-"<>"-"<>$stepstring<>"--ellmax-"<>ToString@(ellnow);
];

setNewParamList[];
jacob=computeJacobianMatrix[];
newFisher=jacobianTransform[fisher,jacob];
newFisherFile=$resultsDirectoryWL<>newfisherID<>".txt";
Export[newFisherFile,newFisher,"Table",  "FieldSeparators" -> " "];
newErrTabMarg=errorsTable[$newParamFiducial,oneSigmaErrors[newFisher],$newParamLabel];
errorsMargTableFile=$resultsDirectoryWL<>newfisherID<>".png";
Export[errorsMargTableFile,newErrTabMarg,"PNG"];
newEllipses=GraphicsGrid[ArrayReshape[Table[plotConfRegion[newFisher,tt[[1]],tt[[2]],$newParamFiducial,$newParamLabel,$blueishTones,2.],{tt,proper2Subset[Range@Length@$newParamList]}],{4,4}]];
newEllipsesFile=$resultsDirectoryWL<>newfisherID<>"--Ellipses.png";
Export[newEllipsesFile,newEllipses,"PNG"];
Switch[StringMatchQ["Exp*"][$nameOfModel],
True,
setHubbleParamList[];
jacobHubble=computeHubbleJacobianMatrix[];
HubbleFisher=jacobianTransform[fisher,jacobHubble];
HubbleFisherFile=$resultsDirectoryWL<>hubblefisherID<>".txt";
Export[HubbleFisherFile,HubbleFisher,"Table",  "FieldSeparators" -> " "];
HubbleErrTabMarg=errorsTable[$HubbleParamFiducial,oneSigmaErrors[HubbleFisher],$HubbleParamLabel];
errorsMargTableFile=$resultsDirectoryWL<>hubblefisherID<>".png";
Export[errorsMargTableFile,HubbleErrTabMarg,"PNG"];




setCPLParamList[];
jacobCPL=computeCPLJacobianMatrix[];
CPLFisher=jacobianTransform[fisher,jacobCPL];
CPLFisherFile=$resultsDirectoryWL<>cplfisherID<>".txt";
Export[CPLFisherFile,CPLFisher,"Table",  "FieldSeparators" -> " "];
CPLErrTabMarg=errorsTable[$CPLParamFiducial,oneSigmaErrors[CPLFisher],$CPLParamLabel];
errorsMargTableFile=$resultsDirectoryWL<>cplfisherID<>".png";
Export[errorsMargTableFile,CPLErrTabMarg,"PNG"];,
False,
None;
];
];


(*Dynamic[{"param_a:",Global`proga, "param_b:",Global`progb, "ell:", Global`progell, "time in s:", Global`progt}]*)


Dynamic[{"param_a:",Global`proga, "param_b:",Global`progb, "ell:", Global`progell, "time in s:", Global`progt}]
Do[
ClearAll[$fisherMatrix];
ClearAll[$fisherMatrixdamp];

initialdate=DateString[];
Print["Initial Date: ", initialdate];

ellnow=$elllist[[index]];
Print["ell_max=",ellnow];
setWLBinsSpecifications[ellmax->ellnow];

Print["Computing Fisher Matrix for all parameters..."];
SetOptions[FisherWLSimple, kDamping->False];
totaltime=First@AbsoluteTiming@($fisherMatrix=FisherWLMatrix[]);
Print[ToString@(ellnow)<>"--fisher matrix WL finished"];
Print["Total time spent in sec: "<>ToString[totaltime]];
fisherID=($fisherMatWLResultID<>ToString@(ellnow));
Export[$resultsDirectoryWL<>fisherID<>".txt",$fisherMatrix,"Table"];

errsTable=errorsTable[$paramfidus,oneSigmaErrors[$fisherMatrix],$paramlabels];
errsTabFile=$resultsDirectoryWL<>"errorsTable-margpars-"<>fisherID<>".png";
Export[errsTabFile,errsTable];
derivedParams[$fisherMatrix, ellnow, False];


SetOptions[FisherWLSimple, kDamping->($kmaxHard)];
Print["Computing Fisher Matrix for all parameters..."];
totaltime=First@AbsoluteTiming@($fisherMatrixdamp=FisherWLMatrix[]);
Print[ToString@(ellnow)<>"--fisher matrix WL finished"];
Print["Total time spent in sec: "<>ToString[totaltime]];
fisherIDdamp=($fisherMatWLResultID<>ToString@(ellnow))<>"_damped";
Export[$resultsDirectoryWL<>fisherIDdamp<>".txt",$fisherMatrixdamp,"Table"];

errsTable=errorsTable[$paramfidus,oneSigmaErrors[$fisherMatrixdamp],$paramlabels];
errsTabFile=$resultsDirectoryWL<>"errorsTable-primaryParsMarg-"<>fisherIDdamp<>".png";
Export[errsTabFile,errsTable];
derivedParams[$fisherMatrixdamp, ellnow, True];

Print["Files exported to: "<>$resultsDirectoryWL];

finalDate=DateString[];
Print["Final Date: ", finalDate]
,
{index,1,Length@$elllist}]


$paramoptions
