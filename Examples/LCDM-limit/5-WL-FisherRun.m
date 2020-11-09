(* ::Package:: *)

(* ::Subsection:: *)
(*Initialization and Specifications*)


$resultsDirectoryWL=mkDirectory[notebookdir<>"/Results-WL-"<>$dateshort<>"-"<>$stepstring<>"/"]


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


myInterrupt[]


$IAswitch


initializeLensKernels[]


(* ::Subsection:: *)
(*Compute Fisher Matrix*)


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

Print["Files exported to: "<>$resultsDirectoryWL];

finalDate=DateString[];
Print["Final Date: ", finalDate]
,
{index,1,Length@$elllist}]


$paramoptions
