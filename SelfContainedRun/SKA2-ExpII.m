(* ::Package:: *)

(* ::Section:: *)
(*Initialization Cells*)


(* ::Subsection:: *)
(*Set paths to local folders and directories*)


SetOptions[$Output, FormatType->OutputForm];
notebookdir=If[aa=DirectoryName@$InputFileName; aa=="", NotebookDirectory[], aa]


SetDirectory[notebookdir]


$globalpackageparallel=True


Get["Initialization.m"]


(* ::Section:: *)
(*Run type specifications*)


(*$MinPrecision=7;*)
$debugPlots=False;
$debugPrint=False;


$dateshort=DateString["YearShort"]<>DateString["MonthNameShort"]<>"d"<>DateString["DayShort"]

$nameOfModel="Exp II"

Get["0-Runtype.m"]
(* 
epsizstepdenominator=10000;
$epsilonzstep=N@(1/epsizstepdenominator);   epsilon step for numerical derivatives of z-dependent parameters lnDA and lnH for the moment only *)


(* ::Section:: *)
(*Run type strings*)


loadInputFilesNotebook="2-LoadInputFiles-SteM"<>".m"


FileExistsQ[notebookdir<>loadInputFilesNotebook]


inputfilename=StringReplace[loadInputFilesNotebook,{".m"->""}]


setPScosmoOptions[lcdmBool->False,linearBool->True, kdependentGrowth->False];


(* ::Section:: *)
(*Cosmological Parameters, choose fiducial values.*)


Get["1-Parameters.m"]


(* ::Section:: *)
(*Load cosmological functions from files*)


$deleteInterpolationFiles=False


Get[loadInputFilesNotebook]


(*Get["2-LoadInputFiles-SteM.m"]*)


(* ::Section:: *)
(*Load survey specifications and set other survey parameters*)


$pscosmoopts


SetOptions[$Output, FormatType->OutputForm];
notebookdir=If[aa=DirectoryName@$InputFileName; aa=="", NotebookDirectory[], aa]
SetDirectory[notebookdir]
$surveychosen="SKA2";
Get["3-Specifications.m"]


zaverage@$zbinGlobal


zmean=$zmeansurvey


(* ::Section:: *)
(*Fisher Matrix Run*)


(*myInterrupt[]*)


Get["4-FisherRun-GC.m"]
Get["5-FisherRun-WL.m"]


(* ::Section:: *)
(*Send results by mail*)


finaldate=DateString[]


sessionTime=SessionTime[]


(*sendMail["Fisher Matrix computation: fisherID: "<>$fisherMatResultID,"Fisher matrix computation finished: \n Computing Fisher time in seconds: "<>ToString@Total[totaltime]<>
"\n initial notebook date: "<>ToString[999]<>"\n final notebook date: "<>ToString[finaldate]<>"\n Kernel session time in sec: "<>
ToString@sessionTime<>"\n\n\n"<>"\n Positivity of matrix: "<>ToString[posmat]<>
"\n Largest Eigenvalue of matrix: "<>ToString[FortranForm@maxeig]<>"\nFully Marginalized errors on shape parameters: "<>"\n"<>ToString[$paramlabels]<>"\n"<>ToString[shapeErrors]<>"\n",
"casas@thphys.uni-heidelberg.de",{($resultsDirectoryGC<>$fisherMatResultID<>".txt"), errTabFile }]*)


(*CloseKernels[];
Pause[5];*)


Quit[];
