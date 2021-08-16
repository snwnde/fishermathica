(* ::Package:: *)

BeginPackage["GenerationParameters`"]
findInputFilesDir::usage="Give the input files directory."
paramsFolderFunc::usage="Give the input file folders per parameter and epsilon."
$folderpars::usage="List of input file parameters."

epsstrfun::usage="Give string of epsilon."
$epslist::usage="List of epsions."
$epslistPM::usage="List of epsilons (minus and plus)."
epsstrRule::usage="Epsilon string expression rule."
plusminstr::usage="plus and munius string + and -"

Needs["CosmologyFunctions`"];
Needs["UsefulTools`"];
Needs["SurveySpecifications`"];
Needs["FunctionApproximations`"];
Needs["Quintessence`"];

EndPackage[]

BeginPackage["GenerationParameters`", {"CosmologyFunctions`","UsefulTools`","SurveySpecifications`","FunctionApproximations`","Quintessence`","FisherTools`"}]
Begin["`Private`"]
(* DirectoryName@$InputFileName
notebookdir=If[aa=DirectoryName@$InputFileName; aa=="", NotebookDirectory[], aa]
SetDirectory[notebookdir]; *)
machine=Global`machine;

DirectoryName@$InputFileName
geneParamPacDir=If[aa=DirectoryName@$InputFileName; aa=="", NotebookDirectory[], aa]
SetDirectory[geneParamPacDir];

initialdate=DateString[]

(*$MinPrecision=7;*)
$debugPlots=False;
$debugPrint=False;


dateshort=DateString["YearShort"]<>DateString["MonthNameShort"]<>"d"<>DateString["DayShort"]


parametertestIndex=0;
$Mnuchoice=0.15;


Options[dObsPowerDZpi]


(*$epsilonstep=0.003;   (* epsilon step for numerical derivatives of shape parameters, this is connected to the +/-epsilon input files *)
*)
$folderpars={"fiducial"}~Join~$paramlabels
plusminstr={"mn","pl"};
epsstrfun[eps_]:=scientificfortranform[eps, "E", "p", "eps_"]

$epslist = {0.00625, 0.0125, 0.01875, 0.025, 0.0375, 0.05, 0.10}; 
$epslistPM = Sort[(-1*$epslist)~Join~$epslist~Join~{0.}]
epsstrRule=Thread[Rule[$epslist, (epsstrfun[#]&/@$epslist)]]

paramsFolderFunc[param_, epsilon_]:=Block[{pm, string, parstring, epsstr, epsi}, 
                                          pm=Which[Sign[epsilon]==-1, "_"<>plusminstr[[1]], Sign[epsilon]==1, "_"<>plusminstr[[2]], Sign[epsilon]==0, ""];
                                          epsi=If[Chop[epsilon, 10^-6]==0., 0, epsilon];

                                          If[param=="gamma" && Abs[epsi]>=0.10, Return[None]];
                                          (* If[param=="gamma" && epsi==-0.10, Return[None]]; *)

                                          If[param=="fiducial" && epsi!=0, Return[None]];
                                          If[param!="fiducial" && epsi==0, Return[None]];

                                          If[param=="fiducial", pm=""];
                                          epsstr=epsstrfun[Abs[epsi]];
										                      string=(param<>pm<>"_"<>epsstr);						
                                          Return[string]
                                          ];




findInputFilesDir[modelname_]:=Switch[modelname,
"Exp II",
Return[geneParamPacDir<>"../InputFiles/AlphaAttractorExpII/"<>"/"],
"Exp I",
Return[geneParamPacDir<>"../InputFiles/AlphaAttractorExpI/"<>"/"],
"lcdm",
Return[geneParamPacDir<>"../InputFiles/LCDM/"<>"/"]
]



(* ::Section::Closed:: *)
(*End of package*)


End[]
EndPackage[]
