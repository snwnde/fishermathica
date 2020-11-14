(* ::Package:: *)

(* ::Section:: *)
(*Initialization Cells*)


(* ::Subsection:: *)
(*Set paths to local folders and directories*)


Switch[ToString@RunThrough["whoami",""],
"santiago",
cosmomathicadir="/home/santiago/CosmoProjects/cosmomathica/";
fishertoolsdir="/home/santiago/CosmoProjects/fishermathica/";,
"scasas",
cosmomathicadir="/local/home/scasas/CosmoProjects/cosmomathica/";
fishertoolsdir="/local/home/scasas/CosmoProjects/fishermathica/";
];



With[{p=$Path},ParallelEvaluate[$Path=p]];


AppendTo[$Path, cosmomathicadir];
AppendTo[$Path, fishertoolsdir];


(* ::Subsection:: *)
(*Load packages *)


Get["cosmomathica.m"];
Get["CosmologyFunctions.m"];
Get["FisherTools.m"];
Get["UsefulTools.m"];
Get["SurveySpecifications.m"];
Get["WeakLensingFisherTools.m"];
Get["GenerateInputInfo.m"];
(*Get["sendMailPackage.m"];*)


(* ::Subsection:: *)
(*Load packages in parallel for parallel evaluations*)


If[SameQ[$globalpackageparallel,True],
Print["Evaluating packages in parallel..."];
With[{p=$Path},ParallelEvaluate[$Path=p]];
ParallelEvaluate[Needs["cosmomathica`interface`","cosmomathica.m"]];
ParallelNeeds["CosmologyFunctions`"];
ParallelNeeds["FisherTools`"];
SetOptions[ParallelTable,
DistributedContexts->{"CosmologyFunctions`","FisherTools`", "WeakLensingFisherTools`","Global`"}]; 
 (*This can be set in a package by changing the variable $DistributedContexts*)
DistributeDefinitions[memo];
SetSharedFunction[valueOnMain,setValueOnMain];
Off[NIntegrate::slwcon];
ParallelTable[Off[NIntegrate::slwcon],{i,$KernelCount}];
Off[InterpolatingFunction::dmval];
ParallelEvaluate[Off[InterpolatingFunction::dmval]];
,
Print["Not evaluating packages in parallel. Reload package after setting $globalpackageparallel=True, to evaluate packages in parallel kernels."]
]

