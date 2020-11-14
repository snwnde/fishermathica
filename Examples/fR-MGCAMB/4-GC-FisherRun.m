(* ::Package:: *)

(* ::Section:: *)
(*Fisher matrix blocks and structure*)


$resultsDirectoryGC=mkDirectory[notebookdir<>"/Results-GCsp-"<>"-"<>$stepstring<>"/"]


$zdependderivatives="numerical";


epsizstepdenominator=10000;
$epsilonzstep=N@(1/epsizstepdenominator);   (* epsilon step for numerical derivatives of z-dependent parameters lnDA and lnH for the moment only *)


setNumericalEpsilon[$epsilonstep, epsilonStepForZdependentParameters->$epsilonzstep]


$stencilpoints=5;    (* number of points for the stencil of the numerical derivative of z-dependent quantities *)


SetOptions[dObsPowerDZpi, stencilpoints->$stencilpoints]


$onlyzdeppars=False;


zdepVariablesVector={d8,d1};


FisherBlockBuilder[zdepVariablesVector, 
fisherBlockDerivativeMethod->"FullNumerical", fbsigma8Variables->True, onlyZdependentParameters->$onlyzdeppars]


Print["z-dependent Vector: "]


Print[$zdependDerivVector];


(*GCfisherint$precision=11;
GCfisherint$accuracy=11;
GCfisherint$maxrecursion=1000;*)


(*SetOptions[NIntegrate, PrecisionGoal->Automatic, AccuracyGoal->GCfisherint$accuracy, MinRecursion->2, 
MaxRecursion->GCfisherint$maxrecursion, Method->{Automatic,"SymbolicProcessing"->0},  WorkingPrecision->Floor[$MachinePrecision]]*)


$externalPowerSpectrumInterpolatingFunction


kbinning="kbin_"
interpolstr=interpolParticularString
zderivativesstr=$APoptString<>"_"<>$zdependderivatives<>"_lnHlnDalnfs8lnbs8_"<>(ToString@$stencilpoints)<>"pt_sten"
epszstr=replaceStringPoint[epsizstepdenominator, "epsizd"]
epsshstr=replaceStringPoint[$epsilonstep, "epsish"]
kmaxstrid=replaceStringPoint[NumberForm[$kmaxHard,{3,2}], "kmax"]


$paramoptions


(*myInterrupt[]*)


(* ::Section:: *)
(*Compute derivatives to files*)


$computeDerivatives=False


$epsilonzstep


$zdependDerivVector


$zdependderivatives


$zdependDerivVector


$parampositions


(*krangetest=Exp@Range[Log[hubbleToday[]*$krange[[1]]], Log[hubbleToday[]*$krange[[-1]]], 0.12];
krangetest={0.1,0.2};
If[$computeDerivatives==True,

mucases={0.3};
zderivativesstr="z-k-mu-"<>"P_obs-"<>$APoptString<>"_"<>$zdependderivatives<>"_lnHlnDalnfs8lnbs8Ps";
indi=DateString[];
Print["case: ",  zderivativesstr];
expofilename=$resultsDirectoryGC<>"derivativesFile-krange-"<>zderivativesstr<>".txt";
derivsTable[zderivativesstr]={};

(*Do[
setNumericalEpsilon[$epsilonstep, epsilonStepForZdependentParameters->stepii];
(*Print["krange Length (of fiducial file)"];
Print[Length[$krange]];*)
Print[" step in z-dep: ", $epsilonzstep];
stepstr=replaceStringPoint[$epsilonzstep, "epsz"];
(*
pos1=1;
pos2=-1;
Print[Length[$krange[[pos1;;pos2]]]];
Print@($krange[[pos1]]);
Print@($krange[[pos2]]);*)
*)
Do[
Do[ 
zind=(zaverage@$zbinGlobal)[[zi]];
(*iind=9+5*(zi-1);
rang=Range[4]~Join~Range[iind-4,iind];*)
temptab=ParallelTable[{zind, kki, muu, 
dObsPowerDpi[zind,kki, muu, 1,dlnPdpDerivative->True],
dObsPowerDpi[zind,kki, muu, 2,dlnPdpDerivative->True],
dObsPowerDpi[zind,kki, muu, 3,dlnPdpDerivative->True],
dObsPowerDpi[zind,kki, muu, 4,dlnPdpDerivative->True],
dObsPowerDpi[zind,kki, muu, 5,dlnPdpDerivative->True],
dObsPowerDpi[zind,kki, muu, 6,dlnPdpDerivative->True],
d8[zind,kki,muu],
d1[zind,kki,muu]
},
{kki,krangetest} ];
Print["z value: ", zind];
Print["mu value: ", muu];
derivsTable[zderivativesstr]=AppendTo[derivsTable[zderivativesstr],temptab];
,
{muu, mucases}
];
,
{zi,1,Length@(zaverage@$zbinGlobal)}
];
(*,
{stepii, {10.0^(-5),2*10.0^(-4)}(*{10.0^(-5), 5*10.0^(-5), 1*10.0^(-4), 5*10.0^(-4), 1*10.0^(-3), 5*10.0^(-3), 1*10.0^(-2), 2*10.0^(-2), 5*10.0^(-2), 1*10.0^(-1)}*)}
];*)
endi=DateString[];
Print["time passed: "];
Print[DateDifference[indi,endi, "Minutes"]];
Print["Exporting:  ", expofilename];
Export[expofilename, 
Flatten[derivsTable[zderivativesstr],1], "Table"];
derivfile[zderivativesstr]=Import[expofilename,"Table"];
];*)


zaverage@$zbinGlobal


$paramoptions


(*compaTabs[tab_, tabbench_]:=100*(tab-tabbench)/tabbench*)


(*santiObsCosmDerivs=Table[{dObsPowerDpi[1.4, 0.1, 0.3, pp, dlnPdpDerivative->True], dObsPowerDpi[1.4,0.2, 0.3, pp, dlnPdpDerivative->True]},{pp,1,6}] *)


(*domenicoObsCosmDerivs={{1.6757 ,2.59094},{1.0439, -2.48047},
{-3.35459, -2.91231},
{-0.02749, 0.66669},
{0.10913, 0.11233},
{-0.003, 0.026}
}*)


(*compaTabs[santiObsCosmDerivs, domenicoObsCosmDerivs]*)


(*d8[0.9,0.22,0.44]*)


(*?dObsPowerDpi*)


(*LogLinearPlot[dObsPowerDpi[1, kkki, 0., 1,dlnPdpDerivative->True], {kkki, 0.001, 1}, PlotRange\[Rule]Full]
LogLinearPlot[dObsPowerDpi[1, kkki, 0., 3,dlnPdpDerivative->True], {kkki, 0.001, 1}, PlotRange\[Rule]Full]*)


(*dObsdPi$1=Table[{klk, dObsPowerDpi[1, klk, 0., 1, dlnPdpDerivative->True]}, {klk, logarithmicDivisions[{0.001,1},100]}];
dObsdPi$2=Table[{klk, dObsPowerDpi[1, klk, 0., 2, dlnPdpDerivative->True]}, {klk, logarithmicDivisions[{0.001,1},100]}];
dObsdPi$3=Table[{klk, dObsPowerDpi[1, klk, 0., 3, dlnPdpDerivative->True]}, {klk, logarithmicDivisions[{0.001,1},100]}];*)


(*Export["dObsdPi$3.txt",dObsdPi$3]
Export["dObsdPi$2.txt",dObsdPi$2]
Export["dObsdPi$1.txt",dObsdPi$1]*)


(*MYdObsdPi$1 = Import["MYdObsdPi$1.txt", "Table"];
MYdObsdPi$2 = Import["MYdObsdPi$2.txt", "Table"];
MYdObsdPi$3 = Import["MYdObsdPi$3.txt", "Table"];*)


(*MYdObsdPi$1//Dimensions*)


(*ListLogLinearPlot[{dObsdPi$1,MYdObsdPi$1} , Joined\[Rule]True, PlotLegends\[Rule]Automatic]*)


?powerSpectrum


covMatGC[zv_, kv_, muv_]:=Block[{pref=2*(2 Pi)^3, veffs, pobs, cov},
veffs= volumeEffective[zv, kv, muv]*volumeSurvey[zv, $zbinGlobal];
pobs=observedPowerSpectrum[zv,kv,muv];
cov = (pref/veffs)*pobs^2 * (1/kv)^3;
Return[cov]
]


covMatGC[0.9,0.001, 0.99]


?volumeSurvey


ztestaa=0.9;
mutestaa=0.1;


(*LogLogPlot[{observedPowerSpectrum[ztestaa,kiki, mutestaa]-Sqrt[covMatGC[ztestaa,kiki, mutestaa]], observedPowerSpectrum[ztestaa,kiki, mutestaa],
observedPowerSpectrum[ztestaa,kiki, mutestaa]+Sqrt[covMatGC[ztestaa,kiki, mutestaa]]}, {kiki, 0.001, 2.}, PlotLegends->Automatic, PlotRange->Full, PlotStyle->{Dashed, Thick, Dashed}]*)


Options[volumeEffective]


LogLogPlot[{kiki*observedPowerSpectrum[ztestaa,kiki, mutestaa]}, {kiki, 0.01, 0.8}, PlotLegends->Automatic(*, PlotRange\[Rule]{60000, 65000}*)]


(*LogLinearPlot[{kiki*Sqrt[covMatGC[ztestaa,kiki, 0.9]]}, {kiki, 0.001, 0.8}, PlotLegends->Automatic(*, PlotRange\[Rule]{60000, 65000}*)]*)


(*LogLogPlot[{Sqrt[covMatGC[ztestaa,kiki, 0.9]]/observedPowerSpectrum[ztestaa,kiki, 0.9]}, {kiki, 0.001, 3}, PlotLegends->Automatic(*, PlotRange\[Rule]{60000, 65000}*)]*)





(* ::Section:: *)
(*Set Integration Options*)


Options[FisherIntegration]


Options[NIntegrate]


(** Original settings:   SetOptions[NIntegrate, Method->{Automatic,"SymbolicProcessing"->0}, AccuracyGoal->Infinity, PrecisionGoal->Automatic]  **)
(*SetOptions[NIntegrate, Method->{Automatic,"SymbolicProcessing"->0}, AccuracyGoal->12, PrecisionGoal->Automatic, MaxRecursion->1000, MinRecursion\[Rule]2]*)


SetOptions[FisherIntegration,functionNIntegrate->NIntegrate]


(*Options[NIntegrateInterpolatingFunction]*)


(*(zaverage@$zbinGlobal)[[1;;3]]*)


(*zderivativesstr*)


(*$stencilpoints*)


(*$zdependderivatives*)


(*$resultsDirectoryGC*)


(*lnDAderivInterpol["krange-5pt"]=Interpolation[{{#1,#2,#3,#4},#5}&@@@derivfile[zderivativesstr], InterpolationOrder->3, Method->"Spline"];*)


(*lnHderivInterpol["krange-5pt"]=Interpolation[{{#1,#2,#3,#4},#6}&@@@derivfile[zderivativesstr], InterpolationOrder->3, Method->"Spline"];*)


(*lnHderivInterpol[3]*)


(*Export[$resultsDirectoryGC<>"interpolation-lnH-finerK-smallsteps-numerical-5pts.mx",lnHderivInterpol["finer-5pt"], "MX"]*)


(*Export[$resultsDirectoryGC<>"interpolation-lnDa-finerK-smallsteps-numerical-5pts.mx",lnDAderivInterpol["finer-5pt"], "MX"]*)


(*If[$justDerivatives==True, myInterrupt[]]*)


(*$paramoptions*)


(*$zdependDerivPositions*)


(*$numZindependentPars*)


(*$numZdependentPars*)


(*zdepstring=StringJoin[(ToString[#])&/@zdepVariablesVector]*)


(*Print[zdepstring]*)


(*$fisherCosmoFullChainBlock[0.95,k,mu]*)


?$fisherCosmoFullChainBlock


$paramoptions


(* ::Section:: *)
(*Run the GC Fisher Matrix Calculation*)


zaverage@$zbinGlobal


(*FisherMatrixGC[1., 1,2,$fisherCosmoFullChainBlock[1.0, k, mu]]*)


$fisherMatResultID="IST_GC_Fisher-WP6_Casas_NLrecipe-vT5-"<>kbinning<>"-"<>interpolstr<>"-"<>zderivativesstr<>"-"<>epszstr<>"-"<>epsshstr<>"-"<>kmaxstrid<>"-"<>$stepstring


Print["$fisherMatResultID:  "<>$fisherMatResultID]


?$fisherCosmoFullChainBlock


($fisherCosmoParsBlock[#,k,mu]&/@(zaverage@$zbinGlobal);)//AbsoluteTiming


$parampositions


MatrixPlot[$fisherCosmoParsBlock[1.2,0.1,0.99]]


(*Interrupt[]*)


(*FisherIntTest[$fisherCosmoParsBlock, 1.0, cosmoBlockPkderivatives\[Rule]True]*)


FisherMatrixGC[1.0,1,1,$fisherCosmoParsBlock[1.0,k,mu][[1,1]]]


If[$justAnalysys==False,
totaltime=FisherMatrixGCCalculation[$fisherCosmoParsBlock, $fisherDiagonalZTermsBlock, $fisherCrossTermBlock, 
resultsDirectory->$resultsDirectoryGC,clearFisherIntegrationCache->False, fisherMatrixExportNameID->$fisherMatResultID, 
storeTemporaryResults->True, compressTemporaryResults->False, createProgressDialogMonitor->False]
]


(* ::Input:: *)
(**)


Total@totaltime


(Total@totaltime)/60


$resultsDirectoryGC


DirectoryQ[$resultsDirectoryGC]


FisherFinalAssembly[exportFisherMatrix->True,compressExportMatrix->False, exportParametersUsed->True]


(*If[$newrecipelinear==True,
sp={5,6};
lsp=Length@sp;
lspi=-(lsp+1);
complementlin=Complement[Range[$numTotalEntries], sp];
$fisherMatrix=$fisherMatrix[[complementlin,complementlin]];
Export[$resultsDirectoryGC<>$fisherMatResultID<>"-linrecipe.txt",$fisherMatrix,"Table"];
$numZindependentPars=$numZindependentPars-lsp;
$paramfidus=$paramfidus[[1;;lspi]];
$paramlabels=$paramlabels[[1;;lspi]];
]*)


posmat=PositiveDefiniteMatrixQ[$fisherMatrix]


Print["Positive definite matrix: "<>ToString[posmat]]


$numZindependentPars


$numZdependentPars


$paramoptions


$fisherMatrix[[1;;$numZindependentPars,1;;$numZindependentPars]]


oneSigmaErrors[$fisherMatrix[[1;;$numZindependentPars,1;;$numZindependentPars]]]


shapeErrors=oneSigmaErrors[$fisherMatrix][[1;;$numZindependentPars]]


$paramfidus


$paramlabels


errTab=errorsTableTxt[$paramfidus, shapeErrors, $paramlabels]


new={1.4,2.6,4.9,0.29,1.2,0.91}


old={1.29,2.20,5.13,0.25,1.38,2.09}


Abs[(1-old/new)]*100


$paramoptions


$numZdependentPars


zdeperrorslist=oneSigmaErrors[$fisherMatrix][[$numZindependentPars+1;;-1]];


$zdependDerivPositions


errorszdeptab=Table[{NumberForm[(zaverage@$zbinGlobal)[[ii]],{3,2}]}~Join~(NumberForm[#,{6,5}]&/@(zdeperrorslist[[1+$numZdependentPars*(ii-1);;$numZdependentPars+$numZdependentPars(ii-1)]])), 
{ii, 1, Length@(zaverage@$zbinGlobal)}]


exporterrorszdeptab=Prepend[errorszdeptab, {"z", "lnbs8", "Ps"}]


PrependTo[exporterrorszdeptab,(NumberForm[#,{6,5}]&/@shapeErrors) ]


PrependTo[exporterrorszdeptab,(NumberForm[#,{6,5}]&/@$paramlabels) ]


$resultsDirectoryGC


errTabFile=$resultsDirectoryGC<>"marg-shape-errors--"<>$fisherMatResultID<>".txt";
allerrsTabFile=$resultsDirectoryGC<>"All-shape-zdep-errors--"<>$fisherMatResultID<>".txt";


Export[allerrsTabFile,exporterrorszdeptab, "Table"]


Export[errTabFile,errTab, "Table"]


nuisanceParsList={indexdlnbs8, indexdPs};


fisherFullNumDerivMargBias=fullFisherPostTransformation[$fisherMatrix,exportTableName->$resultsDirectoryGC<>$fisherMatResultID<>"-margnuisance", marginalizedElementsIndex->nuisanceParsList,
performReduceOperations->"Marginalize", repairPositivity->False]


shapeErrors


oneSigmaErrors[fisherFullNumDerivMargBias]
