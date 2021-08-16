(* ::Package:: *)

(* ::Section:: *)
(*Fisher matrix blocks and structure*)


$resultsDirectoryGC=mkDirectory[notebookdir<>"/Results-GC-"<>$dateshort<>"-"<>$stepstring<>"-"<>$surveychosen<>"/"]

(*Print["d1[1.8,0.1,0.88]="<>ToString[d1[1.8,0.1,0.88]]];*)


$zdependderivatives


kbinning="kbin_"<>fileinputkbinning
interpolstr="interp_"<>interpInputTabs<>"_ord_"<>ToString[fiIO]
zderivativesstr=$APoptString<>"_"<>$zdependderivatives<>"_lnHlnDalnfs8lnbs8_"<>(ToString@$stencilpoints)<>"pt_sten"
epszstr=replaceStringPoint[epsizstepdenominator, "epsizd"]
epsshstr=replaceStringPoint[$epsilonstep, "epsish"]
kmaxstrid=replaceStringPoint[NumberForm[$kmaxHard,{3,2}], "kmax"]





(* ::Section:: *)
(*Compute derivatives to files*)


powerSpectrum[0.5, 0.1]


observedPowerSpectrum[0.8, 0.1, 0.9]


fGrowthRate[0.5]


Growth[0.4]





$computeDerivatives=False


$epsilonzstep


(*krangetest=Exp@Range[Log[hubbleToday[]*$krange[[1]]], Log[hubbleToday[]*$krange[[-1]]], 0.12];*)


(*krangetest=$krange;*)


$zdependDerivVector


$zdependderivatives


mucases={0.5,1.0};
zderivativesstr="z-k-mu-"<>"P_obs-"<>$APoptString<>"_"<>$zdependderivatives<>"_lnHlnDalnfs8lnbs8Ps";
indi=DateString[];
Print["case: ",  zderivativesstr];
expofilename=$resultsDirectoryGC<>"derivativesFile-krange-"<>zderivativesstr<>".txt";
derivsTable[zderivativesstr]={};
If[$computeDerivatives==True,
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
d4[zind,kki,muu], 
d5[zind,kki,muu],
d7[zind,kki,muu],
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
];


expofilename


Options[FisherIntegration]


Options[NIntegrate]


SetOptions[NIntegrate, Method->{Automatic,"SymbolicProcessing"->0}, AccuracyGoal->Infinity, PrecisionGoal->Automatic]


SetOptions[FisherIntegration,functionNIntegrate->NIntegrateInterpolatingFunction]



(* ::Section:: *)
(*Run the GC Fisher Matrix Calculation*)


$fisherMatResultID="GC_Fisher-"<>kbinning<>"-"<>interpolstr<>"-"<>zderivativesstr<>"-"<>epszstr<>"-"<>epsshstr<>"-"<>kmaxstrid<>"-"<>$stepstring


Print["$fisherMatResultID:  "<>$fisherMatResultID]


(*Interrupt[]*)


Options@FisherMatrixGCCalculation


If[$justAnalysys==False,
totaltime=FisherMatrixGCCalculation[$fisherCosmoParsBlock, $fisherDiagonalZTermsBlock , $fisherCrossTermBlock, 
resultsDirectory->$resultsDirectoryGC,clearFisherIntegrationCache->False, fisherMatrixExportNameID->$fisherMatResultID, 
storeTemporaryResults->True, compressTemporaryResults->False, createProgressDialogMonitor->False]
]


(* ::Input:: *)
(**)


(Total@totaltime)/60


$resultsDirectoryGC


DirectoryQ[$resultsDirectoryGC]


FisherFinalAssembly[exportFisherMatrix->True,compressExportMatrix->False, exportParametersUsed->True]


If[$newrecipelinear==True,
sp={5,6};
lsp=Length@sp;
lspi=-(lsp+1);
complementlin=Complement[Range[$numTotalEntries], sp];
$fisherMatrix=$fisherMatrix[[complementlin,complementlin]];
Export[$resultsDirectoryGC<>$fisherMatResultID<>"-linrecipe.txt",$fisherMatrix,"Table"];
$numZindependentPars=$numZindependentPars-lsp;
$paramfidus=$paramfidus[[1;;lspi]];
$paramlabels=$paramlabels[[1;;lspi]];
]


posmat=PositiveDefiniteMatrixQ[$fisherMatrix]


Print["Positive definite matrix: "]


Print[posmat]


$numZindependentPars


$numZdependentPars


$paramoptions


shapeErrors=oneSigmaErrors[$fisherMatrix][[1;;$numZindependentPars]]


errTab=errorsTableTxt[$paramfidus, shapeErrors, $paramlabels]


$numZdependentPars


zdeperrorslist=oneSigmaErrors[$fisherMatrix][[$numZindependentPars+1;;-1]];


errorszdeptab=Table[{NumberForm[(zaverage@$zbinGlobal)[[ii]],{3,2}]}~Join~(NumberForm[#,{6,5}]&/@(zdeperrorslist[[1+$numZdependentPars*(ii-1);;$numZdependentPars+$numZdependentPars(ii-1)]])), 
{ii, 1, Length@(zaverage@$zbinGlobal)}]


exporterrorszdeptab=Prepend[errorszdeptab, {"z", "lnD", "lnH", "lnfs8", "lnbs8", "Ps"}]


PrependTo[exporterrorszdeptab,(NumberForm[#,{6,5}]&/@shapeErrors) ]


PrependTo[exporterrorszdeptab,(NumberForm[#,{6,5}]&/@$paramlabels) ]


errTabFile=$resultsDirectoryGC<>"marg-shape-errors--"<>$fisherMatResultID<>".txt";
allerrsTabFile=$resultsDirectoryGC<>"All-shape-zdep-errors--"<>$fisherMatResultID<>".txt";


Export[allerrsTabFile,exporterrorszdeptab, "Table"]


Export[errTabFile,errTab, "Table"]


Eigenvalues@$fisherMatrix


fisherFullNumDerivMargBias=fullFisherPostTransformation[$fisherMatrix,exportTableName->$resultsDirectoryGC<>$fisherMatResultID<>"-"<>"-margbias", marginalizedElementsIndex->{indexdlnb, indexdPs},
performReduceOperations->"Marginalize",repairPositivity->False]


errTabMarg=errorsTable[$paramfidus,oneSigmaErrors[fisherFullNumDerivMargBias],parlabels]
errorsMargTableFile=$resultsDirectoryGC<>"errorsMargTable-"<>$fisherMatResultID<>".png";
Export[errorsMargTableFile,errTabMarg,"PNG"]

plotConfRegion[fish_,a_,b_,par_,label_,colorslist_,thick_]:=Graphics[draw123Ellipses[marginalized2Matrix[fish,a,b],a,b,par,colorslist,lineThickness->thick],Frame->True,FrameLabel->styFunc[label[[{a,b}]]],AspectRatio->1.,FrameTicksStyle->Directive[Black,16],ImageSize->{400,400}]
gridConfidenceRegions=GraphicsGrid[ArrayReshape[Table[plotConfRegion[fisherFullNumDerivMargBias,tt[[1]],tt[[2]],$paramfidus,$paramlabels,$blueishTones,2.],{tt,proper2Subset[Range@Length@$paramoptions]}],{4,4}]]

ellipsesFile=$resultsDirectoryGC<>"ellipses-"<>$fisherMatResultID<>".png";
Export[ellipsesFile,gridConfidenceRegions,"PNG"]


(* ::Section:: *)
(*Other Stuff*)


setNewParamList[]


jacob=computeJacobianMatrix[]


newFisher=jacobianTransform[fisherFullNumDerivMargBias,jacob];


newFisherFile=$resultsDirectoryGC<>"newFisher"<>$stepstring<>".txt";


Export[newFisherFile,newFisher,"Table",  "FieldSeparators" -> " "]


newErrTabMarg=errorsTable[$newParamFiducial,oneSigmaErrors[newFisher],$newParamLabel]
errorsMargTableFile=$resultsDirectoryGC<>"newErrorsMargTable.png";
Export[errorsMargTableFile,newErrTabMarg,"PNG"];


newEllipses=GraphicsGrid[ArrayReshape[Table[plotConfRegion[newFisher,tt[[1]],tt[[2]],$newParamFiducial,$newParamLabel,$blueishTones,2.],{tt,proper2Subset[Range@Length@$newParamList]}],{4,4}]]
newEllipsesFile=$resultsDirectoryGC<>"newEllipses.png";
Export[newEllipsesFile,newEllipses,"PNG"]


Switch[StringMatchQ["Exp*"][$nameOfModel],
True,
setHubbleParamList[];
jacobHubble=computeHubbleJacobianMatrix[];
HubbleFisher=jacobianTransform[fisherFullNumDerivMargBias,jacobHubble];

HubbleFisherFile=$resultsDirectoryGC<>"HubbleFisher"<>$stepstring<>".txt";
Export[HubbleFisherFile,HubbleFisher,"Table",  "FieldSeparators" -> " "];

HubbleErrTabMarg=errorsTable[$HubbleParamFiducial,oneSigmaErrors[HubbleFisher],$HubbleParamLabel];
errorsMargTableFile=$resultsDirectoryGC<>"HubbleErrorsMargTable.png";
Export[errorsMargTableFile,HubbleErrTabMarg,"PNG"];

setCPLParamList[];
jacobCPL=computeCPLJacobianMatrix[];
CPLFisher=jacobianTransform[fisherFullNumDerivMargBias,jacobCPL];

CPLFisherFile=$resultsDirectoryGC<>"CPLFisher"<>$stepstring<>".txt";
Export[CPLFisherFile,CPLFisher,"Table",  "FieldSeparators" -> " "];

CPLErrTabMarg=errorsTable[$CPLParamFiducial,oneSigmaErrors[CPLFisher],$CPLParamLabel];
errorsMargTableFile=$resultsDirectoryGC<>"CPLErrorsMargTable.png";
Export[errorsMargTableFile,CPLErrTabMarg,"PNG"];,
False,
None;
];
