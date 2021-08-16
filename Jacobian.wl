(* ::Package:: *)

BeginPackage["Jacobian`"]


(* ::Section::Closed:: *)
(*Contextualization*)


w0CPL::usage="CPL parameter w0"
waCPL::usage="CPL parameter wa"
M2::usage="Quintessential parameter M2"
r::usage="Tensor-to-scalar ratio r"
w0dpi::usage="Derivative with respect to cosmological parameters of CPL parameter w0"
wadpi::usage="Derivative with respect to cosmological parameters of CPL parameter wa"
M2dpi::usage="Derivative with respect to cosmological parameters of model parameter M^2"
rdpi::usage="Derivative with respect to cosmological parameters of tensor-to-scalar ratio r"
hubbledpi::usage="Derivative with respect to cosmological parameters of hubble rate today"


setNewParamList::usage="Set new parameter list, fiducial values and labels"
$newParamList::usage="New parameter list"
$newParamFiducial::usage="New parameter fiducial values"
$newParamLabel::usage="New parameter labels"


setCPLParamList::usage="Set CPL parameter list, fiducial values and labels"
$CPLParamList::usage="CPL parameter list"
$CPLParamFiducial::usage="CPL parameter fiducial values"
$CPLParamLabel::usage="CPL parameter labels"


setHubbleParamList::usage="Set Hubble parameter list, fiducial values and labels"
$HubbleParamList::usage="Hubble parameter list"
$HubbleParamFiducial::usage="Hubble parameter fiducial values"
$HubbleParamLabel::usage="Hubble parameter labels"


computeJacobianMatrix::usage="Jacobian matrix computer"
computeCPLJacobianMatrix::usage="CPL Jacobian matrix computer"
computeHubbleJacobianMatrix::usage="Hubble Jacobian matrix computer"


Needs["CosmologyFunctions`"];
Needs["UsefulTools`"];
Needs["SurveySpecifications`"];
Needs["FunctionApproximations`"];
Needs["Quintessence`"];
EndPackage[]


BeginPackage["Jacobian`", {"CosmologyFunctions`","UsefulTools`","SurveySpecifications`","FunctionApproximations`","Quintessence`","FisherTools`"}]
Begin["`Private`"]


(* ::Section::Closed:: *)
(*CPL Parameters*)


Options[w0CPL]=$paramoptions;
w0CPL[opts:OptionsPattern[]]:=Module[{paropts},
paropts=complementParamValues[{opts},w0CPL,returnList->"Full", filterCosmoPars->True];
wDE[NfOfz[0],paropts]
];


Options[waCPL]=$paramoptions;
waCPL[opts:OptionsPattern[]]:=Module[{paropts,waValue},
paropts=complementParamValues[{opts},waCPL,returnList->"Full", filterCosmoPars->True];
waValue=D[wDE[NfOfz[z],paropts],z]/.z->0
];


(* ::Section::Closed:: *)
(*Parameter derivatives*)


Options[M2dpi]={setEpsilon->$epsilonstep,transformDerivativeFunction->"None"}~Join~$paramoptions;

M2dpi[indexA_?IntegerQ,deropts:OptionsPattern[]]:=If[indexA<=Length[$paramoptions],
Block[{extpowopts,optplus,optminus,plusfunc,minusfunc,parfids,epsil=OptionValue[setEpsilon],step,transf=OptionValue[transformDerivativeFunction]},
{step,optplus,optminus,parfids}=numericalParamsDerivative[M2dpi,{deropts},List[indexA],epsil, transformDerivativeFunction->transf];
plusfunc=M2Gen[optplus];
minusfunc=M2Gen[optminus];
(plusfunc-minusfunc)/(step)],
Print["Option index "<>ToString[indexA]<>" not valid"];Abort[]
]


Options[rdpi]={setEpsilon->$epsilonstep,transformDerivativeFunction->"None"}~Join~$paramoptions;

rdpi[indexA_?IntegerQ,deropts:OptionsPattern[]]:=If[indexA<=Length[$paramoptions],
Block[{extpowopts,optplus,optminus,plusfunc,minusfunc,parfids,epsil=OptionValue[setEpsilon],step,transf=OptionValue[transformDerivativeFunction]},
{step,optplus,optminus,parfids}=numericalParamsDerivative[M2dpi,{deropts},List[indexA],epsil, transformDerivativeFunction->transf];
plusfunc=rGen[optplus];
minusfunc=rGen[optminus];
(plusfunc-minusfunc)/(step)],
Print["Option index "<>ToString[indexA]<>" not valid"];Abort[]
]


Options[w0dpi]={setEpsilon->$epsilonstep,transformDerivativeFunction->"None"}~Join~$paramoptions;

w0dpi[indexA_?IntegerQ,deropts:OptionsPattern[]]:=If[indexA<=Length[$paramoptions],
Block[{extpowopts,optplus,optminus,plusfunc,minusfunc,parfids,epsil=OptionValue[setEpsilon],step,transf=OptionValue[transformDerivativeFunction]},
{step,optplus,optminus,parfids}=numericalParamsDerivative[w0dpi,{deropts},List[indexA],epsil, transformDerivativeFunction->transf];
plusfunc=w0CPL[optplus];
minusfunc=w0CPL[optminus];
(plusfunc-minusfunc)/(step)],
Print["Option index "<>ToString[indexA]<>" not valid"];Abort[]
]


Options[wadpi]={setEpsilon->$epsilonstep,transformDerivativeFunction->"None"}~Join~$paramoptions;

wadpi[indexA_?IntegerQ,deropts:OptionsPattern[]]:=If[indexA<=Length[$paramoptions],
Block[{extpowopts,optplus,optminus,plusfunc,minusfunc,parfids,epsil=OptionValue[setEpsilon],step,transf=OptionValue[transformDerivativeFunction]},
{step,optplus,optminus,parfids}=numericalParamsDerivative[wadpi,{deropts},List[indexA],epsil, transformDerivativeFunction->transf];
plusfunc=waCPL[optplus];
minusfunc=waCPL[optminus];
(plusfunc-minusfunc)/(step)],
Print["Option index "<>ToString[indexA]<>" not valid"];Abort[]
]


Options[hubbledpi]={setEpsilon->$epsilonstep,transformDerivativeFunction->"None"}~Join~$paramoptions;

hubbledpi[indexA_?IntegerQ,deropts:OptionsPattern[]]:=If[indexA<=Length[$paramoptions],
Block[{extpowopts,optplus,optminus,plusfunc,minusfunc,parfids,epsil=OptionValue[setEpsilon],step,transf=OptionValue[transformDerivativeFunction]},
{step,optplus,optminus,parfids}=numericalParamsDerivative[hubbledpi,{deropts},List[indexA],epsil, transformDerivativeFunction->transf];
plusfunc=hubbleToday[optplus];
minusfunc=hubbleToday[optminus];
(plusfunc-minusfunc)/(step)],
Print["Option index "<>ToString[indexA]<>" not valid"];Abort[]
]


(* ::Section:: *)
(*Jacobian matrix builder*)


(* ::Subsection::Closed:: *)
(*For M^2 and r*)


Options[setNewParamList]=$paramoptions;
setNewParamList[opts:OptionsPattern[]]:=Module[{paropts},
paropts=complementParamValues[{opts},setNewParamList,returnList->"Full", filterCosmoPars->True];

$newParamList=Module[{nsSelect,AsSelect,nsPos,AsPos,newlist},
nsSelect=KeySelect[paropts,StringMatchQ[ToString[#],"*"<>SymbolName[ns]<>"*"]&];
AsSelect=KeySelect[paropts,StringMatchQ[ToString[#],"*"<>SymbolName[As]<>"*"]&];
nsPos=First@First@Position[paropts,First@Keys@nsSelect];
AsPos=First@First@Position[paropts,First@Keys@AsSelect];
newlist=paropts;
newlist[[nsPos]]=M2E10->M2Gen[paropts]*10^10;
newlist[[AsPos]]=r->rGen[paropts];
newlist];
$newParamFiducial=$newParamList[[All,2]];
$newParamLabel=SymbolName/@$newParamList[[All,1]];
];

Options[computeJacobianMatrix]=$paramoptions;
computeJacobianMatrix[opts:OptionsPattern[]]:=Module[{paropts,newTOoldJacobianMatrix,oldTOnewJacobianMatrix},
paropts=complementParamValues[{opts},computeJacobianMatrix,returnList->"Full", filterCosmoPars->True];

newTOoldJacobianMatrix=Module[{M2Pos,rPos,nsSelect,AsSelect},
nsSelect=KeySelect[paropts,StringMatchQ[ToString[#],"*"<>SymbolName[ns]<>"*"]&];
AsSelect=KeySelect[paropts,StringMatchQ[ToString[#],"*"<>SymbolName[As]<>"*"]&];
M2Pos=First@First@Position[paropts,First@Keys@nsSelect];
rPos=First@First@Position[paropts,First@Keys@AsSelect];

Table[KroneckerDelta[i,j](1-KroneckerDelta[i,M2Pos]-KroneckerDelta[i,rPos])
+10^10*M2dpi[j]*KroneckerDelta[i,M2Pos]+rdpi[j]*KroneckerDelta[i,rPos]
,{i,1,Length@paropts},{j,1,Length@paropts}]
];

oldTOnewJacobianMatrix=Inverse[newTOoldJacobianMatrix]
];


(* ::Subsection:: *)
(*For w0 and wa*)


Options[setCPLParamList]=$paramoptions;
setCPLParamList[opts:OptionsPattern[]]:=Module[{paropts},
paropts=complementParamValues[{opts},setCPLParamList,returnList->"Full", filterCosmoPars->True];

$CPLParamList=Module[{OmegamSelect,gammaSelect,OmegamPos,gammaPos,newlist},
OmegamSelect=KeySelect[paropts,StringMatchQ[ToString[#],"*"<>SymbolName[Omegam]<>"*"]&];
gammaSelect=If[Association[]!=KeySelect[paropts,StringMatchQ[ToString[#],"*"<>SymbolName[gamma]<>"*"]&],KeySelect[paropts,StringMatchQ[ToString[#],"*"<>SymbolName[gamma]<>"*"]&],KeySelect[paropts,StringMatchQ[ToString[#],"*"<>SymbolName[hubble]<>"*"]&]];
OmegamPos=First@First@Position[paropts,First@Keys@OmegamSelect];
gammaPos=First@First@Position[paropts,First@Keys@gammaSelect];
newlist=paropts;
newlist[[OmegamPos]]=w0->w0CPL[paropts];
newlist[[gammaPos]]=wa->waCPL[paropts];
newlist];
$CPLParamFiducial=$CPLParamList[[All,2]];
$CPLParamLabel=SymbolName/@$CPLParamList[[All,1]];
];

Options[computeCPLJacobianMatrix]=$paramoptions;
computeCPLJacobianMatrix[opts:OptionsPattern[]]:=Module[{paropts,newTOoldJacobianMatrix,oldTOnewJacobianMatrix},
paropts=complementParamValues[{opts},computeCPLJacobianMatrix,returnList->"Full", filterCosmoPars->True];

newTOoldJacobianMatrix=Module[{w0Pos,waPos,OmegamSelect,gammaSelect},
OmegamSelect=KeySelect[paropts,StringMatchQ[ToString[#],"*"<>SymbolName[Omegam]<>"*"]&];
gammaSelect=If[Association[]!=KeySelect[paropts,StringMatchQ[ToString[#],"*"<>SymbolName[gamma]<>"*"]&],KeySelect[paropts,StringMatchQ[ToString[#],"*"<>SymbolName[gamma]<>"*"]&],KeySelect[paropts,StringMatchQ[ToString[#],"*"<>SymbolName[hubble]<>"*"]&]];
w0Pos=First@First@Position[paropts,First@Keys@OmegamSelect];
waPos=First@First@Position[paropts,First@Keys@gammaSelect];

Table[KroneckerDelta[i,j](1-KroneckerDelta[i,w0Pos]-KroneckerDelta[i,waPos])
+w0dpi[j]*KroneckerDelta[i,w0Pos]+wadpi[j]*KroneckerDelta[i,waPos]
,{i,1,Length@paropts},{j,1,Length@paropts}]
];

oldTOnewJacobianMatrix=Inverse[newTOoldJacobianMatrix]
];


(* ::Subsection::Closed:: *)
(*For hubble*)


Options[setHubbleParamList]=$paramoptions;
setHubbleParamList[opts:OptionsPattern[]]:=Module[{paropts},
paropts=complementParamValues[{opts},setHubbleParamList,returnList->"Full", filterCosmoPars->True];

$HubbleParamList=Module[{gammaSelect,gammaPos,varphiFPos,newlist},
gammaSelect=KeySelect[paropts,StringMatchQ[ToString[#],"*"<>SymbolName[gamma]<>"*"]&];
gammaPos=First@First@Position[paropts,First@Keys@gammaSelect];
newlist=paropts;
newlist[[gammaPos]]=hubble->hubbleToday[paropts];
newlist];
$HubbleParamFiducial=$HubbleParamList[[All,2]];
$HubbleParamLabel=SymbolName/@$HubbleParamList[[All,1]];
];

Options[computeHubbleJacobianMatrix]=$paramoptions;
computeHubbleJacobianMatrix[opts:OptionsPattern[]]:=Module[{paropts,newTOoldJacobianMatrix,oldTOnewJacobianMatrix},
paropts=complementParamValues[{opts},computeHubbleJacobianMatrix,returnList->"Full", filterCosmoPars->True];

newTOoldJacobianMatrix=Module[{hubblePos,gammaSelect},
gammaSelect=KeySelect[paropts,StringMatchQ[ToString[#],"*"<>SymbolName[gamma]<>"*"]&];
hubblePos=First@First@Position[paropts,First@Keys@gammaSelect];

Table[KroneckerDelta[i,j](1-KroneckerDelta[i,hubblePos])
+hubbledpi[j]*KroneckerDelta[i,hubblePos],{i,1,Length@paropts},{j,1,Length@paropts}]
];

oldTOnewJacobianMatrix=Inverse[newTOoldJacobianMatrix]
];


(* ::Section::Closed:: *)
(*Miscellaneous*)


(*ListPlot@Table[{xx,w0dpi[6,setEpsilon\[Rule]xx]},{xx,0.0001,0.005,0.0001}]*)


(*ListPlot@Table[{xx,wadpi[1,setEpsilon\[Rule]xx]},{xx,0.0001,0.005,0.0001}]*)


Options@paramDerivative





(* ::Section:: *)
(*End of package*)


End[]
EndPackage[]
