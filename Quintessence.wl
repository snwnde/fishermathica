(* ::Package:: *)

BeginPackage["Quintessence`"]


(* ::Section:: *)
(*Contextualization*)


V::usage="Quintessence potential. V[phi_,opts:OptionsPattern[]]"
alphaGen::usage="Gives the alpha-attractor model parameter alpha"
gammaGen::usage="Gives the alpha-attractor model parameter gamma"
M2Gen::usage="Gives the alpha-attractor model parameter M^2"
LambdaGen::usage="Gives the alpha-attractor model parameter Lambda"
varphiFGen::usage="Gives the initial value of varphi when solving the system"
$h0computed::usage="Dimensionless h0 computed recursively"
saveH0computed::usage="Save the dimensionless h0 computed recursively avoiding repetitive calculation"
curlyPhi::usage="The quintessential field varphi"


rGen::usage="Gives the tensor-to-scalar ratio r"


weff::usage="Effective EOS parameter"
wDE::usage="Dark energy EOS parameter"


alpha::usage="Set the quintessential protected variable alpha"
gamma::usage="Set the quintessential protected variable gamma"
Lambda::usage="Set the quintessential protected variable Lambda"
varphiF::usage="Set the quintessential protected variable varphiF"


$nameOfModel::usage="Model name"
$unitsScalingFactor::usage="Units change factor"
setQuintessenceFixedValues::usage="Ser quintessence fixed values"
$gammaFixed::usage="Fixed value of gamma"
$alphaFixed::usage="Fixed value of alpha"
$varphiFFixed::usage="Fixed value of varphiF"


$scalarParametersProtectedList::usage="gives the list of available quintessential parameters, whose names are protected and are available for their use in a Fisher forecast"


Needs["CosmologyFunctions`"];
Needs["UsefulTools`"];
Needs["SurveySpecifications`"];
Needs["FunctionApproximations`"];
EndPackage[]


Protect[alpha,gamma,Lambda,varphiF]


BeginPackage["Quintessence`", {"CosmologyFunctions`","UsefulTools`","SurveySpecifications`","FunctionApproximations`"}]
Begin["`Private`"]
$unitsScalingFactor=3.08567758*10^22/10^3/(5.39116*10^-44)/Sqrt[8Pi];(*Mpc*sec/km/tp/Sqrt[8Pi]*)


(* ::Section::Closed:: *)
(*Some parameter transformations*)


$scalarParametersProtectedList={alpha,gamma,Lambda}


gammaGen::noparam="Symbol gamma or a function of it, is not part of the $paramoptions parameters"
LambdaGen::noparam="Symbol Lambda or a function of it, is not part of the $paramoptions parameters"
alphaGen::noparam="Symbol alpha or a function of it, is not part of the $paramoptions parameters"
varphiFGen::usage="Gives the initial value of varphi when solving the system"
M2Gen::noparamNs="There is no parameter ns"
rGen::noparamNs="There is no parameter ns"


Options[alphaGen]=$paramoptions~Join~{transformedParameter->True};
alphaGen[opts:OptionsPattern[]]:=Block[{paropts, alphaValue, rulepar},
paropts=complementParamValues[{opts},alphaGen,returnList->"Full", filterCosmoPars->True];
rulepar=FilterRules[paropts,_?(StringMatchQ[ToString[#],"*"<>SymbolName[alpha]<>"*"]&)];
alphaValue=Which[
  SameQ[rulepar, {}],
    (*Message[alphaGen::noparam];
    Return[$Failed],*)
    $alphaFixed,
  OptionValue[transformedParameter]==True && UnsameQ[First@rulepar[[1]], alpha],
    setORreadTransformedParameters[alpha, paropts, convertTransformedParameter->True],
  (OptionValue[transformedParameter]==False && FilterRules[paropts,alpha]=!={}),
    alpha/.FilterRules[paropts,alpha],
  (OptionValue[transformedParameter]==True && SameQ[First@rulepar[[1]], alpha]),
    SetOptions[alphaGen, transformedParameter->False];
    alpha/.FilterRules[paropts,alpha]
  ]
]

Options[gammaGen]=$paramoptions~Join~{transformedParameter->True};
gammaGen[opts:OptionsPattern[]]:=Block[{paropts, gammaValue, rulepar},
paropts=complementParamValues[{opts},gammaGen,returnList->"Full", filterCosmoPars->True];
rulepar=FilterRules[paropts,_?(StringMatchQ[ToString[#],"*"<>SymbolName[gamma]<>"*"]&)];
gammaValue=Which[
  SameQ[rulepar, {}],
    (*Message[gammaGen::noparam];
    Return[$Failed],*)
    $gammaFixed,
  OptionValue[transformedParameter]==True && UnsameQ[First@rulepar[[1]], gamma],
    setORreadTransformedParameters[gamma, paropts, convertTransformedParameter->True],
  (OptionValue[transformedParameter]==False && FilterRules[paropts,gamma]=!={}),
    gamma/.FilterRules[paropts,gamma],
  (OptionValue[transformedParameter]==True && SameQ[First@rulepar[[1]], gamma]),
    SetOptions[gammaGen, transformedParameter->False];
    gamma/.FilterRules[paropts,gamma]
  ]
]

Options[M2Gen]=$paramoptions;
M2Gen[opts:OptionsPattern[]]:=Block[{paropts, M2Value,AsValue,NValue,nsValue,alphaValue},
paropts=complementParamValues[{opts},M2Gen,returnList->"Full", filterCosmoPars->True];
AsValue=scalarAmplitude[paropts];
nsValue=If[FilterRules[paropts,ns]!={},OptionValue[ns],Message[M2Gen::noparamNs]];
NValue=2/(1-nsValue);
alphaValue=alphaGen[paropts];
M2Value=144Pi^2*alphaValue*NValue/(2NValue-3alphaValue)^3 * AsValue
]

Options[LambdaGen]=$paramoptions~Join~{transformedParameter->True};
LambdaGen[opts:OptionsPattern[]]:=Block[{paropts, LambdaValue, rulepar},
paropts=complementParamValues[{opts},LambdaGen,returnList->"Full", filterCosmoPars->True];
rulepar=FilterRules[paropts,_?(StringMatchQ[ToString[#],"*"<>SymbolName[Lambda]<>"*"]&)];
LambdaValue=Which[
  SameQ[rulepar, {}],
    Message[LambdaGen::noparam];
    Return[$Failed],
  OptionValue[transformedParameter]==True && UnsameQ[First@rulepar[[1]], Lambda],
    setORreadTransformedParameters[Lambda, paropts, convertTransformedParameter->True],
  (OptionValue[transformedParameter]==False && FilterRules[paropts,Lambda]=!={}),
    Lambda/.FilterRules[paropts,Lambda],
  (OptionValue[transformedParameter]==True && SameQ[First@rulepar[[1]], Lambda]),
    SetOptions[LambdaGen, transformedParameter->False];
    Lambda/.FilterRules[paropts,Lambda]
  ]
]

(*Options[varphiFGen]=$paramoptions;
varphiFGen[opts:OptionsPattern[]]:=-35;*)

Options[varphiFGen]=$paramoptions~Join~{transformedParameter->True};
varphiFGen[opts:OptionsPattern[]]:=Block[{paropts,varphiFValue,rulepar},
paropts=complementParamValues[{opts},varphiFGen,returnList->"Full", filterCosmoPars->True];
rulepar=FilterRules[paropts,_?(StringMatchQ[ToString[#],"*"<>SymbolName[varphiF]<>"*"]&)];
varphiFValue=Which[
  SameQ[rulepar, {}],
    (*Message[varphiFGen::noparam];
    Return[$Failed],*)
    $varphiFFixed,
  OptionValue[transformedParameter]==True && UnsameQ[First@rulepar[[1]], varphiF],
    setORreadTransformedParameters[varphiF, paropts, convertTransformedParameter->True],
  (OptionValue[transformedParameter]==False && FilterRules[paropts,varphiF]=!={}),
    varphiF/.FilterRules[paropts,varphiF],
  (OptionValue[transformedParameter]==True && SameQ[First@rulepar[[1]], varphiF]),
    SetOptions[varphiFGen, transformedParameter->False];
    varphiF/.FilterRules[paropts,varphiF]
  ]
]


Options[rGen]=$paramoptions;
rGen[opts:OptionsPattern[]]:=Block[{paropts, rValue,AsValue,NValue,nsValue,alphaValue},
paropts=complementParamValues[{opts},rGen,returnList->"Full", filterCosmoPars->True];
AsValue=scalarAmplitude[paropts];
nsValue=If[FilterRules[paropts,ns]!={},OptionValue[ns],Message[rGen::noparamNs]];
NValue=2/(1-nsValue);
alphaValue=alphaGen[paropts];
rValue=12*alphaValue/NValue^2
]


(* ::Section::Closed:: *)
(*Alpha - attractor model functions*)


Options[V]=$paramoptions;
V[phi_,opts:OptionsPattern[]]:=
Module[{OmegaM0,hubref,OmegaR0,paropts,flatness,alpha,gamma,M2,Lambda},
paropts=complementParamValues[{opts},V,returnList->"Complement", filterCosmoPars->True];
OmegaM0=OmegaM0Today[paropts];
OmegaR0=OmegaR0Today[paropts];
hubref=$hubblereference*100/$unitsScalingFactor;
alpha=alphaGen[paropts];
gamma=gammaGen[paropts];
Which[
StringMatchQ[$nameOfModel,"Exp* I"],
M2=M2Gen[paropts];
M2*Exp[gamma*(Tanh[phi/Sqrt[6*alpha]]-1)],
StringMatchQ[$nameOfModel,"lCDM"]||StringMatchQ[$nameOfModel,"lcdm"],
(1-OmegaM0-OmegaR0)*3hubref^2,
StringMatchQ[$nameOfModel,"Exp* II"],
M2=M2Gen[paropts];
M2*Exp[-2gamma]*( Exp[gamma (Tanh[phi/Sqrt[6 alpha]]+1)]-1),
StringMatchQ[$nameOfModel,"linear*"],
Lambda=LambdaGen[paropts];
gamma*Sqrt[6 alpha](Tanh[phi/Sqrt[6 alpha]]+1)+Lambda
]
];


Options[reducedV]=$paramoptions;
reducedV[phi_,opts:OptionsPattern[]]:=Block[{paropts},
paropts=complementParamValues[{opts},reducedV,returnList->"Complement", filterCosmoPars->True];
V[phi,opts]/($h0computed^2)];

dReducedV[phi_,opts:OptionsPattern[]]:=\!\(
\*SubscriptBox[\(\[PartialD]\), \(phi\)]\(reducedV[phi, opts]\)\);
dCurlyPhi[xx_,opts:OptionsPattern[]]:=D[curlyPhi[x,opts],x]/.x-> xx;


Options[ee]=$paramoptions;
ee[eFolds_,phi_,opts:OptionsPattern[]]:=Block[{OmegaM0,OmegaR0,paropts},
paropts=complementParamValues[{opts},ee,returnList->"Complement", filterCosmoPars->True];
OmegaM0=OmegaM0Today[paropts];
OmegaR0=OmegaR0Today[paropts];
(reducedV[phi[eFolds],opts]/3+OmegaM0*Exp[-3eFolds]+OmegaR0*Exp[-4eFolds])/(1-1/6 phi'[eFolds]^2)
];

Options[epsilon]=$paramoptions;
epsilon[eFolds_,phi_,opts:OptionsPattern[]]:=Block[{OmegaM0,OmegaR0,paropts},
paropts=complementParamValues[{opts},epsilon,returnList->"Complement", filterCosmoPars->True];
OmegaM0=OmegaM0Today[paropts];
OmegaR0=OmegaR0Today[paropts];
1/2*(phi'[eFolds]^2+(3OmegaM0*Exp[-3eFolds]+4OmegaR0*Exp[-4eFolds])/ee[eFolds,phi,opts])
];


(* ::Section:: *)
(*Equation solver*)


Options[findH0computed]=$paramoptions~Join~{doIteration->True};
findH0computed[opts:OptionsPattern[]]:=Block[{OmegaM0,paropts,OmegaR0},
paropts=complementParamValues[{opts},findH0computed,returnList->"Complement", filterCosmoPars->True];
OmegaM0=OmegaM0Today[paropts];
OmegaR0=OmegaR0Today[paropts];

If[OptionValue[doIteration],
(iteration[paropts];
findH0computed[paropts,doIteration->False]),

If[NumericQ[$h0computed],
Sqrt[V[curlyPhi[0,paropts,doIteration->False],paropts]/(3*(1-1/6*dCurlyPhi[0,paropts,doIteration->False]^2-OmegaM0-OmegaR0))],
$hubblereference*100/$unitsScalingFactor]
]
];


Options[iteration]=$paramoptions;
iteration[opts:OptionsPattern[]]:=Block[{OmegaM0,OmegaR0,paropts,tole=10^-10},
paropts=complementParamValues[{opts},iteration,returnList->"Complement", filterCosmoPars->True];
OmegaM0=OmegaM0Today[paropts];
OmegaR0=OmegaR0Today[paropts];

Do[
$h0computed = findH0computed[paropts,doIteration->False];
If[Abs[dimensionlessHubbleWrtN[0,paropts,doIteration->False]-1]<tole,Return[Null]],
Infinity]
];


Options[saveH0computed]=$paramoptions;
saveH0computed[opts:OptionsPattern[]]:=saveH0computed[opts]=findH0computed[opts,doIteration->True];


Options[dimensionlessHubbleWrtN]=$paramoptions~Join~{doIteration->True};
dimensionlessHubbleWrtN[eFolds_,opts:OptionsPattern[]]:=Block[{OmegaM0,OmegaR0,paropts},
paropts=complementParamValues[{opts},dimensionlessHubbleWrtN,returnList->"Complement", filterCosmoPars->True];
OmegaM0=OmegaM0Today[paropts];
OmegaR0=OmegaR0Today[paropts];
If[OptionValue[doIteration],
($h0computed=saveH0computed[paropts];
dimensionlessHubbleWrtN[eFolds,paropts,doIteration-> False]),

Sqrt[(reducedV[curlyPhi[eFolds,paropts,doIteration->False],paropts]/3+OmegaM0*Exp[-3eFolds]+OmegaR0*Exp[-4eFolds])/(1-1/6 dCurlyPhi[eFolds,paropts,doIteration->False]^2)]
]
];


Options[curlyPhi]=$paramoptions~Join~{doIteration->True};

curlyPhi[eFolds_,opts:OptionsPattern[]]:=
Module[{paropts,varphiF},
paropts=complementParamValues[{opts},curlyPhi,returnList->"Complement", filterCosmoPars->True];
varphiF=varphiFGen[paropts];

If[OptionValue[doIteration],
($h0computed=saveH0computed[paropts];
curlyPhi[eFolds,paropts,doIteration-> False]),

(Module[{eqnPhi,iniPhi,iniPhiPrime,phiii,phii,xx},
eqnPhi= (phii'')[xx]+(3-epsilon[xx,phii,paropts]) Derivative[1][phii][xx]+dReducedV[phii[xx],paropts]/ee[xx,phii,paropts]==0;
iniPhi=phii[-15]==varphiF;
iniPhiPrime=Derivative[1][phii][-15]==0;
phiii=NDSolveValue[{eqnPhi,iniPhi, iniPhiPrime},phii,{xx,-15,10}];

phiii[eFolds]])
]
];


Options[slowrollWrtN]=$paramoptions~Join~{doIteration->True};
slowrollWrtN[eFolds_,opts:OptionsPattern[]]:=Block[{OmegaM0,OmegaR0,paropts},
paropts=complementParamValues[{opts},slowrollWrtN,returnList->"Complement", filterCosmoPars->True];
OmegaM0=OmegaM0Today[paropts];
OmegaR0=OmegaR0Today[paropts];
If[OptionValue[doIteration],
($h0computed=saveH0computed[paropts];
slowrollWrtN[eFolds,paropts,doIteration-> False]),
1/2*(dCurlyPhi[eFolds,paropts,doIteration->False]^2+(3OmegaM0*Exp[-3eFolds]+4OmegaR0*Exp[-4eFolds])/dimensionlessHubbleWrtN[eFolds,paropts,doIteration->False]^2)
]
];


(* ::Section:: *)
(*Overwriting hubbleToday*)


hubbleToday[opts:OptionsPattern[]]:=Block[{paropts, hubblevalue},
paropts=complementParamValues[{opts},hubbleToday,returnList->"Full", filterCosmoPars->True];
hubblevalue=If[FilterRules[paropts,hubble]!={},OptionValue[hubble],saveH0computed[paropts]*$unitsScalingFactor/100]
];


(* ::Section:: *)
(*Set fixed values*)


Options[setQuintessenceFixedValues]={gammaFixed->$gammaFixed,alphaFixed->$alphaFixed,varphiFFixed->$varphiFFixed};
setQuintessenceFixedValues[opts:OptionsPattern[]]:=Module[{},
$gammaFixed=OptionValue[gammaFixed];
$alphaFixed=OptionValue[alphaFixed];
$varphiFFixed=OptionValue[varphiFFixed];
];


(* ::Section::Closed:: *)
(*Equation of state*)


Options[weff]=$paramoptions;
weff[eFolds_,opts:OptionsPattern[]]:=Module[{paropts,phi},
paropts=complementParamValues[{opts},weff,returnList->"Full", filterCosmoPars->True];
phi=curlyPhi[#,paropts]&;
-1+2slowrollWrtN[eFolds,paropts]/3
];


Options[wDE]=$paramoptions;
wDE[eFolds_,opts:OptionsPattern[]]:=Module[{paropts,redshift,hub,potential,phi,dphi,shorthand},
paropts=complementParamValues[{opts},wDE,returnList->"Full", filterCosmoPars->True];
redshift=zOfNf[eFolds];
phi=curlyPhi[#,paropts]&;
dphi=dCurlyPhi[#,paropts]&;
hub=dimensionlessHubbleWrtN[eFolds,paropts]*saveH0computed[paropts];
potential=V[phi[eFolds],paropts];
shorthand=1/2*dphi[eFolds]^2*hub^2;
(shorthand-potential)/(shorthand+potential)
];


(* ::Section::Closed:: *)
(*End of package*)


End[]

EndPackage[]
