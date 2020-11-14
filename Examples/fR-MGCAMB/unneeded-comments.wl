(* ::Package:: *)

(*fiducialPkzNL = externalMatrixListImport["Pnlzk","fiducial",0.];*)


(*fiducialPkzL = externalMatrixListImport["Plzk","fiducial",0.];*)


(*(fiducialPkzL//Dimensions)*)


(*ListLogLogPlot[{Transpose[{$krangePowerSpectrumFiducial, fiducialPkzL[[1,All]]}], 
Transpose[{$krangePowerSpectrumFiducial, fiducialPkzNL[[1,All]]}]}]*)


(*obsnam="Pnlzk";
zval=0.;
paramval="Om";
epsilonvalue=$epsilonstep;
fiducialobs=externalObservableFunc[obsnam,"fiducial", 0.];
plusobs=externalObservableFunc[obsnam,paramval, epsilonvalue];
minusobs=externalObservableFunc[obsnam,paramval, -epsilonvalue];
LogLinearPlot[{(plusobs[zval,kk]-fiducialobs[zval,kk])/fiducialobs[zval,kk], 
(minusobs[zval,kk]-fiducialobs[zval,kk])/fiducialobs[zval,kk]} , {kk, 10^-4, 5}];
Manipulate[
kval=kki;
ppTable=observableParamDependency[obsnam,paramval,zval,kval];
{der,rmse,lm,cutpartab}=steMderiv[ppTable];
twoptder=twoPointDerivative[obsnam,paramval,zval,kval, epsilonvalue, True];
{der4,rmse4,lm4,cutpartab4}=steMderiv[ppTable,4];
Show[Plot[{1+x*der,1+x*der4,1+x*twoptder},{x,-0.1,0.1}, PlotStyle->{Red,Orange,Purple}, PlotRange->{0.99,1.01} ],
Plot[{lm4[x]},{x,-0.1,0.1}, PlotStyle->{Orange},PlotRange->Full],
ListPlot[{ppTable,cutpartab},PlotLabel->"k="<>ToString@kval,PlotStyle->{Red,Orange,Purple}],
Frame->True]
,{{kki,0.8}, 0.001,1.,0.001}]*)


(*twoptder*fiducialobs[zval,0.8]


der*fiducialobs[zval,0.8]


der4*fiducialobs[zval,0.8]
*)


(*
obsnam2="Hz";
zval2=0.3;
paramval2="Om";
epsilonvalue2=0.05;
fiducialobs2=externalObservableFunc[obsnam2,"fiducial", 0.];
plusobs2=externalObservableFunc[obsnam2,paramval2, epsilonvalue2];
minusobs2=externalObservableFunc[obsnam2,paramval2, -epsilonvalue2];
Plot[{(plusobs2[zz]-fiducialobs2[zz])/fiducialobs2[zz], 
(minusobs2[zz]-fiducialobs2[zz])/fiducialobs2[zz]} , {zz, 0.001, 2.2}]


Manipulate[
zval=zzi;
ppTable2=observableParamDependency[obsnam2,paramval2,zval2,kval2];
{der2,rmse2,lm2,cutpartab2}=steMderiv[ppTable2];
twoptder2=twoPointDerivative[obsnam2,paramval2,zval2,kval2, epsilonvalue2, True];
Show[ListPlot[{ppTable2,cutpartab2},PlotRange->Full,PlotLabel->"k="<>ToString@kval],Plot[lm2[x],{x,-0.1,0.1},PlotRange->Full],
Plot[1+x*twoptder2,{x,-0.1,0.1}, PlotStyle->Red, PlotRange->Full ],Frame->True]
,{zzi, 0.001,2.,0.1}]*)



(*LogLogPlot[{externalObservableFunc["Plzk","logfR0", 0.1][1.0,kk]/externalObservableFunc["Plzk","logfR0", 0.1][0.,kk], 
externalObservableFunc["Plzk","logfR0", 0.1][1.0,kk]/externalObservableFunc["Plzk","fiducial", 0.][1.0,kk],
{externalObservableFunc["Plzk","logfR0", 0.05][0.,kk]/externalObservableFunc["Plzk","logfR0", 0.1][0.5,kk]}} , {kk, 10^-4, 5}, PlotLegends->Automatic]*)
(*myInterrupt[]*)


(*parnam="h";
fiduval=$parameterFiducials[parnam];
paptab=observableParamDependency["Pnlzk",parnam,0.,0.9];
range={(1+$epslistPM[[1]])*$parameterFiducials[parnam],(1+$epslistPM[[-1]])*$parameterFiducials[parnam]};
{relderiv,rmse,lml,cutpartab}=steMderiv[paptab,fiduval,1];
{relderiv,rmse,lml4,cutpartab}=steMderiv[paptab,fiduval,4];
Show[ListPlot[paptab],Plot[{lml[x],lml4[x]},{x,range[[1]],range[[2]]}]]*)


(*observableDerivativeZK["Pnlzk","h",0.,0.9, derivativeType->"2pt"]*)


(*observableDerivativeZK["Pnlzk","h",0.,0.9, derivativeType->"SteM1"]*)


(*observableDerivativeZK["Pnlzk","h",0.,0.9, derivativeType->"SteM4"]*)


(*
kkkk=$krangePowerSpectrumFiducial[[402]]


zzzz=$zrangePowerSpectrum[[50]]


observableDerivativeZK["Pnlzk","Om",zzzz,kkkk, derivativeType->"2pt"]


observableDerivativeZK["Pnlzk","Om",zzzz,kkkk, derivativeType->"Stem1"]


observableDerivativeZK["Pnlzk","Om",zzzz,kkkk, derivativeType->"Stem4"]


observableDerivativeZK["fgzk","Om",zzzz,kkkk, derivativeType->"2pt"]


observableDerivativeZK["fgzk","Om",zzzz,kkkk, derivativeType->"Stem1"]


observableDerivativeZK["fgzk","Om",zzzz,kkkk, derivativeType->"Stem4"]


observableDerivativeZK["Hz","Om",zzzz,kkkk, derivativeType->"2pt"]


observableDerivativeZK["Hz","Om",zzzz,kkkk, derivativeType->"Stem1"]


observableDerivativeZK["Hz","Om",zzzz,kkkk, derivativeType->"Stem4"]


$paramlabels[[1;;1]]


$krangePowerSpectrumFiducial;


$observablesNames


$observablesNames[[1;;2]]*)


(*zeto=1.0;
Table[
paramstr=pp;
obsstr="Pnlzk";
LogLinearPlot[{Power[externalObservableDerivative2ptFunc[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak],2]
,Power[externalObservableDerivativeStem1Func[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak],2],
Power[externalObservableDerivativeStem4Func[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak],2]}, {kak, 1*^-4,10}, PlotLegends->Automatic,PlotLabel->pp],
{pp,$paramlabels}]*)


(*pk2ptMatPlot=GraphicsGrid[Table[
paramstr=pp;
obsstr="Plzk";
LogLinearPlot[{Identity[externalObservableDerivative2ptFunc[obsstr,pi][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]
*Identity[externalObservableDerivative2ptFunc[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]}, {kak, 1*^-4,10}, PlotLegends->Automatic,PlotLabel->pi<>" x "<>pp,
ImageSize->200],
{pp,$paramlabels},{pi,$paramlabels}
], Frame->All]
Export["pkL_derivatives2pt_MatPlot.pdf",pk2ptMatPlot,"PDF"]
pk2ptMatPlot=GraphicsGrid[Table[
paramstr=pp;
obsstr="Pnlzk";
LogLinearPlot[{Identity[externalObservableDerivative2ptFunc[obsstr,pi][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]
*Identity[externalObservableDerivative2ptFunc[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]}, {kak, 1*^-4,10}, PlotLegends->Automatic,PlotLabel->pi<>" x "<>pp,
ImageSize->200],
{pp,$paramlabels},{pi,$paramlabels}
], Frame->All]
Export["pkNL_derivatives2pt_MatPlot.pdf",pk2ptMatPlot,"PDF"]*)


(*pkSteM1MatPlot=GraphicsGrid[Table[
paramstr=pp;
obsstr="Plzk";
LogLinearPlot[{Identity[externalObservableDerivativeStem1Func[obsstr,pi][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]
*Identity[externalObservableDerivativeStem1Func[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]}, {kak, 1*^-4,10}, PlotLegends->Automatic,PlotLabel->pi<>" x "<>pp,
ImageSize->200],
{pp,$paramlabels},{pi,$paramlabels}
], Frame->All]
Export["pkL_derivativesSteM1_MatPlot.pdf",pkSteM1MatPlot,"PDF"]
pkSteM1MatPlot=GraphicsGrid[Table[
paramstr=pp;
obsstr="Pnlzk";
LogLinearPlot[{Identity[externalObservableDerivativeStem1Func[obsstr,pi][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]
*Identity[externalObservableDerivativeStem1Func[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]}, {kak, 1*^-4,10}, PlotLegends->Automatic,PlotLabel->pi<>" x "<>pp,
ImageSize->200],
{pp,$paramlabels},{pi,$paramlabels}
], Frame->All]
Export["pkNL_derivativesSteM1_MatPlot.pdf",pkSteM1MatPlot,"PDF"]*)


(*pkSteM4MatPlot=GraphicsGrid[Table[
paramstr=pp;
obsstr="Plzk";
LogLinearPlot[{Identity[externalObservableDerivativeStem4Func[obsstr,pi][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]
*Identity[externalObservableDerivativeStem4Func[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]}, {kak, 1*^-4,10}, PlotLegends->Automatic,PlotLabel->pi<>" x "<>pp,
ImageSize->200],
{pp,$paramlabels},{pi,$paramlabels}
], Frame->All]
Export["pkL_derivativesSteM4_MatPlot.pdf",pkSteM4MatPlot,"PDF"]
pkSteM4MatPlot=GraphicsGrid[Table[
paramstr=pp;
obsstr="Pnlzk";
LogLinearPlot[{Identity[externalObservableDerivativeStem4Func[obsstr,pi][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]
*Identity[externalObservableDerivativeStem4Func[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]}, {kak, 1*^-4,10}, PlotLegends->Automatic,PlotLabel->pi<>" x "<>pp,
ImageSize->200],
{pp,$paramlabels},{pi,$paramlabels}
], Frame->All]
Export["pkNL_derivativesSteM4_MatPlot.pdf",pkSteM4MatPlot,"PDF"]*)


Hubble[1.0, lcdmBool->True, externalFile->False]


Hubble[1.0]


Hubble[0., physicalUnits->True]


Hubble[0., lcdmBool->True, externalFile->False, physicalUnits->True]


Hubble[0., hubble->0.68, physicalUnits->True]


ztesta=0.


(Hubble[ztesta, hubble->0.67*(1+0.01), lcdmBool->True, externalFile->False, physicalUnits->False]-Hubble[ztesta, hubble->0.67*(1-0.01), lcdmBool->True, externalFile->False, physicalUnits->False])/(2*0.01*0.67)


(Hubble[ztesta, hubble->0.67*(1+0.01), physicalUnits->False]-Hubble[ztesta, hubble->0.67*(1-0.01), physicalUnits->False])/(2*0.01*0.67)


(*paropts=complementParamValues[{Omegam->0.33},externalPowerSpectrumFunction,returnList->"Complement",filterCosmoPars->True]


First@First@paropts


First@First@$paramoptions


First@complementParamValues[{paropts},externalFunctionTaylor,returnList->"Fiducials",filterCosmoPars->True]


$paramoptions


Hubble[0., Omegam->0.33]


setExternalCosmoInterpolatingFunctions[externalPowerSpectrumInput->extPowerSpecFunc, externalGrowthRateInput ->extGrowthRateFunc , 
externalSigma8ofZInput->extSigma8FuncOfZ, externalPowerSpectrumNoWiggleInput->extPowerNWSpecFunc, externalGrowthInput ->extGrowthFunc]

*)


(*ClearAll[externalFunctionTaylorTest];

Options[externalFunctionTaylorTest]=$paramoptions~Join~{functionInterpolatedInLogLog->False,derivativeInterpolatedInLogLin->False,
  externalFunction->True,externalDerivativeFunction->True, interpolatedArguments->"k", kdependence->True, externalSpectraType -> "pk"};

externalFunctionTaylorTest[zkarg__,opts:OptionsPattern[]]:=Module[{zz, kk, kdep, par,eps,paropts,parfidu,filename,
  fidu,return,intpDerlinlog,intpFunloglog,funcy=(#)&,funcx=(#)&,funcyprime=(#)&,funcxprime=(#)&,fiduval,parval,Delta,tolerance=10^(-6)
(*hard coded tolerance to distinguish between fiducial and value away from fiducial*), intpargs,
  extFunction,extFunctionDerivative,
  externalFuncTemplate, externalFuncDerivTemplate, parname, parlab, parsrules},
  
  parsrules=Thread[Rule[$paramnames,$paramlabels]];
  extFunction=OptionValue[externalFunction];
  extFunctionDerivative=OptionValue[externalDerivativeFunction];
  If[BooleanQ[extFunction] || BooleanQ[extFunctionDerivative], Message[externalFunctionTaylor::extfuncs]; Abort[]];

  kdep=OptionValue[kdependence];

  {zz,kk}=twoArgumentsCheck[zkarg,kdep];
  kk=ReleaseHold[kk];

  paropts=complementParamValues[{opts},externalFunctionTaylor,returnList->"Complement",filterCosmoPars->True];
  (*fidu=getCosmoParameterFileName[externalFunctionTaylor,paropts, getFiducialFileName->True];*)
  intpFunloglog=OptionValue[functionInterpolatedInLogLog];
  intpDerlinlog=OptionValue[derivativeInterpolatedInLogLin];
  If[intpFunloglog==True,funcy=(10^#)&;funcx=Log10[#]&;];
  If[intpDerlinlog==True,funcyprime=(#)&;funcxprime=Log10[#]&;];

  (*If[StringMatchQ[fidu,"Fiducial*",IgnoreCase->True]==False,Message[getCosmoParameterFileName::fiduwarn,fidu]];*)
  debugPrint["paropts: "<>ToString@paropts,1];

  If[Length@paropts>=1,
    paropts=First@paropts;
    parfidu=First@complementParamValues[{paropts},externalFunctionTaylor,returnList->"Fiducials",filterCosmoPars->True];
    fiduval=Last@parfidu;
    parval=Last@paropts;
    (*filename=getCosmoParameterFileName[externalFunctionTaylor,paropts,inputNumericalDerivatives->True];*)
    parname=First@paropts;
    Delta=Chop[(parval-fiduval),tolerance];
    , (*Else*)
    parname=First@First@$paramoptions; (*get a parameter name that is not "fiducial", to avoid finding a wrong name in derivative file"*)
    Delta=0;
  ];
    
    parlab=parname/.parsrules;
   
  intpargs=OptionValue[interpolatedArguments];

  Which[
    intpargs=="zk",
    externalFuncTemplate[zaa_,kaa_]:=extFunction[zaa,funcx@kaa];
    externalFuncDerivTemplate[zaa_,kaa_]:=extFunctionDerivative[parlab][zaa,funcxprime@kaa];
    ,
    intpargs=="k",
    externalFuncTemplate[zaa_,kaa_]:=extFunction[[zValTozIndex[zaa]]][funcx@kaa];
    externalFuncDerivTemplate[zaa_,kaa_]:=extFunctionDerivative[parlab][[zValTozIndex[zaa]]][funcxprime@kaa];
    ,
    intpargs=="z",
    externalFuncTemplate[zaa_]:=extFunction[zaa];
    externalFuncDerivTemplate[zaa_]:=extFunctionDerivative[parlab][zaa];
  ];

  debugPrint["fidu val: "<>ToString@fiduval];
  debugPrint["param val: "<>ToString@parval];
  debugPrint["param label: "<>ToString@parlab];
  debugPrint["Delta: "<>ToString@Delta];
  debugPrint["Delta: "<>ToString@Delta];

  return=funcy@(externalFuncTemplate[zz,kk])+Delta*funcyprime@(externalFuncDerivTemplate[zz,kk]);

  Return[return]
]*)


Options[powerSpectrum]


$biasInterpFunc[0.1]


$kminRef


Options[NIntegrate]


NIntegrate[powerSpectrum[0.5,kk],{kk,$kminRef,0.2}]//AbsoluteTiming


NIntegrate[powerSpectrum[0.5,kk],{kk,$kminRef,0.2}]//AbsoluteTiming


NIntegrate[observedPowerSpectrum[0.5,kk, 0.9],{kk,$kminRef,0.2}]//AbsoluteTiming


logkrange=logarithmicDivisions[{$kminRef, 0.2}, 250];
krange=Table[kk,{kk, $kminRef, 0.2, 0.001}];
Dimensions[krange]


(pp=observedPowerSpectrum[0.5,#, 0.88]&/@krange);//AbsoluteTiming


(pp2=observedPowerSpectrum[0.5,#, 0.88]&/@logkrange);//AbsoluteTiming





ttab1=Transpose[{krange, pp}];
ttab2=Transpose[{logkrange, pp2}];


integratetrapz[ttab1]//AbsoluteTiming


integratetrapz[ttab2]//AbsoluteTiming


$kmaxRef


kmaxChoice[1.0]


powerSpectrum[0.,0.8]


externalPowerSpectrumFunction[0.,200., externalSpectraType->"pk"]


externalPowerSpectrumFunction[0.,200., externalSpectraType->"pknw"]


pthetathetaMoments[0.1,1]


sigmavNLNew[0.1,0.2,0.9]


integratetrapz


observedPowerSpectrum[0.2,0.4,0.8]


observedPowerSpectrum[0.9,0.4,0.8]//AbsoluteTiming


kminn=$kminRef


kmaxx=kmaxChoice[1.0]


posis=Position[$krangePowerSpectrumFiducial,_?(kminn<=#<kmaxx&)];


karray=Extract[$krangePowerSpectrumFiducial, posis];


karray//Dimensions


fishArra=$fisherCosmoParsBlock[1.0,#,0.1]&/@karray


fishArra//Dimensions


$paramoptions


fishInteg=Transpose[{karray,fishArra[[All,1,1]]}];


fishIntegFunc[zz_?NumericQ,muu_?NumericQ, a_?IntegerQ, b_?IntegerQ]:=($fisherCosmoParsBlock[zz,#,muu][[a,b]])&/@karray


observedPowerSpectrum[1.2,0.1,0.9]//AbsoluteTiming


0.06*450


powerSpectrum[1.2,0.1]//AbsoluteTiming


dObsPowerDpi[1.2,0.1,0.44, 3, dlnPdpDerivative->True]//AbsoluteTiming


($fisherCosmoParsBlock[1.2,#,0.9][[1,1]]&/@karray)//AbsoluteTiming


fishIntegFunc[1.2,0.05, 1, 1]//AbsoluteTiming


fishIntegFunc[1.2,0.7, 1, 1]//AbsoluteTiming;


fishInteg;


fishInteg;


ListLogLinearPlot[fishInteg]


(Differences[#1].MovingAverage[#2, 2] & @@ Transpose[fishInteg])//AbsoluteTiming


interpo=Interpolation[fishInteg, InterpolationOrder->1]


NIntegrate[interpo[xx], {xx,kminn,kmaxx}]//AbsoluteTiming


NIntegrate[(($fisherCosmoParsBlock[1.0,kk,0.1][[1,1]])),{kk,kminn,kmaxx}]


NIntegrate[(($fisherCosmoParsBlock[1.0,kk,muu])),{kk,kminn,kmaxx},{muu,0.,1.}]


NIntegrate[(($fisherCosmoParsBlock[1.0,kk,muu]))*volumeEffective[1.0,kk,muu]*kk^2*dampingTerm[1.0,kk,muu],{kk,kminn,kmaxx},{muu,0.,1.}]*volumeSurvey[1.0,$zbinGlobal]*(2/(8 Pi^2))



