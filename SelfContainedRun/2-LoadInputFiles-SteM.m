(* ::Package:: *)

SetOptions[$Output, FormatType->OutputForm];
nbdir=If[aa=DirectoryName@$InputFileName; aa=="", NotebookDirectory[], aa]
Get[nbdir<>"/GenerateInputFiles/"<>"GenerationParameters.m"]

inputDataDir=findInputFilesDir[$nameOfModel];
Print["Directory of input files exists: ", DirectoryQ[inputDataDir]];

$EHuDeWiggling=False;

$paramfidus


$parameterFiducialsRule=Thread[Rule[$paramlabels, $paramfidus]]


$parameterFiducials[pp_]:=Block[{pv}, pv=(pp/.$parameterFiducialsRule); Return[pv]]


$parameterFiducials["s8"]


epsstrfun[0]


$epslist = {0.00625, 0.0125, 0.01875, 0.025, 0.0375, 0.05, 0.10}; 


$epslistPM = Sort[(-1*$epslist)~Join~$epslist~Join~{0.}]


$epsilonstep=0.0125


setNumericalEpsilon[$epsilonstep, epsilonStepForZdependentParameters->$epsilonzstep]


(*epsstrfun[#]&/@$epslist*)


(*epsstrRule=Thread[Rule[$epslist, (epsstrfun[#]&/@$epslist)]]*)


(*paramsFolderFunc[param_, epsilon_]:=Block[{pm, string, parstring, epsstr, epsi}, 
                                          pm=Which[Sign[epsilon]==-1, "_"<>plusminstr[[1]], Sign[epsilon]==1, "_"<>plusminstr[[2]], Sign[epsilon]==0, ""];
                                          epsi=If[Chop[epsilon, 10^-6]==0., 0, epsilon];
                                          If[param=="fiducial" && epsi!=0, Return[None]];
                                          If[param!="fiducial" && epsi==0, Return[None]];
                                          If[param=="fiducial", pm=""];
                                          epsstr=epsstrfun[Abs[epsi]];
										  string=(param<>pm<>"_"<>epsstr);						
                                          Return[string]
                                          ];*)


paramsFolderFunc["h", 0.00625]


paramsFolderFunc["fiducial", 0.555]


paramsFolderFunc["gamma", 0]


ClearAll[paramsFolderTable]


$folderpars


paramsFolderList = DeleteCases[Flatten[Table[paramsFolderFunc[parr, epps], {parr, $folderpars}, {epps, $epslistPM}]],None];


paramsFolderList


paramsFolderList//Length


Length@$epslistPM * Length@$folderpars - (Length@$epslistPM-1)*1 - Length@$paramlabels


$paramsDirectoryNames=paramsFolderList;


inputDataDir


({#,DirectoryQ@(inputDataDir<>#)})&/@$paramsDirectoryNames


AllTrue[Flatten[#], TrueQ]&@((DirectoryQ@(inputDataDir<>#))&/@$paramsDirectoryNames)


(*?$paramsDirectoryNames*)


$paramsDirectoryNames[[1]]


$ignorepars={"AIA","etaIA","betaIA"};
(*ignoreparamsfunc[str_, ignolist_:ignorepars]:=If[MemberQ[paramsFolderTable[ignolist], str], Return[$paramsDirectoryNames[[1]] ], Return[str]]*)


(*(***(*order the params directory names according to the order of $paramoptions *)**)*)


$paramoptions


Print["Names of input Directories : ", $paramsDirectoryNames];


funcsInputNames={"background_Hz","D_Growth-z","f_GrowthRate-z","Plin-zk","sigma8-z"}


$observablesNames={"Hz","Dgz","fgz","Plzk","s8z"};


(*$zrangeFile="z_values_list.txt"*)


(*$krangeFile="k_values_list.txt"*)


inputDataDir


frontendversion=replaceStringPoint[$VersionNumber, "math"];  (* Use the mathematica version number, in case of compatibility issues *)


interpolationDirectory=mkDirectory[inputDataDir<>"/InterpolationFunctions-"<>frontendversion<>"/"]


$deleteInterpolationFiles


If[$deleteInterpolationFiles==True,
Print[" Deleting interpolation files..."];
interpfiles=FileNames[interpolationDirectory<>"*.mc"];
Print["Number of filenames: ", Length@interpfiles];
DeleteFile[#]&/@interpfiles;
Print["Remaining files in folder: "];
FileNames[interpolationDirectory<>"*.mc"],
(*ELSE*)

Print[" Reading interpolation files..."];
Print["Number of filenames: "];
Length@FileNames[interpolationDirectory<>"*.mc"]
]


interpInputTabs="LinLin";


SetOptions[Interpolation, InterpolationOrder->interpOrd, Method->interpMeth]


intpmeth=(FilterRules[Options[Interpolation], Method])[[1,2]]


interpolParticularString="IntpMeth_"<>intpmeth<>"-IntpOrd_"<>(ToString@interpOrd)<>"-IntpFunc_"<>interpInputTabs


interpInputTabs


Switch[interpInputTabs
,
"Log10Log10",
intpFuncLogBool=True;
funcx=Log10;
funcy=Log10;
,
"LinLin",
intpFuncLogBool=False;
funcx=Identity;
funcy=Identity;
,
_,
intpFuncLogBool=False;
funcx=Identity;
funcy=Identity;
];


MemberQ[$ignorepars,"h"]


MemberQ[$ignorepars,"betaIA"]


$paramsDirectoryNames


funcsInputNames


SetOptions[fGrowthRate,solveDifferentialEquation->False]


$rulesfuncsInputNames=Thread[Rule[$observablesNames,funcsInputNames]]


ignoreParameters[parameps_]:=Block[{ignolist=$ignorepars, return}, 
                    return=If[MemberQ[ignolist, StringSplit[parameps, "_"][[1]]]==True, 
                                                       paramsFolderFunc["fiducial", 0],
                                                       parameps]
                                                       ]


ignoreParameters["fiducial_eps_0"]


ignoreParameters["betaIA_pl_eps_1p0E-1"]


ignoreParameters["AIA_pl_eps_1p0E-1"]


ClearAll[externalInputFile]


externalInputFile[par_, observable_, format_:".txt"]:=Block[{para,obsstr},
para=ignoreParameters[par];
If[MemberQ[$observablesNames,observable]==False, Print["Observable string is not supported"]; Abort[]];
obsstr=observable/.$rulesfuncsInputNames;
((ToString[para]<>"/"<>obsstr<>format))
]


kzRangeFile[par_, grid_String, format_:".txt"]:=Block[{para, gfile},
para=ignoreParameters[par];
gfile=grid<>"_values_list";
((ToString[para]<>"/"<>gfile<>format))
]


externalInputFile["ns_mn_eps_1p9E-2", "Plzk"]


kzRangeFile["ns_mn_eps_1p9E-2", "z"]


externalInputFile["AIA_pl_eps_1p0E-1", "Plzk"]


$observablesNames


AllTrue[Flatten[#], TrueQ]&@(FileExistsQ@(inputDataDir<>externalInputFile[#, "Plzk"])&/@$paramsDirectoryNames)


Table[(AllTrue[Flatten[#], TrueQ]&@(FileExistsQ@(inputDataDir<>externalInputFile[#, ii])&/@$paramsDirectoryNames)),{ii,$observablesNames}]


paramsFolderFunc["fiducial", 0.]==$paramsDirectoryNames[[1]]


$zrangePowerSpectrum = Import[(inputDataDir<>kzRangeFile[paramsFolderFunc["fiducial", 0.], "z"]), "List"];


$krangePowerSpectrumFiducial=Import[(inputDataDir<>kzRangeFile[paramsFolderFunc["fiducial", 0.], "k"]), "List"];


Print["Length of input k-range: ", Length@$krangePowerSpectrumFiducial]


Print["Length of input z-range: ", Length@$zrangePowerSpectrum]


ClearAll[externalMatrixListImport]


externalMatrixListImport[obsname_,param_,epsilon_]:=externalMatrixListImport[obsname,param,epsilon]=Block[{folderpar, tab, dim},
folderpar=paramsFolderFunc[param, epsilon];
tab=Import[(inputDataDir<>externalInputFile[folderpar, obsname]), "Table"];
dim=Dimensions[tab];
If[dim[[2]]==1,
tab=Flatten[tab,1];  (*If imported Table is just a list, flatten to have a 1-d List*)
];
Return[tab]]

$folderpars






ClearAll[externalObsMatrixList,externalObsTable,externalObservableFunc]








ClearAll[externalObservableFunc]
ClearAll[externalObservableDerivative2ptFunc,externalObservableDerivativeStem1Func,externalObservableDerivativeStem4Func]


$paramlabels


$observablesNames


(*myInterrupt[]*)





externalObservableDerivative2ptFunc["Hz","OmegaM"][0.2]


externalObservableDerivative2ptFunc["Plzk","OmegaM"][0.1,0.3]


lala=$zeroFunction


lala[0.1,0.3]


$paramlabels





importBools={}
Do[
Print["----- Observable: ", obse];
Print["Importing fiducials"];
interpolationNameFiducial=interpolationDirectory<>"interpolatingFunction_Fiducial_"<>obse<>".mc";
AppendTo[importBools,Quiet@Check[externalObservableFunc[obse,"fiducial", 0.]=Uncompress@Import[interpolationNameFiducial,"String"]; True, False]];
Do[
Print["Parameter: ", par];
Print["Importing 2pt derivatives"];
interpolationNameParam2=interpolationDirectory<>"interpolatingFunction_"<>obse<>"_Derivative_2pt_"<>par<>".mc";
AppendTo[importBools,Quiet@Check[externalObservableDerivative2ptFunc[obse,par]=Uncompress@Import[interpolationNameParam2,"String"]; True, 
If[MemberQ[$ignorepars, par]==True,
Print["Parameter ", par, " will be ignored and derivative taken to be zero"];
externalObservableDerivative2ptFunc[obse,par]=$zeroFunction;
True,
(*ELSE*)
Print["Not found: ", interpolationNameParam2];
False]]];

If[StringMatchQ[obse,"*zk*"]==True,
Print["Importing SteM derivatives"];
(*interpolationNameParam1=interpolationDirectory<>"interpolatingFunction_"<>obse<>"_Derivative_Stem1_"<>par<>".mc";
AppendTo[importBools,Quiet@Check[externalObservableDerivativeStem1Func[obse,par]=Uncompress@Import[interpolationNameParam1,"String"]; True, Print["Not found: ", interpolationNameParam1]; False]];*)
interpolationNameParam4=interpolationDirectory<>"interpolatingFunction_"<>obse<>"_Derivative_Stem4_"<>par<>".mc";
AppendTo[importBools,Quiet@Check[externalObservableDerivativeStem4Func[obse,par]=Uncompress@Import[interpolationNameParam4,"String"]; True, 
If[MemberQ[$ignorepars, par]==True,
Print["Parameter ", par, " will be ignored and derivative taken to be zero"];
externalObservableDerivativeStem4Func[obse,par]=$zeroFunction;
True,
(*ELSE*)
Print["Not found: ", interpolationNameParam2];
False]]]];
,   
{par, $paramlabels} ]
,
{obse,  $observablesNames}]


importBools


AllTrue[importBools, TrueQ]


If[AllTrue[importBools, TrueQ]==True,
Print["All interpolation files found"];
$importNumericalTables=False,
(*Else*)
$importNumericalTables=True;];

Print["Importing Numerical Tables for observables: " , $importNumericalTables]


If[(*True*)$importNumericalTables,
numfiles=0;
Print["-------Memory in Use: ", memoryUseMB[]];
Do[
inidate=DateString[];
Print["**Param: ", par];
Do[
folderpar = paramsFolderFunc[par, epsi];
If[SameQ[folderpar, None],
(*Parameter and epsilon combination does not exist, continue to next step of the loop*)
Continue[]
];
Print["****Epsilon: ", epsi];
Print["Pre-table: ------Memory in Use: ", memoryUseMB[]];

krange[par,epsi]=Import[(inputDataDir<>kzRangeFile[folderpar, "k"]), "List"];
Do[
Print["Observable being read: ", (obsname/.$rulesfuncsInputNames)];
externalObsMatrixList[obsname,par,epsi]=externalMatrixListImport[obsname,par,epsi];

If[StringMatchQ[obsname, "*zk*"]==True,
Print["Scale-dependence and z-dependence"];
externalObsTable[obsname, par, epsi]=Flatten[Table[{{$zrangePowerSpectrum[[zzi]],krange[par,epsi][[kki]]},externalObsMatrixList[obsname,par,epsi][[zzi,kki]]},
{zzi,1,Length@$zrangePowerSpectrum},{kki,1,Length@krange[par,epsi]}],1];
, (*Else: Only z-dependence*)
Print["Only z-dependence"];
externalObsTable[obsname, par, epsi]=Table[{$zrangePowerSpectrum[[zzi]],externalObsMatrixList[obsname,par,epsi][[zzi]]},
{zzi,1,Length@$zrangePowerSpectrum}];
];

externalObservableFunc[obsname, par, epsi]=Interpolation[externalObsTable[obsname, par, epsi]];


Print["------Memory in Use: ", memoryUseMB[]];
(*extGrowthRateFunc[par]=Interpolation[extGrowthRateTable[par]];*)
ClearAll[externalObsTable];
ClearAll[externalObsMatrixList];
numfiles += 1;
,
{obsname,$observablesNames}]
,
{epsi, $epslistPM}
];
finidate=DateString[];
timeElapsed[inidate, finidate];
,
{par,$folderpars}];
Print["Observable Files Read: "<>ToString@numfiles];
]


(*myInterrupt[]*)


memoryUseMB[]


(*externalObservableFunc["Plzk", "h", -0.1]*)


LogLogPlot[{externalObservableFunc["Plzk", "fiducial", 0.][0.,kk] , externalObservableFunc["Plzk", "fiducial", 0.][0.,kk]}, {kk, 10^-4, 5}]


obsername="Plzk"
LogLogPlot[{externalObservableFunc[obsername,"fiducial", 0.][1.0,kk]} , {kk, 10^-4, 5}]


obsername="Plzk"
LogLogPlot[Evaluate@Table[externalObservableFunc[obsername,"fiducial", 0.][zz,kki],{kki,{0.001,0.01,0.1,0.5,1.0,5.0}}] , {zz, 10^-4, 2.5}]



observableParamDependency[obsname_,par_,zzi_,kki_,epsilonlist_:$epslistPM]:=Block[{varfun, fiducialfun, paramFuncTab, var, fiduval},
            If[StringMatchQ[obsname, "*zk*"]==True,
            (*Scale-dependence and z-dependence*)
            varfun = externalObservableFunc[obsname,par, #][zzi,kki]&;            
            fiducialfun=externalObservableFunc[obsname,"fiducial", 0.][zzi,kki];
            ,
            varfun = externalObservableFunc[obsname,par, #][zzi]&;            
            fiducialfun=externalObservableFunc[obsname,"fiducial", 0.][zzi];
            ];
            fiduval=$parameterFiducials[par];
            paramFuncTab = Table[{fiduval*(1+ee), 
                            If[ee==0.,
                            var=fiducialfun,
			                var=varfun[ee] ];
            var/fiducialfun}, {ee, epsilonlist}];
            Return[paramFuncTab]]


ClearAll[steMderiv]


ClearAll[twoPointDerivative]


$epsilonstep


$epslist


twoPointDerivative[obsname_,param_,zzi_,kki_, epsi_:0.00625, relative_:True]:=Block[{ee, twopoint, relder, retu, varfun,fiducialfun, fiduval}, 
If[StringMatchQ[obsname, "*zk*"]==True,
            (*Scale-dependence and z-dependence*)
            varfun = externalObservableFunc[obsname,param, #][zzi,kki]&;            
            fiducialfun=externalObservableFunc[obsname,"fiducial", 0.][zzi,kki];
            ,
            varfun = externalObservableFunc[obsname,param, #][zzi]&;            
            fiducialfun=externalObservableFunc[obsname,"fiducial", 0.][zzi];
            ];
fiduval=$parameterFiducials[param];
twopoint=(varfun[epsi]-varfun[-epsi])/(2*epsi*fiduval);
relder=twopoint/fiducialfun;
retu=If[relative==True,
relder,
twopoint];
Return[retu]]

steMderiv[partable_,x0_, order_:1]:=Block[{lm,rmse,relativederiv,halflen, cutptab, tolerance=0.005, polyst},
            If[order<1, Print["Polynomial Order must be larger or equal 1"];Abort[]];
             polyst=Table[x^ii,{ii,Range[0,order]}];
             halflen=Range[Floor[Length[partable]/2]];
             Do[
                cutptab=partable[[ii;;-ii]];
				lm=LinearModelFit[cutptab, polyst,x];
				rmse=Sqrt[Total[(lm["FitResiduals"])^2]];
				relativederiv=(D[lm[x],x]/.x->x0);
				If[rmse <= tolerance, 
				Break[]]
				,
				{ii,halflen}];
			Return[{relativederiv, rmse, lm, cutptab}]
			]



ClearAll[observableDerivativeZK]


Options[observableDerivativeZK]={derivativeType->"SteM1", twopointepsilon->$epsilonstep};


derivativeType::usage="Option for observableDerivativeZK. Options: '2pt'->Two point derivative, 'Stem1'->SteM derivative fitted with order 1, 'Stem4': SteM derivative fitted with order 4.";


twopointepsilon::usage="Option for observableDerivativeZK. Value of the epsilon step used in the 2-point derivatives.";


observableDerivativeZK[obsname_,parname_,zzi_,kki_, opts:OptionsPattern[]]:=Block[{fiduobs,fiduval,pptab,rmse=-1.,lml,cutpartab,
																	twoptder, epsvalue,
                                                                    zgrid=$zrangePowerSpectrum, zlenhalf,
                                                                    kgrid=$krangePowerSpectrumFiducial,klenhalf,kindex,
                                                                    dertype, relderiv, deriv,
                                                                    derivElem},
             epsvalue=OptionValue[twopointepsilon];                                                       
            dertype=OptionValue[derivativeType];
            If[MemberQ[{"SteM1","SteM4","2pt"},dertype]==False, Abort[]];
            klenhalf=Floor[(Length@kgrid)/2];
            zlenhalf=Floor[(Length@zgrid)/2];
             If[StringMatchQ[obsname, "*zk*"]==True,
            (*Scale-dependence and z-dependence*)
            fiduobs=externalObservableFunc[obsname,"fiducial", 0.][zzi,kki];
            kindex=kki;
            ,
            (*Only z-dependence*)
            fiduobs=externalObservableFunc[obsname,"fiducial", 0.][zzi];
            kindex = kgrid[[klenhalf]] (*have a fake kki*)
            ];
            Which[dertype=="2pt",
            relderiv=twoPointDerivative[obsname,parname,zzi,kki, epsvalue, True];
            ,
            dertype=="SteM1",
            pptab=observableParamDependency[obsname,parname,zzi,kki];
            fiduval=$parameterFiducials[parname];
            {relderiv,rmse,lml,cutpartab}=steMderiv[pptab, fiduval];
            ,
            dertype=="SteM4",
            pptab=observableParamDependency[obsname,parname,zzi,kki];
            fiduval=$parameterFiducials[parname];
            {relderiv,rmse,lml,cutpartab}=steMderiv[pptab,fiduval,4];
            ];
           If[zzi==zgrid[[zlenhalf]] && kindex==kgrid[[klenhalf]],
				Print["z=", zzi];
				Print["k=", kki];
				Print["RMSE residual: ", rmse];  (*Will print -1 for 2pt deriv*)
				Print["Relative "<>dertype<>" tangent at fiducial: ", relderiv];
				];
			 If[StringMatchQ[obsname, "*zk*"]==True,
            (*Scale-dependence and z-dependence*)
			 derivElem = {{zzi,kki},relderiv*fiduobs} ;
			 ,
			 (*Only z-dependence*)
			 derivElem = {zzi,relderiv*fiduobs} ;
			];
			Return[derivElem]
			];


$epslist



$paramlabels


If[$importNumericalTables,  (*Only execute this cell if interpolating Functions not found above*)
globaliniT=DateString[];
paramlist=$paramlabels;(*[[1;;1]];*)
(*observablesNamesZ={"Hz","s8z"};
observablesNamesZK={"Dgzk","fgzk","Plzk","Pnlzk"};;*)
Do[
Print["----- Observable: ", obse];
ii=DateString[];
Do[
Print["Parameter: ", par];
dataDerivatives[par]=Reap[
zgrid=$zrangePowerSpectrum;
kgrid=$krangePowerSpectrumFiducial;
Do[
If[MemberQ[$zrangePowerSpectrum[[1;;-1;;10]], zzii]==True,
Print["Redshift z=",zzii]];
 Do[
     Sow[observableDerivativeZK[obse,par,zzii,kkii, derivativeType->"2pt"],"2pt"];
     Sow[observableDerivativeZK[obse,par,zzii,kkii, derivativeType->"SteM1"],"SteM1"];
      Sow[observableDerivativeZK[obse,par,zzii,kkii, derivativeType->"SteM4"],"SteM4"];
      If[StringMatchQ[obse, "*zk*"]==False,
      Break[]]; (*Continue with next z*)
  ,   {kkii,  kgrid}]           
  ,   {zzii, zgrid}]
  ];
  externalObservableDerivative2ptFunc[obse,par]=Interpolation[dataDerivatives[par][[2,1]], InterpolationOrder->1, Method->"Spline"];
  externalObservableDerivativeStem1Func[obse,par]=Interpolation[dataDerivatives[par][[2,2]], InterpolationOrder->1, Method->"Spline"];
  externalObservableDerivativeStem4Func[obse,par]=Interpolation[dataDerivatives[par][[2,3]], InterpolationOrder->1, Method->"Spline"];
  ClearAll[dataDerivatives];
  ,   {par, paramlist} ];
  Print["------Memory in Use: ", memoryUseMB[]];
  ff=DateString[];
  timeElapsed[ii,ff];
  ,
  {obse,  $observablesNames}];
Print["***Total Time Elapsed: "];
timeElapsed[globaliniT,DateString[]];

fufuRule=Thread[Rule[{externalObservableDerivative2ptFunc,externalObservableDerivativeStem1Func,externalObservableDerivativeStem4Func}, {"_Derivative_2pt_","_Derivative_Stem1_","_Derivative_Stem4_"}]];
]


If[$importNumericalTables,
Print["*** EXPORTING Interpolating Functions ***"];
Do[
Print["----- Observable: ", obse];
interpolationNameFiducial=interpolationDirectory<>"interpolatingFunction_Fiducial_"<>obse<>".mc";
If[StringMatchQ[obse,"*zk*"]==True,
table=Flatten[Table[{{zzi,kki},externalObservableFunc[obse,"fiducial", 0.][zzi,kki]},
{zzi,$zrangePowerSpectrum},{kki,$krangePowerSpectrumFiducial}],1];
,
table=Table[{zzi,externalObservableFunc[obse,"fiducial", 0.][zzi]},
{zzi,$zrangePowerSpectrum}];
];
Export[interpolationNameFiducial,Compress[Interpolation[table]],"String"];
Do[
Print["Parameter: ", par];
Do[
Print["Derivative: ", ToString@fufu];
If[StringMatchQ[obse,"*zk*"]==True,
If[Quiet@Check[fufu[obse,par][$zrangePowerSpectrum[[1]],$krangePowerSpectrumFiducial[[1]] ], False]==False,
Print[ToString@fufu," for observable ", obse, "not exported"];
Break[];
];
table=Flatten[Table[{{zzi,kki},fufu[obse,par][zzi,kki]},
{zzi,$zrangePowerSpectrum},{kki,$krangePowerSpectrumFiducial}],1];
, (*ELSE, z-dependence only *)
If[Quiet@Check[fufu[obse,par][$zrangePowerSpectrum[[1]] ], False]==False,
Print[ToString@fufu," for observable ", obse, "not exported"];
Break[];
];
table=Table[{zzi,fufu[obse,par][zzi]},
{zzi,$zrangePowerSpectrum}];
];
derivname=(fufu/.fufuRule);
interpolationNameParam="interpolatingFunction_"<>obse<>derivname<>par<>".mc";
Export[interpolationDirectory<>interpolationNameParam,Compress[Interpolation[table]],"String"];
,
{fufu,{externalObservableDerivative2ptFunc,externalObservableDerivativeStem4Func, externalObservableDerivativeStem1Func}}]
(*interpolationNameParam1="interpolatingFunction_"<>obse<>"_Derivative_Stem1_"<>par<>".mc";
Export[interpolationDirectory<>interpolationNameParam1,Compress[externalObservableDerivativeStem1Func[obse,par]],"String"];*)
(*interpolationNameParam4="interpolatingFunction_"<>obse<>"_Derivative_Stem4_"<>par<>".mc";
Export[interpolationDirectory<>interpolationNameParam4,Compress[externalObservableDerivativeStem4Func[obse,par]],"String"];*)
,   
{par, paramlist} ]
,
{obse,  $observablesNames}];
Print["*** Interpolating Functions Exported***"];
]


$zrangePowerSpectrum


$krangePowerSpectrumFiducial;


externalObservableDerivative2ptFunc["Hz","OmegaM"][0.4]


externalObservableDerivativeStem4Func["Plzk","OmegaM"]


Table[
LogLinearPlot[{Identity[externalObservableDerivative2ptFunc["Hz",pp][zaz]]}, {zaz, 0.001,2.4}, PlotLegends->Automatic,PlotLabel->pp, PlotRange->Full], {pp,$paramlabels}]


Table[
obsstr="fgz";
LogLinearPlot[{Identity[externalObservableDerivative2ptFunc[obsstr,pp][zaz]],
Identity[externalObservableDerivativeStem4Func[obsstr,pp][zaz]],Identity[externalObservableDerivativeStem4Func[obsstr,pp][zaz]]}, {zaz, 0.001,2.4}, PlotLegends->Automatic,PlotLabel->pp, PlotRange->Full], {pp,$paramlabels}]


legends={"2pt", "SteM"};


zeto=1.0;
Table[
paramstr=pp;
obsstr="Plzk";
LogLinearPlot[{Identity[externalObservableDerivative2ptFunc[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]],
Identity[externalObservableDerivativeStem4Func[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]}, {kak, 1*^-4,10}, PlotLegends->legends,PlotLabel->pp, PlotRange->Full],
{pp,$paramlabels}]



zeto=1.0;
Table[
paramstr=pp;
obsstr="Plzk";
LogLinearPlot[{Identity[externalObservableDerivative2ptFunc[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]],
Identity[externalObservableDerivativeStem4Func[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]}, {kak, 1*^-4,10}, PlotLegends->legends,PlotLabel->pp,PlotRange->Full],
{pp,$paramlabels}]


legends={"2pt", "SteM1", "SteM4"};


zeto=1.5;
Table[
paramstr=pp;
obsstr="Plzk";
LogLinearPlot[{Identity[externalObservableDerivative2ptFunc[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]],
Identity[externalObservableDerivativeStem1Func[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]],
Identity[externalObservableDerivativeStem4Func[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]}, {kak, 1*^-4,10}, PlotLegends->legends,PlotLabel->pp,PlotRange->Full],
{pp,$paramlabels}]


externalObservableDerivative2ptFunc["Plzk","AIA"][zeto,0.1]


externalObservableDerivativeStem4Func["Plzk","AIA"][zeto,0.1]


zeto=1.0;
Table[
paramstr=pp;
obsstr="Plzk";
LogLinearPlot[{Identity[externalObservableDerivative2ptFunc[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]],
Identity[externalObservableDerivativeStem1Func[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]],
Identity[externalObservableDerivativeStem4Func[obsstr,paramstr][zeto,kak]/externalObservableFunc[obsstr,"fiducial", 0.][zeto,kak]]}, {kak, 1*^-4,10}, PlotLegends->legends,PlotLabel->pp],
{pp,$paramlabels}]



$paramlabels


$pks8ratio


?$pks8ratio


If[$pks8ratio==0, False, True]


intpFuncLogBool


Print["Power spectrum P(k)/sigma8 (yes=1, no=0): ", $pks8ratio];
setPowerSpectrumRelatedSpecs[pks8Ratio->If[$pks8ratio==0, False, True]]

Print[intpFuncLogBool];
(*Options@setPScosmoOptions;
setPScosmoOptions[kdependentGrowth->True,lcdmBool->False,linearBool->False];*)

SetOptions[externalPowerSpectrumFunction, functionInterpolatedInLogLog->intpFuncLogBool];


funco=(externalObservableDerivative2ptFunc["Hz",#]&)


funco["OmegaM"][0.1]


(*myInterrupt[]*)


setCosmologyFixedValues[internalHubbleUnits->"1/Mpc", externalHubbleUnits->"1/Mpc",
internalDistanceUnits->"Mpc",
internalPkUnits->"1/Mpc", externalPkUnits->"h/Mpc"]


Options[setCosmologyFixedValues]


$externalPkUnits


$internalDistanceUnits


externalObservableFunc["Hz","fiducial", 0.][0.]*$lightspeed


$lightspeed


setExternalCosmoInterpolatingFunctions[externalHubbleInput->(externalObservableFunc["Hz","fiducial", 0.]),externalHubbleDerivativesInput->(externalObservableDerivative2ptFunc["Hz",#]&)]


(*Options[Hubble]*)


(*Hubble[1.0, lcdmBool->True, externalFile->False]*)


(*Hubble[1.0]*)


(*Hubble[0., physicalUnits->True]*)


(*Hubble[0., lcdmBool->True, externalFile->False, physicalUnits->True]*)


(*Hubble[0., hubble->0.68, physicalUnits->True]*)


ztesta=0.


(*(Hubble[ztesta, hubble->0.67*(1+0.01), lcdmBool->True, externalFile->False, physicalUnits->False]-Hubble[ztesta, hubble->0.67*(1-0.01), lcdmBool->True, externalFile->False, physicalUnits->False])/(2*0.01*0.67)*)


Hubble[0.]


Hubble[0., physicalUnits->True]


$externalPkUnits


$internalPkUnits


externalObservableDerivative2ptFunc["Hz","gamma"][ztesta]


$paramoptions


Options[setExternalCosmoInterpolatingFunctions]


$observablesNames


setExternalCosmoInterpolatingFunctions[externalGrowthInput->(externalObservableFunc["Dgz","fiducial", 0.]),externalGrowthDerivativesInput->(externalObservableDerivativeStem4Func["Dgz",#]&),
externalSigma8ofZInput->(externalObservableFunc["s8z","fiducial", 0.]), externalSigma8ofZDerivativesInput->(externalObservableDerivative2ptFunc["s8z",#]&),
externalGrowthRateInput->(externalObservableFunc["fgz","fiducial", 0.]), externalGrowthRateDerivativesInput->(externalObservableDerivativeStem4Func["fgz",#]&)]


Options[setExternalCosmoInterpolatingFunctions]


$kdependentGrowth


$fkfix


Options[fGrowthRate]


Options[externalGrowthRateFunction]


setExternalCosmoInterpolatingFunctions[externalPowerSpectrumInput->(externalObservableFunc["Plzk","fiducial", 0.]),
externalPowerSpectrumDerivativesInput->(externalObservableDerivativeStem4Func["Plzk",#]&)]


powerSpectrum[0.1,0.5, Omegam->0.33]


$zrangePowerSpectrum


setKandZMinMaxFixedValues[kmaxIntegrations->Max[$krangePowerSpectrumFiducial], kminIntegrations->Min[$krangePowerSpectrumFiducial], kmaxInterpolations->Max[$krangePowerSpectrumFiducial],
zmaxIntegrations->($externalHubbleInterpolatingFunction["Domain"]//Max)]


$zMaxIntegral

