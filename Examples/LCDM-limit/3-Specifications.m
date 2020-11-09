(* ::Package:: *)

(* ::Section:: *)
(*Load survey specifications for GCsp*)


$specsDirectory=fishertoolsdir<>"/Euclid+SKA1+SKA2+DESI-specifications/";


If[DirectoryQ[$specsDirectory], Print["Specifications directory exists"]]


setGCspectroSurveySpecifications[$surveychosen]


setGalaxySurveySpecs[zMeanSurvey->$zmeansurvey]


Print["Test numeric output from specs: nz[zbinstart]="<>ToString[$ndensGalZFuncGlobal[(zaverage@$zbinGlobal)[[1]]]]<>",  bz[zbinend]="<>ToString[$biasInterpFunc[(zaverage@$zbinGlobal)[[-1]]]] ]


If[$debugPlots==False, Plot[{$ndensGalZFuncGlobal[zz] ,$biasInterpFunc[zz]/1000}, {zz,(zaverage@$zbinGlobal)[[1]],(zaverage@$zbinGlobal)[[-1]]}, PlotLegends->Automatic]]


(* ::Section:: *)
(*Load survey specifications and set other survey parameters*)


ellmxlin=1500;
ellmxnonlin=5000;


$elllist={ellmxlin,ellmxnonlin}


$invertkmxscales=False;


$runtype


$photozparameters


setWLBinsSpecifications[ellmin->10,ellmax->ellmxnonlin, numberOfellBins->100, minZinSurvey->0.001, maxZinSurvey->2.5, numberOfZBins->10]


$surveychosen


setWeakLensingSurveySpecifications[$surveychosen, activateIntrinsicAlignments->$includeIAterms]


Check[existsFile[$specsDirectory<>"scaledmeanlum-E2Sa.dat"]; 
$luminosityFileInterpolation=Interpolation[Import[$specsDirectory<>"scaledmeanlum-E2Sa.dat", "Table"], InterpolationOrder->1];
Print["IA specifications loaded"], 
Print["IA specification file not found"]
]
