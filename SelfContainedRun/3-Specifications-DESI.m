(* ::Package:: *)

(* ::Section:: *)
(*Load survey specifications and set other survey parameters*)


$paramoptions


specsDir=fishertoolsdir<>"/Euclid+SKA1+SKA2+DESI-specifications/";


If[DirectoryQ[specsDir], Print["Specifications directory exists"]]

$surveychosen="DESI"

specsFiles["DESI-ndens"]="DESI-ELG-dN_dzddeg-zbins-deltaz_0p1.txt";  
specsFiles["DESI-bias"]="DESI-ELG-fiducial_bias_factor.txt";


If[FileExistsQ[specsDir<>specsFiles["DESI-ndens"]], Print["DESI specifications for nofz: file exists"]]
If[FileExistsQ[specsDir<>specsFiles["DESI-bias"]], Print["DESI specifications for bias: file exists"]]





readSurveySpecsFromFiles[numberDensityZbinsFile->specsFiles["DESI-ndens"],
biasSpecsZbinsFile->specsFiles["DESI-bias"],setGalaxySurveySpecifications->True,
specificationsDirectory->specsDir,surveyName->"DESI", numberDensityFileFormat->"CentralZinBin-dndOmegadz", 
biasFileFormat->"BiasFactor-Growth"]


100.0/1700


$zbinGlobal


$zmeansurvey


setGalaxySurveySpecs[zMeanSurvey->1.2]


$zmeansurvey


zb1=(zaverage@$zbinGlobal)[[1]]


zb2=(zaverage@$zbinGlobal)[[-1]]


Print["Test numeric output from specs: nz[zbinstart]="<>ToString[$ndensGalZFuncGlobal[zb1]]<>",  bz[zbinend]="<>ToString[$biasInterpFunc[zb2]] ]

If[$debugPlots==False, Plot[{$ndensGalZFuncGlobal[zz] ,$biasInterpFunc[zz]/1000}, {zz,zb1,zb2}, PlotLegends->Automatic]]


ellmxlin=1500;
ellmxnonlin=5000;


$elllist={ellmxlin,ellmxnonlin}


$invertkmxscales=False;


$runtype


$photozparameters


setWLBinsSpecifications[ellmin->10,ellmax->ellmxnonlin, numberOfellBins->100, minZinSurvey->0.001, maxZinSurvey->2.5, numberOfZBins->10]


$surveychosen


setWeakLensingSurveySpecifications[$surveychosen, activateIntrinsicAlignments->$includeIAterms]


$luminosityFileInterpolation=Interpolation[Import[specsDir<>"scaledmeanlum-E2Sa.dat", "Table"], InterpolationOrder->1];
