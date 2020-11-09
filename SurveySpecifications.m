(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: SurveySpecifications *)
(* :Context: SurveySpecifications` *)
(* :Author: casas *)
(* :Date: 2016-09-08 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2016 casas *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["SurveySpecifications`"]
(* Exported symbols added here with SymbolName::usage *)


Needs["CosmologyFunctions`"]
Needs["UsefulTools`"]
Needs["FisherTools`"]

EndPackage[]


BeginPackage["SurveySpecifications`", {"FisherTools`", "CosmologyFunctions`", "UsefulTools`"}]


biasOfzFuncSKAFit::usage="fitting function for bias in SKA. biasOfzSKAfit[z, cvector]"

dndzFuncSKAFit::usage=" fitting function for dndz in SKA. dndzFuncSKAFit[z, cvector]"

ndensOfzFuncSKAFit::usage="number density function for SKA in units of [1/distanceAngular^3]. Form: ndensOfzFuncSKAFit[z_Real, cvect_List]."

setNumbDensBiasSKA::usage="sets the number density and the bias for a SKA survey, given the redshift bins and the vector of c-parameters. Form: setNumbDensBiasSKA[zbins_List, cvect_List]"
importSKAspecs::usage="Imports the SKA specifications from a file containing a row with the z-bin centers, a row with c-parameters
and a row specifying the volume of the survey with the pattern 'volume = some_number'. The optional parameter deltaz, specifies the width of each redshift bin, default: 0.1.
Form: importSKAspecs[fileTable_, deltaz_:0.1]"

zbinsSurvey::usage="Option for different specifications functions, that sets the redshift bins to be used. Default: $zbinGlobal"

biasInterpFuncFromGrowth::usage="compute the fiducial bias at the survey redshift bins center, by dividing the given biasfact through the Growth[z].
Options: $paramoptions, zbinsSurvey->$zbinGlobal.
Form: biasInterpFuncFromGrowth[biasfact_, options]"

readSurveySpecsFromFiles::usage="Reads from two text files the specifications for the
number density of galaxies and the fiducial bias as a function of redshift. All input given in terms of function options.
See: Options[readSurveySpecsFromFiles] for a list of options."

surveyName::usage="Option for readSurveySpecsFromFiles. It specifies the name of a survey so that the specification files are read correctly."


numberDensityFileFormat::usage="Option for readSurveySpecsFromFiles, specifies format of input file for ndens(z).
 File formats:
  CentralZinBin-dndOmegadz:  Consists of two columns, central redshift in bin and corresponding dndOmegadz. A line with the format volume = value or area = value of the survey is also accepted.
  CentralZinBin-FittingParameters: Consists of two rows, central redshift in bins and a row of fit parameters to construct ndens and bias. And a third line with the format volume = value or area = value of the survey."

biasFileFormat::usage="Option for readSurveySpecsFromFiles, specifies format of input file for b(z).
 File formats:
 CentralZinBin-biasFiducial: Consists of two columns, central redshift in bin and corresponding fiducial value of bias.
 BiasFactor-Growth: One line specifying bias_factor = value. The bias function is then bias[z]=bias_factor/Growth[z].
 Other: Other format, obtained from ndens file or from constants. "


setBiasFunctionFromFile::usage="Option for readSurveySpecsFromFiles, if set to True,
bias function  b(z) is read and set from files"
setNdensFunctionFromFile::usage="Option for readSurveySpecsFromFiles, if set to True,
number density function  n(z) is read and set from files"


setGalaxySurveySpecifications::usage="Option for readSurveySpecsFromFiles.
If set to True, readSurveySpecsFromFiles calls the function setGalaxySurveySpecs and passes to it the option: totalSkyAreaSurvey.
If all other specs are read correctly from files, then with this Option, there is no need to call setGalaxySurveySpecs separately."

specificationsDirectory::usage="Option for readSurveySpecsFromFiles:
Specify the directory in which the survey specifications are stored."
numberDensityZbinsFile::usage="Option for readSurveySpecsFromFiles:
Specify a file that contains specifications of number density as a function of z."
biasSpecsZbinsFile::usage="Option for readSurveySpecsFromFiles:
Specify a file that contains specifications of fiducial bias as a function of z."
dzBinWidthInFile::usage="Option for readSurveySpecsFromFiles: Specify widths of the redshift bins for the input file."


$ndensfile::usage="File name containing specifications for galaxy number densities as a function of z"
$biasfile::usage="File name containing specifications for fiducial bias as a function of z"
$specsDirectory::usage="Directory where specification files for surveys are stored."

$spectroSpecFiles::usage="Container for the specification filenames";
setGCspectroSurveySpecifications::usage="Function to set automatically specifications from files for specific Surveys."

specsSemiValidation::usage="Function that validates if an interpolating function depending on redshift returns a numeric value
when evaluated at a the center of a random redshift bin. Form: specsSemiValidation[func_InterpolatingFunction]"

Begin["`Private`"]

$biasfile="SurveyName-fiducial_bias-zbins.dat"  (*some default example of the filename of a txt file specifying the fiducial bias*)
$ndensfile="SurveyName-fiducial_dNdzddeg-zbins.dat"  (*some default example of the filename of a txt file specifying the fiducial dN/dz/ddeg number density per redshift bin per degree*)


biasOfzFuncSKAFit[zet_, cvect_] :=
    With[{c4 = cvect[[4]], c5 = cvect[[5]]}, c4 Exp[c5 zet]];

dndzFuncSKAFit[zet_, cvect_]:=With[{c1=cvect[[1]], c2=cvect[[2]], c3=cvect[[3]]},10^c1 * zet^c2*Exp[-c3*zet]];

ndensOfzFuncSKAFit[z_, cvect_ ] :=With[{degg=Pi/180},
  dndzFuncSKAFit[z,
    cvect]*(1/(degg^2 *((1 + z) distanceAngular[z])^2
      *(1/Hubble[z])  )
  )
] (*this will have the unis of 1/distanceAngular, see that function for more detail*)

setNumbDensBiasSKA[zbins_, cparamsvector_]:=Block[{nofzTable, bofzTable, zbinsav, interpfuncs, ndensinterp, biasinterp},
  zbinsav = zaverage[zbins];
  nofzTable=Table[{zz, ndensOfzFuncSKAFit[zz,cparamsvector]}, {zz, zbinsav}];
  $dndOmdz = nofzTable[[All, 2]];
  bofzTable=Table[{zz, biasOfzFuncSKAFit[zz,cparamsvector]}, {zz, zbinsav}];
  ndensinterp = Interpolation[nofzTable];
  biasinterp = Interpolation[bofzTable];
  interpfuncs = {ndensinterp, biasinterp};
  Return[interpfuncs]
];

importSKAspecs::wrongspecs="Wrong specifications file for SKA. It must contain a row with zbins centers, a row with c parameters and a row containing the word area or volume with an equal sign and a number"
importSKAspecs[fileTable_, deltaz_:0.1]:=Block[{cvect, zbins, volline, area, third, test},

  zbins = zAveragetoBins[fileTable[[1]], deltaz];
  cvect = fileTable[[2]];
  third = fileTable[[3]];
  volline = parselineWithNameEqualValue[third];
  test = StringMatchQ[volline[[1]], ("*vol*"|"*area*"), IgnoreCase->True];
  area = volline[[2]];
  Which[VectorQ[zbins,NumericQ]==False,
    debugPrint["zbins not numeric"];
    Message[importSKAspecs::wrongspecs],
    VectorQ[cvect,NumericQ]==False,
    debugPrint["cvector not numeric"];
    Message[importSKAspecs::wrongspecs],
    test==False,
    debugPrint["volume or area not specified in file"];
    Message[importSKAspecs::wrongspecs],
    NumberQ[area]==False,
    debugPrint["area value not numeric"];
    Message[importSKAspecs::wrongspecs]
    ];

    Return[{zbins, cvect, area}]
]

setNumbDens[ndenslist_List, zbinssurv_List]:=Block[{par, dndens3, zbinsSurvey, ndensfunct},
    dndens3=Flatten[ndenslist,1];
    (*zvaluesSurvey=zAveragetoBins[ndenstable[[All,1]],deltazbin];*)
    $dndOmdz=dndens3;
    zbinsSurvey=zbinssurv;
    ndensfunct = ndensVolumeInterpFunc[dndens3,zbinsSurvey];
    Return[ndensfunct]
  ]


Options[biasInterpFuncFromGrowth]=$paramoptions~Join~{zbinsSurvey->$zbinGlobal}

biasInterpFuncFromGrowth[biasfact_, opts:OptionsPattern[]]:=Block[{bintp, zbinscent=zaverage@OptionValue[zbinsSurvey], paropts=FilterRules[{opts}, $paramoptions]},
  bintp = Interpolation[Table[{zz, biasfact/Growth[zz,paropts]}, {zz,zbinscent}]];
  Return[bintp]
]

specsSemiValidation[func_InterpolatingFunction]:=Block[{zr,rb,zn, ll, tn},
  zr=zaverage@$zbinGlobal;
  ll = Length@zr;
  rb = RandomInteger[{1,ll}];
  Check[zn=func[zr[[rb]]], Return[False]];
  tn=NumericQ[zn];
  Return[tn]
]


zAveragetoBins[zlist_,dz_]:=Block[{dzs=dz/2,zbins},
zbins=N@DeleteDuplicates@Union[Rationalize[zlist+dzs],Rationalize[zlist-dzs]]
];

$specsDirectory=None;

Options[readSurveySpecsFromFiles]={numberDensityZbinsFile->$ndensfile, numberDensityFileFormat->"CentralZinBin-dndOmegadz",
  biasSpecsZbinsFile->$biasfile, biasFileFormat->"CentralZinBin-biasFiducial", dzBinWidthInFile->0.1, setBiasFunctionFromFile->True,
  setNdensFunctionFromFile->True, setGalaxySurveySpecifications->True,
  totalSkyAreaSurvey->$areaSurvey, specificationsDirectory->$specsDirectory, surveyName->None}


readSurveySpecsFromFiles::setbias=messageFormatText["Set bias function b(z) from input file succesfully.", FontColor->Green];
readSurveySpecsFromFiles::setndens=messageFormatText["Set number density function n(z) from input file succesfully.", FontColor->Green];


readSurveySpecsFromFiles::nosetbias=messageFormatText["Set bias function b(z) from input file was not succesful.", FontColor->Red];
readSurveySpecsFromFiles::nosetndens=messageFormatText["Set number density function n(z) from input file was not succesful.", FontColor->Red];

readSurveySpecsFromFiles::formatndens=messageFormatText["File formats:
  CentralZinBin-dndOmegadz:  Consists of two columns, central redshift in bin and corresponding dndOmegadz. A line with the format volume = value or area = value of the survey is also accepted.
  CentralZinBin-FittingParameters: Consists of two rows, central redshift in bins and a row of fit parameters to construct ndens and bias. And a third line with the format volume = value or area = value of the survey.",
  FontColor->Red];
readSurveySpecsFromFiles::formatbias=messageFormatText["File formats:
 CentralZinBin-biasFiducial: Consists of two columns, central redshift in bin and corresponding fiducial value of bias.
 BiasFactor-Growth: One line specifying bias_factor = value. The bias function is then bias[z]=bias_factor/Growth[z].
 FittingParameters: Bias computed with fitting parameters and fitting function.
 Other: Other format, obtained from ndens file or from constants.", FontColor->Red];




readSurveySpecsFromFiles::areapass=messageFormatText["setGalaxySurveySpecifications is set to True, but no compatible totalSkyAreaSurvey number was passed or read from files.
 Try again or set $areaSurvey", FontColor->Red]

readSurveySpecsFromFiles::filewarn="Warning: No valid string for the file name of the `1` specification file."

readSurveySpecsFromFiles::surveywarn="Warning: No valid name of the survey specified, if it is a different survey and the format of files coincides with the existing ones, it will work.
See numberDensityFileFormat biasFileFormat usage messages for more information."

readSurveySpecsFromFiles[opts:OptionsPattern[]]:=Block[{nfileform,bfileform,dzbinfile, ndenstable, biastable, dndens3,
  zvaluesSurvey, zvaluesBias, biasfidu, biasZs, biasfunc, biaslist, ndensfunc, ndenslist,
  surveynam=OptionValue[surveyName], vol, dndztab, biastab, cvecto, biasfact, biform, ndform},

  If[SameQ[OptionValue[specificationsDirectory], None],
    $specsDirectory=NotebookDirectory[]<>"/";
    ,
    $specsDirectory=OptionValue[specificationsDirectory];
  ];

  dzbinfile=OptionValue[dzBinWidthInFile];

  (*If[SameQ[surveynam, None],
    Message[readSurveySpecsFromFiles::surveywarn]
    ,
    {ndform, biform} = matchSurveyNamesWithFileFormats[surveynam];
    SetOptions[readSurveySpecsFromFiles, biasFileFormat->biform];
    SetOptions[readSurveySpecsFromFiles, numberDensityFileFormat->ndform];
    SetOptions[readSurveySpecsFromFiles, surveyName->surveynam];
  ];*)

  SetOptions[readSurveySpecsFromFiles,specificationsDirectory->$specsDirectory]; (*option set here, to avoid FrontEnd Problem*)
  (*Print[$specsDirectory];*)
  If[ StringQ[OptionValue[numberDensityZbinsFile]] == True,
    $ndensfile=$specsDirectory<>OptionValue[numberDensityZbinsFile];
    existsFile[$ndensfile],
    Message[readSurveySpecsFromFiles::filewarn, "ndens"];
    Interrupt[]
  ];
  If[ StringQ[OptionValue[biasSpecsZbinsFile]] == True,
    $biasfile=$specsDirectory<>OptionValue[biasSpecsZbinsFile];,
    Message[readSurveySpecsFromFiles::filewarn, "bias"];
    Interrupt[]
  ];



  If[OptionValue[setGalaxySurveySpecifications]==True && OptionValue[setNdensFunctionFromFile]==True,
    Which[
      OptionValue[numberDensityFileFormat]=="CentralZinBin-FittingParameters",
        existsFile[$ndensfile];
        ndenstable=Import[$ndensfile,"Table"];
        {zvaluesSurvey, cvecto, vol} = importSKAspecs[ndenstable,dzbinfile];
        {ndensfunc, biasfunc} = setNumbDensBiasSKA[zvaluesSurvey, cvecto];
        $areaSurvey = vol;
      ,
      OptionValue[numberDensityFileFormat]=="CentralZinBin-dndOmegadz",
        existsFile[$ndensfile];
        ndenstable = importNumericTxtTable[$ndensfile];
        ndenslist = ndenstable[[All,2]];
        dzbinfile=Abs[ ndenstable[[2,1]] - ndenstable[[1,1]] ];
        zvaluesSurvey=zAveragetoBins[ndenstable[[All,1]], dzbinfile]; (*FIX THIS*)
        ndensfunc = setNumbDens[ndenslist, zvaluesSurvey];
        vol = importValueOfParameterFromTxtFile[$ndensfile, ("*vol*" | "*area*")];
        $areaSurvey = vol;
        ,
      OptionValue[numberDensityFileFormat]=="LimitsZBin-dndz-bias",
      existsFile[$ndensfile];
      ndenstable = importNumericTxtTable[$ndensfile];
      ndenslist = ndenstable[[All,4]];
      zvaluesSurvey=Union[ndenstable[[All,1]], ndenstable[[All,3]] ];
      ndensfunc = setNumbDens[ndenslist, zvaluesSurvey];
      vol = importValueOfParameterFromTxtFile[$ndensfile, ("*vol*" | "*area*")];
      $areaSurvey = vol;
      ,
      True,
      Message[readSurveySpecsFromFiles::formatndens];
      Message[readSurveySpecsFromFiles::surveywarn];
    ];

    $zbinGlobal=zvaluesSurvey;
    $dzBinWidth=dzbinfile;
    $ndensGalZFuncGlobal = ndensfunc;
    If[specsSemiValidation[$ndensGalZFuncGlobal] === True,
      Message[readSurveySpecsFromFiles::setndens];
      ,
      Message[readSurveySpecsFromFiles::nosetndens];
      Interrupt[];
    ];
    SetOptions[setGalaxySurveySpecs, surveyRedshiftBins->$zbinGlobal, ndensFunction->$ndensGalZFuncGlobal,
        dzBinWidth->$dzBinWidth, dnumbdOmegadz->$dndOmdz, customNdensFunction->None];
  ];


  If[OptionValue[setGalaxySurveySpecifications]==True && OptionValue[setBiasFunctionFromFile]==True,
    Which[
      OptionValue[biasFileFormat]=="CentralZinBin-biasFiducial",
      existsFile[$biasfile];
      biastable=importNumericTxtTable[$biasfile];
      biasfidu=biastable[[All,2]];
      biasZs = biastable[[All,1]];
      biasfunc = Interpolation[Transpose[{biasZs,biasfidu}]];
      ,
      OptionValue[biasFileFormat]=="LimitsZBin-dndz-bias",
      existsFile[$biasfile];
      biastable=importNumericTxtTable[$biasfile];
      biasfidu=biastable[[All,5]];
      biasZs = biastable[[All,2]];
      biasfunc = Interpolation[Transpose[{biasZs,biasfidu}]];
      ,
      OptionValue[biasFileFormat]=="FittingParameters",
      debugPrint["bias was set together with ndens file"]
      ,
      OptionValue[biasFileFormat]=="BiasFactor-Growth",
      existsFile[$biasfile];
      biasfact = importValueOfParameterFromTxtFile[$biasfile, ("*bias*")];
      biasfunc = biasInterpFuncFromGrowth[biasfact];
      ,
      OptionValue[biasFileFormat]=="Other",
      biaslist = {#, 1.0}&/@(zaverage@$zbinGlobal);  (*bias is just a constant 1, as a test*)
      biasfunc = Interpolation[biaslist];
      ,
      True,
      Message[readSurveySpecsFromFiles::formatbias];
      Interrupt[]
    ];
    $biasInterpFunc = biasfunc;
    If[specsSemiValidation[$biasInterpFunc] === True,
      Message[readSurveySpecsFromFiles::setbias];,
      Message[readSurveySpecsFromFiles::nosetbias];
      Interrupt[];
    ];
    SetOptions[setGalaxySurveySpecs, biasFunction->$biasInterpFunc, customBiasFunction->None];
  ];


  If[OptionValue[setGalaxySurveySpecifications]==False,
    Return[{ndensfunc,biasfunc}]
  ];

  If[OptionValue[setGalaxySurveySpecifications]==True,
    If[NumberQ[$areaSurvey]==False,
      Message[readSurveySpecsFromFiles::areapass],
      SetOptions[readSurveySpecsFromFiles, totalSkyAreaSurvey->$areaSurvey];
      setGalaxySurveySpecs[totalSkyAreaSurvey->$areaSurvey];
    ];
  ]

];

specsFiles["Euclid-ndens"]="Euclid_3-zmin-zm-zmax-dN_dzddeg-fiducial_bias.dat";
specsFiles["Euclid-bias"]="Euclid_3-zmin-zm-zmax-dN_dzddeg-fiducial_bias.dat";
specsFiles["DESI-ndens"]="DESI-ELG-dN_dzddeg-zbins-deltaz_0p1.txt";
specsFiles["DESI-bias"]="DESI-ELG-fiducial_bias_factor.txt";
specsFiles["SKA2-ndens"]="SKA2-zbincenter-cparams_fit.txt";
specsFiles["SKA1-ndens"]="SKA1-zbincenter-cparams_fit.txt";
specsFiles["SKA2-bias"]="SKA2-zbincenter-cparams_fit.txt";
specsFiles["SKA1-bias"]="SKA1-zbincenter-cparams_fit.txt";
$spectroSpecFiles=specsFiles;


setGCspectroSurveySpecifications[chosenSurvey_]:=Block[{nzstr, bzstr, nzfile, bfile, nzformat, bzformat},
  nzstr = chosenSurvey<>"-ndens";
  bzstr = chosenSurvey<>"-bias";
  nzfile = specsFiles[nzstr];
  bzfile = specsFiles[bzstr];

If[FileExistsQ[$specsDirectory<>nzfile],
Print["specifications for number density n(z): file exists"]];
If[FileExistsQ[$specsDirectory<>bzfile],
Print["specifications for galaxy bias b(z): file exists"]];

  Which[chosenSurvey=="Euclid",
  nzformat="LimitsZBin-dndz-bias";
  bzformat="LimitsZBin-dndz-bias";
,
chosenSurvey=="DESI",
  nzformat="CentralZinBin-dndOmegadz";
  bzformat="BiasFactor-Growth";
,
chosenSurvey=="SKA2",
  nzformat="CentralZinBin-FittingParameters";
  bzformat="FittingParameters"
,
chosenSurvey=="SKA1",
  nzformat="CentralZinBin-FittingParameters";
  bzformat="FittingParameters"
  ];

  readSurveySpecsFromFiles[numberDensityZbinsFile->nzfile,
biasSpecsZbinsFile->bzfile,setGalaxySurveySpecifications->True,
specificationsDirectory->$specsDirectory,surveyName->chosenSurvey,
numberDensityFileFormat->nzformat,
biasFileFormat->bzformat];

Print[styleFunction[chosenSurvey<>" GC spectroscopic Specifications loaded ****",25]];
]


matchSurveyNamesWithFileFormats[survnam_String]:=Block[{nam, bform, ndform},
  Which[
    StringMatchQ[surveynam, "*SKA*"],
    ndform="CentralZinBin-FittingParameters";
    bform="FittingParameters";
    ,
    StringMatchQ[surveynam, "*Euclid*"],
    ndform = "CentralZinBin-dndOmegadz";
    bform = "CentralZinBin-biasFiducial";
    ,
    StringMatchQ[surveynam, "*DESI*"],
    ndform = "CentralZinBin-dndOmegadz";
    bform = "BiasFactor-Growth";
  ];
  Return[{ndform, bform}]
]


End[] (* `Private` *)

EndPackage[]
