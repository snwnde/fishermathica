(* ::Package:: *)

BeginPackage["GenerateInput`"]

Unprotect["GenerateInput`*", "GenerateInput`*`*"]  (*Unprotect and Protect definitions for using ClearAll in the case of developing code*)

ClearAll["GenerateInput`*", "GenerateInput`*`*"]




haloversion::usage="An option for CAMB. Integer number that chooses the Halofit version used in the non-linear matter power spectrum.
1: Original, 2: Bird, 3: Peacock, 4: Takahashi, 5: Mead, 6: Halo model, 7: Casarini, 9:Winther";

$haloversionRule::usage="Correspondance between haloversion number and name of the Halofit prescription."

$parameterInfoRules::usage="Rule containing all the information of the parameters sent to CAMB."
computeEpsilonParameterFiles::usage="Function: computeEpsilonParameterFiles[listofzrange_,epsi_,opts].
Function that computes and outputs input files in the common input format, for a specified listofzrange,
containing a list of output redshifts and for a specific epsilon value varying the fiducial
parameters in $paramoptions."
outputSpectraFormat::usage="Option for computeEpsilonParameterFiles.
1: power spectra as a table of many z and k values, standard for IST:F in WL and XC.
2: one file per z, with {k, P(k), sigma8}, standard for IST:F GCsp.
3: Common Input for all probes, matrices for P(z,k). New standard for all future forecasts."

computeEpsilons::usage="If set to True, computes variation of the power spectra with respect to the cosmological parameters."
computeFiducial::usage="If set to True, computes the fiducial power spectrum."
referenceAmplitude::usage="Reference primordial amplitude for the power spectrum. Usually $Asreference."
referenceSigma8::usage="Reference sigma8 (at z=0) amplitude for the power spectrum. Usually $sigma8reference."

exportDirectory::usage="Directory where all files are exported to. Default: './' "
extraNumericCAMBparams::usage="List of precision and numerical parameters to be passed to CAMB."
exportDetailedInfo::usage="If set to True, prints a file with detailed information on the cosmological and numerical parameters, file names and time stamps."
exportPowerSpectra::usage="If set to True, exports the power spectra files. Default: True. Set to False to produce only information files."
constantNeutrinoParams::usage="List of parameters for the case of constant neutrino properties passed to CAMB."
parameterIndices::usage="List of integers, corresponding to the positions of parameters in $paramoptions for which to generate the power spectra variations."
fiducialFolderName::usage="Name of the folder containing the fiducial values.
Default: `fiducial_eps_0`"

Needs["cosmomathica`interface`"];
Needs["CosmologyFunctions`"];
Needs["UsefulTools`"];
EndPackage[]


BeginPackage["GenerateInput`", {"CosmologyFunctions`","UsefulTools`", "cosmomathica`interface`"}]




Begin["`Private`"]


$haloversionRule={1->"Original", 2->"Bird", 3->"Peacock", 4->"Takahashi", 5->"Mead",
6->"HaloModel", 7->"Casarini", 9->"Winther"}

cambbasis[optip_, Asvar_]:=({OmegaBaryon0Today[optip]*(hubbleToday[optip]^2),(OmegaCDM0Today[optip])*(hubbleToday[optip]^2),
(hubbleToday[optip]),10^9*Asvar});
printcamb[cblist_]:=({"Corresponding physical base parameters: "<>
"{omegab, "<>" omegac, "<>"h, "<>"10^9 As} :> "
<>ToString[{FortranForm@cblist[[1]],FortranForm@cblist[[2]],FortranForm@cblist[[3]],FortranForm@cblist[[4]]}]});

standbasis[optip_, Asvar_]:=({OmegaBaryon0Today[optip],(OmegaCDM0Today[optip]),
(hubbleToday[optip]),10^9*Asvar});
printstand[cblist_]:=({"Corresponding non-physical base parameters: "<>
"{Omegab, "<>" Omegac, "<>"h, "<>"10^9 As} :> "
<>ToString[{FortranForm@cblist[[1]],FortranForm@cblist[[2]],FortranForm@cblist[[3]],FortranForm@cblist[[4]]}]});



ClearAll[computeEpsilonParameterFiles];

Options[computeEpsilonParameterFiles]={computeEpsilons->True,
  computeFiducial->True, referenceAmplitude:>$Asreference, TransferKmax->50,
HalofitVersion->4,
exportDirectory->"./", extraNumericCAMBparams->False,
exportDetailedInfo->True, exportPowerSpectra->True, constantNeutrinoParams->False,
parameterIndices:>Range[Length[$paramoptions]], outputSpectraFormat->1,
referenceSigma8:>$sigma8reference, outputFileExtension->".txt",
fiducialFolderName->"fiducial_eps_0"};

(* ::Code::Bold:: *)
(**)


computeEpsilonParameterFiles[listofzrange_,epsi_,opts:OptionsPattern[]]:=Block[{compueps,
  compufid,epsilonstep,halovers,Tkmax, numcambopt, neutrinopts,
  dirnames,
  optplus, optminus, step, Asref, AsrescaledList,
  Asoptplus, Asoptminus,
  lenlistofz,
  expodir, plusmin, detailedInfo={},
  expoinfo,Omeganuu, paramindexlist,
  outext,
  outZtab,outKtab,outZtabPlus,outKtabPlus,outZtabMinus,outKtabMinus,
  growthTab,fgrowthrateTab,growthTabPlus,fgrowthrateTabPlus,growthTabMinus,fgrowthrateTabMinus,
  HofZtab,s8ofZtab,HofZtabPlus,s8ofZtabPlus,HofZtabMinus,s8ofZtabMinus,
  psLinTab,psNonLinTab, psLinTabPlus, psNonLinTabPlus,
  psLinTabMinus, psNonLinTabMinus,
  fiduDir, inidate, findate, calctime
},
inidate=DateString[];
fiduDir=OptionValue[fiducialFolderName];
epsilonstep=epsi;
outext=OptionValue[outputFileExtension];
halovers=OptionValue[HalofitVersion];
Tkmax=OptionValue[TransferKmax];
lenlistofz=Length[listofzrange];
compufid=OptionValue[computeFiducial];
compueps=OptionValue[computeEpsilons];
Asref=OptionValue[referenceAmplitude];
s8ref=OptionValue[referenceSigma8];
$parameterInfoRules={"Asref"->Asref,
"s8ref"->s8ref,
 "haloversion"->halovers,
 "Tkmax"->Tkmax
};
AppendTo[$parameterInfoRules,
"ListOfzRange"->listofzrange];
SetOptions[LCDMCAMBPsPre, {TransferKmax->Tkmax, HalofitVersion->halovers}];
AsrescaledList={};
expodir=OptionValue[exportDirectory];
expoinfo=OptionValue[exportDetailedInfo];
mkDirectory[expodir];
plusmin={"pl","mn"};
numcambopt=OptionValue[extraNumericCAMBparams];
If[UnsameQ[numcambopt, False],
SetOptions[LCDMCAMBPsPre, numcambopt];
];
AppendTo[$parameterInfoRules, "NumericalCambParameters"->numcambopt];
neutrinopts=OptionValue[constantNeutrinoParams];
If[UnsameQ[neutrinopts, False],
SetOptions[LCDMCAMBPsPre, neutrinopts];
Omeganuu=(OmegaNu/.neutrinopts);
];
AppendTo[$parameterInfoRules, "NeutrinoCambParameters"->neutrinopts];
$Asreference=Asref;
paramindexlist=OptionValue[parameterIndices];

If[compufid==True,
  {outZtab,outKtab,HofZtab,
    s8ofZtab,growthTab,fgrowthrateTab,
    psLinTab, psNonLinTab
    }=LCDMCAMBPsPre[listofzrange, checkConsistency->False];
  mkDirectory[expodir<>fiduDir<>"/"];
  Export[expodir<>fiduDir<>"/"<>"z_values_list"<>outext, outZtab, "Table"];
  Export[expodir<>fiduDir<>"/"<>"k_values_list"<>outext, outKtab, "Table"];
  Export[expodir<>fiduDir<>"/"<>"background_Hz"<>outext, HofZtab, "Table"];
  Export[expodir<>fiduDir<>"/"<>"sigma8-z"<>outext, s8ofZtab, "Table"];
  Export[expodir<>fiduDir<>"/"<>"D_Growth-zk"<>outext, growthTab, "Table"];
  Export[expodir<>fiduDir<>"/"<>"f_GrowthRate-zk"<>outext, fgrowthrateTab, "Table"];
  Export[expodir<>fiduDir<>"/"<>"Plin-zk"<>outext, psLinTab, "Table"];
  Export[expodir<>fiduDir<>"/"<>"Pnonlin-zk"<>outext, psNonLinTab, "Table"];
];

If[expoinfo==True,
detailedInfo=printDetailedInfo["fiducial", detailedInfo,
expodir, $parameterInfoRules];
];

If[compueps==False, Print["No epsilon files computed."]; Return[False] ];

Do[
   Print["Current parameter: "];
   Print[$paramoptions[[indi]][[1]] ];
   Print["Relative epsilon value : "<>ToString[epsilonstep]];
   $debugPrint=False;

   {optplus,optminus,step}=numericalDerivativeStep[indi,epsilonstep];
   If[UnsameQ[$paramoptions[[indi]][[1]] , sigma8],
     Asoptplus=LCDMCAMBPsPre[0.,  {optplus}~Join~{checkConsistency->True, returnRescaledAs->True}];
     Asoptminus=LCDMCAMBPsPre[0., {optminus}~Join~{checkConsistency->True, returnRescaledAs->True}];
     AppendTo[$parameterInfoRules, "epsilonstep"->epsilonstep];
     AppendTo[$parameterInfoRules, "optplus"->optplus];
     AppendTo[$parameterInfoRules, "Asoptplus"->Asoptplus];
     AppendTo[$parameterInfoRules, "optminus"->optminus];
     AppendTo[$parameterInfoRules, "Asoptminus"->Asoptminus];
     If[expoinfo==True,
        detailedInfo=printDetailedInfo["epsilon", detailedInfo,
        "False", $parameterInfoRules];
      ];
     ,(*Else, SameQ[$paramoptions[[indi]][[1]] , sigma8] *)
     Asoptplus=LCDMCAMBPsPre[0., {optplus}~Join~{checkConsistency->False, returnRescaledAs->True}];
     Asoptminus=LCDMCAMBPsPre[0., {optminus}~Join~{checkConsistency->False, returnRescaledAs->True}];
     AppendTo[$parameterInfoRules, "epsilonstep"->epsilonstep];
     AppendTo[$parameterInfoRules, "optplus"->optplus];
     AppendTo[$parameterInfoRules, "Asoptplus"->Asoptplus];
     AppendTo[$parameterInfoRules, "optminus"->optminus];
     AppendTo[$parameterInfoRules, "Asoptminus"->Asoptminus];
     If[expoinfo==True,
     detailedInfo=printDetailedInfo["epsilon", detailedInfo,
     expodir, $parameterInfoRules];
     ];
   ];

  If[UnsameQ[$paramoptions[[indi]][[1]], sigma8] , $Asreference=Asoptplus; , $Asreference=Asref;];
  {outZtabPlus,outKtabPlus,HofZtabPlus,
    s8ofZtabPlus,growthTabPlus,fgrowthrateTabPlus,
    psLinTabPlus, psNonLinTabPlus}=LCDMCAMBPsPre[listofzrange, {optplus}~Join~{checkConsistency->False}];

  If[UnsameQ[$paramoptions[[indi]][[1]], sigma8], $Asreference=Asoptminus; , $Asreference=Asref;];
  {outZtabMinus,outKtabMinus,HofZtabMinus,
    s8ofZtabMinus,growthTabMinus,fgrowthrateTabMinus,
    psLinTabMinus, psNonLinTabMinus}=LCDMCAMBPsPre[listofzrange, {optminus}~Join~{checkConsistency->False}];

  dirnames=Flatten[((Table[(#<>"_"<>plm<>"_"<>scientificfortranform[epsilonstep, "E", "p", "eps_"]), {plm,plusmin}])&/@$paramlabels)];
  expodirplus=expodir<>dirnames[[2*indi-1]];
  expodirminus=expodir<>dirnames[[2*indi]];
  mkDirectory[expodirplus];
  mkDirectory[expodirminus];

  Print["Exporting plus epsilon files..."];
  Export[expodirplus<>"/"<>"z_values_list"<>outext,          outZtabPlus, "Table"];
  Export[expodirplus<>"/"<>"k_values_list"<>outext,          outKtabPlus, "Table"];
  Export[expodirplus<>"/"<>"background_Hz"<>outext,          HofZtabPlus, "Table"];
  Export[expodirplus<>"/"<>"sigma8-z"<>outext,              s8ofZtabPlus, "Table"];
  Export[expodirplus<>"/"<>"D_Growth-zk"<>outext,          growthTabPlus, "Table"];
  Export[expodirplus<>"/"<>"f_GrowthRate-zk"<>outext, fgrowthrateTabPlus, "Table"];
  Export[expodirplus<>"/"<>"Plin-zk"<>outext,               psLinTabPlus, "Table"];
  Export[expodirplus<>"/"<>"Pnonlin-zk"<>outext,         psNonLinTabPlus, "Table"];

  Print["Exporting minus epsilon files..."];
  Export[expodirminus<>"/"<>"z_values_list"<>outext,          outZtabMinus, "Table"];
  Export[expodirminus<>"/"<>"k_values_list"<>outext,          outKtabMinus, "Table"];
  Export[expodirminus<>"/"<>"background_Hz"<>outext,          HofZtabMinus, "Table"];
  Export[expodirminus<>"/"<>"sigma8-z"<>outext,              s8ofZtabMinus, "Table"];
  Export[expodirminus<>"/"<>"D_Growth-zk"<>outext,          growthTabMinus, "Table"];
  Export[expodirminus<>"/"<>"f_GrowthRate-zk"<>outext, fgrowthrateTabMinus, "Table"];
  Export[expodirminus<>"/"<>"Plin-zk"<>outext,               psLinTabMinus, "Table"];
  Export[expodirminus<>"/"<>"Pnonlin-zk"<>outext,         psNonLinTabMinus, "Table"];
  ,
  {indi,paramindexlist}
  ];
findate=DateString[];
calctime = timeElapsed[inidate,findate];
Return[calctime]
];


printDetailedInfo[infoType_, detailedInfo_, expodir_:"False",
parameterInfoRules_]:=Block[{deta=detailedInfo,
  pAsref="Asref"/.parameterInfoRules,
  ps8ref="s8ref"/.parameterInfoRules,
  phalovers="haloversion"/.parameterInfoRules,
  pTkmax="Tkmax"/.parameterInfoRules,
  pnumcambopt="NumericalCambParameters"/.parameterInfoRules,
  pneutrinopts="NeutrinoCambParameters"/.parameterInfoRules,
  plistofzrange="ListOfzRange"/.parameterInfoRules,
  expofile,
  pepsilonstep,
  poptplus,
  poptminus,
  pAsoptplus,
  pAsoptminus},
AppendTo[deta, {"Initial date: "<>ToString[DateString[]]}];
AppendTo[deta, {"Fiducial parameters: "<>ToString[$paramoptions]}];
AppendTo[deta, {"Reference Amplitude As: "<>ToString[FortranForm@pAsref]}];
AppendTo[deta, {"Reference sigma8: "<>ToString[ps8ref]}];
AppendTo[deta, printstand[standbasis[{}, pAsref]]];
AppendTo[deta, printcamb[cambbasis[{}, pAsref]]];
AppendTo[deta, {"Halofit version in CAMB: "<>ToString[phalovers]}];
AppendTo[deta, {"Maximum k value in CAMB (h/Mpc): "<>ToString[pTkmax]}];
AppendTo[deta, {"Numeric CAMB parameters : "<>ToString[pnumcambopt]}];
AppendTo[deta, {"Constant Neutrino CAMB parameters : "<>ToString[pneutrinopts]}];
AppendTo[deta, {"List of z points: "<>ToString[plistofzrange]}];
If[
infoType=="fiducial",
infoname="Fiducial-detailedInfo";
expofile=expodir<>infoname<>".txt";
If[expodir != "False",
Print["Exporting file: "<>expofile];
Export[expofile, deta, "Table"];
];
];
If[
infoType=="epsilon",
pAsoptminus="Asoptminus"/.parameterInfoRules;
pAsoptplus="Asoptplus"/.parameterInfoRules;
poptplus="optplus"/.parameterInfoRules;
poptminus="optminus"/.parameterInfoRules;
pepsilonstep="epsilonstep"/.parameterInfoRules;
(** Adding information to detailedInfo table**)

   AppendTo[deta, {"***Parameter Variations:  ***"}];
   AppendTo[deta, {"Relative epsilon value : "<>ToString[pepsilonstep]}];
   AppendTo[deta, {"Fiducial reference sigma8 : "<>ToString[(ps8ref)]}];
   AppendTo[deta, {"Fiducial reference As : "<>ToString[(pAsref)]}];
   AppendTo[deta, {"++Parameter plus change: "<>ToString[poptplus]}];
   AppendTo[deta, {"Corresponding As value: "<>ToString[FortranForm[pAsoptplus]]}];
   AppendTo[deta, printcamb[cambbasis[poptplus, pAsoptplus]]];
   AppendTo[deta, printstand[standbasis[poptplus, pAsoptplus]]];
   AppendTo[deta, {"--Parameter minus change: "<>ToString[poptminus]}];
   AppendTo[deta, {"Corresponding As value: "<>ToString[FortranForm[pAsoptminus]]}];
   AppendTo[deta, printcamb[cambbasis[poptminus, pAsoptminus]]];
   AppendTo[deta, printstand[standbasis[poptminus, pAsoptminus]]];

If[expodir != "False",
infoname=scientificfortranform[pepsilonstep, "E", "p", "eps_"]<>"-detailedInfo";
expofile=expodir<>infoname<>".txt";
Print["Exporting file: "<>expofile];
Export[expofile, deta, "Table"];
];
];
Return[deta]
];

End[]
EndPackage[]

(*
If[expoinfo==True,
AppendTo[detailedInfo, {"***Parameter Variations:  ***"}];
AppendTo[detailedInfo, {"++Parameter plus change: "<>ToString[optplus]}];
AppendTo[detailedInfo, {"Corresponding sigma8(z="<>ToString[listofzrange[[zzi]]]<>"): "<>ToString[sigma8valofzplus]}];
AppendTo[detailedInfo, {"Corresponding sigma8(z=0): "<>ToString[sigma8val0plus]}];
AppendTo[detailedInfo, {"--Parameter minus change: "<>ToString[optminus]}];
AppendTo[detailedInfo, {"Corresponding sigma8(z="<>ToString[listofzrange[[zzi]]]<>"): "<>ToString[sigma8valofzminus]}];
AppendTo[detailedInfo, {"Corresponding sigma8(z=0): "<>ToString[sigma8val0minus]}];
Print["Exporting info file..."];
AppendTo[detailedInfo, {"Parameter variation generation finished. Date: "<>ToString[DateString[]]}];
Export[expodir<>scientificfortranform[epsilonstep, "E", "p", "eps_"]<>"-detailedInfo.txt", detailedInfo, "Table"]
];
*)
