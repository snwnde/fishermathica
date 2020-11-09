(* ::Package:: *)

ClearAll[computeEpsilonParameterFiles]


Options[computeEpsilonParameterFiles]={computeEpsilons->True, computeFiducial->True, referenceAmplitude:>$Asreference, TransferKmax->50,
HalofitVersion->4,
exportDirectory->"./", extraNumericCAMBparams->False, exportDetailedInfo->True, exportPowerSpectra->True, constantNeutrinoParams->False,
parameterIndices:>Range[Length[$paramoptions]], outputSpectraFormat->1, referenceSigma8:>$sigma8reference, outputFileExtension->".txt"};


haloversion::usage="An option for CAMB. Integer number that chooses the Halofit version used in the non-linear matter power spectrum.
1: Original, 2: Bird, 3: Peacock, 4: Takahashi, 5: Mead, 6: Halo model, 7: Casarini, 9:Winther";


haloversionRule={1->"Original", 2->"Bird", 3->"Peacock", 4->"Takahashi", 5->"Mead", 6->"HaloModel", 7->"Casarini", 9->"Winther"}


outputSpectraFormat::usage="Option for computeEpsilonParameterFiles. 1: power spectra as a table of many z and k values, useful for WL.
2: one file per z, with {k, P(k), sigma8}, useful for GC."

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


(* ::Code::Bold:: *)
(**)


computeEpsilonParameterFiles[listofzrange_,epsi_,opts:OptionsPattern[]]:=Block[{compueps,compufid,epsilonstep,halovers,Tkmax, numcambopt, neutrinopts,
zpkTableFidu, zpkTablePlus, zpkTableMinus, dirnames,
optplus, optminus, step, Asref, AsrescaledList,Asoptplus, Asoptminus,pklistfidu,pklistplus,pklistminus, lenlistofz,
expodir,expospectra, plusmin, detailedInfo={},
expoinfo, wmegac, wmegab, wbwcfunc,wopt, cambbasis, printcamb,printstand,standbasis, Omeganuu, paramindexlist, outfmt,
  lastindx, pkintpfidu,log10kvalsintpfidu,sigma8valofz,s8list, sigma8val0, kvalsintpfidu,pklistfiduovers8sq, pkintp,
pkintpplus,log10kvalsintpplus,sigma8valofzplus, sigma8val0plus, pklistplusovers8sq,s8listplus,
pkintpminus,log10kvalsintpminus,sigma8valofzminus, sigma8val0minus, pklistminusovers8sq,s8listminus, s8ref, outext,
  outZtab,outKtab,outZtabPlus,outKtabPlus,outZtabMinus,outKtabMinus,
  growthTab,fgrowthrateTab,growthTabPlus,fgrowthrateTabPlus,
  growthTabMinus,fgrowthrateTabMinus,
  allpertsTab, allpkTab
},
epsilonstep=epsi;
outfmt=OptionValue[outputSpectraFormat];
outext=OptionValue[outputFileExtension];
halovers=OptionValue[HalofitVersion];
Tkmax=OptionValue[TransferKmax];
lenlistofz=Length[listofzrange];
compufid=OptionValue[computeFiducial];
compueps=OptionValue[computeEpsilons];
Asref=OptionValue[referenceAmplitude];
s8ref=OptionValue[referenceSigma8];
AsrescaledList={};
expodir=OptionValue[exportDirectory];
expospectra=OptionValue[exportPowerSpectra];
expoinfo=OptionValue[exportDetailedInfo];
mkDirectory[expodir];
plusmin={"pl","mn"};
numcambopt=OptionValue[extraNumericCAMBparams];
If[UnsameQ[numcambopt, False],
SetOptions[LCDMCAMBPsPre, numcambopt];
];
neutrinopts=OptionValue[constantNeutrinoParams];
If[UnsameQ[neutrinopts, False],
SetOptions[LCDMCAMBPsPre, neutrinopts];
Omeganuu=(OmegaNu/.neutrinopts);
];
$Asreference=Asref;
paramindexlist=OptionValue[parameterIndices];
(*wmegab=OmegaBaryon0Today[]*(hubbleToday[]^2);
wmegac=OmegaCDM0Today[]*(hubbleToday[]^2);*)
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

If[compufid==True,
Which[
outfmt==1,
  Do[
  If[FileExistsQ[expodir<>"pkz-Fiducial"<>outext], Print["Fiducial file already exists"]; Break[]];
  $debugPrint=False;
  Print["z values: ", Reverse[listofzrange[[zzi]]] ];
  If[zzi==Reverse[Range[lenlistofz]][[1]], zpkTableFidu={}; $debugPrint=True; Print["Printing debug values for LCDMCAMBPsPre"]];
  pklistfidu=LCDMCAMBPsPre[listofzrange[[zzi]], returnInterpolated->False ,checkConsistency->True, TransferKmax->Tkmax, HalofitVersion->halovers];
  Asref=$Asreference;  (* double check in case passed Asref is not the correct one. The checkConsistency function inside LCDMCAMBPsPre,
  overwrites the $Asreference with the correct one. This is then passed to Asref.*)
  pklistfidu=Flatten[pklistfidu,1];
  zpkTableFidu=zpkTableFidu~Join~pklistfidu;
  If[zzi==Reverse[Range[lenlistofz]][[-1]],
  Export[expodir<>"pkz-Fiducial"<>outext, zpkTableFidu, "Table"];
  ];
  ,
  {zzi,Reverse[Range[lenlistofz]]}
  ];
,
outfmt==2,
  {outZtab,outKtab,growthTab,fgrowthrateTab}=LCDMCAMBPsPre[listofzrange, returnGrowthRate->True, returnInterpolated->False, TransferKmax->Tkmax, HalofitVersion->halovers, checkConsistency->False];
  mkDirectory[expodir<>"fiducial/"];
  Export[expodir<>"z_values_list"<>outext, outZtab, "Table"];
  Export[expodir<>"k_values_list"<>outext, outKtab, "Table"];
  Export[expodir<>"fiducial/"<>"D_Growth_zk-fiducial"<>outext, growthTab, "Table"];
  Export[expodir<>"fiducial/"<>"f_GrowthRate_zk-fiducial"<>outext, fgrowthrateTab, "Table"];
  Do[
  If[FileExistsQ[expodir<>"fiducial/"<>"Pk-fiducial-z_"<>IntegerString[zzi-1,10,3]<>outext], Print["Fiducial file already exists"]; Break[]];
  If[zzi==1,$debugPrint=True,$debugPrint=False];
  Print["z value: ", listofzrange[[zzi]] ];
  {pkintpfidu,log10kvalsintpfidu,sigma8valofz,sigma8val0}=LCDMCAMBPsPre[listofzrange[[zzi]],
  returnInterpolated->True, returnSigma8->True, returnGrowthRate->False, TransferKmax->Tkmax, HalofitVersion->halovers, checkConsistency->False];
  Export[expodir<>"log_k_sampling"<>outext, log10kvalsintpfidu, "Table"];
  kvalsintpfidu=10^(log10kvalsintpfidu);
  pklistfidu=10^(Table[pkintpfidu[kkl],{kkl,log10kvalsintpfidu}]);
  pklistfiduovers8sq=pklistfidu/(sigma8valofz^2);
  s8list=Table[sigma8valofz,{kkl,1,Length@log10kvalsintpfidu}];
  zpkTableFidu=Transpose[{kvalsintpfidu,pklistfiduovers8sq,s8list}];
  PrependTo[zpkTableFidu, {"#  k [h/Mpc]          Pk/s8^2 [Mpc/h]^3          s8"}];
  mkDirectory[expodir<>"fiducial/"];
  Export[expodir<>"fiducial/"<>"Pk-fiducial-z_"<>IntegerString[zzi-1,10,3]<>outext, zpkTableFidu, "Table"];
    ,
  {zzi,1,lenlistofz}
  ];
];
];

If[expoinfo==True,
AppendTo[detailedInfo, {"Initial date: "<>ToString[DateString[]]}];
AppendTo[detailedInfo, {"Fiducial parameters: "<>ToString[$paramoptions]}];
AppendTo[detailedInfo, {"Reference Amplitude As: "<>ToString[FortranForm@Asref]}];
AppendTo[detailedInfo, {"Reference sigma8: "<>ToString[s8ref]}];
AppendTo[detailedInfo, printstand[standbasis[{}, Asref]]];
AppendTo[detailedInfo, printcamb[cambbasis[{}, Asref]]];
AppendTo[detailedInfo, {"Halofit version in CAMB: "<>ToString[halovers]}];
AppendTo[detailedInfo, {"Maximum k value in CAMB (h/Mpc): "<>ToString[Tkmax]}];
AppendTo[detailedInfo, {"Numeric CAMB parameters : "<>ToString[numcambopt]}];
AppendTo[detailedInfo, {"Constant Neutrino CAMB parameters : "<>ToString[neutrinopts]}];
AppendTo[detailedInfo, {"List of z points: "<>ToString[listofzrange]}];
If[FileExistsQ[expodir<>"Fiducial"<>"-detailedInfo.txt"],
Print["Fiducial detailedInfo file already exists"];
,
Print["Printing Fiducial detailedInfo file"];
Export[expodir<>"Fiducial"<>"-detailedInfo.txt", detailedInfo, "Table"]
];
];


If[compueps==False, Print["No epsilon files computed."]; Return[False] ];

Do[
Do[
Print["Current parameter: "];
Print[$paramoptions[[indi]][[1]] ];
Print["Relative epsilon value : "<>ToString[epsilonstep]];
$debugPrint=False;

Which[outfmt==1,
If[zzi==Reverse[Range[lenlistofz]][[1]],
   zpkTablePlus={}; zpkTableMinus={}; $debugPrint=True;

   {optplus,optminus,step}=numericalDerivativeStep[indi,epsilonstep];
   If[UnsameQ[$paramoptions[[indi]][[1]] , sigma8],
   Asoptplus=LCDMCAMBPsPre[0., {optplus}~Join~{returnInterpolated->False ,checkConsistency->True, TransferKmax->Tkmax, HalofitVersion->halovers, returnRescaledAs->True}];
   Asoptminus=LCDMCAMBPsPre[0., {optminus}~Join~{returnInterpolated->False ,checkConsistency->True, TransferKmax->Tkmax, HalofitVersion->halovers, returnRescaledAs->True}];
   AppendTo[AsrescaledList, {Asoptplus,Asoptminus}];
   AppendTo[detailedInfo, {"***Parameter Variations:  ***"}];
   AppendTo[detailedInfo, {"Relative epsilon value : "<>ToString[epsilonstep]}];
   AppendTo[detailedInfo, {"Constant sigma8 : "<>ToString[(sigma8/.$paramoptions)]}];
   AppendTo[detailedInfo, {"++Parameter plus change: "<>ToString[optplus]}];
   AppendTo[detailedInfo, {"Corresponding As value: "<>ToString[FortranForm[Asoptplus]]}];
   AppendTo[detailedInfo, printcamb[cambbasis[optplus, Asoptplus]]];
   AppendTo[detailedInfo, printstand[standbasis[optplus, Asoptplus]]];
   AppendTo[detailedInfo, {"--Parameter minus change: "<>ToString[optminus]}];
   AppendTo[detailedInfo, {"Corresponding As value: "<>ToString[FortranForm[Asoptminus]]}];
   AppendTo[detailedInfo, printcamb[cambbasis[optminus, Asoptminus]]];
   AppendTo[detailedInfo, printstand[standbasis[optminus, Asoptminus]]];
   Print[AsrescaledList];
   ,(*Else, SameQ[$paramoptions[[indi]][[1]] , sigma8] *)
   Asoptplus=LCDMCAMBPsPre[0., {optplus}~Join~{returnInterpolated->False ,checkConsistency->False, TransferKmax->Tkmax, HalofitVersion->halovers, returnRescaledAs->True}];
   Asoptminus=LCDMCAMBPsPre[0., {optminus}~Join~{returnInterpolated->False ,checkConsistency->False, TransferKmax->Tkmax, HalofitVersion->halovers, returnRescaledAs->True}];
   AppendTo[detailedInfo, {"++Parameter plus change: "<>ToString[optplus]}];
   AppendTo[detailedInfo, {"As plus: "<>ToString[FortranForm@Asoptplus]}];
   AppendTo[detailedInfo, printcamb[cambbasis[optplus, Asoptplus]]];
   AppendTo[detailedInfo, printstand[standbasis[optplus, Asoptplus]]];
   AppendTo[detailedInfo, {"--Parameter minus change: "<>ToString[optminus]}];
   AppendTo[detailedInfo, {"As minus: "<>ToString[FortranForm@Asoptminus]}];
   AppendTo[detailedInfo, printcamb[cambbasis[optminus, Asoptminus]]];
   AppendTo[detailedInfo, printstand[standbasis[optminus, Asoptminus]]];
   ];
];
(*Print["z values: ", Reverse[listofzrange[[zzi]]] ];*)
If[UnsameQ[$paramoptions[[indi]][[1]], sigma8] , $Asreference=Asoptplus; , $Asreference=Asref;];
pklistplus=LCDMCAMBPsPre[listofzrange[[zzi]], {optplus}~Join~{returnInterpolated->False ,checkConsistency->False, TransferKmax->Tkmax,
HalofitVersion->halovers}];
If[UnsameQ[$paramoptions[[indi]][[1]], sigma8], $Asreference=Asoptminus; , $Asreference=Asref;];
pklistminus=LCDMCAMBPsPre[listofzrange[[zzi]], {optminus}~Join~{returnInterpolated->False ,checkConsistency->False, TransferKmax->Tkmax,
HalofitVersion->halovers}];
pklistplus=Flatten[pklistplus,1];
pklistminus=Flatten[pklistminus,1];
zpkTablePlus=zpkTablePlus~Join~pklistplus;
zpkTableMinus=zpkTableMinus~Join~pklistminus;
If[zzi==Reverse[Range[lenlistofz]][[-1]],
If[expospectra==True,
dirnames=Flatten[((Table[(#<>"_"<>plm<>"_"<>scientificfortranform[epsilonstep, "E", "p", "eps_"]), {plm,plusmin}])&/@$paramlabels)];
mkDirectory[expodir<>dirnames[[2*indi-1]]];  mkDirectory[expodir<>dirnames[[2*indi]]];
Print["Exporting plus file..."];
Export[expodir<>dirnames[[2*indi-1]]<>"/pkz-"<>dirnames[[2*indi-1]]<>outext, zpkTablePlus, "Table"];
Print["Exporting minus file..."];
Export[expodir<>dirnames[[2*indi]]<>"/pkz-"<>dirnames[[2*indi]]<>outext, zpkTableMinus, "Table"];
];
If[expoinfo==True,
Print["Exporting info file..."];
AppendTo[detailedInfo, {"Parameter variation generation finished. Date: "<>ToString[DateString[]]}];
Export[expodir<>scientificfortranform[epsilonstep, "E", "p", "eps_"]<>"-detailedInfo.txt", detailedInfo, "Table"]
];
$Asreference=Asref;
];

,
outfmt==2,

  lastindx=Reverse[Range[lenlistofz]][[1]];

If[zzi==lastindx,$debugPrint=True,$debugPrint=False];
Print["z value: ", listofzrange[[zzi]] ];
{optplus,optminus,step}=numericalDerivativeStep[indi,epsilonstep];

{pkintpplus,log10kvalsintpplus,sigma8valofzplus, sigma8val0plus}=LCDMCAMBPsPre[listofzrange[[zzi]], {optplus}~Join~{
returnInterpolated->True, returnSigma8->True, returnGrowthRate->False, TransferKmax->Tkmax, HalofitVersion->halovers, checkConsistency->False}];
{pkintpminus,log10kvalsintpminus,sigma8valofzminus, sigma8val0minus}=LCDMCAMBPsPre[listofzrange[[zzi]], {optminus}~Join~{
returnInterpolated->True, returnSigma8->True, returnGrowthRate->False, TransferKmax->Tkmax, HalofitVersion->halovers, checkConsistency->False}];

If[zzi==lastindx,
  Print["Exporting Growth and GrowthRate at eps value"];
  {outZtabPlus,outKtabPlus,growthTabPlus,fgrowthrateTabPlus}=LCDMCAMBPsPre[listofzrange,{optplus}~Join~{returnGrowthRate->True,
  returnInterpolated->False, TransferKmax->Tkmax, HalofitVersion->halovers, checkConsistency->False}];
  {outZtabMinus,outKtabMinus,growthTabMinus,fgrowthrateTabMinus}=LCDMCAMBPsPre[listofzrange,{optminus}~Join~{returnGrowthRate->True,
    returnInterpolated->False, TransferKmax->Tkmax, HalofitVersion->halovers, checkConsistency->False}];
];

log10kvalsintpfidu=If[FileExistsQ[expodir<>"log_k_sampling"<>outext],
Import[expodir<>"log_k_sampling"<>outext, "Table"]
,
Print["fiducial k-sampling file not available, using k-sampling of plus step"];
log10kvalsintpplus
];

kvalsintpfidu=10^(log10kvalsintpfidu);

pklistplus=10^(Table[pkintpplus[kkl],{kkl,log10kvalsintpfidu}]);
pklistminus=10^(Table[pkintpminus[kkl],{kkl,log10kvalsintpfidu}]);

pklistplusovers8sq=pklistplus/(sigma8valofzplus^2);
pklistminusovers8sq=pklistminus/(sigma8valofzminus^2);

s8listplus=Table[sigma8valofzplus,{kkl,1,Length@log10kvalsintpfidu}];
s8listminus=Table[sigma8valofzminus,{kkl,1,Length@log10kvalsintpfidu}];

kvalsintpfidu = N@Flatten[kvalsintpfidu,1];
pklistplusovers8sq = N@Flatten[pklistplusovers8sq,1];
pklistminusovers8sq = N@Flatten[pklistminusovers8sq,1];  (*temporary fix for formatting*)

zpkTablePlus=Transpose[{kvalsintpfidu,pklistplusovers8sq,s8listplus}];
zpkTableMinus=Transpose[{kvalsintpfidu,pklistminusovers8sq,s8listminus}];

PrependTo[zpkTablePlus, {"#  k [h/Mpc]          Pk/s8^2 [Mpc/h]^3          s8"}];
PrependTo[zpkTableMinus, {"#  k [h/Mpc]          Pk/s8^2 [Mpc/h]^3          s8"}];

dirnames=Flatten[((Table[(#<>"_"<>plm<>"_"<>scientificfortranform[epsilonstep, "E", "p", "eps_"]), {plm,plusmin}])&/@$paramlabels)];
mkDirectory[expodir<>dirnames[[2*indi-1]]];
mkDirectory[expodir<>dirnames[[2*indi]]];

Print["Exporting plus file for Pk..."];
Export[expodir<>dirnames[[2*indi-1]]<>"/Pk-"<>dirnames[[2*indi-1]]<>"-z_"<>IntegerString[zzi-1,10,3]<>outext, zpkTablePlus, "Table"];

If[zzi==lastindx,
  Print["Exporting plus files for Growth..."];
  Export[expodir<>dirnames[[2*indi-1]]<>"/z_values_list"<>outext, outZtabPlus, "Table"];
  Export[expodir<>dirnames[[2*indi-1]]<>"/k_values_list"<>outext, outKtabPlus, "Table"];
  Export[expodir<>dirnames[[2*indi-1]]<>"/D_Growth_zk-"<>dirnames[[2*indi-1]]<>outext, growthTabPlus, "Table"];
  Export[expodir<>dirnames[[2*indi-1]]<>"/f_GrowthRate_zk-"<>dirnames[[2*indi-1]]<>outext, fgrowthrateTabPlus, "Table"];

  Print["Exporting minus files for Growth..."];
  Export[expodir<>dirnames[[2*indi]]<>"/z_values_list"<>outext, outZtabMinus, "Table"];
  Export[expodir<>dirnames[[2*indi]]<>"/k_values_list"<>outext, outKtabMinus, "Table"];
  Export[expodir<>dirnames[[2*indi]]<>"/D_Growth_zk-"<>dirnames[[2*indi]]<>outext, growthTabMinus, "Table"];
  Export[expodir<>dirnames[[2*indi]]<>"/f_GrowthRate_zk-"<>dirnames[[2*indi]]<>outext, fgrowthrateTabMinus, "Table"];
];

Print["Exporting minus file for Pk..."];
Export[expodir<>dirnames[[2*indi]]<>"/Pk-"<>dirnames[[2*indi]]<>"-z_"<>IntegerString[zzi-1,10,3]<>outext, zpkTableMinus, "Table"];




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

];

,
{zzi,Reverse[Range[lenlistofz]]}
];
,
{indi,paramindexlist}
]
];


$debugPlots=False;
$justAnalysys=False;
