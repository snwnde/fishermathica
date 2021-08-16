(* ::Package:: *)

(* ::Section:: *)
(*Initialization*)


DirectoryName@$InputFileName
notebookdir=If[aa=DirectoryName@$InputFileName; aa=="", NotebookDirectory[], aa]
SetDirectory[notebookdir];

$globalpackageparallel=True
Get[notebookdir<>"/../Initialization.m"]
$nameOfModel=$nomDeModele
Get[notebookdir<>"GenerationParameters.m"]
setPScosmoOptions[lcdmBool->True,linearBool->True, kdependentGrowth->False];


(* ::Section:: *)
(*Cosmological Parameters, choose fiducial values.*)


parametertestIndex=0;

Get[notebookdir<>"/../0-Runtype.m"]
Get[notebookdir<>"/../1-Parameters.m"]



(* {OmegaK0Today[],OmegaM0Today[],OmegaCDM0Today[],OmegaBaryon0Today[],OmegaDE0Today[],1-OmegaDE0Today[]-OmegaM0Today[],$OmegaR,hubbleToday[]} *)


(* ::Section::Bold:: *)
(*Generate Many Files*)


$filesdir=findInputFilesDir[$nameOfModel];

zlist=Flatten@Import[$filesdir<>"z_values_list.txt","Table"];
klist=Flatten@Import[$filesdir<>"k_values_list.txt","Table"];


$zrange=zlist;
(*$zrange={1.0, 1.2, 1.4, 1.65}*)
lenlistofz=Length[$zrange]
$krange=klist;
(* $krange=Select[klist,#<2&]; *)
lenlistofk=Length[$krange]



(* 
paramsFolderTable[parslist_]:=Flatten[Table[pp<>"_"<>pm<>"_"<>epsstr, {pp,parslist},{pm,plusminstr}],1];


paramsFolderTable[folderpars]


$paramsDirectoryNames={"fiducial"}~Join~paramsFolderTable[folderpars] *)


(* ::Subsection:: *)
(*Produce Files*)


OmegaBaryon0Today[]


powerSpectrum[0.1,0.2]


writeInputeFiles[]:=Switch[StringMatchQ["Exp*"][$nameOfModel],
True,

Module[{omegaM,omegaB,nS,AS109,Gamma,dir},
Do[
omegaM=OmegaM0Today[];omegaB=OmegaBaryon0Today[];nS=nsref;AS109=as109ref;Gamma=gammaref;
Which[
StringMatchQ[ii,"*Ome*M*"],omegaM=OmegaM0Today[]*(1+epsilonstep),
StringMatchQ[ii,"*Ome*B*"],omegaB=OmegaBaryon0Today[]*(1+epsilonstep),
StringMatchQ[ii,"*ns*"],nS=nsref*(1+epsilonstep),
StringMatchQ[ii,"*As*"],AS109=as109ref*(1+epsilonstep),
StringMatchQ[ii,"*gamma*"],Gamma=gammaref*(1+epsilonstep),
StringMatchQ[ii,"fiducial"],None
];

folderpar = paramsFolderFunc[ii, epsilonstep];
If[SameQ[folderpar, None],
(*Parameter and epsilon combination does not exist, continue to next step of the loop*)
Continue[]
];

dir=$filesdir<>paramsFolderFunc[ii,epsilonstep];
mkDirectory[dir];
(*If[StringMatchQ["Exp II*"][$nameOfModel],gammaLowerBoundary=121.2;gammaUpperBoundary=135;,gammaLowerBoundary=0;gammaUpperBoundary=200];*)
If[StringMatchQ["Exp II*"][$nameOfModel],gammaLowerBoundary=0.9495*gammaref;gammaUpperBoundary=1.058*gammaref;,gammaLowerBoundary=0;gammaUpperBoundary=200];
If[Gamma<gammaLowerBoundary||Gamma>gammaUpperBoundary,Print["gamma is too small or too big. Skip file generation."],

  If[FileExistsQ[dir<>"/background_Hz.txt"], Print["File already exists in "<>paramsFolderFunc[ii,epsilonstep]],
  HzTable=Table[hubbleToday[Omegam->omegaM,Omegab->omegaB,ns->nS,As109->AS109,gamma->Gamma]*dimensionlessHubbleWrtN[NfOfz[zz]]*100/$lightspeed,{zz,$zrange}];
  Export[dir<>"/background_Hz.txt", HzTable, "Table"];
  ];

  If[FileExistsQ[dir<>"/Plin-zk.txt"], Print["File already exists in "<>paramsFolderFunc[ii,epsilonstep]],
  zpkTable=Table[powerSpectrum[zz,kk,Omegam->omegaM,Omegab->omegaB,ns->nS,As109->AS109,gamma->Gamma],{zz,$zrange},{kk,$krange}];
  Export[dir<>"/Plin-zk.txt", zpkTable, "Table"];
  ];

  If[FileExistsQ[dir<>"/D_Growth-z.txt"], Print["File already exists in "<>paramsFolderFunc[ii,epsilonstep]],
  dzkTable=Table[Growth[zz,Omegam->omegaM,Omegab->omegaB,ns->nS,As109->AS109,gamma->Gamma],{zz,$zrange}];
  Export[dir<>"/D_Growth-z.txt", dzkTable, "Table"];
  ];
 
  If[FileExistsQ[dir<>"/f_GrowthRate-z.txt"], Print["File already exists in "<>paramsFolderFunc[ii,epsilonstep]],
  fzTable=Table[fGrowthRate[zz,Omegam->omegaM,Omegab->omegaB,ns->nS,As109->AS109,gamma->Gamma],{zz,$zrange}];
  Export[dir<>"/f_GrowthRate-z.txt", fzTable, "Table"];
  ];
  
  If[FileExistsQ[dir<>"/sigma8-z.txt"], Print["File already exists in "<>paramsFolderFunc[ii,epsilonstep]],
  sigma8zTable=Table[sigma8ofZ[zz,Omegam->omegaM,Omegab->omegaB,ns->nS,As109->AS109,gamma->Gamma],{zz,$zrange}];
  Export[dir<>"/sigma8-z.txt", sigma8zTable, "Table"];
  ];

  If[FileExistsQ[dir<>"/z_values_list.txt"], Print["File already exists in "<>paramsFolderFunc[ii,epsilonstep]],
  Export[dir<>"/z_values_list.txt", zlist, "Table"];
  ];

  If[FileExistsQ[dir<>"/k_values_list.txt"], Print["File already exists in "<>paramsFolderFunc[ii,epsilonstep]],
  Export[dir<>"/k_values_list.txt", klist, "Table"];
  ];
  ],


{ii,$folderpars},{epsilonstep,$epslistPM}]
];,

False,

Module[{omegaM,omegaB,nS,AS109,h,dir},
Do[
omegaM=OmegaM0Today[];omegaB=OmegaBaryon0Today[];nS=nsref;AS109=as109ref;h=hubbleref;
Which[
StringMatchQ[ii,"*Ome*M*"],omegaM=OmegaM0Today[]*(1+epsilonstep),
StringMatchQ[ii,"*Ome*B*"],omegaB=OmegaBaryon0Today[]*(1+epsilonstep),
StringMatchQ[ii,"*ns*"],nS=nsref*(1+epsilonstep),
StringMatchQ[ii,"*As*"],AS109=as109ref*(1+epsilonstep),
StringMatchQ[ii,"*hubble*"],h=hubbleref*(1+epsilonstep),
StringMatchQ[ii,"fiducial"],None
];

folderpar = paramsFolderFunc[ii, epsilonstep];
If[SameQ[folderpar, None],
(*Parameter and epsilon combination doesr not exist, continue to next step of the loop*)
Continue[]
];

dir=$filesdir<>paramsFolderFunc[ii,epsilonstep];
mkDirectory[dir];

If[FileExistsQ[dir<>"/background_Hz.txt"], Print["File already exists in "<>paramsFolderFunc[ii,epsilonstep]],
  HzTable=Table[hubbleToday[Omegam->omegaM,Omegab->omegaB,ns->nS,As109->AS109,hubble->h]*dimensionlessHubbleWrtN[NfOfz[zz]]*100/$lightspeed,{zz,$zrange}];
  Export[dir<>"/background_Hz.txt", HzTable, "Table"];
  ];

If[FileExistsQ[dir<>"/Plin-zk.txt"], Print["File already exists in "<>paramsFolderFunc[ii,epsilonstep]],
  zpkTable=Table[powerSpectrum[zz,kk,Omegam->omegaM,Omegab->omegaB,ns->nS,As109->AS109,hubble->h],{zz,$zrange},{kk,$krange}];
  Export[dir<>"/Plin-zk.txt", zpkTable, "Table"];
  ];

If[FileExistsQ[dir<>"/D_Growth-z.txt"], Print["File already exists in "<>paramsFolderFunc[ii,epsilonstep]],
  dzkTable=Table[Growth[zz,Omegam->omegaM,Omegab->omegaB,ns->nS,As109->AS109,hubble->h],{zz,$zrange}];
  Export[dir<>"/D_Growth-z.txt", dzkTable, "Table"];
  ];
 
 If[FileExistsQ[dir<>"/f_GrowthRate-z.txt"], Print["File already exists in "<>paramsFolderFunc[ii,epsilonstep]],
  fzTable=Table[fGrowthRate[zz,Omegam->omegaM,Omegab->omegaB,ns->nS,As109->AS109,hubble->h],{zz,$zrange}];
  Export[dir<>"/f_GrowthRate-z.txt", fzTable, "Table"];
  ];
  
  If[FileExistsQ[dir<>"/sigma8-z.txt"], Print["File already exists in "<>paramsFolderFunc[ii,epsilonstep]],
  sigma8zTable=Table[sigma8ofZ[zz,Omegam->omegaM,Omegab->omegaB,ns->nS,As109->AS109,hubble->h],{zz,$zrange}];
  Export[dir<>"/sigma8-z.txt", sigma8zTable, "Table"];
  ];

  If[FileExistsQ[dir<>"/z_values_list.txt"], Print["File already exists in "<>paramsFolderFunc[ii,epsilonstep]],
  Export[dir<>"/z_values_list.txt", zlist, "Table"];
  ];

  If[FileExistsQ[dir<>"/k_values_list.txt"], Print["File already exists in "<>paramsFolderFunc[ii,epsilonstep]],
  Export[dir<>"/k_values_list.txt", klist, "Table"];
  ];,
  
 {ii,$folderpars},{epsilonstep,$epslistPM}]
];


]


generateInputFiles[]:=Module[{AlreadyGenerate=True},

Do[

folderpar = paramsFolderFunc[ii, epsilonstep];
If[SameQ[folderpar, None],
(*Parameter and epsilon combination does not exist, continue to next step of the loop*)
Continue[]
];

dir=$filesdir<>paramsFolderFunc[ii,epsilonstep];
If[FileExistsQ[dir<>"/background_Hz.txt"], None, AlreadyGenerate=False;writeInputeFiles[]];
If[FileExistsQ[dir<>"/Plin-zk.txt"], None, AlreadyGenerate=False;writeInputeFiles[]];
If[FileExistsQ[dir<>"/D_Growth-z.txt"], None, AlreadyGenerate=False;writeInputeFiles[]];
If[FileExistsQ[dir<>"/f_GrowthRate-z.txt"], None, AlreadyGenerate=False;writeInputeFiles[]];
If[FileExistsQ[dir<>"/sigma8-z.txt"], None, AlreadyGenerate=False;writeInputeFiles[]];
If[FileExistsQ[dir<>"/z_values_list.txt"], None, AlreadyGenerate=False;writeInputeFiles[]];
If[FileExistsQ[dir<>"/k_values_list.txt"], None, AlreadyGenerate=False;writeInputeFiles[]];,

{ii,$folderpars},{epsilonstep,$epslistPM}];



If[AlreadyGenerate==True,Print["The input files were already generated."]]
];

generateInputFiles[]



