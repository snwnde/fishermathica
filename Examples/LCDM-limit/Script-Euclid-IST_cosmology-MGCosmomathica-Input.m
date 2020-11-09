(* ::Package:: *)

(* ::Section:: *)
(*Initialization Cells*)


(* ::Subsection:: *)
(*Set paths to local folders and directories*)


notebookdir=If[aa=DirectoryName@$InputFileName; aa=="", NotebookDirectory[], aa]


SetDirectory[notebookdir]


$globalpackageparallel=True


Get["0-Initialization.m"]


(* ::Section:: *)
(*Run type specifications*)


$stepstring="-IST_cosmo-SteM4-nohunits-yesIA-sameeps-morezinte-";   (* String identifying this run and configuration of Fisher matrix analysis*)


$dateshort=DateString["ISODate"]    (*Date used in Results folder name*)


$debugPlots=False;      (* If set to True, intermediate plots will be produced. Needs to be set to False when running in terminal. *)
$debugPrint=False;      (* If set to True, intermediate quantities inside package functions will be printed. Set to True only when debugging a very specific function *)


$justAnalysys=False;   (*If set to True, import Fishers and run postprocessing analysis only *)


inputsfolder=notebookdir<>"../"<>"Input/";


inputDataDir=inputsfolder<>"MGCosmomathica-CommonInputFiles-IST_cosmology-LCDM/"


(* ::Section:: *)
(*Cosmological Parameters, choose fiducial values.*)


setPScosmoOptions[lcdmBool->False,linearBool->False, kdependentGrowth->True];


Get["1-Parameters.m"]


(* ::Section:: *)
(*Load cosmological functions from files*)


DirectoryQ[inputDataDir]


myInterrupt[]


Get["2-LoadInputFiles-SteM.m"]


(* ::Section:: *)
(*Load survey specifications and set other survey parameters*)


myInterrupt[]


$surveychosen="Euclid"


Get["3-Specifications.m"]


$zintepoints


(* ::Section:: *)
(*Fisher Matrix Run*)


myInterrupt[]


$newrecipelinear


(*Get["4-GC-FisherRun.m"]*)


myInterrupt[]


Get["5-WL-FisherRun.m"]



