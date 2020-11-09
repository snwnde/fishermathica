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


$stepstring="-WP6-HS5-SteM4-int1-";   (* String identifying this run and configuration of Fisher matrix analysis*)


$dateshort=DateString["ISODate"]    (*Date used in Results folder name*)


$debugPlots=False;      (* If set to True, intermediate plots will be produced. Needs to be set to False when running in terminal. *)
$debugPrint=False;      (* If set to True, intermediate quantities inside package functions will be printed. Set to True only when debugging a very specific function *)


$justAnalysys=False;   (*If set to True, import Fishers and run postprocessing analysis only *)


inputsfolder=notebookdir<>"../"<>"Input/";


inputDataDir=inputsfolder<>"WP6-fofR-commonInput-Cosmomathica-InputFiles-baseline_cosmology-DirectPk--Mnu_0p06-HS5/"


(* ::Section:: *)
(*Cosmological Parameters, choose fiducial values.*)


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


(*Get["4-GC-FisherRun.m"]*)


myInterrupt[]


Get["5-WL-FisherRun.m"]



