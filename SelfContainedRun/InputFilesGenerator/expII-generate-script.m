(* ::Package:: *)

DirectoryName@$InputFileName
notebookdir=If[aa=DirectoryName@$InputFileName; aa=="", NotebookDirectory[], aa]
SetDirectory[notebookdir];

$nomDeModele="Exp II"

Get[notebookdir<>"Generate-InputFiles.m"]



