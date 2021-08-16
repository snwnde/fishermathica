DirectoryName@$InputFileName
notebookdir=If[aa=DirectoryName@$InputFileName; aa=="", NotebookDirectory[], aa]
SetDirectory[notebookdir];

$nomDeModele="lcdm"

Get[notebookdir<>"Generate-InputFiles.m"]