(* ::Package:: *)

BeginPackage["UsefulTools`"]


extractRuleSide::usage="extractRuleSide[rulelist_,side_] extracts from a rulelist, the rhs side or the lhs side, format: lhs or rhs"

myColorData::usage="my favorite Color Data for discrete values, ColorData=63"

styTexFunc::usage="styTexFunc[text_,size_:Large], Style a text string in TeX mode, using default size Large"

styFunc::usage="Style a text string using default size 22, Bold, and 'Times' fontfamily as default optional arguments.
Form: styFunc[text_, size:22, w:Bold, family:'Times']:"


replaceStringPoint::usage="replaceStringPoint[val_,prestr_String,sep_String:'_',pstr_String:p]: Replaces the decimal point in val,
with the default letter p and prepends to val, the string 'prestr', with the default separator '_'"

parselineWithNameEqualValue::usage="given a string it looks for a pattern of a parameter equal to a digit value. It returns a list with the parameter name and the numeric value.
The regular expression to be matched is given by the internal variable 'regexval'"

importNumericTxtTable::usage="Import a table or matrix from an ASCII text file, neglecting rows or columns which contain headers with strings that are not numeric.
Form: importNumericTxtTable@fileName"
numericImport::usage="Import a table or matrix from an ASCII text file,
neglecting rows or columns which contain headers with strings that are not numeric.
Form: numericImport@fileName"

importValueOfParameterFromTxtFile::usage="Import from a text file only the value of the parameter corresponding to the parampattern string pattern.
 Options: IgnoreCase->True.  Form: importValueOfParameterFromTxtFile[file_String,parampattern_, opts:OptionsPattern[]]"


mkDirectory::usage="Create a directory and return the directory name. It works even if directory exists already. If directory cannot be created,
it returns an error message."

removetrailinghyphens::usage="Remove trailing hyphens from a string."


numfill::usage="numfill[n,ffillString(optional)], Pad a list of numbers or a number with the fillString, so that the final number is n digits long. Fill string is by default '0'"

myInterrupt::usage="myInterrupt[mess_:'Session '], Interrupts evaluation and outputs a message with the total session time"

rescaleApplyFunctionTable::usage="Form: rescaleApplyFunctionTable[funcy1_,funcx2_:(1&), funcx1_:(1#&)], applies pure functions
to the columns of a 2dim table. If table is formed by elements {x,y}, it returns {funcx1@x, funcx2@x*funcy1@y)}.
Default return: {x,funcy1@y}"

applyFunctionToTable::usage="Form: applyFunctionToTable[funcx1_:(1#&),funcy1_:(1#&), funcx2_:(1&), funcy2_:(1&)],
applies pure functions to the x and y columns of a 2-dimensional table.
Default: funcx1 and funcy1 are the identity function, while funcx2 and funcy2 are the number 1.
{x,y}->{funcx1@x * funcy2@y, funcy1@y * funcx2@x}"

ratioOfLists::usage="Form: ratioOfLists[list1,list2] or ratioOfLists[{list1,list2}].
Calculates for two lists list1 and list2, which have equal 'x' values and different 'y' values, the ratio of the 'y' values, keeping the 'x' intact."


logarithmicDivisions::usage="Get n uniformly spaced divisions in Logspace  between x and y.
Form:logarithmicDivisions[{x,y},n]"

tentox::usage="Form: tentox=(10^#)&; Pure function that computes 10^x."

removeDuplicateXY::usage="Remove from a table of X-Y pairs, all the rows in which X is repeated."

twoArgumentsCheck::usage="Form: twoArgumentsCheck[txarg__,bool_]. Checks if the arguments passed as txarg are composed of one or two arguments,
and returns the arguments separately."

timeElapsed::usage="timeElapsed[initialDate_, finalDate_].
Compute elapsed time in hours, minutes, seconds, between two DateTime strings."

commonMatrixProperties::usage="Prints common information about a Matrix."

removeDownValues::usage="Removes DownValues of a function f specified with arguments matching patterns. Form: removeDownValues[p:f_[___]].
Used for example to remove memoized values of a function, without removing the function itself.
With pattern matching one can remove just certain DownValues as needed.
To remove the function definition completely, give f with arguments matching every possible pattern"


parallelmemo::usage="Form: parallelmemo[f_Symbol][x__]. Takes a function f, with a variable number of arguments and options x__ and memoizes it on the main
and the parallel kernels."
memo::usage="Form: memo[___]. Auxiliary function for parallelmemo.
Accepts function and arguments and memoizes them."
valueOnMain::usage="Form: valueOnMain[f_,x__]. Auxiliary function for parallelmemo. Saves value of function and arguments on main kernel."
setValueOnMain::usage="Form: setValueOnMain[f_,x__,fx_]. Auxiliary function for parallelmemo. Assigns UpValues to f."
clearParallelCache::usage="Form: clearParallelCache[f_Symbol]. Function that clears the parallel cache of the function 'f' created by the function parallelmemo."

memoryUseMB::usage="Returns current memory used in storing Mathematica objects, in MB."

messageFormatText::usage="messageFormatText[text_String,options:OptionsPattern[Style]]:
This function formats the text with StyleBox in order to be used for generating formatted Message[] messages.
Accepts all the options from Style."


$blueishTones::usage="Triplet of blueish colors for contour plots"
$redishTones::usage="Triplet of redish colors for contour plots"
$yellowishTones::usage="Triplet of yellowish colors for contour plots"

$linecolorsty1::usage="small, small dashed line and greenish color style for plots."
$linecolorsty2::usage="small, medium dashed line and blueish color style for plots."
$linecolorsty3::usage="medium, medium dashed line and orangeish color style for plots."
$linecolorsty4::usage="large, medium dashed line and pinkish color style for plots."


numberFormat::usage="Option for exportSimpleFormatTable. Default: NumberForm. Specify an output form for display, like NumberForm or ScientificForm."
formArgument::usage="Option for exportSimpleFormatTable. Default: {3,2}. Argument of the specified numberFormat function."
fortranOutput::usage="Option for exportSimpleFormatTable. Default: False. Specify if a fortran output in scientific notation is wanted."
texOutput::usage="Option for exportSimpleFormatTable. Specify if output to LaTeX is wanted. Output is written as an array in display math mode."
exportSimpleFormatTable::usage="exportSimpleFormatTable[filename, table or several arguments (each containing a column vector), options].
 This function exports a table into a text or tex file, using the specified formatting in the options."
scientificfortranform::usage="Function that converts any decimal number to scientific form and the to Fortran Form.
Arguments: scientificfortranform[num_, expstr_String:'e', pstr_String:'p', prestr_String:'', sep_String:''].
expstr: string for denoting the exponent,
pstr: string for denoting the decimal point,
prestr: string for labeling the number, sep: string to separate label and number."


writeTabletoTeXArray::usage="writeTabletoTeXArray[table, filename, options]. This function exports a numeric table into an array in LaTeX in display math mode.
It also writes a predefined preamble that allows the tex file to be compiled. Among other options, it can highlight certain values of the entries with a specified color."
columnAlignment::usage="Option for writeTabletoTeXArray. Options: 'c', 'l' or 'r'. Sets the alignment of each column in the array."
highlightColorString::usage="Option for writeTabletoTeXArray. Specify a valid LaTeX color for highlighting certain entries of the array."
highlightFunction::usage="Option for writeTabletoTeXArray. Specify a pure function containing a conditional test. Each entry of the array will be checked acording to this
 test function. If the output is True, the entry will be highlighted using the color specified in highlightColorString."


printConditionalInformation::usage="A redefinition of the print function in which a bool variable passed as option, decides if the text passed is to be printed or not. Options: printBool->True (Default)"

existsFile::usage="Check if filename exists as a file, otherwise prints message and Abort. Form: existsFile[filename]"

existsFile::nofile="Text file `1` not found. Aborting..."

$debugPrint::usage="Boolean variable that activates Prints in debugging mode"

debugPrint::usage="Use this instead of Print when printing variables inside package fucntions for debugging.
This avoids problems with printing inside Plots. Only prints if $debugPrint is set to True."

EndPackage[]


BeginPackage["UsefulTools`"]
Begin["`Private`"]

myColorData=ColorData[63,"ColorList"];

extractRuleSide[rulelist_,side_]:=Block[{hs=side},Cases[rulelist,HoldPattern[lhs_->rhs_]->hs]]


rescaleApplyFunctionTable[funcy1_,funcx2_:(1&), funcx1_:(1#&)]:={funcx1@First[#],(funcx2@First[#]*funcy1@Last[#])}&;

applyFunctionToTable[funcx1_:(1#&),funcy1_:(1#&), funcx2_:(1&), funcy2_:(1&)]:=
Block[{table},table=(#/.{{x_,y_}->{(funcx1@x)*(funcy2@y),(funcy1@y)*(funcx2@x)}}); table]&

ratioOfLists[list1_, list2_]:=Transpose[{list1[[All,1]],list2[[All,2]]/list1[[All,2]]}];
ratioOfLists[listoflist_List]:=With[{l1=listoflist[[1]],l2=listoflist[[2]]}, ratioOfLists[l1, l2]];



existsFile[filename_]:=If[(FileExistsQ[filename] && DirectoryQ[filename]==False ), debugPrint["File"<>ToString@filename<>"  exists"], Message[existsFile::nofile, filename]; Interrupt[]];

$debugPrint=False;

debugPrint[var__, switch_:0]:=Block[{test=$debugPrint},
  If[test==True && switch==0,
    Print[var];]]

toDirective[{ps__}|ps__]:=Flatten[Directive@@Flatten[{#}]]&/@{ps}

styleJoin[style_,base_]:=Module[{ps,n},ps=toDirective/@{PlotStyle/.Options[base],style};
ps=ps/.Automatic:>Sequence[];
n=LCM@@Length/@ps;
MapThread[Join,PadRight[#,n,#]&/@ps]]

pp={Plot, ListLogLogPlot,LogLogPlot, LogLinearPlot, LogPlot};

Unprotect/@pp;

(#[a__,b:OptionsPattern[]]:=Block[{$alsoPS=True,sh},sh=Cases[{b},("MessagesHead"->hd_):>hd,{-2},1]/.{{z_}:>z,{}->#};
With[{new=styleJoin[OptionValue[PlotStyle],sh]},#[a,PlotStyle->new,b]]]/;!TrueQ[$alsoPS];
DownValues[#]=RotateRight[DownValues@#];(*fix for versions 9 and 10*))&/@pp;

SetOptions[pp,PlotStyle->myColorData];

styTexFunc[text_,size_:Large]:=Text[Style[ToExpression[text, TeXForm, HoldForm], size]]

styFunc[text_, size_:22, w_:Bold, family_:"Times"]:=Style[text, size, w, FontFamily->family]
SetAttributes[styFunc, Listable]


replaceStringPoint[val_,prestr_String,sep_String:"_",pstr_String:"p"] := Block[{repval},
repval=StringReplace[ToString[val], "." -> pstr];
prestr <> sep <> repval]

SetAttributes[replaceStringPoint, Listable]

ClearAll[scientificfortranform]
scientificfortranform[num_Integer, expstr_String:"e", pstr_String:".", prestr_String:"", sep_String:""]:=(prestr <> sep <> ToString[num]);
scientificfortranform[num_Real, expstr_String:"e", pstr_String:".", prestr_String:"", sep_String:"", sigdigts_Integer:2, decdigts_Integer:1]:=replaceStringPoint[ToString[NumberForm[NumberForm[num,
  NumberFormat->(#1&),ExponentFunction->(#&)],{sigdigts,decdigts}]], prestr,sep,pstr]<>expstr<>ToString[NumberForm[num,
  NumberFormat->(#3&),ExponentFunction->(#&)]];
scientificfortranform[rul_Rule, opt___]:=Rule[rul[[1]], scientificfortranform[rul[[2]]], opt ];
SetAttributes[scientificfortranform, Listable]


scientificfortranform[num_String, opt___]:=ToString[num];

scientificfortranform[num_, opt___]:=ToString[num];

SetAttributes[scientificfortranform, Listable]

numericImport = (Select[Import[#, "Table"], VectorQ[#, NumericQ] &]) &;

importNumericTxtTable=Select[Import[#, "Table"], VectorQ[#, NumericQ]&]&;


regexval = "(\\w+)\\s*=\\s*(\\d+\\.?\\d*)"
(*This regular expression matches a word character (letter,
  digit or _),followed by any amount of whitespace,
  followed by an equal sign,followed by any amount of white space,
  followed by one or more digits,
  followed by one or zero occurences of a dot,
  followed by zero or more digits*)

parselineWithNameEqualValue[strline_String] :=
    Block[{pars},
      pars = StringCases[
        strline,(*this one works if the passed strline is indeed just a string*)
        RegularExpression[regexval] -> {"$1", "$2"}];
      Return[Flatten[pars]]
    ]


parselineWithNameEqualValue[strline_List]:=Block[{pars, strl},
  strl = StringJoin@(ToString[#]&/@strline);  (*this one works if the passed strline is a list of strings and/or numbers *)
  pars = StringCases[strl,
    RegularExpression[regexval] -> {"$1", "$2"}];
  pars = Flatten[pars];
  pars = {ToString@pars[[1]], ToExpression@pars[[2]]};
  Return[pars]
]

Options[importValueOfParameterFromTxtFile]={IgnoreCase->True};

importValueOfParameterFromTxtFile::patt="Either pattern `1` not found or this pattern is not followed by an equal sign and a numerical value in the text.";

importValueOfParameterFromTxtFile[file_String,parampattern_,opts:OptionsPattern[]]:=Block[{par=parampattern,st,val, filt},
  filt=Select[ReadList[file,String],StringMatchQ[#,par,opts]&];
  st = Select[filt,!StringMatchQ[#,RegularExpression[".*#.*"]]&]; (*filter lines containing comments with #*)
  If[st=={},Message[importValueOfParameterFromTxtFile::patt, par];Return[Null]];
  val=parselineWithNameEqualValue[st[[1]]];
  If[val=={},Message[importValueOfParameterFromTxtFile::patt, par];Return[Null]];
  val = ToExpression[val[[2]]]
]


mkDirectory[str_?StringQ]:=Module[{dir},
  If[SameQ[DirectoryQ[str],True],dir=str;
  Return[dir],dir=Check[CreateDirectory[str],
  Print["Creation of directory "<>str<>" Failed"];
  $Failed];
  Return[dir]]
]



numfill[n_,fill_String:"0"]:=Function[{s},StringJoin@@PadLeft[Characters@ToString@s,n,fill]];

$startTime=SessionTime[];

myInterrupt[mess_:"Session "]:=Module[{sess=SessionTime[],diff,str},str=If[mess,ToString[mess]];
diff=sess-$startTime;
Print[DateString[]];
Print[mess<>"Time elapsed since initialization: ", diff];
FrontEndExecute[FrontEndToken["FindEvaluatingCell"]];
Interrupt[]]

tentox=(10^#)&;

logarithmicDivisions[{x_?Positive,y_?Positive},n_Integer/;n>1]:= (x (y/x)^Range[0,1,1/(n-1)])


removeDuplicateXY[tableXY_List]:=Map[First,GatherBy[tableXY,First],1]



twoArgumentsCheck[txarg__,bool_]:=Block[{len,outlist,tv,xv},
  len=Length@List@txarg;
  If[bool==True && len==2,
    tv=First@List@txarg; xv=Last@List@txarg,
    tv=First@List@txarg; xv=Hold[Sequence@@{}]];
  outlist={tv,xv}
]


commonMatrixProperties[mat_?MatrixQ]:=Block[{dim,sym,posd,diag},dim=ToString[Dimensions[mat]];
sym=ToString@SymmetricMatrixQ[mat];
posd=ToString@PositiveDefiniteMatrixQ[mat];
(*diag=ToString@DiagonalizableMatrixQ[mat];*)
Print["Dimensions: "<>dim<>"  Symmetric: "<>sym<>"  Positive definite: "<>posd(*<>"  Diagonalizable: "<>diag*)]];




$blueishTones={ColorData["Crayola"]["Indigo"],Cyan,ColorData["Crayola"]["TealBlue"]};
$redishTones={ColorData["Crayola"]["Scarlet"],ColorData["Crayola"]["Salmon"],ColorData["Crayola"]["Magenta"]};
$yellowishTones={ColorData["Crayola"]["Sunglow"],ColorData["Crayola"]["YellowOrange"],ColorData["Crayola"]["Shadow"]};

$linecolorsty1={Directive[RGBColor[0.243137,0.752941,0.556863],Dashing[{Small,Small}]],RGBColor[0.313725,0.752941,0.290196],RGBColor[0.129412,0.482353,0.101961]};
$linecolorsty2 = {Directive[RGBColor[0.494118,0.698039,1.],Dashing[{Small,Medium}]],RGBColor[0.117647,0.498039,1.],RGBColor[0.109804,0.0431373,0.498039]};
$linecolorsty3 = {Directive[RGBColor[1.,0.733333,0.109804],Dashing[{Large,Medium}]],RGBColor[0.92549,0.65098,0.101961],RGBColor[1.,0.490196,0.0745098]};
$linecolorsty4 = {Directive[RGBColor[1.,0.25098,0.988235],Dashing[{Large,Medium}]],RGBColor[0.894118,0.0901961,1.],RGBColor[0.85098,0.215686,0.576471]};



matrixConditionNumber[mat_?MatrixQ]:=Module[{condnumb, error},
  condnumb=LinearAlgebra`MatrixConditionNumber[mat];
  Print["Condition Number: ", cond, " You might lose "<>ToString[Log10[condnumb]]<>"digits of accuracy"];
  If[condnumb>(1/$MachineEpsilon), Print["Matrix is ill-conditioned. You might want to check your algorithms."]]
];




Options[exportSimpleFormatTable]={numberFormat->NumberForm,formArgument->{3,2},fortranOutput->False, texOutput->False};
exportSimpleFormatTable[filename_?StringQ,vectorsortab__,opts:OptionsPattern[{exportSimpleFormatTable,writeTabletoTeXArray}]]:=Module[{len,list,
  formfu,formarg,numfun,
  formattedtabl, tabl, optlat},
list=List@vectorsortab;
len=Length@list;
Which[len>1,
tabl=Thread[list],
len==1,
tabl=vectorsortab];
numfun=OptionValue[numberFormat];
formarg=OptionValue[formArgument];
formfu=Function[numfun[#,formarg]];
If[OptionValue[fortranOutput]==True,formfu=Function[OutputForm@ScientificForm[#,NumberFormat->(Row[{#1,"e",If[ToString[#3]=="","0",#3]}]&)]]
];
formattedtabl=MapAt[formfu,tabl,{All,All}];
If[OptionValue[texOutput]==False,
Export[filename,formattedtabl,"Table"];
,
optlat=FilterRules[{opts},{ Options[writeTabletoTeXArray]}];
writeTabletoTeXArray[formattedtabl,filename,optlat];
]
];




Options[writeTabletoTeXArray]={columnAlignment->"c", highlightColorString->"red", highlightFunction->((If[#>0, True])&)};

writeTabletoTeXArray[tabl_List,fileName_String,opts:OptionsPattern[]]:=Block[{line1, lineall, centstr, lencol=Dimensions[tabl][[2]],lenrow=Dimensions[tabl][[1]], file, entryij, posi, funchi, funccol,
histr=""},
posi=OptionValue[columnAlignment];
funccol=OptionValue[highlightColorString];
funccol="\\color{"<>ToString[funccol]<>"}";
funchi = OptionValue[highlightFunction];
file=OpenWrite[fileName];
centstr=Fold[StringJoin[posi,ToString[#]]&,"",ConstantArray[posi,lencol]];
WriteLine[file, "\\documentclass[]{article}"];
WriteLine[file, "\\usepackage{color}"];
WriteLine[file,"\\begin{document}"];
WriteLine[file,"\\["];
WriteLine[file, "\\left("];
WriteLine[file,"\\begin{array}{"<>centstr<>"}"];
lineall="";
Do[
line1="";
Do[
entryij=tabl[[ii,jj]];
(*Print[entryij];
Print[funchi@(ToExpression@ToString@entryij)];*)
If[SameQ[funchi@(ToExpression@ToString@entryij), True],
histr=(funccol<>" "),
histr=""
];
If[jj<lencol,
line1=line1<>histr<>ToString@(tabl[[ii,jj]])<>" "<>"&"<>" ";
,
(*Print[jj];*)
line1=line1<>histr<>ToString@(tabl[[ii,jj]])<>" "<>"\\\\";
(*Print[line1];*)
];
,{jj,1,lencol}
];
WriteLine[file, line1];
lineall=lineall<>line1;
,
{ii,1,lenrow}
];
Print[lineall];
Print["Printed to file: "<>ToString[fileName]];
WriteLine[file,"\\end{array}"];
WriteLine[file,"\\right)"];
WriteLine[file,"\\]"];
WriteLine[file,"\\end{document}"];
Close[file];
];


removetrailinghyphens=StringReplace[#, {RegularExpression["-/$"]->"/",
  RegularExpression["-$"]->"", RegularExpression["^-"]->"", RegularExpression["-\\."]->"."}]&;

Options[printConditionalInformation]={printBool->True};

printConditionalInformation[text_, opts:OptionsPattern[]]:=Block[
  {pbool=OptionValue[printBool]},
  If[pbool==True,
   Print[text]
  ];
]

(*  from mathematica stackexchange, modfied to do actual memoization and to accept any number of arguments. ***)



parallelmemo[f_Symbol][x__]:=With[{result=memo[f,x]},
  If[result===Null,
    With[{value=valueOnMain[f,x]},
      If[value===Null,
        f/:memo[f,x]=setValueOnMain[f,x,f[x]]
        ,
        f[x]=value;
        f/:memo[f,x]=value]
    ],
    f[x]=result;
    result]
]

memo[___]:=Null
(*DistributeDefinitions[memo];*)

valueOnMain[f_,x__]:=memo[f,x]

setValueOnMain[f_,x__,fx_]:=f/:memo[f,x]=fx

(*SetSharedFunction[valueOnMain,setValueOnMain]*)

memoryUseMB[]:=N[MemoryInUse[]/1024/1024];

clearParallelCache[f_Symbol]:=(UpValues[f]={};ParallelEvaluate[UpValues[f]={}];)

(***  from mathematica stackexchange, modfied to do actual memoization and to accept any number of arguments.*)


SetAttributes[removeDownValues,HoldAllComplete];
removeDownValues[p:f_[___]]:=DownValues[f]=DeleteCases[DownValues[f,Sort->False],HoldPattern[Verbatim[HoldPattern][p]:>_]];


messageFormatText[text_String,optis:OptionsPattern[Style]]:=ToString[DisplayForm[StyleBox[text,optis]],StandardForm];

timeElapsed[initialDate_, finalDate_, prettyPrint_:True]:=Block[{ss,mm,hh, secs, mins},
secs = UnixTime[finalDate]-UnixTime[initialDate];
mins = Floor[secs,60]/60;
ss = Mod[secs,60];
mm = Mod[mins, 60];
hh = Floor[mins, 60]/60;
If[prettyPrint,
Print["Elapsed Time: ", "HH: ",hh, " MM: ",mm, " SS: ",ss];
];
Return[secs]
];
(*Example plot*)
(*ListLogLogPlot[(nlPSList["YA1"][replaceStringPoint["z", ztoplot]][replaceStringPoint["kref", #]])&/@koutvals,
PlotLegends -> Placed[LineLegend["YA1-knl="<>ToString[#]&/@koutvals, LegendFunction -> "Frame",
LegendLayout -> "Column", LegendMarkers -> Automatic, LegendMarkerSize -> Large], {0.8, 0.85}], Frame -> True, ImageSize -> 1200, Joined -> True,
PlotMarkers -> Automatic, FrameLabel -> {styTexFunc["k"], styTexFunc["P(k)"]}, PlotRange -> {{0.001, 2}, Automatic},
Epilog -> {Inset[LogLogPlot[Evaluate[(Interpolation[nlPSList["YA1"][replaceStringPoint["z", ztoplot]][replaceStringPoint["kref", #]]][kk]/Interpolation[nlPSList["YA1"][replaceStringPoint["z", ztoplot]]["kref-0p001"]][kk])&/@koutvals],
{kk, 0.001, 1.5},PlotRange\[Rule]All, ImageSize -> 150, Frame -> True, FrameLabel -> {Style["k", Bold, 14], Style["Ratio to knl=0.001", Bold, 14]}], Scaled[{0.5, 0.3}], Right, Scaled[0.5]] ,
Inset[Framed[Style["YA1 varying knl at z="<>ToString@ztoplot, 20,FontFamily\[Rule]"Helvetica"], Background -> White], Scaled[{0.65, 0.2}]]}]*)


End[]
EndPackage[]
