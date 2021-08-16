(* ::Package:: *)

Switch[$nameOfModel,
"Exp II",
gammaref=127.645140`9;
setQuintessenceFixedValues[alphaFixed->7/3,varphiFFixed->-10];
(*gammaref=120.664765`9;
setQuintessenceFixedValues[alphaFixed->7/3,varphiFFixed->-35];*)
(*gammaref=124.257929`9;
setQuintessenceFixedValues[alphaFixed->2/3,varphiFFixed->-10];*)
$stepstring="-AlphaAttractor-ExpII";,
"lcdm",
setQuintessenceFixedValues[alphaFixed->7/3,varphiFFixed->-10];
$stepstring="-lcdm";,
"Exp I",
(*gammaref=126.254076`9;
setQuintessenceFixedValues[alphaFixed->1/3,varphiFFixed->-35];*)
gammaref=127.871242`9;
setQuintessenceFixedValues[alphaFixed->7/3,varphiFFixed->-10];
$stepstring="-AlphaAttractor-ExpI";
];





$justAnalysys=False;
$justDerivatives=False;
$computeDerivatives=True;
$deleteInterpolationFiles=False;


(*fiIO=3; (* interpolation order for input files*)

interpmeth="Spline";  (* interpolation method for input files*)*)

fileinputkbinning="logint50";  (* input CAMB files to be used *)

interpInputTabs="LinLin";     (* interpolation of Tables, either in LinLin or LogLog *)  (*LogLog has a yet unsolved problem for numerical derivs *)

(* $newrecipelinear=True;  nonlinear recipe if False, reduces to linear recipe if True *)

kmaxLinear=0.20;     (*maximum k scale in Fisher *)

$zdependderivatives="numerical";  (*numerical or analytical derivatives of z-dependent variables*)

(* $epsilonstep=0.003;   (* epsilon step for numerical derivatives of shape parameters, this is connected to the +/-epsilon input files *)

epsizstepdenominator=10000;
$epsilonzstep=N@(1/epsizstepdenominator);   epsilon step for numerical derivatives of z-dependent parameters lnDA and lnH for the moment only *)

$stencilpoints=5;    (* number of points for the stencil of the numerical derivative of z-dependent quantities *)

(*$stepstring="-vNL_1";   (* string identifying this specific type of Fisher code runs*)*)

frontendversion=replaceStringPoint[$VersionNumber, "math"];
