(* ::Package:: *)

(* Weak Lensing Fisher Matrix Package *)
(* This package belongs to the FisherMatrixTools package and needs Cosmology.m and FisherMatrix.m  *)

BeginPackage["WeakLensingFisherTools`"]



(**** Global input, numerical, physical or survey parameters  ****)

$ellmin::usage="Minimum multipole ell considered in the WL Fisher Matrix sum."
$ellmax::usage="Maximum multipole ell considered in the WL Fisher Matrix sum."
$deltaell::usage="Spacing in multipole ell when summing over all multipoles in the WL Fisher Matrix sum."
$logellrange::usage="Range (bins) of log10 multipoles (Log10[ell]) where Fisher matrix is evaluated. "
$ellrange::usage="Range (bins) of multipoles where Fisher matrix is evaluated. "
$ellbinscenters::usage="Centers of the bins in multipoles where Fisher matrix is evaluated. "
$galaxyDensArcminSq::usage="Galaxy density per arcmin squared. WL specification."
$numZbinsWL::usage="Number of correlated redshift bins for WL."
$numellbins::usage="Number of multipole bins (ell) for WL Fisher matrix."


$zminSurvey::usage="Minimum redshift for the z-bins of the survey"
$zmaxSurvey::usage="Maximum redshift for the z-bins of the survey"

$zintepoints::usage="Number of zpoints used for interpolation of z-dependent quantities of C_ells"
$photozparameters::usage="Vector of photometric redshift error parameters for the survey."


$photozparametersAnaly::usage="Default: $photozparametersAnaly. Analytical parameter names."
$photozParsRule::usage="Default: $photozParsRule. Rule matching analytical names and numbers."


$numZbinsWL::usage="Number of bins for the tomographic Weak Lensing survey. Default: $numZbinsWL=10."

$customZbinsWL::usage="Default: $customZbinsWL={}. Variable that stores a custom redshift-binning of the WL functions.
Used in case that the survey provides a specific survey or for testing different binnings."

$biasParamsGCph::usage="Names of bias parameters for the GCph probe. One independent parameter per redshift bin."

$zmeanWL::usage="Mean redshift of the galaxy number density distribution. Survey specification parameter."

$phZerror::usage="Photometric redshift error for the specified survey. Only used in the simple Gaussian photo-z distribution."

$intrinsicShear::usage="Sets the intrinsic shear specification for the given survey. Default: $intrinsicShear=0.22"

$areaSurveyWL::usage="Area in sq. deg. covered by the survey. Used to calculate $fsky. Survey specification parameter."

(**** Functions for the Weak Lensing analysis  ****)


setWLBinsSpecifications::usage="setWLBinsSpecifications sets the specifications in multipoles and redshift bins for a WL Fisher Matrix."

setWeakLensingSurveySpecifications::usage="Form: setWeakLensingSurveySpecifications[survey_String].
Set the survey specifications of the provided WL survey string. Hard coded at the moment."

activateIntrinsicAlignments::usage="Option for setWeakLensingSurveySpecifications.
If set to True, the IA term is included in the shear shear Cij."

customWLredshiftBins::usage="Default: customWLredshiftBins->False. Option for setWeakLensingSurveySpecifications. If set to True,
then a custom biinning of redshift bins is read from the variable $customZbinsWL."

nGalaxy::usage="Function of redshift z and mean redshift $z0WL, specifying the distribution of galaxies in redshift"
nGalaxyArea::usage="Integral over the redshift distribution of galaxies."

nGalaxyNormalized::usage="Form: nGalaxyNormalized[z_]. Computes the normalized n(z) galaxy distribution, by dividing
the raw n(z) by the $nGalaxyNormalization."

calculateZbins::usage="Form: calculateZbins[distr,count,a,b].
Takes a distribution 'distr' of number of galaxies as a function of z and returns 'count' number of equally populated redshift bins in the range a<z<b ."

photoZprobabilityPerfect::usage="Form: photoZprobabilityPerfect[z_, bin_?IntegerQ]. Probability distribution of photometric redshifts, for the
case of a perfect determination (Heaviside step function)."

normalizationOfconvoluted$nGalaxyBins$PerfectPhotoZ::usage="Calculates the normalization of the convolution of the galaxy number density and
the perfect photo-z distribution over the range in z of the survey for the bin 'bin'.
Form: normalizationOfconvoluted$nGalaxyBins$PerfectPhotoZ[bin_]"


convoluted$nGalaxyBins$perfectPhotoZ::usage="Form: convoluted$nGalaxyBins$perfectPhotoZ[z_, bin_].
The convolution between the distribution of galaxy number density and the photometric redshift probability in z, for the case
of a perfect photo-z determination."

photoZprobabilityGaussian::usage="Form: photoZprobabilityGaussian[z_,zp_, sigz_],
the functional form of the photo-z probability distribution function."

photozprobGaussianWithErrors::usage="Form: photozprobGaussianWithErrors[z_,zp_, errvect_],
the functional form of the photo-z probability distribution function with double exponential and error parameters.
Accepts an error vector of 7 parameters."

convoluted$nGalaxyBins$Gaussian::usage="Form: convoluted$nGalaxyBins$Gaussian[z, bin, opts].
Calculates the convolution between the photo-z probability distribution and the number density of galaxies.
Accepts a redshift value z, the bin number 'bin' and Options."


convniz::usage="Form: convniz: Analytic integral of convoluted n(z)."
integralniz::usage="Form: integralniz. Analytic integral."


photometricZprobabilityDistribution::usage="Option for convoluted$nGalaxyBins$Gaussian to specify
which photo-z probability distribution will be convoluted with the galaxy number density distribution."

photoZerrorsVector::usage="Option for convoluted$nGalaxyBins$Gaussian to specify the
 photo-z errors for the photo-z distribution. Default: $photozparameters"


galaxyNorm$convoluted$nGalaxyBins$Gaussian::usage="Form: galaxyNorm$convoluted$nGalaxyBins$Gaussian[bin, opts].
Normalization of the convoluted galaxy number density."


convoluted$nGalaxyBins$GaussianNormalized::usage="Form: convoluted$nGalaxyBins$GaussianNormalized[z_, bin_, opts].
Normalized convolution of the galaxy number density with the photometric redshift error distribution."


galaxySurfaceDensity::binmethod="At the moment, only equipopulated tomographic bins are implemented for the n_i."
galaxySurfaceDensity::usage="Number of galaxies per tomographic bin per surface area. For output in sterradians, use default option: dimensionless->False."
tomographicBinMethod::usage="Option for galaxySurfaceDensity. At the moment only the option 'equipopulated' is implemented."


windowA::usage="Form: windowA[bin,opts]. Computes the integral of the convoluted$nGalaxyBins$GaussianNormalized, using NDSolveValue. First term of the window function."
windowB::usage="Form: windowB[bin,opts]. Computes the integral of the convoluted$nGalaxyBins$GaussianNormalized over
the comoving distance, using NDSolveValue. Second term of the window function."

windowAND::usage="Form: windowAND[bin,zi,opts]. Equivalent to windowA[bin,opts][zi]."
windowBND::usage="Form: windowBND[bin,zi,opts]. Equivalent to windowB[bin,opts][zi]."

windowAInt::usage="Form: windowA[bin,opts]. Computes the integral of the convoluted$nGalaxyBins$GaussianNormalized. First term of the window function."
windowBInt::usage="Form: windowB[bin,opts]. Computes the integral of the convoluted$nGalaxyBins$GaussianNormalized over
the comoving distance. Second term of the window function."

windowTildeGG::usage="Form: windowTildeGG[z,bin,opts]. Computes the window function from the two integrated parts windowA and windowB."

windowAFunction::usage="Default: windowAFunction->windowAND. Option for windowTildeGG to
get the function that will compute the window integrals."

windowBFunction::usage="Default: windowBFunction->windowBND. Option for windowTildeGG to
get the function that will compute the window integrals."

windowGammaGamma::usage="Form: windowGammaGamma[z,bin,opts]. Computes the window function for the shear from the windowTildeGG function and the cosmological parameters."

constantBias::usage="Option for biasGalGalParam. If set to True, the bias at each bin is just a constant nuisance parameter."
windowTildeGalGal::usage="The window function for the photometric galaxy clustering."
biasGalGalParam::usage="Function specifying how the bias function for GCph is calculated as a function of paramters."
biasGalGal::usage="Function specifying how the bias is calculated as a function of redshift."

windowIA::usage="The window function for the photometric weak lensing Intrinsic Alignment terms."
windowLL::usage="The window function for the photometric total (shear) weak lensing term."

kofell::usage="Form: kofell[z_,ell,opts]. Computes for a given multipole ell, the corresponding k mode in the Limber approximation, as a function of redshift z."
ellofk::usage="Form: ellofk[z_,ell,opts]. Computes for a given k mode, the corresponding multipole ell in the Limber approximation, as a function of redshift."

pijIntegrandGG::usage="Form: pijIntegrandGG[z,ell,i,j,opts]. The integrand of the Cij function for the shear (gamma,gamma).
Contains the Kernel functions with its window functions and the Limber power spectrum, using the standard window gamma functions.
For GR+LCDM cases."

pijIntegrandGI::usage="Form: pijIntegrandGI[z,ell,i,j,opts]. The integrand of the Cij function for the (gamma,IA) component.
Contains the Kernel functions with its window functions and the Limber power spectrum, using the standard window gamma functions.
For GR+LCDM cases."

pijIntegrandII::usage="Form: pijIntegrandII[z,ell,i,j,opts]. The integrand of the Cij function for the (IA,IA) component.
Contains the Kernel functions with its window functions and the Limber power spectrum, using the standard window gamma functions.
For GR+LCDM cases."

pijIntegrand2::usage="Form: pijIntegrand2[z,ell,i,j,opts]. The integrand of the Cij function for the shear (gamma,gamma).
Contains the Kernel functions with its window functions and the Limber power spectrum. It also contains a Sigma lensing function for non-LCDM cases."

pijIntegrandGalGal::usage="Function that calculates the integrand for the Cls of galaxy-galaxy photometric GC."

curlyFIA::usage="Form: curlyFIA[zz_, opts]. The model for the empirically motivated extension of the nonlinear alignment model. Enters into the power spectra of the
 (gamma, IA) and (IA,IA) terms."

$luminosityFileInterpolation::usage="Variable to store the interpolation function of the scaled mean luminosity density."



lensKernelGG::usage="Form: lensKernelGG[z,ii,jj, opts]. Computes the lensing Kernel function for gamma-gamma in LCDM."
lensKernelGI::usage="Form: lensKernelGI[z,ii,jj, opts]. Computes the lensing Kernel function for gamma-IA in LCDM."
lensKernelII::usage="Form: lensKernelII[z,ii,jj, opts]. Computes the lensing Kernel function for IA-IA in LCDM."

lensKernelGGTab::usage="Table created for interpolating lensKernelGG in redshift."
lensKernelGGInterpol::usage="Function for interpolating lensKernelGG in redshift."
lensKernelGITab::usage="Table created for interpolating lensKernelGI in redshift."
lensKernelGIInterpol::usage="Function for interpolating lensKernelGI in redshift."
lensKernelIITab::usage="Table created for interpolating lensKernelII in redshift."
lensKernelIIInterpol::usage="Function for interpolating lensKernelII in redshift."

initializeLensKernels::usage="Function that initializes the Lensing Kernel interpolating functions for all redshift bins.
Saves time in the computation of Cells and the Fisher Matrix."

pLimberInterpol::usage="function that interpolates the pLimber."
pLimber::usage="Limber approximation power spectrum for ell multipoles, including MG Sigma."
pijIntegrandGGInterpol::usage="Function to interpolate pijIntegrandGG"
pijIntegrandGGInterpolFunc::usage="Interpolating function wrapper (to have z as standard arg) for pijIntegrandGG"
pijIntegrandIIInterpolFunc::usage="Interpolating function wrapper for pijIntegrandII"
pijIntegrandIIInterpol::usage="Interpolating function for pijIntegrandII"
pijIntegrandGIInterpol::usage="function to Interpolate pijIntegrandGI"
pijIntegrandGIInterpolFunc::usage="Interpolating function wrapper for pijIntegrandGI"


lensKernelGalGal::usage="Function that computes the kernel for the galaxy-galaxy term in photometric GC."
galaxygalaxyIntegrand::usage="Option for CijGalGal to choose the integrand for the gal-gal Cij."

(*numericalDerivativeStep::usage="Form: numericalDerivativeStep[index_, epsi_].
Computes the steps and the parameter options for a simple numerical derivative with 3point stencil. Accepts the index of the parameter corresponding to its
position in $paramoptions and the epsilon value epsi."
*)

dlensKernGGdp::usage="Form: dlensKernGGdp[zz,i,j,alpha, opts]. Derivative of the lensing gamma-gamma Kernel function."

dPLimberdp::usage="Form: dPLimberdp[zz,ell,alpha, opts]. Derivative of the Limber power spectrum w.r.t. the cosmological parameters.
Options: setEpsilon->$epsilonstep."


Cij::usage="Form: Cij[ell,i,j,opts]. Cij function for the total lensing spectra, including IA terms.
Accepts $paramoptions as options. "
CijGG::usage="Form: CijGG[ell,i,j,opts]. Cij function for the lensing gamma-gamma component. Accepts $paramoptions as options. "
CijGI::usage="Form: CijGI[ell,i,j,opts]. Cij function for the lensing gamma-IA correlation. Accepts $paramoptions as options. "
CijII::usage="Form: CijII[ell,i,j,opts]. Cij function for the lensing IA-IA component. Accepts $paramoptions as options. "
CijGalGal::usage="Function that calculates the Cls of galaxy-galaxy photometric GC."

shearshearIntegrand::usage="Option for Cij, specifying the integrand of the shear-shear Cij tomographic spectra."
shearIAIntegrand::usage="Option for Cij, specifying the integrand of the shear-IA Cij tomographic spectra."
IAIAIntegrand::usage="Option for Cij, specifying the integrand of the IA-IA Cij tomographic spectra."

integrationOptions::usage="Default: integrationOptions->{AccuracyGoal->8,
  Method->{'GlobalAdaptive','MaxErrorIncreases'->10000, 'SymbolicProcessing'->False}}.
Option for Cij and other functions that use NIntegrate as integration method. For full automatic integration, use integrationOptions->{}."

dCijdp::usage="Form: dCijdp[ell,i,j,alpha]. Computes the derivative of the Cij function w.r.t. the cosmological parameters in $paramoptions, using
an epsilon step of $epsilonstep."

dCijdpIntegral::usage="Form: dCijdpIntegral[ell,i,j,alpha]. Computes the integral in redshift (from $zminSurvey to $zmaxSurvey) over
the derivative of the integrand of Cij (pijIntegrand2) w.r.t. the cosmological parameters sepcified in $paramoptions.
It uses the lensing Kenrel gamma-gamma function and the Limber power spectrum."

noiseMatrix::usage="Form: noiseMatrix[i_,j_]. Specifies the noise term that goes into the lensing Covariance matrix."

kDamping::usage="Option for kScaleErrorDamping and FisherWL calculations. After this scale in k, the covariance matrix is exponentially damped."
tolerance::usage="Option for kScaleErrorDamping. If the damping factor is smaller than 'tolerance', it is set to zero, to avoid numerical errors."
dampingExponent::usage="Option for kScaleErrorDamping. Exponent 'n' inside the exponential function for a damping in k: Exp[-k(ell,z)^n/kDamping^n]."
kScaleErrorDamping::usage="Function 'kScaleErrorDamping[ell, i, opts]' to damp exponentially the covariance matrix after a specific scale in 'kDamping' as a function of redshift bin 'i' and multipole 'ell'. Options: {kDamping, dampingExponent, tolerance}."
kScaleErrorDampingMatrix::usage="Matrix fof the kScaleErrorDamping function for all redshift bins for a WL observable. Options: {kDamping}."


inversecovariance::usage="Form: inversecovariance[ell]. Computes the inverse of the covariance Matrix for each multipole ell.
Options: cijFunction->Cij : The function that computes the Cij's.
         numberOfZBins->numbins : The number of bins of the lensing survey, that correspond to the dimensions of the covariance matrix."
cijFunction::usage="Option for inverseCovariance. Default: cijFunction->Cij. Specifies which function to use for computing the Cij's.
Can be used to specify external input Cij functions."
numberOfZBins::usage="Option for Fisher WL functions.
It specifies the number of bins of the lensing survey, that correspond to the dimensions of the covariance matrix. Default: numberOfZBins->$numZbinsWL"

numberOfellBins::usage="Option for Fisher WL functions.
It specifies the number of bins in multipole ell of the Fisher matrix sum. Default: numberOfellBins->$numellbins"


Jmat::usage="Form: Jmat[ell,alpha]. Matrix with components 'ij' of derivatives of the Cij function w.r.t. parameter alpha.
Options: numberOfZBins->numbins : The number of bins of the lensing survey, that correspond to the dimensions of the covariance matrix.
         dcijdpFunction->dCijdpIntegral : This option specifies the function that computes the derivatives of the Cij's. It can be used to provide
external interpolated input  of Cij derivatives."
dcijdpFunction::usage="Option for WL Fisher functions. It specifies which function should be used to compute the derivative
of the Cij's. It can be used to provide external interpolated input  of Cij derivatives. Default: dcijdpFunction->dCijdpIntegral."

FisherWL::usage="Form: FisherWL[alpha,beta,opts]. Function that computes the WL Fisher matrix components alpha and beta.
Options:
ellmin->$ellmin : Minimum multipole ell to sum over.
ellmax->$ellmax : Maximum multipole ell to sum over.
Deltaell->$deltaell : Logarithmic (base10) spacing for the sum over ell.
cijFunction->Cij : This option specifies the function that computes the Cij's.
dcijdpFunction->dCijdpIntegral : This option specifies the function that computes the derivatives of the Cij's. It can be used to provide
external interpolated input  of Cij derivatives."

FisherWLSimple::usage="A more efficient form of FisherWL. Check documentation for FisherWL."


FisherWLMatrix::usage="FisherWLMatrix[opts]. Function that evaluates FisherWLSimple into a Matrix
by computing only the symmetric elements. Options: keepParameters."
keepParameters::usage="Option for FisherWLMatrix.
Default: 'All'
If set to a List of integers, the Fisher Matrix will be computed only for the parameters
at the positions in $paramnames specified by the list."

falseCij::usage="Form: falseCij. Function with Pause to test parallelization."

ellmin::usage="Option for Fisher WL functions. Minimum multipole ell to sum over. Default: ellmin->$ellmin"
ellmax::usage="Option for Fisher WL functions. Maximum multipole ell to sum over. Default: ellmax->$ellmax"
Deltaell::usage="Option for Fisher WL functions. Default: Deltaell->$deltaell. Logarithmic (base10) spacing for the sum over ell."
minZinSurvey::usage="Option for Fisher WL functions. Minimum redshift z of the tomographic survey. Default: minZinSurvey->$zminSurvey."
maxZinSurvey::usage="Option for Fisher WL functions. Maximum redshift z of the tomographic survey. Default: maxZinSurvey->$zmaxSurvey."

(**** Global variables calculated at run-time  ****)

$zbinsEquiPopu::usage="List with the redshift bin limits of the equally populated bins in the range $zminSurvey<z<$zmaxSurvey.
Number of bins: $numZbinsWL."
$nGalaxyNormalization::usage="Normalization factor of the galaxy number density redshift distribution.
Computed as: $galaxyDensArcminSq/nGalaxyArea[$zmaxSurvey]."
$z0WL::usage="Parameter used in WL functions which is equal to the mean redshift $zmeanWL of the n(z) distribution, divided by Sqrt[2]."
$fsky::usage="Fraction of the sky covered by the survey. Calculated from $areaSurveyWL."
$alphanz::usage="Exponent in the perfect photometric n(z) function."

(*IST baselinewl notation*)
$aIA::usage="Protected name of parameter A_IA for intrinsic alignment."
$eIA::usage="Protected name of parameter eta_IA for intrinsic alignment."
$bIA::usage="Protected name of parameter beta_IA for intrinsic alignment."
$cIA::usage="Name of parameter C_IA for intrinsic alignment."
(*IST baselinewl notation*)
$IAswitch::usage="Switch that sets if Intrinsic Alignment is present."

$aIAfidu::usage="Fiducial of parameter A_IA for intrinsic alignment."
$eIAfidu::usage="Fiducial of parameter eta_IA for intrinsic alignment."
$bIAfidu::usage="Fiducial of parameter beta_IA for intrinsic alignment."

clearFisherWLMemoizedQuantities::usage="Function that clears the memoized values of WL functions. Clears in main and parallel Kernels."


Needs["FisherTools`"];
Needs["CosmologyFunctions`"];
Needs["UsefulTools`"];
Needs["FunctionApproximations`"];
EndPackage[]

BeginPackage["WeakLensingFisherTools`", {"FisherTools`","CosmologyFunctions`","UsefulTools`","FunctionApproximations`"}]

Protect[$aIA, $eIA, $bIA];

Begin["`Private`"]




$photozparameters={1.0,0.0,0.05,1.0,0.1,0.05,0.10};
$photozparametersAnaly={cb,zb,sb,c0,z0,s0,fo};
$photozParsRule=Thread[Rule[$photozparametersAnaly,$photozparameters]];
(*surveyString=chosenSurvey;*)
$fsky:=$areaSurveyWL/(4 Pi*(180/Pi)^2);
$deltaell=0.01;(*logarithm Log10 spacing between $ellmin and $ellmax*)
$ellmin=10;
$ellmax=5000;
$numellbins=100;
(*ellmax0nonlin=5000;(*kmax=0.2;*)
ellmax0lin=1000;(*kmax=0.2;*) *)
$zminSurvey=0.01;
$zmaxSurvey=3.0;
$numZbinsWL=10;
$customZbinsWL={};

$IAswitch=0;
$aIAfidu=1.72;
$bIAfidu=2.17;
$eIAfidu=-0.41;

$zintepoints=300;
zmaxInt=$zmaxSurvey+0.1; (*for limits in integrals*) (* end of $zrangeGlobal  *)
zminInt=10^-6; (*for limits in integrals*)


Options[setWLBinsSpecifications]={ellmin:>$ellmin, ellmax:>$ellmax, Deltaell:>$deltaell, minZinSurvey:>$zminSurvey,
  maxZinSurvey:>$zmaxSurvey, numberOfZBins:>$numZbinsWL, numberOfellBins:>$numellbins}

setWLBinsSpecifications[opts:OptionsPattern[]]:=Block[{checkvect, delell, ellbins, aia},

  $ellmin=OptionValue[ellmin];
  $ellmax=OptionValue[ellmax];
  delell=OptionValue[Deltaell];
  ellbins=OptionValue[numberOfellBins];
  Which[NumberQ[ellbins]==True,
    $numellbins=ellbins;
    $deltaell=N@((Log[10,$ellmax]-Log[10,$ellmin])/$numellbins);
    ,
    NumberQ[delell]==True,
    $deltaell=delell
  ];
  $logellrange=Range[Log10[$ellmin], Log10[$ellmax], $deltaell];
  $ellbinscenters=tentox@(MovingAverage[$logellrange,2]);
  $ellrange=tentox@$logellrange;
  $numZbinsWL=OptionValue[numberOfZBins];
  $zminSurvey=OptionValue[minZinSurvey];
  $zmaxSurvey=OptionValue[maxZinSurvey];
]

Options[setWeakLensingSurveySpecifications]={customWLredshiftBins->False, activateIntrinsicAlignments->True}

setWeakLensingSurveySpecifications::binNum="Number of bins in $numZbinsWL does not match size of $zbinsEquiPopu"

setWeakLensingSurveySpecifications[survey_String, opts:OptionsPattern[]]:=Block[{sets,
  custwlzbins=OptionValue[customWLredshiftBins],
  boolIA=OptionValue[activateIntrinsicAlignments]},
  chosenSurvey=survey;
  Which[StringMatchQ[chosenSurvey,"*Euclid*"],
    $zmeanWL=0.9;(*z_mean over median*)
    $z0WL=$zmeanWL/Sqrt[2];
    $z0pre=$z0WL;
    $alphanz=3/2;
    $phZerror=.05;(*photometric redshift error*)
    $intrinsicShear=.22;
    $cIA=0.0134;  (*IA nonlinear algnment model*)
    $galaxyDensArcminSq=30;(*euclid spec*)
    $areaSurveyWL=15000;
    $customZbinsWL={0.001, 0.418, 0.560, 0.678, 0.789, 0.900 , 1.019, 1.155, 1.324, 1.576, 2.50};
    Print[styleFunction[chosenSurvey<>" Specifications loaded ****",25]];
    ,
    StringMatchQ[chosenSurvey,"*SKA2*"],
    $zmeanWL=1.6;(*z_mean over median*)
    $z0WL=$zmeanWL/Sqrt[2];
    $z0pre=$z0WL;
    $alphanz=3/2;
    $phZerror=.05;(*photometric redshift error*)
    $cIA=0.0134;  (*IA nonlinear algnment model*)
    $intrinsicShear=0.3;
    $galaxyDensArcminSq=10;(*SKA2 spec*)
    $areaSurveyWL=30940;
    Print[styleFunction[chosenSurvey<>" Specifications loaded ****",25]
    ];
    ,
    StringMatchQ[chosenSurvey,"*SKA1*"],
  (*SKA1 advanced,not early one*)
    $zmeanWL=1.0;(*z_mean over median*)
    $z0WL=$zmeanWL/Sqrt[2];
    $z0pre=$z0WL;
    $alphanz=3/2;
    $phZerror=.05;(*photometric redshift error*)
    $intrinsicShear=0.3;
    $cIA=0.0134;  (*IA nonlinear algnment model*)
    $galaxyDensArcminSq=2.7;(*SKA1 spec*)
    $areaSurveyWL=5000;
    Print[styleFunction[chosenSurvey<>" Specifications loaded ****",25]];
    ,
    StringMatchQ[chosenSurvey,"*LSST*"],
    $zmeanWL=1.03;(*z_mean over median*)
    $z0WL=0.11;
    $z0pre=1;
    $alphanz=0.68;
    $phZerror=.05;(*photometric redshift error*)
    $cIA=0.0134;  (*IA nonlinear algnment model*)
    $intrinsicShear=0.26;
    $galaxyDensArcminSq=27;
    $areaSurveyWL=14300;
    Print[styleFunction[chosenSurvey<>" Specifications loaded ****",25]
    ];
  ];
  integralniz=Integrate[photozprobGaussianWithErrors[za,zap,$photozparametersAnaly]*nGalaxy[za,z00],{zap,zill,zull}];
  If[
    SameQ[custwlzbins, True] && VectorQ[$customZbinsWL, NumericQ],
    $zbinsEquiPopu=$customZbinsWL;
    ,
    $zbinsEquiPopu=calculateZbins[nGalaxy[#, $z0WL]&,$numZbinsWL,$zminSurvey,$zmaxSurvey];
    ];
  (*Consistency check*)
     If[$numZbinsWL != Length[$zbinsEquiPopu]-1,
       Message[setWeakLensingSurveySpecifications::binNum];
       Abort[] ];

  $biasParamsGCph = (ToExpression["b"<>ToString[#]]&/@Range[$numZbinsWL]);

  If[boolIA==True,
     $IAswitch=1;
     ,
     $IAswitch=0
  ];
 ]





nGalaxy[z_, z0_] := (z/$z0pre)^2*Exp[-(z/z0)^($alphanz)];

(*nGalaxy[z_]:=1.50*Exp[-((z-0.70)/0.32)^2]+0.20*Exp[-((z-1.20)/0.46)^2]  (*CFHTLens specs*) *)


nGalaxyAntidev[zmi_,zma_]:=nGalaxyAntidev[zmi,zma]=Block[{zzz},
NDSolve[{y'[zzz]==nGalaxy[zzz, $z0WL],y[0]==0},y,{zzz,zmi,zma}][[1,1,2]]
];

nGalaxyArea[z_?NumericQ]:=nGalaxyAntidev[$zminSurvey,$zmaxSurvey][z];

nGalaxyNormalized[z_]:=With[{ngalnormalization=$galaxyDensArcminSq/nGalaxyArea[$zmaxSurvey]},
  $nGalaxyNormalization=ngalnormalization; ngalnormalization*nGalaxy[z, $z0WL] ]


(*Compute z-bins such that they divide redshift space into intervals with equal area under the nGalaxy curve*)
calculateZbins[distr_,count_,a_?NumericQ,b_?NumericQ] :=calculateZbins[distr,count,a,b]=Module[{bins={a},area,tt},
  integral[x1_?NumericQ,x2_?NumericQ] := NIntegrate[distr[x],{x,x1,x2}];
  area = integral[a,b]/count;
  Do[
    tt = FindRoot[integral[bins[[i]],xx]==area,{xx,bins[[i]]*1.1+.1}][[1,2]];
    bins=Append[bins,tt],{i,1,count}];
  bins];


photoZprobabilityPerfect[z_, bin_?IntegerQ]:=Block[{zl=$zbinsEquiPopu[[bin]], zu=$zbinsEquiPopu[[bin+1]]},HeavisideTheta[z-zl]HeavisideTheta[zu-z]]

normalizationOfconvoluted$nGalaxyBins$PerfectPhotoZ[bin_]:=normalizationOfconvoluted$nGalaxyBins$PerfectPhotoZ[bin]=Block[{zi},
  Integrate[photoZprobabilityPerfect[zi,1]nGalaxyNormalized[zi], {zi,$zminSurvey,$zmaxSurvey}]]


convoluted$nGalaxyBins$perfectPhotoZ[z_, bin_]:=
    convoluted$nGalaxyBins$perfectPhotoZ[z, bin]=
        photoZprobabilityPerfect[z,bin]nGalaxyNormalized[z]/(normalizationOfconvoluted$nGalaxyBins$PerfectPhotoZ[bin])


photoZprobabilityGaussian[z_,zp_, sigz_]:=Exp[-(zp-z)^2/(2 sigz)]

photozprobGaussianWithErrors[z_,zp_, errvect_]:=With[{cb=errvect[[1]],zb=errvect[[2]],
  sigmab=errvect[[3]], c0=errvect[[4]], z0=errvect[[5]], sigma0=errvect[[6]],fout=errvect[[7]]},
  ((1-fout)/(Sqrt[2 Pi]*sigmab*(1+z)))*Exp[-(1/2)((z-cb zp-zb)/(sigmab(1+z)))^2]+
      ((fout)/(Sqrt[2 Pi]*sigma0*(1+z)))*Exp[-(1/2)((z-c0 zp-z0)/(sigma0(1+z)))^2]]



(*convniz[zva_,zvi_,zvu_]:=Block[{zw},(integralniz/.($photozParsRule~Join~{z00->$z0WL}))/.{Global`za->zva, Global`zill->zvi, Global`zull->zvu}];*)
convniz[zva_,zvi_,zvu_]:=Block[{zw},(integralniz/.($photozParsRule~Join~{z00->$z0WL}))/.{za->zva, zill->zvi, zull->zvu}];


Options[convoluted$nGalaxyBins$Gaussian]={photometricZprobabilityDistribution->photozprobGaussianWithErrors, photoZerrorsVector->$photozparameters};

convoluted$nGalaxyBins$Gaussian[z_?NumericQ, bin_, opts:OptionsPattern[]]:=
    convoluted$nGalaxyBins$Gaussian[z, bin, opts]=Block[{photoZpdf, photoerrvec,zl,zu},
  photoZpdf=OptionValue[photometricZprobabilityDistribution];
  zl=$zbinsEquiPopu[[bin]];
  zu=$zbinsEquiPopu[[bin+1]];
  (*photoerrvec=OptionValue[photoZerrorsVector];
  NIntegrate[photoZpdf[z,zp, photoerrvec]nGalaxy[zp,$z0WL],{zp,zl,zu}]*)
   convniz[z,zl,zu]
    ]


galaxyNorm$convoluted$nGalaxyBins$Gaussian[bin_, opts:OptionsPattern[{convoluted$nGalaxyBins$Gaussian}]]:=galaxyNorm$convoluted$nGalaxyBins$Gaussian[bin, opts]=Block[{zl},
  NIntegrate[convoluted$nGalaxyBins$Gaussian[zz, bin, opts],{zz,$zminSurvey,$zmaxSurvey}, AccuracyGoal->5,
  Method->{"GlobalAdaptive",(*"MaxErrorIncreases"\[Rule]10000,*)"SymbolicProcessing"->False}]]


convoluted$nGalaxyBins$GaussianNormalized[z_, bin_, opts:OptionsPattern[{convoluted$nGalaxyBins$Gaussian}]]:=
    convoluted$nGalaxyBins$Gaussian[z, bin, opts]/galaxyNorm$convoluted$nGalaxyBins$Gaussian[bin, opts]



Options[galaxySurfaceDensity]={dimensionless->False, tomographicBinMethod->"equipopulated"}
galaxySurfaceDensity[bin_]:=Block[{nis, sterradians},
  sterradians = 3600*(180/Pi)^2;
  If[OptionValue[dimensionless]==True,
    sterradians=1;
  ];
  If[OptionValue[tomographicBinMethod]!="equipopulated",
    Message[galaxySurfaceDensity::binmethod];
  ];
  nis=(($galaxyDensArcminSq/$numZbinsWL)*sterradians);
  Return[nis]
]

windowA[bin_,opts:OptionsPattern[$paramoptions]]:=windowA[bin,opts]=Block[{zzz},
    NDSolveValue[{convoluted$nGalaxyBins$GaussianNormalized[zzz,bin]==f'[zzz],f[zmaxInt]==0},f,{zzz,zminInt,zmaxInt},MaxStepSize->.002]
    ];

windowAND[bin_,zi_?NumericQ, opts:OptionsPattern[$paramoptions]]:=(-windowA[bin,opts][zi]);(*minus sign comes from integration limits*)

windowAInt[bin_,zi_?NumericQ, opts:OptionsPattern[$paramoptions]]:=windowAInt[bin,zi,opts]=
    Block[{zzz}, NIntegrate[convoluted$nGalaxyBins$GaussianNormalized[zzz,bin],{zz,zi,$zmaxSurvey}]];


windowB[bin_,opts:OptionsPattern[$paramoptions]]:=windowB[bin,opts]=Block[{zzz},
    NDSolveValue[{convoluted$nGalaxyBins$GaussianNormalized[zzz,bin]/(comovingDistance[zzz, {opts}~Join~{dimensionless->True}])==f'[zzz],
      f[zmaxInt]==0 },f,{zzz,zminInt,zmaxInt},MaxStepSize->.002]];

windowBND[bin_,zi_?NumericQ, opts:OptionsPattern[$paramoptions]]:=(-windowB[bin,opts][zi]);(*minus sign comes from integration limits*)

windowBInt[bin_,zi_?NumericQ, opts:OptionsPattern[$paramoptions]]:=windowBInt[bin,zi,opts]=Block[{zz},
    NIntegrate[convoluted$nGalaxyBins$GaussianNormalized[zz,bin]/(comovingDistance[zz, {opts}~Join~{dimensionless->True}]),
      {zz,zi,$zmaxSurvey}]];


Options[windowTildeGG]=$paramoptions~Join~{windowAFunction->windowAND, windowBFunction->windowBND};

windowTildeGG[z_,bin_,opts:OptionsPattern[]]:=windowTildeGG[z,bin,opts]=Block[{windA=OptionValue[windowAFunction], windB=OptionValue[windowBFunction],
   cosmopts},
    cosmopts=complementParamValues[{opts},windowTildeGG,returnList->"Full",filterCosmoPars->True];
    (windA[bin,z,cosmopts]-comovingDistance[z, {cosmopts}~Join~{dimensionless->True}]*windB[bin,z,cosmopts])
];

windowGammaGamma[z_,bin_,opts:OptionsPattern[$paramoptions]]:=windowGammaGamma[z,bin,opts]=Block[{compunts,units,hval},
      compunts=compatibleUnits[$internalDistanceUnits, invertUnits->True];
      hval=hubbleToday[opts];
      units=externalUnitsConverter[compunts,hval, physicalUnits->False, inverseDistanceUnits->True, internalUnits->True];
      (units)^2*(3/2)*(1+z)*OmegaM[0.,opts]*comovingDistance[z,{opts}~Join~{dimensionless->False}]*windowTildeGG[z,bin,opts]
];


Options[windowTildeGalGal]=$paramoptions~Join~{};

windowTildeGalGal[z_,bin_,opts:OptionsPattern[]]:=windowTildeGalGal[z,bin,opts]=Block[{cosmopts, hval, units, compunts, ret},
  compunts=compatibleUnits[$internalDistanceUnits, invertUnits->True];
  hval=hubbleToday[opts];
  units=externalUnitsConverter[compunts,hval, physicalUnits->False, inverseDistanceUnits->True, internalUnits->True];
  cosmopts=complementParamValues[{opts},windowTildeGalGal,returnList->"Full",filterCosmoPars->True];
  ret = biasGalGalParam[z,bin, cosmopts] * (convoluted$nGalaxyBins$GaussianNormalized[z, bin]) * DimensionlessHubble[z,cosmopts];

  Return[ret]
]

biasGalGal[z_,bin_]:=Block[{par, zbar},
  zbar = Mean[{$zbinsEquiPopu[[bin]], $zbinsEquiPopu[[bin+1]]}];
  Return[Sqrt[1+zbar]]
];


Options[biasGalGalParam]=$paramoptions~Join~{constantBias->True};
biasGalGalParam[z_,bin_, opts:OptionsPattern[]]:=Block[{par, parname, parval, parfidu, fidus, biopts, fixval, bivalue, tolerance=10^-3},
  fidus=complementParamValues[{opts},biasGalGalParam,returnList->"Fiducials"];
  biopts=complementParamValues[{opts},biasGalGalParam,returnList->"Complement"];
  fixval = biasGalGal[z,bin];
  par = $biasParamsGCph[[bin]];
  parfidu = par/.fidus;
  debugPrint[biopts, 1];
  debugPrint[par, 1];
  If[biopts=={}, (* Default parameters passed, return bias function*)
    debugPrint["default", 1];
    Return[fixval]
  ];
  If[MemberQ[$biasParamsGCph, biopts[[1]][[1]]  ]==False, (* non-bias parameter passed, return  bias function*)
    debugPrint["non bias param given", 1];
    Return[fixval]
  ];
  If[SameQ[biopts[[1]][[1]],par] == True,  (*Correct bias parameter passed for the current bin, return value of bias param.*)
    debugPrint["fiducial given", 1];
    Return[biopts[[1]][[2]]]
     , (*Else*)
    debugPrint["non matching bias param given", 1];  (*Incorrect bias parameter passed for the current bin, return 0.*)
    Return[0.]
  ]
  (* Print["Other bias models not implemented yet."];*)
 ]


Options[windowIA]=$paramoptions;

windowIA[z_,ii_, opts:OptionsPattern[]]:=windowIA[z,ii, opts]=
    Block[{compunts, units, hval},
      compunts=compatibleUnits[$internalPkUnits, invertUnits->False];
      hval=hubbleToday[opts];
      units=externalUnitsConverter[compunts,hval, physicalUnits->False, inverseDistanceUnits->True, internalUnits->True];
      (units) * (DimensionlessHubble[z,opts])*
          (convoluted$nGalaxyBins$GaussianNormalized[z,ii])
    ];


Options[windowLL]=$paramoptions;

windowLL[z_,ii_, opts:OptionsPattern[]]:=windowLL[z,ii, opts]=Block[{ret, OmegaM0},
      OmegaM0=OmegaM0Today[opts];
      ret = windowGammaGamma[z,ii] + (OmegaM0/Growth[z,opts])*(curlyFIA[z, opts])*windowIA[z,ii,opts];
      Return[ret]
]


kofell[z_?NumericQ,ell_?NumericQ,opts:OptionsPattern[$paramoptions]]:=(ell+0.5)/(comovingDistance[z,{opts}~Join~{dimensionless->False}])
(*to fix: need to pass custom units to comovingDistance, such that it is compatible with the units of k*)
ellofk[z_?NumericQ,k_?NumericQ,opts:OptionsPattern[$paramoptions]]:=(k*(comovingDistance[z,{opts}~Join~{dimensionless->False}])-0.5)


Options[lensKernelGG]=$paramoptions~Join~{dimensionless->True};

lensKernelGG[z_,ii_,jj_, opts:OptionsPattern[]]:=lensKernelGG[z,ii,jj, opts]=If[ii<jj, lensKernelGG[z,jj,ii, opts],
    Block[{compunts, units, hval, dimless},
      dimless = OptionValue[dimensionless];
      compunts=compatibleUnits[$internalPkUnits, invertUnits->False];
      hval=hubbleToday[opts];
      units=If[dimless,
      1,
      externalUnitsConverter[compunts,hval, physicalUnits->False, inverseDistanceUnits->True, internalUnits->True];
      ];
      (units)^3* ((3/2)*(1+z)*OmegaM0Today[opts])^2 * (1/DimensionlessHubble[z,opts])*windowTildeGG[z,ii,opts]*windowTildeGG[z,jj,opts]
]
];


Options[lensKernelGI]=$paramoptions~Join~{dimensionless->True};

lensKernelGI[z_,ii_,jj_, opts:OptionsPattern[]]:=lensKernelGI[z,ii,jj, opts]=If[ii<jj, lensKernelGI[z,jj,ii, opts],
    Block[{compunts, units, hval, dimless},
      dimless = OptionValue[dimensionless];
      compunts=compatibleUnits[$internalPkUnits, invertUnits->False];
      hval=hubbleToday[opts];
      units=If[dimless,
      1,
      externalUnitsConverter[compunts,hval, physicalUnits->False, inverseDistanceUnits->True, internalUnits->True];
      ];
      (units)^3* ((3/2)*(1+z)*OmegaM0Today[opts]) * (1/comovingDistance[z, {opts}~Join~{dimensionless->True}])*
          (convoluted$nGalaxyBins$GaussianNormalized[z,ii]*windowTildeGG[z,jj,opts] + convoluted$nGalaxyBins$GaussianNormalized[z,jj]*windowTildeGG[z,ii,opts])
]
];


Options[lensKernelII]=$paramoptions~Join~{dimensionless->True};

lensKernelII[z_,ii_,jj_, opts:OptionsPattern[]]:=lensKernelII[z,ii,jj, opts]=If[ii<jj, lensKernelII[z,jj,ii, opts],
    Block[{compunts, units, hval, dimless},
      dimless = OptionValue[dimensionless];
      compunts=compatibleUnits[$internalPkUnits, invertUnits->False];
      hval=hubbleToday[opts];
      units=If[dimless,
      1,
      externalUnitsConverter[compunts,hval, physicalUnits->False, inverseDistanceUnits->True, internalUnits->True];
      ];
      (units)^3 * (DimensionlessHubble[z,opts]/(comovingDistance[z, {opts}~Join~{dimensionless->True}])^2)*
          (convoluted$nGalaxyBins$GaussianNormalized[z,ii] * convoluted$nGalaxyBins$GaussianNormalized[z,jj])
]
];


Options[lensKernelGalGal]=$paramoptions~Join~{dimensionless->True};

lensKernelGalGal[z_,ii_,jj_, opts:OptionsPattern[]]:=lensKernelGalGal[z,ii,jj, opts]=If[ii<jj, lensKernelGalGal[z,jj,ii, opts],
    Block[{compunts, units, hval, dimless},
      dimless = OptionValue[dimensionless];
      compunts=compatibleUnits[$internalPkUnits, invertUnits->False];
      hval=hubbleToday[opts];
      units=If[dimless,
      1,
      externalUnitsConverter[compunts,hval, physicalUnits->False, inverseDistanceUnits->True, internalUnits->True];
      ];
      (units)^3 * (1/(DimensionlessHubble[z,opts]*comovingDistance[z, {opts}~Join~{dimensionless->True}]^2 ) ) * windowTildeGalGal[z,ii,opts]*windowTildeGalGal[z,jj,opts]
]
];


Options[lensKernelGGTab]=$paramoptions~Join~{dimensionless->True};
Options[lensKernelGGInterpol]=$paramoptions~Join~{dimensionless->True};
lensKernelGGTab[i_,j_, opts:OptionsPattern[]]:=lensKernelGGTab[i,j,opts]=Block[{zz},
Table[{zz,lensKernelGG[zz,i,j, opts]},{zz, $zminSurvey,$zmaxSurvey,($zmaxSurvey-$zminSurvey)/$zintepoints}]];
lensKernelGGInterpol[i_,j_, opts:OptionsPattern[]]:=lensKernelGGInterpol[i,j,opts]=Interpolation[lensKernelGGTab[i,j,opts]]

Options[lensKernelGITab]=$paramoptions~Join~{dimensionless->True};
Options[lensKernelGIInterpol]=$paramoptions~Join~{dimensionless->True};
lensKernelGITab[i_,j_, opts:OptionsPattern[]]:=lensKernelGITab[i,j,opts]=Block[{zz},
Table[{zz,lensKernelGI[zz,i,j, opts]},{zz, $zminSurvey,$zmaxSurvey,($zmaxSurvey-$zminSurvey)/$zintepoints}]];
lensKernelGIInterpol[i_,j_, opts:OptionsPattern[]]:=lensKernelGIInterpol[i,j,opts]=Interpolation[lensKernelGITab[i,j,opts]]

Options[lensKernelIITab]=$paramoptions~Join~{dimensionless->True};
Options[lensKernelIIInterpol]=$paramoptions~Join~{dimensionless->True};
lensKernelIITab[i_,j_, opts:OptionsPattern[]]:=lensKernelIITab[i,j,opts]=Block[{zz},
Table[{zz,lensKernelII[zz,i,j, opts]},{zz, $zminSurvey,$zmaxSurvey,($zmaxSurvey-$zminSurvey)/$zintepoints}]];
lensKernelIIInterpol[i_,j_, opts:OptionsPattern[]]:=lensKernelIIInterpol[i,j,opts]=Interpolation[lensKernelIITab[i,j,opts]]

(*Options[pijIntegrand2]=$paramoptions~Join~{kdamping->0.5};

pijIntegrand2[z_?NumericQ,ell_?NumericQ,i_,j_,opts:OptionsPattern[]]:=pijIntegrand2[z,ell,i,j,opts]= Block[{kda=OptionValue[kdamping], opts},
  opts=complementParamValues[{opts},pijIntegrand2,returnList->"Full",filterCosmoPars->True];
  ($H0^3 /(1+z)^6) (3.0/2.0)^2 windowTildeGG[z,i,opts]*windowTildeGG[z,j,opts] OmegaM[z,opts]^2  DimensionlessHubble[z,opts]^3  *
      (SigmaFunction[z, kofell[z,ell,opts], opts])^2*powerSpectrum[z,kofell[z,ell,opts],opts]
];
*)

initializeLensKernels[]:=Block[{iniT, finiT, tt},
Print["Initializing Lensing Kernels for all redshift bins..."];
iniT=DateString[];
Print["...lensKernelGG..."]
Table[lensKernelGGInterpol[b1,b2,{}],{b1,1,$numZbinsWL},{b2,1,$numZbinsWL}];
Table[lensKernelGGInterpol[b1,b2],{b1,1,$numZbinsWL},{b2,1,$numZbinsWL}];
If[$IAswitch==1,
Print["...lensKernelGI..."]
Table[lensKernelGIInterpol[b1,b2,{}],{b1,1,$numZbinsWL},{b2,1,$numZbinsWL}];
Table[lensKernelGIInterpol[b1,b2],{b1,1,$numZbinsWL},{b2,1,$numZbinsWL}];
Print["...lensKernelII..."]
Table[lensKernelIIInterpol[b1,b2,{}],{b1,1,$numZbinsWL},{b2,1,$numZbinsWL}];
Table[lensKernelIIInterpol[b1,b2],{b1,1,$numZbinsWL},{b2,1,$numZbinsWL}];
];
finiT=DateString[];
tt=timeElapsed[iniT,finiT, False];
Print["Time spent in seconds: ", tt];
];



Options[pLimber]=$paramoptions;
pLimber[zz_, ell_, opts:OptionsPattern[]]:=pLimber[zz,ell,opts]=Block[{kofl},
kofl=kofell[zz,ell,opts];
powerSpectrum[zz,kofl,opts]*(SigmaLensing[zz,kofl,opts])^2
];

Options[pLimberInterpol]=$paramoptions;
pLimberInterpol[ell_, opts:OptionsPattern[]]:=pLimberInterpol[ell,opts]=Block[{zz},
Interpolation[Table[{zz,
pLimber[zz,ell,opts]}, {zz, $zminSurvey,$zmaxSurvey,($zmaxSurvey-$zminSurvey)/$zintepoints}]]];

Options[pijIntegrandGG]=$paramoptions~Join~{dimensionless->True};
pijIntegrandGG[z_?NumericQ,ell_?NumericQ,i_,j_,opts:OptionsPattern[]]:=
    pijIntegrandGG[z,ell,i,j,opts]=Block[{ret},
        lensKernelGG[z,i,j,opts]*pLimber[z, ell, opts]
  ];

Options[pijIntegrandGGInterpol]=$paramoptions~Join~{dimensionless->True};
pijIntegrandGGInterpol[ell_,i_,j_,opts:OptionsPattern[]]:=pijIntegrandGGInterpol[ell,i,j,opts]=Block[{zz},
Interpolation[Table[{zz,
lensKernelGGInterpol[i,j,opts][zz]*pLimberInterpol[ell,opts][zz]},
{zz, $zminSurvey,$zmaxSurvey,($zmaxSurvey-$zminSurvey)/$zintepoints}]]];

Options[pijIntegrandGGInterpolFunc]=$paramoptions~Join~{dimensionless->True};
pijIntegrandGGInterpolFunc[z_,ell_,i_,j_,opts:OptionsPattern[]]:=pijIntegrandGGInterpol[ell,i,j,opts][z]


Options[pijIntegrandGI]=$paramoptions~Join~{dimensionless->True};
pijIntegrandGI[z_?NumericQ,ell_?NumericQ,i_,j_,opts:OptionsPattern[]]:=
    pijIntegrandGI[z,ell,i,j,opts]=Block[{ret, OmegaM0},
      OmegaM0=OmegaM0Today[opts];
      (OmegaM0/Growth[z,kofell[z,ell,opts],opts])*curlyFIA[z, opts]*lensKernelGI[z,i,j,opts]*pLimber[z,ell,opts]
    ];

Options[DIAInterpol]=$paramoptions;
DIAInterpol[ell_, opts:OptionsPattern[]]:=DIAInterpol[ell, opts]=Block[{ret, OmegaM0, zz},
  OmegaM0 = OmegaM0Today[opts];
  Interpolation[
    Table[{zz, ((OmegaM0/Growth[zz, kofell[zz, ell, opts], opts])*curlyFIA[zz, opts])},
         {zz, $zminSurvey, $zmaxSurvey, ($zmaxSurvey-$zminSurvey)/$zintepoints}]
         ]
    ];


Options[pijIntegrandGIInterpol]=$paramoptions~Join~{dimensionless->True};
pijIntegrandGIInterpol[ell_?NumericQ,i_,j_,opts:OptionsPattern[]]:=
    pijIntegrandGIInterpol[ell,i,j,opts]=Block[{ret, OmegaM0,zz},
      OmegaM0=OmegaM0Today[opts];
      Interpolation[
      Table[{zz,(DIAInterpol[ell,opts][zz]*lensKernelGIInterpol[i,j,opts][zz]*pLimberInterpol[ell,opts][zz])},
      {zz, $zminSurvey,$zmaxSurvey,($zmaxSurvey-$zminSurvey)/$zintepoints}]
      ]
    ];

Options[pijIntegrandGIInterpolFunc]=$paramoptions~Join~{dimensionless->True};
pijIntegrandGIInterpolFunc[z_,ell_,i_,j_,opts:OptionsPattern[]]:=pijIntegrandGIInterpol[ell,i,j,opts][z]

Options[pijIntegrandII]=$paramoptions~Join~{dimensionless->True};
pijIntegrandII[z_?NumericQ,ell_?NumericQ,i_,j_,opts:OptionsPattern[]]:=
    pijIntegrandII[z,ell,i,j,opts]=Block[{ret, OmegaM0},
      OmegaM0=OmegaM0Today[opts];
      (OmegaM0/Growth[z,kofell[z,ell,opts],opts])^2*(curlyFIA[z, opts])^2*lensKernelII[z,i,j,opts]*pLimber[z,ell,opts]
    ];


Options[pijIntegrandIIInterpol]=$paramoptions~Join~{dimensionless->True};
pijIntegrandIIInterpol[ell_?NumericQ,i_,j_,opts:OptionsPattern[]]:=
    pijIntegrandII[ell,i,j,opts]=Block[{ret, OmegaM0,zz},
      OmegaM0=OmegaM0Today[opts];
      Interpolation[Table[{zz,
        ((DIAInterpol[ell,opts][zz])^2*lensKernelIIInterpol[i,j,opts][zz]*pLimberInterpol[ell,opts][zz])},
      {zz, $zminSurvey,$zmaxSurvey,($zmaxSurvey-$zminSurvey)/$zintepoints}]]
    ];

Options[pijIntegrandIIInterpolFunc]=$paramoptions~Join~{dimensionless->True};
pijIntegrandIIInterpolFunc[z_,ell_,i_,j_,opts:OptionsPattern[]]:=pijIntegrandIIInterpol[ell,i,j,opts][z]




Options[curlyFIA]=$paramoptions~Join~{modelIA->"e", luminosityFunction->$luminosityFileInterpolation};

curlyFIA[zz_, opts:OptionsPattern[]]:=Block[{switch=1, lumo, lumofunc, fia, aia, eia, bia, model, optos, cia=$cIA, sign=-1},
  model=OptionValue[modelIA];
  optos=complementParamValues[{opts},curlyFIA,returnList->"Full"];
  eia=$eIA/.optos;
  bia=$bIA/.optos;
  aia=$aIA/.optos;
  If[NumericQ[aia]==False, aia=$aIAfidu];
  If[NumericQ[bia]==False, bia=$bIAfidu];
  If[NumericQ[eia]==False, eia=$eIAfidu];
  lumofunc=OptionValue[luminosityFunction];
  lumo=lumofunc[zz];
  fia=Which[model=="c",
     1
     ,
     model=="d",
     (1+zz)^eia
     ,
     model=="e",
     (1+zz)^eia * lumo^bia
   ];
  debugPrint[model];
  debugPrint[lumo];
  debugPrint[fia];

  sign*aia*cia*fia
];


Options[pijIntegrandGalGal]=$paramoptions~Join~{dimensionless->True};
pijIntegrandGalGal[z_?NumericQ,ell_?NumericQ,i_,j_,opts:OptionsPattern[]]:=
    pijIntegrandGalGal[z,ell,i,j,opts]=Block[{ret},
      lensKernelGalGal[z,i,j,opts]*powerSpectrum[z,kofell[z,ell,opts],opts]
    ];



Options[dlensKernGGdp]={setEpsilon->$epsilonstep, dimensionless->True};

dlensKernGGdp[zz_?NumericQ,i_,j_,alpha_, opts:OptionsPattern[]]:=dlensKernGGdp[zz,i,j,alpha, opts]=Block[{opt1,opt2, step, epsi},
  epsi=OptionValue[setEpsilon];
  {opt1,opt2,step}=numericalDerivativeStep[alpha, epsi];
  (lensKernelGG[zz,i,j,opt1~Join~{opts}]-lensKernelGG[zz,i,j,opt2~Join~{opts}])/step
]



Options[dPLimberdp]={setEpsilon->$epsilonstep};

dPLimberdp[zz_?NumericQ,ell_?NumericQ,alpha_, opts:OptionsPattern[]]:=dPLimberdp[zz,ell,alpha, opts]=Block[{opt1,opt2, step, epsi},
  epsi=OptionValue[setEpsilon];
  {opt1,opt2,step}=numericalDerivativeStep[alpha, epsi];
  (powerSpectrum[zz,kofell[zz,ell,opt1],opt1]-powerSpectrum[zz,kofell[zz,ell,opt2],opt2])/step
]


unitsCij[opts:OptionsPattern[]]:=Block[{compunts,hval,units},
compunts=compatibleUnits[$internalPkUnits, invertUnits->False];
hval=hubbleToday[opts];
units=externalUnitsConverter[compunts,hval, physicalUnits->False, inverseDistanceUnits->True, internalUnits->True];
units=units^3;
Return[units]
];


Options[CijGG]=$paramoptions~Join~{shearshearIntegrand->pijIntegrandGG}

CijGG[ell_?NumericQ,i_,j_,opts:OptionsPattern[{CijGG,NIntegrate}]]:=CijGG[ell,i,j,opts]=If[i<j,CijGG[ell,j,i,opts],
  Block[{pIntegrandFuncGG=OptionValue[shearshearIntegrand],
    intopts, popts, copts, hval, units},
    popts=complementParamValues[{opts},CijGG,returnList->"Complement", filterCosmoPars->True];
    units = unitsCij[popts];
    units*NIntegrate[pIntegrandFuncGG[zaz,ell,i,j,popts],{zaz,$zminSurvey,$zmaxSurvey}]
  ]
];


Options[CijGI]=$paramoptions~Join~{shearIAIntegrand->pijIntegrandGI}

CijGI[ell_?NumericQ,i_,j_,opts:OptionsPattern[{CijGI,NIntegrate}]]:=CijGI[ell,i,j,opts]=If[i<j,CijGI[ell,j,i,opts],
  Block[{pIntegrandFuncGI=OptionValue[shearIAIntegrand],
    intopts, popts, copts, units},
    popts=complementParamValues[{opts},CijGI,returnList->"Complement", filterCosmoPars->True];
    units = unitsCij[popts];
    units*NIntegrate[pIntegrandFuncGI[zaz,ell,i,j,popts],{zaz,$zminSurvey,$zmaxSurvey}]
  ]
];


Options[CijII]=$paramoptions~Join~{IAIAIntegrand->pijIntegrandII}

CijII[ell_?NumericQ,i_,j_,opts:OptionsPattern[{CijII,NIntegrate}]]:=CijII[ell,i,j,opts]=If[i<j,CijII[ell,j,i,opts],
  Block[{pIntegrandFuncII=OptionValue[IAIAIntegrand],
    intopts, popts, copts, units},
    popts=complementParamValues[{opts},CijII,returnList->"Complement", filterCosmoPars->True];
    units = unitsCij[popts];
    units*NIntegrate[pIntegrandFuncII[zaz,ell,i,j,popts],{zaz,$zminSurvey,$zmaxSurvey}]
  ]
];


Options[CijGalGal]=$paramoptions~Join~{galaxygalaxyIntegrand->pijIntegrandGalGal}

CijGalGal[ell_?NumericQ,i_,j_,opts:OptionsPattern[{CijGalGal,NIntegrate}]]:=CijGalGal[ell,i,j,opts]=If[i<j,CijGalGal[ell,j,i,opts],
  Block[{pIntegrandFuncGalGal=OptionValue[galaxygalaxyIntegrand],
    intopts, popts, copts, units},
    popts=complementParamValues[{opts},CijGalGal,returnList->"Complement", filterCosmoPars->True];
    units = unitsCij[popts];
    units*NIntegrate[pIntegrandFuncGalGal[zaz,ell,i,j,popts],{zaz,$zminSurvey,$zmaxSurvey}]
  ]
];

Options[Cij]=$paramoptions
Cij[ell_?NumericQ,i_,j_,opts:OptionsPattern[{Cij,NIntegrate}]]:=Cij[ell,i,j,opts]=If[i<j,Cij[ell,j,i,opts],
  Block[{CijTotal,CijIAs,
     popts, copts},
    popts=complementParamValues[{opts},Cij,returnList->"Complement", filterCosmoPars->True];
    CijIAs=If[$IAswitch!=0,
      CijGI[ell,i,j,popts]+CijII[ell,i,j,popts],
      0];
    CijTotal=CijGG[ell,i,j,popts]+CijIAs;
    Return[CijTotal]
  ]
];



Options[dCijdp]=$paramoptions~Join~{cijFunction->Cij};

dCijdp[ell_,i_,j_,alpha_, opts:OptionsPattern[]]:=dCijdp[ell,i,j,alpha, opts]=Block[
  {opt1,opt2, CijFun, step},
  CijFun=OptionValue[cijFunction];
  {opt1,opt2,step}=numericalDerivativeStep[alpha, $epsilonstep];
  (CijFun[ell,i,j,opt1]-CijFun[ell,i,j,opt2])/step
];


dCijdpIntegral[ell_?NumericQ,i_,j_,alpha_]:=dCijdpIntegral[ell,i,j,alpha]=Block[{opt1,opt2, step},
  {opt1,opt2,step}=numericalDerivativeStep[alpha, $epsilonstep];
  NIntegrate[
    (dlensKernGGdp[zzz,i,j,alpha]*powerSpectrum[zzz,kofell[zzz,ell]]+dPLimberdp[zzz,ell,alpha]*lensKernelGG[zzz,i,j]),
    {zzz,$zminSurvey,$zmaxSurvey}]
];




noiseMatrix[i_,j_]:=noiseMatrix[i,j]=KroneckerDelta[i,j]*($intrinsicShear^2/(($galaxyDensArcminSq/$numZbinsWL)*3600*(180/Pi)^2))



Options[inversecovariance]=$paramoptions~Join~{cijFunction->Cij, numberOfZBins->$numZbinsWL};

SetSharedFunction[invcov];
ParallelEvaluate[invcov[Cijfunc_,ell_,i_,j_]:=invcov[Cijfunc, ell,i,j]=If[i<j, invcov[Cijfunc, ell,j,i],
  Cijfunc[ell,i,j]+noiseMatrix[i,j]
], First[Kernels[]]]

inversecovariance[ell_,opts:OptionsPattern[]]:=inversecovariance[ell,opts]=Block[{ret, CijFunc, Nbins=OptionValue[numberOfZBins]},
  CijFunc=OptionValue[cijFunction];
  ParallelTable[Table[invcov[CijFunc,ell, i, j], {j, 1, i}], {i, Nbins}];
  ret=Inverse[Table[invcov[CijFunc,ell, i, j],{i,Nbins},{j,Nbins}]];
  Return[ret]
]


Options[kScaleErrorDamping]={kDamping -> 1.5, dampingExponent->4, tolerance->10^(-15)};
kScaleErrorDamping[ell_?NumericQ, i_,opts : OptionsPattern[]] :=
  Block[{kda = OptionValue[kDamping], n=OptionValue[dampingExponent], to=OptionValue[tolerance]},
  Chop[Exp[-N[kofell[$zbinsEquiPopu[[i]], ell], $MachinePrecision]^n / kda^n], to]  ];

kScaleErrorDampingMatrix[ell_?NumericQ,
   opts : OptionsPattern[{kDamping -> 1.5}]] := DiagonalMatrix[Table[kScaleErrorDamping[ell, ii, opts], {ii, $numZbinsWL}]];


SetSharedFunction[jmat];
ParallelEvaluate[jmat[dCijfunc_,ell_,i_,j_,alph_]:=jmat[dCijfunc,ell,i,j,alph]=If[i<j, jmat[dCijfunc,ell,j,i,alph],
   dCijfunc[ell,i,j,alph]
], First[Kernels[]]]

Options[Jmat]=$paramoptions~Join~{dcijdpFunction->dCijdp, numberOfZBins->$numZbinsWL};

Jmat[ell_?NumericQ,alpha_,opts:OptionsPattern[]]:=Jmat[ell,alpha, opts]=Block[{ret, dCijdpFunc, Nbins=OptionValue[numberOfZBins]},
  dCijdpFunc=OptionValue[dcijdpFunction];
  ParallelTable[Table[jmat[dCijdpFunc,ell, i, j, alpha], {j, 1, i}], {i, Nbins}];
  ret=Table[jmat[dCijdpFunc,ell, i, j, alpha],{i,Nbins},{j,Nbins}];
  Return[ret]
];

Global`prog1=Global`prog2=Global`prog3=0;

Options[FisherWL]=$paramoptions~Join~{ellmin:>$ellmin,ellmax:>$ellmax,Deltaell:>$deltaell, cijFunction->Cij, dcijdpFunction->dCijdp};

FisherWL[alpha_,beta_,opts:OptionsPattern[]]:=FisherWL[alpha,beta,opts]=If[alpha<beta,FisherWL[beta,alpha,opts],
  Block[{CijFunc, dCijdpFunc, lmin=OptionValue[ellmin], lmax=OptionValue[ellmax], dell=OptionValue[Deltaell]},
    Global`prog1=Global`prog2=Global`prog3=0;
    SetSharedVariable[Global`prog1=Global`prog2=Global`prog3];
    $logellrange;
    $ellbinscenters;
    dCijdpFunc=OptionValue[dcijdpFunction];
    CijFunc=OptionValue[cijFunction];
    Sum[With[{ell=$ellbinscenters[[i]] , delell=(tentox@($logellrange[[i+1]])-tentox@($logellrange[[i]]))},
      Global`prog1++;
      ($fsky/2)*(2 ell+1)*delell*
          Tr[Jmat[ell,alpha, dcijdpFunction->dCijdpFunc].inversecovariance[ell, cijFunction->CijFunc].Jmat[ell,beta, dcijdpFunction->dCijdpFunc].inversecovariance[ell,cijFunction->CijFunc]]],
      {i,Length@$ellbinscenters-1}]]
];

Options[FisherWLSimple]=$paramoptions~Join~{cijFunction->Cij,dcijdpFunction->dCijdp, kDamping->False};

Global`progell=0;
Global`proga=0;
Global`progb=0;
Global`progt=0;

FisherWLSimple[alpha_,beta_,opts:OptionsPattern[]]:=Block[
{CijFunc, dCijdpFunc, kdamp, errmat, fishentry},
CijFunc=OptionValue[cijFunction];
dCijdpFunc=OptionValue[dcijdpFunction];
kdamp=OptionValue[kDamping];

If[alpha==1 && beta==1,
  If[UnsameQ[kdamp, False] && NumberQ[kdamp],
     Print["Applying covariance matrix damping at scale: "<>ToString[kdamp]]]];
(*Repeat $ellbinscenters calculation here, in case $logellrange or $ellmax, $ellmin are changed outside*)
$ellbinscenters=tentox@(MovingAverage[$logellrange,2]);
fishentry=Sum[Block[{
                    ell=$ellbinscenters[[i]],
                    delell=(tentox@($logellrange[[i+1]])-tentox@($logellrange[[i]])),
                    alphader, betader, invcov, cov,init,finit},
Global`progell=ell;
init=DateString[];
cov=Table[CijFunc[ell,b1,b2]+noiseMatrix[b1,b2], {b1,1,$numZbinsWL},{b2,1,$numZbinsWL}];
alphader=Table[dCijdpFunc[ell,b1,b2,alpha], {b1,1,$numZbinsWL},{b2,1,$numZbinsWL}];
betader=Table[dCijdpFunc[ell,b1,b2,beta], {b1,1,$numZbinsWL},{b2,1,$numZbinsWL}];
invcov=Inverse[cov];
finit=DateString[];
Global`progt=timeElapsed[init,finit,False];
If[UnsameQ[kdamp, False] && NumberQ[kdamp],
    errmat=kScaleErrorDampingMatrix[ell, kDamping->kdamp];
    invcov=errmat.invcov.Transpose[errmat];
];
($fsky/2)*(2 ell+1)*delell*Tr[alphader.invcov.betader.invcov] ],
{i,Length@$ellbinscenters-1} ];
Return[fishentry];
];


Options[FisherWLMatrix] = {keepParameters -> All};
FisherWLMatrix[opts : OptionsPattern[]]:=Block[
{fishMat, fishEval, i, j, rangepar, parlist, part, size},
part = OptionValue[keepParameters];
rangepar = Range@Length@$paramnames;
parlist = Part[rangepar, part];

  fishEval[a_, b_]:=fishEval[a, b]=fishEval[b, a]=FisherWLSimple[a, b]; (*memoization*)

  fishMat = Table[Global`proga=i; Global`progb=j; fishEval[i, j], {i, parlist}, {j, parlist}];
  Return[fishMat];
];



clearFisherWLMemoizedQuantities[opts:OptionsPattern[]]:=Module[{var},


  ParallelTable[removeDownValues[nGalaxyAntidev[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[nGalaxyAntidev[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[calculateZbins[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[calculateZbins[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[normalizationOfconvoluted$nGalaxyBins$PerfectPhotoZ[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[normalizationOfconvoluted$nGalaxyBins$PerfectPhotoZ[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[convoluted$nGalaxyBins$perfectPhotoZ[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[convoluted$nGalaxyBins$perfectPhotoZ[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[convoluted$nGalaxyBins$Gaussian[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[convoluted$nGalaxyBins$Gaussian[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[galaxyNorm$convoluted$nGalaxyBins$Gaussian[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[galaxyNorm$convoluted$nGalaxyBins$Gaussian[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[convoluted$nGalaxyBins$GaussianNormalized[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[convoluted$nGalaxyBins$GaussianNormalized[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[windowA[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[windowA[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[windowAInt[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[windowAInt[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[windowB[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[windowB[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[windowBInt[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[windowBInt[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[windowTildeGG[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[windowTildeGG[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[windowGammaGamma[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[windowGammaGamma[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[lensKernelGG[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[lensKernelGG[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[lensKernelGI[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[lensKernelGI[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[lensKernelII[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[lensKernelII[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[pijIntegrand2[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[pijIntegrand2[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[pijIntegrandGG[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[pijIntegrandGG[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[pijIntegrandGI[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[pijIntegrandGI[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[pijIntegrandII[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[pijIntegrandII[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[dlensKernGGdp[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[dlensKernGGdp[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[dPLimberdp[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[dPLimberdp[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[CijGG[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[CijGG[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[CijGI[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[CijGI[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[CijII[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[CijII[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[Cij[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[Cij[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[dCijdp[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[dCijdp[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[dCijdpIntegral[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[dCijdpIntegral[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[noiseMatrix[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[noiseMatrix[Except[_Pattern | _PatternTest],___]];

  removeDownValues[invcov[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[inversecovariance[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[inversecovariance[Except[_Pattern | _PatternTest],___]];

  removeDownValues[jmat[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[Jmat[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[Jmat[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[FisherWL[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[FisherWL[Except[_Pattern | _PatternTest],___]];

  ParallelTable[removeDownValues[FisherWLSimple[Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[FisherWLSimple[Except[_Pattern | _PatternTest],___]];

  ]

End[]
EndPackage[]
