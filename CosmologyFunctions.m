(* ::Package:: *)

BeginPackage["CosmologyFunctions`"]

Unprotect["CosmologyFunctions`*", "CosmologyFunctions`*`*"]  (*Unprotect and Protect definitions for using ClearAll in the case of developing code*)

ClearAll["CosmologyFunctions`*", "CosmologyFunctions`*`*"]


setParameterOptions::usage="set the default values of cosmological parameters by passing a list of fiducial values and a list of parameter names
setParameterOptions[paramNames_,paramFiducials_]"

setPScosmoOptions::usage="set the Boolean options for specifying the power spectrum and cosmology.
Default Options:
kdependentGrowth->False: Set to True if the Growth or GrowthRate depend on modes k.
lcdmBool->True: Set to False if a different cosmological model is going to be imported from external files.
linearBool->True : Set to False if the default power spectra used for calculations, should be non-linear power spectra.
"


setPowerSpectrumRelatedSpecs::usage="Set specifications that are related to the calculation of the observed power spectrum.
This are usually not changed very often. See Options[setPowerSpectrumRelatedSpecs] for a list of available specifications"

complementParamValues::usage=" complementParamValues[extrapars_,functionsym_,
opts:globalPars->$paramoptions,returnList->Complement,filterCosmoPars->False]
Finds the complement of the set of parameter options passed as extra arguments extrapars and
the default options for the function functionsym.
If returnList=Full, it returns the full list of parameters with the extrapars replacing the positions of the defaultpars.
If returnList=Default, it returns the default options for the function functionsym.
If filterCosmoPars->True it filters the final list, accepting only options defined in the local globalPars option
(usually cosmological parameters defined by $paramoptions)."

setCosmologyFixedValues::usage="set global cosmological variables, the gamma parameter of the growth rate and H0 in units of h/Mpc"
setTermsSwitches::usage="sets the options switches for different contributions, setTermsSwitches[betaswitch,dampSwitch],
 betaSwitch shuts down cosmological information in R.S.Distortions,
dampSwitch turns off the Seo&Eisenstein nonlinear damping term in the observed power spectrum"

setKandZMinMaxFixedValues::usage="setKandZMinMaxFixedValues[options],
Options:
Type ?kminIntegrations to see usage message.
Type ?kmaxIntegrations to see usage message.
Type ?kmaxInterpolations to see usage message.
Type ?zmaxInterpolations to see usage message.
Type ?zmaxIntegrations to see usage message.
Type ?fixGrowthScalek to see usage message.
"

wDEzFunc::usage="calculate w(z) from a parametrization of w0 and wa Form: wDEzFunc[z_number, w0_number, wa_number]=function"
curvatureSin::usage="Gives the Sinus function depending on the curvature value and sign. For Abs[ok]<=0.001 it returns just the argument x.
Form: curvatureSin[x,ok]= either (Sinh[x], Sin[x], x)"
OmegaM0Today::usage="Outputs the value of OmegaM today. Valid options: $paramoptions"
OmegaBaryon0Today::usage="Outputs the value of OmegaB (Baryons) today. Valid options: $paramoptions"
OmegaCDM0Today::usage="Outputs the value of OmegaC (Cold Dark Matter) today. Valid options: $paramoptions"
OmegaDE0Today::usage="Outputs the value of OmegaDE (Dark Energy or Lambda) today. Valid options: $paramoptions"
hubbleToday::usage="Outputs the dimensionless value of h (Hubble) for the given set of cosmological parameters."
OmegaK0Today::usage="Outputs the value of the curvature energy fraction OmegaK today. Computed from the deviation from flatness."


wEoSparameters::usage="Returns a list with {w0,wa} in the CPL parametrization, if they are defined as parameters in $paramoptions.
Otherwise default is: {-1,0}"
parametersModifiedGravity::usage=" Returns the value for a modified gravity parameter of name 'parname'.
If not found within $paramoptions, returns 0.
Form: parametersModifiedGravity[parname_,opts:OptionsPattern[]]"

inverseTransformation::usage="Option for parametersModifiedGravity.
Apply inverse transformation to the extra MG parameter specified in $paramoptions.
Used to vary transform(param) but passing param to Cosmomathica CAMB.
The function passed to this option should be the inverse of the transform function, written in pure form.
For example:
transform = Log10[#]& ;
inverseTransform = (10^#)& ;"


scalarAmplitude::usage="Returns the value of the scalar amplitude of fluctuations As. Options:
transformedParameter\[Rule]True : If set to True, the parameter As is obtained by the inverse function of the  transformed parameter.
See function setORreadTransformedParameters for details."
transformedParameter::usage="Option for scalarAmplitude. Set to True by default. If set to True, it is assumed that $paramoptions contains a transformed parameter of As"


setORreadTransformedParameters::usage="If a fundamental cosmological parameter is transformed with a symbolic function,
setReadTransformedParameters finds the inverse given the function and reads the actual bare value of the parameter.
It stores for each of the parameters in $paramoptions, the function that transforms the parameters and its inverse,
$funcsOfParameters[index] and  $inverseFuncsOfParameters[posi] respectively.
The Symbol of the transformed parameter 'par' must be a string specifying the symbolic function f[par]
or the Symbol  of the transformed parameter 'par', must contain the string 'par' in it and the function f[par]
must be passed separately.
Options:
generateInverseFunction -> False: If set to True, the inverse functions are calculated.
parameterToTransform -> None: Give here a symbol which is then used to match the parameter inside $paramoptions.
convertTransformedParameter -> False: If set to True, the parameter given by parameterToTransform is obtained using $inverseFuncsOfParameters[posi].
explicitTransformFunction->None: If function is not given from parameter string in $paramoptions, then pass here the function that transforms the parameter.

Warning: Certain symbolic functions have more than one inverse or need extra conditions,
please check the validity of the transform before proceding.
"


generateInverseFunction::usage="Option for setORreadTransformedParameters"
parameterToTransform::usage="Option for setORreadTransformedParameters"
convertTransformedParameter::usage="Option for setORreadTransformedParameters"
explicitTransformFunction::usage="Option for setORreadTransformedParameters"


externalUnitsConverter::usage="Converts the value of external and internal Hubble and Distance functions into units of h/Mpc (to be compatible with a power spectrum in Mpc^3/h^3).
Form: externalUnitsConverter[unitstr_String,hval_].
Options: $paramoptions, inverseDistanceUnits->False, physicalUnits->False.
For Distance units set the option inverseDistanceUnits to False.
If option physicalUnits is set to True the Hubble function is returned in [km/s/Mpc].
Accepted values for Hubble: 'h/Mpc', '1/Mpc', 'km/s/Mpc'
Accepted values for distanceAngular: 'Mpc', 'Mpc/h'
Internal code units for power spectrum k-modes is always [h/Mpc]"

compatibleUnits::usage=" Two possible Forms: compatibleUnits[unitstr_String]; compatibleUnits[sourceunit_String,targetunit_String, hval_Real].
If called in the first form, and option: invertUnits->True, returns the 'inverse' of the unit string. Example: 'Mpc/h' -> 'h/Mpc'; If option invertUnits->False, returns the given units.
If called in the second form, returns the numeric value needed to convert the sourceunit into the targetunit, usually dependent on hval.
If option inverseDistanceUnits->True, it performs this operation for units that are inverse distances."
testUnitString::usage="Form: testUnitString[ax_]. Test if the expression 'ax' corresponds to one of the valid internal units of the code and returns
the corresponding string."


Hubble::usage="Gives the value of the dimensionful Hubble parameter for other model different from LCDM if lcdmBool->False"
LCDMDimensionlessHubble::usage="Gives the value of the dimensionless Hubble parameter, at redshift z, accepts w0-wa parametrization for DE"

DimensionlessHubble::usage="Gives the value of the dimensionless Hubble parameter, at redshift z.
If lcdmBool->False, gives the DimensionlessHubble parameter from an external loaded function"

OmegaMLCDM::usage="Gives the value of OmegaM as a function of redshift z in LCDM"
OmegaM::usage="Gives the value of OmegaM in a general cosmology as a function of redshift, accepts the option lcdmBool->True for LCDM"
LCDMGrowth::usage="gives the LCDM growth factor as a function of z defined through OmegaM^gamma.
Gamma can be set globally using setCosmologyFixedValues"
Growth::usage="gives the Growth of a general cosmology as a function of z, accepts the option lcdmBool->True for LCDM"
fGrowthRate::usage="gives the growth rate f of matter perturbations for LCDM or a general cosmology. Accepts the option lcdmBool->True for LCDM"



globalPars::usage="Option for complementParamValues"
returnList::usage="Option for complementParamValues"
filterCosmoPars::usage="Option for complementParamValues"
debugOptions::usage="Option for complementParamValues"


setBiasFunction::usage="sets the bias interpolation function when given a list of redshifts and a list of bias values
Accepts options for Interpolation"
betaRSDfunction::usage="gives the redshift space distortion beta function as a function of z. rsdbeta=fGrowthRate/bias"
distanceAngular::usage="gives the angular diameter distance as a function of z"
comovingDistance::usage="gives the comoving distance as a function of z"
distanceSoundHorizon::usage="gives the sound horizon at CMB time as a function of the decoupling redshift"
dimensionless::usage="Default: dimensionless->True. Option for comovingDistance and other functions. If set to True, the dimensionless quantity is returned."

dampingTerm::usage="dampingTerm[z,k,mu,opts] gives the nonlinear damping term correction to the observed power spectrum.
Needs the setting of the variables $dampingGrowthNorm, $dampingSigma0Value and $dampSwitch"


volumeZbin::usage="Gives the volume contained in a z bin given by z1<z<z2.
Form: volumeZbin[z1,z2,opts]. Uses the angular diameter distance which has to be in the correct units of h/Mpc"


sigmaZ::usage="sets the redshift error of spectroscopic measurements"
sigmaR::usage="sets the spectroscopic redshift error as a function of redshift"
sigmaV::usage="sets the error on redshift due to nonlinear velocity dispersions. Needs $sigmaPecVel to be set in km/s."
errorZ::usage="redshift error as a function of scale k, angle mu and redshift z with optional cosmological parameters given by opts
Form: errorZ[k,mu,z,opts]"

windowk::usage="Window function to integrate over in calculation of sigma8 for the power spectrum"
sigma8Function::usage="Calculates sigma8 for a power spectrum interpolation function at a radius scale scal
Form: sigma8Function[ps, scal, opts] Takes opts for NIntegrate"
normalizationFactor::usage="normalizationFactor[ps_,sigma8_] calculates sigma8^2/sigma8Function[ps,$radiusScale]
$radiusScale set by default to 8 "

sigma8ofZ::usage="sigma8ofZ[z, opts]:Calculates the function sigma8(z) by integrating the powerspectrum over the windowk function at a radius of $radiusScale.
Accepts options for parameters: $paramoptions and $pscosmoopts and options for NIntegrate."

topHatWindowR::usage="Radius scale for the top hat window function used to calculate sigma8ofZ. Default value: $radiusScale in Mpc/h."

rescaleSigma8fromAs::usage="rescaleSigma8fromAs[s8ref,asvalue]: Given a reference value for sigma8 and a modified value for the scalar amplitude As,
outputs the correspondingly rescaled sigma8."
rescaleAsfromSigma8::usage="rescaleAsfromSigma8[s8value, s8ref, asref]: Given a reference value for sigma8 and As, and a modified value for the amplitude sigma8,
outputs the correspondingly rescaled As."

LCDMCAMBPsPre::usage="Gives LogLog power spectrum from CAMB LCDM, Needs cosmomathica package. Options: $paramoptions and linearBool->True"
LcdmCambPk::usage=" LcdmCambPk[zred, k, sigma8, opts], Gives power spectrum interpolation function from CAMB LCDM, Needs cosmomathica package. Options: $paramoptions and linearBool->True"
returnInterpolated::usage="Option for LCDMCAMBPsPre. If Set to True, returns the table of k, Pklin and Pknonlin values.
Otherwise it returns an interpolating function in (k,P(k))."
observablesCAMB::usage="Function to extract observables from CAMB cosmomathica data array.";
returnObservable::usage="Option for observablesCAMB. Observables can be: s8fgDz, Hubble, Growth, GrowthRate."
checkConsistency::usage="Option for LCDMCAMBPsPre. If set to True, the code is run in such a way that
only one parameter changes at a time. Useful for sigma8."
returnRescaledAs::usage="Option for LCDMCAMBPsPre. If Set to True, returns the rescaled amplitude As, in order to match the wanted sigma8.
Does not return power spectra."
returnSigma8::usage="Option for LCDMCAMBPsPre. If Set to True, returns the sigma8 value, together with
the spectra"
returnGrowthRate::usage="Option for LCDMCAMBPsPre. If Set to True, returns  a table with f(z), D(z) and sigma8(z) values, provided LCDMCAMBPsPre is evaluated with a list of redshifts."
returnDeprecated::usage="Option for LCDMCAMBPsPre.
If Set to True, returns the power spectrum in a deprecated format, used for the
IST:F WL format=1 in GenerateInput.m."
activateModifiedGravity::usage="Option for LCDMCAMBPsPre.
If set to False, MGCAMB options are set to default, by passing an empty set {}.
If set to True, MGCAMB options have to be passed manually, using optionsMG.
Predefined models that can be set by passing a string:
'fofR-HS'
'mueta-param'
"
optionsMG::usage="Option for LCDMCAMBPsPre. Used for passing directly MGCAMB flags to Cosmomathica CAMB.
Active only if 'activateModifiedGravity' is set to True."

LCDMTrHuPre::usage="LCDMTrHuPre[sigma8val_,opts:], Gives power spectrum from Eisenstein&Hu Transfer functions,
using the Growth function of this package. opts: cosmological parameters and linearBool->True. If linearBool->False,
nonlinear PS is calculated with the function HalofitCorrection"

powerSpectrum::usage="powerSpectrum[z_,k_,opts]. Gives the power spectrum for a specified cosmology and a calculation method.
Options=$paramoptions, $pscosmoopts and the extra options:
spectrumMethod\[Rule](Transfer or CAMB) calculation with Eisenstein&Hu Transfer functions or CAMB from cosmomathica
sigma8reference\[Rule]0.8 (reference sigma8 normalization for PS, just in case it is not specified as a cosmological parameter)"

dlnPdkDerivativeFunction::usage=" dlnPdkDerivativeFunction[zred_,kr_,opts], Gives the numerical Logarithm derivative of P(k) wrt k, as a function of k and z.
Accepts all options that are accepted by powerSpectrum"


pthetathetaInt::usage="Integral of the theta-theta power spectrum divided by k^2 (linear theory velocity dispersion), which simplifies to an integral over P_matter(k). Accepts the same options as powerSpectrum."
sigmapNL::usage="Function computing the sigmapnl parameter, entering the Lorentzian factor of the FoG.  Accepts the same options as powerSpectrum.
To match the recipe, this actually computes sigmap(mean)/sigma8(mean). If the option computeNLevolution==True, it computes the z-dependent sigmapnl/sigma8. Only available for the $fbsigma8 recipe nat the moment."
sigmavNL::usage="Function computing the sigmavnl parameter, entering the BAO damping term g_mu. Accepts the same options as powerSpectrum."


ignoreSigmaPVCosmoDependence::usage="Option for FingersOfGod and BAOdamping,
if set to True, sigmapNLNew and sigmavNLNew ignore the explicit dependence on cosmo parameters.
The only dependence enters then through kAP and muAP."
FingersOfGod::usage="FingersOfGod[zz_, kk_, muu_, opts:OptionsPattern[powerSpectrum]],
Fingers of God term for GCsp. Depending on the value of $kdependentGrowth, there are two
different approaches."
BAOdamping::usage="BAOdamping[zz_, kk_, muu_, opts:OptionsPattern[powerSpectrum]],
BAO damping term for GCsp. Depending on the value of $kdependentGrowth,
there are two different approaches."


epsilonChangeParameter::usage="Form: epsilonChangeParameter[par_, epsi_:$epsilonstep, paropts_:$paramoptions].
Gives value of parameter par, changed by an epsilon epsi, according to the fiducials in paropts."

SigmaLensing::usage="Function that computes the Sigma lensing function, altering the Weyl potential, which is the one appearing in the Poisson equation for relativistic particles (sum of Psi and Phi).
For LCDM and GR, this is identically equal to 1."

externalSigmaLensingFunction::usage="External function for the SigmaLensing function."
$externalSigmaLensingInterpolatingFunction::usage="Interpolating function for the external SigmaLensing function."
$externalSigmaLensingDerivativesInterpolatingFunction::usage="Interpolating function for the
derivatives of the external SigmaLensing function."
theoreticalMGFunction::usage="Option for SigmaLensing. Specifies a theoretical Sigma function defined for modified gravity. Default: False"
theoreticalKdependence::usage="Option for SigmaLensing. Specifies if theoreticalMGFunction is k-dependent. Default: False"

$vd4::usage="Function to vary numerically the z-dependent function corresponding to d4. Default: $vd4[z]:=1"
$vd5::usage="Function to vary numerically the z-dependent function corresponding to d5. Default: $vd5[z]:=1"
$vd7::usage="Function to vary numerically the z-dependent function corresponding to d7. Default: $vd7[z]:=1"
$vd8::usage="Function to vary numerically the z-dependent function corresponding to d8. Default: $vd8[z]:=1"


RAPfunction::usage=" RAPfunction[zred_,mur_,opts:] : Gives the R value of the Alcock Pazcynski effect
as a function of redshift and mu. Options: $paramoptions"


numericalDerivativeStep::usage="Form: numericalDerivativeStep[index_, epsi_].
Computes the steps and the parameter options for a simple numerical derivative with 3point stencil. Accepts the index of the parameter corresponding to its
position in $paramoptions and the epsilon value epsi."

variationOfzdependentQuantities::usage="Function that activates the variation of z-dependent quantities numerically. It accepts any member of $zdependDerivVector."
clearZdependentVariations::usage="Function that clears all the variation of z-dependent quantities, for example $vd4."

varylogquantities::usage="Option for variationOfzdependentQuantities, if set to True, the epsilon is applied to the natural log of the corresponding quantity."

kmuAlcockPaczynski::usage="kmuAlcockPaczynski[zr_,kr_,mur_,opts:] :
Gives a list {k,mu} with the modified k and mu values due to the Alcock-Paczynski effect as a function of redshift "

observedPowerSpectrum::usage=" observedPowerSpectrum[zred_,k_,mu_,opts:] :
Gives the observed power spectrum as a function of redshift z, scale k and observing angle mu.
It includes the raw P(k), bias, redshift space distortions, photometric redshift error, shot noise,
and the Alcock-Paczynski effect on k,mu and the volume.
Options: $paramoptions, APeffect->True, and all other options available for the powerSpectrum function"

$biasInterpFunc::usage="interpolation function for the fiducial bias"

$PowerExtraShot::usage="Set some function of z to account for some extra shot noise term in the observed power spectrum"


$shaParVarInZdep::usage="Global setting parameter option for observedPowerSpectrum and related quantities.
Meaning: shape parameters variation in z-dependent quantities.
If set to False, z-dependent quantities do not depend on the shape parameters.
If set to True, all z-dependent quantities also depend on shape parameters."

zdepFunctionsCosmoVariation::usage="Option for observedPowerSpectrum, that sets the variable $shaParVarInZdep."

$paramoptions::usage="list of rules between protected parameter names and parameter fiducial values"
$paramfidus::usage="list of fiducial values for parameters"
$paramnames::usage="list of protected parameter names for fiducials"
$paramlabels::usage="list of parameter label names for fiducials (i.e. use in plots)"
$parampositions::usage="list of rules between parameter names and parameter position in vector of parameters"
$ncosmopars::usage="number of fixed cosmological fiducial parameters"
$pslinearopt::usage="Boolean option  specifying linear (True) or nonlinear (False) power spectrum calculation"
$pscosmoopts::usage="Boolean options for cosmology and power spectrum, usually: lcdmBool->True and linearBool->True"
$modcosmoopt::usage="Boolean option specifying LCDM (True) or non-LCDM (False) cosmological functions"
$partesta::usage="debug variable"


$modelName::usage="Representative name of model chosen, either linear or nonlinear and LCDM or external custom cosmology"


$runtype::usage="Type of run either linear or nonlinear and LCDM or external custom cosmology"


$kminRef::usage="Minimum k value where power spectrum is defined or to be used in k integrals"
$kmaxRef::usage="Maximum reference k value for integrals or other calculations in k. Used in sigma8Function"
$kMaxInt::usage="Maximum allowed k for interpolations and computations in k. Used in dlnP/dlnk"
$zMaxIntp::usage="sets the maximum z until which to interpolate internal power spectrum functions."
$zMaxIntegral::usage="sets the maximum z for upper value of integration of Hubble distances and related quantities."

kminIntegrations::usage="Option for setKandZMinMaxFixedValues to set the mimimum k in integrations of some functions. Default: $kminRef"
kmaxIntegrations::usage="Option for setKandZMinMaxFixedValues to set the maximum k in integrations of some functions. Default: $kmaxRef"
kmaxInterpolations::usage="Option for setKandZMinMaxFixedValues to set the maximum k in interpolations of some functions. Default: $kMaxInt"
zmaxInterpolations::usage="Option for setKandZMinMaxFixedValues to set the maximum z in interpolations of power spectra functions. Default: $zMaxIntp"
zmaxIntegrations::usage="Option for setKandZMinMaxFixedValues to set the maximum z in integrations of distances and related quantities. Default: $zMaxIntegral"
fixGrowthScalek::usage="Option for setKandZMinMaxFixedValues and the density growth functions (Growth, fGrowthRate), that fixes a scale k for the k-dependent growth of perturbations."


$fkfix::usage="Value of a fixed scale k, for the growth of density perturbations in the case they are k-dependent. Default value: False, which means that full k-dependence is taken into account. "


spectrumMethod::usage="Option for power spectrum, 'Transfer', 'CAMB', 'CosmicEmulator' or 'ExternalFile'. Default: Transfer"

sigma8reference::usage="Option for power spectrum, when sigma8 is not specified as a parameter. Default: 0.8"
Asreference::usage="Option for power spectrum, when As is not specified as a parameter. Default: 2.1265*10^-9"

APeffect::usage="Option for observed power spectrum. Activates or deactivates the Alcock-Pazcynski effect, Default: False"

varyZdependentparameter::usage="Option for observed power spectrum. Default: False. If a member of $zdependDerivVector is passed, then
 it varies numerically with an $epsilonzstep the corresponding function."

lcdmBool::usage="Option for the package. If set to True, cosmological functions are calculated for LCDM using internal functions.
If set to false, cosmological functions are taken from loaded input files. Which functions are taken is specified
in the externalCosmoFunctions options"

linearBool::usage="Option for the package. If set to True, only linear quantities are used for the power spectrum calculation.
If set to false, nonlinear quantities are used for the power spectrum calculation."

$lightspeed::usage="Exact speed of light in units of km/s"

$H0::usage="value of H0 reference in units of h/Mpc"

$paramsDirectoryNames::usage="List of the names of directories and files for external input of cosmological functions,
corresponding to the each cosmological parameter varied in relative +/- epsilon.
pl = par*(1+eps)
mn = par*(1-eps)
Exception: when par=0, then use absolute epsilons.
Ordering of parameters should be the same as for the list $paramoptions with the minus term preceding the plus term. The fiducial directory name should be the first in the list."

internalHubbleUnits::usage="Option for setCosmologyFixedValues that sets the units to be used internally in the LCDM code for the Hubble function. Desirably 'h/Mpc'"
internalDistanceUnits::usage="Option for setCosmologyFixedValues that sets the units to be used internally in the LCDM code for the Distance function. Desirably 'Mpc/h'"


$externalH0units::usage="This parameter has to be set to the corresponding units of the external input files: 'h/Mpc' or '1/Mpc'."
$externalDistanceUnits::usage="This parameter has to be set to the corresponding units of the external input files: 'h/Mpc' or '1/Mpc'."
$internalH0units::usage="This parameter has to be set to the wanted units of the internal Hubble function: 'h/Mpc', '1/Mpc' or 'km/s/Mpc.
Internally the code matches the units to those of the power spectrum in order to calculate survey volumes and densities correctly.
However, for derivatives it can have other different units."
$internalDistanceUnits::usage="This parameter has to be set to the wanted units of the internal Distance function: 'Mpc/h' or 'Mpc'.
Internally the code matches the units to those of the power spectrum in order to calculate survey volumes and densities correctly.
However, for derivatives it can have other different units."

externalHubbleUnits::usage="Option for setCosmologyFixedValues. String specifying the units of the external input file containing H(z).
Possible units: 'h/Mpc', '1/Mpc', 'km/s/Mpc'.
Internally in the code, H(z) is transformed to units of h/Mpc."

externalDistanceUnits::usage="Option for setCosmologyFixedValues. String specifying the units of the external input file containing D(z) (distance).
Possible units: 'Mpc/h', 'Mpc'.
Internally in the code, distanceAngular(z) is transformed to units of Mpc/h, for compatibility with Fisher volume number density."

internalUnits::usage="Option for externalUnitsConverter in H(z). If set to True, calculates the needed units to multiply in front of the dimensionless Hubble function."
inverseDistanceUnits::usage="Option for externalUnitsConverter and compatibleUnits. If set to True, returns units inverse to the distance units."
invertUnits::usage="Option for compatibleUnits. If set to True, compatibleUnits returns the 'inverse' string of the supplied string. Otherwise it returns the input string."

externalPkUnits::usage="Option for setCosmologyFixedValues, sets the external units of k for the external power spectrum."
internalPkUnits::usage="Option for setCosmologyFixedValues, sets the internal units of k for the internal power spectrum."

$externalPkUnits::usage="Variable that sets the external units of k for the external power spectrum."
$internalPkUnits::usage="Variable that sets the internal units of k for the internal power spectrum."

pkUnitsConverter::usage="Function that sets the correct units for k and P(k). A string with the desired units is passed as a variable and the value of h."
units::usage="Option for powerSpectrum, to set the units of k and correspondingly of P(k) for the powerSpectrum.
Default: $internalPkUnits='h/Mpc'. It can also take the value '1/Mpc'."

fGammaFit::usage="Option for setCosmologyFixedValues. Gamma value appearing in the expression Omega^gamma that fits the LCDM growth rate f.
Default value fGammaFit=$fGamma=5.0/9.0"

$fGamma::usage="Gamma value appearing in the expression Omega^gamma that fits the LCDM growth rate f. Default value $fGamma=5.0/9.0"

hubbleReference::usage="Option for setCosmologyFixedValues. In case h is a derived parameter, this sets the default value for h,
when all other parameters are at its fiducial value and cannot be obtained from them.
In this case it should be set at the beginning of all calculations"

physicalUnits::usage="Option for Hubble function, default set to False.
If set to true, it will use the previously set units to convert the result to km/s/Mpc"


OmegaExtraMatterSpecies::usage="Option for setCosmologyFixedValues that sets a value of a matter component in the Hubble function,
which is itself not a parameter of the fiducial model. E.g. Neutrinos."
OmegaRadiation::usage="Option for setCosmologyFixedValues that sets a fixed value of a radiation component in the Hubble function."


$Omegabfixed::usage="Default: $Omegabfixed=0. Reference value for a fixed Omegab (baryons), when Omegab is not considered a varying parameter in $paramoptions."

OmegabReference::usage="Default: OmegabReference->$Omegabfixed. Option for setCosmologyFixedValues. This parameter has to be set to some value if one does not want to use Omegab as a varying cosmological parameter in $paramoptions."

getCosmoParameterFileName::usage="Get the correct file from the external cosmological input files, corresponding to the given numerical cosmological parameter"

inputNumericalDerivatives::usage="Option for getCosmoParameterFileName and externalPowerSpectrumFunction for activating the input of numerical derivatives from external files."
getNonFiducialFileName::usage="Option for getCosmoParameterFileName. If set to True, getCosmoParameterFileName returns a random file name which is not the fiducial one."
externalHubbleFunction::usage="Gives the value of the Hubble function computed from external input files. Form: externalHubbleFunction[z, opts]"

externalDistanceFunction::usage="Gives the value of the Distance function computed from external input files. Form: externalDistanceFunction[z, opts]"


shiftedDomain::usage="Option for external cosmological background functions. If the external function is not defined at z=0,
then the first entry of the file is taken as the initial value zdom and the function is evaluated at z=z+zdom."


externalGrowthRateFunction::usage="Gives the growth rate f, from an external input file, for any value of z and (optionally)
k and the cosmological parameters as options. For a k-dependent Growth, use the option kdependentGrowth->True.
Form: externalGrowthRateFunction[z, k(optional),opts]"

externalGrowthFunction::usage="Gives the growth function G, from an external input file, for any value of z and (optionally)
k and the cosmological parameters as options. For a k-dependent Growth, use the option kdependentGrowth->True. Form:
externalGrowthFunction[z, k(optional),opts]"

normalizeGrowth::usage="Option for external Growth function. If True, Growth is normalized using input files:
externalGrowthFunction[z]=$externalGrowthInterpolatingFunction[z]/$externalGrowthInterpolatingFunction[0]."

extGrowthDz::usage="Gives the growth function G, from the integration of the growth rate f(z,k) from an external file. (k is optional)
Form: extGrowthDz[zk__,opts]"

solveDifferentialEquation::usage="Option for fGrowthRate. If set to True, function f(z) is calculated from solving differential equation.
Includes support for the CPL parametrization w0-wa."

zOfNf::usage="Function that converts N-foldings Nf into redshift z."
NfOfz::usage="Function that converts redshift z into N-foldings Nf."
aOfz::usage="Function that converts redshift z into scale factor a."
zOfa::usage="Function that converts scale factor a into redshift z."


externalSigma8ofZFunction::usage="Function that gives the cosmological sigma8(z) function, from an external input file."

externalPowerSpectrumFunction::usage="Gives the power spectrum function from an external file at the respective z and k value.
The list of ranges in z for the external files has to be set correctly using the option filesRangeInZ or setting the variable  $zrangePowerSpectrum.
If the z value given as argument to externalPowerSpectrumFunction does not correspond to a value in $zrangePowerSpectrum, then 'interpolation' using the Growth function is used."


filesRangeInZ::usage="Option for setting the list of values in z at which the external power spectrum files are defined."
interpolatedInLogLog::usage= "Option for externalPowerSpectrumFunction. If set to True,
then it is because the power spectrum input files were interpolated in Log10,Log10 space.
Therefore a convertion is applied to get P(k)."

functionInterpolatedInLogLog::usage="Option for externalFunctionTaylor. If set to True,
the fiducial function (f(x_0)) has been interpolated in Log10,Log10 space (x,y)."

derivativeInterpolatedInLogLin::usage="Option for externalFunctionTaylor. If set to True,
the derivative of the function (f'(x_0)) has been interpolated in Log10,Linear space (x,y)."


externalFunctionTaylor::extfuncs="In function externalFunctionTaylor, the externalFunction or its derivative are not set to a function."
externalFunctionTaylor::usage="Form: externalFunctionTaylor[zkarg__,opts:OptionsPattern[]]. Function to calculate the value of a function around
the fiducial point in parameter space, given a function and the derivative of that function w.r.t. the parameters.
The function can depend on z or on both z and k. Check options and the usage messages for more information."

powerFit::usage="Tool to compute a linear fit on the logarithm of a function. Used for extrapolation purposes."
functionExtrapolationFit::usage="Function used to extrapolate k-dependent functions. Using a power law fit (linear fit in log-log space)."

externalFunction::usage="Option for externalFunctionTaylor. This option provides the Taylor function with the externalFunction for the specified quantity."
externalDerivativeFunction::usage="Option for externalFunctionTaylor. This option provides the Taylor function with the externalDerivativeFunction
for the derivative of the specified quantity."
interpolatedArguments::usage="Option for externalFunctionTaylor. Possible values:
'z':  The function is interpolated in z and depends only on z.
'k':  The function is interpolated in k, while there might be function for each redshift z.
'zk': The function is interpolated in both z and k."
kdependence::usage="Option for externalFunctionTaylor. If set to True, the passed externalFunction is both z and k dependent."

externalSpectraType::usage="Option for external power spectrum function, specifying which kind of spectra are read from files."

getDomain::usage="Option for different external interpolated functions to show the valid
domain instead of producing output."

externalPowerSpectrumFunction::zrange="z-range of external Power Spectrum files not set properly. Returning: Null."


$zrangePowerSpectrum::usage="List of values in z at which the external power spectrum files are defined."
zValTozIndex::usage="Function to convert from z-value to index in z-range of specified files."
zIndexTozVal::usage="Function to convert from z-index to z-value in z-range of specified files."

setExternalCosmoInterpolatingFunctions::usage="Set the corresponding interpolating functions or lists of interpolating functions for
the external cosmological input files needed.
Each cosmological interpolation function is indexed according to the parameter names set by $paramsDirectoryNames.
Each cosmological function is set using the corresponding option. If the option is not set, the function is set to 'False'
If an interpolation file is given for a cosmological function (Hubble, distanceAngular, Growth, fGrowthRate or powerSpectrum),
the option 'externalFile' is set to True in the corresponding global cosmological function"

externalHubbleInput::usage="Option for setExternalCosmoInterpolatingFunctions"
externalDistanceInput::usage="Option for setExternalCosmoInterpolatingFunctions"
externalGrowthInput::usage="Option for setExternalCosmoInterpolatingFunctions"
externalGrowthRateInput::usage="Option for setExternalCosmoInterpolatingFunctions"
externalPowerSpectrumInput::usage="Option for setExternalCosmoInterpolatingFunctions"
externalPowerSpectrumDerivativesInput::usage="Option for setExternalCosmoInterpolatingFunctions"
externalSigma8ofZInput::usage="Option for setExternalCosmoInterpolatingFunctions"
externalPowerSpectrumNoWiggleInput::usage="Option for setExternalCosmoInterpolatingFunctions"
externalPowerSpectrumNoWiggleDerivativesInput::usage="Option for setExternalCosmoInterpolatingFunctions"

parameterDirectoriesNames::usage="Option for the function getCosmoParameterFileName, sets the correct order of the directories which contain the external input files ffor the corresponding cosmological parameters.
Default: $paramsDirectoryNames"


kdependentGrowth::usage="Option for setPSCosmoOptions. If True, Growth and GrowthRate are k-dependent either from internal functions or from external input."
$kdependentGrowth::usage="Value of the kdependentGrowth Option. If set to True, Growth and Growth Rate are calculated as k and z-dependent functions."


growthDerivative::usage="Option for the fgrowthRate Function. If set to True, fGrowthRate is calculated from the derivative of the function Growth."
growthRateIntegral::usage="Option for the Growth function. If set to True, Growth is obtained from integrating the fGrowthRate in time."


$externalHubbleInterpolatingFunction::usage="external input interpolating function for the Hubble function"
externalHubbleDerivativesInput::usage="Option for setExternalCosmoInterpolatingFunctions. If set to True, the externalHubble function accepts derivatives."
$externalHubbleDerivativesInterpolatingFunction::usage="Function that stores the interpolating function for the external derivative of the Hubble function."

$externalDistanceInterpolatingFunction::usage="external input interpolating function for the Angular Diameter Distance function"
externalDistanceDerivativesInput::usage="Option for setExternalCosmoInterpolatingFunctions. If set to True, the externalDistance function accepts derivatives."
$externalDistanceDerivativesInterpolatingFunction::usage="Function that stores the interpolating function for the external derivative of the angular diameter distance function."

$externalGrowthInterpolatingFunction::usage="external input interpolating function for the linear Growth function"
$externalGrowthDerivativesInterpolatingFunction::usage="external input interpolating function for the derivatives of the linear Growth function"

$externalGrowthRateInterpolatingFunction::usage="external input interpolating function for the linear Growth Rate function"
$externalGrowthRateDerivativesInterpolatingFunction::usage="external input interpolating function for the derivatives of the Growth Rate f function"

$externalPowerSpectrumInterpolatingFunction::usage="external input interpolating function for the Power Spectrum function"
$externalPowerSpectrumDerivativesInterpolatingFunction::usage="external input interpolating function for the derivatives of the Power Spectrum
with respect to the cosmological parameters."
$externalPowerSpectrumNoWiggleInterpolatingFunction::usage="external input interpolating function for the non-wiggle power Spectrum function"
$externalPowerSpectrumNoWiggleDerivativesInterpolatingFunction::usage="external input interpolating function for the non-wiggle power Spectrum function"

$externalSigma8ofZInterpolatingFunction::usage="external input interpolating function for the linear sigma8(z) function"
$externalSigma8ofZDerivativesInterpolatingFunction::usage="external input interpolating function for the derivatives of the linear sigma8(z) function"

$externalfGrowthRateSigma8ofZInterpolatingFunction::usage="external input interpolating function for the linear fsigma8(z) function"
$externalfGrowthRateSigma8ofZDerivativesInterpolatingFunction::usage="external input interpolating function for the derivatives of the linear fsigma8(z) function"

externalGrowthDerivativesInput::usage="Option for setExternalCosmoInterpolatingFunctions. If set to True, the externalGrowth function accepts derivatives."
externalGrowthRateDerivativesInput::usage="Option for setExternalCosmoInterpolatingFunctions. If set to True, the externalfGrowthRate function accepts derivatives."
externalSigma8ofZDerivativesInput::usage="Option for setExternalCosmoInterpolatingFunctions. If set to True, the externalSigma8ofZ function accepts derivatives."

externalfGrowthRateSigma8ofZInput::usage="Option for setExternalCosmoInterpolatingFunctions. If set to True, the externalfGrowthRateSigma8ofZ function accepts external function."
externalfGrowthRateSigma8ofZDerivativesInput::usage="Option for setExternalCosmoInterpolatingFunctions. If set to True, the externalfGrowthRateSigma8ofZ function accepts derivatives."

externalSigmaLensingInput::usage="Option for setExternalCosmoInterpolatingFunctions. If set to True, the externalSigmaLensing function accepts external function."
externalSigmaLensingDerivativesInput::usage="Option for setExternalCosmoInterpolatingFunctions. If set to True, the externalSigmaLensing function accepts external derivatives."

externalfGrowthRateSigma8ofZFunction::usage="Function that accepts external interpolating functions and passes it to fGrowthRateSigma8ofZ"
fGrowthRateSigma8ofZ::usage="Function fGrowthRateSigma8ofZ, fsigma8(z)."

externalFile::usage="Option for importing external cosmology input files"


$dampingGrowthNorm::usage="The normalization of the Growth function that goes into the Eisenstein Lagrangian damping"
$dampingSigma0Value::usage="The Sigma0 value (phenomenological) that goes into the Eisenstein Lagrangian damping"
$radiusScale::usage="The default radius of 8 Mpc/h for the calculation of the normalization sigma8. Under a change of units, it is automatically changes in the sigma8 function."
$dampingSwitch::usage="If set to zero, turns off the Eisenstein Lagrangian damping on the observed power spectrum"
$dzSpectError::usage="Error in redshift dz in spectroscopic measurements, due to relative velocities"
$sigmaPecVel::usage="Damping induced in redshift space by nonlinear peculiar velocity dispersion in km/s"
$APeffectBool::usage="Activates or deactivates the Alcock-Pazcynski effect. Default: False"
$APswitch::usage="The value of this variable is coupled to the value of $APeffectBool.
If set to 0, the AP effect is turned off completely.
If set to 1, it is taken into account in the derivatives dlnP/dnlH, dlnP/dlnDa. Default: 0"
$APoptString::usage="The value of this string is coupled to the value of $APeffectBool.
Shows either: \"yes-AP\" or \"no-AP\" depending if the AP effect is turned on or off."
$pks8ratio::usage="Switch, if set to 1 then the external powerSpectrum from input files is the ratio P(z,k)/s8(z)^2. If set to zero, it is P(z,k). The powerSpectrum[] function returns always the true P(z,k)"

$fbsigmaBool::usage="If set to True, the variables fs8 and bs8 are used in the calculation of the observed power spectrum derivatives."

radius8MpcScale::usage="Option for the default radius of 8 Mpc/h for the calculation of the normalization sigma8.
Default value: $radiusScale=8"
dampingGrowthNorm::usage="Option for the normalization of the Growth function that goes into the Seo&Eisenstein Lagrangian damping.
Default value: $dampingGrowthNorm=0.758."
dampingSigma0::usage="Option for the Sigma0 value (phenomenological) that goes into the Seo&Eisenstein Lagrangian damping.
Default value: $dampingSigma0Value=11."
dampingSwitch::usage="Option for turning off the Eisenstein Lagrangian damping on the observed power spectrum. If set to zero, damping is turned off.
Default value: $dampingSwitch=0"
dzSpectroscopicError::usage="Option for setting error in redshift dz induced by spectroscopic errors, set number as a relative error (not percent). Default value: $dzSpectError=0."
dzRelativeVelocityError::usage="Option for setting error in redshift dz induced by spectroscopic errors (previously called: relative velocities, left here temporarily for compatibility issues),
set number as a relative error (not percent). Default value: $dzSpectError=0."
peculiarVelocityDispersionDamping::usage="Option for setting damping in redshift space induced by nonlinear peculiar velocity dispersion in km/s. Default value: $sigmaPecVel=0."
extraPowerShotNoiseFunction::usage="Option for including an extra shot noise for the observed power spectrum. Should be a function of z and k."
pks8Ratio::usage="Option for the power spectrum. If set to 1, then the power spectrum read from files is P(k,z)/s8(z)^2. Default: set to 0."



hubble::usage="Protected name of cosmological parameter h"
OmegaDE::usage="Protected name of cosmological parameter OmegaDE or OmegaL"
Omegam::usage="Protected name of cosmological parameter Omegam"
Omegab::usage="Protected name of cosmological parameter Omegab"
Omegac::usage="Protected name of cosmological parameter Omegac"
omegac::usage="Protected name of cosmological parameter omegac"
omegab::usage="Protected name of cosmological parameter omegab"
omegam::usage="Protected name of cosmological parameter omegam"

ns::usage="Protected name of cosmological parameter ns"
As::usage="Protected name of cosmological parameter As (Notice that there are several conventions for As (logAs, 10^9As, ...)"
sigma8::usage="Protected name of cosmological parameter sigma8 (normalization of the power spectrum)"
w0::usage="Protected name of cosmological parameter w0"
wa::usage="Protected name of cosmological parameter wa"
tau::usage="Protected name of cosmological parameter tau (optical depth)"

sigmapnl::usage="Protected name of cosmological parameter sigmapnl, used in the non-linear recipe. Enters the Lorentzian factor for the FoG."
sigmavnl::usage="Protected name of cosmological parameter sigmavnl, used in the non-linear recipe. Enters the BAO damping factor g_mu."

$Omeganu::usage="Global name of cosmological parameter Omeganu, containing a value of Omega for massive neutrinos or extra matter species"
$OmegaR::usage="Global name of cosmological parameter OmegaR, containing a value of Omega for radiation species"


$hubblereference::usage="Reference value for h (Hubble_0/100 today), when it is needed but it cannot be obtained from $paramoptions."
$sigma8reference::usage="Reference value for sigma_8 (at z=0), when it is needed but it cannot be obtained from $paramoptions."
$Asreference::usage="Reference value for A_s, when it is needed but it cannot be obtained from $paramoptions."


$funcsOfParameters::usage="Array of size Length@$paramoptions, that contains the functions that transform the fundamental parameters in $paramoptions."
$inverseFuncsOfParameters::usage="Array of size Length@$paramoptions, that contains the inverse functions that transform back the parameters in $paramoptions to fundamental parameters."


$cosmoParametersProtectedList::usage="gives the list of available cosmological parameters, whose names are protected and are available for their use in a Fisher forecast"

clearCosmologyMemoizedQuantities::usage="Clears all functions of cosmology, which have memoized values. Useful for computing many fiducials in a loop."


withCodeAfter::usage="Trick for protecting memoized functions"


$epsilonzstep::usage="Global value of epsilon step for numerical derivatives of z-dependent functions"


$debugPlots::usage="Boolean variable that activates Plots in debugging mode"

$zeroFunction::usage="Function equivalent to zero, used as a tool in certain specific cases."


Needs["cosmomathica`interface`"];
Needs["NumericalCalculus`"];
Needs["UsefulTools`"];
EndPackage[]


BeginPackage["CosmologyFunctions`", {"cosmomathica`interface`", "NumericalCalculus`","UsefulTools`"}]






Protect[hubble, OmegaDE, Omegam, Omegab, Omegac, omegac, omegab, omegam,
  omeganu, Omeganu, ns, As, sigma8, w0, wa, tau, sigmapnl, sigmavnl]


Begin["`Private`"]




$debugPlots=False;


$zeroFunction[x__] := 0

$cosmoParametersProtectedList={hubble, OmegaDE, Omegam, Omegab, Omegac,
  omegac, omegab, omegam, omeganu, Omeganu, ns, As, sigma8, w0, wa, tau, sigmapnl, sigmavnl}


$lightspeed=299792.458; (*N[QuantityMagnitude@UnitConvert[Quantity["SpeedOfLight"], "km/s"],9];*)


setParameterOptions[paramNames_,paramFiducials_,paramLabels_]:=Module[{chkli},
chkli=Head[paramNames];
If[chkli!=List, Print["Wrong parameter format passed"]];
$paramoptions=Thread[paramNames->paramFiducials];
$paramnames=paramNames;
$paramfidus=paramFiducials;
$parampositions=Thread[#->First@First@Position[$paramnames, #]]&/@$paramnames;
$ncosmopars=Length@$paramnames;
$paramlabels=paramLabels;
  Table[$funcsOfParameters[i] = None, {i,Length@$paramoptions}];
  Table[$inverseFuncsOfParameters[i] = None, {i,Length@$paramoptions}];
];

Options[setPScosmoOptions]={kdependentGrowth->False, lcdmBool->True, linearBool->True};

setPScosmoOptions[opts:OptionsPattern[]]:=Module[{chkli},
$pscosmoopts={lcdmBool->OptionValue[lcdmBool], linearBool->OptionValue[linearBool]};
SetOptions[setPScosmoOptions, $pscosmoopts];
$pslinearopt=FilterRules[$pscosmoopts,_?(StringMatchQ[SymbolName[#], "lin*"] &)];
$modcosmoopt=FilterRules[$pscosmoopts,_?(StringMatchQ[SymbolName[#], "lcdm*"] &)];
$kdependentGrowth=OptionValue[kdependentGrowth];

If[(lcdmBool/.$pscosmoopts) ,
If[(linearBool/.$pscosmoopts),$runtype="linear-LCDM",$runtype="nonlinear-LCDM"],
If[(linearBool/.$pscosmoopts),$runtype="linear-external-Cosmo",$runtype="nonlinear-external-Cosmo"]];  (*this should be fixed, so that it is clear that lcdm, linear is chosen without having to set options*)
$modelName="-"<>$runtype;           (*model name for file output save*)
Print["This run $modelname = ", $modelName];
];





Options[setORreadTransformedParameters]=$paramoptions~Join~{
generateInverseFunction->False, convertTransformedParameter->False, explicitTransformFunction->None}

setORreadTransformedParameters::symbolarg="Function setORreadTransformedParameters only accepts arguments with head Symbol"
setORreadTransformedParameters::param="Parameter `1` is not a String that contains a Symbolic function of a parameter or explicitFunction was not provided."
setORreadTransformedParameters::wrongopts="Options for transformation of a parameter are wrong or inconsistent, please check."

setORreadTransformedParameters[___]:=Message[setORreadTransformedParameters::symbolarg]

setORreadTransformedParameters[parsymb_Symbol,opts:OptionsPattern[]]:=Block[{paropts,invetra,parval,func,funcpass, parout,
  rulepar, readpar, readbool, posi, xs},
  readbool = OptionValue[convertTransformedParameter];
  invetra = OptionValue[generateInverseFunction];
  paropts=complementParamValues[{opts},setORreadTransformedParameters,returnList->"Full", filterCosmoPars->True];
  (*Print[paropts];*)

  rulepar = FilterRules[paropts,_?(StringMatchQ[ToString[#],"*"<>SymbolName[parsymb]<>"*"]&)];
  posi=First@First@(Position[paropts,First@rulepar]);
  parval=(First@rulepar)[[2]];

  If[SameQ[((First@rulepar)[[1]]), parsymb]==True,
    Return[parval]  (*In case the user forgot to unset the option transformedParameter, teh function still returns the correct value*)
  ];
  If[(invetra==True || (readbool==True && SameQ[$inverseFuncsOfParameters[posi], None]) ),
    Unset[$funcsOfParameters[posi]];
    Unset[$inverseFuncsOfParameters[posi]];
    Which[
      SameQ[OptionValue[explicitTransformFunction], None] && StringQ[((First@rulepar)[[1]])],
      func = ToExpression[((First@rulepar)[[1]])]/.{parsymb->xi};
      $funcsOfParameters[posi][xi_]=func;
      $inverseFuncsOfParameters[posi] = InverseFunction[$funcsOfParameters[posi]];
      ,
      UnsameQ[OptionValue[explicitTransformFunction], None],
      func = OptionValue[explicitTransformFunction];
      $funcsOfParameters[posi]=func;
      $inverseFuncsOfParameters[posi] = InverseFunction[$funcsOfParameters[posi]];
      ,
      True,
      Message[setORreadTransformedParameters::param, ((First@rulepar)[[1]])]
    ];
  ];
  If[readbool==False,
    Return[$inverseFuncsOfParameters[posi][Global`x]//Simplify[#,Element[Global`x,Reals]]&]
  ];
  If[readbool==True,
    Return[$inverseFuncsOfParameters[posi][parval]]
  ];
]

$dampingGrowthNorm=0.758;
$dampingSigma0Value=11.0;
$radiusScale=8;
$dampingSwitch=0;
$dzSpectError=0;
$sigmaPecVel=0;
$APeffectBool;
$zerofunc[x__]:=0;
$PowerExtraShot[zred_,kref_]:=0*zred*kref;
$APoptString="no-AP"
dzRelativeVelocityError = dzSpectroscopicError; (*duplicate option name, left here to avoid compatibility issues*)
$pks8ratio=0;
$shaParVarInZdep=False;

Options[setPowerSpectrumRelatedSpecs]={radius8MpcScale->$radiusScale,
  dampingGrowthNorm->$dampingGrowthNorm, dampingSigma0->$dampingSigma0Value,
dampingSwitch->$dampingSwitch, dzSpectroscopicError->$dzSpectError,
peculiarVelocityDispersionDamping->$sigmaPecVel,
APeffect->$APeffectBool, extraPowerShotNoiseFunction->$PowerExtraShot,
sigma8reference->$sigma8reference, Asreference->$Asreference, pks8Ratio->False,
zdepFunctionsCosmoVariation->$shaParVarInZdep}

setPowerSpectrumRelatedSpecs[opts:OptionsPattern[]]:=Block[{check(*performing consistency checks on options would be a good idea*),
pks8bool},
$dampingGrowthNorm=OptionValue[dampingGrowthNorm];
$dampingSigma0Value=OptionValue[dampingSigma0];
$radiusScale=OptionValue[radius8MpcScale];
$dampingSwitch=OptionValue[dampingSwitch];
$dzSpectError=OptionValue[dzSpectroscopicError];
$sigmaPecVel=OptionValue[peculiarVelocityDispersionDamping];
$sigma8reference=OptionValue[sigma8reference];
pks8bool=OptionValue[pks8Ratio];
$shaParVarInZdep=OptionValue[zdepFunctionsCosmoVariation];
$pks8ratio=If[pks8bool==True, 1, 0];
$APeffectBool=OptionValue[APeffect];
If[$APeffectBool==False, $APswitch=0; $APoptString="no-AP", $APswitch=1; $APoptString="yes-AP"];
$PowerExtraShot=OptionValue[extraPowerShotNoiseFunction];
SetOptions[setPowerSpectrumRelatedSpecs, radius8MpcScale->$radiusScale, dampingGrowthNorm->$dampingGrowthNorm, dampingSigma0->$dampingSigma0Value,
dampingSwitch->$dampingSwitch, dzSpectroscopicError->$dzSpectError,
peculiarVelocityDispersionDamping->$sigmaPecVel,
APeffect->$APeffectBool, extraPowerShotNoiseFunction->$PowerExtraShot,
sigma8reference->$sigma8reference,Asreference->$Asreference,
  pks8Ratio->pks8bool, zdepFunctionsCosmoVariation->$shaParVarInZdep
];
]

Options[complementParamValues]={globalPars->$paramoptions,returnList->"Complement",filterCosmoPars->False, debugOptions->False}
complementParamValues[extrapars_,functionsym_,opts:OptionsPattern[]]:=Block[{parext,pars,globpars,compars,posix,funcopts,filt},
funcopts=Options[functionsym];
parext=FilterRules[extrapars,funcopts];
compars=Complement[parext,funcopts];
globpars=OptionValue[globalPars];
filt=OptionValue[filterCosmoPars];
If[OptionValue[debugOptions]==True,
debugPrint[funcopts, 1];
debugPrint[parext, 1];
debugPrint[compars, 1];
];
Switch[OptionValue[returnList],
"Default"
  ,
pars=funcopts
  ,
"Complement"
  ,
pars=If[parext=={},{},compars],
"Fiducials"
  ,
If[parext=={},
pars={},
posix=Flatten[Position[funcopts[[All,1]],#]&/@compars[[All,1]]];
pars=funcopts[[posix]]
];
  ,
"Full"
  ,
  If[parext=={},
  pars=funcopts,
  posix=Flatten[Position[funcopts[[All,1]],#]&/@compars[[All,1]]];
  pars=ReplacePart[funcopts,Thread[posix->compars]]
  ]
];
If[filt,
pars=FilterRules[pars,globpars],
pars]
]


$fGamma=5.0/9.0;
(*$hubblereference=0.7;*)
(*$sigma8reference=0.8;*)
$Omeganu=0;
$H0 = 100/$lightspeed; (* for reference: 100/299792.458, in units of h/Mpc*)
$OmegaR = 0;
$Omegabfixed = 0; (*value of Omegab in case it wants to be fixed and not considered a varying cosmo parameter*)
(*$externalPkUnits="h/Mpc";
$internalPkUnits="h/Mpc";
$externalH0units="1/Mpc";
$externalDistanceUnits="Mpc";
$internalH0units="h/Mpc";
$internalDistanceUnits="Mpc/h";
*)

Options[setCosmologyFixedValues]={externalHubbleUnits:>"h/Mpc", externalDistanceUnits:>"Mpc/h", internalHubbleUnits:>"h/Mpc",
internalDistanceUnits:>"Mpc/h", internalPkUnits:>"h/Mpc", externalPkUnits->"h/Mpc",
fGammaFit->$fGamma, hubbleReference->$hubblereference, OmegaExtraMatterSpecies->$Omeganu, sigma8reference->$sigma8reference,
  Asreference->$Asreference,
OmegaRadiation->$OmegaR, OmegabReference->$Omegabfixed}

setCosmologyFixedValues[opts:OptionsPattern[]]:=Module[{gamma, H0u, D0u, intH0u, intD0u, pkU, intpkU},
$fGamma=OptionValue[fGammaFit];
H0u=OptionValue[externalHubbleUnits];
D0u=OptionValue[externalDistanceUnits];
intH0u=OptionValue[internalHubbleUnits];
intD0u=OptionValue[internalDistanceUnits];
pkU=OptionValue[externalPkUnits];
intpkU=OptionValue[internalPkUnits];
$Omeganu=If[SameQ[OptionValue[OmegaExtraMatterSpecies], None] || OptionValue[OmegaExtraMatterSpecies]==0,
            0, OptionValue[OmegaExtraMatterSpecies]];
$OmegaR = OptionValue[OmegaRadiation];
$fGamma=getgamma[$paramoptions];
$hubblereference=OptionValue[hubbleReference]; (*Reference h, when h is a derived parameter and is not included in $paramoptions*)
$sigma8reference=OptionValue[sigma8reference]; (*Reference sigma8, when sigma8 is a derived parameter and is not included in $paramoptions*)
$Asreference=OptionValue[Asreference]; (*Reference sigma8, when sigma8 is a derived parameter and is not included in $paramoptions*)
$Omegabfixed=OptionValue[OmegabReference];     (*Reference Omegab, when Omegab is a derived parameter and is not included in $paramoptions*)
$externalH0units=H0u;
$externalDistanceUnits=D0u;
$internalH0units=intH0u;  (*internal LCDM H(z) is usually in units of h/Mpc*)
$internalDistanceUnits=intD0u;(*internal LCDM D(z) is usually in units of Mpc/h*)
$externalPkUnits=pkU;
$internalPkUnits=intpkU;
  SetOptions[setCosmologyFixedValues, externalHubbleUnits:>$externalH0units, externalDistanceUnits:>$externalDistanceUnits,
    internalHubbleUnits:>$internalH0units, internalDistanceUnits:>$internalDistanceUnits,
    internalPkUnits:>$internalPkUnits, externalPkUnits:>$externalPkUnits,
    fGammaFit->$fGamma, hubbleReference->$hubblereference, OmegabReference->$Omegabfixed,
    OmegaExtraMatterSpecies->$Omeganu, OmegaRadiation->$OmegaR, sigma8reference->$sigma8reference, Asreference->$Asreference
  ];
];



$kminRef=N[2Pi/800];
$kmaxRef=1.5;
$kMaxInt=10.0;
$zMaxIntp=5.0;
$zMaxIntegral=3000;
$fkfix=False;


Options[setKandZMinMaxFixedValues]={kminIntegrations->$kminRef, kmaxIntegrations->$kmaxRef, kmaxInterpolations->$kMaxInt,
zmaxInterpolations->$zMaxIntp, zmaxIntegrations->$zMaxIntegral, fixGrowthScalek->$fkfix}


setKandZMinMaxFixedValues[opts:OptionsPattern[]]:=Module[{checks, hvalue, unitsfactor},
  hvalue=hubbleToday[];
  unitsfactor=compatibleUnits[$externalPkUnits,$internalPkUnits, hvalue, inverseDistanceUnits->True];
$kminRef=OptionValue[kminIntegrations]*unitsfactor;
$kmaxRef=OptionValue[kmaxIntegrations]*unitsfactor;
$kMaxInt=OptionValue[kmaxInterpolations]*unitsfactor;
$zMaxIntp=OptionValue[zmaxInterpolations];
$zMaxIntegral=OptionValue[zmaxIntegrations];
$fkfix=OptionValue[fixGrowthScalek];
SetOptions[setKandZMinMaxFixedValues, kminIntegrations->$kminRef, kmaxIntegrations->$kmaxRef, kmaxInterpolations->$kMaxInt,
zmaxInterpolations->$zMaxIntp, zmaxIntegrations->$zMaxIntegral, fixGrowthScalek->$fkfix];
]


Options[setExternalCosmoInterpolatingFunctions]={
  externalHubbleInput->False,
  externalHubbleDerivativesInput->False,
  externalDistanceInput->False,
  externalDistanceDerivativesInput->False,
  externalGrowthInput->False,
  externalGrowthDerivativesInput->False,
  externalGrowthRateInput->False,
  externalGrowthRateDerivativesInput->False,
  externalSigma8ofZInput->False,
  externalSigma8ofZDerivativesInput->False,
  externalfGrowthRateSigma8ofZInput->False,
  externalfGrowthRateSigma8ofZDerivativesInput->False,
  externalPowerSpectrumInput->False,
  externalPowerSpectrumDerivativesInput->False,
  externalPowerSpectrumNoWiggleInput->False,
  externalPowerSpectrumNoWiggleDerivativesInput->False,
  externalSigmaLensingInput->False,
  externalSigmaLensingDerivativesInput->False};


setExternalCosmoInterpolatingFunctions[opts:OptionsPattern[]]:=Block[{
  distopt=OptionValue[externalDistanceInput],
  distderivsopt=OptionValue[externalDistanceDerivativesInput],
  hubopt=OptionValue[externalHubbleInput],
  hubderivsopt=OptionValue[externalHubbleDerivativesInput],
  growopt=OptionValue[externalGrowthInput],
  growderivsopt=OptionValue[externalGrowthDerivativesInput],
  growrateopt=OptionValue[externalGrowthRateInput],
  growratederivsopt=OptionValue[externalGrowthRateDerivativesInput],
  s8ofzopt=OptionValue[externalSigma8ofZInput],
  s8ofzderivsopt=OptionValue[externalSigma8ofZDerivativesInput],
  fs8ofzopt=OptionValue[externalfGrowthRateSigma8ofZInput],
  fs8ofzderivsopt=OptionValue[externalfGrowthRateSigma8ofZDerivativesInput],
  psopt=OptionValue[externalPowerSpectrumInput],
  psderivsopt=OptionValue[externalPowerSpectrumDerivativesInput],
  pknwopt=OptionValue[externalPowerSpectrumNoWiggleInput],
  sigopt=OptionValue[externalSigmaLensingInput],
  sigderivsopt=OptionValue[externalSigmaLensingDerivativesInput]
},

If[UnsameQ[hubopt,False],
  $externalHubbleInterpolatingFunction=hubopt;
  SetOptions[setExternalCosmoInterpolatingFunctions,externalHubbleInput->hubopt];
  SetOptions[externalHubbleFunction, inputNumericalDerivatives->False];
  SetOptions[Hubble,externalFile->True];
];
If[UnsameQ[hubderivsopt,False],
  $externalHubbleDerivativesInterpolatingFunction=hubderivsopt;
  SetOptions[setExternalCosmoInterpolatingFunctions,externalHubbleDerivativesInput->hubderivsopt];
  SetOptions[externalHubbleFunction, inputNumericalDerivatives->True];
  SetOptions[Hubble,externalFile->True];
];
If[UnsameQ[distopt,False],
  $externalDistanceInterpolatingFunction=distopt;
  SetOptions[setExternalCosmoInterpolatingFunctions,externalDistanceInput->distopt];
  SetOptions[externalDistanceFunction, inputNumericalDerivatives->False];
  SetOptions[distanceAngular,externalFile->True];
  SetOptions[comovingDistance,externalFile->True];
];
If[UnsameQ[distderivsopt,False],
  $externalDistanceDerivativesInterpolatingFunction=distderivsopt;
  SetOptions[setExternalCosmoInterpolatingFunctions,externalDistanceDerivativesInput->distderivsopt];
  SetOptions[externalDistanceFunction, inputNumericalDerivatives->True];
  SetOptions[distanceAngular,externalFile->True];
  SetOptions[comovingDistance,externalFile->True];
];
If[UnsameQ[growopt,False],
  $externalGrowthInterpolatingFunction=growopt;
  SetOptions[setExternalCosmoInterpolatingFunctions,externalGrowthInput->growopt];
  SetOptions[Growth,externalFile->True];
  SetOptions[Growth,growthRateIntegral->False];
  ,
  SetOptions[Growth,growthRateIntegral->True]
];
If[UnsameQ[growderivsopt,False],
  $externalGrowthDerivativesInterpolatingFunction=growderivsopt;
  SetOptions[setExternalCosmoInterpolatingFunctions,externalGrowthDerivativesInput->growderivsopt];
  SetOptions[Growth,externalFile->True];
  SetOptions[Growth,growthRateIntegral->False];
  SetOptions[externalGrowthFunction, inputNumericalDerivatives->True];
];
If[UnsameQ[growrateopt,False],
  $externalGrowthRateInterpolatingFunction=growrateopt;
  SetOptions[setExternalCosmoInterpolatingFunctions,externalGrowthRateInput->growrateopt];
  SetOptions[fGrowthRate,externalFile->True];
  SetOptions[fGrowthRate,growthDerivative->False];
  ,
  SetOptions[fGrowthRate,growthDerivative->True]
];
If[UnsameQ[growratederivsopt,False],
  $externalGrowthRateDerivativesInterpolatingFunction=growratederivsopt;
  SetOptions[setExternalCosmoInterpolatingFunctions,externalGrowthRateDerivativesInput->growratederivsopt];
  SetOptions[fGrowthRate,externalFile->True];
  SetOptions[fGrowthRate,growthDerivative->False];
  SetOptions[externalGrowthRateFunction, inputNumericalDerivatives->True];
  ];
If[UnsameQ[s8ofzopt,False],
  $externalSigma8ofZInterpolatingFunction=s8ofzopt;
  SetOptions[setExternalCosmoInterpolatingFunctions, externalSigma8ofZInput->s8ofzopt];
  SetOptions[sigma8ofZ,externalFile->True];
];
If[UnsameQ[s8ofzderivsopt,False],
  $externalSigma8ofZDerivativesInterpolatingFunction=s8ofzderivsopt;
  SetOptions[setExternalCosmoInterpolatingFunctions, externalSigma8ofZDerivativesInput->s8ofzderivsopt];
  SetOptions[sigma8ofZ,externalFile->True];
  SetOptions[externalSigma8ofZFunction, inputNumericalDerivatives->True];
];
If[UnsameQ[fs8ofzopt,False],
  $externalfGrowthRateSigma8ofZInterpolatingFunction=fs8ofzopt;
  SetOptions[setExternalCosmoInterpolatingFunctions, externalfGrowthRateSigma8ofZInput->fs8ofzopt];
  SetOptions[fGrowthRateSigma8ofZ,externalFile->True];
];
If[UnsameQ[fs8ofzderivsopt,False],
  $externalfGrowthRateSigma8ofZDerivativesInterpolatingFunction=fs8ofzderivsopt;
  SetOptions[setExternalCosmoInterpolatingFunctions, externalfGrowthRateSigma8ofZDerivativesInput->fs8ofzderivsopt];
  SetOptions[fGrowthRateSigma8ofZ,externalFile->True];
  SetOptions[externalfGrowthRateSigma8ofZFunction, inputNumericalDerivatives->True];
];
If[UnsameQ[sigopt,False],
  $externalSigmaLensingInterpolatingFunction=sigopt;
  SetOptions[setExternalCosmoInterpolatingFunctions, externalSigmaLensingInput->sigopt];
  SetOptions[SigmaLensing,externalFile->True];
];
If[UnsameQ[sigderivsopt,False],
  $externalSigmaLensingDerivativesInterpolatingFunction=sigderivsopt;
  SetOptions[setExternalCosmoInterpolatingFunctions, externalSigmaLensingDerivativesInput->sigderivsopt];
  SetOptions[SigmaLensing, externalFile->True];
  SetOptions[externalSigmaLensingFunction, inputNumericalDerivatives->True];
];
If[UnsameQ[psopt,False],
  $externalPowerSpectrumInterpolatingFunction=psopt;
  SetOptions[setExternalCosmoInterpolatingFunctions,externalPowerSpectrumInput->psopt];
  SetOptions[powerSpectrum, externalFile->True, spectrumMethod->"ExternalFile"];
  SetOptions[observedPowerSpectrum, externalFile->True, spectrumMethod->"ExternalFile"];
  SetOptions[externalPowerSpectrumFunction, inputNumericalDerivatives -> False, externalFunction->True];
];
If[UnsameQ[psderivsopt,False],
  $externalPowerSpectrumDerivativesInterpolatingFunction=psderivsopt;
  SetOptions[setExternalCosmoInterpolatingFunctions, externalPowerSpectrumDerivativesInput->psderivsopt];
  SetOptions[powerSpectrum, externalFile->True, spectrumMethod->"ExternalFile"];
  SetOptions[observedPowerSpectrum, externalFile->True, spectrumMethod->"ExternalFile"];
  SetOptions[externalPowerSpectrumFunction, inputNumericalDerivatives -> True,
    externalFunction->psopt, externalDerivativeFunction->psderivsopt];
];
If[UnsameQ[pknwopt,False],
  $externalPowerSpectrumNoWiggleInterpolatingFunction=pknwopt;
  SetOptions[setExternalCosmoInterpolatingFunctions,externalPowerSpectrumNoWiggleInput->pknwopt];
  SetOptions[externalPowerSpectrumFunction,externalPowerSpectrumNoWiggleInput->pknwopt];
  SetOptions[powerSpectrum, externalFile->True, spectrumMethod->"ExternalFile"];
  SetOptions[observedPowerSpectrum, externalFile->True, spectrumMethod->"ExternalFile"]; (* ExternalFile option here, since ExternalNoWiggle is just a specific sub-case for the non-linear recipe *)
  (*SetOptions[externalPowerSpectrumFunction, inputNumericalDerivatives-> False, externalFunction->True]; No need for this option, external NoWiggle spectrum is only used jointly with standard Pk *)
];
If[UnsameQ[pknwderivsopt,False],
  $externalPowerSpectrumNoWiggleDerivativesInterpolatingFunction=pknwderivsopt;
  SetOptions[setExternalCosmoInterpolatingFunctions,externalPowerSpectrumNoWiggleDerivativesInput->pknwopt];
  SetOptions[externalPowerSpectrumFunction,externalPowerSpectrumNoWiggleDerivativesInput->pknwderivsopt];
  SetOptions[powerSpectrum, externalFile->True, spectrumMethod->"ExternalFile"];
  SetOptions[observedPowerSpectrum, externalFile->True, spectrumMethod->"ExternalFile"]; (* ExternalFile option here, since ExternalNoWiggle is just a specific sub-case for the non-linear recipe *)
  (*SetOptions[externalPowerSpectrumFunction, inputNumericalDerivatives-> False, externalFunction->True]; No need for this option, external NoWiggle spectrum is only used jointly with standard Pk *)
];
]


getCosmoParameterFileName::fiduwarn="Warning: first element `1` in $paramsDirectoryNames does not contain the string 'fiducial'"

Options[getCosmoParameterFileName]=$paramoptions~Join~$modcosmoopt~Join~{parameterDirectoriesNames->$paramsDirectoryNames,
inputNumericalDerivatives->False, getNonFiducialFileName->False, getFiducialFileName->False}

getCosmoParameterFileName[func_,paroptions_,funcopts:OptionsPattern[]]:=Block[{ran,pardirs,parextraopt, parfidu, fiduval, parextval,
diffratio, tole, filestring="None", index,numderivsin},
  tole=10*^-6; (*tolerance should be smaller than the epsilon of the file parameters, but not too small to avoid problems*)
parextraopt=complementParamValues[{paroptions},func,returnList->"Complement", filterCosmoPars->True];
pardirs=OptionValue[parameterDirectoriesNames];
numderivsin=OptionValue[inputNumericalDerivatives];

  If[StringMatchQ[pardirs[[1]],"Fiducial*",IgnoreCase->True]==False,Message[getCosmoParameterFileName::fiduwarn, pardirs[[1]] ]];

  If[OptionValue[getFiducialFileName]==True,
    filestring=pardirs[[1]];
    Return[filestring]
  ];

  If[OptionValue[getNonFiducialFileName]==True,
  ran=RandomInteger[{1, (Length@$paramsDirectoryNames - 1)}];
  index=Flatten[Position[!StringMatchQ[#, "Fiducial",IgnoreCase -> True] & /@ pardirs, True]][[ran]];
  filestring=pardirs[[index]];
  Return[filestring]
];
(*debugPrint[pardirs];*)
Which[numderivsin==False,
  If[Length@parextraopt==0
    ,
     filestring=pardirs[[1]];
    ,
    parextraopt=First@parextraopt;
    parfidu=First@complementParamValues[{paroptions},func,returnList->"Fiducials", filterCosmoPars->True];
    parextval=Last@parextraopt;
    fiduval=Last@parfidu;
    index=Position[$paramoptions,First@parextraopt][[1,1]];
    If[fiduval != 0,
        diffratio=(parextval/fiduval);
        ,
        diffratio = 1 + (parextval);  (*for absolute epsilons used when fiducial param == 0)*)
        ];
    tole=Abs[tole];
    debugPrint["diffratio: "<> ToString[diffratio], 1];
     Which[Chop[Abs[1-diffratio],tole]==0,
      filestring=pardirs[[1]];,
      diffratio < 1,  (*corresponds to (1-eps)*)
      filestring=pardirs[[2*index]];,
      diffratio > 1,   (*corresponds to (1+eps)*)
      filestring=pardirs[[2*index+1]];
     ];
  ];
  ,
  numderivsin==True,
  parextraopt=First@parextraopt;
  index=Position[$paramoptions,First@parextraopt][[1,1]];
  filestring=pardirs[[1+index]];
];
  debugPrint["filename: "<>filestring];
  Return[filestring]
]


aOfz[z_]:=1/(1+z);

zOfa[a_]:=(1/a)-1;

NfOfz[z_]:=Log[aOfz[z]];

zOfNf[nf_]:=Exp[-nf]-1;



wDEzFunc[ztt_, w0_, wa_]=Assuming[ztt\[Element]Reals && ztt>0,Integrate[3(1+w0+wa*zx/(1+zx))/(1+zx),{zx,0,ztt}]];

Options[wEoSparameters]=$paramoptions;

wEoSparameters[opts:OptionsPattern[]]:=Block[{return, paropts, w0r, war},
paropts=complementParamValues[{opts},wEoSparameters,returnList->"Full", filterCosmoPars->True];
If[FilterRules[paropts,w0]!={},w0r=OptionValue[w0],w0r=-1];
If[FilterRules[paropts,wa]!={},war=OptionValue[wa],war=0];
return={w0r,war}
]

Options[parametersModifiedGravity]:=$paramoptions~Join~{inverseTransformation->(Identity[#]&)};
parametersModifiedGravity[parname_,opts:OptionsPattern[]]:=Block[{return, paropts, par1, par1T, inv},
inv=OptionValue[inverseTransformation];
paropts=complementParamValues[{opts},parametersModifiedGravity,returnList->"Full", filterCosmoPars->True];
If[FilterRules[paropts,parname]!={},
  par1=OptionValue[parname],
  par1=0
];
par1T=inv@par1;  (*Apply inverse transformation. Used for example to vary log(param) but passing param to Cosmomathica CAMB*)
Return[par1T]
]


curvatureSin[x_,ok_]:=Which[ok>.001,Sinh[x],ok<-.001,Sin[x],Abs[ok]<=.001,x];

curvatureFunc=Function[{xi,ki},Limit[Sinh[Sqrt[-k]*xi]/Sqrt[-k],k->ki]];

Options[unity]=$paramoptions;
unity[opts:OptionsPattern[]]:=Block[{OmR=$OmegaR},
  1 - OmR   (* extra matter species like Omeganu are taken into account in OmegaCDM and OmegaBaryon*)
]; (*this function returns 1 minus the fixed parameters contributing to energy density. In principle accepts options, but doesn't use them*)


Options[OmegaBaryon0Today]=$paramoptions;
OmegaBaryon0Today[opts:OptionsPattern[]]:=Block[{OmegaB0,OmegaC0,OmegaDE0,paropts, hubref, uno},
paropts=complementParamValues[{opts},OmegaBaryon0Today,returnList->"Full", filterCosmoPars->True];
hubref=hubbleToday[paropts];
uno = unity[];
OmegaB0=Which[
  FilterRules[paropts,Omegab]!={},OptionValue[Omegab],
  Length[FilterRules[paropts,{omegab}]]==1,OptionValue[omegab]/(hubref^2),
  Length[FilterRules[paropts,{Omegam,Omegac}]]==2,(OptionValue[Omegam]-OptionValue[Omegac]-$Omeganu),
  Length[FilterRules[paropts,{Omegac,OmegaDE}]]==2,(uno-OptionValue[Omegac]-OptionValue[OmegaDE]),
  Length[FilterRules[paropts,{omegam,omegac,hubble}]]==3,(OptionValue[omegam]-OptionValue[omegac])/(hubref^2)-$Omeganu,
  Length[FilterRules[paropts,{omegac,omegab,OmegaDE}]]==3,(OptionValue[omegac]*(uno-OptionValue[OmegaDE])/(OptionValue[omegac]+OptionValue[omegab])),
  FilterRules[paropts,Omegab]=={} && FilterRules[paropts,Omegac]=={} && (FilterRules[paropts,Omegam]!={} || FilterRules[paropts,omegam]!={}) , $Omegabfixed
  ]
];

Options[OmegaCDM0Today]=$paramoptions;
OmegaCDM0Today[opts:OptionsPattern[]]:=Block[{OmegaC0,paropts, hubref, OmegaB0, OmegaDE0, uno},
paropts=complementParamValues[{opts},OmegaCDM0Today,returnList->"Full", filterCosmoPars->True];
uno = unity[];
hubref=hubbleToday[paropts];
OmegaB0=OmegaBaryon0Today[paropts];
OmegaC0=Which[
  FilterRules[paropts,Omegac]!={},OptionValue[Omegac],
  Length[FilterRules[paropts,{omegac}]]==1,OptionValue[omegac]/(hubref^2),
  Length[FilterRules[paropts,{Omegam,Omegab}]]==2,(OptionValue[Omegam]-OmegaB0-$Omeganu),
  Length[FilterRules[paropts,{Omegab,OmegaDE}]]==2,(uno-OmegaB0-OptionValue[OmegaDE]),
  Length[FilterRules[paropts,{omegab,OmegaDE}]]==2,(uno-OmegaB0-OptionValue[OmegaDE]),
  Length[FilterRules[paropts,{omegam,omegab,hubble}]]==3,(OptionValue[omegam]-OptionValue[omegab])/(hubref^2)-$Omeganu,
  Length[FilterRules[paropts,{omegac,omegab,OmegaDE}]]==3,(OptionValue[omegac]*(uno-OptionValue[OmegaDE])/(OptionValue[omegac]+OptionValue[omegab])),
  OmegaB0!=0 && FilterRules[paropts,Omegac]=={} && (FilterRules[paropts,Omegam]!={}), (OptionValue[Omegam]-OmegaB0-$Omeganu),
  OmegaB0!=0 && FilterRules[paropts,Omegac]=={} && (FilterRules[paropts,omegam]!={}), (OptionValue[omegam]/hubref^2-OmegaB0-$Omeganu),
  OmegaB0==0 && FilterRules[paropts,Omegac]=={} && FilterRules[paropts,Omegam]!={}, OptionValue[Omegam]
  ]
];

Options[OmegaDE0Today]=$paramoptions;
OmegaDE0Today[opts:OptionsPattern[]]:=Block[{OmegaC0,paropts, hubref, OmegaB0, OmegaM0,OmegaDE0, uno},
paropts=complementParamValues[{opts},OmegaDE0Today,returnList->"Full", filterCosmoPars->True];
uno = unity[];
hubref=hubbleToday[paropts];
OmegaM0=OmegaM0Today[paropts];
OmegaDE0=Which[
  FilterRules[paropts,OmegaDE]!={},OptionValue[OmegaDE],
  True,
  uno-OmegaM0
  ]
];

Options[OmegaM0Today]=$paramoptions;
OmegaM0Today[opts:OptionsPattern[]]:=Block[{OmegaM0,hubref,OmegaDE0, OmegaC0, OmegaB0,paropts, uno},
paropts=complementParamValues[{opts},OmegaM0Today,returnList->"Full", filterCosmoPars->True];
uno = unity[];
OmegaB0=OmegaBaryon0Today[paropts];
hubref=hubbleToday[paropts];
OmegaC0=OmegaCDM0Today[paropts];
(*If[OmegaM0!=(OmegaB0+OmegaC0+$Omeganu), Print["Omegam0 does not match Omegab0, Omegac0 and Omeganu0"; Abort[]]];*)
OmegaM0=Which[
  FilterRules[paropts,Omegam]!={}, OptionValue[Omegam],
  Length[FilterRules[paropts,{omegam, hubble}]]==2,(OptionValue[omegam])/(hubref^2),
  Length[FilterRules[paropts,{OmegaDE}]]==1, uno-OptionValue[OmegaDE],
  True,
  OmegaB0+OmegaC0
  ]
]


Options[hubbleToday]=$paramoptions;
hubbleToday[opts:OptionsPattern[]]:=Block[{paropts, hubblevalue, uno},
paropts=complementParamValues[{opts},hubbleToday,returnList->"Full", filterCosmoPars->True];
uno = unity[];
hubblevalue=Which[
  FilterRules[paropts,hubble]!={},OptionValue[hubble],
  Length[FilterRules[paropts,{omegam,OmegaDE}]]==2,Sqrt[OptionValue[omegam]/(uno-OptionValue[OmegaDE])],
  Length[FilterRules[paropts,{omegac,omegab,OmegaDE}]]==3,Sqrt[(OptionValue[omegab]+OptionValue[omegac])/(uno-OptionValue[OmegaDE])],
  FilterRules[paropts,hubble]=={},$hubblereference
  ]
];



Options[OmegaK0Today]=$paramoptions;

OmegaK0Today[opts:OptionsPattern[]]:=
Block[{OmegaDE0,hubref,OmegaC0, OmegaB0, w0ref,waref, paropts, flatness, tole=0.0001, curv},
paropts=complementParamValues[{opts},OmegaK0Today,returnList->"Full", filterCosmoPars->True];
OmegaDE0=OmegaDE0Today[paropts];
hubref=hubbleToday[paropts];
OmegaB0=OmegaBaryon0Today[paropts];
OmegaC0=OmegaCDM0Today[paropts];

flatness = $OmegaR + OmegaB0 + OmegaC0 + $Omeganu + OmegaDE0;
curv=Chop[1-flatness, tole];

Return[curv]
];




LCDMDimensionlessHubble::flatness="Warning: The sum of the energy density fractions (Omega_i) at z=0 is not 1."

Options[LCDMDimensionlessHubble]=$paramoptions
LCDMDimensionlessHubble[z_,opts:OptionsPattern[]]:=LCDMDimensionlessHubble[z,opts]=
Block[{OmegaDE0,hubref,OmegaC0, OmegaB0,Omegak, w0ref,waref, paropts, flatness},
paropts=complementParamValues[{opts},LCDMDimensionlessHubble,returnList->"Complement", filterCosmoPars->True];
OmegaDE0=OmegaDE0Today[paropts];
hubref=hubbleToday[paropts];
OmegaB0=OmegaBaryon0Today[paropts];
OmegaC0=OmegaCDM0Today[paropts];
{w0ref,waref}=wEoSparameters[paropts];
Omegak=OmegaK0Today[paropts];

(*If[ UnsameQ[Chop[flatness-1, 10^-8], 0], Message[LCDMDimensionlessHubble::flatness] ];*)

Sqrt[$OmegaR*(1+z)^4 + OmegaB0*(1+z)^3 + OmegaC0*(1+z)^3 + $Omeganu*(1+z)^3 + OmegaDE0*Exp[wDEzFunc[z,w0ref,waref]]+ (Omegak)*(1+z)^2]
];

(*Example on how to protect a memoized function*)


(*ClearAll[withCodeAfter];*)
SetAttributes[withCodeAfter, HoldRest];
withCodeAfter[before_, after_] := (after; before);


(*Unprotect[LCDMDimensionlessHubble];
ClearAll[LCDMDimensionlessHubble];*)


(*call_LCDMDimensionlessHubble /; MemberQ[Attributes[LCDMDimensionlessHubble], Protected] :=
   withCodeAfter[Unprotect[LCDMDimensionlessHubble]; call, Protect[LCDMDimensionlessHubble]]*)


(*Protect[LCDMDimensionlessHubble]*)


(*End of Example on how to protect a memoized function*)


testUnitString::invunit="The unit expression `1` is not a valid unit for this quantity. Aborting.";
testUnitString[ax_]:=Block[{h=Global`h,Mpc=Global`Mpc,km=Global`km,s=Global`s},
  Which[
  ax===1/Mpc,
  "1/Mpc",
  ax===Mpc/h,
  "Mpc/h",
  ax===h/Mpc,
  "h/Mpc",
  ax===Mpc,
  "Mpc",
  ax===((km/s)/Mpc),
  "km/s/Mpc",
  ax===(Mpc/(km/s)),
  "Mpc/(km/s)",
  True,
  Message[testUnitString::invunit,ax];
  Abort[]
  ]
]


Options[compatibleUnits]=$paramoptions~Join~{invertUnits->True, inverseDistanceUnits->False, physicalUnits->False};

compatibleUnits[unitstr_String, opts:OptionsPattern[]]:=Block[{compunit, expu},
  expu=ToExpression[unitstr];
  If[OptionValue[invertUnits]==True,
    expu = 1/expu;
  ];
  compunit=testUnitString[expu];
  Return[compunit]
];


compatibleUnits[sourceunit_String,targetunit_String, opts:OptionsPattern[]]:=Block[{invopt=OptionValue[invertUnits], sunit, unitsfact, hval, paropts},
       paropts=complementParamValues[{opts},Hubble,returnList->"Full", filterCosmoPars->True];
       hval=hubbleToday[paropts];
       If[invopt==True,
         sunit=compatibleUnits[sourceunit, invertUnits->invopt],
         sunit=sourceunit
         ];
  ];

compatibleUnits[sourceunit_String,targetunit_String, hval_Real, opts:OptionsPattern[]]:=Block[{invunit, unitsfact},
  Which[
    OptionValue[inverseDistanceUnits]==False,
    Switch[targetunit,
      "Mpc",unitsfact=Which[sourceunit=="Mpc",1,sourceunit=="Mpc/h",(1/hval)],
      "Mpc/h",unitsfact=Which[sourceunit=="Mpc",(hval), sourceunit=="Mpc/h", 1]
    ],
    OptionValue[inverseDistanceUnits]==True,
    Switch[targetunit,
      "h/Mpc",unitsfact=Which[sourceunit=="1/Mpc", (1/hval), sourceunit=="h/Mpc", 1, sourceunit=="km/s/Mpc", (1/($lightspeed*hval))];,
      "1/Mpc",unitsfact=Which[sourceunit=="1/Mpc", 1, sourceunit=="h/Mpc", (hval), sourceunit=="km/s/Mpc", (1/$lightspeed)];
    ];
  ];
  Return[unitsfact]
]




Options[externalUnitsConverter]=$paramoptions~Join~{inverseDistanceUnits->False, physicalUnits->False, internalUnits->True};

externalUnitsConverter[unitstr_String,hval_,opts:OptionsPattern[]]:=Block[{pars, unitsopt, unitsfact},
unitsopt=unitstr;
debugPrint["inside 1"];
debugPrint[OptionValue[internalUnits]];
debugPrint[OptionValue[inverseDistanceUnits]];
Which[
  OptionValue[internalUnits]==True,
  Which[
    OptionValue[inverseDistanceUnits]==True,
      If[OptionValue[physicalUnits]==True,unitsopt="km/s/Mpc"];
      debugPrint["inside hubble"];
      Switch[unitsopt, (*Outputs internal Hubble function in [h/Mpc], [1/Mpc] or [km/s/Mpc]*)
        "h/Mpc", unitsfact=(100/$lightspeed),
        "1/Mpc", unitsfact=(100/$lightspeed)*hval,
        "km/s/Mpc",unitsfact=100*hval
      ],
    OptionValue[inverseDistanceUnits]==False,
      Switch[unitsopt, (*Outputs internal Distance function in [Mpc] or [Mpc/h]*)
       "Mpc",unitsfact=$lightspeed/(100*hval),
       "Mpc/h",unitsfact=$lightspeed/100;
      ]
  ],
  OptionValue[internalUnits]==False,
      debugPrint["inside external"];

      Which[
        OptionValue[inverseDistanceUnits]==False,
          Switch[unitsopt, (*Converts external units of Distance function into $internalDistanceUnits*)
            "Mpc",unitsfact=Which[$internalDistanceUnits=="Mpc",1,$internalDistanceUnits=="Mpc/h",(hval)],
            "Mpc/h",unitsfact=Which[$internalDistanceUnits=="Mpc",(1/hval),$internalDistanceUnits=="Mpc/h",1]
          ],
        OptionValue[inverseDistanceUnits]==True,
          Switch[unitsopt, (*Converts external units of Hubble function into $internalH0Units*)
            "h/Mpc",unitsfact=Which[$internalH0units=="1/Mpc",hval,$internalH0units=="h/Mpc",1];,
            "1/Mpc",unitsfact=Which[$internalH0units=="1/Mpc",1,$internalH0units=="h/Mpc",(1/hval)];,
            "km/s/Mpc",unitsfact=Which[$internalH0units=="1/Mpc",(1/$lightspeed),$internalH0units=="h/Mpc",(1/($lightspeed*hval))]
          ];
          If[OptionValue[physicalUnits]==True,
            unitsfact=unitsfact*Which[$internalH0units=="1/Mpc",($lightspeed),$internalH0units=="h/Mpc",($lightspeed*hval)]
          ];
        ]
    ];
 Return[unitsfact]
]




ClearAll[externalFunctionTaylor];

Options[externalFunctionTaylor]=$paramoptions~Join~{functionInterpolatedInLogLog->False,derivativeInterpolatedInLogLin->False,
  externalFunction->True,externalDerivativeFunction->True, interpolatedArguments->"k", kdependence->True, externalSpectraType -> "pk"};

externalFunctionTaylor::extfuncs="Option `1` was not properly set to a Function."

externalFunctionTaylor[zkarg__,opts:OptionsPattern[]]:=Module[{zz, kk, kdep, par,eps,paropts,parfidu,filename,
  fidu,return,intpDerlinlog,intpFunloglog,funcy=(#)&,funcx=(#)&,funcyprime=(#)&,funcxprime=(#)&,
  fiduval,parval,Delta,tolerance=10^(-6),
(*hard coded tolerance to distinguish between fiducial and value away from fiducial*)
  intpargs,
  extFunction,extFunctionDerivative,
  externalFuncTemplate, externalFuncDerivTemplate, parname, parlab, parsrules,
  externalFuncResult, externalFuncDerivResult,
  maxkdomain},

  parsrules=Thread[Rule[$paramnames,$paramlabels]];
  extFunction=OptionValue[externalFunction];
  extFunctionDerivative=OptionValue[externalDerivativeFunction];
  Which[BooleanQ[extFunction],
     Message[externalFunctionTaylor::extfuncs, "extFunction"]; Abort[]
     ,
     BooleanQ[extFunctionDerivative],
     Message[externalFunctionTaylor::extfuncs, "extFunctionDerivative"]; Abort[]
     ];

  kdep=OptionValue[kdependence];

  {zz,kk}=twoArgumentsCheck[zkarg,kdep];
  kk=ReleaseHold[kk];

  paropts=complementParamValues[{opts},externalFunctionTaylor,returnList->"Complement",filterCosmoPars->True];
  (*fidu=getCosmoParameterFileName[externalFunctionTaylor,paropts, getFiducialFileName->True];*)
  intpFunloglog=OptionValue[functionInterpolatedInLogLog];
  intpDerlinlog=OptionValue[derivativeInterpolatedInLogLin];
  If[intpFunloglog==True,funcy=(10^#)&;funcx=Log10[#]&;];
  If[intpDerlinlog==True,funcyprime=(#)&;funcxprime=Log10[#]&;];

  (*If[StringMatchQ[fidu,"Fiducial*",IgnoreCase->True]==False,Message[getCosmoParameterFileName::fiduwarn,fidu]];*)
  debugPrint["paropts: "<>ToString@paropts,1];

  If[Length@paropts>=1,
    paropts=First@paropts;
    parfidu=First@complementParamValues[{paropts},externalFunctionTaylor,returnList->"Fiducials",filterCosmoPars->True];
    fiduval=Last@parfidu;
    parval=Last@paropts;
    (*filename=getCosmoParameterFileName[externalFunctionTaylor,paropts,inputNumericalDerivatives->True];*)
    parname=First@paropts;
    Delta=Chop[(parval-fiduval),tolerance];
    , (*Else*)
    parname=First@First@$paramoptions; (*get a parameter name that is not "fiducial", to avoid finding a wrong name in derivative file"*)
    Delta=0;
  ];

    parlab=parname/.parsrules;

  intpargs=OptionValue[interpolatedArguments];

  Which[
    intpargs=="zk",
    maxkdomain = Max@extFunction["Domain"][[2]];
    externalFuncResult=functionExtrapolationFit[kk, extFunction[zz,#]&, maxkdomain];
    externalFuncDerivResult=If[SameQ[extFunctionDerivative[parlab], $zeroFunction]==True,
    0,(*for the case in which the derivative is set to be zero at the input*)
    (*ELSE*)
    maxkdomain = Max@extFunctionDerivative[parlab]["Domain"][[2]];
    functionExtrapolationFit[kk, extFunctionDerivative[parlab][zz,#]&, maxkdomain]
    ];
    ,
    intpargs=="k",
    (*TODO: Implement extrapolation Fit here*)
    externalFuncResult=extFunction[[zValTozIndex[zaa]]][funcx@kaa];
    externalFuncDerivResult=If[SameQ[extFunctionDerivative[parlab][[zValTozIndex[zaa]]], $zeroFunction]==True,
    0,(*for the case in which the derivative is set to be zero at the input*)
    (*ELSE*)
    extFunctionDerivative[parlab][[zValTozIndex[zaa]]][funcxprime@kaa]
    ];
    ,
    intpargs=="z",
    externalFuncResult=extFunction[zz];
    externalFuncDerivResult=extFunctionDerivative[parlab][zz];
  ];

  debugPrint["fidu val: "<>ToString@fiduval];
  debugPrint["param val: "<>ToString@parval];
  debugPrint["param label: "<>ToString@parlab];
  debugPrint["Delta: "<>ToString@Delta];

  return=funcy@(externalFuncResult)+Delta*funcyprime@(externalFuncDerivResult);

  Return[return]
]


powerFit[Funcaofz_, doma_]:=Block[{domrangmin=0.50, domrangstep=20,
  pktab,lmfit,retu, klist, fittedPk,xa},
klist=Range[doma*domrangmin,doma,doma/domrangstep];
pktab=Table[{kkk,Funcaofz[kkk]}, {kkk,klist}];
lmfit=NonlinearModelFit[pktab, a*xa^b, {a,b},xa];
(*fittedPk[xa]= ((Exp[lmfit[xa]]/.xa->Log[xa]));  for the previous fit in log-log space*)
fittedPk[xa]= lmfit[xa];
Return[fittedPk[xa]]
]

functionExtrapolationFit[kako_, Funcaofz_, maxdomain_]:=Block[{extra,kval,nomcoeff,retu,doma, tole=0.98},
If[kako > tole*maxdomain,
retu= ((powerFit[Funcaofz, maxdomain])/.(xa->kako))
, (*ELSE*)
retu=Funcaofz[kako]
];
Return[retu]
]



Options[externalHubbleFunction]=$paramoptions~Join~$modcosmoopt~Join~{shiftedDomain->False, inputNumericalDerivatives->False};

externalHubbleFunction[zz_, opts:OptionsPattern[]]:=Block[{cosmopts,filename,zdom=0, inpunumder, exthub},

  inpunumder=OptionValue[inputNumericalDerivatives];

  cosmopts=complementParamValues[{opts}, externalHubbleFunction, returnList->"Complement", filterCosmoPars->True];

  If[OptionValue[shiftedDomain]==True, zdom=$externalHubbleInterpolatingFunction[filename]["Domain"][[1,1]] ];


If[inpunumder==False,
  filename=getCosmoParameterFileName[externalHubbleFunction,{opts}];
  exthub=$externalHubbleInterpolatingFunction[filename][zz+zdom];
  ,
  exthub=(externalFunctionTaylor[zz, cosmopts~Join~{
                        externalFunction->$externalHubbleInterpolatingFunction,
                        externalDerivativeFunction->$externalHubbleDerivativesInterpolatingFunction,
                        interpolatedArguments->"z", kdependence->False}])
  ];
Return[exthub]
]


Options[Hubble]=$paramoptions~Join~$modcosmoopt~Join~{physicalUnits->False, externalFile->False}

Hubble[z_?NumericQ,opts:OptionsPattern[]]:=Block[{lcdmOpt,paropts, parcosmopts, hvalue, Hubfunc, units, conversionfactor,physunits},
parcosmopts=complementParamValues[{opts},Hubble,returnList->"Complement", filterCosmoPars->True];
lcdmOpt=OptionValue[lcdmBool];  (*this should be fixed, so that it is clear that lcdm, linear is chosen without having to set options*)
debugPrint["Hubble lcdmopt ", lcdmOpt];
hvalue=hubbleToday[parcosmopts];
physunits=OptionValue[physicalUnits];  (*If True, returns H(z) in units of km/s/Mpc*)
Which[lcdmOpt==True || OptionValue[externalFile]==False,
units=externalUnitsConverter[$internalH0units,hvalue, physicalUnits->physunits, inverseDistanceUnits->True, internalUnits->True];
units*LCDMDimensionlessHubble[z,parcosmopts],
lcdmOpt==False && OptionValue[externalFile]==True,
units=externalUnitsConverter[$externalH0units,hvalue, physicalUnits->physunits, inverseDistanceUnits->True, internalUnits->False];
units*externalHubbleFunction[z, parcosmopts]
]];

Options[DimensionlessHubble]=$paramoptions~Join~$modcosmoopt

DimensionlessHubble[zred_?NumericQ,opts:OptionsPattern[]]:=Block[{paropts},
paropts=complementParamValues[{opts},DimensionlessHubble,returnList->"Full"];
debugPrint["Dimensionlesshubble paropts ", paropts];
Hubble[zred,paropts]/Hubble[0.,paropts]];



Options[ndsolDistance]=$paramoptions~Join~$modcosmoopt;
ndsolDistance[opts:OptionsPattern[]]:=ndsolDistance[opts]=NDSolveValue[{yd'[zx]==(1.0/DimensionlessHubble[zx,opts]),yd[0]==0.},yd,{zx,0.,$zMaxIntegral}];


Options[comovingDistance]=$paramoptions~Join~$modcosmoopt~Join~{externalFile->False, dimensionless->True}


comovingDistance[zred_?NumericQ,opts:OptionsPattern[]]:=Block[{zt=zred, paropts,parcosmopts, lcdmOpt, hvalue, hubfidpar, unitsfact,unitsfactext,unitsfinal, comdis},
  parcosmopts=complementParamValues[{opts},comovingDistance,returnList->"Complement", filterCosmoPars->True];
  paropts=complementParamValues[{opts},comovingDistance,returnList->"Complement"];
  lcdmOpt=OptionValue[lcdmBool];
  hvalue=hubbleToday[parcosmopts];
  hubfidpar=hubbleToday[];
  unitsfact=externalUnitsConverter[$internalDistanceUnits,hvalue, internalUnits->True, inverseDistanceUnits->False, physicalUnits->False];
  Which[
      (OptionValue[externalFile]==False),
        comdis=ndsolDistance[paropts][zt];
        unitsfinal=If[OptionValue[dimensionless]==True, 1, unitsfact];
    ,
      (OptionValue[externalFile]==True),
        If[OptionValue[dimensionless]==True,
          unitsfactext=externalUnitsConverter[$externalDistanceUnits,hubfidpar, internalUnits->False, inverseDistanceUnits->False, physicalUnits->False];
          unitsfact=externalUnitsConverter[$internalDistanceUnits,hubfidpar, internalUnits->True, inverseDistanceUnits->False, physicalUnits->False];
          unitsfinal=1/(unitsfact*unitsfactext);
          If[MemberQ[Flatten[List@@@parcosmopts , 1], hubble]==True,
            parcosmopts=complementParamValues[{parcosmopts},comovingDistance,returnList->"Fiducials", filterCosmoPars->True];
            ];
          ,
          unitsfactext=externalUnitsConverter[$externalDistanceUnits,hvalue, internalUnits->False, inverseDistanceUnits->False, physicalUnits->False];
          unitsfinal=unitsfactext;
          ];
        comdis=(1+zt)*externalDistanceFunction[zt, parcosmopts];

  ];
  Return[unitsfinal*comdis]
];

Options[distanceAngular]=$paramoptions~Join~$modcosmoopt~Join~{externalFile->False}

(*removed memoization of: distanceAngular[zred,opts], this allows to change units, by changing $internalDistanceUnits or $externalDistanceUnits*)

distanceAngular[zred_?NumericQ,opts:OptionsPattern[]]:=Block[{zt=zred, paropts,parcosmopts, lcdmOpt, hvalue, unitsfact, comdist, Omegak},
parcosmopts=complementParamValues[{opts},distanceAngular,returnList->"Complement", filterCosmoPars->True];
paropts=complementParamValues[{opts},distanceAngular,returnList->"Complement"];
lcdmOpt=OptionValue[lcdmBool];
debugPrint["distAng lcdmopt ", lcdmOpt];
debugPrint["distAng paropts ", paropts];
hvalue=hubbleToday[parcosmopts];
Omegak=OmegaK0Today[parcosmopts];
Which[(lcdmOpt==True || OptionValue[externalFile]==False),
unitsfact=externalUnitsConverter[$internalDistanceUnits,hvalue, internalUnits->True,
  inverseDistanceUnits->False, physicalUnits->False];


comdist=If[Omegak==0,
ndsolDistance[paropts][zt]
,
Sin[Sqrt[-Omegak]*ndsolDistance[paropts][zt]]/Sqrt[-Omegak]
];

unitsfact*(1/(1+zt))*comdist

,

lcdmOpt==False && OptionValue[externalFile]==True,
unitsfact=externalUnitsConverter[$externalDistanceUnits,hvalue, internalUnits->False, inverseDistanceUnits->False, physicalUnits->False];
unitsfact*externalDistanceFunction[zt, parcosmopts]
]  (*The volume, the density and the power spectrum, should have compatible units.
    *)
];


Options[externalDistanceFunction]=$paramoptions~Join~$modcosmoopt~Join~{shiftedDomain->False, inputNumericalDerivatives->False};

externalDistanceFunction[zz_, opts:OptionsPattern[]]:=Block[{cosmopts,filename,zdom=0, inpunumder, exthub},

  inpunumder=OptionValue[inputNumericalDerivatives];

  cosmopts=complementParamValues[{opts}, externalDistanceFunction, returnList->"Complement", filterCosmoPars->True];

  If[OptionValue[shiftedDomain]==True, zdom=$externalDistanceInterpolatingFunction[filename]["Domain"][[1,1]] ];


  If[inpunumder==False,
    filename=getCosmoParameterFileName[externalDistanceFunction,{opts}];
    exthub=$externalDistanceInterpolatingFunction[filename][zz+zdom];
    ,
    exthub=(externalFunctionTaylor[zz, cosmopts~Join~{
      externalFunction->$externalDistanceInterpolatingFunction,
      externalDerivativeFunction->$externalDistanceDerivativesInterpolatingFunction,
      interpolatedArguments->"z", kdependence->False}])
  ];
  Return[exthub]
]



Options[ndsolSoundHorizon]=$paramoptions~Join~$modcosmoopt;
ndsolSoundHorizon[opts:OptionsPattern[]]:=ndsolSoundHorizon[opts]=With[{omb=(OmegaBaryon0Today[opts]*hubbleToday[opts]^2)},
  NDSolveValue[{yd'[zx]==((1.0/Sqrt[3.0*(1+30000.0*omb*(1/(1+zx)))])/DimensionlessHubble[zx,opts]),yd[0]==0.},yd,{zx,0.,500.}]
];

Options[distanceSoundHorizon]=$paramoptions~Join~$modcosmoopt~Join~{externalFile->False}

distanceSoundHorizon[zred_?NumericQ,opts:OptionsPattern[]]:=Block[{zt=zred, paropts,parcosmopts, lcdmOpt, hvalue, unitsfact},
  parcosmopts=complementParamValues[{opts},distanceSoundHorizon,returnList->"Complement", filterCosmoPars->True];
  paropts=complementParamValues[{opts},distanceSoundHorizon,returnList->"Complement"];
  lcdmOpt=OptionValue[lcdmBool];
  (*debugPrint["distAng parcosmopts ", parcosmopts];*)
  hvalue=hubbleToday[parcosmopts];
  Which[(lcdmOpt==True || OptionValue[externalFile]==False),
    unitsfact=externalUnitsConverter[$internalDistanceUnits,hvalue, internalUnits->True, inverseDistanceUnits->False, physicalUnits->False];
    unitsfact*ndsolSoundHorizon[paropts][zt]
  ]
];

Options[volumeZbin]=$paramoptions~Join~$modcosmoopt

(*removed memoization*)



volumeZbin[z1_,z2_,opts:OptionsPattern[]]:=Block[{paropts, intunits, return},
  intunits=$internalDistanceUnits;
  $internalDistanceUnits=compatibleUnits[$internalPkUnits];
  paropts=complementParamValues[{opts},volumeZbin,returnList->"Complement"];
  return=4(Pi/3)*((1+z2)^3*distanceAngular[z2,paropts]^3-(1+z1)^3*distanceAngular[z1,paropts]^3)/(4*Pi*(180/Pi)^2);
  $internalDistanceUnits=intunits;
  Return[return];
]
  (*distance Angular has the units set by $internalDistanceUnits or $externalDistanceUnits *)
  (* the volume inside a z-bin, which goes into the survey volume, has to be compatible with the units of the power spectrum for the Fisher matrix *)
  (*volumeZbin is a volume per square steradian*)


Options[OmegaMLCDM]=$paramoptions
OmegaMLCDM[z_?NumericQ,opts:OptionsPattern[]]:=OmegaMLCDM[z,opts]=Block[{OmegaM0, paropts},
paropts=complementParamValues[{opts},OmegaMLCDM,returnList->"Complement", filterCosmoPars->True];
OmegaM0=OmegaM0Today[paropts];
OmegaM0*(1+z)^3/((LCDMDimensionlessHubble[z,paropts])^2 )];


Options[OmegaM]=$paramoptions~Join~$modcosmoopt
OmegaM[z_?NumericQ,opts:OptionsPattern[]]:=OmegaM[z,opts]=Block[{OmegaM0, parcosmopts, paropts},
parcosmopts=complementParamValues[{opts},OmegaM,returnList->"Complement", filterCosmoPars->True];
paropts=complementParamValues[{opts},OmegaM,returnList->"Complement"];
OmegaM0=OmegaM0Today[paropts];
OmegaM0*(1+z)^3/DimensionlessHubble[z,paropts]^2];



scalarAmplitude::noparam="Symbol As or a function of it, is not part of the $paramoptions parameters"


Options[scalarAmplitude]=$paramoptions~Join~{transformedParameter->True};

scalarAmplitude[opts:OptionsPattern[]]:=Block[{paropts, Asvalue, rulepar},
paropts=complementParamValues[{opts},scalarAmplitude,returnList->"Full", filterCosmoPars->True];
rulepar=FilterRules[paropts,_?(StringMatchQ[ToString[#],"*"<>SymbolName[As]<>"*"]&)];
Asvalue=Which[
  SameQ[rulepar, {}],
    Message[scalarAmplitude::noparam];
    Return[$Failed],
  OptionValue[transformedParameter]==True && UnsameQ[First@rulepar[[1]], As],
    setORreadTransformedParameters[As, paropts, convertTransformedParameter->True],
  (OptionValue[transformedParameter]==False && FilterRules[paropts,As]=!={}),
    As/.FilterRules[paropts,As],
  (OptionValue[transformedParameter]==True && SameQ[First@rulepar[[1]], As]),
    SetOptions[scalarAmplitude, transformedParameter->False];
    As/.FilterRules[paropts,As]
  ]
]



Options[extdGNDSolv]=$paramoptions~Join~$modcosmoopt


extdGNDSolv[kvalu_,zerovalu_?NumericQ,endvalu_?NumericQ,opts:OptionsPattern[]]:=extdGNDSolv[kvalu,zerovalu,endvalu,opts]=Block[
  {kv=ReleaseHold[kvalu], znd, yy},

  debugPrint[zerovalu];
  debugPrint[endvalu];

  NDSolveValue[{-externalGrowthRateFunction[znd,kv,opts]/(1+znd)==yy'[znd],yy[zerovalu]==0.},yy,{znd,zerovalu,endvalu}]
];


Options[extGrowthDz]=$paramoptions~Join~$modcosmoopt~Join~{kdependentGrowth->$kdependentGrowth}
extGrowthDz[zk__, fulldomain_,opts:OptionsPattern[]]:=extGrowthDz[zk,fulldomain,opts]=Block[{kdep, zet, paropts, dom, growdz, zeroval, endval, zval, kvoi, kval, dompos, dummya, dummyb},
paropts=complementParamValues[{opts},extGrowthDz,returnList->"Complement", filterCosmoPars->True];
kdep=OptionValue[kdependentGrowth];
{zval, kvoi}= twoArgumentsCheck[zk,kdep];
kval=ReleaseHold[kvoi];
dom=Last@fulldomain;
zeroval=dom[[1]];
endval=dom[[2]];
growdz[kval]=extdGNDSolv[kvoi,zeroval,endval,paropts];
If[zval<zeroval, zval=zeroval]; (*To avoid getting out of the domain of the interpolation function*)
Exp[growdz[kval][zval]]
]

getgamma[opts_]:=Block[{gam, gamrule},
  gamrule=FilterRules[opts,gamma];
  gam=If[gamrule!={} ,
    gamma/.gamrule
    ,
    $fGamma
  ];
  Return[gam]
]

Options[logGdz]=$paramoptions
logGdz[opts_]:=logGdz[opts]=Block[{zn,za,yy, gam, methodsolve=OptionValue[fGrowthRate, solveDifferentialEquation]},
  gam = getgamma[opts];
  If[methodsolve==True,
    NDSolveValue[{yy'[za]==-fGrowthRate[za, opts]/(1+za),yy[0.]==0},yy,{za,0.,100.}],
    NDSolveValue[{-OmegaMLCDM[zn,opts]^gam/(1+zn)==yy'[zn],yy[0.]==0.},yy,{zn,0.,100.}]
  ]
];


(*remove memoization of LCDMGrowth, logGdz is already memoized*)
Options[LCDMGrowth]=$paramoptions
LCDMGrowth[z_?NumericQ,opts:OptionsPattern[]]:=Block[{paropts},
  paropts=complementParamValues[{opts},LCDMGrowth,returnList->"Complement", filterCosmoPars->True];
  Exp[logGdz[paropts][z]]
]


Options[Growth]=$paramoptions~Join~$modcosmoopt~Join~{growthRateIntegral->False, externalFile->False,
  kdependentGrowth:>$kdependentGrowth,
getDomain->False, fixGrowthScalek:>$fkfix};

Growth[zk__,opts:OptionsPattern[]]:=Block[{lcdmOpt,paropts, kdep, zval,
  kval, kuf, puf, kscale, return, domain, dummya, dummyb, pargropts},
paropts=complementParamValues[{opts},Growth,returnList->"Complement", filterCosmoPars->True];
pargropts=complementParamValues[{opts},externalGrowthFunction,returnList->"Complement"];
lcdmOpt=OptionValue[lcdmBool];
kdep=OptionValue[kdependentGrowth];
{zval, kval}= twoArgumentsCheck[zk,kdep];
kval=ReleaseHold[kval];
If[kdep && NumberQ[$fkfix], kval=$fkfix];
Which[lcdmOpt==True && kdep==False || (lcdmOpt==False && kdep==False && OptionValue[externalFile]==False),
return=LCDMGrowth[zval,paropts],
(OptionValue[growthRateIntegral]==True) || (OptionValue[externalFile]==False),
domain=Head[externalGrowthRateFunction[dummya,dummyb,pargropts]]["Domain"];
return=extGrowthDz[zval,kval,domain, paropts],
(OptionValue[externalFile]==True && lcdmOpt==False) || (OptionValue[externalFile]==True && OptionValue[growthRateIntegral]==False),
domain=Head[externalGrowthFunction[dummya,dummyb,pargropts]]["Domain"];
If[OptionValue[getDomain]==True, Return[domain]];
{kuf, puf} = pkUnitsConverter[$internalPkUnits, hubbleToday[paropts], externalFile->True];
kscale=kuf*kval;
return=externalGrowthFunction[zval,kscale,paropts]
]];


Options[externalGrowthFunction]=$paramoptions~Join~$modcosmoopt~Join~{
  kdependentGrowth->$kdependentGrowth,
  shiftedDomain->False, normalizeGrowth->False,
  inputNumericalDerivatives->False}


externalGrowthFunction[zetk__,opts:OptionsPattern[]]:=Block[{filename, len, kdep, zval,
  kval, zdom=0, normfact=1, inpunumder, returnfunc, cosmopts},
inpunumder=OptionValue[inputNumericalDerivatives];
cosmopts=complementParamValues[{opts},Growth,returnList->"Complement", filterCosmoPars->True];
filename=getCosmoParameterFileName[externalGrowthFunction,{opts}];
kdep=OptionValue[kdependentGrowth];
{zval, kval}= twoArgumentsCheck[zetk,kdep];
kval=ReleaseHold[kval];
If[inpunumder==True,
returnfunc=(externalFunctionTaylor[zval,kval, cosmopts~Join~{
                        externalFunction->$externalGrowthInterpolatingFunction,
                        externalDerivativeFunction->$externalGrowthDerivativesInterpolatingFunction,
                        interpolatedArguments->"zk", kdependence->kdep}])
, (*ELSE*)
If[OptionValue[shiftedDomain]==True && kdep==False, zdom=$externalGrowthInterpolatingFunction[filename]["Domain"][[1,1]]];
If[OptionValue[normalizeGrowth]==True,
returnfunc=$externalGrowthInterpolatingFunction[filename][zval+zdom, kval]/$externalGrowthInterpolatingFunction[filename][0.+zdom, kval],
returnfunc=$externalGrowthInterpolatingFunction[filename][zval+zdom, kval]
]
];
Return[returnfunc]
]


Options[externalGrowthRateFunction]=$paramoptions~Join~$modcosmoopt~Join~{
  kdependentGrowth->$kdependentGrowth,
  shiftedDomain->False, inputNumericalDerivatives->False}


externalGrowthRateFunction[zk__,opts:OptionsPattern[]]:=Block[{filename,zval,kval,kdep,dompos,
  kscale, kuf, puf,zdom=0, returnfunc, cosmopts},
filename=getCosmoParameterFileName[externalGrowthRateFunction,{opts}];
inpunumder=OptionValue[inputNumericalDerivatives];
kdep=OptionValue[kdependentGrowth];
cosmopts=complementParamValues[{opts},Growth,returnList->"Complement", filterCosmoPars->True];
{zval, kval}= twoArgumentsCheck[zk,kdep];
kval=ReleaseHold[kval];

If[inpunumder==True,
returnfunc=(externalFunctionTaylor[zval,kval, cosmopts~Join~{
                        externalFunction->$externalGrowthRateInterpolatingFunction,
                        externalDerivativeFunction->$externalGrowthRateDerivativesInterpolatingFunction,
                        interpolatedArguments->"zk", kdependence->kdep}])
, (*ELSE*)
If[OptionValue[shiftedDomain]==True && kdep==False,
zdom=$externalGrowthRateInterpolatingFunction[filename]["Domain"][[1,1]]
];
returnfunc=$externalGrowthRateInterpolatingFunction[filename][zval+zdom, kval]
];

Return[returnfunc]
]

Options[fGrowthSolveND]=$paramoptions
Tini = NfOfz[100.];
Tfin = NfOfz[0.];
(*debugPrint["paropts in fGrowthSolveND", paropts];*)
fGrowthSolveND[optis_]:=fGrowthSolveND[optis]=NDSolveValue[
	{f'[t]+f[t]^2+(2+D[LCDMDimensionlessHubble[zOfNf[t],optis],t] /LCDMDimensionlessHubble[zOfNf[t],optis])f[t]==3/2 OmegaMLCDM[zOfNf[t], optis],
f[Tini]==OmegaMLCDM[zOfNf[Tini],optis]^$fGamma},
f,{t,Tini-2,Tfin}];


Options[fGrowthRate]=$paramoptions~Join~$modcosmoopt~Join~{growthDerivative->False,
  externalFile->False, kdependentGrowth:>$kdependentGrowth,
fixGrowthScalek:>$fkfix, fGammaFit->$fGamma, solveDifferentialEquation -> False}
(*added memoization of fGrowthRate*)
fGrowthRate[zk__,opts:OptionsPattern[]]:=fGrowthRate[zk,opts]=Block[{lcdmOpt,paropts, fromDeriv,kdep, zval,
  kval, kscale, kuf, puf, return, gamma, solveND},
paropts=complementParamValues[{opts},fGrowthRate,returnList->"Complement", filterCosmoPars->True];
lcdmOpt=OptionValue[lcdmBool];
kdep=OptionValue[kdependentGrowth];
If[(FilterRules[$paramoptions,w0]!={} || FilterRules[$paramoptions,wa]!={}),
solveND=True;
SetOptions[fGrowthRate, solveDifferentialEquation -> True],
solveND=OptionValue[solveDifferentialEquation]];
{zval, kval}= twoArgumentsCheck[zk,kdep];
kval=ReleaseHold[kval];
gamma=getgamma[paropts];
fromDeriv=OptionValue[growthDerivative];
debugPrint["z "<>(ToString@zval)];
(*debugPrint["k "<>(ToString@kval)];*)
If[kdep && NumberQ[$fkfix], kval=$fkfix];
Which[
lcdmOpt==True && kdep==False && solveND==False,
(*Print["Omegam^gamma"];*)
OmegaMLCDM[zval,paropts]^gamma
,
solveND == True && OptionValue[externalFile]==False && kdep==False,
(*Print["fGrowthNDSolve"];*)
fGrowthSolveND[paropts][NfOfz[zval]]
,
(OptionValue[externalFile]==False) || (fromDeriv==True),
-(1+zval)/Growth[zval,kval,paropts]*(D[Growth[zx,kval,paropts],zx]/.{zx->zval})(*careful, this has to be tested with k-dep Growth*)
,
(lcdmOpt==False && OptionValue[externalFile]==True) || (fromDeriv==False && OptionValue[externalFile]==True),
{kuf, puf} = pkUnitsConverter[$internalPkUnits, hubbleToday[paropts], externalFile->True];
kscale=kuf*kval;
externalGrowthRateFunction[zval,kscale,paropts]
]];


Options[setBiasFunction]=InterpolationOrder->1
setBiasFunction[redshiftsbins_,biasvals_,opts:OptionsPattern[Interpolation]]:=Block[{intorder,meanredshifts,lista},
meanredshifts=MovingAverage[redshiftsbins,2];
$biasInterpFunc=Interpolation[Transpose[{meanredshifts,biasvals}],opts,Sequence@@Options@setBiasFunction];
  Return[$biasInterpFunc]
];


Options[betaRSDfunction]=$paramoptions~Join~$modcosmoopt~Join~{kdependentGrowth:>$kdependentGrowth}
(*remove memoization*)
betaRSDfunction[zk__,opts:OptionsPattern[]]:=Block[{kval,zval,kdep},
kdep=OptionValue[kdependentGrowth];
{zval, kval}= twoArgumentsCheck[zk,kdep];
kval=ReleaseHold[kval];
fGrowthRate[zval,kval,opts]/$biasInterpFunc[zval]]


Options[sigmaV]=$paramoptions~Join~$modcosmoopt
sigmaV[zref_,opts:OptionsPattern[]]:=Block[{unitsfactor,return,paropts, hvalue},
  paropts=complementParamValues[{opts},sigmaV, returnList->"Full", filterCosmoPars->True];
  hvalue=hubbleToday[paropts];
  unitsfactor=compatibleUnits["1/Mpc",$internalPkUnits, hvalue, inverseDistanceUnits->True];
  return=unitsfactor*($sigmaPecVel*(1+zref))/Hubble[zref, physicalUnits->True];
  Return[return]
(* sigmaV,r = sigmaV[km/s] /(H(z)[km/s/Mpc]*a(z)) units: 1/Mpc *)
(* $sigmaPecVel nonlinear peculiar velocity dispersion in km/s*)
(*remove opts from Hubble, this Hubble function should be just the fiducial one in units of km/s/Mpc*)
]

sigmaZ[zref_]:=Abs[$dzSpectError](1+zref);

Options[sigmaR]=$paramoptions~Join~$modcosmoopt;
(*sigmaR[zref_,opts:OptionsPattern[]]:=Module[{ret, intH0units},
  intH0units=$internalH0units;  (*save the global H0 units to a temporary file, then use the compatible units and reset them to the previous one*)
  $internalH0units = compatibleUnits[$internalPkUnits, invertUnits->False];
  ret=sigmaZ[zref]/(Hubble[zref,opts]*$vd5[zref]);
  $internalH0units = intH0units;
  Return[ret]
];
*)
sigmaR[zref_,opts:OptionsPattern[]]:=Module[{ret, paropts, hvalue, compunits, intH0units, unitsfact},
  paropts=complementParamValues[{opts},sigmaR, returnList->"Full", filterCosmoPars->True];
  hvalue=hubbleToday[paropts];
  ret=sigmaZ[zref]/(Hubble[zref,opts]*$vd5[zref]);
  compunits=compatibleUnits[$internalPkUnits, invertUnits->False];
  unitsfact=compatibleUnits[$internalH0units, compunits, hvalue, inverseDistanceUnits->True];
  Return[(1/unitsfact)*ret]
];


Options[errorZ]=$paramoptions~Join~$modcosmoopt
errorZ[k_,mu_,zref_,opts:OptionsPattern[]]:=Exp[-k^2*mu^2*(sigmaR[zref,opts]^2+sigmaV[zref,opts]^2)];
(*errorZ should not change when cosmological parameters are changed, if it is evaluated in the observed power spectrum*)

Options[dampingTerm]=$paramoptions~Join~$modcosmoopt~Join~{growthDerivative->False, externalFile->False, kdependentGrowth:>$kdependentGrowth}

dampingTerm[zref_,k_,mu_,opts:OptionsPattern[]]:=Block[{damp=$dampingSwitch},
If[damp==0, 1,
Exp[-(k*$dampingGrowthNorm*Growth[zref,k,opts])^2*($dampingSigma0Value^2+mu^2*(($dampingSigma0Value*(1+fGrowthRate[zref,k,opts]))^2-$dampingSigma0Value^2))]
]];


windowk[x_] := 3*x^-3*(Sin[x] - x*Cos[x]);

Options[sigma8Function]={MaxRecursion->15,PrecisionGoal->6,AccuracyGoal->6}
sigma8Function[ps_, scal_?NumericQ,opts:OptionsPattern[NIntegrate]]:=(1/(2Pi^2))*NIntegrate[Exp[3 lk]ps[Exp@lk]*windowk[Exp[lk]*scal]^2,
{lk, Log[$kminRef],Log[$kmaxRef]},opts];

normalizationFactor[ps_,sig8_]:=(sig8^2)/sigma8Function[ps, $radiusScale];

rescaleSigma8fromAs::noparam="Parameter options do not contain As or any function of it."

rescaleSigma8fromAs[s8ref_?NumberQ, asvalue_?NumberQ]:=rescaleSigma8fromAs[s8ref,asvalue]=Block[{s8val, asref, fiduopts},
fiduopts=$paramoptions;
If[SameQ[FilterRules[{fiduopts},_?(StringMatchQ[ToString[#],"*"<>SymbolName[As]<>"*"]&)], {}],
Message[rescaleSigma8fromAs::noparam];
Return[$Failed]
,
asref=scalarAmplitude[fiduopts];
s8val = s8ref*Sqrt[asvalue/asref]
]
];



rescaleAsfromSigma8[s8val_?NumberQ, s8ref_?NumberQ, asref_?NumberQ]:=Block[{asval, fiduopts},
    asval = (s8val^2/s8ref^2)*asref
];


Options[externalSigma8ofZFunction]=$paramoptions~Join~$modcosmoopt~Join~{shiftedDomain->False, inputNumericalDerivatives->False};

externalSigma8ofZFunction[zred_,opts:OptionsPattern[]]:=Block[{ zdom=0, cosmopts, inpunumder, exts8z, filename},
  inpunumder=OptionValue[inputNumericalDerivatives];
  cosmopts=complementParamValues[{opts}, externalSigma8ofZFunction, returnList->"Complement", filterCosmoPars->True];
  If[OptionValue[shiftedDomain]==True, zdom=$externalSigma8ofZInterpolatingFunction[filename]["Domain"][[1,1]]];

  If[inpunumder==False,
    filename=getCosmoParameterFileName[externalSigma8ofZFunction,{opts}];
    exts8z=$externalSigma8ofZInterpolatingFunction[filename][zred+zdom];
    ,
    exts8z=(externalFunctionTaylor[zred, cosmopts~Join~{
      externalFunction->$externalSigma8ofZInterpolatingFunction,
      externalDerivativeFunction->$externalSigma8ofZDerivativesInterpolatingFunction,
      interpolatedArguments->"z", kdependence->False}])
  ];
  Return[exts8z]
]

Options[sigma8ofZ]={MaxRecursion->15,PrecisionGoal->6,AccuracyGoal->6, topHatWindowR->$radiusScale,
spectrumMethod->"Transfer", sigma8reference->$sigma8reference, externalFile->False}~Join~$paramoptions~Join~$pscosmoopts

sigma8ofZ[zet_?NumericQ,opts:OptionsPattern[{NIntegrate, sigma8ofZ}]]:=sigma8ofZ[zet,opts]=Block[{kk=Global`k,
paropts=complementParamValues[{opts},sigma8ofZ, returnList->"Complement", filterCosmoPars->True],
  powopts=FilterRules[complementParamValues[{opts},sigma8ofZ, returnList->"Complement"], Options@powerSpectrum],
scale=OptionValue[topHatWindowR], lcdmOpt=OptionValue[lcdmBool], intopts, s8res, unitconv},
  intopts=FilterRules[Options[sigma8ofZ], Options[NIntegrate]];
  unitconv = (pkUnitsConverter[$internalPkUnits, (hubble/.$paramoptions )])[[1]];
  scale = scale*unitconv;
  Which[
    (lcdmOpt==True && OptionValue[externalFile]==False) &&  $sigma8reference>=0.1, (*only do this if sigma8 has a reasonable value*)
     s8res=$sigma8reference*Growth[zet]
    ,
    (lcdmOpt==False && OptionValue[externalFile]==False),
       s8res = Sqrt[(1/(2Pi^2))*NIntegrate[Exp[3 kk]*powerSpectrum[zet,Exp[kk], powopts]*windowk[Exp[kk]*scale]^2, {kk, Log@$kminRef,Log@$kmaxRef}(**)]]
    ,
    OptionValue[externalFile]==True,
    s8res = externalSigma8ofZFunction[zet,paropts]
  ];
  Return[s8res]
];


Options[fGrowthRateSigma8ofZ]={sigma8reference->$sigma8reference, externalFile->False, kdependentGrowth:>$kdependentGrowth,
  fixGrowthScalek:>$fkfix}~Join~$paramoptions~Join~$pscosmoopts;

fGrowthRateSigma8ofZ[zk__?NumericQ,opts:OptionsPattern[]]:=fGrowthRateSigma8ofZ[zk,opts]=Block[{kk=Global`k,
  paropts=complementParamValues[{opts},fGrowthRateSigma8ofZ, returnList->"Complement", filterCosmoPars->True],
  lcdmOpt=OptionValue[lcdmBool], intopts, s8res, kdep=OptionValue[kdependentGrowth], fres, fs8res, zval, kval},
  {zval, kval}= twoArgumentsCheck[zk,kdep];
  kval=ReleaseHold[kval];
  If[kdep && NumberQ[$fkfix], kval=$fkfix];
  Which[
    lcdmOpt==True || OptionValue[externalFile]==False,
    s8res = sigma8ofZ[zval,opts];
    fres=fGrowthRate[zval,kval, paropts];
    fs8res=fres*s8res;
    ,
    lcdmOpt==False || OptionValue[externalFile]==True,
    fs8res = externalfGrowthRateSigma8ofZFunction[zval,paropts]
  ];
  Return[fs8res]
];


Options[externalfGrowthRateSigma8ofZFunction]=$paramoptions~Join~$modcosmoopt~Join~{shiftedDomain->False, inputNumericalDerivatives->False};

externalfGrowthRateSigma8ofZFunction[zred_,opts:OptionsPattern[]]:=Block[{ zdom=0, cosmopts, inpunumder, extfs8z, filename},
  inpunumder=OptionValue[inputNumericalDerivatives];
  cosmopts=complementParamValues[{opts}, externalfGrowthRateSigma8ofZFunction, returnList->"Complement", filterCosmoPars->True];
  If[OptionValue[shiftedDomain]==True, zdom=$externalfGrowthRateSigma8ofZInterpolatingFunction[filename]["Domain"][[1,1]] ];

  If[inpunumder==False,
    filename=getCosmoParameterFileName[externalfGrowthRateSigma8ofZFunction,{opts}];
    extfs8z=$externalSigma8ofZInterpolatingFunction[filename][zred+zdom];
    ,
    extfs8z=(externalFunctionTaylor[zred, cosmopts~Join~{
      externalFunction->$externalfGrowthRateSigma8ofZInterpolatingFunction,
      externalDerivativeFunction->$externalfGrowthRateSigma8ofZDerivativesInterpolatingFunction,
      interpolatedArguments->"z", kdependence->False}])
  ];
  Return[extfs8z]
]


Options[SigmaLensing]={sigma8reference->$sigma8reference, externalFile->False,
       theoreticalMGFunction->False, theoreticalKdependence->False}~Join~$paramoptions~Join~$pscosmoopts;

(*Remove memoization:   SigmaLensing[zk__,opts:OptionsPattern[]]:=SigmaLensing[zk,opts]=Block[{kk=Global`k,*)
SigmaLensing[zk__,opts:OptionsPattern[]]:=Block[{kk=Global`k,
  paropts=complementParamValues[{opts},SigmaLensing, returnList->"Complement", filterCosmoPars->True],
  lcdmOpt=OptionValue[lcdmBool], intopts, Sres, zval, kval,
  MGfunc=OptionValue[theoreticalMGFunction],
  kdep=OptionValue[theoreticalKdependence]},
  {zval, kval}= twoArgumentsCheck[zk,True];
  kval=ReleaseHold[kval];
  Which[
    lcdmOpt==True && OptionValue[externalFile]==False,
    Sres=1;
    ,
    lcdmOpt==False && OptionValue[externalFile]==True,
    Sres = externalSigmaLensingFunction[zval,kval,paropts]
    ,
    lcdmOpt==False && UnsameQ[MGfunc,False],
    Sres = If[kdep,
              MGfunc[zval,kval,paropts],
              MGfunc[zval,paropts]];
    ,
    True (*If all above fails*),
    Sres=1
    ];
  Return[Sres]
];



Options[externalSigmaLensingFunction]=$paramoptions~Join~$modcosmoopt~Join~{shiftedDomain->False, inputNumericalDerivatives->False};

externalSigmaLensingFunction[zred_,kk_,opts:OptionsPattern[]]:=Block[{ zdom=0, cosmopts, inpunumder, extSigma, filename},
  inpunumder=OptionValue[inputNumericalDerivatives];
  cosmopts=complementParamValues[{opts}, externalSigmaLensingFunction, returnList->"Complement", filterCosmoPars->True];
  If[OptionValue[shiftedDomain]==True, zdom=$externalSigmaLensingInterpolatingFunction[filename]["Domain"][[1,1]] ];

  If[inpunumder==False,
    filename=getCosmoParameterFileName[externalSigmaLensingFunction,{opts}];
    extSigma=$externalSigmaLensingInterpolatingFunction[filename][zred+zdom,kk];
    ,
    extSigma=(externalFunctionTaylor[zred,kk, cosmopts~Join~{
      externalFunction->$externalSigmaLensingInterpolatingFunction,
      externalDerivativeFunction->$externalSigmaLensingDerivativesInterpolatingFunction,
      interpolatedArguments->"k"}])
  ];
  Return[extSigma]
]


Options[LCDMCAMBPsPre]=$paramoptions~Join~$pslinearopt~Join~{returnInterpolated->False,
  checkConsistency->True, returnRescaledAs->False, returnSigma8->False, returnGrowthRate->False,
  returnDeprecated->False}~Join~{
  HalofitVersion->4,TransferKmax->5, AccuracyBoost->1,TransferHighPrecision->False,TransferKperLogInt->50,
  OmegaNu->0.,MasslessNeutrinos->3.046, MassiveNeutrinos->0.,NuMassFractions->{1},
  activateModifiedGravity->False,
optionsMG->{}};

LCDMCAMBPsPre::numericarg="Input z argument is not a number or list of numbers."
LCDMCAMBPsPre::neutrino="Input Omeganu does not correspond to $Omeganu from setCosmologyFixedValues."


LCDMCAMBPsPre[zred_,opts:OptionsPattern[]]:=LCDMCAMBPsPre[zred,opts]=Block[
{consistencycheck,tole=10^-5,retInterp,cambdata,w0ref,linearString,
  halofitver,accboost,trhighpr,tklogint, numericopts,
  psLog,pslKvalues,fR0,mgopts={},e11par,e22par,
  pslPvalues,psTemp,psFull,norm,OM,OB,OC,OK,ODE,nsval,h,
  waref,posi,linearB,pardefopts,exopts,paropts, pargrowopts, Aspar, Asopt, Asbool, s8ref,
  s8val, cambs8, psTableLin, psTableNonLin, psfullTab, returnQuant, zets, zoutind,zoutintp, zii,zind,zlist,
  tkmax, retas, rets8, s8z,
  Onu, mlessnu, msivenu, numassfrac, numassdeg, neutrinoopts,
  retfDG,redzs, cambs8z,cambfz,cambDz,pertstable, returnAll=False, allQuantsTab
},
retInterp=OptionValue[returnInterpolated];
retas=OptionValue[returnRescaledAs];
rets8=OptionValue[returnSigma8];
retfDG = OptionValue[returnGrowthRate];
retDepr = OptionValue[returnDeprecated];

If[retInterp==True && retfDG==True,
  returnAll=True
  allQuantsTab={}
];

consistencycheck=OptionValue[checkConsistency];


paropts=complementParamValues[{opts}, LCDMCAMBPsPre, returnList->"Full", filterCosmoPars->True];
pargrowopts=complementParamValues[{opts},Growth, returnList->"Complement"];
h=hubbleToday[paropts];
OM=OmegaM0Today[paropts];
OB=OmegaBaryon0Today[paropts];
OC=OmegaCDM0Today[paropts];
Onu=OptionValue[OmegaNu];
If[$Omeganu != Onu, Message[LCDMCAMBPsPre::neutrino]; Abort[] ];
OK=OmegaK0Today[paropts];

Aspar=Quiet[scalarAmplitude[paropts]];
Asbool=(UnsameQ[Aspar,$Failed]&&NumericQ[Aspar]);
If[Asbool,
  Asopt=Aspar;
  s8val=s8ref=$sigma8reference;
  ,
  s8val=(sigma8/.paropts);
  s8ref=(sigma8/.$paramoptions);
  Asopt=rescaleAsfromSigma8[s8val, s8ref, $Asreference];
];
nsval=OptionValue[ns];
linearB=OptionValue[linearBool];
{w0ref,waref}=wEoSparameters[paropts];
halofitver=OptionValue[HalofitVersion];
Which[
  SameQ[OptionValue[activateModifiedGravity],False],
    mgopts={},
  OptionValue[activateModifiedGravity]=="fofR-HS",
    SetOptions[parametersModifiedGravity, inverseTransformation->((10^#)&)];
    fR0 = parametersModifiedGravity[logfR0, paropts];
    mgopts={MGflag->3, QSAflag->4, FR0 -> fR0, FRn->1, DEmodel->0};
    SetOptions[LCDMCAMBPsPre, optionsMG->mgopts];
  ,
  OptionValue[activateModifiedGravity]=="mueta-param",
    e11par = parametersModifiedGravity[E11, paropts];
    e22par = parametersModifiedGravity[E22, paropts];
    mgopts={MGflag->1, pureMGflag->1, mugammapar->2, DEmodel->0, E11->e11par, E22->e22par};
    SetOptions[LCDMCAMBPsPre, optionsMG->mgopts];
  ,
  SameQ[OptionValue[activateModifiedGravity], True],
    mgopts=OptionValue[optionsMG];
];

tkmax=OptionValue[TransferKmax];
accboost=OptionValue[AccuracyBoost];
trhighpr=OptionValue[TransferHighPrecision];
tklogint=OptionValue[TransferKperLogInt];
numericopts={AccuracyBoost->accboost,TransferHighPrecision->trhighpr,
  TransferKmax->tkmax, TransferKperLogInt->tklogint};

mlessnu=OptionValue[MasslessNeutrinos];    (*float, in standard case: massless + massive = 3.046*)
msivenu=OptionValue[MassiveNeutrinos];    (*integer*)
numassfrac=OptionValue[NuMassFractions];  (*length of this list, specifies number of neutrino eigenstates, list consists of fractions adding to 1*)
(*numassdeg=OptionValue[NuMassDegeneracies]; (*list of individual degeneracies for each neutrino *) not available at the moment in CAMB wrapper*)
neutrinoopts={OmegaNu->Onu,MasslessNeutrinos->mlessnu, MassiveNeutrinos->msivenu,NuMassFractions->numassfrac};

If[linearB==True,linearString="PSlinear",linearString="PSnonlinear"];

zets=Which[(Length@zred)==0 && Chop[zred, tole]==0,
  zoutind=zoutintp=-1;
  {0}
  ,
  (Length@zred)==0 && Chop[zred, tole]!=0,
  zoutind=zoutintp=-2;
  {0, zred}
  ,
  (Length@zred)!=0 && Chop[zred[[1]], tole ]!=0,
  zoutind=Range[(Length@zred)];
  zoutintp=-2;
  Prepend[zred,0]
  ,
  (Length@zred)!=0 && Chop[zred[[1]], tole ]==0,
  zoutind=Range[(Length@zred)];
  zoutintp=-1;
  zred
  ,
  (NumericQ@zred)==False || VectorQ[zred, NumericQ]==False,
  Message[LCDMCAMBPsPre::numericarg];
  Abort[]
  ,
  True,
  zoutind=-1;
  {0}
];


debugPrint["Numeric CAMB options : "<>ToString[scientificfortranform[numericopts]], 1];
debugPrint["Constant Neutrino CAMB options : "<>ToString[scientificfortranform[neutrinoopts]], 1];
debugPrint["Modified Gravity CAMB options : "<>ToString[scientificfortranform[mgopts]], 1];

debugPrint["Parameters : "<>ToString[paropts], 1];

cambdata=CAMB[OC,OB,h,
  TransferRedshifts->zets,w0ppf->w0ref, wappf->waref,
  ScalarSpectralIndex->{nsval}, ScalarPowerAmplitude->{Asopt},
  HalofitVersion->halofitver, NonLinear->"pk", numericopts~Join~neutrinoopts~Join~mgopts];

cambs8=((CAMB["sigma8(z)"] /. cambdata)[[-1]]); (*sigma8 at z=0*)


debugPrint["zlist: "<>ToString@(CAMB["redshifts"] /. cambdata), 1];

debugPrint["Values first run:  Omegac: "<>(ToString@OC)<>" , Omegab: "<>(ToString@OB)<>
    " , h: "<>(ToString@h)<>
    " , ns: "<>(ToString@nsval)<>" , As: "<>(ToString@FortranForm[Asopt])
    <>" , sigma8: "<>(ToString@cambs8), 1];


If[rets8==True && retas==False && retInterp==False && consistencycheck==False,
  Return[cambs8]  (*sigma8 at z=0*)
];


If[consistencycheck==True (*&& Asbool==False*),
If[s8val==s8ref, (*if they are equal, sigma has to remain constant when changing other parameters*)
  If[Chop[cambs8-s8val,tole]!=0,
    Asopt=rescaleAsfromSigma8[s8val,cambs8, $Asreference];
    debugPrint["second CAMB run...", 1];
    debugPrint["Rescaled As: "<>(ToString@FortranForm[Asopt]), 1];
    cambdata=CAMB[OC,OB,h,
      TransferRedshifts->zets,w0ppf->w0ref, wappf->waref,ScalarSpectralIndex->{nsval},
      ScalarPowerAmplitude->{Asopt},
      HalofitVersion->halofitver, NonLinear->"pk", numericopts~Join~neutrinoopts~Join~mgopts];
    cambs8=((CAMB["sigma8(z)"] /. cambdata)[[-1]]);
    debugPrint["CAMB sigma8: "<>(ToString@FortranForm[cambs8]), 1];
    debugPrint["Values second run:  Omegac: "<>(ToString@OC)<>" , Omegab: "<>(ToString@OB)<>
        " , h: "<>(ToString@h)<>
        " , ns: "<>(ToString@nsval)<>" , As: "<>(ToString@FortranForm[Asopt])
        <>" , sigma8: "<>(ToString@cambs8), 1];
    If[paropts==$paramoptions, $Asreference = Asopt];
    If[Chop[cambs8-s8val,tole]!=0, Print["Rescaling of As failed, check settings."]; Abort[] ];
    ];
  ];
];



debugPrint["CAMB: As-s8 : "<>ToString[FortranForm[Asopt]]<>" , "<>ToString[cambs8] , 1];


If[rets8==True,
  s8z = ((CAMB["sigma8(z)"] /. cambdata)[[zoutind]]);  (*sigma8 at the specified redshift*)
];

If[retas==True, Return[Asopt]];

If[retfDG==True,
  pertstable = observablesCAMB[cambdata, returnObservable->"Growth", returnInterpolated->False];
   debugPrint["dimensions pert", Dimensions@pertstable, 0];
  Return[pertstable]
];


Which[
  retInterp==True && Length@zoutind == 0,
    pslKvalues=Log10[(CAMB[linearString]/.cambdata)[[zoutintp]][[All,1]]];
    pslPvalues=Log10[(CAMB[linearString]/.cambdata)[[zoutintp]][[All,2]]];
    psLog=Interpolation[Transpose[{pslKvalues,pslPvalues}],InterpolationOrder->3,Method->"Spline"];
    If[rets8==False, returnQuant=psLog, returnQuant={psLog,pslKvalues,s8z, cambs8}];
    If[returnAll==True, AppendTo[allQuantsTab, returnQuant]];
    ,
  Length@zoutind == 0,
    pslKvalues=(CAMB["PSlinear"]/.cambdata)[[zoutind]][[All,1]];
    psTableLin=(CAMB["PSlinear"]/.cambdata)[[zoutind]][[All,2]];
    psTableNonLin=(CAMB["PSnonlinear"]/.cambdata)[[zoutind]][[All,2]];
    psfullTab=Transpose[{pslKvalues,psTableLin,psTableNonLin}];
    returnQuant=psfullTab;
  ,
  Length@zoutind !=0  && retDepr==True,
  psfullTab={};
  debugPrint["Indices of returned redshifts: "<>ToString[zoutind], 1];
  Do[
    pslKvalues=(CAMB["PSlinear"]/.cambdata)[[zii]][[All,1]];
    psTableLin=(CAMB["PSlinear"]/.cambdata)[[zii]][[All,2]];
    psTableNonLin=(CAMB["PSnonlinear"]/.cambdata)[[zii]][[All,2]];
    zind=(CAMB["redshifts"]/.cambdata)[[zii]];
    zlist=ConstantArray[zind, Length@pslKvalues];
    AppendTo[psfullTab,Transpose[{zlist,pslKvalues,psTableLin,psTableNonLin}]];
  ,{zii,zoutind}];
  returnQuant=psfullTab;
  ,
  (*Else, with Default options and an array of z-values*)
  Length@zoutind !=0,
  returnQuant=observablesCAMB[cambdata, returnObservable->"All"];
  ];


Return[returnQuant]
];

(**Needs to be fixed according to new observablesCAMB**)
(*
Options[LcdmCambPk]=$paramoptions~Join~$pslinearopt~Join~{HalofitVersion->4,
checkConsistency->True};

LcdmCambPk[zred_,kk_,opts:OptionsPattern[{LCDMCAMBPsPre, LcdmCambPk}]]:=
10^(LCDMCAMBPsPre[zred,{opts}~Join~{returnInterpolated->True}][Log10[kk]]);
*)



Options[observablesCAMB]={returnObservable->"All", returnInterpolated->False};

observablesCAMB[cambdata_, opts:OptionsPattern[]]:=Block[{return, retInt, obs, transfer,
  redz, revZrang, kvals, krang,
  growth, growthrate, s8growth, s8growthZtab, s8growthrate, s8growthrateZtab,
  s8ofz, s8ofZtab,
  HofZ, HofZtab,
  zdepPertsTab, kzdepTab, zzi, growthTab, growthrateTab, gInterp, fgFunc,
  fgInterp, outZtab, outKtab},
  obs = OptionValue[returnObservable];
  retInt = OptionValue[returnInterpolated];
  redz = CAMB["redshifts"] /. cambdata;
  (*redshifts from CAMB are ordered from larger to smaller, with 0. at the end*)
  transfer = (CAMB["Transfer"]/.cambdata);
  kvals = transfer[[1, All,1]];

  krang = Range[Length[kvals]];
  revZrang = Reverse[Range[Length[redz]]];

  outKtab = Table[kvals[[kki]] , {kki, krang}];
  outZtab = Table[redz[[zzi]]  ,  {zzi, revZrang}];

  psKvalues=(CAMB["PSlinear"]/.cambdata)[[1]][[All,1]];
  If[Total[Chop[outKtab-psKvalues]]!=0,
   Print["** Warning: The k-values from Transfer and from Pk are not equal. Check for consistency"]
  ];

  s8growth[zind_]:=(CAMB["D+(z)"]/.cambdata)[[zind]];
  s8growthrate[zind_]:=(CAMB["f(z)"]/.cambdata)[[zind]];
  s8ofz[zind_]:=(CAMB["sigma8(z)"]/.cambdata)[[zind]];
  zdepPertsTab = Table[{s8ofz[zzi], s8growth[zzi], s8growthrate[zzi]}, {zzi, revZrang}];
  s8ofZtab = zdepPertsTab[[All,1]];
  s8growthZtab = zdepPertsTab[[All,2]];
  s8growthrateZtab = zdepPertsTab[[All,3]];

  Hofz[zind_]:=(CAMB["H(z)"]/.cambdata)[[zind]];
  HofZtab = Table[Hofz[zzi], {zzi, revZrang}];

  debugPrint["obs, retInt: ", obs, retInt, 0];
  growth[zind_, kind_]:=transfer[[zind, kind, 7]]/(transfer[[-1,kind,7]]) ; (*D(z), normalized to 1 at z=0 *)
  growthrate[zind_, kind_]:=transfer[[zind, kind, 14]]/(transfer[[1,kind,14]]) ; (*D(z), normalized to 1 at z=0 *)
  growthTab =Flatten[Table[{{redz[[zzi]], kvals[[kki]]}, growth[zzi, kki]}, {zzi, revZrang}, {kki, krang }], 1];
  gInterp = Interpolation[growthTab, InterpolationOrder->3, Method->"Spline"];

  growthrateTab =Flatten[Table[{{redz[[zzi]], kvals[[kki]]}, growthrate[zzi, kki ]}, {zzi, revZrang}, {kki, krang }], 1];
  fgInterp = Interpolation[growthrateTab, InterpolationOrder->3, Method->"Spline"];

  growthTab = ArrayReshape[growthTab[[All,2]], {Length@revZrang, Length@krang} ];
  growthrateTab = ArrayReshape[growthrateTab[[All,2]], {Length@revZrang, Length@krang} ];
  debugPrint["dimensions: ",
             Dimensions@outZtab, Dimensions@outKtab,
             Dimensions@growthTab, Dimensions@growthrateTab, 1];


  psTableLin=(CAMB["PSlinear"]/.cambdata);
  psTableNonLin=(CAMB["PSnonlinear"]/.cambdata);
  psLin[zind_, kind_]:=psTableLin[[zind,kind]][[2]]; (*first column is k, second is Pk*)
  psNonLin[zind_, kind_]:=psTableNonLin[[zind,kind]][[2]]; (*first column is k, second is Pk*)
  psLinTab =Flatten[Table[{{redz[[zzi]], kvals[[kki]]}, psLin[zzi, kki]}, {zzi, revZrang}, {kki, krang }], 1];
  psNonLinTab =Flatten[Table[{{redz[[zzi]], kvals[[kki]]}, psNonLin[zzi, kki]}, {zzi, revZrang}, {kki, krang }], 1];

  psLinInterp = Interpolation[psLinTab, InterpolationOrder->3, Method->"Spline"];
  psNonLinInterp = Interpolation[psNonLinTab, InterpolationOrder->3, Method->"Spline"];

  psLinTab = ArrayReshape[psLinTab[[All,2]], {Length@revZrang, Length@krang} ];
  psNonLinTab = ArrayReshape[psNonLinTab[[All,2]], {Length@revZrang, Length@krang} ];

  Which[
    obs=="s8fgDz",
    return=Transpose[{outZtab, s8ofZtab, s8growthZtab, s8growthrateZtab}];
    ,
    obs=="Hubble",
    return=Transpose[{outZtab, HofZtab}];
    ,
    obs=="Growth" || obs=="GrowthRate",
    return = {outZtab, outKtab, growthTab, growthrateTab};
    ,
    obs=="powerSpectrum",
    return = {outZtab, outKtab, psLinTab, psNonLinTab};
    ,
    obs=="All",
    return={outZtab, outKtab, HofZtab, s8ofZtab, growthTab, growthrateTab, psLinTab, psNonLinTab};
  ];
Return[return]
];






Options[LCDMTrHuPre]=$paramoptions~Join~$pscosmoopts~Join~{growthRateIntegral->False,
  externalFile->False, kdependentGrowth:>$kdependentGrowth};

LCDMTrHuPre[sigma8val_,opts:OptionsPattern[]]:=LCDMTrHuPre[sigma8val,opts]=Block[
{tfdata,ps,fullps,psInt,psnl,norm,zrange,OM,OB,h,nspar,paropts,pargrowopts, linearB, linps, klist, klistpad, trlist, trlistpad,
duma, dumb, domaink, domainz, domain, Aspar, sigma8par},
paropts=complementParamValues[{opts},LCDMTrHuPre, returnList->"Complement", filterCosmoPars->True];
pargrowopts=complementParamValues[{opts},Growth, returnList->"Complement"];
h=hubbleToday[paropts];
nspar=OptionValue[ns];
OM=OmegaM0Today[paropts];
OB=OmegaBaryon0Today[paropts];
Aspar=Quiet[scalarAmplitude[paropts], {scalarAmplitude::noparam}];
sigma8par=sigma8val;
zrange=Range[0., $zMaxIntp,.2];
tfdata=Transfer[OM,OB/OM , 2.726, h];
debugPrint[{OM,OB/OM , 2.726, h}];
klist=Transfer["kvalues"]/.tfdata;
trlist=Transfer["full"]/.tfdata;
If[OptionValue@lcdmBool==False,
domain=Growth[duma,dumb,pargrowopts~Join~{getDomain->True}];
domainz=Last@domain; (*z-Domain must be always in second place*)
zrange=Range[domainz[[1]], domainz[[2]], 0.2];
If[OptionValue[kdependentGrowth]==True,
domaink=domain[[1]];
duma=LengthWhile[klist,#<domaink[[1]]&];
dumb=LengthWhile[klist,#<domaink[[2]]&];
klist=klist[[duma;;dumb]];
trlist=trlist[[duma;;dumb]];
debugPrint[klist];
debugPrint[Length@klist,Length@trlist];
]];
ps=Transpose[{klist, klist^nspar * trlist^2}];
debugPrint[ps];
psInt=Interpolation@ps;

If[UnsameQ[Aspar,$Failed]&&NumericQ[Aspar],
sigma8par=rescaleSigma8fromAs[sigma8val,Aspar];
];

norm=normalizationFactor[psInt,sigma8par];
(*If Growth is only z-dependent, this still works,
since Growth can accept z and/or k as parameters and it all depends of the option $kdependentGrowth*)

debugPrint["norm: ", norm];

linps=Table[{klist[[kk]],(klist[[kk]]^nspar*trlist[[kk]]^2*norm*(Growth[zz,klist[[kk]],pargrowopts])^2)},{zz,zrange},{kk,1,Length@klist}];

fullps=Flatten[Table[{{zrange[[zi]],klist[[kk]]},(linps[[zi,kk,2]])},{zi,1, Length@zrange},{kk,1,Length@klist}],1];

linearB=OptionValue[linearBool];
debugPrint[fullps[[1]], fullps[[-1]]];

If[linearB,
Interpolation[fullps, InterpolationOrder->2,Method->"Spline"],
psnl=Table[{zrange[[zi]],HaloFitCorrection[linps[[zi]], OM,1-OM]},{zi,1, Length@zrange}];
psnl=Table[{{row[[1]],#[[1]]},#[[2]]}&/@row[[2]],{row,psnl}];
Interpolation@Flatten[psnl,1]
]
];

LCDMCosmicEmuPre::invalpars="Parameter `1` specified with value `2` is not allowed by CosmicEmulator"

Options[LCDMCosmicEmuPre]=$paramoptions~Join~$pslinearopt

LCDMCosmicEmuPre[zred_,sigma8ref_,opts:OptionsPattern[]]:=
LCDMCosmicEmuPre[zred,sigma8ref,opts]=Block[
{emudata,atime,zlist,w0ref,psTemp,psKfunc,norm,oM,oB,nsval,h,waref,posi,pardefopts,exopts,paropts,s8val,sigma8par,Aspar},
paropts=complementParamValues[{opts},LCDMCosmicEmuPre,returnList->"Full", filterCosmoPars->True];
h=hubbleToday[paropts];
oM=OmegaM0Today[paropts]*h^2;
oB=OmegaBaryon0Today[paropts]*h^2;
s8val=sigma8ref;
Aspar=Quiet[scalarAmplitude[paropts]];
If[UnsameQ[Aspar,$Failed]&&NumericQ[Aspar],
s8val=rescaleSigma8fromAs[s8val,Aspar];
];
nsval=OptionValue[ns];

{w0ref,waref}=wEoSparameters[paropts];
If[NumericQ@waref && waref!=0, Message[LCDMCosmicEmuPre::invalpars, ToString[waref], waref]];

atime={aOfz[zred]};

emudata=FrankenEmu[oM,oB,h,s8val,nsval,w0ref,atime];

psTemp=First@((rescaleApplyFunctionTable[h^3*# &, 1&, (1.0/h)*# &] /@ (#)) & /@ (FrankenEmu["pk"] /. emudata)); (*Rescales emuPk to units of h/Mpc*)

psKfunc=Interpolation[psTemp, InterpolationOrder->3, Method->"Spline"]


];

powerSpectrum::extopts="Inconsistent options chosen, lcdmBool=False, but externalFile has not been defined. Switching to lcdm=True."

Options[powerSpectrum]=$paramoptions~Join~$pscosmoopts~Join~{spectrumMethod->"Transfer", sigma8reference->$sigma8reference,
  externalFile->False, kdependentGrowth:>$kdependentGrowth, units:>$internalPkUnits}

powerSpectrum[z_?NumericQ,k_,opts:OptionsPattern[]]:=Block[{paropts, psmethod, extfile, s8, hvalue, unitsPk,
  lcdmOpt, parcosmopts, kuf, puf, pk},

paropts=complementParamValues[{opts},powerSpectrum, returnList->"Complement"];
parcosmopts=complementParamValues[{opts},powerSpectrum, returnList->"Full", filterCosmoPars->True];
psmethod=OptionValue[spectrumMethod];
extfile=OptionValue[externalFile];
lcdmOpt=OptionValue[lcdmBool];  (*lcdmBool option, controls if lcdm internal power spectra are being used*)

If[
  (extfile==False && lcdmOpt==False),
  lcdmOpt=True;
  psmethod="Transfer";
  If[paropts=={}, paropts={(lcdmBool->lcdmOpt)}]; (*forces paropts to have a value, when changing to lcdmBool->True*)
  Message[powerSpectrum::extopts]
];
paropts = paropts /. {(externalFile->x_) -> (externalFile -> extfile),
  (spectrumMethod -> x_)->(spectrumMethod->psmethod),
  (lcdmBool->x_)->(lcdmBool->lcdmOpt)};  (*replaces paropts, with the modified values. *)
s8 = FilterRules[parcosmopts, sigma8];
s8 = If[s8=={},OptionValue[sigma8reference],s8[[1,2]]];
hvalue=hubbleToday[parcosmopts];
unitsPk=OptionValue[units];
debugPrint[paropts];
{kuf,puf}=pkUnitsConverter[unitsPk,hvalue, externalFile->extfile];
Which[
  psmethod=="CAMB",
    pk=puf*LcdmCambPk[z,k*kuf,FilterRules[paropts,Options@LcdmCambPk]]
  ,
  psmethod=="Transfer",
    pk=puf*LCDMTrHuPre[s8,FilterRules[paropts,Options@LCDMTrHuPre]][z,k*kuf]
  ,
  psmethod=="CosmicEmulator",
    pk=puf*LCDMCosmicEmuPre[z,s8,FilterRules[paropts,Options@LCDMCosmicEmuPre]][k*kuf]
  ,
  psmethod=="ExternalFile",
    pk=puf*externalPowerSpectrumFunction[z,k*kuf,(FilterRules[paropts,Options@externalPowerSpectrumFunction])~Join~({externalSpectraType->"pk"})];
    pk=If[$pks8ratio==1, (*external pk is given as pk/s8^2 *) pk*(sigma8ofZ[z,FilterRules[paropts, Options@sigma8ofZ]])^2, (*external pk is matter power spectrum pk*) pk]
  ,
  psmethod=="ExternalNoWiggle",
    pk=puf*externalPowerSpectrumFunction[z,k*kuf,(FilterRules[paropts,Options@externalPowerSpectrumFunction])~Join~({externalSpectraType->"pknw"})];
    pk=If[$pks8ratio==1, (*external pk is given as pk/s8^2 *) pk*(sigma8ofZ[z,FilterRules[paropts, Options@sigma8ofZ]])^2, (*external pk is matter power spectrum pk*) pk]
  ];
  Return[pk]
];


Options[pkUnitsConverter]=$paramoptions~Join~$pscosmoopts~Join~{externalFile->False}

pkUnitsConverter[unitsstr_String, hval_, opts:OptionsPattern[]]:=Block[
  {unit,kfact, extfile, pfact, h},
  unit=unitsstr;
  h=hval;
  extfile=OptionValue[externalFile];
  If[extfile==False,
    Which[
      unit=="h/Mpc",  (*cosmomathica power spectrum functions are all obtained in "h/Mpc"*)
      kfact=1;
      pfact=1;
      ,
      unit=="1/Mpc",   (*cosmomathica power spectrum functions are all obtained in "h/Mpc"*)
      kfact=(1/h);
      pfact=(1/h^3);
      ]
  ];
  If[extfile==True,
    Which[
      unit==$externalPkUnits,
      kfact=1;
      pfact=1;
      ,
      unit!=$externalPkUnits && $externalPkUnits=="h/Mpc",
      kfact=(1/h);   (*inverse of intuitive k change, since the P function has to be evaluated in the units it was interpolated on*)
      pfact=(1/h^3);
      ,
      unit!=$externalPkUnits && $externalPkUnits=="1/Mpc",
      kfact=(h);  (*inverse of intuitive k change, since the P function has to be evaluated in the units it was interpolated on*)
      pfact=h^3;
    ]
  ];
  {kfact,pfact}
]



zValTozIndex[zval_]:=Block[
  {zindex, zrang=$zrangePowerSpectrum},zindex=0;
  If[zval < Sort[zrang][[1]], Return[zindex=1] ]; (*return the first index, in case the value lies below the first element.*)
  While[zindex+1<=Length@zrang && zrang[[zindex+1]]<= zval, zindex++];
  zindex];


zIndexTozVal[zind_?IntegerQ]:=Module[
  {zval, zrang=$zrangePowerSpectrum},
  zval=$zrangePowerSpectrum[[zind]]
  ];

externalPowerSpectrumFunction::noexternalfile="Attempted to output power spectrum from external file,
but external file was not defined in setExternalCosmoInterpolatingFunctions. Returning: Null"


externalPowerSpectrumFunction::nozrange="Range of redshifts defined for external power spectrum is not a numeric list or an accepted keyword.
 Returning: Null"

Options[externalPowerSpectrumFunction]=$paramoptions~Join~$pscosmoopts~Join~{filesRangeInZ:>$zrangePowerSpectrum,
  getDomain->False, functionInterpolatedInLogLog->False, derivativeInterpolatedInLogLin->False,
  inputNumericalDerivatives->False,
externalFunction->True, externalDerivativeFunction->False, externalSpectraType->"pk",
externalPowerSpectrumNoWiggleInput->False, externalPowerSpectrumNoWiggleDerivativesInput->False};

externalPowerSpectrumFunction[zet_?NumericQ, kps_?NumericQ, opts:OptionsPattern[]]:=Block[{filename,zindex,
zval, paropts , cosmopts, dompos, zrang, grfact, inpunumder, intpFuncloglog, intpDerloglog, funcy=((1#)&), funcx=((1#)&),
extFunc, extDerFunc, extPk, externalSpectraFunc, externalSpectraDerivFunc, typeOpt, extwigg,
returnPk},

  paropts=complementParamValues[{opts},Growth,returnList->"Complement"];
  cosmopts=complementParamValues[{opts},externalPowerSpectrumFunction,returnList->"Complement", filterCosmoPars->True];
  inpunumder=OptionValue[inputNumericalDerivatives];
  intpFuncloglog=OptionValue[functionInterpolatedInLogLog];
  intpDerloglin=OptionValue[derivativeInterpolatedInLogLin];
  typeOpt = OptionValue[externalSpectraType];
  extFunc=OptionValue[externalFunction];
  extwigg=OptionValue[externalPowerSpectrumNoWiggleInput];
  extDerFunc=OptionValue[externalDerivativeFunction];

  Which[typeOpt=="pk",
    externalSpectraFunc=$externalPowerSpectrumInterpolatingFunction;
    externalSpectraDerivFunc = $externalPowerSpectrumDerivativesInterpolatingFunction;
    ,
    typeOpt=="pknw"
    ,
    If[SameQ[extwigg,False],
       Return[$zerofunc[zet,kps]],
       externalSpectraFunc=$externalPowerSpectrumNoWiggleInterpolatingFunction;
       externalSpectraDerivFunc=$externalPowerSpectrumNoWiggleDerivativesInterpolatingFunction];
    (*TO DO: add derivatives here, for the case in which derivatives are provided.*)
    ,
    True,
    externalSpectraFunc=$externalPowerSpectrumInterpolatingFunction;
    externalSpectraDerivFunc = $externalPowerSpectrumDerivativesInterpolatingFunction;
  ];

  debugPrint["cosmopts: "<>ToString@scientificfortranform[cosmopts],1];

  If[OptionValue[setExternalCosmoInterpolatingFunctions,externalPowerSpectrumInput]==False,
    Message[externalPowerSpectrumFunction::noexternalfile];
    Return[0]
  ];

  zrang=OptionValue[filesRangeInZ];
  (*$zrangePowerSpectrum = zrang;*)
  If[intpFuncloglog==True,
    funcy=(10^#)&;
    funcx=Log10[#]&;
  ];

 Which[
  (ListQ[zrang]==True && VectorQ[zrang, NumericQ]==False),
    Message[externalPowerSpectrumFunction::zrangefalse];
    Message[externalPowerSpectrumFunction::nozrange];
    Return[Null]
  ,
  (VectorQ[zrang, NumericQ]==True) && (inpunumder==False),
      zindex=zValTozIndex[zet];
      zval= zrang[[zindex]];
      filename=getCosmoParameterFileName[externalPowerSpectrumFunction,cosmopts];
      If[OptionValue[getDomain]==True,
        getCosmoParameterFileName[externalFunctionTaylor,cosmopts, getFiducialFileName->True]; (*get domain of fiducial file*)
        Return[ funcy@(First@externalSpectraFunc[filename][[zindex]]["Domain"]) ]
        ];
      extPk=funcy@(externalSpectraFunc[filename][[zindex]][funcx@kps]);
    ,
    (VectorQ[zrang, NumericQ]==True) && (inpunumder==True),
        extPk=(externalFunctionTaylor[zet, kps, cosmopts~Join~{functionInterpolatedInLogLog->intpFuncloglog,
               derivativeInterpolatedInLogLin->intpDerloglin,
               externalFunction->externalSpectraFunc,
               externalDerivativeFunction->externalSpectraDerivFunc,
               interpolatedArguments->"zk", kdependence->True}]);
      ];
   Return[extPk]
 ]







numericalDerivativeStep[index_, epsi_]:=Block[{opt1, opt2,step, retlist},
  If[$paramoptions[[index]][[2]]==0.,
    opt1=$paramoptions[[index]]/.pp_?NumericQ->pp+epsi;
    opt2=$paramoptions[[index]]/.pp_?NumericQ->pp-epsi;
    step=2 epsi
    ,
    opt1=$paramoptions[[index]]/.pp_?NumericQ->pp*(1+epsi);
    opt2=$paramoptions[[index]]/.pp_?NumericQ->pp*(1-epsi);
    step=2 epsi*$paramfidus[[index]]
  ];
  retlist={opt1, opt2, step};
  Return[retlist]
]


$vd4[zr_]:=1;
$vd5[zr_]:=1;
$vd7[zr_]:=1;
$vd8[zr_]:=1;


variationOfzdependentQuantities::wrongsign="Wrong function assigned for sign. It has to be either Plus or Minus.";

Options[variationOfzdependentQuantities]={setEpsilon:>$epsilonzstep, varylogquantities->True};

variationOfzdependentQuantities[dvara_, dsign_, epsicoeff_, deropts:OptionsPattern[] ]:=Module[{step, sign, vlog},
  step=OptionValue[setEpsilon]*epsicoeff;
  vlog=OptionValue[varylogquantities];
  sign=If[((dsign==Plus || dsign==Minus)),
    dsign
    (* , Message[variationOfzdependentQuantities::wrongsign]; Abort[] *) ];
  Switch[dvara
    ,
    (*option d5 = Hubble*)
    FisherTools`d5
    ,
    clearZdependentVariations[];
    ClearAll[$vd5];
    $vd5[zr_]:=If[vlog==True, Hubble[zr]^(sign@step), 1+(sign@step) ]; (*variation of logH -> (1+epsilon) logH = log[H^(1+epsilon)]*)
    (*$vd5[zr_]:=Hubble[zr]^(sign@step) (*variation of logH -> (1+epsilon) logH = log[H^(1+epsilon)]*)*)
    ,
    (*option d4 = distance*)
    FisherTools`d4
    ,
    clearZdependentVariations[];
    ClearAll[$vd4];
    $vd4[zr_]:=If[vlog==True, distanceAngular[zr]^(sign@step), 1+(sign@step) ]; (*variation of logH -> (1+epsilon) logD = log[D^(1+epsilon)]*)
    (*$vd4[zr_]:=distanceAngular[zr]^(sign@step)  (*variation of logH -> (1+epsilon) logD = log[D^(1+epsilon)]*)*)
    ,
    (*option d7 = fsigma8*)
    FisherTools`d7
    ,
    clearZdependentVariations[];
    ClearAll[$vd7];
    $vd7[zr_]:=If[vlog==True, fGrowthRateSigma8ofZ[zr]^(sign@step), 1+(sign@step) ];
    ,
    (*option d8 = bsigma8*)
    FisherTools`d8
    ,
    clearZdependentVariations[];
    ClearAll[$vd8];
    $vd8[zr_]:=If[vlog==True, ($biasInterpFunc[zr]sigma8ofZ[zr])^(sign@step), 1+(sign@step) ];
    ,
    _
    ,
    clearZdependentVariations[]
  ];
]

clearZdependentVariations[]:=(
  ClearAll[$vd4];
  ClearAll[$vd5];
  ClearAll[$vd7];
  ClearAll[$vd8];
  $vd4[zr_]:=1;
  $vd5[zr_]:=1;
  $vd7[zr_]:=1;
  $vd8[zr_]:=1;
);


(*$shaParVarInZdep=False;*)


Options[RAPfunction]=$paramoptions~Join~{varyZdependentparameter->{False,Plus,1},
 zdepFunctionsCosmoVariation:>$shaParVarInZdep}
RAPfunction[zred_,mur_,opts:OptionsPattern[]]:=Module[{muref,Href,dref,dvaropt,paroptis, fidoptis,
  zoptis, Rref, qpar, qper, brack},
(*paropts=complementParamValues[{opts},RAPfunction];  Is anyway just called by one function and only has $paramoptions as OptionValues*)
  paroptis=complementParamValues[{opts},RAPfunction, returnList->"Complement", filterCosmoPars->True];
  fidoptis=complementParamValues[{opts},kmuAlcockPaczynski, returnList->"Default", filterCosmoPars->True];


  zoptis=If[$shaParVarInZdep==False, fidoptis, paroptis];

  dvaropt=OptionValue[varyZdependentparameter];
  (*Print[First@dvaropt, Last@dvaropt];*)
  variationOfzdependentQuantities[dvaropt[[1]], dvaropt[[2]], dvaropt[[3]] ];

  qpar=Hubble[zred]/($vd5[zred]*Hubble[zred,zoptis]);
  qper=($vd4[zred]*distanceAngular[zred, zoptis])/(distanceAngular[zred]);
  brack=(1+mur^2*(qper^2/qpar^2 - 1));
  Rref = (1/qper)*Sqrt[brack];
  Return[Rref]
  (*Print["vd5 ",$vd5[zred]];*)
  (*Print["vd4 ",$vd4[zred]];*)
  (* Rref=Sqrt[((($vd5[zred]*Hubble[zred,paroptis])*($vd4[zred]*distanceAngular[zred,paroptis])*mur)^2-(Hubble[zred]*distanceAngular[zred])^2 * (mur^2-1))/
      ((Hubble[zred]*($vd4[zred]*distanceAngular[zred, paroptis]))^2)] *)
]

Options[kmuAlcockPaczynski]=$paramoptions~Join~{varyZdependentparameter->{False,Plus,1},
zdepFunctionsCosmoVariation:>$shaParVarInZdep
};

kmuAlcockPaczynski[zr_,kr_,mur_,opts:OptionsPattern[]]:=Module[{Rap,kap,muap,paroptis,zoptis,
  fidoptis, dvaropt, hubepsval=1, toler=10^(-7) },
(*paropts=complementParamValues[{opts},kmuAlcockPaczynski];*)
  dvaropt=OptionValue[varyZdependentparameter];
  paroptis=complementParamValues[{opts},kmuAlcockPaczynski, returnList->"Full", filterCosmoPars->True];
  fidoptis=complementParamValues[{opts},kmuAlcockPaczynski, returnList->"Default", filterCosmoPars->True];

  zoptis=If[$shaParVarInZdep==False, fidoptis, paroptis];

  variationOfzdependentQuantities[dvaropt[[1]], dvaropt[[2]], dvaropt[[3]] ];

  Rap=RAPfunction[zr,mur,opts];

  If[$internalPkUnits=="1/Mpc",
    hubepsval=(hubble/.paroptis)/(hubble/.fidoptis);  (*this gives 1+/-epsilon *)
    hubepsval = Chop[hubepsval, toler]; (*smooth out very small values of epsilon*)
  ];

  debugPrint["Rap : "<>ToString[Rap],1];
  debugPrint["heps : "<>ToString[hubepsval],1];
  kap=Rap*kr*hubepsval;
  muap=((Hubble[zr,zoptis]*$vd5[zr])/Hubble[zr])*(mur/Rap);
  debugPrint["H :"<>ToString[((Hubble[zr,zoptis]*$vd5[zr])/Hubble[zr])],1];

  Return[{kap,muap}]
]



Options[observedPowerSpectrum]=$paramoptions~Join~$pscosmoopts~Join~{spectrumMethod->"Transfer", sigma8reference->$sigma8reference,
  APeffect->$APeffectBool,kdependentGrowth:>$kdependentGrowth, externalFile->False, varyZdependentparameter->{False,Plus,1},
  zdepFunctionsCosmoVariation->False }
$shaParVarInZdep
observedPowerSpectrum[zred_?NumericQ,k_,mu_,opts:OptionsPattern[]]:=Block[{paroptis, fidoptis, zoptis, linearRec, kAP,muAP,effAP,baoterm, powspec, powopts, kapopts, kaiserterm, s8z, s8zpw, dvar, dvarsign, dvarcoeff, lorentzFoG, gmudamping, powspecDW, powspecNW},
paroptis=complementParamValues[{opts},observedPowerSpectrum, returnList->"Complement", filterCosmoPars->True];
fidoptis=complementParamValues[{opts},observedPowerSpectrum, returnList->"Fiducials", filterCosmoPars->True];
powopts=complementParamValues[{opts},powerSpectrum, returnList->"Complement"];
kapopts=complementParamValues[{opts},kmuAlcockPaczynski, returnList->"Complement"];
effAP=OptionValue[APeffect];
linearRec=(linearBool/.$pscosmoopts);
dvar=(OptionValue[varyZdependentparameter])[[1]];
dvarsign=(OptionValue[varyZdependentparameter])[[2]];
dvarcoeff=(OptionValue[varyZdependentparameter])[[3]];

zoptis=If[$shaParVarInZdep==False, fidoptis, paroptis];

variationOfzdependentQuantities[dvar, dvarsign, dvarcoeff];


debigPrint["Print: d4, d5: "];
debugPrint[$vd4[zred], 0];
debugPrint[$vd5[zred], 0];

{kAP,muAP}=If[effAP,kmuAlcockPaczynski[zred,k,mu,kapopts],{k,mu}];

debugPrint["kap, muap: "<>ToString@kAP<>ToString@muAP, 0];

baoterm=((Hubble[zred,zoptis]*$vd5[zred])*distanceAngular[zred]^2)/(Hubble[zred]*(distanceAngular[zred,zoptis]*$vd4[zred])^2);

debugPrint["baoterm: "<>ToString@baoterm, 0];

kaiserterm=1;
powspec=1;
s8z=1;  (*set dummy values here to avoid nasty debugging*)
lorentzFoG=1;
gmudamping=0; (*when zero, return to standard linear recipe*)

Which[
  $fbsigmaBool==True,

  s8z = sigma8ofZ[zred,fidoptis];  (*sigma8(z) which enters z-dependent variables, independent of shape parameters*)
  s8zpw = sigma8ofZ[zred,paroptis];  (*sigma8(z) that comes with the power spectrum, changes as a function of shape parameters*)

  kaiserterm = ( (($biasInterpFunc[zred]*s8z)*$vd8[zred]) + (fGrowthRateSigma8ofZ[zred,kAP, zoptis]*$vd7[zred] )*muAP^2)^2;

  If[linearRec==False,
    lorentzFoG = FingersOfGod[zred, kAP, muAP, powopts];
    gmudamping = BAOdamping[zred, kAP, muAP, powopts];    (*cosmological parameters are varied here only if external input files exist for those parameters*)
    powspecNW=powerSpectrum[zred,kAP,powopts~Join~{spectrumMethod->"ExternalNoWiggle"}]/(s8zpw^2);
  ];
  powspec=powerSpectrum[zred,kAP,powopts]/(s8zpw^2);
  ,
  $fbsigmaBool==False,
  kaiserterm = $biasInterpFunc[zred]^2*(1+betaRSDfunction[zred,kAP,zoptis]*muAP^2)^2;
  powspec=powerSpectrum[zred,kAP,powopts];
];

debugPrint["s8z: "<>ToString[s8z], 0];
debugPrint["kaiser: "<>ToString[kaiserterm], 0];
debugPrint["lorentz: "<>ToString[lorentzFoG], 0];
debugPrint["gmu: "<>ToString[gmudamping], 0];
debugPrint["pnw: "<>ToString[powspecNW], 0];
debugPrint["pow: "<>ToString[powspec], 0];


powspecDW = powspec;
If[linearRec==False && $fbsigmaBool==True, (*for the moment NL recipe valid only in fbsigma8 case*)
  (*non-linear recipe, damping of wiggles in BAO oscillations*)
  powspecDW = powspec * Exp[- gmudamping * kAP^2] + powspecNW * ( 1-Exp[- gmudamping*kAP^2] );
];


debugPrint["pdw: "<>ToString[powspecDW], 0];

baoterm*powspecDW*lorentzFoG*kaiserterm*errorZ[kAP,muAP,zred,zoptis]
+$PowerExtraShot[zred, k]

];


Options[FingersOfGod]={ignoreSigmaPVCosmoDependence->True}
FingersOfGod[zz_, kk_, muu_, opts:OptionsPattern[{powerSpectrum, FingersOfGod}]]:=Block[{lorentzfog,brack, paroptis,
  fidoptis, zoptis, sgopts, ignS=OptionValue[ignoreSigmaPVCosmoDependence]},
paroptis=complementParamValues[{opts},powerSpectrum, returnList->"Complement", filterCosmoPars->True];
fidoptis=complementParamValues[{opts},powerSpectrum, returnList->"Fiducials", filterCosmoPars->True];
zoptis=If[$shaParVarInZdep==False, fidoptis, paroptis];
sgopts=If[ignS, fidoptis, paroptis];

If[$kdependentGrowth==True,

lorentzfog = 1/( 1 +  (kk * muu * sigmapNLNew[zz,sgopts])^2 )
,
lorentzfog = 1/( 1 +  (kk * muu * fGrowthRateSigma8ofZ[zz,kk,zz]*$vd7[zz] * sigmapNL[zz,paroptis])^2 )
];

Return[lorentzfog]
];



Options[BAOdamping]={ignoreSigmaPVCosmoDependence->True}
BAOdamping[zz_, kk_, muu_, opts:OptionsPattern[{powerSpectrum, BAOdamping}]]:=Block[{gmu, paroptis, fidoptis,
  zoptis, sgopts, ignS=OptionValue[ignoreSigmaPVCosmoDependence]},
paroptis=complementParamValues[{opts},powerSpectrum, returnList->"Complement", filterCosmoPars->True];
fidoptis=complementParamValues[{opts},powerSpectrum, returnList->"Fiducials", filterCosmoPars->True];
zoptis=If[$shaParVarInZdep==False, fidoptis, paroptis];
sgopts=If[ignS, fidoptis, paroptis];

If[$kdependentGrowth==True,

gmu = sigmavNLNew[zz,kk,muu,sgopts]^2;
,
gmu = sigmavNL[zz, paroptis]^2 * ( 1 - muu^2 + muu^2 * (1 + fGrowthRate[zz, kk, zoptis])^2 );
];

Return[gmu]
];


pthetathetaMoments[zz_, mom_?IntegerQ, opts:OptionsPattern[powerSpectrum]]:=pthetathetaMoments[zz,mom, opts]=Block[{paroptis, ptt, fmom},
paroptis=complementParamValues[{opts},powerSpectrum, returnList->"Complement", filterCosmoPars->True];
 If[$kdependentGrowth==True,
  fmom=(fGrowthRate[zz, #, paroptis]^mom &);
  ,
  fmom=(fGrowthRate[zz, paroptis]^mom &);
  ];
	ptt=(1/(6 Pi^2)) * NIntegrate[powerSpectrum[zz, kk, opts]*fmom[kk], {kk, $kminRef, $kmaxRef}];
	Return[ptt]
];



sigmavNLNew[zz_, kk_, muu_, opts:OptionsPattern[powerSpectrum]]:=Block[{sv, f0, f1, f2},
 f0 = pthetathetaMoments[zz,0, opts];
 f1 = pthetathetaMoments[zz,1, opts];
 f2 = pthetathetaMoments[zz,2, opts];
 sv=Sqrt[f1+ 2*muu^2*f1 + muu^2*f2];
 Return[sv]
];

sigmapNLNew[zz_, opts:OptionsPattern[powerSpectrum]]:=sigmapNLNew[zz, opts]=Block[{sp},
 sp=Sqrt[pthetathetaMoments[zz, 1, opts]];
 Return[sp]
];



pthetathetaInt[zz_, opts:OptionsPattern[powerSpectrum]]:=pthetathetaInt[zz, opts]=Block[{optis},
  (1/(6 Pi^2)) * NIntegrate[powerSpectrum[zz, kk, opts], {kk, $kminRef, $kmaxRef}]
];


Options[sigmapNL]={computeNLevolution->False};
(*This is actually sigmap/sigma8, in order to match the fbsigma8 recipe*)
sigmapNL[zz_, opts:OptionsPattern[{powerSpectrum, sigmapNL}]]:=Block[{paroptis, pkoptis, sig8norm=1, zeva, sigp, zmean=FisherTools`$zmeansurvey, integ},
  paroptis=complementParamValues[{opts},powerSpectrum, returnList->"Full", filterCosmoPars->True];
  pkoptis=complementParamValues[{opts},powerSpectrum, returnList->"Full"];
  If[(linearBool/opts)==True, sigmapnl=0; Return[sigmapnl]];
  If[OptionValue[computeNLevolution]==True,
    zeva=zz
    ,
    zeva=zmean];
  If[FilterRules[paroptis,sigmapnl]!={},
    sigp=OptionValue[sigmapnl]
    ,
    integ = pthetathetaInt[zeva, pkoptis];
    sigp= Sqrt[integ]/sigma8ofZ[zeva];
    ];
  Return[sigp]
];

Options[sigmavNL]={computeNLevolution->False, kdependentGrowth:>$kdependentGrowth};
sigmavNL[zk__, opts:OptionsPattern[{powerSpectrum, sigmavNL}]]:=Block[{paroptis, pkoptis, sigvmean,
  zmean=FisherTools`$zmeansurvey, integ, kdep, zval, kval},
  paroptis=complementParamValues[{opts},powerSpectrum, returnList->"Full", filterCosmoPars->True];
  pkoptis=complementParamValues[{opts},powerSpectrum, returnList->"Full"];

  kdep=OptionValue[kdependentGrowth];
  {zval, kval}= twoArgumentsCheck[zk,kdep];
  kval=ReleaseHold[kval];
  If[kdep && NumberQ[$fkfix], kval=$fkfix];

  If[(linearBool/.opts)==True, sigmavnl=0; Return[sigmavnl]];
  If[FilterRules[paroptis,sigmavnl]!={},
    sigvmean=OptionValue[sigmavnl]
    ,
    integ = pthetathetaInt[zmean, pkoptis];
    sigvmean= Sqrt[integ]
  ];
  sigvmean * (Growth[zval]/Growth[zmean])
];

Options[dlnPdkDerivativeFunction]=$paramoptions~Join~$pscosmoopts~Join~{spectrumMethod->"Transfer", sigma8reference->$sigma8reference,
  externalFile->False, kdependentGrowth:>$kdependentGrowth, derivativeMethod->"D"}

dlnPdkDerivativeFunction[zred_?NumericQ,kr_?NumericQ,opts:OptionsPattern[]]:=dlnPdkDerivativeFunction[zred,kr,opts]=Block[{kk,parsfullopts,
parcosmoopts,dlnP,dlnPdkz,extFileBool,kdom, derivm},
parcosmoopts=FilterRules[{opts},Options[powerSpectrum]];
(D[Log[powerSpectrum[zred,kk, parcosmoopts]],kk])/.kk->kr
(*derivm=OptionValue[derivativeMethod];
Which[derivm=="ND",
ND[Log[powerSpectrum[zred,kk, parcosmoopts]],kk,kr, Terms->10], (*Terms=10 assures convergence with D in some tests*)
derivm=="D",
(D[Log[powerSpectrum[zred,kk, parcosmoopts]],kk])/.kk->kr,
True,
(D[Log[powerSpectrum[zred,kk, parcosmoopts]],kk])/.kk->kr
]*)
];


epsilonChangeParameter[par_, epsi_:$epsilonstep, paropts_:$paramoptions]:=(par->(par/.paropts)*(1+epsi));

(*
memofuncs={OmegaM[_Real,___],
OmegaMLCDM[_Real,___],
ndsolDistance[_List],
LCDMDimensionlessHubble[_Real],
extdGNDSolv[_,_Real,_Real,___] ,
extGrowthDz[__Real,_List,___] ,
logGdz[_List] ,
fGrowthSolveND[_List],
rescaleSigma8fromAs[_Real,_Real],
sigma8ofZ[_Real,___],
LCDMCAMBPsPre[_Real,_,___] ,
LCDMTrHuPre[_Real,___] ,
LCDMCosmicEmuPre[_Real,sigma8ref_,___],
powerSpectrum[_Real,k_,___],
RAPfunction[_Real,_,___] ,
kmuAlcockPaczynski[_Real,_,_,___] ,
observedPowerSpectrum[_Real,_,_,___],
dlnPdkDerivativeFunction[_Real,_Real,___]
}




*)

clearCosmologyMemoizedQuantities[opts:OptionsPattern[]]:=Module[{var},

  ParallelTable[removeDownValues[OmegaM[Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[OmegaM[Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[OmegaMLCDM[Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[OmegaMLCDM[Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[ndsolDistance[Except[_Pattern | _PatternTest]]],{i, $KernelCount}];
  removeDownValues[ndsolDistance[Except[_Pattern | _PatternTest]]];
  ParallelTable[removeDownValues[LCDMDimensionlessHubble[Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[LCDMDimensionlessHubble[Except[_Pattern | _PatternTest], ___]];
  ParallelTable[removeDownValues[extdGNDSolv[_,Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[extdGNDSolv[_,Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[logGdz[Except[_Pattern | _PatternTest]]],{i, $KernelCount}];
  removeDownValues[logGdz[Except[_Pattern | _PatternTest]]];
  ParallelTable[removeDownValues[fGrowthSolveND[Except[_Pattern | _PatternTest]]],{i, $KernelCount}];
  removeDownValues[fGrowthSolveND[Except[_Pattern | _PatternTest]]];
  ParallelTable[removeDownValues[fGrowthRate[Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[fGrowthRate[Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[rescaleSigma8fromAs[Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest]]],{i, $KernelCount}];
  removeDownValues[rescaleSigma8fromAs[Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest]]];
  ParallelTable[removeDownValues[sigma8ofZ[Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[sigma8ofZ[Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[LCDMCAMBPsPre[Except[_Pattern | _PatternTest],_,___]],{i, $KernelCount}];
  removeDownValues[LCDMCAMBPsPre[Except[_Pattern | _PatternTest],_,___]];
  ParallelTable[removeDownValues[LCDMTrHuPre[Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[LCDMTrHuPre[Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[LCDMCosmicEmuPre[Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[LCDMCosmicEmuPre[Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[powerSpectrum[Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[powerSpectrum[Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[RAPfunction[Except[_Pattern | _PatternTest],_,___]],{i, $KernelCount}];
  removeDownValues[RAPfunction[Except[_Pattern | _PatternTest],_,___]];
  ParallelTable[removeDownValues[kmuAlcockPaczynski[Except[_Pattern | _PatternTest],_,_,___]],{i, $KernelCount}];
  removeDownValues[kmuAlcockPaczynski[Except[_Pattern | _PatternTest],_,_,___]];
  ParallelTable[removeDownValues[observedPowerSpectrum[Except[_Pattern | _PatternTest],_,_,___]],{i, $KernelCount}];
  removeDownValues[observedPowerSpectrum[Except[_Pattern | _PatternTest],_,_,___]];
  ParallelTable[removeDownValues[dlnPdkDerivativeFunction[Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[dlnPdkDerivativeFunction[Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[Transfer[Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest]]], {i, $KernelCount}];
  removeDownValues[Transfer[Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest]]];
  ParallelTable[removeDownValues[CAMB[Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[CAMB[Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[FrankenEmu[Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest]]], {i, $KernelCount}];
  removeDownValues[FrankenEmu[Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest]]];
]

(*Protect[APeffect,
externalDistanceInput,
externalGrowthInput,
externalGrowthRateInput,
externalHubbleInput,
externalHubbleUnits,
externalPowerSpectrumInput,
fGammaFit,
hubbleReference,
parameterDirectoriesNames,
physicalUnits,
sigma8reference,
spectrumMethod]  (*Protect the symbols used as Options for Package Functions*)*)


Protect[$lightspeed]  (*Protect inmutable constants*)


(*(*Protect definitions of Package functions*)
Protect[
complementParamValues,
curvatureSin,
dampingTerm,
(*DimensionlessHubble,*)
(*distanceAngular,  memoized functions need a special treatment*)
(*dlnPdkDerivativeFunction*)
(*errorZ*)
externalDistanceFunction,
externalHubbleFunction,
externalHubbleUnits,
(*fGrowthRate,*)
getCosmoParameterFileName,
(*Growth,
Hubble,
hubbleToday,
kmuAlcockPaczynski,
LcdmCambPk,
LCDMCAMBPsPre,
LCDMDimensionlessHubble,
LCDMGrowth,
LCDMTrHuPre,
normalizationFactor,
observedPowerSpectrum,
OmegaBaryon0Today,
OmegaM0Today,
OmegaMLCDM,
powerSpectrum,
RAPfunction,
betaRSDfunction,*)
setBiasFunction,
setCosmologyFixedValues,
setExternalCosmoInterpolatingFunctions,
setKandZMinMaxFixedValues,
setParameterOptions,
setPScosmoOptions,
setTermsSwitches
(*sigma8Function,
sigmaR,
sigmaZ,
wDEzFunc,
windowk*)]*)


End[]
EndPackage[]
