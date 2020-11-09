(* ::Package:: *)

BeginPackage["FisherTools`"]

Unprotect["FisherTools`*", "FisherTools`*`*"]  (*Unprotect and Protect definitions for using ClearAll in the case of developing code*)

ClearAll["FisherTools`*", "FisherTools`*`*"]


FisherTools::usage="FisherTools package: This package provides functionality for calculating and analyzing a galaxy clustering Fisher Matrix."

setNumericalEpsilon::usage="set the global value of epsilon for the h step in numerical derivatives.
Value can be set as an option individually for all derivative functions using setEpsilon->epsilon"

epsilonStepForZdependentParameters::usage="Option for setNumericalEpsilon. This sets a different epsilon step for the redshift(z)-dependent parameters."

setFisherDimensions::usage="set dimensions of Fisher Matrix for package functions,
nc=number of z-independent cosmological parameters,
nz=number of z-depdendent cosmological parameter functions,
nzbi=number of z-bins.
Form:setfisherDimensions[nc_,nz_,nzbi_]=Null"

setFisherKmaxValues::usage="setFisherKmaxValues[kmaxhard, kmaxmethod, maxKsigma8:(optional), opts]:
kmaxhard: Sets the maximum k allowed generally, saved as the variable: $kmaxHard:  Can be also set with the option: maxFisherKCut.
kmaxmethod: Sets the method for integration limits in k. Can be also set with the option: fisherKCutMethod.
maxKsigma8: Optional, sets the maximum sigma8 used to calculate the maximum kmax limit for $kmaxMethod 1 and 2.
Can be also set with the option: sigma8KCutValue.
For a display of available kmax calculation methods, type: ?$kmaxMethod"

setFisherParameterFiducials::usage="setFisherParameterFiducials[paramNames_,paramFiducials_]
sets the global default list of rules between parameter names and parameter fiducial values
Not necessary to set if CosmologyFunctions package is loaded and variables are set there"

setGalaxySurveySpecs::usage="
Sets the global lists of number of galaxies for each redshift bin and other specifications for the survey.
See Options[setGalaxySurveySpecs] for the list of possible specifications.
For getting the explanation of a specification option 'spec', enter: ?spec"






$solidAngleInDegrees::usage="Solid angle 4Pi in degrees"





$ZBiasList::usage="List of redshifts and its corresponding fiducial bias"

$zbinGlobal::usage="List with the values of the intervals for the redshift bins of the survey"
$zmeansurvey::usage="Mean redshift of the redshift bins of the survey"
$dzBinWidth::usage="Width of each redshift bin for the galaxy survey"
$areaSurvey::usage="Total sky area in sq deg covered by the survey"
$dndOmdz::usage="number density of galaxies per solid angle in sq. degrees and per redshift bin"

$photoZerr::usage="photometric redshift error of survey"
$intrShear::usage="Intrinsic shear in lensing power spectrum"
$fSky::usage="Covered sky fraction by the survey"
$galDensArcMin2::usage="galaxy density per square arc minute"
$z0GalDist::usage="median galaxy distance for number density of lensing galaxies"

$betaswitch::usage="If set to 0, the cosmological information from the RSD beta functions is turned off. Default: 1"



surveyRedshiftBins::usage="Option for setGalaxySurveySpecs. List with the values of the intervals for the redshift bins of the survey"
dzBinWidth::usage="Option for setGalaxySurveySpecs. Width of each redshift bin for the galaxy survey"
totalSkyAreaSurvey::usage="Option for setGalaxySurveySpecs. Total sky area in sq deg covered by the survey"
zMeanSurvey::usage="Option for setGalaxySurveySpecs. Mean redshift of Survey. If not specified (None) it is calculated from the Mean
of the bins"
photometricRedshiftError::usage="Option for setGalaxySurveySpecs. Photometric redshift error of survey. DefaultValue: $photoZerr=0.05"
lensingIntrinsicShear::usage="Option for setGalaxySurveySpecs. Intrinsic shear in lensing power spectrum"
coveredSkyFraction::usage="Option for setGalaxySurveySpecs. Covered sky fraction by the survey"
galaxyDensityPerArcMinSq::usage="Option for setGalaxySurveySpecs. Galaxy density per square arc minute"
lensingGalaxyDistance::usage="Option for setGalaxySurveySpecs. Median galaxy distance for number density of lensing galaxies"
rsdBetaSwitch::usage="Option for setGalaxySurveySpecs. Turns off the cosmological information of the RSD beta function in the derivatives of P"
dnumbdOmegadz::usage="Option for setGalaxySurveySpecs. Sets the number density of galaxies per solid angle in sq. degrees and per redshift bin"
APeffectSwitch::usage="Option for setGalaxySurveySpecs. Turns off the AP effect in the derivatives of P for the calculation of the Fisher Matrix"


complementOptionValues::usage="Form: complementOptionValues[customOptsList_,internalFunctOptsList_,externalFunctOptsList_, options: returnList->'CustomComplementOptions', filterCosmoPars->'All'].
For a list of custom options, a list of options of an internal function and a list of options for an external function, it returns the following, according to the string set in the option returnList:
CustomFullOptions: gives all the custom options and the default options accepted by the internal function, which are NOT options of the external function given.
CustomComplementOptions: gives only the custom options accepted by the internal function, which are NOT options of the external function given.
ExternalComplementOptions: gives the custom options belonging to the set of options of the external function, which are set to a different value than the defaults of the external function.
ExternalUnsetOptions: gives the custom options belonging to the set of options of the external function, which were not previously defined as options for the internal function.
DefaultOptions: gives back the default unchanged options of the internal function.
Option filterCosmoPars->:
'All': returns the above options output including cosmo parameters.
'In':  only returns from the above options output, the ones which are cosmo parameters.
'Out': returns the above options output, filtering out all cosmo parameters."

$compOptValuesStrList::usage="Option strings for the function complementOptionValues.
{CustomFullOptions,CustomComplementOptions,ExternalComplementOptions,ExternalUnsetOptions,DefaultOptions}";

complementOptionValues::wrstr="Wrong string for option returnList, switching to CustomOptions"

zaverage::usage = "When given a list of initial bin positions, gives back a list of the average quantity in each bin
Form: zaverage[list] = list"
zBinExtendLims::usage="Extend the end of the redshift bin limit by optional value extopt (default 0.2) and extend first bin value to 0"
zbinStart::usage = "Gives back the starting value of the list of bins
Form: zbinStart[list] = number"
zbinEnd::usage = "Gives back the final value of the list of bins
Form: zbinEnd[list] = number"
zAveragetoBins::usage="Gives back the list of redshift bin intervals after being provided with the central redshift of each bin
and the bin width. Form: zAveragetoBins[zlist_,dz_]"
nbins::usage = "Gives back the number of bins in a list of binned values
Form: nbins[list] = number"
nearestBinLimts::usage = "Gives back the lower and upper limits of the bin that contains the point pz.
If pz is smaller than smallest value of list_zbins, it uses 0 as smallest bin limit.
If pz is bigger than biggest number in list_zbins, it returns a list that goes from the last element of list_zbins to pz.
Form: nearestBinLims[pz_, list_zbins] = list of two elements"
centralZbin::usage = "Gives the central value of the bin, when the input is the bin index.
Form: centralZbin[index_number,list_zbins] = number"



fixElements::usage = "fix parameter in Fisher Matrix, by removing row and column associated with the respective index.
Form: fixElements[matrix,list] = matrix"
marginalizeElements::usage = "Marginalize parameters of the Fisher matrix associated with the list of indices given as input.
Form: marginalizeElements[matrix,list] = matrix"
oneSigmaErrors::usage = "1-sigma fully marginalized errors of the Fisher Matrix.
Form: onesigmaErrors[matrix] = list"
oneSigmaMaximizedErrors::usage = "1-sigma fully maximized errors of the Fisher Matrix.
Form: oneSigmaMaximizedErrors[matrix] = list"
insertColumn::usage = "insert vector as a column at position pos.
Form: insertColumn[matrix,vector,pos] = matrix"
insertRow::usage = "insert vector as a row at position pos.
Form: insertRow[matrix,vector,pos] = matrix"
jacobianTransform::usage ="Transform a Matrix by using a Jacobian Matrix: J^t.M.J
Form: jacobianTransform[matrix,jacobianmatrix] = matrix"
toPercentage::usage = "transform error to percentage.
Form: toPercentage[number] = string"
errorsTable::usage ="Generate a table of absolute and relative errors for each fiducial parameter of the Fisher Matrix.
The default row headers are: fiducial, abs.error, rel.error. This can be changed by the option RowHeaders.
 Form: errorsTable[vectorFiducials, vectorOneSigmaErrors, vectorParameterLabels, RowHeaders->vectorStringTitles] = TableForm"

errorsTableTxt::usage ="Generate a table of absolute and relative errors for each fiducial parameter of the Fisher Matrix.
The default row headers are: fiducial, abs.error, rel.error. This can be changed by the option RowHeaders.
 Form: errorsTable[vectorFiducials, vectorOneSigmaErrors, vectorParameterLabels, RowHeaders->vectorStringTitles] = TableForm"


errorsTableBinsPars::usage="Generate a table of absolute and relative errors for redshift dependent parameters of the Fisher Matrix.
Accepts a list with the central values of the bins, the fiducial lists and error lists.
The last accepted list contains the headers labels of the  table.
Form: errorsTableBinsPars[binsList_,fidusTable_,binerrorsTable_,parlabList_,opts:]"

marginalized2Matrix::usage = "Marginalizes Fisher Matrix over every parameter except the ones specified by indices a,b.
Form: marginalized2Matrix[matrix, a_index, b_index] = 2x2 Matrix"

fixed2Matrix::usage=" Fixes Fisher Matrix at every parameter except the ones specified by indices a,b.
Form: fixed2Matrix[matrix_,a_index,b_index] = 2x2 Matrix"

placeBlocks::usage = "Add contribution block matrices to Fisher Matrix by using blocks.
Form: placeBlocks[am,cm,dm,pos,size] = matrix,
Dimensions constraint: Width[cm]==Width[dm] && Length[cm]==Length[am],
am=The upper left corner square matrix,
cm=Cross term rectangular block matrix,
dm=Diagonal square block matrix,
pos=Block position where to set the blocks cm and dm:
dm is set in the diagonal at (pos,pos),
cm is set in the block position (1,pos),
Transpose[cm] is set in the block position (pos,1),
size: Final desired matrix size"

(*confidence contour tools*)
ellipseAngle::usage = "Calculate inclination angle of ellipse w.r.t. the vertical axis, for the confidence region of a 2x2 fisher matrix F_2.
Form: ellipseAngle[matrix]=angle in radians"
drawEllipse::usage="Ellipse given by the Sqrt of the Eigenvalues of the Matrix and rotated
by the angle given by the function ellipseAngle for the parameters specified by indices a and b.
scale parameter scales the area for the different confidence regions.
Form: drawEllipse[matrix, a_index, b_index, fiducialParameter_list, scale] = Graphic"

draw123Ellipses::usage="draws 1-,2- and 3-sigma confidence regions given by the ellipses of the function
drawEllipse using the scalings 1.51, 2.49 and 3.44 respectively. Colors for each respective line are passed as a list
Form: draw123Ellipses[matrix, a_index, b_index, fiducialParameter_list, colors_List] = {Graphic1,Graphic2,Graphic3}
Options (option->default):
lineThickness->0.05  (thickness of ellipses in primter points)
contourSigmas->{1,2,3}  (list of sigma error contours to be plotted)
"

lineThickness::usage="Option for draw123Ellipses to set the absolute thickness of the lines in printer points."

contourSigmas::usage="Option for draw123Ellipses. Default: {1,2,3}. List of sigma contours to be plotted by function draw123Ellipses.
List must be of maximum length = 3, containing any of the numbers between 1 and 3."

styleFunction::usage="Form: styleFunction[text_, size_,opts:]. Accepts all the options of Style. Default: {FontFamily\[Rule]'Times', Bold}"

proper2Subset::usage="Gives the indices of all possible non-repeating subsets of size two from a list of parameters.
Form: proper2Subset[list]=nested list
Tip: Use paramslist[[#]]&/@proper2Subset[paramslist] to obtain the values of the subsets of paramlist"

$blueTones::usage="set of 3 default blue tones for ellipses"



plotConfidenceRegionsGrid::usage="plotConfidenceRegionsGrid[fishmats_List,opts:OptionsPattern[{Graphics, plotConfidenceRegionsGrid, GraphicsGrid}]]:
Takes as input a list with Fisher matrices. Plots the confidence regions for the specified parameters (all parameters by default).
Options with defaults:
{parameterFiducials->{$paramfidus},parameterLabels->{$paramlabels},labelStyleFunction->styFunc,
lineThickness->2.0, contourSigmas->{1,2,3},
colorLists->{$blueTones, $redishTones, $yellowishTones}, legendLabels->{'Matrix 1', 'Matrix 2', 'Matrix 3'},legendPositioning->{1-0.15,1-0.08},whichParameters->{_,_},whichContours->marginalized2Matrix,
Frame->True,AspectRatio->1.,FrameTicksStyle->Directive[Black,16], Spacings->{Scaled[.05],Scaled[.005]},
ImageSize->{400,400}, PlotRangePadding->Scaled[0.05], PlotRange->{Automatic, Automatic}, PlotRangeClipping->True, ImageMargins->10, frameLabelPositions->{{True,False},{True,False}}}, extraEpilog->{}"


parameterFiducials::usage="Option for plotConfidenceRegionsGrid. Default: $paramfidus. Can be a list of several lists of fiducials, one for each Fisher matrix."
parameterLabels::usage="Option for plotConfidenceRegionsGrid. Default: $paramlabels. Specifies the labels for the parameters to be plotted."
labelStyleFunction::usage="Option for plotConfidenceRegionsGrid. Default: styFunc. Specifies a style function for the labels of the plot."
colorLists::usage="Option for plotConfidenceRegionsGrid. Default: $blueTones. List of lists containing three colors for the 1,2 and 3-sigma contours.
There can be one list for each Fisher Matrix."
legendLabels::usage="Option for plotConfidenceRegionsGrid. Default:{}. List of strings, specifying a label for each Fisher Matrix."
extraEpilog::usage="Option for plotConfidenceRegionsGrid. Default:{}. Inset[{'string'}, Scaled[{x,y}]], specifying an extra text going into the Plot figure. Any other object accepted by Epilog could also work."
legendPositioning::usage="Option for plotConfidenceRegionsGrid. Default:{1-0.15,1-0.08}. Specifies the positioning of the legend labels on each subplot in relative units, where {1,1} is the upper right corner."
frameLabelPositions::usage="Option for plotConfidenceRegionsGrid. Default: {{True, False}, {True, False}}. Specifies if Frame labels should be drawn at the corresponding positions {{letf, right}, {bottom, top}}."
whichParameters::usage="Option for plotConfidenceRegionsGrid. Default:{_,_}. Specifies a pattern rule for the parameters to be plotted.
Example: {4,3} only plots 1 ellipse, corresponding to the parameter combination 4-3. {1,_} only plots combinations where parameter 1 appears."
whichContours::usage="Option for plotConfidenceRegionsGrid. Default: marginalized2Matrix. Specifies which kind of contours will be plotted.
 Marginalizing over everything except two parameters: marginalized2Matrix or Fixing over everything except two parameters: fixed2Matrix."




fixPositiveDefiniteMatrix::usage="Fixes the positivity of a matrix by reconstructing the
matrix using its singular values as eigenvalues in the eigensystem. Only to be used for symmetric and real matrices.
Form: fixPositiveDefiniteMatrix[matrix]=matrix "

extractCorner::usage="Extracts first (num x num) elements of the upper left corner of a matrix
Outputs warning if Matrix is not symmetric.
Form: extractCorner[matrix,num]=matrix"

extractErrorAt::usage="Extracts fully marginalized errors at each redshift bin for z-dependent (default) or z-independent parameter
corresponding to the index given.
Form: extractErrorAt[errorlist_,index_,zDependent->True]=list"

parameterPositionInBins::usage="Gives indices of the position inside the Fisher Matrix of the parameter given by index,
index is the number of the parameter among the z-depdendent ones.
Form: parameterPositionInBins[number]=list"


functionOutputImportCheck::usage="Checks if dataFileName exists and imports it as a table of values.
If the file is not found, then the function calcFunc is evaluated with the argument argument_ and then it exports
its results as a table to the file dataFileName.
Form: functionOutputImportCheck[calcFunc_,argument_,dataFileName_]=list (and exporting or importing of files)"


volumeSurvey::usage="gives the volume of a survey in the central value zmid of the redshift bins given by zbino
volumeSurvey[zmid,zbino]=number. Needs to specify globally the covered sky fraction of the survey with Global`volumeEuclid"

totalGal::usage="Global variable counting the total galaxy"

totalGalaxies::usage=" totalGalaxies[numb_,index_] counts total number of Galaxies in survey,
when given number numb and returns cumulative total galaxies in bin index"

PobsSymbolic::usage="PobsSymbolic[kfunc_,mufunc_,Dfunc_,Hfunc_,Dfuncref_,Hfuncref_,Gfunc_,bfunc_,betafunc_,P0func_,Pshotfunc_]
 gives the symbolic form of the observed power spectrum"

RfuncSymbolic::usage="RfuncSymbolic[muref_,Dfunc_,Hfunc_,Dfuncref_,Hfuncref_]
 gives the symbolic R function that relates k and mu fiducial to a different cosmology"

mufuncSymbolic::usage="mufuncSymbolic[muref_,Dfunc_,Hfunc_,Dfuncref_,Hfuncref_]
gives the symbolic mu function that relates mu and mu in the reference cosmology "

kfuncSymbolic::usage="kfuncSymbolic[kref_,muref_,Dfunc_,Hfunc_,Dfuncref_,Hfuncref_]
gives the symbolic k function relating k and k in the reference cosmology"

dmudlogH::usage="derivative of mu function wrt log(Hubble)"
dmudlogD::usage="derivative of mu function wrt log(DistanceAngular)"
dkdlogH::usage="derivative of k function wrt log(Hubble)"
dkdlogD::usage="derivative of k function wrt log(DistanceAngular)"

numericalParamsDerivative::usage="
numericalParamsDerivative[defaultopts_,extraopts_,indexlist_,epsilonValue_,opts]
computes the step, and the parameter options for numerical derivatives of Cosmology dependent functions.
defaultopts=Global parameter options for the function,
extraopts=local parameter options for the function to override the global ones,
indexlist=List giving the index of parameters where to take derivatives from (corresponding to $paramoptions order),
epsilonValue= epsilon step for numerical derivative,
options:
derivativeOrder->First, Second or SecondMixed (passed as strings),
transformDerivativeFunction-> (any function known symbolically by Mathematica, such as Log or Sin)"

epsilonStepForFiducialZero::usage="Option for numericalParamsDerivative. If set to True, the epsilon step for parameters whose fiducial is zero is taken from $epsilonstepfidu0
and it is in principle different from $epsilonstep."

derivativeOrder::usage="Option for numericalParamsDerivative: order of derivatives to be computed"


ndensVolumeDepTable::usage="Gives the number density of galaxies for a parameter dependent volume.
Form ndensVolumeCalcPars[dn,zbino,opts] dn:number count in bin,zbino:redshift bins,
Options: dz=redshift bin width (default 0.1) and cosmological parameters given by $paramoptions"

ndensVolumeInterpFunc::usage="ndensVolumeInterpFunc[dn_,zbino_,opts] calculates an interpolation function for the number density
using ndensVolumeDepTable and the arguments dn: number count, zbino:redshift bins. Options are:
The cosmological parameter options $paramoptions and all options available for Interpolation.
InterpolationOrder set to 1 by default."

dNdensDparTable::usage="dNdensDparTable[dn_,zbino_,alpha_,opts] calculates the derivative of the function ndens wrt the parameter alpha,
 for a number count dn, redshift bins zbino and options which are the cosmological parameters and setEpsilon"

dNdensDparInterpFunc::usage="dNdensDparInterpFunc[dn_,zbino_,alpha_,opts]
returns the interpolation function of the derivative of the function ndens wrt the parameter alpha,
 for a number count dn, redshift bins zbino and options which are the cosmological parameters and
options for the function Interpolation"

ndensLensFunc::usage="Gives the number density function of redshift for Weak Lensing forecasts.
Uses the previously defined $z0GalDist as a parameter"


galaxyNumberIntegrated::usage="Integrated number of galaxies in a redshift shell z2-z1. galaxyNumberIntegrated[z1_,z2_].
This assumes fiducial cosmology and $ndensGalZFuncGlobal as the number density n(z) function"
ndensReBinning::usage="ndensReBinning[newbin_, opts]: From the usual number density function n(z) ($ndensGalZFuncGlobal) and the bins $zbinGlobal,
it calculates a new n(z) in new bins given by newbin.
It uses the integrated number counts of galaxies as a model-independent observable. Options:
resetSurveySpecs: Default=True. If True, calculate and then overwrite the previous settings $ndensGalZFuncGlobal and $zbinGlobal in setGalaxySurveySpecs.
returnInterpolatedInZ: Default=True. If True, return an interpolating function n(z), otherwise return only a list with new n(z) values for the specified z-bins in newbin.
InterpolationOrder: Default=1. Order of the Interpolation function."


returnInterpolatedInZ::usage="Option for ndensReBinning, If True, return an interpolating function n(z), otherwise return only a list with new n(z) values for the specified z-bins in newbin."
resetSurveySpecs::usage="Option for ndensReBinning, If True, calculate and then overwrite the previous settings $ndensGalZFuncGlobal and $zbinGlobal in setGalaxySurveySpecs."



sigma8FunctionForZK::usage="Accepts an interpolation function for powerspectrum[z,k]
and calculates sigma8Function[powS,Pi/2/kvals] on a table of redshifts and k values"

dPowerDpi::usage="Gives the first derivative of the power spectrum wrt the parameter alpha.
Options: setEpsilon->$epsilonstep, transformDerivativeFunction->None, dlnPdpDerivative->False, $paramoptions, $pscosmoopts, spectrumMethod->'Transfer'
Form:dPowerDpi[zred_,k_,alpha_,opts:]"

dObsPowerDpi::usage="Gives the first derivative of the observed power spectrum wrt the parameter alpha.
Options: setEpsilon->$epsilonstep, transformDerivativeFunction->None, dlnPdpDerivative->False, $paramoptions, $pscosmoopts, spectrumMethod->'Transfer', APeffect->True
Form:dObsPowerDpi[zred_,k_,mu_,alpha_,opts:]"

dObsPowerDZpi::usage="Varies the observedPower spectrum with respect to a z-dependent function in a numerical way. Accepts a member of $zdependDerivVector."

d2PowerD2pii::usage="Gives the second derivative of the power spectrum wrt the parameter alpha.
Options: setEpsilon->$epsilonstep, transformDerivativeFunction->None, dlnPdpDerivative->False, $paramoptions, $pscosmoopts, spectrumMethod->'Transfer'
Form:d2PowerD2pii[zred_,k_,alpha_,opts:]"

d2ObsPowerD2pii::usage="Gives the second derivative of the observed power spectrum wrt the parameter alpha.
Options: setEpsilon->$epsilonstep, transformDerivativeFunction->None, dlnPdpDerivative->False, $paramoptions, $pscosmoopts, spectrumMethod->'Transfer', APeffect->True
Form:d2ObsPowerD2pii[zred_,k_,alpha_,opts:]"

d2PowerD2pij::usage="Gives the second mixed derivative of the power spectrum wrt the parameter alpha and beta.
Options: setEpsilon->$epsilonstep, transformDerivativeFunction->None, $paramoptions, $pscosmoopts, spectrumMethod->'Transfer'
Form: d2PowerD2pij[zred_,k_,alpha_,beta_,deropts:]"

d2ObsPowerD2pij::usage="Gives the second mixed derivative of the observed power spectrum wrt the parameter alpha and beta.
Options: setEpsilon->$epsilonstep, transformDerivativeFunction->None, $paramoptions, $pscosmoopts, spectrumMethod->'Transfer', APeffect->True
Form: d2ObsPowerD2pij[zred_,k_,alpha_,beta_,deropts:]"

dBackgroundDpi::usage=" Gives the derivative of the background function, bfunction at z=zred, for the parameter name param.
Options: setEpsilon->$epsilonstep, dLogDerivative\[Rule]True, $paramoptions.
Form: dBackgroundDpi[bfunction_,zred_,param_,optpars]"


dVolumeSurveyDpi::usage="Computes the derivative of the survey volume (volumeSurvey) function with respect to one cosmological parameter.
dVolumeSurveyDpi[zmid_,alpha__Integer,options]. zmid: Redshift at which the derivative is calculated.
alpha: Integer representing the index of the cosmological parameter to derive with. (see $paramoptions).
options->defaults:  1. setEpsilon->$epsilonstep. Sets the epsilon used in the numerical derivative.
					2. zBins->$zbinGlobal. Sets the bins used in the redshift slices of volumeSurvey.
					3. dlnPdpDerivative\[Rule]True. If True, logarithmic derivative is performed."


zBins::usage="Option for certain functions that sets the actual redshift bins used in the calculation.
For normal usage it should be set to the default: $zbinGlobal."


dPobsdmu::usage="dPobsdmu[zr_,k_,mu_]: Derivative dlnPobs/dmu for the observed power spectrum."
$fsigmabBool::usage="Switch between expressions for the derivative of Pobs wrt mu.
If set to True, the derivative dPobsdmu[zr_,k_,mu_] is calculated in terms of fsigma8(z) and bsigma8(z).
If set to False, the derivative dPobsdmu[zr_,k_,mu_] is calculated in terms of betaRSD(z)."


d1::usage=" Log derivative of the observed power spectrum wrt Pshot  (extra shot noise): dlnPobs/dPshot"
d2::usage=" Log derivative of the observed power spectrum wrt Log bias: dlnPobs/dlnb"
d3::usage=" Log derivative of the observed power spectrum wrt Log f (growth rate): dlnPobs/dlnf"
d4::usage=" Log derivative of the observed power spectrum wrt Log D (distance): dlnPobs/dlnD"
d5::usage=" Log derivative of the observed power spectrum wrt Log H (Hubble): dlnPobs/dlnH"
d6::usage=" Log derivative of the observed power spectrum wrt Log G (growth): dlnPobs/dlnG"
d7::usage=" Log derivative of the observed power spectrum wrt Log fsigma8(z) (growth rate times sigma8(z)): dlnPobs/dln(fsigma8(z))"
d8::usage=" Log derivative of the observed power spectrum wrt Log bsigma8(z) (bias times sigma8(z)): dlnPobs/dln(bsigma8(z))"
d9::usage=" Log derivative of the observed power spectrum wrt Log sigmaV(z)^2 (square of the velocity dispersion): dlnPobs/dln(sigmaV(z)^2)"


indexdPs::usage="Index of variable Pshot in Fisher Matrix."
indexdlnb::usage="Index of variable Log bias in Fisher Matrix."
indexdlnf::usage="Index of variable Log f (growth rate) in Fisher Matrix."
indexdlnD::usage="Index of variable Log distance angular in Fisher Matrix."
indexdlnH::usage="Index of variable Log Hubble in Fisher Matrix."
indexdlnG::usage="Index of variable Log G (growth) in Fisher Matrix."
indexdlnfs8::usage="Index of variable fsigma8(z) (growth rate times sigma8(z)) in Fisher Matrix."
indexdlnbs8::usage="Index of variable bsigma8(z) (bias times sigma8(z)) in Fisher Matrix."
indexdlnsigmav::usage="Index of variable sigmaV(z)^2 (square of the velocity dispersion in Fisher Matrix."

Protect[indexdPs,indexdlnb,indexdlnf, indexdlnD, indexdlnH, indexdlnG, indexdlnfs8, indexdlnbs8, indexdlnsigmav]; (*Protected names of indices*)

(*variationOfzdependentQuantities::usage="Function that activates the variation of z-dependent quantities numerically. It accepts any member of $zdependDerivVector."
*)

setFisherDerivativesVector::usage="Sets the vector of the z-dependent intermediate derivative variables $zdependDerivVector in the Fisher Matrix and assigns their position to a Rule $zdependDerivPositions for later use"
$zdependDerivVector::usage="Vector of the z-dependent intermediate derivative variables for the Fisher Matrix. Contains a list of any of the d1-d6 derivatives in the wanted order"
$zdependDerivPositions::usage="Rule assigning one of the indices {indexdP,indexdlnb,indexdlnf, indexdlnD, indexdlnH, indexdlnG} to the corresponding position in $zdependDerivVector"
$zdependDerivRulesFunctions::usage="Vector of rules of the z-dependent intermediate derivative variables for the Fisher Matrix. Contains a list of rules matching the symbols to its corresponding
z-dependent functions"
$zdependderivatives::usage="Switch for the derivatives of the z-dependent functions. Options: 'analytical', 'numerical'. "

$epsilonstep::usage="Global value of epsilon step for numerical derivatives"
$epsilonstepfidu0::usage="Global value of epsilon step for numerical derivatives of parameters whose fiducial values are 0."

$numbGalaxiesGlob::usage="Global list with number of galaxies for each redshift bin"


$galaxiesCountInBins::usage="List with number count of galaxies for each redshift bin,
computed from specifications of dz, dndOmegadz and area of survey"


$numZindependentPars::usage="Number of z-independent cosmological parameters in the Fisher matrix. Default: Length[$paramoptions]"
$numZdependentPars::usage="Number of z-dependent cosmological function parameters in the Fisher matrix. Default: Length[$zdependDerivVector]"
$numZbins::usage="Number of z-bins used in the calculation of the Fisher  Matrix. Default: nbins[$zbinGlobal]"
$numTotalEntries::usage="Number of total entries of the final Fisher Matrix, default: $numTotalEntries=$numZdependentPars*$numZbins+$numZindependentPars"


$ndensGalZFuncGlobal::usage="The global galaxy number density interpolation function used in the Fisher Matrix analysis"

ndensFunction::usage="Option for setGalaxySurveySpecs, shows the global ndens function n(z) ($ndensGalZFuncGlobal) used in the Fisher run."
galaxyCountsInBins::usage="Option for setGalaxySurveySpecs, shows the number of galaxies per bin derived from $ndensGalZFuncGlobal in the Fisher run."
biasFunction::usage="Option for setGalaxySurveySpecs, shows the global bias function b(z) ($biasInterpFunc) used in the Fisher run."

customNdensFunction::usage="Option for setGalaxySurveySpecs, to override the number density function n(z) read from specification files"
customBiasFunction::usage="Option for setGalaxySurveySpecs, to override the bias function b(z) read from specification files"


$relativePkNoise::usage="Function of k and z, which determines an extra noise in the power spectrum,
entering an effective ndens in the covariance"



maxFisherKCut::usage="Option for setFisherKmaxValues: Sets $kmaxHard."
fisherKCutMethod::usage="Option for setFisherKmaxValues: Sets $kmaxMethod.
method for maximum integration limit in Fisher Matrix integration in k. Options:
1: kmax=$kmaxHard is chosen always
2: kmax= Chosen from kmaxInterpolationFunction[z] if smaller than $kmaxHard
3: kmax= Chosen from kmaxInterpolationFunction[z] always"
sigma8KCutValue::usage="Option for setFisherKmaxValues: Sets $maxKsigma8."


$kmaxHard::usage="Maxmimum integration limit in k for Fisher Matrix"

$kmaxMethod::usage="method for maximum integration limit in Fisher Matrix integration in k. Options:
1: kmax=$kmaxHard is chosen always
2: kmax= Chosen from kmaxInterpolationFunction[z] if smaller than $kmaxHard
3: kmax= Chosen from kmaxInterpolationFunction[z] always"

$kmaxInterpFunc::usage="Function which gives the maximum kmax value as a function of redshift, for the k integration of the Fisher Matrix"

$maxKsigma8::usage=" Default=0.35.
Maximum sigma8 which is used to calculate the allowed kmax as a function of redshift in the Fisher Matrix integration."

kmaxZFuncTable::usage="Form: kmaxZFuncTable[sigma8FuncZK_,opts].
Calculates a table of maximum kmax values as a function of redshift for the fisher Matrix Integration.
It uses a constraint on the maximum allowed normalization of P(k), $maxKsigma8,
to calculate the value kmax at each z, by finding the root of the function Sqrt[sigma8(k,z)]==Sqrt[$maxKsigma8].
Other options: rootIniValK=Pi/40.0(Starting point for finding root), zTableSteps=0.1(Steps in z to create final table),
zbinValues=$zbinGlobal(zbins used to calculate the final table)"

kmaxChoice::usage="Function that chooses the maximum limit of integration in k for the Fisher Matrix, according to $kmaxMethod.
For a list of available options, type: ?$kmaxMethod.
Options for kmaxChoice: kmaxHardValue (Default: $kmaxHard), kmaxChoiceMethod: (Default: $kmaxMethod),
kmaxInterpFunction: (Default: $kmaxInterpFunc)"



sigmaPkNoise::usage="Defines an extra error in Pk entering an effective ndens in the covariance,
by multiplying $relativePkNoise with the observedPowerSpectrum. Accepts all the options for observedPowerSpectrum"


ndensEffective::usage="calculate the effective ndens used in the Fisher Matrix calculation
by including or not extra noise affecting the power spectrum. Form: ndensEffective[zref,k,mu,(Options not available yet)]
This noise function has to be set manually using the global function $relativePkNoise[zz,kk].
Options: Deactivate this extra noise by setting: extraErrorNoise->False (Default)."

volumeEffective::usage="calculate the 'effective volume' used in the Fisher Matrix calculation,
using the ndens effective and the observed power spectrum"


exponentVolumeEffective::usage="Option for volumeEffective. This parameter specifies the power of the bracket inside volumeEffective.
Default case: 2, for higher-order tensor purposes this can be set to 3 or 1"

numeratorPobsSwitch::usage="Option for volumeEffective. If this is set to true, volume effective is calculated as (nP/(nP+1)^2.
If set to False, the expression becomes (n/(1+nP))^2.
When set to True, make sure that the derivatives of P wrt cosmological parameters are logarithmic."

FisherIntegration::usage=" Form: FisherIntegration[zref_,fisherBlock_, opts:]. Integrates in k and mu, a block matrix fisherBlock of Fisher integrands evaluated at the redshift zref.
Accepts all available options for NIntegrate with these as default: MaxRecursion->100,PrecisionGoal->7,AccuracyGoal->7,Method->{'GlobalAdaptive','MaxErrorIncreases'->15000}.
Other options:
lower integration k limit: kminIntLimit->$kminRef
zbin limits for calculating survey volume: zbins->$zbinGlobal"

kminIntLimit::usage="Option for FisherIntegration: Lower intergation limit in k for Fisher integrals"

FisherMatrixGC::usage="Defines the integration for the Galaxy Clsutering Fisher Matrix and takes care of the symmetry properties
Options: symmetric->True
Form: FisherMatrixGC[zref_?NumericQ,alpha_,beta_, fisherBlock_,opts:]
Only evaluates if zref is a numeric quantity, to prevent error messages."



FisherMatrixGCCalculation::usage="Calculates the full fisher matrix for galaxy clustering power spectrum in the specified redshift bins.
Form: FisherMatrixGCCalculation[fishCosmoBlock_,fishZBlockDiagonal_:None,fishZBlockCross_:None,fish3Tensor_:None,fish4Tensor_:None,opts].
The arguments are the fisher blocks containing derivatives of the power spectrum P. The z-dependent blocks (diagonal and cross terms) and the blocks containing higher derivatives of P,
are optional arguments.
Options:
zDependentIntermediateBlocks->True (if set to True, calculation of z-dependent fisher blocks is performed)
tensorFisherBlocks->False   (if set to True, calculation of higher order tensor blocks is performed)
resultsDirectory:>$resultsDirectoryGC (global setting for the Notebook, path to location where fisher matrix results are to be saved)
storeTemporaryResults->False  (if set to True, a temporary directory at $resultsDirectoryGC<>'temporaryFiles/' is created where blocks are stored for each redshift bin in the .mc format)
clearFisherIntegrationCache->False (if set to True, all the previous results cached (memoized)
in FisherMatrixGC and FisherTensor are cleared in order to start a new calculation without the need of redifining all parameters,
(This has to be tested more carefully since there might be other memoized values in cosmological functions that are also to be cleared in some cases))"


$vectorCosmo::usage="Vector function of z,k,mu containing derivatives of Pobs wrt the fundamental cosmological parameters"
$vectorZFunctions::usage="Vector function of z,k,mu containing derivatives of Pobs wrt z-dependent cosmological functions (such as bias, Pshot, fGrowth, ...)"
$vectorFullChainSum::usage="Vector function of z,k,mu containing derivatives of Pobs wrt to fundamental cosmological parameters
and wrt z-dependent cosmological functions (usually which do not depend on cosmological parameters), using the chain rule"
FisherBlockBuilder::incop="The length `1` of the zDependent Vector is incompatible with this Fisher calculation method.";
$fisherCosmoFullChainBlock::usage="Matrix built by the Kronecker Product of $vectorFullChainSum with itself, setting to zero
the entries where products of different z-dependent cosmological functions appear."
$fisherDiagonalZTermsBlock::usage="Matrix built by the Kronecker Product of $vectorZFunctions with itself, which will be placed at the diagonal of the Fisher Matrix."
$fisherCrossTermBlock::usage="Matrix built by the Kronecker Product of $vectorCosmo with $vectorZFunctions, which correspond to the cross-correlation terms of the Fisher Matrix."
$fisherCosmoParsBlock::usage="Matrix built by the Kronecker Product of $vectorCosmo with itself, which correspond to the fundamental cosmological parameters in the Fisher Matrix."



logFirstDerivativeOfP::usage="Option for FisherBlockBuilder, controls if the first derivative of P(k) should be a Log derivative, Default=True."

$s8zderivswitch::usage="Switch that turns off (when set to zero) the derivative of Sigma8ofZ in the ChainRule method of FisherBlockBuilder."

FisherBlockBuilder::usage="Given a vector of z-dependent functions, it creates the needed Fisher blocks for a Fisher analysis.
It accepts three different methods:  FullNumerical, ChainRule or Jacobian.
The blocks that are stored are: $fisherDiagonalZTermsBlock, $fisherCrossTermBlock, $fisherCosmoParsBlock or $fisherCosmoFullChainBlock.
The computed vectors are: $vectorCosmo, $vectorZFunctions or $vectorFullChainSum.
For more info on the vectors or blocks use ?$vectorCosmo.
The vectors can also be handed as custom vectors computed outside, by using the options: customVectorCosmoPars, customVectorZFunctions, customVectorFullChainSum.
The option: logFirstDerivativeOfP, controls if the first derivative of P(k) should be a Log derivative, Default=True.
The option: printInformation (Default: True) can be set to False, if no information about the blocks and the output is needed."



printInformation::usage="Option for FisherBlockBuilder. If set to True, information about size of blocks, method used and parameters needed for FisherMatrixGCCalculation are printed. (Default: True)"
customVectorCosmoPar::usage="Option for FisherBlockBuilder. Specify a vector that has been constructed for the cosmo parameters Fisher Block."
customVectorZFunction::usage="Option for FisherBlockBuilder. Specify a vector that has been constructed for the z-dependent function parameters Fisher Block."
customVectorFullChainSum::usage="Option for FisherBlockBuilder. Specify a vector that has been constructed for the full chain rule first derivative of Pobs,
which will be used to build the Fisher Matrix."
fbsigma8Variables::usage="If set to True, derivatives of Pobs are calculated wrt to the variables fsigma8 and bsigma8"
onlyZdependentParameters::usage="If set to True, and the method is ChainRule, then only a Fisher matrix of the z-Dependent Parameters is constructed"

checkSurveySpecsConsistency::usage="Option for FisherMatrixGCCalculation: Before starting calculation, do coonsistency checks on survey specifications"
zDependentIntermediateBlocks::usage="Option for FisherMatrixGCCalculation:  If set to True, calculation of z-dependent fisher blocks is performed.
Used for debugging only, option gets overriden by fisherBlockDerivativeMethod."
tensorFisherBlocks::usage="Option for FisherMatrixGCCalculation: if set to True, calculation of higher order tensor blocks is performed"
resultsDirectory::usage="Option for FisherMatrixGCCalculation: global setting for the Notebook, path to location where fisher matrix results are to be saved."
storeTemporaryResults::usage="Option for FisherMatrixGCCalculation: if set to True, a temporary directory at $resultsDirectoryGC<>'temporaryFiles/' is created
where blocks are stored for each redshift bin in the .mc format"
compressTemporaryResults::usage="Option for FisherMatrixGCCalculation: If storeTemporaryResults is set to True, then with this option
you can choose if the temporary results are stored as compressed .mc files or standard ascii .txt files."
createProgressDialogMonitor::usage="Option for FisherMatrixGCCalculation: If set to True, a Pop-up Dialog will appear, showing the progress of the Fisher calculation."
clearFisherIntegrationCache::usage="Option for FisherMatrixGCCalculation:  if set to True, all the previous results cached (memoized) in FisherMatrixGC and FisherTensor
are cleared in order to start a new calculation without the need of redifining all parameters."
cosmoBlockPkmuDependence::usage="Option for FisherMatrixGCCalculation: If set to True, the dP/dparam cosmo block of the Fisher matrix, depends on z, k and mu.
(Used if P is the observed power spectrum Pobs in this block). If set to False, this block depends only on z, k. (used when P is the raw P(k))
Used for debugging only, option gets overriden by fisherBlockDerivativeMethod."
fisherMatrixExportNameID::usage="Option for FisherMatrixGCCalculation:  String indicating an ID for the result of the Fisher matrix calculation.
This string is used to create files and directories and export quantities."
fisherBlockDerivativeMethod::usage="Option for FisherMatrixGCCalculation: One of three options: FullNumerical, ChainRule, Jacobian.
A fourth option \"Debug\" is available for debugging purposes."
parallelEvaluation::usage="Option for FisherMatrixGCCalculation: Default:True, If set to False, there is no parallelization and ParallelTables are evaluated as Tables."
fisherOrder3TensorBlock::usage="Option for FisherMatrixGCCalculation: The third order tensor integrand S can be given as an option to the FisherMatrixGCCalculation."
fisherOrder4TensorBlock::usage="Option for FisherMatrixGCCalculation: The fourth order tensor integrand Q can be given as an option to the FisherMatrixGCCalculation."
fisherOrder111TensorBlock::usage="Option for FisherMatrixGCCalculation: The third order tensor integrand J can be given as an option to the FisherMatrixGCCalculation."

$fisherMethodCalcBlock::usage="One of three options: FullNumerical, ChainRule, Jacobian. Variable which shows the Fisher method used."


$cosmoblock::usage="Array of Matrix blocks containing the cosmo Fisher blocks"
$crossZblock::usage="Array of Matrix blocks containing the Fisher z-dependent cross terms blocks"
$diagZblock::usage="Array of Matrix blocks containing the Fisher z-dependent terms diagonal blocks"
$Stensorblock::usage="Array of Matrix blocks containing the Fisher 3rd order tensor blocks"
$Qtensorblock::usage="Array of Matrix blocks containing the Fisher 4th order tensor blocks"
$Jtensorblock::usage="Array of Matrix blocks containing the Fisher 3rd order tensor blocks"

InitializeFisherMatrixGCCalculation::usage="Function that wraps other FisherGC calculating functions."
zDependentVariablesVector::usage="Vector of z-dependent variables used in the Fisher matrix."

FisherFinalAssembly::usage="This function uses the results from FisherMatrixGCCalculation, which are: $cosmoblock, $crossZblock and $diagZblock
to construct the final full Fisher Matrix depending on the chosen method specified by $fisherMethodCalcBlock.
It produces the final Fisher Matrix stored in $fisherMatrix and
a reduced Fisher Matrix block $fisherMatrixCosmoPars, corresponding to just the entries of the fundamental cosmological parameters.
Options: $paramoptions, tensorFisherBlocks->False, fisherBlockDerivativeMethod->$fisherMethodCalcBlock,
resultsDirectory:>$resultsDirectoryGC, storeTemporaryResults->False,
fisherMatrixExportNameID->$fisherMatResultID, exportFisherMatrix->True, compressExportMatrix->True"


exportParametersUsed::usage="Option for FisherFinalAssembly. If set to True, a vector containing
{$paramoptions, $paramnames, $paramfidus, $paramlabels} is exported with the $fisherMatResultID name into the $resultsDirectoryGC."


exportFisherMatrix::usage="Option for FisherFinalAssembly. If set to True, the resulting $fisherMatrix will be exported to a file."
compressExportMatrix::usage="Option for FisherFinalAssembly. If exportFisherMatrix->True, the resulting $fisherMatrix will be exported to a file.
If compressExportMatrix->True, the matrix will be exported as a compressed *.mc file.
If compressExportMatrix->False, the matrix will be exported as a table *.txt file."


$fisherMatrix::usage="Final Fisher Matrix computed by FisherFinalAssembly"
$fisherMatrixCosmoPars::usage="Reduced Fisher Matrix block computed by FisherFinalAssembly. Corresponds to just the entries of the fundamental cosmological parameters."


FisherIntTest::usage="Function to test integration"


symmetricMatrix::usage="Option for FisherMatrixGC, specifying if the matrix is to be regarded as symmetric to avoid repeating calculations"
tensorOrder::usage="Option for FisherTensor, specifying the order of the Fisher Tensors"


FisherTensor::usage="Computes higher orders of the Fisher matrix, by constructing tensors of third or fourth order. Options:
symmetric->[True, False], tensorOrder->[3, 4].
Form: FisherTensor[zref,alpha,beta,gamma,delta,fisherBlock,opts]"

FisherTensorS::usage="FisherTensorS[zref_,alpha_,beta_,gamma_,fisherBlockab_,opts]. Function that calculates the Fisher order three tensor formed by P,alpha * P,{beta gamma}"
FisherTensorJ::usage="FisherTensorJ[zref_,a_,b_,c_,fisherBlockab_,opts:]. Function that calculates the Fisher order three tensor formed by P,a * P,b * P,c"
FisherTensorQ::usage="FisherTensorQ[zref_,a_,b_,c_,d_,fisherBlockab_,opts]. Function that calculates the Fisher order four tensor formed by P,ab * P,bc"
prefactor::usage="Option for FisherTensor(S,J,Q). Specifies a prefactor that is multiplied in front of the integral."


fullFisherPostTransformation::usage="Performs transformations on the computed Fisher Matrix.
Form: fullFisherPostTransformation[fishermat_,options]
Options include:
repairPositivity\[Rule]False,
fisherMatrixDimensions\[Rule]{$numZindependentPars,$numZdependentPars,$numZbins},
marginalizedElementsIndex\[Rule]{False},exportTableName\[Rule]'fisherDefaultName', fixedElementsIndex\[Rule]{False},\[IndentingNewLine]performReduceOperations\[Rule]'None', zBinsJacobian\[Rule]False, cosmoParamsJacobian\[Rule]False"


repairPositivity::usage="Option for fullFisherPostTransformation. Default: False. If some eigenvalues of the Fisher matrix are negative but small,
this option allows to make them positive again by doing a singular value decomposition. For details see the function: fixPositiveDefiniteMatrix"
fisherMatrixDimensions::usage="Option for fullFisherPostTransformation. Default: {$numZindependentPars,$numZdependentPars,$numZbins}.
Change this List just if you have a non-standard Fisher matrix."
marginalizedElementsIndex::usage="Option for fullFisherPostTransformation. Default: {False}.
List containing one or several elements of {indexdPs, indexdlnb, indexdlnf, indexdlnD, indexdlnH, indexdlnG}.
The variables chosen, will be marginalized over according to the positions in the matrix set by $zdependDerivPositions"
exportTableName::usage="Option for fullFisherPostTransformation. Default: \"fisherDefaultName\". File name where the final Fisher matrix table will be exported.
The extension .txt will be appended automatically."
fixedElementsIndex::usage="Option for fullFisherPostTransformation. Default: {False}.
List containing one or several elements of {indexdPs, indexdlnb, indexdlnf, indexdlnD, indexdlnH, indexdlnG}.
The variables chosen, will be fixed according to the positions in the matrix set by $zdependDerivPositions"
performReduceOperations::usage="Option for fullFisherPostTransformation. Default: \"None\". Options: \"Marginalize\", \"Fix\".
This option sets the operation to be done on the Fisher matrix following the variables set in marginalizedElementsIndex or fixedElementsIndex"
zBinsJacobian::usage="Option for fullFisherPostTransformation. Default: False. Jacobian matrix that transforms z-bin-dependent functions into cosmological parameters.
Be aware of the correct dimensionality."
cosmoParamsJacobian::usage="Option for fullFisherPostTransformation. Default: False. Jacobian matrix that transforms the fundamental cosmological parameters
into a new set of cosmological parameters. Be aware of the correct dimensionality."


dlnPdpDerivative::usage="Option for derivative of the power spectrum"
setEpsilon::usage="Option for derivative of the power spectrum"
stencilpoints::usage="Option for number of stencil points for the numerical derivative."
$stencilpoints::usage="Number of stencil points for the numerical derivative of z-dependent quantities."

transformDerivativeFunction::usage="Option for derivative of the power spectrum"

dLogDerivative::usage="Option, if set True, caclulates  dLog(bckg)/dp derivatives, where bck are background functions."
dlogPdlogZdepDerivative::usage="Option, if set True, calculates dLogP/dLog(zdep) derivatives of the observed power spectrum, where zdep are z-dependent functions."

extraErrorNoise::usage="Option for ndensEffective to shut down the extra relative noise on the power spectrum P(k)."


CreateProgressDialog::usage=" Create a pop-up window with progress bars. Usage: Give an arbitrary number of triplets, where each triplet is: {Dynamic[variable],maximum_of_variable,'name of progress bar'}.
CreateProgressDialog[{Dynamic[var1_],max1_,'name'}]"


$resultsDirectoryGC::usage="Path to directory, where to store the Fisher Matrix results"
$fisherMatResultID::usage="String to identify and name the resulting Fisher matrix of the actual calculation"


$kamu::usage="Global variable containing integration variables Global`k and/or Global`mu."

clearFisherMemoizedQuantities::usage="Clears memoized quantities stored in functions related to Fisher GC calculations."

Needs["CosmologyFunctions`"];
Needs["UsefulTools`"];
Needs["SurveySpecifications`"];
Needs["FunctionApproximations`"];
EndPackage[]


BeginPackage["FisherTools`", {"CosmologyFunctions`","UsefulTools`","SurveySpecifications`","FunctionApproximations`"}]




Begin["`Private`"]


$solidAngleInDegrees=4 Pi*(180/Pi)^2;


setFisherParameterFiducials[paramNames_,paramFiducials_]:=Module[{chkli},
chkli=Head[paramNames];
If[chkli!=List, Print["Wrong parameter format passed"]];
$paramoptions=Thread[paramNames->paramFiducials];
$paramnames=paramNames;
$paramfidus=paramFiducials;];


$compOptValuesStrList={"CustomFullOptions","CustomComplementOptions","ExternalComplementOptions","ExternalUnsetOptions",
"ExternalChangeDefaultOptions","DefaultOptions"};

Options[complementOptionValues]={filterCosmoPars->"All",returnList->"CustomComplementOptions"};

complementOptionValues[extrapars_,funcopts_,externOpts_,opts:OptionsPattern[]]:=
Block[{extraOpts,internpars,cosmopars,parsout,funcInternOpts,posix, extraextOpts,filtcosmo,stri,funcExternOpts,cosmostring, extraUnionOpts},
extraextOpts=FilterRules[extrapars,externOpts]; (*Take options which are accepted by external function,
even if they were not set in options of internal function*)
funcInternOpts=FilterRules[funcopts,Except[externOpts]];
funcExternOpts=FilterRules[funcopts,externOpts];
extraOpts=FilterRules[extrapars,Except[externOpts]]; (*Filter out options from external built in functions*)
extraOpts=FilterRules[extraOpts,funcopts]; (*Filter out not accepted options from the main function*)
extraOpts=Complement[extraOpts,funcopts]; (*Filter out repeated options (default and extra)*)
(*debugPrint[extraextOpts];
debugPrint[extraOpts];*)
cosmopars=$paramoptions;
stri=OptionValue[returnList];
cosmostring=OptionValue[filterCosmoPars];

If[extraOpts=={},
internpars=funcInternOpts,
posix=Flatten[Position[funcInternOpts[[All,1]],#]&/@extraOpts[[All,1]]];
internpars=ReplacePart[funcInternOpts,Thread[posix->extraOpts]]  (*Replace default function options with given extra options*)
];

If[MemberQ[$compOptValuesStrList, stri]==False,Message[complementOptionValues::wrstr];stri="CustomFullOptions"];

Switch[stri,
"CustomFullOptions",
parsout=internpars,
"CustomComplementOptions",
parsout=extraOpts,
"ExternalComplementOptions",
parsout=Complement[extraextOpts, externOpts],
"ExternalChangeDefaultOptions",
extraUnionOpts=Union[funcExternOpts,extraextOpts];
parsout=If[funcExternOpts==extraUnionOpts, extraUnionOpts, Fold[DeleteCases,extraUnionOpts,funcExternOpts]];
parsout=Complement[parsout, externOpts],
"ExternalUnsetOptions",
parsout=Complement[extraextOpts, funcExternOpts],
"DefaultOptions",
parsout=funcopts];

Switch[cosmostring,
"All",
parsout,
"In",
parsout=FilterRules[parsout,cosmopars],
"Out",
parsout=FilterRules[parsout,Except[cosmopars]]
]
]



setFisherDimensions[nc_,nz_,nzbi_]:=(
$numZindependentPars=nc;
$numZdependentPars=nz;
$numZbins=nzbi;
)



$betaswitch=1;



setGalaxySurveySpecs::inconsist=messageFormatText["Careful: $APeffectBool and $APswitch are not defined consistently.", FontColor->Red];


Options[setGalaxySurveySpecs]={surveyRedshiftBins->$zbinGlobal, dzBinWidth->$dzBinWidth, dnumbdOmegadz->$dndOmdz,
totalSkyAreaSurvey->$areaSurvey, zMeanSurvey->None, customNdensFunction->None, ndensFunction->$ndensGalZFuncGlobal, biasFunction->$biasInterpFunc,
customBiasFunction->None, galaxyCountsInBins->$galaxiesCountInBins,
photometricRedshiftError->0.05,
rsdBetaSwitch->$betaswitch, lensingIntrinsicShear->$intrShear, coveredSkyFraction->$fSky, galaxyDensityPerArcMinSq->$galDensArcMin2,
lensingGalaxyDistance->$z0GalDist, APeffectSwitch->$APswitch}

setGalaxySurveySpecs[opts:OptionsPattern[]]:=Block[{checkvect, custndens=OptionValue[customNdensFunction], custbias=OptionValue[customBiasFunction]},

  $zbinGlobal=OptionValue[surveyRedshiftBins];
  If[Length@$zbinGlobal >= 2, $numZbins=nbins[$zbinGlobal]];
  If[SameQ[OptionValue[zMeanSurvey], None], $zmeansurvey=Mean[zaverage@$zbinGlobal], $zmeansurvey=OptionValue[zMeanSurvey] ];
  $dzBinWidth=OptionValue[dzBinWidth];
  $dndOmdz=OptionValue[dnumbdOmegadz];
  $areaSurvey=OptionValue[totalSkyAreaSurvey];
  $galaxiesCountInBins=($dndOmdz*$dzBinWidth)*$areaSurvey;

  Which[
    (UnsameQ[custndens,None] && UnsameQ[custndens,True])
    ,
    $ndensGalZFuncGlobal=custndens;
    $galaxiesCountInBins = ($ndensGalZFuncGlobal[#]*volumeSurvey[#,$zbinGlobal])&/@zaverage[$zbinGlobal];
    Print["Previous Ndens function and galaxyCounts have been overriden with custom number density."];
    Clear[$dndOmdz];
    SetOptions[setGalaxySurveySpecs, customNdensFunction->True, dnumbdOmegadz->$dndOmdz];
    ,
    SameQ[custndens,None]
    ,
    debugPrint["Numb.dens. function taken from input file."];
    ,
    SameQ[custndens,True]
    ,
    Print["Numb.dens. function was already set to custom function."];
  ];

  Which[(UnsameQ[custbias,None] && UnsameQ[custbias,True]),
    $biasInterpFunc=custbias;
    Print["Previous Bias function has been overriden with custom bias."];
    SetOptions[setGalaxySurveySpecs, customBiasFunction->True];,
    SameQ[custbias,None],
    debugPrint["Bias function taken from input file."];,
    SameQ[custbias,True],
    Print["Bias function was already set to custom function."];
  ];

  SetOptions[setGalaxySurveySpecs, ndensFunction->$ndensGalZFuncGlobal, galaxyCountsInBins->$galaxiesCountInBins,
  biasFunction->$biasInterpFunc, totalSkyAreaSurvey->$areaSurvey, surveyRedshiftBins->$zbinGlobal, zMeanSurvey->$zmeansurvey];

  $betaswitch=OptionValue[rsdBetaSwitch];
  If[Xor[$APeffectBool, $APswitch!=0], Message[setGalaxySurveySpecs::inconsist]];
  $photoZerr=OptionValue[photometricRedshiftError];

  $intrShear=OptionValue[lensingIntrinsicShear];
  $fSky=OptionValue[coveredSkyFraction];
  $galDensArcMin2=OptionValue[galaxyDensityPerArcMinSq];
  $z0GalDist=OptionValue[lensingGalaxyDistance];

]


$zDepDerivativesListFull={d1,d2,d3,d4,d5,d6,d7,d8,d9};
$zDepIndicesListFull={indexdPs,indexdlnb,indexdlnf, indexdlnD, indexdlnH, indexdlnG, indexdlnfs8, indexdlnbs8, indexdlnsigmav};

setFisherDerivativesVector::wrongvect="Vector of z-dependent derivative functions should contain elements of {d1,d2,d3,d4,d5,d6,d7,d8,d9} or be an empty vector {}"


setFisherDerivativesVector[Optional[vec_List, {}]]:=Block[{posi,indlist=$zDepIndicesListFull,tr, derivsList=$zDepDerivativesListFull},

$zdependDerivVector=vec;


If[( Total[Count[derivsList,#]&/@$zdependDerivVector]>=1 ) || $zdependDerivVector=={},
$numZdependentPars=Length@$zdependDerivVector;
posi=Position[vec,#]&/@derivsList;
posi=Flatten[(posi/.{}->{{}}),1];
$zdependDerivPositions=Thread[indlist->posi];
,
Message[setFisherDerivativesVector::wrongvect];
myInterrupt[];
];

$numZindependentPars=Length@$paramoptions;
$numTotalEntries = $numZdependentPars*$numZbins+$numZindependentPars;
printConditionalInformation[Style["Number of z-dependent parameter functions: "<>ToString[$numZdependentPars],16],
Style[", Number of z-independent cosmo parameters: "<>ToString[$numZindependentPars],16],
Style[", Number of z bins: "<>ToString[$numZbins],16],
Style["\nTotal entries of Fisher Matrix: "<>ToString[$numTotalEntries],16]];
printConditionalInformation[Style["\nPositions of z-dependent functions in Fisher Matrix: ",16]];
printConditionalInformation[$zdependDerivPositions];
]


$epsilonstepfidu0=$epsilonstep;
$epsilonzstep=$epsilonstep;
$stencilpoints=5;

Options[setNumericalEpsilon]={epsilonStepForFiducialZero:>$epsilonstepfidu0, epsilonStepForZdependentParameters:>$epsilonzstep};
setNumericalEpsilon[epsi_, opts:OptionsPattern[]]:=Module[{eps, epsZ, epsF0},
  epsF0=OptionValue[epsilonStepForFiducialZero];
  epsZ=OptionValue[epsilonStepForZdependentParameters];
  $epsilonstep=epsi;
  $epsilonzstep=epsZ;
  $epsilonstepfidu0=epsF0;
]


setFisherKmaxValues[kmaxhard_?NumberQ,kmaxmethod_?IntegerQ,maxsigma8_?NumberQ]:=($kmaxHard=kmaxhard; $kmaxMethod=kmaxmethod; $maxKsigma8=maxsigma8;)

$kmaxHard=0.2;
$kmaxMethod=1;
$maxKsigma8=0.35;
setFisherKmaxValues::method="Option fisherKCutMethod must be an Integer."

Options[setFisherKmaxValues]={maxFisherKCut->$kmaxHard, fisherKCutMethod->$kmaxMethod, sigma8KCutValue->$maxKsigma8}
setFisherKmaxValues[opts:OptionsPattern[]]:=Block[{maxk=OptionValue[maxFisherKCut],
  method=OptionValue[fisherKCutMethod],maxks8=OptionValue[sigma8KCutValue], hvalue, unitsfactor},
  If[IntegerQ[method]!=True, Message[setFisherKmaxValues::method]];
  hvalue=hubbleToday[];
  unitsfactor=compatibleUnits[$externalPkUnits,$internalPkUnits, hvalue, inverseDistanceUnits->True];
  $kmaxHard=maxk*unitsfactor;
  $kmaxMethod=method;
  $maxKsigma8=maxks8;
  SetOptions[setFisherKmaxValues,maxFisherKCut->$kmaxHard, fisherKCutMethod->$kmaxMethod, sigma8KCutValue->$maxKsigma8];
  Print[Options[setFisherKmaxValues]];
]



zaverage[zbino_]:=MovingAverage[zbino,2];
zbinStart[zbino_]:=zbino[[1]];
zbinEnd[zbino_]:=zbino[[-1]];
zBinExtendLims[zbino_,extopt_:0.2]:=Append[If[zbino[[1]]>0,Prepend[zbino,0.],zbino],zbinEnd[zbino]+extopt];

zAveragetoBins[zlist_,dz_]:=Block[{dzs=dz/2,zbins},
zbins=N@DeleteDuplicates@Union[Rationalize[zlist+dzs],Rationalize[zlist-dzs]]]


nbins[zbino_]:=Length[zbino]-1;
nearestBinLimts[pz_,zbino_]:=Block[{i1},i1=LengthWhile[zbino,#<pz&];
Switch[i1,0,{0,zbino[[1]]},Length[zbino],{zbino[[-1]],pz},_,zbino[[{i1,i1+1}]]]]
centralZbin[index_,zbino_]:=If[index>=1,zaverage[zbino][[index]],0.];

fixElements[mat_,list_]:=Nest[Transpose[#[[Complement[Range[Length@#],list],All]]]&,mat,2];

marginalizeElements[mat_,list_]:=Inverse[fixElements[Inverse@mat,list]];

oneSigmaMaximizedErrors[mat_]:=Sqrt@(1.0/Diagonal@mat);

oneSigmaErrors[mat_]:=Sqrt@Diagonal@Inverse@mat;

insertColumn[mat_,vec_,pos_]:=Transpose@Insert[Transpose@mat,vec,pos];

insertRow[mat_,vec_,pos_]:=Insert[mat,vec,pos];

jacobianTransform[fmat_,jacmat_]:=Transpose[jacmat].fmat.jacmat;

Options[toPercentage]={digitsPrecision->2};
toPercentage[x_, opts:OptionsPattern[]]:=With[{nn=OptionValue[digitsPrecision]},ToString[NumberForm[x*100,nn]]<>"%"];

errorsTable[fid_,fisherrors_,parlab_,opts:OptionsPattern[{RowHeaders->{"fiducial","abs. error","rel. error"}, digitsPrecision->{5,4}}]]:=
Module[{imsum=Im[Total[fisherrors]], mapForm, digi},
  digi=OptionValue[digitsPrecision];
  mapForm=MapAt[NumberForm[#, digi]&,#,{All}]& ;
If[imsum=!=0,Print["Warning: fisher errors contain complex numbers"]];
TableForm[{mapForm[fid],mapForm[Abs@fisherrors],toPercentage/@Abs[Quiet[fisherrors/fid,{Power::infy}]]},TableHeadings->{OptionValue[RowHeaders],parlab}]]

Options[errorsTableTxt]={digitsPrecision->{5,4}};
errorsTableTxt[fid_,fisherrors_,parlab_,opts:OptionsPattern[]]:=
    Module[{imsum=Im[Total[fisherrors]], mapForm, digi},
      digi=OptionValue[digitsPrecision];
      mapForm=MapAt[NumberForm[#, digi]&,#,{All}]& ;
      If[imsum=!=0,Print["Warning: fisher errors contain complex numbers"]];
      {(ToString[#]&/@$paramnames), ToString[#]&/@mapForm[fid], ToString[#]&/@mapForm[Abs@fisherrors],toPercentage/@Abs[Quiet[fisherrors/fid,{Power::infy}]]}
      ]

errorsTableBinsPars[binsList_,fidusTable_,binerrorsTable_,parlabList_,opts:OptionsPattern[]]:=
Module[{imsum=Im[Total[Flatten@binerrorsTable]]},
If[imsum=!=0,Print["Warning: fisher errors contain complex numbers"]];
TableForm[Transpose@({binsList}~Join~Flatten[Table[{fidusTable[[ii]],binerrorsTable[[ii]],
toPercentage/@(binerrorsTable[[ii]]/fidusTable[[ii]])}, {ii,Length@fidusTable}],1]),
TableHeadings->{None,parlabList}]]


extractCorner[mat_,num_] := Module[{},If[SymmetricMatrixQ[mat]==True,
  Print["symmetric matrix"],
  Print["Non-symmetric matrix"]];
mat[[1;;num,1;;num]]]

extractErrorAt[errlist_,indx_,zopt:OptionsPattern[zDependent->True]]:=Module[
  {element},
  If[IntegerQ[indx]==False,
    Print["index for this parameter does not exist"];
    Abort[];,
    If[indx<0||indx>$numZdependentPars,
      Print["index out of range"];
      Abort[];,
      element=parameterPositionInBins[indx,zopt];
      errlist[[element]]
    ]]
]


parameterPositionInBins[indx_,zopt:OptionsPattern[zDependent->True]]:=
    Block[
      {elem},
      If[OptionValue[zDependent]==True,
        elem=Flatten[Table[{$numZindependentPars+$numZdependentPars*ii-($numZdependentPars-indx)},{ii,$numZbins}]],
        elem={indx}]]


placeBlocks::dimension="Blocks given to the function placeBlocks are not compatible in their dimensions to build a Fisher Matrix"

placeBlocks[am_,cm_,dm_,pos_,size_]:=
Block[{crossWide=Dimensions[cm][[2]],crossLength=Dimensions[cm][[1]],diagSize=Length[dm],cornerSize=Length[am],blockPos,sparseSize},
If[(crossWide==diagSize&&crossLength==cornerSize),
blockPos=(pos-1)*crossWide+2;
sparseSize=size-diagSize-cornerSize+2;
ArrayFlatten[Normal[SparseArray[{{1,1}->w[am],{1,blockPos}->w[cm],
{blockPos,1}->w[Transpose@cm],{blockPos,blockPos}->w[dm]},{sparseSize,sparseSize},0.]]/.w->Sequence],
Message[placeBlocks::dimension];
myInterrupt[]]]

ellipseAngle[mat_]:=Module[{evec,eval},
evec=Eigenvectors@mat;
eval=Eigenvalues@mat;
evec=Table[If[(e[[1]]<0 && e[[2]] < 0)||(e[[1]]>0 && e[[2]] < 0), -e,e],{e,evec}];
First[VectorAngle[#,{1,0}]&/@evec]
];

marginalized2Matrix[fishmat_,a_,b_]:=Inverse@marginalizeElements[fishmat,Complement[Range@Length@fishmat,{a,b}]]

fixed2Matrix[fishmat_,a_,b_]:=Inverse@fixElements[fishmat,Complement[Range@Length@fishmat,{a,b}]]

drawEllipse[mat_,a_,b_,par_,scal_]:=Translate[Rotate[Circle[{0,0},scal*Abs@Sqrt@Eigenvalues@mat],ellipseAngle@mat],par[[{a,b}]]]

Options[draw123Ellipses]={lineThickness->0.05, contourSigmas->{1,2,3}}

draw123Ellipses[mat_,a_,b_,par_,colorslist_,opts:OptionsPattern[]] :=Block[{sigsnums,sigs, list,list1, list2, list3},
list={};
sigsnums={1.51,2.49,3.44};
sigs=OptionValue[contourSigmas];
For[ii=1,ii<Length@sigs + 1,ii++,
If[MemberQ[{1,2,3}, sigs[[ii]]],
list1[sigs[[ii]]]={AbsoluteThickness[OptionValue[lineThickness]],colorslist[[ii]],
           drawEllipse[mat,a,b,par,sigsnums[[sigs[[ii]]]]]}
];
];
list=Flatten[(list~Join~(list1[#]))&/@sigs]
]

proper2Subset[flist_]:=DeleteDuplicates[Sort/@Permutations[Range[Length@flist],{2}]];

Options[styleFunction]={FontFamily->"Times",FontWeight->Bold};
styleFunction[text_, size_,opts:OptionsPattern[Style]]:=Style[text, size, FontFamily->"Times",FontWeight->Bold, opts]

$blueTones={ColorData["Crayola"]["Indigo"],Cyan,ColorData["Crayola"]["TealBlue"]};


Options[plotConfidenceRegionsGrid]={parameterFiducials->{$paramfidus},parameterLabels->{$paramlabels},labelStyleFunction->styFunc,
lineThickness->2.0, contourSigmas->{1,2,3},
colorLists->{$blueTones, $redishTones, $yellowishTones}, legendLabels->{"Matrix 1", "Matrix 2", "Matrix 3"},legendPositioning->{1-0.15,1-0.08},whichParameters->{_,_},whichContours->marginalized2Matrix,
Frame->True,AspectRatio->1.,FrameTicksStyle->Directive[Black,16], Spacings->{Scaled[.05],Scaled[.005]},
ImageSize->{400,400}, PlotRangePadding->Scaled[0.05], PlotRange->{Automatic, Automatic}, PlotRangeClipping->True, ImageMargins->10, frameLabelPositions->{{True,False},{True,False}} , extraEpilog->{}}

plotConfidenceRegionsGrid[fishmats_List,opts:OptionsPattern[{Graphics, plotConfidenceRegionsGrid, GraphicsGrid}]]:=Block[{fishMat,fiducialLists,
labels,colors,thick,graphopts, styleF, fishNumber,
fiduLength, contourFunc, paramslist, paramnums, fullset, parchoose,set, fulllength, graphicsopts, gridopts,
legendstylelist, legpos,legends, extraepilog, colsfirst, fishColFids, contsigs, filtopts, framelab, framefunc, frameposi},

fishNumber=Length@fishmats;
fiducialLists = OptionValue[parameterFiducials];
fiduLength=Length@fiducialLists;

If[fiduLength!=fishNumber, fiducialLists=fiducialLists[[1]]&/@Range[fishNumber]];

colors=OptionValue[colorLists];
legends=OptionValue[legendLabels];

If[
Length@colors!=fishNumber, colors=colors[[1;;fishNumber]]
];
If[Length@legends!=fishNumber, legends=legends[[1;;fishNumber]]
];
styleF=OptionValue[labelStyleFunction];
labels=OptionValue[parameterLabels][[1]];
paramnums=Length@labels;
thick=OptionValue[lineThickness];
contourFunc=OptionValue[whichContours];
fullset=proper2Subset[Range@paramnums];
fulllength=Length@fullset;
parchoose=OptionValue[whichParameters];
contsigs=OptionValue[contourSigmas];
set = {};
Do[If[MemberQ[ss,parchoose[[1]]]&&MemberQ[ss,parchoose[[2]]],AppendTo[set,ss]],{ss,fullset}];
(*Print[set];*)
filtopts=FilterRules[{opts}, Options[plotConfidenceRegionsGrid]];
SetOptions[plotConfidenceRegionsGrid,filtopts];
Print[Options[plotConfidenceRegionsGrid]];
graphicsopts=FilterRules[complementOptionValues[Options[plotConfidenceRegionsGrid],Options[Graphics],
{ImageMargins},returnList->"CustomFullOptions"],Options[plotConfidenceRegionsGrid]];
gridopts = FilterRules[complementOptionValues[Options[plotConfidenceRegionsGrid],Options[GraphicsGrid],
{FrameTicksStyle, BoxForm`ItemAspectRatio,DefaultBaseStyle, PlotRangePadding, AspectRatio,ImageSize},returnList->"CustomFullOptions"],
Options[plotConfidenceRegionsGrid]];
Print[graphicsopts];
Print[gridopts];

frameposi=OptionValue[frameLabelPositions];
framefunc[a_,b_,xy_]:=If[frameposi[[a,b]],framelab[xy],""];

extraepilog = OptionValue[extraEpilog];
debugPrint[extraepilog];

colsfirst=#[[1]]&/@colors;
legendstylelist=Style[#1,12,#2]&@@@Thread[List[legends,colsfirst]];
fishColFids=Thread[List[fishmats,colors,fiducialLists]];
legpos=OptionValue[legendPositioning];
If[Length@set==1,
set=set[[1]];
framelab["y"]=styleF@Rotate[labels[[set[[2]]]], -90 Degree];
framelab["x"]=styleF@labels[[set[[1]]]];

((Graphics[#,Frame->True,FrameLabel-> {{framefunc[1,1,"y"],framefunc[1,2,"y"]},{framefunc[2,1,"x"],framefunc[2,2,"x"]}} , Sequence@@graphicsopts,
Epilog->{Inset[Column@legendstylelist,
Scaled[legpos]], extraepilog}])&@(draw123Ellipses[contourFunc[#1,set[[1]],set[[2]]],set[[1]],set[[2]],#3,#2,lineThickness->thick, contourSigmas->contsigs]&@@@fishColFids))
,
GraphicsGrid[ArrayReshape[Table[(Graphics[#,Sequence@@({FrameLabel->styleF[labels[[{tt[[1]],tt[[2]]}]]]}~Join~graphicsopts),
Epilog-> {Inset[Column@legendstylelist,
Scaled[legpos]], extraepilog}]&@(draw123Ellipses[contourFunc[#1,tt[[1]],tt[[2]]],tt[[1]],tt[[2]],#3,#2,lineThickness->thick, contourSigmas->contsigs]&@@@fishColFids)),{tt,set}],
Reverse@FactorInteger[Length[set]][[All,1]]],Sequence@@gridopts]
]
]




fixPositiveDefiniteMatrix[mat_]:=Block[{eigenvals,eigenmat,diagD,newmat},
eigenmat=Eigensystem[mat][[2]];
eigenvals=Eigensystem[mat][[1]];
diagD=Diagonal[SingularValueDecomposition[mat][[2]]];
If[Max[(Abs@eigenvals-Abs@diagD)]>9*10^-6,
Print["Warning, difference between eigenvalues and singular values is too high"];
Print[Abs@eigenvals-Abs@diagD]];
newmat=Inverse[eigenmat].DiagonalMatrix[diagD].eigenmat
]


functionOutputImportCheck[calcFunc_,argument_,dataFileName_]:= Module[{data,calctime,dataFileOut},
If[Check[data=Import[dataFileName, "Table"],error]==error,
Print["File not found during import, calculating now..."];
calctime=AbsoluteTiming[data=calcFunc[argument];];
Print["Calculation time: ",First@calctime];
dataFileOut=Export[dataFileName,data,"Table"];
Print["Exported: ", dataFileOut]];
If[MatrixQ[data,NumericQ]==True,Print["All entries numeric"];data,
Print["Check:Error: Some entries are not numeric"];{}]
]



Options[volumeSurvey]=$paramoptions
volumeSurvey[zmid_, zbino_,opts:OptionsPattern[]]:=
$areaSurvey*volumeZbin[nearestBinLimts[zmid,zbino][[1]],nearestBinLimts[zmid,zbino][[2]],opts];


totalGalaxies[numb_,index_]:=Block[{totGalBin=0},AppendTo[totalGal,numb];totGalBin=Total[totalGal[[1;;index]]];Return[totGalBin]];

PobsSymbolic[kfunc_,mufunc_,Dfunc_,Hfunc_,Dfuncref_,Hfuncref_,Gfunc_,bfunc_,betafunc_,P0func_,Pshotfunc_]:=
Dfuncref^2*Hfunc/(Dfunc^2*Hfuncref) * Gfunc^2*bfunc^2*(1+betafunc* mufunc^2)^2*P0func[kfunc] + Pshotfunc

RfuncSymbolic[muref_,Dfunc_,Hfunc_,Dfuncref_,Hfuncref_]:=
1/(Dfunc*Hfuncref) Sqrt[Dfunc^2 Hfunc^2 muref^2-Dfuncref^2 Hfuncref^2 (muref^2 - 1)]

mufuncSymbolic[muref_,Dfunc_,Hfunc_,Dfuncref_,Hfuncref_]:=(Hfunc muref)/(Hfuncref RfuncSymbolic[muref,Dfunc,Hfunc,Dfuncref,Hfuncref])

kfuncSymbolic[kref_,muref_,Dfunc_,Hfunc_,Dfuncref_,Hfuncref_]:=RfuncSymbolic[muref,Dfunc,Hfunc,Dfuncref,Hfuncref]*kref



Options[ndensVolumeDepTable]=$paramoptions~Join~$modcosmoopt~Join~{APeffect->True};
ndensVolumeDepTable[dn_,zbino_,opts:OptionsPattern[]]:=Module[{dz,ndensTable,ndensList,ndensfunc},
(*dz=$dzBinWidth;*)
  $dzBinWidth=Table[(zbino[[ii+1]]-zbino[[ii]]),{ii,1,(Length[zbino]-1)}];
ndensTable=Table[(dn[[ii]]*$dzBinWidth[[ii]])/volumeZbin[zbino[[ii]],zbino[[ii+1]],opts],{ii,1,Length[dn]}]
(*volumeZbin has units compatible with the powerSpectrum units *)
]


Options[ndensVolumeInterpFunc]=$paramoptions~Join~$modcosmoopt~Join~{InterpolationOrder->1, APeffect->True};

ndensVolumeInterpFunc[dn_,zbino_,opts:OptionsPattern[{ndensVolumeInterpFunc,Interpolation}]]:=
Module[{ndensTable,ndensList,ndensfunc, interpopts,paropts},
paropts=complementOptionValues[{opts},Options[ndensVolumeInterpFunc],Options[Interpolation],returnList->"CustomComplementOptions"];
interpopts=complementOptionValues[{opts},Options[ndensVolumeInterpFunc],Options[Interpolation],returnList->"ExternalChangeDefaultOptions"];
debugPrint[paropts];
ndensTable=ndensVolumeDepTable[dn,zbino,paropts];
ndensList=Thread[List[zaverage[zbino],ndensTable]];
ndensfunc=Interpolation[ndensList,interpopts]
]


galaxyNumberIntegrated[z1_,z2_]:=Block[{galnumbint, fact},
  galnumbint=($areaSurvey/$solidAngleInDegrees)*NIntegrate[(4 * Pi* (1+zz)^2 distanceAngular[zz]^2*D[((1+zz) distanceAngular[zz]),zz])*$ndensGalZFuncGlobal[zz],{zz,z1,z2}]
(*fix: fix units here in case distanceAngular and the power spectrum don't have compatible units. At the moment it is just taken care of in the volumeZbin*)
]


Options[ndensReBinning]={returnInterpolatedInZ->True, InterpolationOrder->1, resetSurveySpecs->True}


ndensReBinning[newbin_, opts:OptionsPattern[]]:=Block[{zf1,zf2, zbino=newbin, ndensaveraged, ndensinterpfunc},
zf1=nearestBinLimts[#,zbino][[1]]&;
zf2=nearestBinLimts[#,zbino][[2]]&;
ndensaveraged=((galaxyNumberIntegrated[zf1[#],zf2[#]]/(volumeZbin[zf1[#],zf2[#]]*$areaSurvey))&/@(zaverage@zbino));
ndensinterpfunc = Interpolation[Thread[List[zaverage@zbino,ndensaveraged]], InterpolationOrder->1];
If[OptionValue[resetSurveySpecs]==True,
setGalaxySurveySpecs[surveyRedshiftBins->zbino, customNdensFunction->ndensinterpfunc]
];
If[OptionValue[returnInterpolatedInZ]==True,
Return[ndensinterpfunc],
Return[ndensaveraged]
]
]


Options[numericalParamsDerivative]={derivativeOrder->"First",transformDerivativeFunction->"None", epsilonStepForFiducialZero->False};

numericalParamsDerivative[functionsym_,extraopts_,indexlist_Integer,epsilonValue_,opts:OptionsPattern[]]:=
    numericalParamsDerivative[functionsym,extraopts,{indexlist},epsilonValue,opts];

numericalParamsDerivative[functionsym_,extraopts_,indexlist_List,epsilonValue_,opts:OptionsPattern[]]:=Block[
{compx,posix,parsInput,fiduValA,fiduValB,optiPlus,optiMinus,
stepA,stepB,alpha,beta,epsi=epsilonValue,epsi0,returnList,
derivativeType=OptionValue[derivativeOrder],
functionForm=ToString[OptionValue[transformDerivativeFunction]],vartransf,variableFunct,defaultopts,
plusepsiA, minusepsiA, epsistepA, plusepsiB, minusepsiB, epsistepB},
defaultopts=complementParamValues[{extraopts},functionsym,returnList->"Default",filterCosmoPars->True];
parsInput=complementParamValues[{extraopts},functionsym,returnList->"Full",filterCosmoPars->True];
If[OptionValue[epsilonStepForFiducialZero]==False,
  $epsilonstepfidu0=epsi;
  epsi0=epsi;
  ,
  epsi0=$epsilonstepfidu0;
];

Switch[derivativeType,
"First",
alpha=indexlist[[1]];
fiduValA=parsInput[[alpha,2]];
If[functionForm=="None",vartransf=1,
variableFunct=Symbol[functionForm];
vartransf=(D[variableFunct[xx],xx]/.(xx->fiduValA));
];
(*Careful, if fiduValA is equal to 0 then a variable transformation might also be singular*)
If[Chop[fiduValA-0.,10^(-6)]==0.,
  plusepsiA=fiduValA+epsi0; minusepsiA=fiduValA-epsi0; epsistepA=epsi0
,
  plusepsiA=fiduValA*(1+epsi); minusepsiA=fiduValA*(1-epsi); epsistepA=epsi*fiduValA
];

optiPlus=Complement[ReplacePart[parsInput,{{alpha,2}->plusepsiA}],defaultopts];
optiMinus=Complement[ReplacePart[parsInput,{{alpha,2}->minusepsiA}],defaultopts];
stepA=2*epsistepA*vartransf;
returnList={stepA,optiPlus,optiMinus,parsInput},

"Second",
alpha=indexlist[[1]];
fiduValA=parsInput[[alpha,2]];

If[Chop[fiduValA-0.,10^(-6)]==0.,
  plusepsiA=fiduValA+epsi0; minusepsiA=fiduValA-epsi0; epsistepA=epsi0
,
  plusepsiA=fiduValA*(1+epsi); minusepsiA=fiduValA*(1-epsi); epsistepA=epsi*fiduValA
];

optiPlus=Complement[ReplacePart[parsInput,{{alpha,2}->plusepsiA}],defaultopts];
optiMinus=Complement[ReplacePart[parsInput,{{alpha,2}->minusepsiA}],defaultopts];
stepA=(epsistepA)^2;
returnList={stepA,optiPlus,optiMinus,parsInput},

"SecondMixed",
{alpha,beta}=indexlist;
fiduValA=parsInput[[alpha,2]];
fiduValB=parsInput[[beta,2]];

If[fiduValA==0., plusepsiA=fiduValA+epsi0; minusepsiA=fiduValA-epsi0; epsistepA=epsi0
,plusepsiA=fiduValA*(1+epsi); minusepsiA=fiduValA*(1-epsi); epsistepA=epsi*fiduValA];
If[fiduValB==0., plusepsiB=fiduValB+epsi0; minusepsiB=fiduValB-epsi0; epsistepB=epsi0
,plusepsiB=fiduValB*(1+epsi); minusepsiB=fiduValB*(1-epsi); epsistepB=epsi*fiduValB];

optiPlus=Complement[ReplacePart[parsInput,{{alpha,2}->plusepsiA,{beta,2}->plusepsiB}],defaultopts];
optiMinus=Complement[ReplacePart[parsInput,{{alpha,2}->minusepsiA,{beta,2}->minusepsiB}],defaultopts];
stepA=epsistepA;
stepB=epsistepB;
returnList={stepA,stepB,optiPlus,optiMinus,parsInput}
]
]




Options[dNdensDparTable]={setEpsilon->$epsilonstep}~Join~$paramoptions

dNdensDparTable[dn_,zbino_,alpha__Integer,deropts:OptionsPattern[]]:=If[alpha<=Length[$paramoptions],Block[
{optplus,optminus,parused,epsil=OptionValue[setEpsilon],
step, plusfunc,minusfunc},
{step,optplus,optminus,parused}=numericalParamsDerivative[dNdensDparTable,{deropts},List[alpha],epsil];
plusfunc=ndensVolumeDepTable[dn,zbino,optplus];
minusfunc=ndensVolumeDepTable[dn,zbino,optminus];
(plusfunc-minusfunc)/(step)
],
Print["Option index "<>ToString[alpha]<>" not valid"];Abort[]]


Options[dNdensDparInterpFunc]=$paramoptions~Join~{InterpolationOrder->1}


dNdensDparInterpFunc[dn_,zbino_,alpha_,opts:OptionsPattern[{dNdensDparInterpFunc,Interpolation}]]:=
Block[{tableN,dndp,interpopts,paropts},
interpopts=complementOptionValues[{opts},Options[dNdensDparInterpFunc],Options[Interpolation],returnList->"ExternalComplementOptions"];
paropts=complementOptionValues[{opts},Options[dNdensDparInterpFunc],Options[Interpolation],returnList->"CustomComplementOptions"];
tableN=Thread[List[zaverage[zbino],dNdensDparTable[dn,zbino,alpha,paropts]]];
dndp=Interpolation[tableN,interpopts]
]


ndensLensFunc[zz_]:=zz^2Exp[-(zz/($z0GalDist))^(3/2)];


dmudlogH[mu_]:=-mu*(mu^2-1);
dmudlogD[mu_]:=-mu*(mu^2-1);
dkdlogH[k_,mu_]:=k*mu^2;
dkdlogD[k_,mu_]:=k*(mu^2-1);


d1[zr_,k_,mu_]:=1.0/observedPowerSpectrum[zr,k,mu];  (*dlnPobs/dPshot*)

d2[zr_,k_,mu_]:=2.0/($biasInterpFunc[zr] + $biasInterpFunc[zr]*betaRSDfunction[zr,k] mu^2); (*dlnPobs/dlnb*)

d3[zr_,k_,mu_]:=(2 betaRSDfunction[zr,k] mu^2)/(1.0+betaRSDfunction[zr,k] mu^2);(*dlnPobs/dlnf*)

$fbsigmaBool=False;

$zdependderivatives="analytical";


dPobsdmu[zr_,k_,mu_]:=If[$fbsigmaBool==False,
($betaswitch*(4*betaRSDfunction[zr,k]*mu)/(1.0+betaRSDfunction[zr,k]*mu^2)),
(4*fGrowthRateSigma8ofZ[zr,k]*mu)/($biasInterpFunc[zr]*sigma8ofZ[zr]+fGrowthRateSigma8ofZ[zr,k]*mu^2)
];
(*dPobs/dmu in the case of using betaRSD as a variable*)

d4[zr_,k_,mu_]:=Which[$zdependderivatives=="analytical",
  -2.0+If[$APswitch==0, 0,
dPobsdmu[zr,k,mu]*dmudlogD[mu]
+dlnPdkDerivativeFunction[zr,k]*dkdlogD[k,mu]
  ]
  ,
  $zdependderivatives=="numerical",
  dObsPowerDZpi[zr,k,mu,d4,setEpsilon->$epsilonzstep, stencilpoints->$stencilpoints]
];(*dlnPobs/dlnD*)

d5[zr_,k_,mu_]:=Which[$zdependderivatives=="analytical",
    1.0+If[$APswitch==0, 0,
         dPobsdmu[zr,k,mu]*dmudlogH[mu]
         +dlnPdkDerivativeFunction[zr,k]*dkdlogH[k,mu]
    ]
  ,
  $zdependderivatives=="numerical",
  dObsPowerDZpi[zr,k,mu,d5,setEpsilon->$epsilonzstep, stencilpoints->$stencilpoints]
    ];(*dlnPobs/dlnH*)

d9[zr_,k_,mu_]:=-k^2*mu^2*(sigmaV[zr]^2); (*dlnPobs/d(ln(sigmaV^2))) *)

d6[zr_,k_,mu_]:=2.0; (*dlnPobs/dlnG*)

d7[zr_,k_,mu_]:=Which[$zdependderivatives=="analytical",
  (2.0 mu^2 * fGrowthRateSigma8ofZ[zr,k])/($biasInterpFunc[zr]*sigma8ofZ[zr]+fGrowthRateSigma8ofZ[zr,k]*mu^2) (*dlnPobs/dlnfsigma8(z)*)
  ,
  $zdependderivatives=="numerical",
  dObsPowerDZpi[zr,k,mu,d7,setEpsilon->$epsilonzstep, stencilpoints->3]
];

d8[zr_,k_,mu_]:=Which[$zdependderivatives=="analytical",
  (2.0 * $biasInterpFunc[zr]*sigma8ofZ[zr])/($biasInterpFunc[zr]*sigma8ofZ[zr]+fGrowthRateSigma8ofZ[zr,k]*mu^2) (*dlnPobs/dlnbsigma8(z)*)
   ,
  $zdependderivatives=="numerical",
  dObsPowerDZpi[zr,k,mu,d8,setEpsilon->$epsilonzstep, stencilpoints->3]
];

$zdependDerivRulesFunctions={d2-> Unevaluated[Log@($biasInterpFunc[xx])], d3->Unevaluated[Log@(fG[xx])] (*fix this*), d4->Unevaluated[Log@(distanceAngular[xx])],
  d5->Unevaluated[Log@(Hubble[xx])], d7->Unevaluated[Log@(fGrowthRateSigma8ofZ[xx])], d8->Unevaluated[Log@($biasInterpFunc[xx]*sigma8ofZ[xx])]
};

$zdependDictionary={d1->"Ps", d2->"lnb", d3->"lnf", d4->"lnDa", d5->"lnH", d6->"lnG", d7->"lnfs8", d8->"lnbs8", d9->"lnsigmav"};


(*
variationOfzdependentQuantities::wrongsign="Wrong function assigned for sign. It has to be either Plus or Minus.";

Options[variationOfzdependentQuantities]={setEpsilon:>$epsilonzstep};

variationOfzdependentQuantities[dvara_, dsign_, deropts:OptionsPattern[] ]:=Module[{step, sign},
  step=OptionValue[setEpsilon];
  sign=If[((dsign==Plus || dsign==Minus) && MemberQ[$zdependDerivVector, dvara]),
    dsign,
    Message[variationOfzdependentQuantities::wrongsign]; Abort[] ];
  Switch[
    dvara,
    d5,
    ClearAll[$vd5];
    $vd5[zr_]:=Hubble[zr]^(sign@step) (*variation of logH -> (1+epsilon) logH = log[H^(1+epsilon)]*)
    ,
    d4,
    ClearAll[$vd4];
    $vd4[zr_]:=distanceAngular[zr]^(sign@step)  (*variation of logH -> (1+epsilon) logD = log[D^(1+epsilon)]*)
    ,
    _,
    ClearAll[$vd4];
    ClearAll[$vd5];
    $vd4[zr_]:=1;
    $vd5[zr_]:=1;
  ];
]
*)

Options[dObsPowerDZpi]={setEpsilon:>$epsilonzstep, stencilpoints->$stencilpoints, dlogPdlogZdepDerivative->False}~Join~$paramoptions~Join~$pscosmoopts~Join~{spectrumMethod->"Transfer", APeffect->$APeffectBool};

dObsPowerDZpi[zred_?NumericQ,k_,mu_,dvar_,deropts:OptionsPattern[]]:=dObsPowerDZpi[zred,k,mu,dvar,deropts]=If[MemberQ[$zdependDerivVector, dvar],
  Module[{extpowopts,plus,minus,parfids,passopts,fiduF,retu,
  epsil=OptionValue[setEpsilon],step, coeff1, coeff2, coeffepsi},
    Which[
      OptionValue[stencilpoints]==3,
      coeff1=1;
      coeff2=0;
      coeffepsi=2;
      ,
      OptionValue[stencilpoints]==5,
      coeff1=8;
      coeff2=-1;
      coeffepsi=12;
    ];
    step = coeffepsi*epsil;
    fiduF=Evaluate@( (dvar/.$zdependDerivRulesFunctions)/.{xx->zred} ); (*needs to be fixed for future cases where a quantity might be k-dependent*)
    (*dlogP/dlogH = H* ((log(P(H*(1+eps)) - log(P(H*(1-eps)) )/ (2*eps*H) )    *)
    If[OptionValue[dlogPdlogZdepDerivative]==False, step=step*fiduF];

    retu =( (  coeff1*(Log[ observedPowerSpectrum[zred,k,mu, varyZdependentparameter->{dvar,Plus,1}] ] -
         Log[ observedPowerSpectrum[zred,k,mu,varyZdependentparameter->{dvar,Minus,1}] ] ) +
        coeff2*(Log[ observedPowerSpectrum[zred,k,mu, varyZdependentparameter->{dvar,Plus,2}] ] -
            Log[ observedPowerSpectrum[zred,k,mu,varyZdependentparameter->{dvar,Minus,2}] ] )
            ) / (step)
          );
    Return[retu]
    ],
Print["Variable to be varied "<>ToString[dvar]<>" not valid"];
Abort[]
  ]
(*
****deprecated function, remove from here, when tested****
sigma8FunctionForZK[power_]:=Block[{powS},Flatten[Table[Table[powS[k_]:=
power[zz,k];{zz, N@xkmax,sigma8Function[powS,Pi/2/xkmax]},{xkmax,N@Pi/200,N@2Pi,N@Pi/6}],{zz,0.,zbinEnd[zBinExtendLims[$zbinGlobal]],0.1}],1]]
*)


Options[dPowerDpi]={setEpsilon->$epsilonstep,
transformDerivativeFunction->"None",dlnPdpDerivative->False}~Join~$paramoptions~Join~$pscosmoopts~Join~{spectrumMethod->"Transfer",
kdependentGrowth->False};

dPowerDpi[zred_?NumericQ,k_,alpha_?IntegerQ,deropts:OptionsPattern[]]:=dPowerDpi[zred,k,alpha,deropts]=If[alpha<=Length[$paramoptions],
Block[{extpowopts,optplus,optminus,parfids,passopts,fiduP,epsil=OptionValue[setEpsilon],step,transf=OptionValue[transformDerivativeFunction]},
{step,optplus,optminus,parfids}=numericalParamsDerivative[dPowerDpi,{deropts},List[alpha],epsil, transformDerivativeFunction->transf];
extpowopts=complementOptionValues[{deropts},Options[dPowerDpi],Options[powerSpectrum], returnList->"ExternalComplementOptions", filterCosmoPars->"Out"];
If[OptionValue[dlnPdpDerivative],fiduP=powerSpectrum[zred,k,parfids~Join~extpowopts],fiduP=1];
(powerSpectrum[zred,k,optplus~Join~extpowopts]-powerSpectrum[zred,k,optminus~Join~extpowopts])/(step*fiduP)],
Print["Option index "<>ToString[alpha]<>" not valid"];Abort[]]


Options[dObsPowerDpi]={setEpsilon->$epsilonstep,transformDerivativeFunction->"None",
dlnPdpDerivative->False}~Join~$paramoptions~Join~$pscosmoopts~Join~{spectrumMethod->"Transfer", APeffect->$APeffectBool};

dObsPowerDpi[zred_?NumericQ,k_,mu_,alpha_?IntegerQ,deropts:OptionsPattern[]]:=dObsPowerDpi[zred,k,mu,alpha,deropts]=If[alpha<=Length[$paramoptions],
Block[{extpowopts,optplus,optminus,parfids,passopts,fiduP,epsil=OptionValue[setEpsilon],step,transf=OptionValue[transformDerivativeFunction]},
{step,optplus,optminus,parfids}=numericalParamsDerivative[dObsPowerDpi,{deropts},List[alpha],epsil, transformDerivativeFunction->transf];
extpowopts=complementOptionValues[{deropts},Options[dObsPowerDpi],Options[observedPowerSpectrum], returnList->"ExternalComplementOptions",
filterCosmoPars->"Out"];
If[OptionValue[dlnPdpDerivative],fiduP=observedPowerSpectrum[zred,k,mu,parfids~Join~extpowopts],fiduP=1];
(observedPowerSpectrum[zred,k,mu,optplus~Join~extpowopts]-observedPowerSpectrum[zred,k,mu,optminus~Join~extpowopts])/(step*fiduP)],
Print["Option index "<>ToString[alpha]<>" not valid"];Abort[]]


Options[d2PowerD2pii]={setEpsilon->$epsilonstep,transformDerivativeFunction->"None",dlnPdpDerivative->False}~Join~$paramoptions~Join~$pscosmoopts~Join~{spectrumMethod->"Transfer"};

d2PowerD2pii[zred_?NumericQ,k_,alpha_?IntegerQ,deropts:OptionsPattern[]]:=d2PowerD2pii[zred,k,alpha,deropts]=If[alpha<=Length[$paramoptions],
Block[{extpowopts,optplus,optminus,parfids,passopts,fiduP,epsil=OptionValue[setEpsilon],step,transf=OptionValue[transformDerivativeFunction]},
{step,optplus,optminus,parfids}=numericalParamsDerivative[d2PowerD2pii,{deropts},List[alpha],epsil,derivativeOrder->"Second"];
extpowopts=complementOptionValues[{deropts},Options[d2PowerD2pii],Options[powerSpectrum], returnList->"ExternalComplementOptions", filterCosmoPars->"Out"];
If[OptionValue[dlnPdpDerivative],fiduP=powerSpectrum[zred,k,parfids~Join~extpowopts],fiduP=1];
(powerSpectrum[zred,k,optplus~Join~extpowopts]+powerSpectrum[zred,k,optminus~Join~extpowopts]-2*powerSpectrum[zred,k,parfids~Join~extpowopts])/(step*fiduP)],
Print["Option index "<>ToString[alpha]<>" not valid"];Abort[]]

Options[d2PowerD2pij]={setEpsilon->$epsilonstep,transformDerivativeFunction->"None"}~Join~$paramoptions~Join~$pscosmoopts~Join~{spectrumMethod->"Transfer"};

d2PowerD2pij[zred_?NumericQ,k_,a_?IntegerQ,b_?IntegerQ,deropts:OptionsPattern[]]:=d2PowerD2pij[zred,k,a,b,deropts]=
Block[{extpowopts,alpha=a,beta=b,optplus,optminus,parfids,passopts,fiduP,epsil=OptionValue[setEpsilon],stepA,stepB,
transf=OptionValue[transformDerivativeFunction]},
If[alpha<=Length[$paramoptions]&&beta<=Length[$paramoptions],
If[alpha==beta,d2PowerD2pii[zred,k,alpha,deropts],
If[alpha<beta,d2PowerD2pij[zred,k,beta,alpha,deropts],
{stepA,stepB,optplus,optminus,parfids}=numericalParamsDerivative[d2PowerD2pij,{deropts},List[alpha,beta],epsil,derivativeOrder->"SecondMixed"];
extpowopts=complementOptionValues[{deropts},Options[d2PowerD2pij],Options[powerSpectrum], returnList->"ExternalComplementOptions", filterCosmoPars->"Out"];
(((powerSpectrum[zred,k,optplus~Join~extpowopts]+powerSpectrum[zred,k,optminus~Join~extpowopts]-2*powerSpectrum[zred,k,parfids~Join~extpowopts])/(2*stepA*stepB))-
(stepA/(2*stepB))*d2PowerD2pii[zred,k,alpha,deropts]-(stepB/(2*stepA))*d2PowerD2pii[zred,k,beta,deropts])
]],Print["Derivative index alpha or beta not valid"];Abort[]]]


Options[d2ObsPowerD2pii]={setEpsilon->$epsilonstep,transformDerivativeFunction->"None",dlnPdpDerivative->False}~Join~$paramoptions~Join~$pscosmoopts~Join~{spectrumMethod->"Transfer",APeffect->$APeffectBool};

d2ObsPowerD2pii[zred_?NumericQ,k_,mu_,alpha_?IntegerQ,deropts:OptionsPattern[]]:=d2ObsPowerD2pii[zred,k,mu,alpha,deropts]=If[alpha<=Length[$paramoptions],
Block[{extpowopts,optplus,optminus,parfids,passopts,fiduP,epsil,step,transf},

transf=OptionValue[transformDerivativeFunction];
epsil=OptionValue[setEpsilon];

{step,optplus,optminus,parfids}=numericalParamsDerivative[d2ObsPowerD2pii,{deropts},List[alpha],epsil,derivativeOrder->"Second"];
extpowopts=complementOptionValues[{deropts},Options[d2ObsPowerD2pii],Options[observedPowerSpectrum], returnList->"ExternalComplementOptions", filterCosmoPars->"Out"];
If[OptionValue[dlnPdpDerivative],fiduP=observedPowerSpectrum[zred,k,mu,parfids~Join~extpowopts],fiduP=1];
(observedPowerSpectrum[zred,k,mu,optplus~Join~extpowopts]+observedPowerSpectrum[zred,k,mu,optminus~Join~extpowopts]-2*observedPowerSpectrum[zred,k,mu,parfids~Join~extpowopts])/(step*fiduP)],
Print["Option index "<>ToString[alpha]<>" not valid"];Abort[]]


Options[d2ObsPowerD2pij]={setEpsilon->$epsilonstep,transformDerivativeFunction->"None"}~Join~$paramoptions~Join~$pscosmoopts~Join~{spectrumMethod->"Transfer",APeffect->$APeffectBool};

d2ObsPowerD2pij[zred_?NumericQ,k_,mu_,a_?IntegerQ,b_?IntegerQ,deropts:OptionsPattern[]]:=d2ObsPowerD2pij[zred,k,mu,a,b,deropts]=
Block[{extpowopts,alpha=a,beta=b,optplus,optminus,parfids,passopts,fiduP,epsil,stepA,stepB,
transf},
  epsil=OptionValue[setEpsilon];
  transf=OptionValue[transformDerivativeFunction];
If[alpha<=Length[$paramoptions]&&beta<=Length[$paramoptions],
If[alpha==beta,d2ObsPowerD2pii[zred,k,mu,alpha,deropts],
If[alpha<beta,d2ObsPowerD2pij[zred,k,mu,beta,alpha,deropts],
{stepA,stepB,optplus,optminus,parfids}=numericalParamsDerivative[d2ObsPowerD2pij,{deropts},List[alpha,beta],epsil,derivativeOrder->"SecondMixed"];
extpowopts=complementOptionValues[{deropts},Options[d2ObsPowerD2pij],Options[observedPowerSpectrum], returnList->"ExternalComplementOptions", filterCosmoPars->"Out"];
(((observedPowerSpectrum[zred,k,mu,optplus~Join~extpowopts]+observedPowerSpectrum[zred,k,mu,optminus~Join~extpowopts]-2*observedPowerSpectrum[zred,k,mu,parfids~Join~extpowopts])/(2*stepA*stepB))-
(stepA/(2*stepB))*d2ObsPowerD2pii[zred,k,mu,alpha,deropts]-(stepB/(2*stepA))*d2ObsPowerD2pii[zred,k,mu,beta,deropts])
]],Print["Derivative index alpha or beta not valid"];Abort[]]]



Options[dBackgroundDpi]={setEpsilon->$epsilonstep,dLogDerivative->True}~Join~$paramoptions
dBackgroundDpi[bfunction_,zredk__?NumericQ,param_,optpars:OptionsPattern[]]:=
dBackgroundDpi[bfunction,zredk,param,optpars]=
Block[{index,optplus,optminus,fiduval,step,parPlus,parMinus,parfids,alphaind,fiduBack, cosmopts,
epsi},
  epsi=OptionValue[setEpsilon];
  alphaind=Position[$paramoptions,param][[1,1]];
  cosmopts=FilterRules[{optpars},$paramoptions];
  {step,optplus,optminus,parfids}=numericalParamsDerivative[dBackgroundDpi,{cosmopts},List[alphaind],epsi];
  If[OptionValue[dLogDerivative],fiduBack=bfunction[zredk,parfids],fiduBack=1];
  (bfunction[zredk,optplus]-bfunction[zredk,optminus])/(step*fiduBack)
]

(*triplet has to be formed by a dynamic variable, a number and a string*)
CreateProgressDialog[triplet__]:=
CreateDialog[Column[{"Calculation Progress:"}~Join~(Row[{ProgressIndicator[#1,{0,#2}]," ",
(#1/#2)," % " <>ToString@#3}]&@@@{triplet})],
WindowSize->Fit,WindowTitle->"Calculation Progress Bars "]


Options[dVolumeSurveyDpi]={setEpsilon->$epsilonstep, zBins->$zbinGlobal, dlnPdpDerivative->True}~Join~$paramoptions;

dVolumeSurveyDpi[zmid_,alpha__Integer,deropts:OptionsPattern[]]:=If[alpha<=Length[$paramoptions],Block[
{optplus,optminus,parused,epsil=OptionValue[setEpsilon],
step, plusfunc,minusfunc, fidu,zbino=OptionValue[zBins]},
{step,optplus,optminus,parused}=numericalParamsDerivative[dVolumeSurveyDpi,{deropts},List[alpha],epsil];
    plusfunc=volumeSurvey[zmid, zbino, optplus];
    minusfunc=volumeSurvey[zmid, zbino, optminus];
If[OptionValue[dlnPdpDerivative],fidu=volumeSurvey[zmid, zbino, parused],fidu=1];
       (plusfunc-minusfunc)/(step*fidu)
    ],
Print["Option index "<>ToString[alpha]<>" not valid"];Abort[]
]


Options[kmaxChoice]={kmaxHardValue:>$kmaxHard, kmaxChoiceMethod:>$kmaxMethod, kmaxInterpFunction->$kmaxInterpFunc};

kmaxChoice[zb_, opts:OptionsPattern[]]:=Block[
{kmaxhard=OptionValue[kmaxHardValue],ktemp, unitsfactor=1 (*units taken care of at each method*), kmaxmeth=OptionValue[kmaxChoiceMethod],kmaxFunc=OptionValue[kmaxInterpFunction],kout, pf, kf},
  Switch[kmaxmeth,
   1,kout=kmaxhard,
   2,ktemp=kmaxFunc[zb];If[ktemp>kmaxhard,kout=kmaxhard,kout=ktemp],
   3,kout=kmaxFunc[zb]
  ];
  kout
  ];

Options[kmaxZFuncTable]={maximumSigma8->$maxKsigma8, rootIniValK->Pi/40.0, zTableSteps->0.1, zbinValues->$zbinGlobal};

kmaxZFuncTable[sigma8FuncZK_,opts:OptionsPattern[]]:=Block[{kTab,maxs8=OptionValue[maximumSigma8],
root0k=OptionValue[rootIniValK],zsteps=OptionValue[zTableSteps], zbins=OptionValue[zbinValues]},
kTab=Table[{zz,FindRoot[Sqrt[sigma8FuncZK[zz,xk]]==Sqrt[maxs8],{xk,root0k}][[1,2]]},{zz,0.,zbinEnd[zBinExtendLims[zbins]],zsteps}]
];


Options[sigmaPkNoise]={};

sigmaPkNoise[zre_,kre_, mure_, opts:OptionsPattern[{sigmaPkNoise,observedPowerSpectrum}]]:=Block[{opti},
opti=OptionValue[opts];
$relativePkNoise[kre,zre]*observedPowerSpectrum[zre,kre,mure,opti]
];

$extraErrorNoise=False;

Options[ndensEffective]={extraErrorNoise->$extraErrorNoise};

ndensEffective[zref_,k_, mu_, opts:OptionsPattern[{ndensEffective,observedPowerSpectrum}]]:=Block[{extnoise, ndensout},
extnoise=OptionValue[extraErrorNoise];
If[extnoise==True,
ndensout=$ndensGalZFuncGlobal[zref]/(1+$ndensGalZFuncGlobal[zref]*sigmaPkNoise[zref,k,mu]),
ndensout=$ndensGalZFuncGlobal[zref]]]

Options[volumeEffective]={extraErrorNoise->$extraErrorNoise, exponentVolumeEffective->2, numeratorPobsSwitch->True};

volumeEffective[zref_,kref_, muref_, opts:OptionsPattern[]]:=Block[{ndens, Pobs, expo, numsw},
expo=OptionValue[exponentVolumeEffective];
numsw=OptionValue[numeratorPobsSwitch];  (*standard set to True, then derivatives of P wrt parameters must be logarithmic*)
ndens=ndensEffective[zref,kref, muref];  (*number density and Pobs must have compatible units.*)
Pobs=observedPowerSpectrum[zref,kref,muref];
(ndens*If[numsw,Pobs,1]/(1+ndens*Pobs))^expo
]



$s8zderivswitch=1;

Options[FisherBlockBuilder]=$paramoptions~Join~{
  tensorFisherBlocks->False,
  fisherBlockDerivativeMethod:>$fisherMethodCalcBlock,
  customVectorCosmoPar->None,
  customVectorZFunction->None,
  customVectorFullChainSum->None,
  logFirstDerivativeOfP->True,
  fbsigma8Variables->False,
  onlyZdependentParameters->False,
  printInformation->True
}


FisherBlockBuilder::wrmethj="Method Jacobian assumes $zdependDerivVector != {}. Use FullNumerical or ChainRule otherwise."


FisherBlockBuilder[zdepVector_,opts:OptionsPattern[]]:=Module[{vectorCosmo, vectorCross, vectorZFunctions,
fishBlockCosmo, fishBlockCross, fishBlockZFunctions, fishBlockFullChain,
method=OptionValue[fisherBlockDerivativeMethod], deltaijVect, deltaij,
custvecCosmo=OptionValue[customVectorCosmoPar],custvecZf=OptionValue[customVectorZFunction],custvecChS=OptionValue[customVectorFullChainSum],
fbs8Opt=OptionValue[fbsigma8Variables], prinf=OptionValue[printInformation]},

SetOptions[printConditionalInformation, printBool->prinf];
setFisherDerivativesVector[zdepVector];
(*this functions sets the dimensions of the final Fisher Matrix (and blocks) and sets the vector $zdependDerivVector.
(*setFisherDerivativesVector works also if $zdependDerivVector={} *) *)
$fisherMethodCalcBlock = method;
If[fbs8Opt==True, $fbsigmaBool=True];
If[$fbsigmaBool==True && $pks8ratio==1, $s8zderivswitch=0]; (*If the fs8 and bs8 variables are used, and the power spectrum is the ratio P(k)/s8^2, then don't do the derivative
of s8(z) w.r.t. the cosmological parameters*)

logFirstDeriv=OptionValue[logFirstDerivativeOfP];

If[UnsameQ[custvecCosmo,None]==True,
$vectorCosmo=custvecCosmo,
    Which[method=="FullNumerical",
       $vectorCosmo[zref_,kref_, muref_]:=((dObsPowerDpi[zref,kref,muref,#,dlnPdpDerivative->logFirstDeriv])&/@Range[$numZindependentPars]),
          method=="Jacobian",
       $vectorCosmo[zref_,kref_]:=((dPowerDpi[zref,kref,#,dlnPdpDerivative->logFirstDeriv])&/@Range[$numZindependentPars]);
    ]
];
If[UnsameQ[custvecZf,None]==True,
$vectorZFunctions=custvecZf,
If[Length@$zdependDerivVector>(Length@$zDepDerivativesListFull), (*{d1,d2,d3,d4,d5,d6,d7,d8,d9}*)
    Message[FisherBlockBuilder::incop, Length@$zdependDerivVector];
    myInterrupt[];
,
If[method=="FullNumerical" || method=="Jacobian",
    $vectorZFunctions[zref_,kref_, muref_]:=#[zref,kref,muref]&/@$zdependDerivVector;
]
]
];

If[UnsameQ[custvecChS,None]==True,
  $vectorFullChainSum=custvecChS,
  If[method=="ChainRule",
  (*deltaijVect=Array[deltaij, Length[$zdependDerivVector]];*)
  (* This deltaij vector removes the cross-terms between elements of the zdependDerivVector.
  Only left here for future tests *)
  $vectorFullChainSum[zr_,kr_,mur_]:=
   (
    (
      (
       If[OptionValue[onlyZdependentParameters]==True,
         setFisherDimensions[0,$numZdependentPars, $numZbins];
         0
         ,
         (dObsPowerDpi[zr,kr, mur, (#/.$parampositions),dlnPdpDerivative->logFirstDeriv])+0
             (*If[$fbsigmaBool==True,
           $s8zderivswitch*(-2)*dBackgroundDpi[sigma8ofZ,zr,#, dLogDerivative->True],
           0
           ]*)
         ]
         +
         If[MemberQ[$zdependDerivVector,d5]==False,
           d5[zr,kr,mur]*dBackgroundDpi[Hubble,zr,#, dLogDerivative->True], 0]
         +
         If[MemberQ[$zdependDerivVector,d4]==False,
           d4[zr,kr,mur]*dBackgroundDpi[distanceAngular,zr,#, dLogDerivative->True], 0]
         +
         If[$fbsigmaBool==False && MemberQ[$zdependDerivVector,d3]==False,
           d3[zr,kr,mur]*dBackgroundDpi[fGrowthRate,zr,#, dLogDerivative->True],  0]
         +
         If[$fbsigmaBool==True && MemberQ[$zdependDerivVector,d7]==False,
           debugPrint[$zdependDerivVector];
           d7[zr,kr,mur]*(dBackgroundDpi[fGrowthRateSigma8ofZ,zr,#, dLogDerivative->True]), 0]
        )&/@($paramnames)
      )~Join~(
        (Sequence@@Thread[Times[Table[$zdependDerivVector[[ii]]@@{zr,kr,mur}, {ii,Length[$zdependDerivVector]}],
                          KroneckerDelta[zr,#]]]
         )&/@(zaverage@$zbinGlobal))
     ) (*This works also if $zdependDerivVector={} *)
   ]
];
Which[method=="Jacobian",
ClearAll[$fisherCosmoParsBlock];
ClearAll[$fisherCrossTermBlock];
ClearAll[$fisherDiagonalZTermsBlock];
If[$zdependDerivVector=={}, Message[FisherBlockBuilder::wrmethj]; myInterrupt[]];
$fisherCosmoParsBlock[zr_,kr_]:=$fisherCosmoParsBlock[zr,kr]=KroneckerProduct[$vectorCosmo[zr,kr],$vectorCosmo[zr,kr]];
$fisherCrossTermBlock[zr_,kr_,mur_]:=$fisherCrossTermBlock[zr,kr,mur]=KroneckerProduct[$vectorCosmo[zr,kr],$vectorZFunctions[zr,kr,mur]];
$fisherDiagonalZTermsBlock[zr_,kr_,mur_]:=$fisherDiagonalZTermsBlock[zr,kr,mur]=KroneckerProduct[$vectorZFunctions[zr,kr,mur],$vectorZFunctions[zr,kr,mur]];
printConditionalInformation["\nFisher method: ", method, " chosen"];
printConditionalInformation[Style["\nPass the variables: $fisherCosmoParsBlock, $fisherDiagonalZTermsBlock and $fisherCrossTermBlock to FisherMatrixGCCalculation", 14, Blue]];
,
method=="FullNumerical",
ClearAll[$fisherCosmoParsBlock];
ClearAll[$fisherCrossTermBlock];
ClearAll[$fisherDiagonalZTermsBlock];
$fisherCosmoParsBlock[zr_,kr_,mur_]:=$fisherCosmoParsBlock[zr,kr,mur]=KroneckerProduct[$vectorCosmo[zr,kr,mur],$vectorCosmo[zr,kr,mur]];
If[UnsameQ[$zdependDerivVector,{}],
$fisherCrossTermBlock[zr_,kr_,mur_]:=$fisherCrossTermBlock[zr,kr,mur]=KroneckerProduct[$vectorCosmo[zr,kr,mur],$vectorZFunctions[zr,kr,mur]];
$fisherDiagonalZTermsBlock[zr_,kr_,mur_]:=$fisherDiagonalZTermsBlock[zr,kr,mur]=KroneckerProduct[$vectorZFunctions[zr,kr,mur],$vectorZFunctions[zr,kr,mur]];
];
printConditionalInformation["\nFisher method: ", method, " chosen"];
If[UnsameQ[$zdependDerivVector,{}],
printConditionalInformation[Style["\nPass the variables: $fisherCosmoParsBlock, $fisherDiagonalZTermsBlock and $fisherCrossTermBlock to FisherMatrixGCCalculation", 14, Blue]],
printConditionalInformation[Style["\nPass the variables: $fisherCosmoParsBlock, None and None to FisherMatrixGCCalculation", 14, Blue]]
];
,
method=="ChainRule",
ClearAll[$fisherCosmoFullChainBlock];
$fisherCosmoFullChainBlock[zr_,kr_,mur_]:=$fisherCosmoFullChainBlock[zr,kr,mur]=
Evaluate[(KroneckerProduct[$vectorFullChainSum[zr,kr,mur],$vectorFullChainSum[zr,kr,mur]](*//.{(deltaij[1]deltaij[2])->0, deltaij[1]->1, deltaij[2]->1}*))];
printConditionalInformation["\nFisher method: ", method, " chosen"];
printConditionalInformation[Style["\nPass the variables: $fisherCosmoFullChainBlock, None and None to FisherMatrixGCCalculation", 14, Blue]]
]
];


Off[NIntegrate::izero];
Options[FisherIntegration]={kminIntLimit:>$kminRef,zBins:>$zbinGlobal}~Join~{exponentVolumeEffective->2,functionNIntegrate->NIntegrate};

FisherIntegration[zref_,fisherBlockab_,opts:OptionsPattern[{FisherIntegration,NIntegrate}]]:=With[{kmaxx=kmaxChoice[zref],
kminn=OptionValue[kminIntLimit],
zbinglob=OptionValue[zBins],kk=Global`k,muu=Global`mu,
intopts=FilterRules[{opts},Options[NIntegrate]],
expopts=OptionValue[exponentVolumeEffective],
accgoal=OptionValue[AccuracyGoal],
precgoal=OptionValue[PrecisionGoal],
funcNInt=OptionValue[functionNIntegrate]},
  debugPrint["accgoal: "<>ToString[accgoal], 1];
  debugPrint["precgoal: "<>ToString[precgoal], 1];

  funcNInt[((fisherBlockab))*volumeEffective[zref,kk,muu,exponentVolumeEffective->expopts]*kk^2*dampingTerm[zref,kk,muu],{kk,kminn,kmaxx},{muu,0.,1.},
    AccuracyGoal->accgoal, PrecisionGoal->precgoal]*volumeSurvey[zref,zbinglob]*(2/(8 Pi^2))
];



ClearAll[FisherIntTest]


Options[FisherIntTest]={cosmoBlockPkderivatives->False};


FisherIntTest[fishBlock_, zm_, opts:OptionsPattern[]]:=Block[{block, pri, cosmokPdep, resu, len, zkk},

cosmokPdep=OptionValue[cosmoBlockPkderivatives];
len=Length[$paramnames];
If[cosmokPdep==True,
  $kamu=Global`k;
  Print["k: $kamu=  ", $kamu];,

  $kamu=Sequence[Global`k, Global`mu];
  Print["k and mu: $kamu=  ", $kamu];
];
resu=ParallelTable[FisherMatrixGC[zm,a,b,fishBlock[zm,$kamu][[a,b]]], {a,1,2},{b,4,5}];
Print[resu]
]


Options[FisherMatrixGC]={symmetricMatrix->True}
FisherMatrixGC[zref_?NumericQ,alpha_,beta_, fisherBlock_,opts:OptionsPattern[]]:=FisherMatrixGC[zref,alpha,beta,fisherBlock,opts]=
If[OptionValue[symmetricMatrix]==True&&alpha<beta,
FisherMatrixGC[zref,beta,alpha,fisherBlock,opts],
FisherIntegration[zref,fisherBlock(*[[alpha,beta]]*)]
];


Options[FisherTensor]={tensorOrder->3};

FisherTensor[zref_,alpha_,beta_,gamma_,
delta_ ,fisherBlock_,opts:OptionsPattern[]]:=FisherTensor[zref,alpha,beta,gamma,delta,fisherBlock,opts]=
Switch[OptionValue[tensorOrder],
3,
If[beta<gamma,
FisherTensor[zref,alpha,gamma,beta,delta,fisherBlock,opts],
FisherIntegration[zref,fisherBlock(*[[alpha,beta,gamma]]*)]
],
4,
If[alpha<beta,
FisherTensor[zref,beta,alpha,gamma,delta,fisherBlock,opts],
If[gamma<delta,
FisherTensor[zref,beta,alpha,delta,gamma,fisherBlock,opts],
FisherIntegration[zref,fisherBlock(*[[alpha,beta,gamma,delta]]*)]
]]]


FisherTensorS::usage="FisherTensorS[zref_,alpha_,beta_,gamma_,fisherBlockab_,opts]. Function that calculates the Fisher order three tensor formed by P,alpha * P,{beta gamma}"
FisherTensorJ::usage="FisherTensorJ[zref_,a_,b_,c_,fisherBlockab_,opts:]. Function that calculates the Fisher order three tensor formed by P,a * P,b * P,c"
FisherTensorQ::usage="FisherTensorQ[zref_,a_,b_,c_,d_,fisherBlockab_,opts]. Function that calculates the Fisher order four tensor formed by P,ab * P,bc"
prefactor::usage="Option for FisherTensor(S,J,Q). Specifies a prefactor that is multiplied in front of the integral."


Options[FisherTensorS]={tensorOrder->3, prefactor->1, exponentVolumeEffective->2};
FisherTensorS[zref_,alpha_,beta_,gamma_,fisherBlockab_,opts:OptionsPattern[]]:=
FisherTensorS[zref,alpha,beta,gamma,fisherBlockab,opts]=
FisherTensorS[zref,alpha,gamma,beta,fisherBlockab,opts]=With[{opti=OptionValue[tensorOrder], pref=OptionValue[prefactor]},
pref*FisherIntegration[zref,fisherBlockab]
]

Options[FisherTensorJ]={tensorOrder->3, prefactor->1, exponentVolumeEffective->3};
FisherTensorJ[zref_,a_,b_,c_,fisherBlockab_,opts:OptionsPattern[]]:=
FisherTensorJ[zref,b,a,c,fisherBlockab,opts]=
FisherTensorJ[zref,b,c,a,fisherBlockab,opts]=
FisherTensorJ[zref,a,c,b,fisherBlockab,opts]=
FisherTensorJ[zref,a,b,c,fisherBlockab,opts]=
FisherTensorJ[zref,c,b,a,fisherBlockab,opts]=
FisherTensorJ[zref,c,a,b,fisherBlockab,opts]=
With[{expi=OptionValue[exponentVolumeEffective], pref=OptionValue[prefactor]},
pref*FisherIntegration[zref,fisherBlockab, exponentVolumeEffective->expi]
]


Options[FisherTensorQ]={tensorOrder->4, prefactor->1, exponentVolumeEffective->2};
FisherTensorQ[zref_,a_,b_,c_,d_,fisherBlockab_,opts:OptionsPattern[]]:=
FisherTensorQ[zref,a,b,c,d,fisherBlockab,opts]=
FisherTensorQ[zref,a,b,d,c,fisherBlockab,opts]=
FisherTensorQ[zref,c,d,a,b,fisherBlockab,opts]=
FisherTensorQ[zref,c,d,b,a,fisherBlockab,opts]=
FisherTensorQ[zref,d,c,a,b,fisherBlockab,opts]=
FisherTensorQ[zref,d,c,b,a,fisherBlockab,opts]=
FisherTensorQ[zref,b,a,d,c,fisherBlockab,opts]=
FisherTensorQ[zref,b,a,c,d,fisherBlockab,opts]=
With[{opti=OptionValue[tensorOrder], pref=OptionValue[prefactor]},
pref*FisherIntegration[zref,fisherBlockab]
]





FisherSum[zref_,alpha_,beta_,fisherBlock_]:=Block[{kmaxx=kmaxChoice[zref],kstep,mustep=0.1},
kstep=(Log[kmaxx]-Log[0.00785398])/30;
Total[Flatten[Table[
(mustep*kstep/2)*((((fisherBlock[[alpha,beta]])/.k->Exp[lk])*volumeEffective[zref,Exp[lk],mu]*Exp[3*lk]*dampingTerm[zref,Exp[lk],mu])+(((fisherBlock[[alpha,beta]])/.k->Exp[lk+kstep])*
volumeEffective[zref,Exp[lk+kstep],mu+mustep]*Exp[3*(lk+kstep)]*dampingTerm[zref,Exp[lk+kstep],mu+mustep])),{lk,Log[0.00785398],Log[kmaxx],kstep},{mu,0,1,mustep}]]]*volumeSurvey[zref,$zbinGlobal]*(2Pi/(2Pi)^3)]


Options[FisherMatrixGCCalculation]=$paramoptions~Join~{
  zDependentIntermediateBlocks->True,
  tensorFisherBlocks->False,
  cosmoBlockPkmuDependence->True,
  fisherBlockDerivativeMethod:>$fisherMethodCalcBlock,
  resultsDirectory:>$resultsDirectoryGC,
  storeTemporaryResults->False,
  compressTemporaryResults->True,
  fisherMatrixExportNameID->"fisher-version-1",
  clearFisherIntegrationCache->False,
  checkSurveySpecsConsistency->True,
  parallelEvaluation->True,
  fisherOrder3TensorBlock->None,
  fisherOrder4TensorBlock->None,
  fisherOrder111TensorBlock->None,
   createProgressDialogMonitor->True}


FisherMatrixGCCalculation::nospecs="Survey specifications are not correctly set. Check: Options[setGalaxySurveySpecs]"
FisherMatrixGCCalculation::warnmeth="Warning: Method chosen for the Fisher Matrix calculation may not be correct or consistent with FisherBlockBuilder function."


FisherMatrixGCCalculation[fishCosmoBlock_,fishZBlockDiagonal_,fishZBlockCross_,opts:OptionsPattern[]]:=Block[{
totaltime, steptime=0, totalGals=0,
time1,time2,time3,time4,time5, time6, tensorOpt, zdepOpt, zmid, kmaxbin,
tempDirOpt, resultdir, tempdir, fisherBlockMethod,
cosmoblock, cosmoblockdims, crossZblock, diagZblock, Stensblock,Qtensblock, fish3Tensor, fish4Tensor, fish111Tensor,
zdepCrossBool,tensorBool,tempdirBool, comprBool, cosmokPdep, zk, fishID, optBlockMethod},
totaltime={0};
Global`blocknum=0;
Global`prog1=Global`prog2=Global`prog3=Global`prog4=Global`prog5=Global`prog6=0;
time1=time2=time3=time4=time5=time6=0;
SetSharedVariable[Global`prog1,Global`prog2,Global`prog3,Global`prog4,Global`prog5, Global`prog6];

optBlockMethod=OptionValue[fisherBlockDerivativeMethod];
If[optBlockMethod != $fisherMethodCalcBlock || StringQ[$fisherMethodCalcBlock]==False, Message[FisherMatrixGCCalculation::warnmeth]; myInterrupt[]];
(*according to the fisherBlockMethod other options can be automatically set, without having to ask for them*)
Switch[optBlockMethod,
"FullNumerical",
fisherBlockMethod=1;
If[($zdependDerivVector=={} || $numZdependentPars==0), zdepCrossBool=False, zdepCrossBool=True];
Print["FullNumerical calculation method, z-dependent parameter functions: ", zdepCrossBool];
cosmokPdep=True;,
"Jacobian",
fisherBlockMethod=2;
zdepCrossBool=True;
Print["Jacobian calculation method"];
cosmokPdep=False;,
"ChainRule",
fisherBlockMethod=3;
zdepCrossBool=False;
Print["ChainRule calculation method"];
cosmokPdep=True;,
"Debug",
fisherBlockMethod=1;
zdepCrossBool=OptionValue[zDependentIntermediateBlocks];
cosmokPdep=OptionValue[cosmoBlockPkmuDependence];,
_,
Print["fisherBlockDerivativeMethod set to a non standard value, proceed with default: FullNumerical"];
fisherBlockMethod=1;
zdepCrossBool=True;
cosmokPdep=True;
];

If[cosmokPdep==True,
  Print["Shape parameters block is k and mu dependent"];
  $kamu=Sequence[Global`k, Global`mu],

  Print["Shape parameters block is only k dependent"];
  $kamu=Global`k
];

If[OptionValue[checkSurveySpecsConsistency]==True,
 If[TrueQ[Length@$galaxiesCountInBins==Length@zaverage@$zbinGlobal==Length[$ndensGalZFuncGlobal[#]&/@zaverage@$zbinGlobal]],
 Print["Survey specifications have consistent dimensions compared to redshift bins"],
 Message[FisherMatrixGCCalculation::nospecs]];
 If[VectorQ[$biasInterpFunc[#]&/@zaverage@$zbinGlobal, NumericQ] && VectorQ[$ndensGalZFuncGlobal[#]&/@zaverage@$zbinGlobal, NumericQ],
 Print["Survey specifications n(z) and b(z) produce numeric output"],
 Message[FisherMatrixGCCalculation::nospecs]]
]

Print["Start Fisher Matrix calculation in redshift bins"];

resultdir=OptionValue[resultsDirectory];
$resultsDirectoryGC=resultdir;

tempdirBool=OptionValue[storeTemporaryResults];
comprBool = OptionValue[compressTemporaryResults];

$fisherMatResultID=OptionValue[fisherMatrixExportNameID];

tensorBool=OptionValue[tensorFisherBlocks];
If[tensorBool==True,
fish3Tensor=OptionValue[fisherOrder3TensorBlock];
fish4Tensor=OptionValue[fisherOrder4TensorBlock];
fish111Tensor=OptionValue[fisherOrder111TensorBlock];
];

If[tempdirBool==True,
tempdir=mkDirectory[resultdir<>"temporaryfiles/"];
tempdir=mkDirectory[resultdir<>"temporaryfiles/"<>$fisherMatResultID<>"/"];
Print["Temporary directory created at:  "<>ToString[tempdir]];
];

If[OptionValue[clearFisherIntegrationCache]==True,
ParallelTable[removeDownValues[FisherMatrixGC[_,_Integer,_Integer,__]], {i, $KernelCount}];
ParallelTable[removeDownValues[FisherTensor[_,_Integer,_Integer, _Integer, _Integer,__]],{i, $KernelCount}];
ClearAll[$cosmoblock];
ClearAll[$crossZblock];
ClearAll[$diagZblock];
ClearAll[$Stensorblock];
ClearAll[$Qtensorblock];
ClearAll[$Jtensorblock];
];

If[fisherBlockMethod==3,
cosmoblockdims = $numZindependentPars+$numZbins*$numZdependentPars;,
cosmoblockdims = $numZindependentPars];

If[OptionValue[createProgressDialogMonitor]==True,
nb=CreateProgressDialog[
{Dynamic[Global`blocknum],$numZbins,"Blocks"},
{Dynamic[Global`prog1],cosmoblockdims^2,"CosmoPars"},
If[zdepCrossBool==True, {Dynamic[Global`prog2],$numZindependentPars*$numZdependentPars,"CrossTerms"}, Sequence@@{}],
If[zdepCrossBool==True, {Dynamic[Global`prog3],$numZdependentPars^2,  "ZdiagonalTerms"}, Sequence@@{}],
If[(tensorBool==True && UnsameQ[fish3Tensor,None]), {Dynamic[Global`prog4],$numZindependentPars^3,  "3-Tensor"}, Sequence@@{}],
If[(tensorBool==True && UnsameQ[fish4Tensor,None]), {Dynamic[Global`prog5],$numZindependentPars^4,  "4-Tensor"}, Sequence@@{}],
If[(tensorBool==True && UnsameQ[fish111Tensor,None]), {Dynamic[Global`prog6],$numZindependentPars^3,  "111-Tensor"}, Sequence@@{}]
];
];

Do[
zmid=centralZbin[indz,$zbinGlobal];
kmaxbin=kmaxChoice[zmid];
Print["Central redshift in bin: ", zmid, "  Maximum k in bin: ",kmaxbin];
totalGals=Total[$galaxiesCountInBins[[1;;indz]]];
Print["Volume in bin: ", volumeSurvey[zmid,$zbinGlobal], "  number density in bin: ", $ndensGalZFuncGlobal[zmid],
"  Number of galaxies: ", $galaxiesCountInBins[[indz]], "  Total galaxies: ", totalGals];

Global`prog1=Global`prog2=Global`prog3=Global`prog4=Global`prog5=Global`prog6=0;

time1=First@AbsoluteTiming[$cosmoblock[indz]=ParallelTable[Global`prog1++;FisherMatrixGC[zmid,a,b,fishCosmoBlock[zmid,$kamu][[a,b]]],{a,cosmoblockdims},{b,cosmoblockdims}];];
If[tempdirBool, Export[tempdir<>"cosmoblock-"<>ToString[indz]<>"-z-"<>ToString[zmid]<>If[comprBool,".mc",".txt"],If[comprBool, Compress@$cosmoblock[indz], $cosmoblock[indz] ],If[comprBool, "String","Table"]]];

If[zdepCrossBool==True && UnsameQ[fishZBlockCross,None],
  time2=First@AbsoluteTiming[$crossZblock[indz]=ParallelTable[Global`prog2++;FisherMatrixGC[zmid,a,b, fishZBlockCross[zmid,Global`k,Global`mu][[a,b]],symmetricMatrix->False],{a,$numZindependentPars},{b,$numZdependentPars}];];
  If[tempdirBool, Export[tempdir<>"crossZblock-"<>ToString[indz]<>"-z-"<>ToString[zmid]<>".mc",Compress[$crossZblock[indz]],"String"]]
];

If[zdepCrossBool==True && UnsameQ[fishZBlockDiagonal,None],
  time3=First@AbsoluteTiming[$diagZblock[indz]=ParallelTable[Global`prog3++;FisherMatrixGC[zmid,a,b,fishZBlockDiagonal[zmid,Global`k,Global`mu][[a,b]]],{a,$numZdependentPars},{b,$numZdependentPars}];];
  If[tempdirBool, Export[tempdir<>"diagonalZblock-"<>ToString[indz]<>"-z-"<>ToString[zmid]<>".mc",Compress[$diagZblock[indz]],"String"]]
];
If[tensorBool==True && UnsameQ[fish3Tensor,None],
  time4 = First@AbsoluteTiming[$Stensorblock[indz]=ParallelTable[Global`prog4++;FisherTensorS[zmid,a,b,c,fish3Tensor[zmid,Global`k, Global`mu][[a,b,c]],tensorOrder->3],{a,$numZindependentPars},{b,$numZindependentPars},{c,$numZindependentPars}];];
  If[tempdirBool, Export[tempdir<>"S-tensor-"<>ToString[indz]<>"-z-"<>ToString[zmid]<>".mc",Compress[$Stensorblock[indz]],"String"]]
];
If[tensorBool==True && UnsameQ[fish4Tensor,None],
  time5 = First@AbsoluteTiming[$Qtensorblock[indz]=ParallelTable[Global`prog5++;FisherTensorQ[zmid,a,b,c,d,fish4Tensor[zmid,Global`k, Global`mu][[a,b,c,d]],tensorOrder->4],{a,$numZindependentPars},{b,$numZindependentPars},{c,$numZindependentPars},{d,$numZindependentPars}];];
  If[tempdirBool, Export[tempdir<>"Q-tensor-"<>ToString[indz]<>"-z-"<>ToString[zmid]<>".mc",Compress[$Qtensorblock[indz]],"String"]]
];
If[tensorBool==True && UnsameQ[fish111Tensor,None],
  time6 = First@AbsoluteTiming[$Jtensorblock[indz]=ParallelTable[Global`prog6++;FisherTensorJ[zmid,a,b,c,fish111Tensor[zmid,Global`k, Global`mu][[a,b,c]]],{a,$numZindependentPars},{b,$numZindependentPars},{c,$numZindependentPars}];];
  If[tempdirBool, Export[tempdir<>"J-tensor-"<>ToString[indz]<>"-z-"<>ToString[zmid]<>".mc",Compress[$Jtensorblock[indz]],"String"]]
];

steptime=time1+time2+time3+time4+time5+time6;
AppendTo[totaltime,steptime];

Print["index:", indz, " -- calculation time: ", " cosmoterms 1: ",time1," crossterms 2: ", time2," zdepterms 3: ", time3,
 If[tensorBool==True, Sequence@@{" Stensor 4: ", time4," Qtensor 5: ", time5," Jtensor 6: ", time6}, "  "], ", time in this step: ", steptime ",  cumulative calculation time: ", Total[totaltime]];

Global`blocknum++;
Print["****----"],
{indz,1,$numZbins}];

Print["Calculation finished, total computing time: ", Total[totaltime]];
Return[totaltime]

]


Options[InitializeFisherMatrixGCCalculation]=$paramoptions~Join~{zDependentVariablesVector->{d2}}~Join~{tensorFisherBlocks->False,
  fisherBlockDerivativeMethod:>$fisherMethodCalcBlock,customVectorCosmoPar->None,customVectorZFunction->None,customVectorFullChainSum->None,
  logFirstDerivativeOfP->True,fbsigma8Variables->False,onlyZdependentParameters->False,printInformation->True}~Join~{zDependentIntermediateBlocks->True,
  tensorFisherBlocks->False,cosmoBlockPkmuDependence->True,fisherBlockDerivativeMethod:>$fisherMethodCalcBlock,
  resultsDirectory:>$resultsDirectoryGC,storeTemporaryResults->False,compressTemporaryResults->True,
  fisherMatrixExportNameID->"fisher-version-1",clearFisherIntegrationCache->False,checkSurveySpecsConsistency->True,parallelEvaluation->True,
  fisherOrder3TensorBlock->None,fisherOrder4TensorBlock->None,fisherOrder111TensorBlock->None,createProgressDialogMonitor->True};

InitializeFisherMatrixGCCalculation[arg___,opts:OptionsPattern[]]:=Module[{zdepvect,method,blocks,fbs8varBool, fisherMatID, totaltime},
zdepvect=OptionValue[zDependentVariablesVector];
setFisherDerivativesVector[zdepvect];
$fisherMethodCalcBlock=OptionValue[fisherBlockDerivativeMethod];
fbs8varBool=OptionValue[fbsigma8Variables];
FisherBlockBuilder[zdepvect, fisherBlockDerivativeMethod->$fisherMethodCalcBlock,fbsigma8Variables->fbs8varBool];
Which[$fisherMethodCalcBlock=="FullNumerical",
blocks={$fisherCosmoParsBlock, $fisherDiagonalZTermsBlock, $fisherCrossTermBlock};
,
$fisherMethodCalcBlock=="ChainRule",
blocks={$fisherCosmoFullChainBlock, None, None};
];
$resultsDirectoryGC=OptionValue[resultsDirectory];
fisherMatID=OptionValue[fisherMatrixExportNameID];
totaltime=FisherMatrixGCCalculation[blocks[[1]], blocks[[2]], blocks[[3]], resultsDirectory:>$resultsDirectoryGC,
  clearFisherIntegrationCache->True, fisherMatrixExportNameID->fisherMatID, storeTemporaryResults->False];
Return[totaltime];
];

Options[FisherFinalAssembly]=$paramoptions~Join~{tensorFisherBlocks->False,
fisherBlockDerivativeMethod->$fisherMethodCalcBlock, resultsDirectory:>$resultsDirectoryGC,
fisherMatrixExportNameID->$fisherMatResultID, exportFisherMatrix->True, compressExportMatrix->True, exportParametersUsed->True}


FisherFinalAssembly::warndirfile="Warning: Exporting of Fisher Matrix will fail.
Either the desired filename `1` is not a string or the results directory `2` is not an existing directory."


FisherFinalAssembly[opts:OptionsPattern[]]:=Module[{zdepRemove, zdepElems, method, finalFisher, cosmoFixedFisher,
resultdir, export, compressopt, fisherID},

zdepElems=Range@Length@$zdependDerivVector;
zdepRemove=Sort[parameterPositionInBins[zdepElems]];

method=OptionValue[fisherBlockDerivativeMethod];
resultdir=OptionValue[resultsDirectory];
export=OptionValue[exportFisherMatrix];
compressopt=OptionValue[compressExportMatrix];
fisherID=OptionValue[fisherMatrixExportNameID];

Which[
(method=="FullNumerical" || method=="Jacobian") && UnsameQ[$zdependDerivVector, {}],
finalFisher=Sum[placeBlocks[$cosmoblock[i],$crossZblock[i],$diagZblock[i],i,$numZbins*$numZdependentPars+$numZindependentPars],{i,$numZbins}];
cosmoFixedFisher=Sum[$cosmoblock[i],{i,$numZbins}];
,
(method=="FullNumerical") && $zdependDerivVector=={},
finalFisher=Sum[$cosmoblock[i],{i,$numZbins}];
cosmoFixedFisher=finalFisher;
,
method=="ChainRule",
finalFisher=Sum[$cosmoblock[i],{i,$numZbins}];
If[$numZindependentPars!=0,
  cosmoFixedFisher=fixElements[finalFisher,zdepRemove];
]
];

If[OptionValue[exportParametersUsed]==True,
  (*DumpSave[resultdir<>"parameters-used--"<>fisherID<>".mx", {$paramoptions, $paramnames, $paramfidus, $paramlabels}];*)
Export[resultdir<>"parameters-used--"<>fisherID<>".txt", {(ToString[#]&/@$paramnames), $paramfidus}, "Table"];
];

$fisherMatrix = finalFisher;
$fisherMatrixCosmoPars = cosmoFixedFisher;

Print["Final Fisher Matrix is symmetric:  ", SymmetricMatrixQ[$fisherMatrix]];
Print["Final Fisher Matrix is postive definite:  ", PositiveDefiniteMatrixQ[$fisherMatrix]];

Print["Final Fisher matrix assembled: $fisherMatrix, Dimensions:  ", Dimensions[$fisherMatrix],
"\nFisher Matrix fund. cosmo. pars. only: $fisherMatrixCosmoPars, Dimensions:  ", Dimensions[$fisherMatrixCosmoPars] ];

Which[export==True && compressopt==True && StringQ[fisherID]==True && DirectoryQ[resultdir]==True,
Print["Compressed $fisherMatrix exported to: "];
Export[resultdir<>fisherID<>".mc",Compress[$fisherMatrix],"String"]
,
export==True && compressopt==False && StringQ[fisherID]==True && DirectoryQ[resultdir]==True,
Print["Text form $fisherMatrix exported to: "];
Export[resultdir<>fisherID<>".txt",$fisherMatrix,"Table"]
,
export==True && (StringQ[fisherID]==False || DirectoryQ[resultdir]==False),
Message[FisherFinalAssembly::warndirfile, fisherID, resultdir];
myInterrupt[];
,
export==False,
Print[Style["Final fisher matrix not exported.
Check options of this function, for exporting the result of the Fisher Matrix calculation: $fisherMatrix", 16, Blue]];
]
]



Options[fullFisherPostTransformation]={repairPositivity->False, fisherMatrixDimensions->{$numZindependentPars,$numZdependentPars,$numZbins},
marginalizedElementsIndex->{False},exportTableName->None, fixedElementsIndex->{False},
performReduceOperations->"None", zBinsJacobian->False, cosmoParamsJacobian->False}

fullFisherPostTransformation[fishermat_,opt:OptionsPattern[]]:=Block[{FisherMat,FisherMat2,FisherMat3, FisherMat4,
projectedZetFish, finalProjectedParsFish, dims, elemMargs,
nrems,margelement,margind, opers, fixelem, fixind, elemRemove,
cJac, zJac, fishname, expfile},
dims=OptionValue[fisherMatrixDimensions];
setFisherDimensions[dims/.{List->Sequence}];
Print["Positivity of original Matrix: ",PositiveDefiniteMatrixQ[fishermat]];
If[OptionValue[repairPositivity],
FisherMat=fixPositiveDefiniteMatrix[fishermat];
Print["Positivity of repaired Matrix: ",PositiveDefiniteMatrixQ[FisherMat]];,
FisherMat=fishermat];
opers=OptionValue[performReduceOperations];
margind = OptionValue[marginalizedElementsIndex];
fixind = OptionValue[fixedElementsIndex];
Switch[opers,
"Marginalize",
  Which[UnsameQ[Intersection[margind,$zDepIndicesListFull],{}],
      margelement=Flatten@(margind/.$zdependDerivPositions);
      elemRemove=Sort[parameterPositionInBins[margelement]];,
      VectorQ[margind,NumericQ],
      elemRemove = margind;
      ];
nrems=Length@Flatten[List[elemRemove]];
Print[nrems , " Elements to marginalize:  ", elemRemove];
FisherMat2=marginalizeElements[FisherMat,elemRemove];
Print["Positivity of marginalized Matrix: ",PositiveDefiniteMatrixQ[FisherMat2]];,
"Fix",
  Which[UnsameQ[Intersection[fixind,$zDepIndicesListFull],{}],
    fixelem=Flatten@(fixind/.$zdependDerivPositions);
    elemRemove=Sort[parameterPositionInBins[fixelem]];,
    VectorQ[fixind,NumericQ],
    elemRemove = fixind;
  ];
nrems=Length@Flatten[List[elemRemove]];
Print[nrems , " Elements to fix:  ", elemRemove];
FisherMat2=fixElements[FisherMat,elemRemove];
Print["Positivity of fixed Matrix: ",PositiveDefiniteMatrixQ[FisherMat2]];,
"None",
FisherMat2=FisherMat;
Print["No reduction of dimensions performed"]
];
(*setFisherDimensions[dims[[1]],dims[[2]]-nrems,dims[[3]]];*)
zJac=OptionValue[zBinsJacobian];
cJac=OptionValue[cosmoParamsJacobian];
If[MatrixQ[zJac]==False,
FisherMat3=FisherMat2,
FisherMat3=jacobianTransform[FisherMat2,zJac];
Print["Positivity of z-projected Matrix: ",PositiveDefiniteMatrixQ[FisherMat3]];
];
If[MatrixQ[cJac]==False,
FisherMat4=FisherMat3,
FisherMat4=jacobianTransform[FisherMat3,cJac];
Print["Positivity of cosmo params transformed Matrix: ",PositiveDefiniteMatrixQ[FisherMat4]];
];
fishname=OptionValue[exportTableName];
If[UnsameQ[fishname, None] && StringQ[fishname],
  expfile=Export[fishname<>".txt",FisherMat4,"Table"];
  Print[expfile],
  Print["Transformed Fisher Matrix not exported"]
];
Return[FisherMat4]
]

clearFisherMemoizedQuantities[opts:OptionsPattern[]]:=Module[{var},

  ParallelTable[removeDownValues[FisherMatrixGC[_,Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],___]], {i, $KernelCount}];
  removeDownValues[FisherMatrixGC[_,Except[_PatternTest],Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[FisherTensor[_,Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[FisherTensor[_,Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[$fisherCosmoParsBlock[Except[_Pattern | _PatternTest],__]],{i, $KernelCount}];
  removeDownValues[$fisherCosmoParsBlock[Except[_Pattern | _PatternTest],__]];
  ParallelTable[removeDownValues[$fisherDiagonalZTermsBlock[Except[_Pattern | _PatternTest],_,_]],{i, $KernelCount}];
  removeDownValues[$fisherDiagonalZTermsBlock[Except[_Pattern | _PatternTest],_,_]];
  ParallelTable[removeDownValues[$fisherCrossTermBlock[Except[_Pattern | _PatternTest],_,_]],{i, $KernelCount}];
  removeDownValues[$fisherCrossTermBlock[Except[_Pattern | _PatternTest],_,_]];
  ParallelTable[removeDownValues[$fisherCosmoFullChainBlock[Except[_Pattern | _PatternTest],_,_]],{i, $KernelCount}];
  removeDownValues[$fisherCosmoFullChainBlock[Except[_Pattern | _PatternTest],_,_]];
  ParallelTable[removeDownValues[dObsPowerDZpi[Except[_Pattern | _PatternTest],_,_,Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[dObsPowerDZpi[Except[_Pattern | _PatternTest],_,_,Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[dObsPowerDpi[Except[_Pattern | _PatternTest],_,_,Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[dObsPowerDpi[Except[_Pattern | _PatternTest],_,_,Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[dPowerDpi[Except[_Pattern | _PatternTest],_,Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[dPowerDpi[Except[_Pattern | _PatternTest],_,Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[d2PowerD2pii[Except[_Pattern | _PatternTest],_,Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[d2PowerD2pii[Except[_Pattern | _PatternTest],_,Except[_Pattern | _PatternTest],___]];
  ParallelTable[removeDownValues[d2PowerD2pij[Except[_Pattern | _PatternTest],_,Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[d2PowerD2pij[Except[_Pattern | _PatternTest],_,Except[_Pattern | _PatternTest],Except[_Pattern | _PatternTest],___]];
   ParallelTable[removeDownValues[d2ObsPowerD2pii[Except[_Pattern | _PatternTest],_,_,Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[d2ObsPowerD2pii[Except[_Pattern | _PatternTest],_,_,Except[_Pattern | _PatternTest],___]];
   ParallelTable[removeDownValues[d2ObsPowerD2pij[Except[_Pattern | _PatternTest],_,_,Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest],___]],{i, $KernelCount}];
  removeDownValues[d2ObsPowerD2pij[Except[_Pattern | _PatternTest],_,_,Except[_Pattern | _PatternTest], Except[_Pattern | _PatternTest],___]];
   ParallelTable[removeDownValues[dBackgroundDpi[_,Except[_Pattern | _PatternTest],_,___]],{i, $KernelCount}];
  removeDownValues[dBackgroundDpi[_,Except[_Pattern | _PatternTest],_,___]];

]



Protect[d1,d2,d3,d4,d5,d6,d7,d8,d9];
Protect[nbins];
End[]
EndPackage[]
