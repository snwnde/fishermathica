
# Fishermathica


A Mathematica code to perform Fisher Matrix forecasts for photometric and spectroscopic galaxy surveys.

This is a fresh commit github to avoid carrying around a lot of commit history and old tracked files. The full commit history (570 commits) can be found on 
https://gitlab.com/santiagocasas/FisherMatrixTools.

## Release 0.1
Proper documentation and end-to-end examples will be added in the next full release.


## Packages documentation

### Package: CosmologyFunctions

`APeffect :: Option for observed power spectrum. Activates or deactivates the Alcock-Pazcynski effect, Default: True


complementParamValues ::  complementParamValues[extrapars_,functionsym_, 
                         opts:globalPars->$paramoptions,returnList->Complement,filterCosmoPars->False] 
                         Finds the complement of the set of parameter options passed as extra arguments extrapars and 
                         the default options for the function functionsym. 
                         If returnList=Full, it returns the full list of parameters with the extrapars replacing the positions of the defaultpars.
                         If returnList=Default, it returns the default options for the function functionsym.
                         If filterCosmoPars->True it filters the final list, accepting only options defined in the local globalPars option 
                         (usually cosmological parameters defined by $paramoptions).
curvatureSin :: Gives the Sinus function depending on the curvature value and sign. For Abs[ok]<=0.001 it returns just the argument x.
                Form: curvatureSin[x,ok]= either (Sinh[x], Sin[x], x)

dampingTerm :: dampingTerm[z,k,mu,opts] gives the nonlinear damping term correction to the observed power spectrum. 
               Needs the setting of the variables $growthNorm, $sigma0Value and $dampSwitch

DimensionlessHubble :: Gives the value of the dimensionless Hubble parameter


distanceAngular :: gives the angular diameter distance as a function of z


dlnPdkDeriv ::  dlnPdkDeriv[zred_,kr_,opts], Gives the numerical Logarithm derivative of P(k) wrt k, as a function of k and z.
               Accepts all options that are accepted by powerSpectrum

errorZ :: redshift error as a function of scale k, angle mu and redshift z with optional cosmological parameters given by opts
          Form: errorZ[k,mu,z,opts]

fGrowthRate :: gives the growth rate f of matter perturbations for LCDM or a general cosmology. Accepts the option lcdmBool->True for LCDM


Growth :: gives the Growth of a general cosmology as a function of z, accepts the option lcdmBool->True for LCDM


Hubble :: Gives the value of the dimensionful Hubble parameter for other model different from LCDM if lcdmBool->False


kmuAlcockPaczynski :: kmuAlcockPaczynski[zr_,kr_,mur_,opts:] : 
                      Gives a list {k,mu} with the modified k and mu values due to the Alcock-Paczynski effect as a function of redshift 

LcdmCambPk ::  LcdmCambPk[zred, k, sigma8, opts], Gives power spectrum interpolation function from CAMB LCDM, Needs cosmomathica package. Options: $paramoptions and linearBool->True


LCDMCAMBPsPre :: Gives LogLog power spectrum from CAMB LCDM, Needs cosmomathica package. Options: $paramoptions and linearBool->True


LCDMGrowth :: gives the LCDM growth factor as a function of z defined through OmegaM^gamma. 
              Gamma can be set globally using setCosmologyFixedValues

LCDMHubble :: Gives the value of the dimensionful LCDM Hubble parameter at redshift z, accepts w0-w1 parametrization for DE


LCDMTrHuPre :: LCDMTrHuPre[sigma8val_,opts:], Gives power spectrum from Eisenstein&Hu Transfer functions, 
               using the Growth function of this package. opts: cosmological parameters and linearBool->True. If linearBool->False, 
               nonlinear PS is calculated with the function HalofitCorrection
normalizationFactor :: normalizationFactor[ps_,sigma8_] calculates sigma8^2/sigma8Function[ps,$radiusScale]
                       $radiusScale set by default to 8 

observedPowerSpectrum ::  observedPowerSpectrum[zred_,k_,mu_,opts:] : 
                         Gives the observed power spectrum as a function of redshift z, scale k and observing angle mu. 
                         It includes the raw P(k), bias, redshift space distortions, photometric redshift error, shot noise, 
                         and the Alcock-Paczynski effect on k,mu and the volume.
                         Options: $paramoptions, APeffect->True, and all other options available for the powerSpectrum function
OmegaBaryon0Today :: Outputs the value of OmegaB today. Valid options: Omegab, Omegac, omegac, omegab, omegam


OmegaM :: Gives the value of OmegaM in a general cosmology as a function of redshift, accepts the option lcdmBool->True for LCDM


OmegaM0Today :: Outputs the value of OmegaM today. Valid options: OmegaDE, Omegam, Omegac, omegac, Omegeab, omegab, omegam


OmegaMLCDM :: Gives the value of OmegaM as a function of redshift z in LCDM


PExtraShot :: Set some function of z to account for some extra shot noise term in the observed power spectrum


powerSpectrum :: powerSpectrum[z_,k_,opts]. Gives the power spectrum for a specified cosmology and a calculation method.
                 Options=$paramoptions, $pscosmoopts and the extra options:
                 spectrumMethod(Transfer or CAMB) calculation with Eisenstein&Hu Transfer functions or CAMB from cosmomathica
                 sigma8reference0.8 (reference sigma8 normalization for PS, just in case it is not specified as a cosmological parameter)
RAPfunction ::  RAPfunction[zred_,mur_,opts:] : Gives the R value of the Alcock Pazcynski effect 
               as a function of redshift and mu. Options: $paramoptions

rsdbetaFunc :: gives the redshift space distortion beta function as a function of z. rsdbeta=fGrowthRate/bias


setBiasFunction :: sets the bias interpolation function when given a list of redshifts and a list of bias values
                   Accepts options for Interpolation

setCosmologyFixedValues :: set global cosmological variables, the gamma parameter of the growth rate and H0 in units of h/Mpc


setKextremeValues :: setKextremeValues[kminr_,kmaxr_,kmaxmaxr_], 
                     kminr:Minimum k where interpolation functions of P are defined,
                     kmaxr:Maximum k value for integrals in sigma8 and other quantities
                     kmaxmaxr: Absolute maximum k value where interpolation functions should be defined
setParameterOptions :: set the default values of cosmological parameters by passing a list of fiducial values and a list of parameter names
                       setParameterOptions[paramNames_,paramFiducials_]

setPScosmoOptions :: set the Boolean options for specifying linear power spectrum and LCDM cosmology.
                     Pass a list of two parameters and a list of two booleans. First argument always lcdmBool, second argument linearBool

setTermsSwitches :: sets the options switches for different contributions, setTermsSwitches[betaswitch,dampSwitch],
                     betaSwitch shuts down cosmological information in R.S.Distortions, 
                    dampSwitch turns off the Seo&Eisenstein nonlinear damping term in the observed power spectrum
sigma8Function :: Calculates sigma8 for a power spectrum interpolation function at a radius scale scal
                  Form: sigma8Function[ps, scal, opts] Takes opts for NIntegrate

sigma8reference :: Option for power spectrum, when sigma8 is not specified as a parameter. Default: 0.8


sigmaR :: sets the redshift error exponent as a function of redshift


sigmaZ :: sets the redshift error of relative velocities as a percentage value


spectrumMethod :: Option for power spectrum, 'Transfer' or 'CAMB'. Default: Transfer


wDEzFunc :: calculate w(z) from a parametrization of w0 and w1 Form: wDEzFunc[z_number, w0_number, w1_number]=function


windowk :: Window function to integrate over in calculation of sigma8 for the power spectrum


$biasInterpFunc :: interpolation function for the fiducial bias


$H0 :: value of H0 in units of h/Mpc


$kMaxMax :: Maximum allowed k for interplations and computations in k. Used in dlnP/dlnk


$kmaxRef :: Maximum reference k value for integrals or other calculations in k. Used in sigma8Function


$kminRef :: Minimum k value where power spectrum is defined or to be used in k integrals


$modcosmoopt :: Boolean option specifying LCDM (True) or non-LCDM (False) cosmological functions


$paramfidus :: list of fiducial values for parameters


$paramnames :: list of parameter names for fiducials


$paramoptions :: list of rules between parameter names and parameter fiducial values


$partesta :: debug variable


$pscosmoopts :: Boolean options for cosmology and power spectrum, usually: lcdmBool->True and linearBool->True


$pslinearopt :: Boolean option  specifying linear (True) or nonlinear (False) power spectrum calculation `



### Package: FisherTools

`centralZbin :: Gives the central value of the bin, when the input is the bin index.
               Form: centralZbin[index_number,list_zbins] = number

complementOptionValues :: Form: complementOptionValues[customOptsList_,internalFunctOptsList_,externalFunctOptsList_, options: returnList->'CustomComplementOptions', filterCosmoPars->'All'].
                          For a list of custom options, a list of options of an internal function and a list of options for an external function, it returns the following, according to the string set in the option returnList:
                          CustomFullOptions: gives all the custom options and the default options accepted by the internal function, which are NOT options of the external function given.
                          CustomComplementOptions: gives only the custom options accepted by the internal function, which are NOT options of the external function given.
                          ExternalComplementOptions: gives the custom options belonging to the set of options of the external function, which are set to a different value than the defaults of the external function. 
                          ExternalUnsetOptions: gives the custom options belonging to the set of options of the external function, which were not previously defined as options for the internal function.
                          DefaultOptions: gives back the default unchanged options of the internal function.
                          Option filterCosmoPars->:
                          'All': returns the above options output including cosmo parameters.
                          'In':  only returns from the above options output, the ones which are cosmo parameters.
                          'Out': returns the above options output, filtering out all cosmo parameters.
CreateProgressDialog ::  Create a pop-up window with progress vars and maximum values max.
                        CreateProgressDialog[Dynamic[var1_],Dynamic[var2_], Dynamic[var3_], Dynamic[var4_],max1_,max2_,max3_,max4_]

d2ObsPowerD2pii :: Gives the second derivative of the observed power spectrum wrt the parameter alpha.
                   Options: setEpsilon->$epsilonstep, transformDerivativeFunction->None, dlnPdpDerivative->False, $paramoptions, $pscosmoopts, spectrumMethod->'Transfer', APeffect->True
                   Form:d2ObsPowerD2pii[zred_,k_,alpha_,opts:]
d2ObsPowerD2pij :: Gives the second mixed derivative of the observed power spectrum wrt the parameter alpha and beta.
                   Options: setEpsilon->$epsilonstep, transformDerivativeFunction->None, $paramoptions, $pscosmoopts, spectrumMethod->'Transfer', APeffect->True
                   Form: d2ObsPowerD2pij[zred_,k_,alpha_,beta_,deropts:]
d2PowerD2pii :: Gives the second derivative of the power spectrum wrt the parameter alpha.
                Options: setEpsilon->$epsilonstep, transformDerivativeFunction->None, dlnPdpDerivative->False, $paramoptions, $pscosmoopts, spectrumMethod->'Transfer'
                Form:d2PowerD2pii[zred_,k_,alpha_,opts:]
d2PowerD2pij :: Gives the second mixed derivative of the power spectrum wrt the parameter alpha and beta.
                Options: setEpsilon->$epsilonstep, transformDerivativeFunction->None, $paramoptions, $pscosmoopts, spectrumMethod->'Transfer'
                Form: d2PowerD2pij[zred_,k_,alpha_,beta_,deropts:]
dBackgroundDpi ::  Gives the derivative of the background function, bfunction at z=zred, for the parameter name param.
                  Options: setEpsilon->$epsilonstep, dLogDerivativeTrue, $paramoptions.
                  Form: dBackgroundDpi[bfunction_,zred_,param_,optpars]
dkdlogD :: derivative of k function wrt log(DistanceAngular)


dkdlogH :: derivative of k function wrt log(Hubble)


dlnPdpDerivative :: Option for derivative of the power spectrum


dmudlogD :: derivative of mu function wrt log(DistanceAngular)


dmudlogH :: derivative of mu function wrt log(Hubble)


dNdensDparInterpFunc :: dNdensDparInterpFunc[dn_,zbino_,alpha_,opts] 
                        returns the interpolation function of the derivative of the function ndens wrt the parameter alpha,
                         for a number count dn, redshift bins zbino and options which are the cosmological parameters and 
                        options for the function Interpolation
dNdensDparTable :: dNdensDparTable[dn_,zbino_,alpha_,opts] calculates the derivative of the function ndens wrt the parameter alpha,
                    for a number count dn, redshift bins zbino and options which are the cosmological parameters and setEpsilon

dObsPowerDpi :: Gives the first derivative of the observed power spectrum wrt the parameter alpha.
                Options: setEpsilon->$epsilonstep, transformDerivativeFunction->None, dlnPdpDerivative->False, $paramoptions, $pscosmoopts, spectrumMethod->'Transfer', APeffect->True
                Form:dObsPowerDpi[zred_,k_,alpha_,opts:]
dPowerDpi :: Gives the first derivative of the power spectrum wrt the parameter alpha.
             Options: setEpsilon->$epsilonstep, transformDerivativeFunction->None, dlnPdpDerivative->False, $paramoptions, $pscosmoopts, spectrumMethod->'Transfer'
             Form:dPowerDpi[zred_,k_,alpha_,opts:]
draw123Ellipses :: draws 1-,2- and 3-sigma confidence regions given by the ellipses of the function
                   drawEllipse using the scalings 1.51, 2.49 and 3.44 respectively. Colors for each respective line are passed as a list
                   Form: draw123Ellipses[matrix, a_index, b_index, fiducialParameter_list, colors_List] = {Graphic1,Graphic2,Graphic3}
drawEllipse :: Ellipse given by the Sqrt of the Eigenvalues of the Matrix and rotated 
               by the angle given by the function ellipseAngle for the parameters specified by indices a and b.
               scale parameter scales the area for the different confidence regions.
               Form: drawEllipse[matrix, a_index, b_index, fiducialParameter_list, scale] = Graphic
ellipseAngle :: Calculate inclination angle of ellipse w.r.t. the vertical axis, for the confidence region of a 2x2 fisher matrix F_2.
                Form: ellipseAngle[matrix]=angle in radians

errorsTable :: Generate a table of absolute and relative errors for each fiducial parameter of the Fisher Matrix. 
               The default row headers are: fiducial, abs.error, rel.error. This can be changed by the option RowHeaders.
                Form: errorsTable[vectorFiducials, vectorOneSigmaErrors, vectorParameterLabels, RowHeaders->vectorStringTitles] = TableForm
errorsTableBinsPars :: Generate a table of absolute and relative errors for redshift dependent parameters of the Fisher Matrix.
                       Accepts a list with the central values of the bins and packed in a resective Table, the fiducial lists and error lists. 
                       The last accepted list contains the headers labels of the  table.
                       Form: errorsTableBinsPars[binsList_,fidusTable_,binerrorsTable_,parlabList_,opts:]
extractCorner :: Extracts first (num x num) elements of the upper left corner of a matrix
                 Outputs warning if Matrix is not symmetric.
                 Form: extractCorner[matrix,num]=matrix
extractErrorAt :: Extracts fully marginalized errors at each redshift bin for z-dependent (default) or z-independent parameter  
                  corresponding to the index given. 
                  Form: extractErrorAt[errorlist_,index_,zDependent->True]=list
extraErrorNoise :: Option for ndensEffective tu shut down the extra relative noise on the power spectrum P(k).


FisherIntegration ::  Form: FisherIntegration[zref_,fisherBlock_, opts:]. Integrates in k and mu, a block matrix fisherBlock of Fisher integrands evaluated at the redshift zref.
                     Accepts all available options for NIntegrate with these as default: MaxRecursion->100,PrecisionGoal->7,AccuracyGoal->7,Method->{'GlobalAdaptive','MaxErrorIncreases'->15000}.
                     Other options: 
                     lower integration k limit: kminLimit->$kminRef
                     zbin limits for calculating survey volume: zbins->$zbinGlobal
FisherMatrixGC :: Defines the integration for the Galaxy Clsutering Fisher Matrix and takes care of the symmetry properties
                  Options: symmetric->True
                  Form: FisherMatrixGC[zref_?NumericQ,alpha_,beta_, fisherBlock_,opts:]
                  Only evaluates if zref is a numeric quantity, to prevent error messages.
fixElements :: fix parameter in Fisher Matrix, by removing row and column associated with the respective index. 
               Form: fixElements[matrix,list] = matrix

fixPositiveDefiniteMatrix :: Fixes the positivity of a matrix by reconstructing the 
                             matrix using its singular values as eigenvalues in the eigensystem. Only to be used for symmetric and real matrices.
                             Form: fixPositiveDefiniteMatrix[matrix]=matrix 
functionOutputImportCheck :: Checks if dataFileName exists and imports it as a table of values. 
                             If the file is not found, then the function calcFunc is evaluated with the argument argument_ and then it exports 
                             its results as a table to the file dataFileName.
                             Form: functionOutputImportCheck[calcFunc_,argument_,dataFileName_]=list (and exporting or importing of files)
insertColumn :: insert vector as a column at position pos. 
                Form: insertColumn[matrix,vector,pos] = matrix

insertRow :: insert vector as a row at position pos.
             Form: insertRow[matrix,vector,pos] = matrix

jacobianTransform :: Transform a Matrix by using a Jacobian Matrix: J .M.J 
                                                                      Form: jacobianTransform[matrix,jacobianmatrix] = matrix

kfuncSymbolic :: kfuncSymbolic[kref_,muref_,Dfunc_,Hfunc_,Dfuncref_,Hfuncref_] 
                 gives the symbolic k function relating k and k in the reference cosmology

kmaxChoice :: Function that chooses the maximum limit of integration in k for the Fisher Matrix, according to $kmaxMethod. 
              For a list of available options, type: ?$kmaxMethod.
              Options for kmaxChoice: kmaxHardValue (Default: $kmaxHard), kmaxChoiceMethod: (Default: $kmaxMethod), 
              kmaxInterpFunction: (Default: $kmaxInterpFunc)
kmaxZFuncTable :: Form: kmaxZFuncTable[sigma8FuncZK_,opts].
                  Calculates a table of maximum kmax values as a function of redshift for the fisher Matrix Integration.
                  It uses a constraint on the maximum allowed normalization of P(k), $maxKsigma8, 
                  to calculate the value kmax at each z, by finding the root of the function Sqrt[sigma8(k,z)]==Sqrt[$maxKsigma8].
                  Other options: rootIniValK=Pi/40.0(Starting point for finding root), zTableSteps=0.1(Steps in z to create final table), 
                  zbinValues=$zbinGlobal(zbins used to calculate the final table)
lineThickness :: Option for draw123Ellipses to set the absolute thickness of the lines in printer points.


marginalized2Matrix :: Marginalizes Fisher Matrix over every parameter except the ones specified by index a,b. 
                       Form: marginalized2Matrix[matrix, a_index, b_index] = 2x2 Matrix

marginalizeElements :: Marginalize parameters of the Fisher matrix associated with the list of indices given as input. 
                       Form: marginalizeElements[matrix,list] = matrix

mufuncSymbolic :: mufuncSymbolic[muref_,Dfunc_,Hfunc_,Dfuncref_,Hfuncref_]
                  gives the symbolic mu function that relates mu and mu in the reference cosmology 

nbins :: Gives back the number of bins in a list of binned values
         Form: nbins[list] = number

ndensEffective :: calculate the effective ndens used in the Fisher Matrix calculation 
                  by including or not extra noise affecting the power spectrum. 
                  This noise function has to be set manually using the global function relativePkNoise[zz,kk].
                  Options: Deactivate this extra noise by setting: extraErrorNoise->False (Default)
ndensLensFunc :: Gives the number density function of redshift for Weak Lensing forecasts.
                 Uses the previously defined $z0GalDist as a parameter

ndensVolumeDepTable :: Gives the number density of galaxies for a parameter dependent volume.
                       Form ndensVolumeCalcPars[dn,zbino,opts] dn:number count in bin,zbino:redshift bins, 
                       Options: dz=redshift bin width (default 0.1) and cosmological parameters given by $paramoptions
ndensVolumeInterpFunc :: ndensVolumeInterpFunc[dn_,zbino_,opts] calculates an interpolation function for the number density 
                         using ndensVolumeDepTable and the arguments dn: number count, zbino:redshift bins. Options are: 
                         The cosmological parameter options $paramoptions and all options available for Interpolation.
                         InterpolationOrder set to 1 by default.
nearestBinLimts :: Gives back the lower and upper limits of the bin that contains the point pz. 
                   If pz is smaller than smallest value of list_zbins, it uses 0 as smallest bin limit. 
                   If pz is bigger than biggest number in list_zbins, it returns a list that goes from the last element of list_zbins to pz.
                   Form: nearestBinLims[pz_, list_zbins] = list of two elements
numericalParamsDerivative :: 
                             numericalParamsDerivative[defaultopts_,extraopts_,indexlist_,epsilonValue_,opts]
                             computes the step, and the parameter options for numerical derivatives of Cosmology dependent functions.
                             defaultopts=Global parameter options for the function,
                             extraopts=local parameter options for the function to override the global ones,
                             indexlist=List giving the index of parameters where to take derivatives from (corresponding to $paramoptions order),
                             epsilonValue= epsilon step for numerical derivative,
                             options: 
                             derivativeOrder->First, Second or SecondMixed (passed as strings),
                             transformDerivativeFunction-> (any function known symbolically by Mathematica, such as Log or Sin)
oneSigmaErrors :: 1-sigma fully marginalized errors of the Fisher Matrix. 
                  Form: onesigmaErrors[matrix] = list

parameterPositionInBins :: Gives indices of the position inside the Fisher Matrix of the parameter given by index, 
                           index is the number of the parameter among the z-depdendent ones.
                           Form: parameterPositionInBins[number]=list
placeBlocks :: Add contribution block matrices to Fisher Matrix by using blocks. 
               Form: placeBlocks[am,cm,dm,pos,size] = matrix,
               Dimensions constraint: Width[cm]==Width[dm] && Length[cm]==Length[am],
               am=The upper left corner square matrix, 
               cm=Cross term rectangular block matrix, 
               dm=Diagonal square block matrix, 
               pos=Block position where to set the blocks cm and dm: 
               dm is set in the diagonal at (pos,pos), 
               cm is set in the block position (1,pos),
               Transpose[cm] is set in the block position (pos,1),
               size: Final desired matrix size
plotConfidenceRegions :: Plots the confidence region ellipses for the fully marginalized 2x2 matrices at indices a, b.
                         With the option fullFisherMatrixFalse, it plots only the ellipses for the given 2x2 Matrix. Then a, b can be omitted from input.
                         Form: plotConfidenceRegions[fish_,par_,a_:1,b_:1,opts:].
                         Accepts Options for Graphics.
                         Other options: parameterLabelsNone,labelStyleFunctionstyleFunction,lineThickness0.05, colorList$blueTones
PobsSymbolic :: PobsSymbolic[kfunc_,mufunc_,Dfunc_,Hfunc_,Dfuncref_,Hfuncref_,Gfunc_,bfunc_,betafunc_,P0func_,Pshotfunc_]
                 gives the symbolic form of the observed power spectrum

printVariable :: print stuff for debug


proper2Subset :: Gives the indices of all possible non-repeating subsets of size two from a list of parameters.
                 Form: proper2Subset[list]=nested list
                 Tip: Use paramslist[[#]]&/@proper2Subset[paramslist] to obtain the values of the subsets of paramlist
relativePkNoise :: Function of k and z, which determines an extra noise in the power spectrum, 
                   entering an effective ndens in the covariance

RfuncSymbolic :: RfuncSymbolic[muref_,Dfunc_,Hfunc_,Dfuncref_,Hfuncref_]
                  gives the symbolic R function that relates k and mu fiducial to a different cosmology

setEpsilon :: Option for derivative of the power spectrum


setFisherDimensions :: set dimensions of Fisher Matrix for package functions, 
                       nc=number of z-independent cosmological parameters,
                       nz=number of z-depdendent cosmological parameter functions,
                       nzbi=number of z-bins.
                       Form:setfisherDimensions[nc_,nz_,nzbi_]=Null
setFisherKmaxValues :: setFisherKmaxValues[kmaxhard, kmaxmethod, maxKsigma8:(optional)]: 
                       Sets the method for integration limits in k and sets the maximum k allowed $kmaxHard.
                       Optional: sets the maximum sigma8 used to calculate the maximum kmax limit for $kmaxMethod 1 and 2.
                       For a display of available kmax calculation methods, type: ?$kmaxMethod
setFisherParameterFiducials :: setFisherParameterFiducials[paramNames_,paramFiducials_] 
                               sets the global default list of rules between parameter names and parameter fiducial values
                               Not necessary to set if CosmologyFunctions package is loaded and variables are set there
setGalaxiesNumberBinsSpecs :: 
                              Sets the global lists of number of galaxies for each redshift bin and other specifications for the survey
                              setGalaxiesNumberBinsSpecs[zbin_?VectorQ,dzBinWidth_,volumeSky_,
                              dzRelVelerr_,photozerr_,intrShear_,frsky_,galdarcm_,z0GalDist]
setNumericalEpsilon :: set the global value of epsilon for the h step in numerical derivatives. 
                       Value can be set as an option individually for all derivative functions using setEpsilon->epsilon

sigma8FunctionForZK :: Accepts an interpolation function for powerspectrum[z,k] 
                       and calculates sigma8Function[powS,Pi/2/kvals] on a table of redshifts and k values

sigmaPkNoise :: Defines an extra error in Pk entering an effective ndens in the covariance, 
                by multiplying relativePkNoise with the observedPowerSpectrum. Accepts all the options for observedPowerSpectrum

styleFunction :: Form: styleFunction[text_, size_,opts:]. Accepts all the options of Style. Default: {FontFamily'Times', Bold}


toPercentage :: transform error to percentage. 
                Form: toPercentage[number] = string

totalGal :: Global variable counting the total galaxy


totalGalaxies ::  totalGalaxies[numb_,index_] counts total number of Galaxies in survey, 
                 when given number numb and returns cumulative total galaxies in bin index

transformDerivativeFunction :: Option for derivative of the power spectrum


volumeEffective :: calculate the 'effective volume' used in the Fisher Matrix calculation, 
                   using the ndens effective and the observed power spectrum

volumeSurvey :: gives the volume of a survey in the central value zmid of the redshift bins given by zbino
                volumeSurvey[zmid,zbino]=number. Needs to specify globally the covered sky fraction of the survey with Global
                
volumeEuclid

volumeZbin :: Gives the volume contained in a z bin given by z1<z<z2. 
              Uses the angular diameter distance and H0 from the CosmologyFunctions package

zaverage :: When given a list of initial bin positions, gives back a list of the average quantity in each bin
            Form: zaverage[list] = list

zAveragetoBins :: Gives back the list of redshift bin intervals after being provided with the central redshift of each bin 
                  and the bin width. Form: zAveragetoBins[zlist_,dz_]

zbinEnd :: Gives back the final value of the list of bins
           Form: zbinEnd[list] = number

zBinExtendLims :: Extend the end of the redshift bin limit by optional value extopt (default 0.2) and extend first bin value to 0


zbinStart :: Gives back the starting value of the list of bins
             Form: zbinStart[list] = number

$blueTones :: set of 3 default blue tones for ellipses


$compOptValuesStrList :: Option strings for the function complementOptionValues.
                         {CustomFullOptions,CustomComplementOptions,ExternalComplementOptions,ExternalUnsetOptions,DefaultOptions}

$epsilonstep :: Global value of epsilon step for numerical derivatives


$kmaxHard :: Maxmimum integration limit in k for Fisher Matrix


$kmaxInterpFunc :: Function which gives the maximum kmax value as a function of redshift, for the k integration of the Fisher Matrix


$kmaxMethod :: method for maximum integration limit in Fisher Matrix integration in k. Options:
               1: kmax=$kmaxHard is chosen always
               2: kmax= Chosen from kmaxInterpolationFunction[z] if smaller than $kmaxHard
               3: kmax= Chosen from kmaxInterpolationFunction[z] always
$maxKsigma8 ::  Default=0.35. 
               Maximum sigma8 which is used to calculate the allowed kmax as a function of redshift in the Fisher Matrix integration.

$ndensGalZFuncGlobal :: The global galaxy number density interpolation function used in the Fisher Matrix analysis


$numbGalaxiesGlob :: Global list with number of galaxies for each redshift bin


$zbinGlobal :: Global list with the values of the intervals for the redshift bins`

