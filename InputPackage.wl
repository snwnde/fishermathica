(* ::Package:: *)

(* ::Code::Initialization:: *)
SGKernel[left_?NonNegative,right_?NonNegative,degree_?NonNegative,derivative_?NonNegative]:=Module[{i,j,k,l,matrix,vector},matrix=Table[(*matrix is symmetric*)l=i+j;
If[l==0,left+right+1,(*Else*)Sum[k^l,{k,-left,right}]],{i,0,degree},{j,0,degree}];
vector=LinearSolve[matrix,MapAt[1&,Table[0,{degree+1}],derivative+1]];
(*vector=Inverse[matrix][[derivative+1]];*)Table[vector.Table[If[i==0,1,k^i],{i,0,degree}],{k,-left,right}]]/;derivative<=degree<=left+right


SGSmooth[list_?VectorQ,window_,degree_,derivative_:0]:=Module[{pairs},pairs=MapThread[List,{Range[Length[list]],list}];
Map[Last,SGSmooth[pairs,window,degree,derivative]]]


SGSmooth[list_?MatrixQ,window_,degree_,derivative_:0]:=Module[{kernel,list1,list2,margin,space,smoothData},(*determine a symmetric margin at the ends of the raw dataset.The window width is split in half to make a symmetric window around a data point of interest*)margin=Floor[window/2];
(*take only the 1st column of data in the list to be smoothed (the independant Values) and extract the data from the list by removing half the window width'i.e.,margin' from the ends of the list*)list1=Take[Map[First,list],{margin+1,Length[list]-margin}];
(*take only the 2nd column of data in the list to be smoothed (the dependent Values) and Map them into list2*)list2=Map[Last,list];
(*get the kernel coefficients for the left and right margins,the degree,and the requested derivative*)kernel=SGKernel[margin,margin,degree,derivative];
(*correlation of the kernel with the list of dependent values*)list2=ListCorrelate[kernel,list2];
(*Data _should_ be equally spaced,but... calculate spacing anyway by getting the minimum of all the differences in the truncated list1,remove the first and last points of list1*)space=Min[Drop[list1,1]-Drop[list1,-1]];
(*condition the dependant values based on spacing and the derivative*)list2=list2*(derivative!/space^derivative);
(*recombine the correlated (x-y) data pairs (that is list1 and list2),put these values back together again to form the smooth data list*)smoothData=Transpose[{list1,list2}]]/;derivative<=degree<=2*Floor[window/2]&&$VersionNumber>=4.0


performSGdewiggling[krange_?ListQ, observable_?ListQ]:=Block[{oorder, wwin, externalKPk, 
externalPk, origKs, cutpos,PkdataCut,extPkInterp, savgolPkCut, newKs, pkdCInterp, SGpkCInterp, ratiosTest, ratspos, intermKs, externalNewPkNW},
oorder=4;
wwin=123;
origKs=krange;  (*observable has the form {{z,k},Pk} *)
externalPk=observable;
cutpos = Position[origKs,_?(origKs[[1]]<=#<5.0&)];
externalKPk=Transpose[{origKs,externalPk}];
PkdataCut = Extract[externalKPk, cutpos];
savgolPkCut = SGSmooth[PkdataCut,wwin,oorder];
newKs = savgolPkCut[[All,1]];
extPkInterp=Interpolation[externalKPk];
pkdCInterp=Interpolation[PkdataCut];
SGpkCInterp = Interpolation[savgolPkCut];
(*ratiosTest = Table[SGpkCInterp[kak]/pkdCInterp[kak], {kak,newKs}]*);
(*ratspos=Position[ratiosTest,_?(0.95<=#<=1.05&)]*);
(*intermKs=Extract[newKs, ratspos]*);
externalNewPkNW=Table[{kki, Which[kki<newKs[[1]],
  extPkInterp[kki],
newKs[[1]]<=kki<=newKs[[-1]],
SGpkCInterp[kki],
kki>newKs[[-1]],
extPkInterp[kki]]},{kki,origKs}];
Return[externalNewPkNW]
];


paramsFolderFunc[param_, epsilon_]:=Block[{pm, string, parstring, epsstr, epsi}, 
                                          pm=Which[Sign[epsilon]==-1, "_"<>plusminstr[[1]], Sign[epsilon]==1, "_"<>plusminstr[[2]], Sign[epsilon]==0, ""];
                                          epsi=If[Chop[epsilon, 10^-6]==0., 0, epsilon];
                                          If[param=="fiducial" && epsi!=0, Return[None]];
                                          If[param!="fiducial" && epsi==0, Return[None]];
                                          If[param=="fiducial", pm=""];
                                          epsstr=epsstrfun[Abs[epsi]];
										  string=(param<>pm<>"_"<>epsstr);						
                                          Return[string]
                                          ];


ignoreParameters[parameps_]:=Block[{ignolist=$ignorepars, return}, 
                    return=If[MemberQ[ignolist, StringSplit[parameps, "_"][[1]]]==True, 
                                                       paramsFolderFunc["fiducial", 0],
                                                       parameps]
                                                       ]


externalInputFile[par_, observable_, format_:".txt"]:=Block[{para,obsstr},
para=ignoreParameters[par];
If[MemberQ[$observablesNames,observable]==False, Print["Observable string is not supported"]; Abort[]];
obsstr=observable/.$rulesfuncsInputNames;
((ToString[para]<>"/"<>obsstr<>format))
]


kzRangeFile[par_, grid_String, format_:".txt"]:=Block[{para, gfile},
para=ignoreParameters[par];
gfile=grid<>"_values_list";
((ToString[para]<>"/"<>gfile<>format))
]


externalMatrixListImport[obsname_,param_,epsilon_]:=externalMatrixListImport[obsname,param,epsilon]=Block[{folderpar, tab, dim},
folderpar=paramsFolderFunc[param, epsilon];
tab=Import[(inputDataDir<>externalInputFile[folderpar, obsname]), "Table"];
dim=Dimensions[tab];
If[dim[[2]]==1,
tab=Flatten[tab,1];  (*If imported Table is just a list, flatten to have a 1-d List*)
];
Return[tab]]


observableParamDependency[obsname_,par_,zzi_,kki_,epsilonlist_:$epslistPM]:=Block[{varfun, fiducialfun, paramFuncTab, var, fiduval},
            If[StringMatchQ[obsname, "*zk*"]==True,
            (*Scale-dependence and z-dependence*)
            varfun = externalObservableFunc[obsname,par, #][zzi,kki]&;            
            fiducialfun=externalObservableFunc[obsname,"fiducial", 0.][zzi,kki];
            ,
            varfun = externalObservableFunc[obsname,par, #][zzi]&;            
            fiducialfun=externalObservableFunc[obsname,"fiducial", 0.][zzi];
            ];
            fiduval=$parameterFiducials[par];
            paramFuncTab = Table[{fiduval*(1+ee), 
                            If[ee==0.,
                            var=fiducialfun,
			                var=varfun[ee] ];
            var/fiducialfun}, {ee, epsilonlist}];
            Return[paramFuncTab]]


twoPointDerivative[obsname_,param_,zzi_,kki_, epsi_:0.00625, relative_:True]:=Block[{ee, twopoint, relder, retu, varfun,fiducialfun, fiduval}, 
If[StringMatchQ[obsname, "*zk*"]==True,
            (*Scale-dependence and z-dependence*)
            varfun = externalObservableFunc[obsname,param, #][zzi,kki]&;            
            fiducialfun=externalObservableFunc[obsname,"fiducial", 0.][zzi,kki];
            ,
            varfun = externalObservableFunc[obsname,param, #][zzi]&;            
            fiducialfun=externalObservableFunc[obsname,"fiducial", 0.][zzi];
            ];
fiduval=$parameterFiducials[param];
twopoint=(varfun[epsi]-varfun[-epsi])/(2*epsi*fiduval);
relder=twopoint/fiducialfun;
retu=If[relative==True,
relder,
twopoint];
Return[retu]]


steMderiv[partable_,x0_, order_:1]:=Block[{lm,rmse,relativederiv,halflen, cutptab, tolerance=0.005, polyst},
            If[order<1, Print["Polynomial Order must be larger or equal 1"];Abort[]];
             polyst=Table[x^ii,{ii,Range[0,order]}];
             halflen=Range[Floor[Length[partable]/2]];
             Do[
                cutptab=partable[[ii;;-ii]];
				lm=LinearModelFit[cutptab, polyst,x];
				rmse=Sqrt[Total[(lm["FitResiduals"])^2]];
				relativederiv=(D[lm[x],x]/.x->x0);
				If[rmse <= tolerance, 
				Break[]]
				,
				{ii,halflen}];
			Return[{relativederiv, rmse, lm, cutptab}]
			]


Options[observableDerivativeZK]={derivativeType->"SteM1", twopointepsilon->$twpointepsilonstep};


derivativeType::usage="Option for observableDerivativeZK. Options: '2pt'->Two point derivative, 'Stem1'->SteM derivative fitted with order 1, 'Stem4': SteM derivative fitted with order 4.";


twopointepsilon::usage="Option for observableDerivativeZK. Value of the epsilon step used in the 2-point derivatives.";


observableDerivativeZK[obsname_,parname_,zzi_,kki_, opts:OptionsPattern[]]:=Block[{fiduobs,fiduval,pptab,rmse=-1.,lml,cutpartab,
																	twoptder, epsvalue,
                                                                    zgrid=$zrangePowerSpectrum, zlenhalf,
                                                                    kgrid=$krangePowerSpectrumFiducial,klenhalf,kindex,
                                                                    dertype, relderiv, deriv,
                                                                    derivElem},
             epsvalue=OptionValue[twopointepsilon];                                                       
            dertype=OptionValue[derivativeType];
            If[MemberQ[{"SteM1","SteM4","2pt"},dertype]==False, Abort[]];
            klenhalf=Floor[(Length@kgrid)/2];
            zlenhalf=Floor[(Length@zgrid)/2];
             If[StringMatchQ[obsname, "*zk*"]==True,
            (*Scale-dependence and z-dependence*)
            fiduobs=externalObservableFunc[obsname,"fiducial", 0.][zzi,kki];
            kindex=kki;
            ,
            (*Only z-dependence*)
            fiduobs=externalObservableFunc[obsname,"fiducial", 0.][zzi];
            kindex = kgrid[[klenhalf]] (*have a fake kki*)
            ];
            Which[dertype=="2pt",
            relderiv=twoPointDerivative[obsname,parname,zzi,kki, epsvalue, True];
            ,
            dertype=="SteM1",
            pptab=observableParamDependency[obsname,parname,zzi,kki];
            fiduval=$parameterFiducials[parname];
            {relderiv,rmse,lml,cutpartab}=steMderiv[pptab, fiduval];
            ,
            dertype=="SteM4",
            pptab=observableParamDependency[obsname,parname,zzi,kki];
            fiduval=$parameterFiducials[parname];
            {relderiv,rmse,lml,cutpartab}=steMderiv[pptab,fiduval,4];
            ];
           If[zzi==zgrid[[zlenhalf]] && kindex==kgrid[[klenhalf]],
				Print["z=", zzi];
				Print["k=", kki];
				Print["RMSE residual: ", rmse];  (*Will print -1 for 2pt deriv*)
				Print["Relative "<>dertype<>" tangent at fiducial: ", relderiv];
				];
			 If[StringMatchQ[obsname, "*zk*"]==True,
            (*Scale-dependence and z-dependence*)
			 derivElem = {{zzi,kki},relderiv*fiduobs} ;
			 ,
			 (*Only z-dependence*)
			 derivElem = {zzi,relderiv*fiduobs} ;
			];
			Return[derivElem]
			];
