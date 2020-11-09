(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: CMBTools *)
(* :Context: CMBTools` *)
(* :Author: santiago *)
(* :Date: 2018-08-08 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2018 santiago *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["CMBTools`"]
(* Exported symbols added here with SymbolName::usage *)

CMBCls::usage="CMB Cls for TT, TB, TE, EE and BB computed from an external Boltzmann code. Array of
{\ell, Cl}, accessed as: CMBCls['XY'], where X and Y represent temperature or polarization components.";

CMBnoise::usage="Noise term for CMB forecasts.";

CMBCovarianceMat::usage="CMB Covariance matrix for all T, E, B, components.";
CMBJMat::usage="CMB dCls/dp array for all T, E, B, components.";

CMBFisher::usage="Fisher matrix calculation for CMB forecasts."  (* TODO: program this*)

$fisherMatrixCMB::usage="Computed Fisher matrix for CMB forecasts.";

$probesNumCMB::usage="Number of CMB probes: {TT, TE, TB, EE, BB}.";

(*probes convention: {1:TT, 2:EE, 3:TE, 4:BB, 5:TB}*)


Begin["`Private`"]

Options[CMBFullCls]=$paramoptions~Join~{cijFunction->ClsCMB, noiseFunction->CMBnoise}
CMBFullCls[ell_, xy_, opts:OptionsPattern[]]:=CMBFullCls[ell, xy, opts]=Block[{arra, ClsFunc=OptionValue[cijFunction], noiCls=OptionValue[noiseFunction]},
  ClsFunc[ell,xy]+ noiCls[ell,xy]
]


Options[CMBcovarianceMat]=$paramoptions~Join~{CMBProbesNum->$probesNumCMB};

CMBcovarianceMat[ell_, opts:OptionsPattern[]]:=CMBcovarianceMat[ell, opts]=Block[{Nprobes=OptionValue[CMBprobesNum], numzero=10^(-11), cc=CMBFullCls, cmbcov},
  cmbcov=ConstantArray[0,{Nprobes,Nprobes}];
  cmbcov[[1,1]]=cc[ell,1]^2 - 2(cc[ell,3]*cc[ell,5])^2/(cc[ell,2]*cc[ell,4]+numzero);
  cmbcov[[2,2]]=cc[ell,2]^2;
  cmbcov[[3,3]]=(cc[ell,3]^2 + cc[ell,1]*cc[ell,2])/2 - cc[ell,2]*(cc[ell,5])^2/(2*cc[ell,4]+numzero);
  cmbcov[[4,4]]=cc[ell,4]^2;
  cmbcov[[5,5]]=(cc[ell,5]^2 + cc[ell,1]*cc[ell,4])/2 - cc[ell,4]*(cc[ell,3])^2/(2*cc[ell,2]+numzero);
  cmbcov[[1,2]]=cc[ell,3]^2;
  cmbcov[[1,4]]=cc[ell,5]^2;
  cmbcov[[1,3]]=cc[ell,1]*cc[ell,3]-(cc[ell,3]cc[ell,5]^2/(cc[ell,4]+numzero));
  cmbcov[[1,5]]=cc[ell,1]*cc[ell,5]-(cc[ell,5]cc[ell,3]^2/(cc[ell,2]+numzero));
  cmbcov[[2,3]]=cc[ell,2]*cc[ell,3];
  cmbcov[[4,5]]=cc[ell,4]*cc[ell,5];
  Return[cmbcov]
]



Options[CMBJmat]=$paramoptions~Join~{dcijdpFunction->dClsdpCMB, CMBProbesNum->$probesNumCMB};

CMBJmat[ell_?NumericQ,alpha_,opts:OptionsPattern[]]:=CMBJmat[ell,alpha, opts]=Block[{ret, dClsdpFunc, Nprobes=OptionValue[CMBprobesNum]},
  dClsdpFunc=OptionValue[dcijdpFunction];
  ret=Table[dClsdpFunc[ell, xy, alpha], {xy, Nprobes}];
  Return[ret]
];



Options[CMBFisher]=$paramoptions~Join~{ellmin:>$ellmin,ellmax:>$ellmax,Deltaell:>$deltaell, cijFunction->ClsCMB, dcijdpFunction->dClsdpCMB};

CMBFisher[alpha_,beta_,opts:OptionsPattern[]]:=FisherWL[alpha,beta,opts]=If[alpha<beta,FisherWL[beta,alpha,opts],
  Block[{ClsFunc, dClsdpFunc, lmin=OptionValue[ellmin], lmax=OptionValue[ellmax], dell=OptionValue[Deltaell]},
    Global`prog1=Global`prog2=Global`prog3=0;
    SetSharedVariable[Global`prog1=Global`prog2=Global`prog3];
    $logellrange;
    $ellbinscenters;
    dClsdpFunc=OptionValue[dcijdpFunction];
    ClsFunc=OptionValue[cijFunction];
    Sum[With[{ell=$ellbinscenters[[i]] , delell=(tentox@($logellrange[[i+1]])-tentox@($logellrange[[i]]))},
      Global`prog1++;
      ($fsky/2)*(2 ell+1)*delell*
          Tr[CMBJmat[ell,alpha, dcijdpFunction->dClsdpFunc].CMBinversecovariance[ell, cijFunction->ClsFunc].CMBJmat[ell,beta, dcijdpFunction->dClsdpFunc].inversecovariance[ell,cijFunction->ClsFunc]]],
      {i,Length@$ellbinscenters-1}]]
];



End[] (* `Private` *)

EndPackage[]