(* ::Package:: *)

dispersion = Function[{k, om}, Block[{Alpha},
	PrintTemporary["Calculating first alpha ..."];	
	PrintTemporary[k];	
	PrintTemporary[om];	
	Alpha = findAlpha[k, om];
	(*Print[MatrixForm[Alpha]]*);
	(*Print["Calculating remain alphas ..."]*);
	findCoeffsAndDet[Alpha, k, om]
 ]];
