(* ::Package:: *)

dispersion = Function[{k, om}, Block[{Alpha},
	(*Print["Calculating first alpha ..."];
	*)
	Alpha = findAlpha[k, om];
	(*Print[MatrixForm[Alpha]]*);
	(*Print["Calculating remain alphas ..."]*);
	lastStep[Alpha, k, om]
 ]];
