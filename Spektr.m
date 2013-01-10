spektr = Function[{k, om}, Block[{Alpha},
	Print["Calculating first alpha ..."];
	Alpha = alphaForOm[k, om];
	Print["Calculating remain alphas ..."];
	lastStep[Alpha, k, om]
 ]];