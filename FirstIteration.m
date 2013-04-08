(* ::Package:: *)

buildSigmas = Function[{U, Fi, x, t}, Block[{
			  Sigma, dSigma, sdSigma, Edge, Ldet, Lcoeff, Di, dDi, ones, sdDi, b,
		  c11, c33, c44, c12, c13, c14, e15, e22, e31, e33,es11, es33,
		  u11, u12, u13, u21, u22, u23, u31, u32, u33,
			u1tt, u2tt, u3tt, fi1, fi2, fi3,		  
s11lin, s12lin, s13lin,
s21lin, s22lin, s23lin ,
s31lin, s32lin, s33lin,
di1lin, di2lin, di3lin,
sigmas, dU, dFi, ddU },
		dU = D[U, { x }];	

	dFi = D[Fi, { x }];
	ddU = D[U[[1;;3]], {t, 2}];
	
	u11 = dU[[1,1]];
	
	u12 = dU[[1,2]];
	u13 = dU[[1,3]];
	u21 = dU[[2,1]];
	u22 = dU[[2,2]];
	u23 = dU[[2,3]];
	u31 = dU[[3,1]];
	u32 = dU[[3,2]];
	u33 = dU[[3,3]];
	
	u1tt = ddU[[1]];
	u2tt = ddU[[2]];
	u3tt = ddU[[3]];
	
	fi1 = dFi[[1]];
	fi2 = dFi[[2]];
	fi3 = dFi[[3]];
	
	c11 = 2.02999999999999969482421875`66.;
	c33 = 2.424`66.;
	c44 = 0.595`66.;
	c12 = 0.5729999999999999237060546875`66.;
	c13 = 0.752`66.;
	c14 = 0.0850000000000000095367431640625`66.;
	e15 = 3.7599999999999997868371792719699442386627197265625`66.;
	e22 = -2.430000000000000159872115546022541821002960205078125`66.;
	e31 = 0.2300000000000000099920072216264088638126850128173828125`66.;
	e33 = 1.3300000000000000710542735760100185871124267578125`66.;
	es11 = 3.9199999999999999049997873677556904647989455270362668670713901519775390625`66.;
	es33 = 2.47000000000000025930105900122727390078036791010163142345845699310302734375`66.;
	
	(*Sigma calc*)
	s11lin = c11 * u11 + c12 * u22 + c13 * u33 + c14 * u32 + c14 * u23 + e22 * fi2 + e31 * fi3;
	s12lin = c14 * u31 + c14 * u13 + (c11-c12) / 2 * u21 + (c11-c12) / 2 * u12 + e22 * fi1;
	s13lin = c44 * u13 + c44 * u31 + c14 * u21 + c14 * u12 + e15 * fi1;
	s21lin = c14 * u31 + c14 * u13 + (c11-c12) / 2 * u21 + (c11-c12) / 2 * u12 + e22 * fi1;
	s22lin = c12 * u11 + c11 * u22 + c13 * u33 - c14 * u32 - c14 * u23 - e22 * fi2 + e31 * fi3;
	s23lin = c44 * u23 + c44 * u32 + c14 * u11 - c14 * u22 + e15 * fi2;
	s31lin = c44 * u13 + c44 * u31 + c14 * u21 +c14 * u12 + e15 * fi1;
	s32lin = c44 * u23 + c44 * u32 + c14 * u11 -c14 * u22 + e15 * fi2;
	s33lin = c13 * u11 + c13 * u22 + c33 * u33 + e33 * fi3;
	di1lin = e15 * u31 + e15 * u13 + e22 * u21 + e22 * u12 - es11 * fi1;
	di2lin = e15 * u32 + e15 * u23 + e22 * u11 - e22 * u22 - es11 * fi2;
	di3lin = e31 * u11 + e31 * u22 + e33 * u33 - es33 * fi3;
	
	Sigma = {
	    { s11lin, s12lin, s13lin },
	    { s21lin, s22lin, s23lin },
	    { s31lin, s32lin, s33lin }
	};
		
	Di = { di1lin, di2lin, di3lin };

	dSigma = Transpose[{
	    D[Sigma[[All, 1]], x[[1]]],
	    D[Sigma[[All, 2]], x[[2]]], 
	    D[Sigma[[All, 3]], x[[3]]]
	}];

	ones = Table[1, {3}];
	sdSigma = dSigma . ones;
	
	dDi = {
	    D[Di[[1]], x[[1]]],
	    D[Di[[2]], x[[2]]],    
	    D[Di[[3]], x[[3]]]
	};
	

	sdDi = dDi . ones;
	Append[(\[Rho] * ddU - sdSigma), sdDi]
]];

findAlpha = Function[{k, \[CapitalOmega]}, Block[{
		Alpha, L, A, U, Fi, deriv, exponenta, x, precision, om, t,
		x1, x2, x3, n2, n3, al, z, A1, A2, A3, B1,
		 B, sigmas, LXdet, LXSolution, rSol, Ldet, Lcoeff, b, \[Rho]
	},

	If[ArrayQ[k], Print[k]];

    n3 = Sqrt[1/2];
	n2 = Sqrt[1/2];
		
	(* Local settings *)
	precision = 7;
	
	(*-----linear task - search alfa----*)
	deriv = Exp[-I * (om * t - k * (n2 * x2 + n3 * x3))];
	exponenta = Exp[al x1] * deriv;

	x = { x1, x2, x3 };
	
	
	
	A = { A1, A2, A3 };
	B = { B1 };
	
	U = A * exponenta;
	Fi = B1 * exponenta;	

	sigmas = buildSigmas[U, Fi, x, t];
	
	L = sigmas /. exponenta -> 1 /. om ^ 2 \[Rho] -> \[CapitalOmega];
(*Print["L: "];
Print[MatrixForm[SetPrecision[L,2]]];
*)	
	{b, Lcoeff} = CoefficientArrays[L, {A1, A2, A3, B1}];
(*Print["Lcoeff: "];
Print[MatrixForm[SetPrecision[Lcoeff,2]]];
*)
	Ldet = Det[Lcoeff];
(*Print["Ldet: "];
Print[MatrixForm[SetPrecision[Ldet,2]]];
*)	
	(*Use canoncial form*)
	LXdet = Collect[Ldet,
	    { al, al ^ 2 , al ^ 3, al ^ 4, al ^ 5, al ^ 6, al ^ 7, al ^ 8 }
	] /. al ^ 2 -> z /. al ^ 4 -> z ^ 2 /. al ^ 6 -> z ^ 3 /. al ^ 8 -> z ^ 4;
(*Print["LXdet: "];
Print[MatrixForm[SetPrecision[LXdet,2]]];
*)
(* ,WorkingPrecision->precision *)
	LXSolution = Solve[LXdet == 0, z];
(*Print["LXSolution: "];
Print[MatrixForm[SetPrecision[LXSolution,2]]];
*)	
	(* Result for I part *)
	rSol = Sqrt[z /. LXSolution];
	Join[rSol, -rSol]
]];
