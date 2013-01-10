lastStep = Function[{Alpha, k, \[CapitalOmega]}, Block[ {
    
    L, A, B, Fi, U, dU, ddU, dFi, deriv, exponenta, x, precision = 7, om, t,
    x1, x2, x3, n2, n3, al, z, A1, A2, A3, B1, fi, u1, u2, u3, u11, u12, u13, u21, u22, u23, u31, u32, u33,
    u1tt, u2tt, u3tt, fi1, fi2, fi3,
    c11, c33, c44, c12, c13, c14, e15, e22, e31, e33,es11, es33,
    S, h = 1,
    al1, al2, al3, al4, al5, al6, al7,al8,
    A11, A12, A13, A14, A15, A16, A17, A18,
    A21, A22, A23, A24, A25, A26, A27, A28,
    A31, A32, A33, A34, A35, A36, A37, A38,
    A41, A42, A43, A44, A45, A46, A47, A48,
    conv,
    
    Sigma, dSigma, sdSigma, Edge, Ldet, Lcoeff, Di, dDi, ones, sdDi, b,
s11lin, s12lin, s13lin,
s21lin, s22lin, s23lin ,
s31lin, s32lin, s33lin,
di1lin, di2lin, di3lin},

(* Local settings *)
	 al1 = Alpha[[1]];
	 al2 = Alpha[[2]];
	 al3 = Alpha[[3]];
	 al4 = Alpha[[4]];
	 al5 = Alpha[[5]];
	 al6 = Alpha[[6]];
	 al7 = Alpha[[7]];
	 al8 = Alpha[[8]];
	
	 (*-----linear task - search alfa----*)
	 deriv = Exp[-I * ( om * t - k1 * (n2*x2 + n3*x3))];
	 
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
	 n3 = Sqrt[0.5];
	 n2 = Sqrt[1/2];
	                                         
	                                       (* II part *)
	
	                                   (*-----------------------*)
	                                   (*---search dependences------*)
	
	                                   (* Need to rename *) 
	
	                                   (*A = Table[Subscript[alpha, i, j], {i, 4}, {j, 8}];*)
	 A = {
	     {A11, A12, A13, A14, A15, A16, A17, A18},
	     {A21, A22, A23, A24, A25, A26, A27, A28},
	     {A31, A32, A33, A34, A35, A36, A37, A38}
	 };
	 B = {A41, A42, A43, A44, A45, A46, A47, A48};
	 U = A . (Exp[Alpha * x1]) * deriv;
	 Fi = B . (Exp[Alpha * x1]) * deriv;
	 dU = D[U, {{x1, x2, x3}}];
	 dFi = D[Fi, {{x1, x2, x3}}];
	 ddU = D[U, {t, 2}];
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
	
	                                                                                        (*TODO: need to convert to matrix form *)
	 s11lin = c11 * u11 + c12 * u22 + c13 * u33 + c14 * u32 + c14 * u23 + e22 * fi2 + e31 * fi3;
	 s12lin = c14 * u31 + c14 * u13 + (c11-c12)/2 * u21 + (c11-c12)/2 * u12 + e22 * fi1;
	 s13lin = c44 * u13 + c44 * u31 + c14 * u21+c14*u12+e15*fi1;
	 s21lin = c14 * u31 + c14 * u13 + (c11-c12)/2 * u21 + (c11-c12)/2*u12+e22*fi1;
	 s22lin = c12 * u11 + c11 * u22 + c13 * u33-c14*u32-c14*u23-e22*fi2+e31*fi3;
	 s23lin = c44 * u23 + c44 * u32 + c14 * u11-c14*u22+e15*fi2;
	 s31lin = c44 * u13 + c44 * u31 + c14 * u21+c14*u12+e15*fi1;
	 s32lin = c44 * u23 + c44 * u32 + c14 * u11-c14*u22+e15*fi2;
	 s33lin = c13 * u11 + c13 * u22 + c33 * u33+e33*fi3;
	 di1lin = e15 * u31 + e15 * u13 + e22 * u21 + e22 * u12-es11*fi1;
	 di2lin = e15 * u32 + e15 * u23 + e22 * u11 - e22 * u22-es11*fi2;
	 di3lin = e31 * u11 + e31 * u22 + e33 * u33 - es33*fi3;
	 Sigma = {
	     { s11lin, s12lin, s13lin },
	     { s21lin, s22lin, s23lin },
	     { s31lin, s32lin, s33lin }
	 };
	 dSigma = Transpose[{
	     D[Sigma[[All, 1]], x1],
	     D[Sigma[[All, 2]], x2], 
	     D[Sigma[[All, 3]], x3]
	 }];
	 ones = Table[1, {3}];
	 sdSigma = dSigma . ones;
	 Di = {di1lin, di2lin, di3lin};
	 dDi = {
	     D[Di[[1]], x1],
	     D[Di[[2]], x2],    
	     D[Di[[3]], x3]
	 };
	 sdDi = dDi . ones;
	 
	 L = Append[(\[Rho] * ddU - sdSigma), sdDi] /.  deriv -> 1 /. om ^ 2 \[Rho] -> \[CapitalOmega];
(*	 Print[L];
*)	 

	Print["Solving Matrix coeffs: "];
	conv = Function[{exp, alpha, vars}, Block[ {res, coeff, ealx1},
		Print["."];
	     {res, coeff} = CoefficientArrays[exp /. Exp[alpha * x1] -> ealx1, ealx1];
	     Solve[{coeff[[1,1]] == 0, coeff[[2,1]] == 0, coeff[[3,1]] == 0}, vars]
	 ]];
	
	 S = conv[L, al1, {A11, A21, A31}];
	 A11 = A11 /. S[[1]];
	 A21 = A21 /. S[[1]];
	 A31 = A31 /. S[[1]];
	 S = conv[L, al2, {A12, A22, A32}];
	 A12 = A12/.S[[1]];
	 A22 = A22/.S[[1]];
	 A32 = A32/.S[[1]];
	 S = conv[L, al3, {A13, A23, A33}];
	 A13 = A13/.S[[1]];
	 A23 = A23/.S[[1]];
	 A33 = A33/.S[[1]];
	 S = conv[L, al4, {A14, A24, A34}];
	 A14 = A14/.S[[1]];
	 A24 = A24/.S[[1]];
	 A34 = A34/.S[[1]];
	 S = conv[L, al5, {A15, A25, A35}];
	 A15 = A15/.S[[1]];
	 A25 = A25/.S[[1]];
	 A35 = A35/.S[[1]];
	 S = conv[L, al6, {A16, A26, A36}];
	 A12 = A12/.S[[1]];
	 A22 = A22/.S[[1]];
	 A32 = A32/.S[[1]];
	 S = conv[L, al7, {A17, A27, A37}];
	 A17 = A17/.S[[1]];
	 A27 = A27/.S[[1]];
	 A37 = A37/.S[[1]];
	 S = conv[L, al8, {A18, A28, A38}];
	 A18 = A18/.S[[1]];
	 A28 = A28/.S[[1]];
	 A38 = A38/.S[[1]];
	
	                                       (*----------Search k and coefficients from boundary conditions-------*)
	 Print["Calculating edge ..."];                                      
	 Edge = U /. deriv -> 1 /. om ^ 2 \[Rho] -> \[CapitalOmega];
	 Edge = Append[Edge, Fi /. deriv -> 1 /. om ^ 2 \[Rho] -> \[CapitalOmega]];
	 Edge = Join[Edge /. x1 -> h, Edge /. x1 -> -h];
	 {b, Lcoeff} = CoefficientArrays[Edge, {A41, A42, A43, A44, A45, A46, A47, A48}];
	 
	 Print["Calculating det of edge ..."];
	 Print[MatrixForm[SetPrecision[Lcoeff,2]]];	                                       
	 Ldet = Det[SetPrecision[Lcoeff,2]]
]];