(* ::Package:: *)


(* ::Section:: *)
(* Region Title *)
lastStep = Function[{Alpha, k, \[CapitalOmega]}, Block[ {
    
    L, A, B, Fi, U, dU, ddU, dFi, deriv, ex, x, precision = 7, om, t,
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
    
    Sigma, sigmas, dSigma, sdSigma, Edge, Ldet, Lcoeff, Di, dDi, ones, sdDi, b,
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
	 
	n3 = Sqrt[1/2];
	n2 = Sqrt[1/2];

	
	 (*-----linear task - search alfa----*)
	 deriv = Exp[-I * ( om * t - k * (n2*x2 + n3*x3))];
	 ex = (Exp[Alpha * x1]);

	 x = { x1, x2, x3 };
	 
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


    U = A . ex * deriv;
    Fi = B . ex * deriv;

    sigmas = buildSigmas[U, Fi, x, t];

	L = sigmas /.  deriv -> 1 /. om ^ 2 \[Rho] -> \[CapitalOmega];

(*	Print["Solving Matrix coeffs: "];
*)
	conv = Function[{exp, alpha, vars}, Block[ {res, coeff, ealx1},
		(*Print["."];*)
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
(*	 Print["Calculating edge ..."];                                      
*)	 Edge = U /. deriv -> 1 /. om ^ 2 \[Rho] -> \[CapitalOmega];
	 Edge = Append[Edge, Fi /. deriv -> 1 /. om ^ 2 \[Rho] -> \[CapitalOmega]];
	 Edge = Join[Edge /. x1 -> h, Edge /. x1 -> -h];
	 {b, Lcoeff} = CoefficientArrays[Edge, {A41, A42, A43, A44, A45, A46, A47, A48}];
	 
(*	 Print["Calculating det of edge ..."];
*)
Ldet = Det[SetPrecision[Lcoeff,2]]

(*d = Det[Table[Subscript[a, i, j], {i, 1, 8}, {j, 1, 8}]];
Do[d = d /. Subscript[a, i, j] ->  Lcoeff[[i, j]], {i, 1, 8}, {j, 1, 8}];
d
*)]];
