(* Mathematica package *)
Clear[findAllRoots]
SyntaxInformation[
   findAllRoots] = {"LocalVariables" -> {"Plot", {2, 2}}, 
   "ArgumentsPattern" -> {_, _, OptionsPattern[]}};
SetAttributes[findAllRoots, HoldAll];

Options[findAllRoots] = 
  Join[{"ShowPlot" -> False, PlotRange -> All}, 
   FilterRules[Options[Plot], Except[PlotRange]]];

findAllRoots[fn_, {l_, lmin_, lmax_}, opts : OptionsPattern[]] := 
 Module[
  {pl, p, x, localFunction, brackets},
  localFunction = ReleaseHold[Hold[fn] /. l :> x];
  If[
   lmin != lmax,
   pl = Plot[localFunction, {x, lmin, lmax}, 
     Evaluate@
      FilterRules[Join[{opts}, Options[findAllRoots]], Options[Plot]]
     ];
   p = Cases[pl, Line[{x__}] :> x, Infinity];
   If[OptionValue["ShowPlot"], 
    Print[Show[pl, PlotLabel -> "Finding roots for this function", 
      ImageSize -> 200, BaseStyle -> {FontSize -> 8}]]],
   p = {}
   ];
  brackets = Map[
    First,
    Select[
     (* This Split trick pretends that two points on 
     the curve are "equal" if the function
     values have _opposite _ sign. Pairs of such 
     sign-changes form the brackets for the subsequent
     FindRoot *)
     Split[p, Sign[Last[#2]] == -Sign[Last[#1]] &],
     Length[#1] == 2 &
     ],
    {2}
    ];
  x /. Apply[FindRoot[localFunction == 0, {x, ##1}] &, 
     brackets, {1}] /. x -> {}
  ]