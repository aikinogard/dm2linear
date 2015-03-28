gaussian2[alpha1_,lmn1_,A_,alpha2_,lmn2_,B_]:=Module[{l1=lmn1[[1]],m1=lmn1[[2]],n1=lmn1[[3]],l2=lmn2[[1]],m2=lmn2[[2]],n2=lmn2[[3]],Ax=A[[1]],Ay=A[[2]],Az=A[[3]],Bx=B[[1]],By=B[[2]],Bz=B[[3]]},Re[(\!\( 
\*SubsuperscriptBox[\(\[Integral]\), \(-\[Infinity]\), \(\[Infinity]\)]\( 
\*SuperscriptBox[\((x - Ax)\), \(l1\)] Exp[\(-alpha1\) 
\*SuperscriptBox[\((x - Ax)\), \(2\)]] 
\*SuperscriptBox[\((x - Bx)\), \(l2\)] Exp[\(-alpha2\) 
\*SuperscriptBox[\((x - Bx)\), \(2\)]] \[DifferentialD]x\)\)) (\!\( 
\*SubsuperscriptBox[\(\[Integral]\), \(-\[Infinity]\), \(\[Infinity]\)]\( 
\*SuperscriptBox[\((y - Ay)\), \(m1\)] Exp[\(-alpha1\) 
\*SuperscriptBox[\((y - Ay)\), \(2\)]] 
\*SuperscriptBox[\((y - By)\), \(m2\)] Exp[\(-alpha2\) 
\*SuperscriptBox[\((y - By)\), \(2\)]] \[DifferentialD]y\)\))(\!\( 
\*SubsuperscriptBox[\(\[Integral]\), \(-\[Infinity]\), \(\[Infinity]\)]\( 
\*SuperscriptBox[\((z - Az)\), \(n1\)] Exp[\(-alpha1\) 
\*SuperscriptBox[\((z - Az)\), \(2\)]] 
\*SuperscriptBox[\((z - Bz)\), \(n2\)] Exp[\(-alpha2\) 
\*SuperscriptBox[\((z - Bz)\), \(2\)]] \[DifferentialD]z\)\))]]

(*Calculate benchmark value: integral of product of two Gaussians*)
(*Example:*)
(*Setting in .bash_profile*)
(*alias mathscript="/Applications/Mathematica.app/Contents/MacOS/MathematicaScript"*)
(*mathscript -script gaussian2.m {0.3,{0,1,1},{1.,2.,3.},0.1,{1,1,0},{1.,1.,1.}}*)
(*Print[gaussian2[0.3,{0,1,1},{1.,2.,3.},0.1,{1,1,0},{1.,1.,1.}]]*)
input = Flatten[ToExpression[Rest[$ScriptCommandLine]]]

If[Length[input]!=14,Print["Error: Please check input! Number of inputs:"];Print[Length[input]];Print[input];Exit[],0]

alpha1=ToExpression[input[[1]]]
order1=ToExpression[{input[[2]],input[[3]],input[[4]]}]
A1=ToExpression[{input[[5]],input[[6]],input[[7]]}]
alpha2=ToExpression[input[[8]]]
order2=ToExpression[{input[[9]],input[[10]],input[[11]]}]
A2=ToExpression[{input[[12]],input[[13]],input[[14]]}]

Print[FortranForm[gaussian2[alpha1,order1,A1,alpha2,order2,A2]]]