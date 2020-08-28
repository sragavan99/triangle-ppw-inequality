(* ::Package:: *)

BeginPackage["Section3`"]

A30int::usage = "Takes p, q, x, y as input and outputs the function that needs to be integrated to obtain A30(p, q).";
B30int::usage = "Takes p, q, x, y as input and outputs the function that needs to be integrated to obtain B30(p, q).";
C30int::usage = "Takes p, q, x, y as input and outputs the function that needs to be integrated to obtain C30(p, q).";
D30int::usage = "Takes p, q, x, y as input and outputs the function that needs to be integrated to obtain D30(p, q).";
E30int::usage = "Takes p, q, x, y as input and outputs the function that needs to be integrated to obtain E30(p, q).";
F30int::usage = "Takes p, q, x, y as input and outputs the function that needs to be integrated to obtain F30(p, q).";

A45int::usage = "Takes p, q, x, y as input and outputs the function that needs to be integrated to obtain A45(p, q).";
B45int::usage = "Takes p, q, x, y as input and outputs the function that needs to be integrated to obtain B45(p, q).";
C45int::usage = "Takes p, q, x, y as input and outputs the function that needs to be integrated to obtain C45(p, q).";
D45int::usage = "Takes p, q, x, y as input and outputs the function that needs to be integrated to obtain D45(p, q).";
E45int::usage = "Takes p, q, x, y as input and outputs the function that needs to be integrated to obtain E45(p, q).";
F45int::usage = "Takes p, q, x, y as input and outputs the function that needs to be integrated to obtain F45(p, q).";

TotalIntegral::usage = "Takes a function and (p, q) as input and outputs the integral of that function over the triangle with vertices at (0, 0), (1, 0), (p, q).";
LeftIntegral::usage = "Takes a function and (p, q) as input and outputs the integral of that function over the triangle with vertices at (0, 0), (1, 0), (p, q), restricting to x < p";
LeftInnerIntegral::usage = "Takes a function and (p, q, x) as input and outputs the integral of that function over the triangle with vertices at (0, 0), (1, 0), (p, q) wrt y, restricting to x < p";
RightIntegral::usage = "Takes a function and (p, q) as input and outputs the integral of that function over the triangle with vertices at (0, 0), (1, 0), (p, q), restricting to x > p";
RightInnerIntegral::usage = "Takes a function and (p, q, x) as input and outputs the integral of that function over the triangle with vertices at (0, 0), (1, 0), (p, q) wrt y, restricting to x > p";


(* Eigenfunctions of the 30-60-90 triangle with vertices at (0, 0), (1/2, 0), (1/2, Sqrt[3]/2) in terms of z and t. *)
Phi301C[z_, t_] := Sin[4z]Sin[2t] - Sin[5z]Sin[t] - Sin[z]Sin[3t]
Phi302C[z_, t_] := Sin[5z]Sin[3t] - Sin[2z]Sin[4t] - Sin[7z]Sin[t]

(* Eigenfunctions of the same triangle in terms of x and y *)
Phi301[x_, y_] := Phi301C[Pi/3 * (2x-1), Pi * (1 - 2y/Sqrt[3])]
Phi302[x_, y_] := Phi302C[Pi/3 * (2x-1), Pi * (1 - 2y/Sqrt[3])]

eval301 = -Laplacian[Phi301[x, y], {x, y}]/Phi301[x, y] // FullSimplify
eval302 = -Laplacian[Phi302[x, y], {x, y}]/Phi302[x, y] // FullSimplify

(* Eigenfunctions of the 45-45-90 triangle with vertices at (0, 0), (1, 0), (0, 1) in terms of x and y *)
Phi451[x_, y_] := Sin[2Pi x]Sin[Pi y] + Sin[Pi x] Sin [2 Pi y]
Phi452[x_, y_] := Sin[3Pi x]Sin[Pi y] - Sin[Pi x] Sin [3 Pi y]

eval451 = -Laplacian[Phi451[x, y], {x, y}]/Phi451[x, y] // FullSimplify
eval452 = -Laplacian[Phi452[x, y], {x, y}]/Phi452[x, y] // FullSimplify

(* Linear transformation to map (0, 0), (p, q), (1, 0) to (1/2, Sqrt[3]/2), (1/2, 0), (0, 0) *)
L30[p_, q_] := TransformationFunction[(\!\(\*
TagBox[GridBox[{
{
RowBox[{
RowBox[{"-", "1"}], "/", "2"}], 
FractionBox[
RowBox[{" ", "p"}], 
RowBox[{"2", "q"}]], 
RowBox[{"1", "/", "2"}]},
{
RowBox[{
RowBox[{"-", 
RowBox[{"Sqrt", "[", "3", "]"}]}], "/", "2"}], 
FractionBox[
RowBox[{
RowBox[{"Sqrt", "[", "3", "]"}], 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", " ", "p"}], ")"}]}], 
RowBox[{"2", "q"}]], 
RowBox[{
RowBox[{"Sqrt", "[", "3", "]"}], "/", "2"}]},
{"0", "0", "1"}
},
AutoDelete->False,
GridBoxDividers->{"Columns" -> {{False}}, "ColumnsIndexed" -> {-2 -> True}, "Rows" -> {{False}}, "RowsIndexed" -> {-2 -> True}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{Automatic}}}],
#& ]\))]

(* Linear transformation to map (0, 0), (p, q), (1, 0) to (1, 0), (0, 0), (0, 1) *)
L45[p_, q_] := TransformationFunction[(\!\(\*
TagBox[GridBox[{
{
RowBox[{"-", "1"}], 
FractionBox[
RowBox[{" ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", "p"}], ")"}]}], "q"], "1"},
{"1", 
RowBox[{"-", 
FractionBox["p", "q"]}], "0"},
{"0", "0", "1"}
},
AutoDelete->False,
GridBoxDividers->{"Columns" -> {{False}}, "ColumnsIndexed" -> {-2 -> True}, "Rows" -> {{False}}, "RowsIndexed" -> {-2 -> True}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{Automatic}}}],
#& ]\))]

(* Defining the integrands for the 30-60-90 triangle *)
A30int[p_, q_, x_, y_] = D[Phi301[L30[p, q][{x, y}][[1]], L30[p, q][{x, y}][[2]]], x]^2 + D[Phi301[L30[p, q][{x, y}][[1]], L30[p, q][{x, y}][[2]]], y]^2 // Simplify;
B30int[p_, q_, x_, y_] = D[Phi301[L30[p, q][{x, y}][[1]], L30[p, q][{x, y}][[2]]], x] *  D[Phi302[L30[p, q][{x, y}][[1]], L30[p, q][{x, y}][[2]]], x]+ D[Phi301[L30[p, q][{x, y}][[1]], L30[p, q][{x, y}][[2]]], y] *  D[Phi302[L30[p, q][{x, y}][[1]], L30[p, q][{x, y}][[2]]], y] // Simplify;
C30int[p_, q_, x_, y_] = D[Phi302[L30[p, q][{x, y}][[1]], L30[p, q][{x, y}][[2]]], x]^2 + D[Phi302[L30[p, q][{x, y}][[1]], L30[p, q][{x, y}][[2]]], y]^2 // Simplify;
D30int[p_, q_, x_, y_] = Phi301[L30[p, q][{x, y}][[1]], L30[p, q][{x, y}][[2]]]^2 // Simplify;
E30int[p_, q_, x_, y_] = Phi301[L30[p, q][{x, y}][[1]], L30[p, q][{x, y}][[2]]] * Phi302[L30[p, q][{x, y}][[1]], L30[p, q][{x, y}][[2]]] // Simplify;
F30int[p_, q_, x_, y_] = Phi302[L30[p, q][{x, y}][[1]], L30[p, q][{x, y}][[2]]]^2 // Simplify;

(* Defining the integrands for the 45-45-90 triangle *)
A45int[p_, q_, x_, y_] = D[Phi451[L45[p, q][{x, y}][[1]], L45[p, q][{x, y}][[2]]], x]^2 + D[Phi451[L45[p, q][{x, y}][[1]], L45[p, q][{x, y}][[2]]], y]^2 // Simplify;
B45int[p_, q_, x_, y_] = D[Phi451[L45[p, q][{x, y}][[1]], L45[p, q][{x, y}][[2]]], x] *  D[Phi452[L45[p, q][{x, y}][[1]], L45[p, q][{x, y}][[2]]], x]+ D[Phi451[L45[p, q][{x, y}][[1]], L45[p, q][{x, y}][[2]]], y] *  D[Phi452[L45[p, q][{x, y}][[1]], L45[p, q][{x, y}][[2]]], y] // Simplify;
C45int[p_, q_, x_, y_] = D[Phi452[L45[p, q][{x, y}][[1]], L45[p, q][{x, y}][[2]]], x]^2 + D[Phi452[L45[p, q][{x, y}][[1]], L45[p, q][{x, y}][[2]]], y]^2 // Simplify;
D45int[p_, q_, x_, y_] = Phi451[L45[p, q][{x, y}][[1]], L45[p, q][{x, y}][[2]]]^2 // Simplify;
E45int[p_, q_, x_, y_] = Phi451[L45[p, q][{x, y}][[1]], L45[p, q][{x, y}][[2]]] * Phi452[L45[p, q][{x, y}][[1]], L45[p, q][{x, y}][[2]]] // Simplify;
F45int[p_, q_, x_, y_] = Phi452[L45[p, q][{x, y}][[1]], L45[p, q][{x, y}][[2]]]^2 // Simplify;

(* Procedure for integration *)
(* We did not directly use these functions in the computations, but this code is included here to illustrate the process we followed for each integration. *)
LeftInnerIntegral[f_, p_, q_, x_] := FullSimplify[Integrate[f[p, q, x, y], {y, 0, q/p * x}, Assumptions->{1>p>0, q>0, x>0}], Assumptions->{1>p>0, q>0, x>0}]
LeftIntegral[f_, p_, q_] := FullSimplify[Integrate[LeftInnerIntegral[f, p, q, x], {x, 0, p}]]
RightInnerIntegral[f_, p_, q_, x_] := FullSimplify[Integrate[f[p, q, x, y], {y, 0, -q/(1-p)*x + q/(1-p)}, Assumptions->{1>p>0, q>0, x>0}], Assumptions->{1>p>0, q>0, x>0}]
RightIntegral[f_, p_, q_] := FullSimplify[Integrate[RightInnerIntegral[f, p, q, x], {x, p, 1}]]
TotalIntegral[f_, p_, q_] := FullSimplify[LeftIntegral[f, p, q] + RightIntegral[f, p, q]]

EndPackage[]


















