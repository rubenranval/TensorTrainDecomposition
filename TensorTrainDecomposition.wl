(* ::Package:: *)

(* Title: TensorTrainDecomposition *)
(* Author: Ruben Ranval *)
(* Summary:
     Decompose a numerical array of any rank into a Tensor Train
     (also known as Matrix Product State) representation. *)
(* Version: 1.0.1 *)

(* Note: Also available as an isntant add-on function on https://resources.wolframcloud.com/FunctionRepository/resources/TensorTrainDecomposition/ *)

BeginPackage["TensorTrainDecomposition`"];


(* ::Section:: *)
(*Usage*)


TensorTrainDecomposition::usage =
  "TensorTrainDecomposition[tensor] decomposes tensor into a list of low-rank cores.";


Begin["`Private`"];


(* ::Section:: *)
(*Options*)


$defaultMaxBondDimension = Infinity;
$defaultTolerance = 0;
$defaultMethod = "SVD";

Options[TensorTrainDecomposition] = {
   "MaxBondDimension" -> $defaultMaxBondDimension,
   Tolerance -> $defaultTolerance,
   Method -> $defaultMethod
};


(* ::Section:: *)
(*Messages*)


TensorTrainDecomposition::notarr  = "Input must be a numeric array.";
TensorTrainDecomposition::badchi  = "MaxBondDimension must be a positive integer or Infinity.";
TensorTrainDecomposition::badeps  = "Tolerance must be a non-negative number.";
TensorTrainDecomposition::badmeth = "Method must be \"SVD\" or \"QR\".";


(* ::Section:: *)
(*Definition*)


TensorTrainDecomposition[tensor_, opts : OptionsPattern[]] :=
  Module[
   {\[Chi]max = OptionValue["MaxBondDimension"], \[CurlyEpsilon] = OptionValue[Tolerance],
    method = OptionValue[Method], dims, d, residual, cores = {},
    r = 1, nk, restDims, U, \[CapitalSigma], V, q, Q, R, rNew, \[Sigma]},

   If[! ArrayQ[tensor, _, NumericQ],
     Message[TensorTrainDecomposition::notarr]; Return[$Failed]];
   If[! (\[Chi]max === Infinity || (IntegerQ[\[Chi]max] && \[Chi]max >= 1)),
     Message[TensorTrainDecomposition::badchi]; Return[$Failed]];
   If[! NumericQ[\[CurlyEpsilon]] || Negative[\[CurlyEpsilon]],
     Message[TensorTrainDecomposition::badeps]; Return[$Failed]];
   If[! MemberQ[{"SVD", "QR"}, method],
     Message[TensorTrainDecomposition::badmeth]; Return[$Failed]];

   dims = Dimensions@tensor;
   d = Length@dims;
   residual = tensor;

   Do[
     nk = dims[[k]];
     restDims = Times @@ Drop[dims, k];
     residual = ArrayReshape[residual, {r*nk, restDims}];

     Switch[method,
       "SVD",
       {U, \[CapitalSigma], V} = SingularValueDecomposition@residual;
       \[Sigma] = Diagonal@\[CapitalSigma];
       rNew = Clip[
          Count[Reverse@Accumulate[Reverse[\[Sigma]^2]], x_ /; x >= \[CurlyEpsilon]^2],
          {1, Min[\[Chi]max, Length@\[Sigma]]}];
       AppendTo[cores, ArrayReshape[U[[All, 1 ;; rNew]], {r, nk, rNew}]];
       (* fixed in 1.0.1: replaced Transpose with ConjugateTranspose for complex numbers *)
       residual = \[CapitalSigma][[1 ;; rNew, 1 ;; rNew]] . ConjugateTranspose[V[[All, 1 ;; rNew]]],

       "QR",
       {q, R} = QRDecomposition[residual];
       Q = ConjugateTranspose[q];
       rNew = Dimensions[R][[1]];
       AppendTo[cores, ArrayReshape[Q, {r, nk, rNew}]];
       residual = R];

     r = rNew, {k, 1, d - 1}];

   AppendTo[cores, ArrayReshape[residual, {r, dims[[d]], 1}]];
   cores];
   


End[];
EndPackage[];
