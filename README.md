# TensorTrainDecomposition

Decompose a numerical array of any rank into a Tensor Train (also known as a Matrix Product State) representation.

> Also published in the [Wolfram Function Repository](https://resources.wolframcloud.com/FunctionRepository/resources/TensorTrainDecomposition/).

## Usage

```mathematica
TensorTrainDecomposition[tensor]
```

Decomposes `tensor` into a list of low-rank cores.

## Details & Options

- `TensorTrainDecomposition` decomposes a numerical array of any rank into a Tensor Train (Matrix Product State) representation.
- The result is a list of rank-3 arrays $\{A_1, A_2, \ldots, A_n\}$ where $n$ is the number of dimensions of the input tensor.
- The bond dimension refers to the size of the internal contracted indices connecting adjacent tensor cores. It effectively dictates the rank and compression level of the Tensor Train.
- Each core $A_k$ has dimensions $\{\chi_{k-1},\, d_k,\, \chi_k\}$ where $d_k$ is the size of the $k^\text{th}$ dimension of the input tensor, and $\chi_k$ is the bond dimension between cores $k$ and $k+1$.
- The first core has $\chi_0 = 1$ and the last core has $\chi_n = 1$.

### Options

| Option | Default | Description |
| --- | --- | --- |
| `Method` | `"SVD"` | decomposition method |
| `"MaxBondDimension"` | `Infinity` | maximum bond dimension $\chi$ at each core |
| `Tolerance` | `0` | truncation threshold $\varepsilon$ |

### Behavior

- With `Method -> "QR"`, the function ignores `"MaxBondDimension"` and `Tolerance` and returns the exact decomposition. `QRDecomposition` does not provide singular values, so truncation is not available.
- With the default settings, `TensorTrainDecomposition` returns an exact decomposition with no truncation.
- Setting only `"MaxBondDimension"` controls memory: bond dimensions are capped at the specified value.
- Setting only `Tolerance` controls accuracy: bond dimensions adapt automatically to meet the error target.
- With `Method -> "SVD"`, singular values are discarded until the sum of the squared discarded values reaches $\varepsilon^2$, where $\varepsilon$ is the value specified by `Tolerance`.

For more mathematical background on this algorithm, see the Wikipedia articles on [Tensor Network](https://en.wikipedia.org/wiki/Tensor_network) and [Matrix Product State](https://en.wikipedia.org/wiki/Matrix_product_state).

## Examples

### Basic Examples

Decompose a simple 2×3×4 tensor into a Tensor Train with exact mathematical precision:

```mathematica
tensor = RandomReal[1, {2, 3, 4}];

TensorTrainDecomposition[tensor]
```

```
{{{{-0.611543, 0.791211}, {-0.791211, -0.611543}}},
 {{{-0.569533, -0.11767, 0.15997, 0.630821},
   {-0.556592, -0.111827, -0.755852, -0.175653},
   {-0.600728, 0.284011, 0.545842, -0.377351}},
  {{-0.0130681, -0.836542, 0.28009, 0.0751993},
   {0.0587628, 0.376598, 0.10921, 0.230513},
   {-0.0365915, -0.226623, 0.121557, -0.608298}}},
 {{{1.06126}, {1.33467}, {1.42006}, {1.40377}},
  {{0.244649}, {0.108796}, {-0.742054}, {0.462269}},
  {{-0.46862}, {0.358826}, {-0.0574899}, {0.0712752}},
  {{-0.110017}, {-0.172284}, {0.0566347}, {0.189685}}}}
```

Inspect the dimensions of the underlying cores forming the Tensor Train:

```mathematica
Dimensions /@ %
```

```
{{1, 2, 2}, {2, 3, 4}, {4, 4, 1}}
```

Note how the cores chain together: the trailing bond dimension of each core matches the leading bond dimension of the next (in this case 2 and 4), and the outer bonds are both 1.

### Options

#### MaxBondDimension

Force a maximum bond dimension (rank limit) to compress a large tensor and save memory:

```mathematica
largeTensor = RandomReal[{-1, 1}, {4, 4, 4, 4}];

TensorTrainDecomposition[largeTensor, "MaxBondDimension" -> 2]
```

```
{{{{0.670903, 0.172305}, {-0.413005, 0.886109}, {-0.482858, -0.406204}, {0.382315, 0.141843}}},
 {{{-0.0772469, 0.368848}, {0.276403, 0.0562738}, {-0.00902274, -0.573394}, {0.616954, -0.373717}},
  {{0.351203, 0.405863}, {0.597212, 0.0667429}, {-0.178279, 0.0973214}, {0.158528, 0.462265}}},
 {{{0.198657, -0.69337}, {-0.218896, 0.293781}, {0.448916, -0.283257}, {0.498525, 0.588681}},
  {{0.502755, 0.0013303}, {0.298684, 0.0497181}, {0.142256, 0.00837742}, {-0.316789, -0.0600596}}},
 {{{0.68475}, {1.69375}, {-1.04008}, {3.00051}}, {{-1.0705}, {-0.3975}, {1.7939}, {1.09051}}}}
```

Notice how the internal connecting ranks are strictly bounded to 2:

```mathematica
Dimensions /@ %
```

```
{{1, 4, 2}, {2, 4, 2}, {2, 4, 2}, {2, 4, 1}}
```

#### Tolerance

Discard singular values dynamically based on a truncation error threshold, naturally isolating the true signal from noise:

```mathematica
noisyTensor = RandomReal[1, {5, 5, 5}];

TensorTrainDecomposition[noisyTensor, Tolerance -> 0.1]
```

```
{{{{-0.393123, -0.707127, -0.545606, 0.217018, -0.0253512},
    {-0.432826, 0.384423, -0.259551, -0.094373, 0.767208},
    {-0.445962, -0.128103, 0.695517, 0.536646, 0.113905},
    {-0.49952, -0.170061, 0.280364, -0.777181, -0.197347},
    {-0.457945, 0.553949, -0.269445, 0.22803, -0.599024}}},
 (* ... *)
 {{{3.23721}, {2.4514}, {2.96524}, {2.84792}, {2.12999}}, ...}}
```

The bond dimensions adapt automatically so that the squared discarded singular values stay within the error target.

#### Method

Use `"QR"` instead of `"SVD"` to perform an exact transformation without truncation (note that `"MaxBondDimension"` and `Tolerance` are ignored):

```mathematica
TensorTrainDecomposition[RandomReal[1, {2, 3, 4}], Method -> "QR"]
```

```
{{{{-0.818112, -0.57506}, {-0.57506, 0.818112}}},
 {{{-0.771848, -0.0757732, 0.220253, -0.247601},
   {-0.161964, -0.469426, -0.73878, -0.414458},
   {-0.370874, -0.0315051, -0.371841, 0.84665}},
  {{0., 0.381762, -0.296661, -0.0361467},
   {-0.0253851, 0.748124, -0.374347, -0.195012},
   {0.489721, -0.259758, -0.198201, 0.103758}}},
 {{{1.42218}, {1.03923}, {0.390375}, {0.46433}}, {{0.}, {0.606159}, {0.930472}, {0.691327}},
  {{0.}, {0.}, {1.09169}, {0.771255}}, {{0.}, {0.}, {0.}, {0.440308}}}}
```

The upper-triangular pattern of zeros in the final core reflects the `R` factor produced by the QR factorization.

### Possible Issues

`TensorTrainDecomposition` respects exact inputs (like integers or fractions), which causes the underlying `SingularValueDecomposition` to attempt symbolic root-finding. For large matrices, this scales exponentially:

```mathematica
exactTensor = RandomInteger[1, {2, 3, 4}];

AbsoluteTiming[
  ttSlow = TensorTrainDecomposition@exactTensor;
]
```

```
{5.68876, Null}
```

Use `N` to convert exact inputs to machine precision for much faster execution:

```mathematica
AbsoluteTiming[
  ttFast = TensorTrainDecomposition@N@exactTensor;
]
```

```
{0.0012639, Null}
```

### Neat Examples

Reconstruct a non-separable pattern from tensor trains of increasing bond dimension to watch a low-rank approximation sharpen into the exact array:

```mathematica
n = 120;
c = (n + 1)/2;
pattern = Table[Sin[Sqrt[(i - c)^2 + (j - c)^2]/4.], {i, n}, {j, n}];
ranks = {1, 3, 5, 20, 30};
recons = Table[
   ResourceFunction["TensorTrainContract"][
    TensorTrainDecomposition[pattern, "MaxBondDimension" -> r]], {r, ranks}];
Row[MapThread[
   ArrayPlot[#1, ColorFunction -> "DeepSeaColors", Frame -> False,
      ImageSize -> 150, PlotLabel -> #2] &,
   {Append[recons, pattern],
    Append[("\[Chi] = " <> ToString[#] &) /@ ranks, "exact"]}]]
```

<img src="/img/reconstruction1.png" width="480">

## Installation

Load the function directly from the Wolfram Function Repository:

```mathematica
ResourceFunction["TensorTrainDecomposition"][tensor]
```

Or load the source from this repository into a session:

```mathematica
Get["path/to/TensorTrainDecomposition.wl"]
```

## License
MIT License

Copyright (c) 2026 Ruben Ranval

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
