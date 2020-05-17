# Fourier Transform in C#

| Algorithm | Limitations                                               |
|:--------- |:-----------                                               |
| DFT       | 1D, 2D, and 3D transformations.                           |
| FFT       | supports N dimensions, each with power of 2 side lengths. |

## Usage:

Compute a simple FFT on data with length power of 2:
```
var data = new Complex[1024];
SetComplexData( data );
data.Fft();
```

Compute a multi-dimensional FFT on cube with sides of length of power of 2:
```
var log2n = 8;
var n = 1 << log2N;
var data = new Complex[n*n*n];
SetComplexData( data );
data.MultiFft( log2n, log2n, log2n );
```

There is an example "clouds" animation (.net core forms app) on a 256 x 256 x 256
using 1/f^2 noise computation in the TrentTobler.Examples.FourierTransform project.
