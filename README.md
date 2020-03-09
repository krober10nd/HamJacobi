# HamJacobi

Solves a Hamilton-Jacobi PDE fast (in seconds) to gradient limit a scalar field defined in 2D or 3D. The input to the solver is packed in column-major order with z being the slowest varying dimension. 

## Compile for MATLAB

This code is designed to be mex'ed using https://github.com/audiofilter/mex-it. From MATLAB enter the following command: 

mex CXXFLAGS="\$CXXFLAGS -std=c++11" FastHJ.cpp

Note: you MAY have to start MATLAB from the terminal (on Linux-like OS) like so:

LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 matlab

...but first try to just use MATLAB without that.

## Compile for Python3 

This code can be called from Python3.x using pybind11 (pip3 install pybind11). Compile the cpp source code as a shared library as such...
c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` FastHJ.cpp -o FastHJ`python3-config --extension-suffix`

# Usage in Python

In python...

import FastHJ 
gradlimited = FastHJ.gradlim(...) 

Type help(FastHJ) to see how to call it.

## Usage in MATLAB

Operate this code from MATLAB by changing the appropriate parts of the code below.

```
dims = [nrows ncols nz]; 
elen = % size of grid cell 
dfdx = % decimal fraction representing smoothness
itmax = % maximum number of iterations to perform 
field = 2D array reshaped in column-major order 

smoothed_field = FastHJ( int32(dims), elen, dfdx, int32(itmax), field);
```
Column-major in MATLAB is readily achieved by using the reshape command like so: 

```
[nrows,ncols]=size(matrix);
vec = reshape(matrix,[numel(matrix),1]); 
```

...and back to the original matrix format like: 

```

matrix = reshape(vec,[nrows,ncols]); 
```

