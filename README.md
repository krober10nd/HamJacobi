# HamJacobi

Solves a Hamilton-Jacobi PDE to gradient limit a scalar field defined in 2 or 3D. The input to the solver is packed in column-major order with z being the slowest varying dimension. 

## Compile

This code is designed to be mex'ed using https://github.com/audiofilter/mex-it. From MATLAB enter the following command: 

mex CXXFLAGS="\$CXXFLAGS -std=c++11" FastHJ.cpp

Note: you may have to start MATLAB from the terminal like so:

LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 matlab

## Usage 

Operate this code from MATLAB by changing the appropriate parts of the code below.

```
dims = [ncols nrows nz]; 
elen = % size of grid cell 
dfdx = % decimal fraction representing smoothness
itmax = % maximum number of iterations to perform 
field = 2D array reshaped in column-major order 

smoothed_field = FastHJ( int32(dims), elen, dfdx, int32(itmax), field);
```


