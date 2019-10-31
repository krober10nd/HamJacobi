# HamJacobi

Solves a Hamilton-Jacobi PDE to gradient limit a scalar field defined in 2 or 3D. The input to the solver is packed in row-major order with z being the slowest varying dimension. 

## Compile

This code is designed to be mex'ed using https://github.com/audiofilter/mex-it

mex CXXFLAGS="\$CXXFLAGS -std=c++11" FastHJ.cpp

## Usage 

Operate this code from MATLAB by changing the appropriate parts of the code below. 

dims = [ny nx nz];
elen = 1.0  % size of grid cell 
dfdx = 0.15 % decimal fraction representing smoothness
itmax = 100 % maximum number of iterations to perform 
field = 2D array reshaped in column-major order 
smoothed_field = FastHJ( int32(dims), elen, dfdx, int32(400), field);


