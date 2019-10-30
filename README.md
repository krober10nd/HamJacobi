# HamJacobi

Solves a Hamilton-Jacobi PDE to gradient limit a scalar field defined in 2 or 3D. The input to the solver is packed in row-major order with z being the slowest varying dimension. 

## Compile

g++ -O3 -std=c++11 fastHJ.cpp -I fastHJ.h -o fastHJ.x
