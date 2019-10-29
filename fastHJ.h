#ifndef FASTHJ_H 
#define FASTHJ_H

#include <vector> 
#include <array> 

double* ffun limgrad( const std::array<int,3> &dims, const double elen, const std::vector<double> &ffun, const double dfdx, const int imax); 

#endif 
 

