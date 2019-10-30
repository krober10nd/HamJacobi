#ifndef FASTHJ_H 
#define FASTHJ_H

#include <vector> 
#include <array> 

std::vector<int> findIndices(const std::vector<int>& A, const int value);

std::vector<double> limgrad( const std::array<int,3> &dims, const double elen, const std::vector<double> &ffun, const double dfdx, const int imax); 

#endif 
 

