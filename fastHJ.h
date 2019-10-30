#ifndef FASTHJ_H 
#define FASTHJ_H

#include <vector> 
#include <array> 

std::vector<int> findIndices(const std::vector<int>& A, const int value);

std::vector<double> limgrad( const std::vector<int> &dims, const double elen, std::vector<double> &ffun, const double dfdx, const int imax); 

bool IsNegative(int i);




#endif 
 

