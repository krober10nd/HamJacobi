/* solve the Hamilton-Jacobi equation to smooth a raster field 
/ kjr, usp, 2019
*/
#include "fastHJ.H"

#include <stdio.h>
#include <vector> 
#include <array>
#include <math.h>

#define EPS 1e-9


// find index in linear time 
std::vector<int> findIndices(const std::vector<int>& A, const int value) {
	// try calling reserve here with some estimated size amount
	std::vector<int> B; 
	for(std::size_t i=0; i<A.size(); i++) {
	   if(A[i] == value) {
	     B.push_back(i); 
	   }
	}
	return B; 
}

// solve the Hamilton-Jacobi equation
std::vector<double> ffun limgrad( const std::vector<int>& dims, const double elen, const std::vector<double> &ffun, const double dfdx, const int imax)  {
 
	int ny = dims[0]; int nx = dims[1]; int nz = dims[2];

        std::array<int,ny*nz*nx> aset;
        aset.fill(0);

        double ftol = std::min_element(ffun.begin(), ffun.end())*std::sqrt(EPS);

        auto rm = array<double,5>;
        rm.fill(0);

        for(int iter=0; iter < itmax; iter++) {

           //------------------------- find "active" nodes this pass
           auto aidx = findIndices( aset, iter); 

	   if(aidx.empty) break; 

	   for(std:size_t :: i; i < aidx.size(); i++) {
	      
              //----- map doubly index to singly indexed 
	      int inod = aidx[i];
	      int ipos = 1 + std::floor((inod-1)/dims[3];
              int jpos = inod - (ipos - 1)*dims[3];  

	      // ---- gather indices using 4 edge stencil 
	      //
	      //
	   };

	}
	return ffun; 
}


//
int main {
   

return 0; 
}
