/* solve the Hamilton-Jacobi equation to smooth a raster field 
   Persson, PO. Engineering with Computers (2006) 22: 95. https://doi.org/10.1007/s00366-006-0014-1 
   kjr, usp, 2019
*/

#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector> 
#include <array>
#include <math.h>
#include <cmath>
#include <fstream>

#include <chrono>

#define EPS 1e-9

bool IsNegative(int i) {return (i < 0);}

// find indices in linear time where A==value
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
std::vector<double> limgrad(const std::vector<int>& dims, const double& elen, const double& dfdx, const int& imax, const std::vector<double> &ffun)  {
 
        int nrows = dims[0]; int ncols = dims[1]; int nz = dims[2];

        std::vector<int> aset(nrows*ncols*nz,-1);

        double ftol = *(std::min_element(ffun.begin(), ffun.end()))*std::sqrt(EPS);

        std::array<int,5> rm;
        rm.fill(0); 
        std::array<int,5> npos;
        npos.fill(0); 

	// allocate output 
	std::vector<double> ffun_s; 
	ffun_s.resize(ffun.size()); 
	ffun_s = ffun; 

        for(int iter=0; iter < imax; iter++) {

           //------------------------- find "active" nodes this pass
           auto aidx = findIndices( aset, iter-1);

           //------------------------- convergence
	   if(aidx.empty()) {
             std::cout << "INFO: Converged in " << iter << " iterations." << std::endl; 
             break; 
	   }

           for(std::size_t i=0; i < aidx.size(); i++) {

              //----- map doubly index to singly indexed 
              int inod = aidx[i];
              int ipos = std::floor((inod-1)/ncols);
              int jpos = inod - ipos*ncols;
    
	          // ---- gather indices using 4 edge stencil 
              npos[0] = inod; 
              npos[1] = (ipos+1)*ncols + jpos;//
              npos[2] = (ipos-1)*ncols + jpos;//nnod of left adj
              npos[3] = ipos*ncols + std::min(jpos+1,nrows);//nnod of above adj
              npos[4] = ipos*ncols + std::max(jpos-1,0);//nnod of below adj

              //----- handle boundary vertex adjs.
              //----- iterator that stores the position of last element 
              auto pend = std::remove_if(npos.begin()+1,npos.end(), IsNegative);

              for(auto p=npos.begin()+1; p!=pend; p++){

                 int nod1 = npos[0];
                 int nod2 = *p; 

                 //----------------- calc. limits about min.-value
                 if (ffun_s[nod1] > ffun_s[nod2]) {

                     double fun1 = ffun_s[nod2] + elen * dfdx ;

                     if (ffun_s[nod1] > fun1+ftol) {
                         ffun_s[nod1] = fun1;
                         aset[nod1] = iter;
                     }
                 
                 } else {

                     double fun2 = ffun_s[nod1] + elen * dfdx ;

                     if (ffun_s[nod2] > fun2+ftol) {
                         ffun_s[nod2] = fun2;
                         aset[nod2] = iter;
                     }
                 }
              }
	       }
	}
	return ffun_s; 
}

// mex it up
// first args are inputs, last arg is output 
void mex_function(const std::vector<int>& dims, const double& elen, const double& dfdx, const int& imax, const std::vector<double>& ffun, std::vector<double>& ffun_s) {

    // this is how you call the function from matlab 
    ffun_s = limgrad( dims, elen, dfdx, imax, ffun);

}

#include "mex-it.h"

