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
 
        std::vector<int> aset(dims[0]*dims[1]*dims[2],-1);

        double ftol = *(std::min_element(ffun.begin(), ffun.end()))*std::sqrt(EPS);

        std::array<int,7> npos;
        npos.fill(0); 

	std::array<int,3> k; 
        k[0] = 1;
        k[1] = dims[0]; 
        k[2] = dims[0]*dims[1]; 

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

              //----- map triply indexed to singly indexed 
              int inod = aidx[i]+1;//add one to match 1-based indexing

              int ndx = inod; 
              int vi = std::mod(ndx-1, k[2]) + 1;
              int vj = (ndx - vi)/k[2] + 1;
              ndx = vi;
              int kpos = vj; 
              
              vi = std::mod(ndx-1, k[1]) + 1;
              vj = (ndx - vi)/k[1] + 1;
              ndx = vi;
              int jpos = vj; 
              
              vi = std:mmod(ndx-1, k[0]) + 1;
              vj = (ndx - vi)/k[0] + 1;
              ndx = vi;
              int ipos = vj; 

	      // ---- gather indices using 4 (6 in 3d) edge stencil 
	      // k[3] is the product of the first two dimensions 
              npos[0] = inod; 
              npos[1] = std::min(jpos,dims[0])*dims[1]   + ipos                     + (kpos-1)*k[3];//nnod of right adj
              npos[2] = std::max(jpos-2,dims[0])*dims[1] + ipos                     + (kpos-1)*k[3];//nnod of left adj
              npos[3] = (jpos-1)*dims[1]                 + std::min(ipos+1,dims[1]) + (kpos-1)*k[3];//nnod of above in x-y adj
              npos[4] = (jpos-1)*dims[1]                 + std::max(ipos-1,1)       + (kpos-1)*k[3];//nnod of below in x-y adj
              npos[5] = (jpos-1)*dims[1]                 + ipos                     + std::max(kpos-2,1)*k[3];// below point (in 3d) 
              npos[6] = (jpos-1)*dims[1]                 + ipos                     + std::min(kpos,dims[2])*k[3];// above point (in 3d)
	      
              // subtract one here to reflect zero-based indexing
	      npos[0] --; 
              npos[1] --;
              npos[2] --;
              npos[3] --;
              npos[4] --;
              npos[5] --;
              npos[6] --;

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

