/* solve the Hamilton-Jacobi equation to smooth a raster field 
 * ref: Per Olof, 
/ kjr, usp, 2019
*/
#include "fastHJ.H"

#include <iostream>
#include <stdio.h>
#include <vector> 
#include <array>
#include <math.h>
#include <cmath>

#include <chrono>


#define EPS 1e-9


// Function to check whether
// the element is odd or not.
bool IsOdd(int i)
{
    return ((i % 2) == 1);
}
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
std::vector<double> limgrad(const std::vector<int>& dims, const double elen, const std::vector<double> &ffun, const double dfdx, const int imax)  {
 
	    int ny = dims[0]; int nx = dims[1]; int nz = dims[2];

        std::vector<int> aset(ny*nx*nz,0);

        double ftol = *(std::min_element(ffun.begin(), ffun.end()))*std::sqrt(EPS);

        std::array<int,5> rm;
        rm.fill(0); 
        std::array<int,5> npos;
        npos.fill(0); 

        for(int iter=0; iter < imax; iter++) {

           //------------------------- find "active" nodes this pass
           auto aidx = findIndices( aset, iter); 

           //------------------------- convergence
	       if(aidx.empty()) break; 

	       for(std::size_t i; i < aidx.size(); i++) {

              //----- map doubly index to singly indexed 
	          int inod = aidx[i];
	          int ipos = 1 + std::floor((inod-1)/dims[2]);
              int jpos = inod - (ipos - 1)*dims[2];  
    
	          // ---- gather indices using 4 edge stencil 
              npos[0] = inod; 
              npos[1] = ipos*ny + jpos;// 
              npos[2] = (ipos-2)*ny + jpos;//nnod of left adj
              npos[3] = (ipos-1)*ny + std::min(jpos+1,ny);//nnod of above adj
              npos[4] = (ipos-1)*ny + std::max(jpos-1,1);//nnod of below adj

              //----- handle boundary vertex adjs.
              //----- iterator that stores the position of last element 
              std :: vector <int>::iterator pend; 
              pend = std::remove_if (npos.begin(),npos.end(), IsOdd);  

              for(std::vector<int>::iterator p=npos.begin(); p != pend; ++p){

                 int nod1 = npos[0];
                 int nod2 = *p; 

                 //----------------- calc. limits about min.-value
                 if (ffun[nod1] > ffun[nod2]) {

                     fun1 = ffun[nod2] + elen * dfdx ;

                     if (ffun[nod1] > fun1+ftol) {
                         ffun[nod1] = fun1;
                         aset[nod1] = iter;
                     }
                 
                 } else {

                     fun2 = ffun[nod1] + elen * dfdx ;

                     if (ffun[nod2] > fun2+ftol) {
                         ffun[nod2] = fun2;
                         aset[nod2] = iter;
                     }
                 }
              }
	       }

	}
	return ffun; 
}


//
int main() {

    std::vector<int> dims = {10,10,10};
    double elen = 1.0;
    std::vector<double> ffun;
    ffun.reserve(1000);
    std::fill(ffun.begin(), ffun.end(),0);
    double dfdx=0.15; 
    int imax=std::sqrt(10000000);
     
    auto begin = std::chrono::steady_clock::now();

    std::vector<double> ffun_s  = limgrad( dims, elen, ffun, dfdx, imax);

    auto end = std::chrono::steady_clock::now();

    std::cout << "Method took " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << " microseconds\n";

    return 0; 
}
