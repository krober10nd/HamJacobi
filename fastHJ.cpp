/* solve the Hamilton-Jacobi equation to smooth a raster field 
   Persson, PO. Engineering with Computers (2006) 22: 95. https://doi.org/10.1007/s00366-006-0014-1 
   kjr, usp, 2019
*/
#include "fastHJ.h"

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
std::vector<double> limgrad(const std::vector<int>& dims, const double elen, std::vector<double> &ffun, const double dfdx, const int imax)  {
 
        int nrows = dims[0]; int ncols = dims[1]; int nz = dims[2];

        std::vector<int> aset(nrows*ncols*nz,-1);

        double ftol = *(std::min_element(ffun.begin(), ffun.end()))*std::sqrt(EPS);

        std::array<int,5> rm;
        rm.fill(0); 
        std::array<int,5> npos;
        npos.fill(0); 

        for(int iter=0; iter < imax; iter++) {

           //------------------------- find "active" nodes this pass
           auto aidx = findIndices( aset, iter-1);

           //------------------------- convergence
	       if(aidx.empty()) break; 

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
                 if (ffun[nod1] > ffun[nod2]) {

                     double fun1 = ffun[nod2] + elen * dfdx ;

                     if (ffun[nod1] > fun1+ftol) {
                         ffun[nod1] = fun1;
                         aset[nod1] = iter;
                     }
                 
                 } else {

                     double fun2 = ffun[nod1] + elen * dfdx ;

                     if (ffun[nod2] > fun2+ftol) {
                         ffun[nod2] = fun2;
                         aset[nod2] = iter;
                     }
                 }
              }
	       }
           std::cout << iter << std::endl;
	}
	return ffun; 
}


//
int main() {

    std::vector<int> dims = {13601,2801,1};
    double elen = 1.25;
    double dfdx=0.15;

    std::vector<double> ffun;

    std::ifstream meshSizes;
    meshSizes.open("/home/keith/HamJacobi/meshsize.txt");
    double d;
    while (meshSizes >> d) //ifstream does text->double conversion
       ffun.push_back(d); // add to vector
    meshSizes.close();

    int imax=std::sqrt(ffun.size());
     
    auto begin = std::chrono::steady_clock::now();

    std::vector<double> ffun_s  = limgrad( dims, elen, ffun, dfdx, imax);

    auto end = std::chrono::steady_clock::now();

    std::cout << "Method took " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << " microseconds\n";

    std::ofstream myfile;
    myfile.open ("meshsize_smoothed.txt");
    for(std::size_t i=0; i<ffun_s.size(); i++)
        myfile << ffun_s[i] << std::endl; 
    myfile.close();

    return 0; 
}
