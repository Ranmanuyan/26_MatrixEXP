#include <iterator>
#include <string>
#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <chrono>

#include <boost/assign/std/vector.hpp> 
#include <boost/random.hpp>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "expm.h"

// using namespace boost::multiprecision; 
using namespace std;
using namespace boost::numeric::ublas;



// the number of states is 3n
// Matrix [3n, 3n]
triangular_matrix<double, upper> structureCon(int n, std::vector<std::vector<double> > vect){
    int states = 3*n;
    triangular_matrix<double, upper> Q(states, states);
    // CTMC row_num
    for (int k = 0; k < n-1; ++k){
        Q(3*k+1, 3*k+2) = vect[k+1][0];
        Q(3*k+1, 3*k+3) = vect[k][1];
        Q(3*k+2, 3*k+4) = vect[k][1];
        Q(3*k+3, 3*k+4) = vect[k+1][0];
    }
    Q(0, 1) = vect[0][0];
    Q(3*n-2, 3*n-1) = vect[n-1][1];
    
    for (int i = 0; i < states; ++i){
        for (int j = i+1; j < states; ++j){
            Q(i,i) += (-1 * Q(i,j));
        }
    }
    triangular_matrix<double, upper> C = project(Q, range (0, states-1), range (0, states-1));    // the submatrix of A specified by the two index ranges r1 and r2

    return C;
}



std::vector<double> cdfCal(triangular_matrix<double, upper> Qstar, boost::numeric::ublas::vector<double> a){
  std::vector<double>cdfVec;
  std::vector<double>xVec;
  int dem = Qstar.size1();

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  
  for (int i = 1; i < 150; ++i){
    // xVec.push_back(i);
    double theta = double(i);
    triangular_matrix<double, upper> mat = theta * Qstar;
    // double cdf = 1- sum(prod(a, expm_pad(mat)));
    double cdf = 1- sum(row(expm_pad(mat), 0));

    cdfVec.push_back(double(cdf));
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time Costs = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[milliseconds]" << std::endl;

  return cdfVec;
}



int main () {
    std::vector<std::vector<double> > vect{
        {0.48, 0.27},{0.41, 0.52},{0.28, 0.17},{0.53, 0.09},{ 0.87, 0.07},{0.61, 0.49},{0.85,0.03},{0.92,0.19 },{ 0.8, 0.18},{0.9, 0.16} 
        // {0.4, 0.7},{0.5, 0.3},{0.8, 0.6},{0.2, 0.1}
    };
    int num_job = 10;


    

    triangular_matrix<double, upper> Q = structureCon(num_job, vect);
    // std::cout << Q << std::endl;
    boost::numeric::ublas::vector<double> alpha(Q.size1(),0);
	alpha[0] = 1;
    // std::cout << alpha << std::endl;

    // std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::vector<double> cdfVec = cdfCal(Q, alpha);
    // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();


    for(int i=0; i<cdfVec.size();i++){
    	std::cout << cdfVec[i] << std::endl;
    }
    // std::cout << "Time Costs = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[milliseconds]" << std::endl;
    return 0;
}

// compile command
// g++ -O3 -std=c++11 -g -I /usr/local/boost_1_73_0 expm_sample.cpp -o test