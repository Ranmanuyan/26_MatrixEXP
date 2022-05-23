#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <iostream>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <chrono>


// using namespace boost::multiprecision; 
using namespace std;
using namespace boost::numeric::ublas;
namespace mp = boost::multiprecision;
using Dec = mp::number<mp::cpp_dec_float<0> >;
using Int = mp::cpp_int;

long double boost_factorial(int num){
    Int n = 1;
    for (Int f = num; f>0; --f)
        n *= f;
    long double fac = (1.0/n.convert_to<Dec>()).convert_to<long double>();
    
    return fac;
}


// the number of states is 3n
// Matrix [3n, 3n]
triangular_matrix<double, upper> structureCon(int n, std::vector<std::vector<double> > vect){
    int states = 3*n;
    triangular_matrix<double, upper> Q(states, states);
    // CTMC row_num
    for (int k = 0; k < n-1; ++k){
        Q(3*k+1, 3*k+2) = vect[k][1];
        Q(3*k+1, 3*k+3) = vect[k+1][0];
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

    return Q;
}

triangular_matrix<double, upper> structurePar(int n, std::vector<std::vector<double> > vect){
    std::vector<int> NPS{2,3};
    int states = 3*n + NPS.size();
    triangular_matrix<double, upper> Q(states, states);
    // CTMC row_num
    for (int k = 0; k < n-1; ++k){
        Q(3*k+1, 3*k+2) = vect[k][1];
        Q(3*k+1, 3*k+3) = vect[k+1][0];
        Q(3*k+2, 3*k+4) = vect[k][1];
        Q(3*k+3, 3*k+4) = vect[k+1][0];
    }
    Q(0, 1) = vect[0][0];
    Q(3*n-2,3*n-1) = vect[n-1][1];

    int count =0;
    for (int i = 3*n-1; i < states-1; ++i){
      Q(i,i+1) = vect[NPS[count]][1];
      count++;      
    }
    
    for (int i = 0; i < states; ++i){
        for (int j = i+1; j < states; ++j){
            Q(i,i) += (-1 * Q(i,j));
        }
    }

    return Q;
}

double qCal(triangular_matrix<double, upper> Q){
	double q = 0;
	for (int i = 0; i < Q.size1(); ++i){
		if(Q(i,i) <q){
			q = Q(i,i);
		}
	}
	q = -1*q;

	return q;

}

triangular_matrix<double, upper> Uniformlization(triangular_matrix<double, upper> Q, double q){
	int dem = Q.size1();
	identity_matrix<double> delta(dem);
	
	triangular_matrix<double, upper> Qstar = Q/q + delta;

	return Qstar;

}

std::vector<double> cdfCal(triangular_matrix<double, upper> Qstar, double q){

	std::vector<double>cdfVec;
	int dem = Qstar.size1();
	identity_matrix<double> delta(dem);
    
    
    for(int x =0; x<90;x++){
    	int x1 = int(x); ///10.0;/10.0
    	matrix<double> delta1 = delta;

    	long double cdf=0.0;
    	for (int i = 1; i < 100; ++i){
    		Dec a = pow(q*x1, i);
            long double b = boost_factorial(i);
            long double c = (a*b).convert_to<long double>();

    		delta1 = prod(delta1, Qstar);
            // long double e = exp(-1* q * x1);
    		cdf = cdf + (exp(-1* q * x1)*c)* delta1(0,dem-1);
    	}
    	cdfVec.push_back(cdf);
    }
    return cdfVec;
}



int main () {
    std::vector<std::vector<double> > vect{
        // {0.7, 0.2},{0.8, 0.2},{0.4, 0.8},{0.5, 0.3},{0.1, 0.7},{0.5, 0.3},{0.3, 0.5},{0.4, 0.8},{0.3, 0.8},{0.7, 0.1}
        {0.48,  0.27},  {0.41,  0.52}  ,  {0.28, 0.17} ,   {0.53,0.09},    {0.87, 0.07}, {0.61,  0.49},  { 0.85,  0.03},    {0.92, 0.19},   { 0.8, 0.18},  {0.9, 0.16}
        // {0.4, 0.7},{0.4, 0.3},{0.8, 0.7},{0.4, 0.3}
    };
    int num_job = 10;

    // std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    triangular_matrix<double, upper> Q = structureCon(num_job, vect);
    std::cout << Q << std::endl;
    double q = qCal(Q);
    // std::cout << q << std::endl;
    triangular_matrix<double, upper> Qstar = Uniformlization(Q, q);
    // std::cout << Qstar << std::endl;
    // std::vector<double> cdfVec = cdfCal(Qstar, q);
    // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();


    // for(int i=0; i<cdfVec.size();i++){
    // 	std::cout << cdfVec[i] << std::endl;
    // }
    // std::cout << "Time Costs = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
    return 0;
}