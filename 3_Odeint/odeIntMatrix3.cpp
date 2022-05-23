/*
* @Author: Lei Liu
* @Date:   2021-11-10 16:14:59
* @Last Modified by:   Lei Liu
* @Last Modified time: 2021-11-16 11:11:08
*/

// https://stackoverflow.com/questions/35656237/dynamic-eigen-vectors-in-boostodeint
// http://boccelliengineering.altervista.org/junk/boost_integration/boost_odeint.html
// https://stackoverflow.com/questions/53855922/pass-parameter-to-boost-integrator


// I slightly changed the ode function, so that I can calculate the first row of the e^(A) matrix
// add this row is the cdf 
// This I think can decrease the calculation burden largely


#include <iostream>
#include <Eigen/Core>
#include <cstdlib>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>
#include <chrono>
#include <Eigen/unsupported/Eigen/MatrixFunctions>
using namespace boost::numeric::odeint;
typedef Eigen::RowVectorXd state_type;
// state_type x(4);


void write_states(const state_type &x, const double t){
    std::cout << t << "\t";
    for (int i = 0; i < x.size(); i++){
        std::cout << *(x.data()+i) << "\t";
    }
    std::cout << std::endl;
}

struct ode{
    const Eigen::MatrixXd A;
    ode( Eigen::MatrixXd A1 ): A(A1)  {}

    void operator()( state_type const& x , state_type& dxdt , double ) const {
        dxdt = x*A;
    }
};



int main(){
    state_type x(4);
    x(0)=1;
    x(1)=0;
    x(2)=0;
    x(3)=0;

    Eigen::MatrixXd A(4,4);
    A(0,0)=-0.128;
    A(0,1)=0.044;   
    A(0,2)=0.012;
    A(0,3)=0.03;
    A(1,0)=0.061;
    A(1,1)=-0.061;  
    A(1,2)=0.0;
    A(1,3)=0.0;
    A(2,0)=0.051;
    A(2,1)=0.0;
    A(2,2)=-0.051;
    A(2,3)=  0.0;
    A(3,0)=0.0;
    A(3,1)=0.0;
    A(3,2)=0.088;
    A(3,3)= -0.088;


    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    typedef runge_kutta_dopri5<state_type, double, state_type, double, vector_space_algebra> stepper;
    ode  ode1(A);
    integrate_adaptive(stepper(), ode1, x, 0.0, 500.0, 500.0,write_states);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    std::cout << "final state vector: " << std::endl << x.sum() << std::endl;

    std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();
    Eigen::MatrixXd mat = (500*A).exp();  
    std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end1 - begin1).count() << "[µs]" << std::endl;


    return 0;
}


// g++ -O3 -std=c++11 -I /usr/local/boost_1_73_0 -g  odeIntMatrix2.cpp -o test