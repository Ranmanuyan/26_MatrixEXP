/*
* @Author: Lei Liu
* @Date:   2021-11-10 16:14:59
* @Last Modified by:   Lei Liu
* @Last Modified time: 2021-11-15 17:42:05
*/

// https://stackoverflow.com/questions/35656237/dynamic-eigen-vectors-in-boostodeint


#include <iostream>
#include <Eigen/Core>
#include <cstdlib>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

using namespace boost::numeric::odeint;
typedef Eigen::VectorXd state_type;

state_type x;
Eigen::MatrixXd A;

void ODE_function (const state_type &x, state_type &dxdt, double){
    dxdt = A * x;
}

void write_states(const state_type &x, const double t){
    std::cout << t << "\t";
    for (int i = 0; i < x.size(); i++)
    {
        std::cout << *(x.data()+i) << "\t";
    }
    std::cout << std::endl;
}

int main(){
    int nr_of_states = 3;

    // std::cout << "How many (random) states would you like to simulate?: ";
    // std::cin >> nr_of_states;
    // std::cout << std::endl;

    x = state_type(nr_of_states);
    A = Eigen::MatrixXd(nr_of_states, nr_of_states);

    srand(365);

    for (int i = 0; i < A.size(); i++){
        *(A.data()+i) = ((double)rand()/(double)RAND_MAX);
    }

    for (int i = 0; i < x.size(); i++){
        *(x.data()) = 0;
    }


    std::cout << x <<std::endl;

    typedef runge_kutta_dopri5<state_type, double, state_type, double, vector_space_algebra> stepper;
    integrate_adaptive(stepper(), ODE_function, x, 0.0, 25.0, 1.0,write_states); //

    std::cout << "final state vector: " << std::endl << x << std::endl;

    return 0;
}


// g++ -O0 -std=c++11 -I /usr/local/boost_1_73_0 -g  odeIntMatrix.cpp -o test