/*
* @Author: Lei Liu
* @Date:   2021-11-10 16:14:59
* @Last Modified by:   Lei Liu
* @Last Modified time: 2021-11-16 10:06:33
*/

// https://stackoverflow.com/questions/35656237/dynamic-eigen-vectors-in-boostodeint
// http://boccelliengineering.altervista.org/junk/boost_integration/boost_odeint.html
// https://stackoverflow.com/questions/53855922/pass-parameter-to-boost-integrator

#include <iostream>
#include <Eigen/Core>
#include <cstdlib>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

using namespace boost::numeric::odeint;
typedef Eigen::VectorXd state_type;
// typedef boost::array< double , 4 > state_type;
state_type x(4);



// void ODE_function (const state_type &x, state_type &dxdt, double){
//     // std::cout << A << std::endl;
//     dxdt = A*x;
// }

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
        dxdt = A*x;
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

    typedef runge_kutta_dopri5<state_type, double, state_type, double, vector_space_algebra> stepper;
    // Stepper stepper, System system, State & start_state, Time start_time, Time end_time, Time dt(not required), Observer observer(write_states,not required)
    ode  ode1(A);
    integrate_adaptive(stepper(), ode1, x, 0.0, 2.0, 1.0); //write_states  

    std::cout << "final state vector: " << std::endl << x(0) << std::endl;

    return 0;
}


// g++ -O0 -std=c++11 -I /usr/local/boost_1_73_0 -g  odeIntMatrix2.cpp -o test