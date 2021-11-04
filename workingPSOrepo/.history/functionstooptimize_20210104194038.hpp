
//#include <iostream>
//#include <vector>
#include <math.h>
#include <unistd.h>

//#define PI 3.141592
# variable_datatype double


#define PI ((variable_datatype)std::acos(-1))


double booth(variable_datatype x, variable_datatype y)
{
    variable_datatype  op1 = x + 2*y - 7,
            op2 = 2*x + y - 5;

    return op1*op1 + op2*op2; //
}

double himmelblau(variable_datatype x, variable_datatype y)
{
    variable_datatype  op1 = x*x + y - 11,
            op2 = x + y*y - 7;

    return op1*op1 + op2*op2;
}

double easom(variable_datatype x, variable_datatype y)
{
    variable_datatype  opcos = -std::cos(x) * std::cos(y),
            xpart = (x-PI)*(x-PI),
            ypart = (y-PI)*(y-PI),
            opexp = std::exp( -( xpart + ypart ) );

    return opcos*opexp;
}

double holdertable(variable_datatype x, variable_datatype y)
{
    variable_datatype  optrigo = std::sin(x) * std::cos(y),
            opexp   = std::exp( std::abs( 1 - std::sqrt( x*x + y*y )/PI ) );

    return std::abs( optrigo*opexp );
}

double ackley(variable_datatype x, variable_datatype y)
{
    variable_datatype  op1 = 20 * std::exp( -0.2 * std::sqrt( 0.5 * (x*x + y*y) ) ),
            op2 = std::exp( 0.5 * ( std::cos( 2*PI*x ) + std::cos( 2*PI*y ) ) );

    return -op1 - op2 + std::exp(1) + 20;
}


////////////////////////////////////////////////
//MY_FUNCITONS BELOW///////////////////////////

double VERY_Complicated(variable_datatype x, variable_datatype y){
    variable_datatype tmp1 = x*x + y*y;
    variable_datatype tmp2 = cos(x / 1) * cos(y / std::sqrt(2));
    variable_datatype temp = std::sin(20. - 20. * std::exp(-0.2 * std::sqrt(0.5 * tmp1)) - std::exp(0.5 * tmp2));
    //std::chrono::sleep_until(std::chrono::system_clock::now() + std::chrono::nanoseconds(100));
    unsigned int microsecond = 1; // 1 sec = 1000000
    usleep(microsecond);//sleeps for 3 second
    return (tmp1 / 4000. - tmp2 + 1) * temp + 23.;
}

double Sphere_f(variable_datatype x, variable_datatype y){
    variable_datatype ret = x*x + y*y;
    return ret;
}

double Schwefel_f(variable_datatype x, variable_datatype y){
    variable_datatype ret = 0;
    variable_datatype temp = x+y;
    ret = x*x + temp*temp;
    return ret;
}

double Rosenbrock_f(variable_datatype x, variable_datatype y){
    variable_datatype ret = 100. * (y + x*x) + (x - 1)*(x - 1);
    return ret;
}

double Schwefel2_f(variable_datatype x, variable_datatype y){
    variable_datatype ret = x * sin(std::sqrt(std::abs(x)) +  y * sin(std::sqrt(std::abs(y))));
    return -ret;
}

double Rastrigin_f(variable_datatype x, variable_datatype y){
    variable_datatype ret = x*x - 10 * cos(2 * PI * x) + y*y - 10 * cos(2 * PI * y);
    return 20 + ret;
    //return 10 + ret;
}

double Ackley_f(variable_datatype x, variable_datatype y){
    variable_datatype tmp1 = x*x + y*y;
    variable_datatype tmp2 = cos(2 * PI * x) + cos(2 * PI * y);
    return 20. - 20. * std::exp(-0.2 * std::sqrt(0.5 * tmp1)) - std::exp(0.5 * tmp2);
}

double Griewank_f(variable_datatype x, variable_datatype y){
    variable_datatype tmp1 = x*x + y*y;
    variable_datatype tmp2 = cos(x / 1) * cos(y / std::sqrt(2));
    return tmp1 / 4000. - tmp2 + 1;
}
/*
///////////////////////////////
template <typename variable_datatype>
double VERY_Complicated(std::vector<variable_datatype> x){
    variable_datatype ret = 0;
    int x_size = x.size();
    for(int i = 0; i < x_size; ++i){
        ret += x[i]*x[i];
    }
    unsigned int microsecond = 1; // 1 sec = 1000000
    usleep(microsecond);//sleeps for 3 second
    return ret;
}


template <typename variable_datatype>
double Sphere_f(std::vector<variable_datatype> x){
    variable_datatype ret = 0;
    int x_size = x.size();
    for(int i = 0; i < x_size; ++i){
        ret += x[i]*x[i];
    }
    return ret;
}

template <typename variable_datatype>
double Schwefel_f(std::vector<variable_datatype> x){
    variable_datatype ret = 0;
    int x_size = x.size();
    for(int i = 0; i < x_size; ++i){
        variable_datatype temp = 0;
        for(int j = 0; j < i; ++j){
            temp += x[j];
        }
        ret += temp * temp;
    }
    return ret;
}

template <typename variable_datatype>
double Rosenbrock_f(std::vector<variable_datatype> x){
    variable_datatype ret = 0;
    int x_size = x.size();
    for(int i = 0; i < x_size-1; ++i){
        ret += 100 * (x[i+1] + x[i]*x[i]) + (x[i] - 1)*(x[i] - 1);
    }
    return ret;
}

template <typename variable_datatype>
double Schwefel2_f(std::vector<variable_datatype> x){
    variable_datatype ret = 0;
    int x_size = x.size();
    for(int i = 0; i < x_size; ++i){
        ret += x[i] * sin(std::sqrt(std::abs(x[i])));
    }
    return -ret;
}

template <typename variable_datatype>
double Rastrigin_f(std::vector<variable_datatype> x){
    variable_datatype ret = 0;
    int x_size = x.size();
    for(int i = 0; i < x_size; ++i){
        ret += x[i]*x[i] - 10 * cos(2 * PI * x[i]);
    }
    return 10 + ret;
}

template <typename variable_datatype>
double Ackley_f(std::vector<variable_datatype> x){
    variable_datatype tmp1 = 0;
    variable_datatype tmp2 = 0;
    int x_size = x.size();
    for(int i = 0; i < x_size; ++i){
        tmp1 += x[i]*x[i];
        tmp2 += cos(2 * PI * x[i]);
    }
    return 20. - 20. * std::exp(-0.2 * std::sqrt(1/x_size * tmp1)) - std::exp(1/x_size * tmp2);
}

template <typename variable_datatype>
double Griewank_f(std::vector<variable_datatype> x){
    variable_datatype tmp1 = 0;
    variable_datatype tmp2 = 1;
    int x_size = x.size();
    for(int i = 0; i < x_size; ++i){
        tmp1 += x[i]*x[i];
        tmp2 *= cos(x[i] / std::sqrt(i+1));
    }
    return tmp1 / 4000. - tmp2 + 1;
} */