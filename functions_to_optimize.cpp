
#include <iostream>
#include <vector>
#include <math.h>

#define PI 3.141592

/*
class Function{
    public:
        double virtual operator()(double x, double y) = 0;

        std::vector<double> get_range(){return _range}
        std::vector<double> get_optimum(){return _optimum}
        double get_value)(){return _value}
    protected:
        std::vector<double> _range; 
        std::vector<double> _optimum;
        double _value;
};
*/

double Sphere_f(double x, double y){
    double ret = x*x + y*y;
    return ret;
}

double Schwefel_f(double x, double y){
    double ret = 0;
    double temp = x+y;
    ret = x*x + temp*temp;
    return ret;
}

double Rosenbrock_f(double x, double y){
    double ret = 100. * (y + x*x) + (x - 1)*(x - 1);
    return ret;
}

double Schwefel2_f(double x, double y){
    double ret = x * sin(std::sqrt(std::abs(x)) +  y * sin(std::sqrt(std::abs(y))));
    return -ret;
}

double Rastrigin_f(double x, double y){
    double ret = x*x - 10 * cos(2 * PI * x) + y*y - 10 * cos(2 * PI * y);
    return 10 + ret;
}

double Ackley_f(double x, double y){
    double tmp1 = x*x + y*y;
    double tmp2 = cos(2 * PI * x) + cos(2 * PI * y);
    return 20. - 20. * std::exp(-0.2 * std::sqrt(0.5 * tmp1)) - std::exp(0.5 * tmp2);
}

double Griewank_f(double x, double y){
    double tmp1 = x*x + y*y;
    double tmp2 = cos(x / 1) * cos(y / std::sqrt(2));
    return tmp1 / 4000. - tmp2 + 1;
}

///////////////////////////////

double Sphere_f(std::vector<double> x){
    double ret = 0;
    int x_size = x.size();
    for(int i = 0; i < x_size; ++i){
        ret += x[i]*x[i];
    }
    return ret;
}

double Schwefel_f(std::vector<double> x){
    double ret = 0;
    int x_size = x.size();
    for(int i = 0; i < x_size; ++i){
        double temp = 0;
        for(int j = 0; j < i; ++j){
            temp += x[j];
        }
        ret += temp * temp;
    }
    return ret;
}

double Rosenbrock_f(std::vector<double> x){
    double ret = 0;
    int x_size = x.size();
    for(int i = 0; i < x_size-1; ++i){
        ret += 100 * (x[i+1] + x[i]*x[i]) + (x[i] - 1)*(x[i] - 1);
    }
    return ret;
}

double Schwefel2_f(std::vector<double> x){
    double ret = 0;
    int x_size = x.size();
    for(int i = 0; i < x_size; ++i){
        ret += x[i] * sin(std::sqrt(std::abs(x[i]));
    }
    return -ret;
}

double Rastrigin_f(std::vector<double> x){
    double ret = 0;
    int x_size = x.size();
    for(int i = 0; i < x_size; ++i){
        ret += x[i]*x[i] - 10 * cos(2 * PI * x[i]);
    }
    return 10 + ret;
}

double Ackley_f(std::vector<double> x){
    double tmp1 = 0;
    double tmp2 = 0;
    int x_size = x.size();
    for(int i = 0; i < x_size; ++i){
        tmp1 += x[i]*x[i];
        tmp2 += cos(2 * PI * x[i]);
    }
    return 20. - 20. * std::exp(-0.2 * std::sqrt(1/x_size * tmp1)) - std::exp(1/x_size * tmp2);
}

double Griewank_f(std::vector<double> x){
    double tmp1 = 0;
    double tmp2 = 1;
    int x_size = x.size();
    for(int i = 0; i < x_size; ++i){
        tmp1 += x[i]*x[i];
        tmp2 *= cos(x[i] / std::sqrt(i+1));
    }
    return tmp1 / 4000. - tmp2 + 1;
}
////////////////

int main(int argc, char *argv[]) {

  std::vector<double> x = [1, 4];

  return 0;
}
