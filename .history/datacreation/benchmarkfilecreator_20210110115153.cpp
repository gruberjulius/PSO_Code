#include<iostream>


int main()
{	
	double a = 0.0;
	for(int i = 0; i < 1600; i++){ //Griewank 0
		std::cout << "0.0" << std::endl; 
	}
	for(int i = 0; i < 1600; i++){ //Ackley 0
		std::cout << "-2.71828.0" << std::endl; 
	}
	for(int i = 0; i < 1600; i++){ //Holdertable 0
		std::cout << "0.0" << std::endl; 
	}
	for(int i = 0; i < 1600; i++){
		std::cout << "0.0" << std::endl; 
	}
	for(int i = 0; i < 1600; i++){
		std::cout << "-100.0" << std::endl; 
	}
	for(int i = 0; i < 1600; i++){
		std::cout << "21.0" << std::endl; 
	}
	return 0;
}

//Rastrigin 0