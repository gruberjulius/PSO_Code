#include<iostream>


int main()
{	
	double a = 0.0;
	unsigned int loopsize = 1700
	for(int i = 0; i < loopsize; i++){ //Griewank 0
		std::cout << "0.0" << std::endl; 
	}
	for(int i = 0; i < loopsize; i++){ //Ackley 0
		std::cout << "-2.71828.0" << std::endl; 
	}
	for(int i = 0; i < loopsize; i++){ //Holdertable 0
		std::cout << "0.0" << std::endl; 
	}
	for(int i = 0; i < loopsize; i++){
		std::cout << "0.0" << std::endl; 
	}
	for(int i = 0; i < loopsize; i++){
		std::cout << "-100.0" << std::endl; 
	}
	for(int i = 0; i < loopsize; i++){
		std::cout << "21.0" << std::endl; 
	}
	return 0;
}

