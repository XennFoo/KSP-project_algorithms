
#include "pch.h"
#include <iostream>
#include <random>
#include <fstream>
#define N 20
#define k 20



int main()
{
	std::ofstream myfile;
	myfile.open("testfile.txt");
	std::random_device rd; 
	std::mt19937 eng(rd());
	std::uniform_int_distribution<> distr(1, 262110);
	myfile << "n  " << N << std::endl;
	myfile << "k  " << k << std::endl;
	for (int n = 0; n < N; ++n) {
		myfile << "p "<< distr(eng) << " " << distr(eng) << std::endl;
	}

}
