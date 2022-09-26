#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <vector>

#include "Function.h"


int main(int argc, char* argv[])
{
		
	// nome file output
	std::string outfilestr;
	std::cout << "Inserisci nome file di output " << std::endl;
	std::cin >> outfilestr;
	
	//apri filestream
	std::ofstream fout(outfilestr);
	
	// imposta la precisione dello stream
	fout.precision(15);
	
	
	std::vector<int> numb{10, 100, 1000, 1000000};
	
	
	// stampa valori
	for (auto i : numb)
	{
		fout << i << " " << fabs(Function::Gauss_integral(i,0)-0.5*sqrt(M_PI)) << " " 
			 << fabs(Function::Gauss_integral(i,1)-0.5*sqrt(M_PI)) << " " 
			 << fabs(Function::Gauss_integral(i,2)-0.5*sqrt(M_PI)) << std::endl; 
		
		
	}
	
	// chiudi filestream 
	fout.close();
	
	
	return 0;
	
}
