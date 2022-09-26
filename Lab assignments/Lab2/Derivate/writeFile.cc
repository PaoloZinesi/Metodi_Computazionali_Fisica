#include <iostream>
#include <string>
#include <fstream>

#include "Function.h"

// funzione per stampare su file il confronto tra f' e g
void writeFile(std::string & outfilestr, Function & f, Function & g)
{
	//apri filestream
	std::ofstream fout(outfilestr);
	
	// imposta la precisione dello stream
	fout.precision(10);
	
	
	// stampa g - f'
	for(int k=0; k<g.getN_(); k++)
	{
		fout << f.getX(k) << " " << (g.getF(k) - f.getF1(0, k)) << " " << (g.getF(k) - f.getF1(1, k)) << " " 
			<< (g.getF(k) - f.getF1(2, k)) << " " << (g.getF(k) - f.getF1(3, k)) << std::endl;
	}
	
	
	// chiudi filestream 
	fout.close();
	
}