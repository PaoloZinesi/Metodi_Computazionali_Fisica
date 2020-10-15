#include <iostream>
#include <string>
#include <fstream>

#include "Function.h"

// funzione per aprire un file e stamparci le informazioni contenute nella funzione f
void writeFile(std::string & outfilestr, Function & f)
{
	//apri filestream
	std::ofstream fout(outfilestr);
	
	// imposta la precisione dello stream
	fout.precision(10);
	
	// stampa valori
	for (int i=0; i<f.getN_(); i++)
	{
		fout << f.getX(i) << " " << f.getF(i) << std::endl;
	}
	
	// chiudi filestream 
	fout.close();
	
}