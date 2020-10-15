#include <iostream>
#include <string>

#include "Function.h"

void writeFile(std::string & outfilestr, Function & f);


int main(int argc, char* argv[])
{
	
	double x_0;
	int n;
	double h;
	
	// inserimento valori
	std::cout << "Inserisci il valore del primo punto della griglia (x_0) " << std::endl;
	std::cin >> x_0;
	std::cout << "Inserisci la distanza di sampling (h) " << std::endl;
	std::cin >> h;
	std::cout << "Inserisci il numero di punti (n) " << std::endl;
	std::cin >> n;
	
	
	// crea l'oggetto e inizializza automaticamente tutta la griglia di valori
	Function sinus(x_0, n, h);
	
	
	// calcola la funzione su tutta la griglia
	sinus.Sin();
	
	
	// nome del file di input da tastiera
	std::string file;
	std::cout << "Inserisci nome file di output" << std::endl;
	std::cin >> file;
	
	// scrivi sul file tutta la funzione
	writeFile(file, sinus);
	
		
	return 0;
	
}