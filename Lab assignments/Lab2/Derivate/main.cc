#include <iostream>
#include <string>

#include "Function.h"


// funzione per stampare su file il confronto tra f e g generiche
void writeFile(std::string & outfilestr, Function & f, Function & g);


int main(int argc, char* argv[])
{
	
	double x_0 = 0;
	int n = 100;
	double h = 0.1;
	
	
	// crea l'oggetto e inizializza la griglia di valori
	Function sinus(x_0, n, h, true);
	
	// calcola la funzione su tutta la griglia
	sinus.Sin();
	
	
	// calcola separatamente i 4 metodi delle derivate
	sinus.Derivative(0);	// O(h) destra
	sinus.Derivative(1);	// O(h) sinistra
	sinus.Derivative(2);	// O(h^2)
	sinus.Derivative(3);	// O(h^4)
	
	
	// funzione coseno per confronto
	Function cosinus(x_0, n, h, false);
	cosinus.Cos();
	
	
	
	// nome del file di output da tastiera
	std::string file;
	std::cout << "Inserisci nome file di output" << std::endl;
	std::cin >> file;
	
	
	// scrivi sul file
	writeFile(file, sinus, cosinus);
	
		
	return 0;
	
}
