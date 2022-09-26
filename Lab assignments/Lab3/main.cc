// output: grafico dell'errore del periodo in relazione all'angolo

#include <iostream>
#include <math.h>
#include <fstream>

struct config
{
	double x;
	double v;
	
};

using namespace std;

constexpr double g = 9.80652;			// accelerazione gravitazionale (m*s^(-2))
constexpr double l = 1.0;				// lunghezza del pendolo (m)


// f(x1,x2) = (y1,y2)
config f (config alpha) {
	
	config beta;
	
	beta.x = alpha.v;
	beta.v = (-g) * sin(alpha.x);
	
	return beta;
}


// theta_0 e ∆t in input, periodo del pendolo in output
double T_reale (double th_iniz, double dt)
{
	// inizializzazione
	config theta{th_iniz, 0.0};
	config theta_new; 		// tenere traccia dei valori precedenti
	
	config Y1, Y2, Y3, Y4;
	double T;
	int num_it = 0;
	
	
	
	do
	{
		Y1 = theta;
		
		Y2.x = theta.x + (f(Y1).x)*(dt/2.0);
		Y2.v = theta.v + (f(Y1).v)*(dt/2.0);
		
		Y3.x = theta.x + (f(Y2).x)*(dt/2.0);
		Y3.v = theta.v + (f(Y2).v)*(dt/2.0);
		
		Y4.x = theta.x + (f(Y3).x)*dt;
		Y4.v = theta.v + (f(Y3).v)*dt;
		
		
		theta_new.x = theta.x + ((f(Y1).x) + 2*(f(Y2).x) + 2*(f(Y3).x) + (f(Y4).x))*(dt/6.0);
		theta_new.v = theta.v + ((f(Y1).v) + 2*(f(Y2).v) + 2*(f(Y3).v) + (f(Y4).v))*(dt/6.0);
		
		
		if (theta.v > 0 && theta_new.v < 0)
		{
			
			T = (theta.v * dt) / (theta.v - theta_new.v) + num_it * dt;
			return T;
			
		}
		
		
		num_it++;
		theta = theta_new;
		
		
	}
	while (num_it < lround(8*M_PI*sqrt(l/g)/dt));
	
	// se va troppo fuori ritorna un tempo negativo (non fisico)
	return -1;
	
}



int main()
{
	
	double dt;							// intervallo temporale
	double dtheta;						// intervallo angolare
	
	
	
	cout << "Inserire valore di ∆t (s)" << endl;
	cin >> dt;
	cout << "Inserire valore di ∆ø (rad)" << endl;
	cin >> dtheta;
	
	
	
	int N = lround(M_PI/(2*dtheta));	// N di sampling dell'angolo

	
	// output dati
	ofstream output("out.txt");
	output.precision(12);


	for(int i=1; i<N; i++)
	{
		output << i*dt << " " << ((T_reale(i*dtheta, dt)/(2*M_PI*sqrt(l/g)))-1)*100.0 << endl;
		
	}
	
	output.close();
	
	
}
