#include "simInfo.h"
#include <iostream>
#include <math.h>

using namespace std;


extern double hx;
extern double hy;


// constructor
simInfo::simInfo() {
}


// destructor
simInfo::~simInfo() {
}


// ask info about simulation via iostream
void simInfo::askInfo() {
	
	double dInfo;
	int iInfo;

	// iteration parameters
	cout << "dt (default 1.0) ";
	cin >> dInfo;
	dt_ = (dInfo ? dInfo : 1.0 );
	
	cout << "nTimeSteps (default 500) ";
	cin >> iInfo;
	nTimeSteps_ = (iInfo ? iInfo : 500 );
	
	cout << "freqPrint (default 5) ";
	cin >> iInfo;
	freqPrint_ = (iInfo ? iInfo : 5 );
	

	// simulation parameters
	cout << "qx in [ " << -M_PI/hx << " , " << M_PI/hx << "] ";
	cin >> qx_;
	
	cout << "qy in [ " << -M_PI/hy << " , " << M_PI/hy << "] ";
	cin >> qy_;
	
	cout << "V0 (soglia ~ " << 0.5/(hx*hy) << ") ";
	cin >> V0_;


  return;

}



