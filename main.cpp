/*
* Discrete Fourier transform (C++)
* by Project Nayuki, 2017. Public domain.
* https://www.nayuki.io/page/how-to-implement-the-discrete-fourier-transform
*/

// Shared definitions
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
using std::size_t;
using std::vector;
using namespace std;
const int MAX_SIZE = 999999;


/*
* (Alternate implementation using only real numbers.)
* Computes the discrete Fourier transform (DFT) of the given complex vector.
* All the array arguments must have the same length.
*/
using std::cos;
using std::sin;
void computeDft(const vector<double> &inreal, const vector<double> &inimag,
	vector<double> &outreal, vector<double> &outimag) {

	size_t n = inreal.size();
	for (size_t k = 0; k < n; k++) {  // For each output element
		double sumreal = 0;
		double sumimag = 0;
		for (size_t t = 0; t < n; t++) {  // For each input element
			double angle = 2 * M_PI * t * k / n;
			sumreal += inreal[t] * cos(angle) + inimag[t] * sin(angle);
			sumimag += -inreal[t] * sin(angle) + inimag[t] * cos(angle);
		}
		outreal[k] = sumreal;
		outimag[k] = sumimag;
	}
}

void main() {
	double * time;
	double * ampl;
	time = new double[MAX_SIZE];
	ampl = new double[MAX_SIZE];

	int i = 0;
	const double PI2 = 6.2832;
	ifstream infile;
	infile.open("indata.txt");
	if (time == nullptr || ampl == nullptr)
		cout << "Error: memory could not be allocated";
	else
	{
		while (!infile.eof())
		{
			infile >> time[i] >> ampl[i];
			//ampls.push_back(temp);
			++i;
		}
		for (size_t i = 0; i < 100; i++)
		{
			cout << time[i] << ' ' << ampl[i] << "\n";
		}
	}

	const vector<double> ampls(ampl, ampl + i -1 );
	const vector<double> amplsimg(i-1,0);
	vector<double> outreal(i-1,0);
	vector<double> outimg(i-1,0);

	computeDft(ampls, amplsimg, outreal, outimg);


	ofstream outfile;
	outfile.open("outdata.txt");
	int p = 0;
	for (vector<double>::iterator it = outreal.begin(); it != outreal.end(); ++it)
	{
		outfile << ampl[p] << "\t\t" << *it << "\n";
		++p;
	}

}


