#ifndef FFT_H
#define FFT_H

#include <complex>
#include <valarray>



const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

void fft(CArray& x);
void ifft(CArray& x);

#endif