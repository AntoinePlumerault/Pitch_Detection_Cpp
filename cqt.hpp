#ifndef CQT_H
#define CQT_H

#include <Eigen/SparseCore>
#include "fft.hpp"

class Cqt
{
public:
	Cqt(float min_freq, float max_freq, size_t n_bins, float fs);
	void run(CArray& x);

private:
	CArray cq_hanning_window(int length);

	size_t _n_bins;
	float _min_freq;
	float _max_freq;
	float _fs;
	float _Q;
	int _fftlen;
	std::vector<CArray> _ker;
};

double* ACF(const double* signal, size_t size);
double* PSD(const double* signal, size_t size);
std::string freq2note(double freq, std::string language);

#endif