#ifndef SPECTRALRECORDER_H
#define SPECTRALRECORDER_H

#include <complex>
#include <valarray>



typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

class SpectralRecorder : public sf::SoundRecorder
{
public:
	SpectralRecorder(size_t octave_min, size_t octave_max, double alpha);
	bool onProcessSamples(const sf::Int16* samples, std::size_t sampleCount);
	size_t getPSDSize();
	CArray getFT();
	double* getPSD();
	double * getCumulativePSD(size_t n);
	std::vector<std::pair<size_t, double>> getPeaks(double threshold, double delta_freq);
	double* getDeviationDistribution(size_t lag_max);
	double getBestFreq();
	std::string getBestNote();

private:
	double _alpha;
	double _delta_freq;
	CArray _buffer;
	size_t _octave_min;
	size_t _octave_max;
	size_t _buffer_i;
	size_t _buffer_size;

};

#endif