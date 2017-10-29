#ifndef NOTERECORDER_HPP
#define NOTERECORDER_HPP

#include <complex>
#include <valarray>



typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

class NoteRecorder : public sf::SoundRecorder
{
public:
	NoteRecorder(double sample_rate, double sample_time, double lag);
	bool onProcessSamples(const sf::Int16* samples, std::size_t sampleCount);
	size_t getPSDSize();
	CArray getFT();
	double * getACF();
	double * getACFPSD();
	double* getPSD();
	std::string NoteRecorder::getBestNoteACF();
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
