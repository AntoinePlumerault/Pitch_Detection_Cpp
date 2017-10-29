#ifndef NOTESRECORDER_HPP
#define NOTESRECORDER_HPP

#include <complex>
#include <valarray>

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

class NotesRecorder : public sf::SoundRecorder
{
public:
	NotesRecorder(double sample_rate, double sample_time, double lag);
	bool onProcessSamples(const sf::Int16* samples, std::size_t sampleCount);
	size_t getPSDSize();
	double * getACF();
	double * getACFPSD();
	std::vector<double> getFreqs(size_t n);
	double* getBufferTest() { return _buffer_test; }
	

private:
	double* _buffer_test;
	double _delta_freq;
	double* _buffer;
	size_t _buffer_i;
	size_t _buffer_size;

};

double* ACF(const double* signal, size_t size);
double* PSD(const double* signal, size_t size);
std::string freq2note(double freq);
#endif