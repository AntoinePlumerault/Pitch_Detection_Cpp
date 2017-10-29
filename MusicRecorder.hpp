#ifndef MUSICRECORDER_HPP
#define MUSICRECORDER_HPP

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

class MusicRecorder : public sf::SoundRecorder
{
public:
	MusicRecorder(double sample_rate, double sample_time, double lag);
	bool onProcessSamples(const sf::Int16* samples, size_t sampleCount);	
	std::vector<double> getFreqs(size_t n);
	double* getBuffer();
	size_t getBufferSize();

private:
	double* _buffer;
	double _delta_freq;
	size_t _buffer_i;
	size_t _buffer_size;

};

double* ACF(const double* signal, size_t size);
double* PSD(const double* signal, size_t size);
std::string freq2note(double freq, std::string language);

#endif