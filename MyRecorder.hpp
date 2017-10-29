#ifndef MYRECORDER_H
#define MYRECORDER_H

#include <complex>
#include <valarray>

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

class MyRecorder : public sf::SoundRecorder
{
public:
	MyRecorder(sf::Time interval, size_t buffer_size, size_t PSD_size);

	/*virtual bool onStart() // optionnelle
	{
		return true;
	}
	*/
	
	bool onProcessSamples(const sf::Int16* samples, std::size_t sampleCount);
	float* getPSD();
	/*
	virtual void onStop() // optionnelle
	{

	}*/
private:
	float _sigma;
	std::size_t _PSD_size;
	float* _PSD;
	CArray _filter;
	CArray _buffer;
	std::size_t _buffer_index;
	std::size_t _buffer_size;
	
};

#endif
