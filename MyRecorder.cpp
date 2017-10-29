#include <SFML/Audio.hpp>
#include <complex>
#include <valarray>
#include "fft.hpp"
#include "MyRecorder.hpp"



typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

double Hann_window(size_t n, size_t i) 
{
	return 0.5*(1 - cos((2 * i) / (n - 1)));
}
/*
double Gaussian(double sigma, double x) 
{
	return 1.0 / (sigma*sqrt(2.0 * PI)) * (exp(-1.0 / 2.0 * pow(x / sigma, 2.0)));
}
*/

MyRecorder::MyRecorder(sf::Time interval, size_t buffer_size, size_t PSD_size) : sf::SoundRecorder() {
	this->setProcessingInterval(interval);
	
	_buffer.resize(buffer_size);
	_buffer_size = buffer_size;

	_PSD = new float[PSD_size];
	_PSD_size = PSD_size;
	
	_sigma = 3.0/6000.0;
	_buffer_index = 0;

	_filter.resize(2*PSD_size);
	for (size_t filter_i = 0; filter_i < PSD_size; ++filter_i)
	{
		Complex complex(Gaussian(_sigma, filter_i), 0.0);
		_filter[filter_i] = complex;
	}
	for (size_t filter_i = PSD_size; filter_i < 2*PSD_size; ++filter_i)
	{
		Complex complex(Gaussian(_sigma, filter_i-2*PSD_size), 0.0);
		_filter[filter_i] = complex;
	}
	ifft(_filter);
}

bool MyRecorder::onProcessSamples(const sf::Int16* samples, std::size_t sampleCount)
{
	// if (sampleCount == 0) { return true; }
	size_t n_iter = _buffer_size - _buffer_index;
	
	for (size_t sample_i = 0; sample_i < std::min(n_iter,sampleCount) ; ++sample_i) 
	{
		Complex complex((double)(samples[sample_i]) / pow(2.,16), 0.0);
		_buffer[_buffer_index] = complex;
		_buffer_index += 1;
	}

	if (_buffer_index == _buffer_size)
	{
		_buffer_index = 0;
		
		for (size_t sample_i = n_iter; sample_i < sampleCount; ++sample_i)
		{
			Complex complex((double)(samples[sample_i]) / pow(2., 16), 0.0);
			_buffer[_buffer_index] = complex;
			_buffer_index += 1;
		}
	}
	return true;
}

float* MyRecorder::getPSD() 
{
	CArray FT; FT.resize(2*_PSD_size);
	
	for (size_t sample_i = 0; sample_i < _buffer_index; ++sample_i)
	{
		FT[sample_i] = _buffer[sample_i] * Hann_window(_buffer_size, sample_i + (_buffer_size - _buffer_index));
	}
	Complex zero(0.0, 0.0);
	for (size_t sample_i = _buffer_index; sample_i < 2 * _PSD_size - (_buffer_size - _buffer_index); ++sample_i)
	{
		FT[sample_i] = zero;
	}
	for (size_t sample_i = 2 * _PSD_size - (_buffer_size - _buffer_index); sample_i < 2 * _PSD_size; ++sample_i)
	{
		size_t buffer_i = sample_i - (2 * _PSD_size - (_buffer_size - _buffer_index)) + _buffer_index;
		FT[sample_i] = _buffer[buffer_i] * Hann_window(_buffer_size, buffer_i - _buffer_index);
	}
	
	for (size_t sample_i = 0; sample_i < 2 * _PSD_size; ++sample_i)
	{
		FT[sample_i] *= _filter[sample_i];
	}
	
	fft(FT);
	//FT = _filter;
	for (unsigned int signal_i = 0; signal_i < _PSD_size; ++signal_i)
	{
		_PSD[signal_i] = (float)norm(FT[signal_i]);// *signal_i;
	}
	return _PSD;
}