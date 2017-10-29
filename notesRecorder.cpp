#include <SFML/Audio.hpp>
#include <complex>
#include <valarray>
#include <iostream>
#include <algorithm>
#include "fft.hpp"
#include "NotesRecorder.hpp"
//#include "notes_en.hpp"
#include "notes_fr.hpp"

std::string freq2note(double freq)
{
	std::string note_name("    ");
	double freq_error = INFINITY;

	if (freq == 0) { return note_name; }

	for (size_t octave = 0; octave < 9; ++octave)
	{
		for (std::map<std::string, double>::iterator note = notes[octave].begin(); note != notes[octave].end(); ++note)
		{
			if (abs(note->second - freq) < freq_error)
			{
				freq_error = abs(note->second - freq);
				note_name = note->first;
			}
		}
	}
	return note_name;
}

double* ACF(const double* signal, size_t size)
{
	double* ACF;
	ACF = new double[size/2];

	double AC0 = 0;
	for (size_t t = 0; t < size; ++t)
		AC0 += signal[t] * signal[t];
	ACF[0] = 1;
	
	AC0 = std::max(AC0, 0.010);  //PARAM
	for (size_t lag = 1; lag < (size / 2); ++lag)
	{
		ACF[lag] = 0;
		for (size_t t = lag; t < size; ++t)
			ACF[lag] += signal[t] * signal[t - lag];
		ACF[lag] /= AC0;
	}
	return ACF;
}

double* PSD(const double* signal, size_t size)
{
	//=============== FFT ===============//
	CArray FT; FT.resize(size);
	for (size_t t = 0; t < size; ++t)
		FT[t] = Complex(signal[t], 0);
	fft(FT);

	//============== POWER ==============//
	double* PSD = new double[size / 2];;
	for (size_t freq_i = 0; freq_i < size / 2; ++freq_i)
		PSD[freq_i] = norm(FT[freq_i])/(size);

	return PSD;
}

NotesRecorder::NotesRecorder(double sample_rate, double sample_time, double lag) : sf::SoundRecorder()
{
	this->setProcessingInterval(sf::Time(sf::microseconds((sf::Int16)(sample_time * 1e6))));
	_delta_freq = 1.0 / sample_time;
	_buffer_i = 0;
	_buffer_size = (size_t)(sample_rate * sample_time);
	//_buffer.resize(_buffer_size);
	_buffer = new double[_buffer_size];
	_buffer_test = new double[_buffer_size/2];
}

bool NotesRecorder::onProcessSamples(const sf::Int16* samples, std::size_t sampleCount)
{
	for (size_t sample_i = 0; sample_i < sampleCount; ++sample_i)
	{
		_buffer[_buffer_i] = (double)(samples[sample_i]) / pow(2.0, 16);
		_buffer_i = (_buffer_i + 1) % _buffer_size;
	}
	return true;
}

size_t NotesRecorder::getPSDSize()
{
	return _buffer_size / 2;
}

double* NotesRecorder::getACF()
{
	double* signal;
	signal = new double[_buffer_size];
	//================ OPTIONAL ? ===============//
	for (size_t sample_i = 0; sample_i < (_buffer_size - _buffer_i); ++sample_i)
		signal[sample_i] = _buffer[sample_i + _buffer_i];

	for (size_t sample_i = (_buffer_size - _buffer_i); sample_i < _buffer_size; ++sample_i)
		signal[sample_i] = _buffer[sample_i - (_buffer_size - _buffer_i)];

	double* acf = ACF(signal, _buffer_size);
	delete[] signal;
	return ACF(_buffer, _buffer_size);
	//return acf;
}

double* NotesRecorder::getACFPSD()
{
	double* ACF = getACF();
	double* ACFPSD = PSD(ACF, _buffer_size/2);
	delete[] ACF;
	return ACFPSD;
}

std::vector<double> NotesRecorder::getFreqs(size_t n)
{
	std::vector<double> freqs;
	double* acf = ACF(_buffer, _buffer_size);
	double* psd = PSD(acf, _buffer_size / 2);
	size_t delta = 3;
	do
	{	
		size_t best_lag = 0;
		for (size_t lag = 0; lag < _buffer_size / 4; ++lag)
		{
			if ((psd[best_lag] < psd[lag]) && (psd[lag] >= 0.5))
				best_lag = lag;
		}
		
		if (best_lag != 0)
		{
			freqs.push_back((best_lag + 2) * _delta_freq);

			for (size_t t = std::max((size_t)0, best_lag - delta); t < std::min(_buffer_size / 2, best_lag + delta); ++t)
					psd[t] = 0;
		}
		else
			freqs.push_back(0);

		n -= 1;
	} while (n > 0);

	delete[] acf;
	delete[] psd;
	return freqs;
}
/*
std::vector<double> NotesRecorder::getFreqs()
{
	std::vector<double> freqs;
	bool stop = false;
	int n = 0;
	double* acf = getACF();
	double* temp = new double[_buffer_size/2];
	
	do
	{	
		n += 1;
		double* psd = PSD(acf, _buffer_size/2);
		size_t best_lag = 0;
		for (size_t lag = 0; lag < _buffer_size/4; ++lag)
		{
			if ((psd[best_lag] < psd[lag]) && (psd[lag] >= 0.5))
				best_lag = lag;
		}	
		delete[] psd;
		if (best_lag != 0)
		{
			freqs.push_back((best_lag + 2) * _delta_freq);
	
			for (size_t t = 0; t < _buffer_size/2; ++t)
			{
				if (t > _buffer_size/2 - best_lag*2)
				{
					temp[t] = -acf[t - best_lag*2] + 1.0*acf[t];
				}
				else
				{			
					temp[t] = acf[t] - 1.0*acf[t + best_lag*2];
				}
			}
			
			for (size_t t = 0; t < _buffer_size / 2; ++t)
			{
				acf[t] = temp[t];
				_buffer_test[t] = acf[t];
			}
		}	
		else
		{
			freqs.push_back(0);
			stop = true;
		}
	} while ( n < 2);
	//delete[] temp_buffer;
	delete[] temp;
	delete[] acf;
	return freqs;
}
*/