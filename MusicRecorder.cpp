#include "MusicRecorder.hpp"
#include "notes_en.hpp"
#include "notes_fr.hpp"


std::string freq2note(double freq, std::string language)
{
	std::map<int, std::map<std::string, double>> notes;
	if (language == "en")
		notes = notes_en;
	else if (language == "fr")
		notes = notes_fr;
		
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
	ACF = new double[size / 2];

	double AC0 = 0;
	for (size_t t = 0; t < size; ++t)
		AC0 += signal[t] * signal[t];
	ACF[0] = 1;

	AC0 = std::max(AC0, 0.50);  //PARAM 0.01
	for (size_t lag = 1; lag < (size / 2); ++lag)
	{
		ACF[lag] = 0;
		for (size_t t = lag; t < size; ++t)
			ACF[lag] += signal[t] * signal[t - lag];
		ACF[lag] /= AC0;
	}
	return ACF;
}

double* PSD(const double* signal, size_t size, CQT Q)
{
	//=============== FFT ===============//
	Complex * sig;
	sig = new Complex[size];
	//CArray FT; FT.resize(size);
	for (size_t t = 0; t < size; ++t)
	{
		//FT[t] = Complex(signal[t], 0);
		sig[t] = Complex(signal[t], 0.0);
		//std::cout << "  " << signal[t] << " <=> " << FT[t];
	}
	//fft(FT);
	Complex * QT = Q.transform(sig);

	//============== POWER ==============//
	double* PSD;

	//PSD = new double[size / 2];;
	//for (size_t freq_i = 0; freq_i < size / 2; ++freq_i)
	//	PSD[freq_i] = norm(FT[freq_i]) / (size);

	PSD = new double[168];;
	for (size_t freq_i = 0; freq_i < 168; ++freq_i)
	{
		PSD[freq_i] = norm(QT[freq_i]) / 168.0;
		//std::cout << freq_i << std::endl;
	}
	delete[] sig;

	return PSD;
}

MusicRecorder::MusicRecorder(double sample_rate, double sample_time, double lag) : sf::SoundRecorder()
{
	this->setProcessingInterval(sf::Time(sf::microseconds((sf::Int16)(sample_time * 1e6))));
	_delta_freq = 1.0 / sample_time;
	_buffer_i = 0;
	_buffer_size = (size_t)(sample_rate * sample_time);
	_buffer = new double[_buffer_size];
	double f_min = 440.0*(pow(2.0, (40-69.0) / 12.0));
	double f_max = 440.0*(pow(2.0, (123-69.0) / 12.0));
	CQT Q(_buffer_size/2, 168, f_min, f_max, sample_rate);
	QT = Q;
}

bool MusicRecorder::onProcessSamples(const sf::Int16* samples, std::size_t sampleCount)
{
	for (size_t sample_i = 0; sample_i < sampleCount; ++sample_i)
	{
		_buffer[_buffer_i] = (double)(samples[sample_i]) / pow(2.0, 16);
		_buffer_i = (_buffer_i + 1) % _buffer_size;
	}
	return true;
}

std::vector<double> MusicRecorder::getFreqs(size_t n)
{
	std::vector<double> freqs;
	// double* acf = ACF(_buffer, _buffer_size);
	// double* psd = PSD(acf, _buffer_size/2, QT);
	double* psd = PSD(_buffer, _buffer_size, QT);
	//std::cout << _buffer_size << std::endl;
	for (size_t i = 0; i < 168; ++i)
	{
		freqs.push_back(psd[i]);
	}
	/*
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
			freqs.push_back((best_lag + 2) * _delta_freq * 2);

			for (size_t t = std::max((size_t)0, best_lag - delta); t < std::min(_buffer_size / 2, best_lag + delta); ++t)
				psd[t] = 0;
		}
		else
			freqs.push_back(0);

		n -= 1;
	} while (n > 0);
	*/

	//delete[] acf;
	delete[] psd;
	return freqs;
}

double* MusicRecorder::getBuffer()
{
	double* buffer;
	buffer = new double[_buffer_size];

	for (size_t t = 0; t < _buffer_i; ++t)
		buffer[_buffer_size - _buffer_i + t] = _buffer[t];

	for (size_t t = _buffer_i; t < _buffer_size; ++t)
		buffer[t - _buffer_i] = _buffer[t];

	return buffer;
}

size_t MusicRecorder::getBufferSize()
{
	return _buffer_size;
}