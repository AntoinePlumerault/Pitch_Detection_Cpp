#include <SFML/Audio.hpp>
#include <complex>
#include <valarray>
#include <iostream>
#include "fft.hpp"
#include "NoteRecorder.hpp"
#include "notes.hpp"



double Plume_window(float alpha, size_t N, size_t n)
{
	return 0.5*(1.0 - cos((2.0 * PI * pow((double)n / (double)(N - 1), 1.0 / alpha))));
}

double Gaussian(double sigma, double x)
{
	return 1.0 / (sigma*sqrt(2.0 * PI)) * (exp(-1.0 / 2.0 * pow(x / sigma, 2.0)));
}


bool compare_second(const std::pair<size_t, float> p1, const std::pair<size_t, float> p2)
{
	return p1.second > p2.second;
}

bool compare_first(const std::pair<size_t, float> p1, const std::pair<size_t, float> p2)
{
	return p1.first > p2.first;
}


NoteRecorder::NoteRecorder(double sample_rate, double sample_time, double lag) : sf::SoundRecorder()
{
	this->setProcessingInterval(sf::Time(sf::microseconds((sf::Int16)(sample_time * 1e6))));
	_buffer_size = (size_t)(sample_rate * sample_time);
	_buffer.resize(_buffer_size); 
	
	_alpha = -log(lag/sample_time) / log(2.0);
	std::cout << _alpha << std::endl;
	_delta_freq = 1.0/sample_time;
	_buffer_i = 0;
}

bool NoteRecorder::onProcessSamples(const sf::Int16* samples, std::size_t sampleCount)
{
	for (size_t sample_i = 0; sample_i < sampleCount; ++sample_i)
	{
		_buffer[_buffer_i] = Complex((double)(samples[sample_i]) / pow(2.0, 16), 0.0);
		_buffer_i = (_buffer_i + 1) % _buffer_size;
	}
	return true;
}

size_t NoteRecorder::getPSDSize()
{
	return _buffer_size / 2;
}

CArray NoteRecorder::getFT()
{
	CArray FT; FT.resize(_buffer_size);

	for (size_t sample_i = 0; sample_i < _buffer_i; ++sample_i)
	{
		FT[sample_i] = _buffer[sample_i] * Plume_window(_alpha, _buffer_size, _buffer_i - sample_i);
	}
	for (size_t sample_i = _buffer_i; sample_i < _buffer_size; ++sample_i)
	{
		FT[sample_i] = _buffer[sample_i] * Plume_window(_alpha, _buffer_size, _buffer_i + _buffer_size - sample_i);
	}
	fft(FT);

	return FT;
}

double* NoteRecorder::getACF()
{
	double* signal; 
	signal = new double[_buffer_size];
	//================ OPTIONAL ===============//

	for (size_t sample_i = 0; sample_i < (_buffer_size - _buffer_i); ++sample_i)
	{
		signal[sample_i] = (_buffer[sample_i + _buffer_i]).real(); // The signal is real
	}
	for (size_t sample_i = (_buffer_size - _buffer_i); sample_i < _buffer_size; ++sample_i)
	{
		signal[sample_i] = (_buffer[sample_i - (_buffer_size - _buffer_i)]).real(); // The signal is real
	}

	//============== CALCUL DE L'ACF =============//
	double* ACF;
	ACF = new double[_buffer_size / 2];
	for (size_t lag = 0; lag < (_buffer_size / 2); ++lag)
	{
		ACF[lag] = 0;
		for (size_t sample_i = 0; sample_i < (_buffer_size - lag); ++sample_i)
		{
			ACF[lag] += signal[sample_i] * signal[sample_i + lag];
		}
	}

	delete[] signal;

	return ACF;
}

double* NoteRecorder::getACFPSD()
{
	CArray FT; FT.resize(_buffer_size / 2);
	double* ACF = getACF();
	//ACF = new double[_buffer_size / 2];
	for (size_t sample_i = 0; sample_i < (_buffer_size / 2); ++sample_i)
	{
		ACF[sample_i] = ACF[sample_i] /** Plume_window(1.0, _buffer_size / 4, sample_i)*/;
	}
	
	for (size_t sample_i = 0; sample_i < (_buffer_size / 2); ++sample_i)
	{
		FT[sample_i] = Complex(ACF[sample_i],0) /* Plume_window(2, _buffer_size/2, sample_i)*/;
	}
	
	fft(FT);
	
	double* PSD = new double[_buffer_size / 4];;

	double max_norm = 0;
	for (size_t signal_i = 0; signal_i < _buffer_size / 4; ++signal_i)
	{
		double norm_i = norm(FT[signal_i]);
		if (norm_i > max_norm) { max_norm = norm_i; }
		PSD[signal_i] = norm_i;
	}
	if (max_norm < 1) {
		max_norm = 10;
	}
	for (size_t signal_i = 0; signal_i < _buffer_size / 4; ++signal_i)
	{
		PSD[signal_i] /= max_norm;
	}
	
	delete[] ACF;
	return PSD;
}

std::string NoteRecorder::getBestNoteACF()
{
	double* PSD = getACFPSD();
	double freq = 0;
	double max_power = 1.0 / 5.0;
	for (size_t lag = 0; lag < _buffer_i / 4; ++lag)
	{
		if (max_power < PSD[lag])
		{
			max_power = PSD[lag];
			freq = (lag+2) * _delta_freq * 2;
		}
	}
	std::cout << freq << " aaa " << _delta_freq << std::endl;

	std::string note_name("NONE");
	double delta_freq = 10000;

	if (freq == 0) { return note_name; }

	for (size_t octave = 0; octave < 9; ++octave)
	{
		for (std::map<std::string, double>::iterator note = notes[octave].begin(); note != notes[octave].end(); ++note)
		{
			if (abs(note->second - freq) < delta_freq)
			{
				delta_freq = abs(note->second - freq);
				note_name = note->first;
			}
		}
	}
	return note_name;
}



double* NoteRecorder::getPSD()
{
	CArray FT = getFT();
	double* PSD = new double[_buffer_size / 2];;

	double max_norm = 0;
	for (size_t signal_i = 0; signal_i < _buffer_size / 2; ++signal_i)
	{
		double norm_i = norm(FT[signal_i]);
		if (norm_i > max_norm) { max_norm = norm_i; }
		PSD[signal_i] = norm_i;
	}
	/*
	for (size_t signal_i = 0; signal_i < _buffer_size / 2; ++signal_i)
	{
		PSD[signal_i] /= max_norm;
	}
	*/
	return PSD;
}

double* NoteRecorder::getCumulativePSD(size_t n) {
	double* CumulativePSD = getPSD();
	for (size_t signal_i = 0; signal_i < _buffer_size / 2; ++signal_i)
	{
		for (size_t k = 2; k <= n; ++k)
		{
			CumulativePSD[signal_i] *= 10 * CumulativePSD[std::min(_buffer_size / 2, k*signal_i)];
		}
	}
	return CumulativePSD;
}

std::vector<std::pair<size_t, double>> NoteRecorder::getPeaks(double threshold, double _freq_diff)
{
	double* PSD = getPSD();
	//double* PSD = getCumulativePSD(2);
	/*================ PEAK EXTRACTION ================*/
	std::vector<std::pair<size_t, double>> list_peaks;
	for (size_t freq_i = 1; freq_i < _buffer_size / 2 - 1; ++freq_i)
	{
		if ((PSD[freq_i] >= PSD[freq_i + 1] && PSD[freq_i] >= PSD[freq_i - 1]) && (PSD[freq_i] > threshold))
		{
			list_peaks.push_back(std::pair<size_t, float>(freq_i, PSD[freq_i]));
		}
	}
	/*================= PEAK SELECTION ================*/
	std::vector<std::pair<size_t, double>> list_best_peaks;
	for (size_t peak_i = 0; peak_i < list_peaks.size(); ++peak_i)
	{
		bool keep_this_peak = true;
		for (size_t best_peak_i = 0; best_peak_i < list_best_peaks.size(); ++best_peak_i)
		{
			if (abs((int)(list_best_peaks[best_peak_i].first - list_peaks[peak_i].first)) < pow(2.0, _octave_min + 1) * _delta_freq)
			{
				if (list_best_peaks[best_peak_i].second <= list_peaks[peak_i].second)
					list_best_peaks.erase(list_best_peaks.begin() + best_peak_i);
				else
					keep_this_peak = false;
				break;
			}
		}
		if (keep_this_peak)
			list_best_peaks.push_back(list_peaks[peak_i]);
	}

	delete[] PSD;

	return list_best_peaks;
}

double* NoteRecorder::getDeviationDistribution(size_t lag_max) {
	size_t sig = 10;

	std::vector<std::pair<size_t, double>> list_peaks = getPeaks(0.05, 80);
	for (size_t i = 0; i < 1; ++i) list_peaks.push_back(std::pair<size_t, double>(0, 0)); // Addition of a peak at the zero frequency
	sort(list_peaks.begin(), list_peaks.end(), compare_first);

	double* deviation_distribution = new double[getPSDSize()];
	for (size_t i = 0; i < getPSDSize(); ++i) { deviation_distribution[i] = 0; }

	for (size_t peak_i = 0; peak_i < list_peaks.size(); ++peak_i)
	{
		for (size_t peak_j = std::min(peak_i + 1, list_peaks.size()); peak_j < std::min(peak_i + lag_max, list_peaks.size()); ++peak_j)
		{
			size_t deviation = (list_peaks[peak_i].first - list_peaks[peak_j].first);
			//================ SMOOTHING ================//
			for (int delta = -2 * (int)sig; delta < 2 * (int)sig; ++delta)
			{
				if (deviation + delta > 0 && deviation + delta < getPSDSize())
					deviation_distribution[deviation + delta] += 2 * sig * Gaussian(sig, delta);
			}
			//===========================================//
		}
	}
	return deviation_distribution;
}

double NoteRecorder::getBestFreq()
{
	double* deviation_distribution = getDeviationDistribution(3);

	std::vector<std::pair<size_t, double>> list_peaks;
	for (size_t freq_i = 1; freq_i < _buffer_size / 2 - 1; ++freq_i)
	{
		if ((deviation_distribution[freq_i] >= deviation_distribution[freq_i + 1] && deviation_distribution[freq_i] >= deviation_distribution[freq_i - 1]) && (deviation_distribution[freq_i] > 0))
		{
			list_peaks.push_back(std::pair<size_t, float>(freq_i, deviation_distribution[freq_i]));
		}
	}
	delete[] deviation_distribution;
	sort(list_peaks.begin(), list_peaks.end(), compare_second);
	if (list_peaks.size() > 0)
		return list_peaks[0].first * _delta_freq;
	else
		return 0;
}
/*
std::string NoteRecorder::getBestNote()
{
	double freq = getBestFreq();
	std::string note_name("NONE");

	if (freq == 0) { return note_name; }

	//std::map<std::string, double>::iterator old_note = notes[0].begin();
	std::string old_note_name = notes[0].begin()->first;
	double old_note_freq = notes[0].begin()->second;
	for (size_t octave = 0; octave < 9; ++octave)
	{
		for (std::map<std::string, double>::iterator note = notes[octave].begin(); note != notes[octave].end(); ++note)
		{
			if (freq <= note->second && freq > old_note_freq)
			{
				std::cout << old_note_freq << "<" << freq << "<" << note->second << "        \r";
				if (note->second - freq < freq - old_note_freq)
					note_name = note->first;
				else
					note_name = old_note_name;
				break;
			}
			old_note_name = note->first;
			old_note_freq = note->second;
		}
		if (note_name != std::string("NONE")) { break; }
	}
	return note_name;
}*/

std::string NoteRecorder::getBestNote()
{
	double freq = getBestFreq();

	std::string note_name("NONE");
	double delta_freq = 10000;

	if (freq == 0) { return note_name; }

	for (size_t octave = 0; octave < 9; ++octave)
	{
		for (std::map<std::string, double>::iterator note = notes[octave].begin(); note != notes[octave].end(); ++note)
		{
			if (abs(note->second - freq) < delta_freq)
			{
				delta_freq = abs(note->second - freq);
				note_name = note->first;
			}
		}
	}
	return note_name;
}