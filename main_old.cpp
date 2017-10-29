#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <iostream>
#include <thread>
#include <chrono>
#include "plot.hpp"
//#include "MyRecorder.hpp"
//#include "SpectralRecorder.hpp"
//#include "NoteRecorder.hpp"
//#include "notesRecorder.hpp"
#include "MusicRecorder.hpp"
#include "fft.hpp"



int main()
{
	if (!sf::SoundBufferRecorder::isAvailable())
	{
		printf("error");
	}

	unsigned int height = 800;
	unsigned int length = 800;
	unsigned int margin = 20;

	sf::RenderWindow window(sf::VideoMode(height, length), "FFT");
	sf::Font font;
	font.loadFromFile("MonospaceTypewriter.ttf");

	double sample_rate = 8000.0;											
	double sample_time = 1.0 / 3.0;
	NotesRecorder recorder(sample_rate, sample_time, 1.0/10.0);
	size_t buffer_size = recorder.getPSDSize();
	std::cout << buffer_size << std::endl;
	recorder.start((size_t)sample_rate);

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}
		
		//===================== DISPLAY =====================//
		window.clear();
		/*
		std::cout << "\r" << recorder.getBestFreq() << "      ";		
		
		sf::Text text(recorder.getBestNote(), font);
		text.setCharacterSize(300);
		text.setStyle(sf::Text::Bold);
		window.draw(text);
		*/
		/*	
		double* ACF = recorder.getACF();
		std::vector<double> X; X.resize(buffer_size / 2.);
		std::vector<double> Y; Y.resize(buffer_size / 2.);
		for (unsigned int signal_i = 0; signal_i < buffer_size / 2; ++signal_i)
		{
			Y[signal_i] = ACF[signal_i] + 1;
			//std::cout << PSD[signal_i] << std::endl;
			X[signal_i] = (double)signal_i;
		}
		*/
		/*
		double* ACFPSD = recorder.getACFPSD();
		std::vector<double> X; X.resize(buffer_size / 2);
		std::vector<double> Y; Y.resize(buffer_size / 2);
		for (unsigned int signal_i = 0; signal_i < buffer_size / 2; ++signal_i)
		{
			Y[signal_i] = ACFPSD[signal_i];
			//std::cout << PSD[signal_i] << std::endl;
			X[signal_i] = (double)signal_i;
		}
		*/
		std::vector<double> freqs = recorder.getFreqs(2);
		for(double freq : freqs)
			std::cout << freq << " - ";
		std::cout << std::endl;
		
		/*
		double* ACF =  recorder.getACF();
		std::vector<double> X; X.resize(buffer_size);
		std::vector<double> Y; Y.resize(buffer_size);
		for (unsigned int signal_i = 0; signal_i < buffer_size; ++signal_i)
		{
		Y[signal_i] = ACF[signal_i];
		X[signal_i] = (double)signal_i;
		}
		delete[] ACF;
		*/

		sf::Text note1(freq2note(recorder.getFreqs(2)[0]) + freq2note(recorder.getFreqs(2)[1]), font);
		note1.setCharacterSize(150);
		note1.setStyle(sf::Text::Bold);
		window.draw(note1);
		
		
		double* test = recorder.getBufferTest();
		std::vector<double> X; X.resize(buffer_size);
		std::vector<double> Y; Y.resize(buffer_size);
		for (unsigned int signal_i = 0; signal_i < buffer_size; ++signal_i)
		{
			Y[signal_i] = test[signal_i];
			//std::cout << PSD[signal_i] << std::endl;
			X[signal_i] = (double)signal_i;
		}
		
		/*
		sf::Text text(recorder.getBestNoteACF(), font);
		text.setCharacterSize(300);
		text.setStyle(sf::Text::Bold);
		window.draw(text);
		*/

		/*
		double* PSD = recorder.getPSD(); // PSD : Power Spectral Density
		std::vector<double> X; X.resize(buffer_size / 2.);
		std::vector<double> Y; Y.resize(buffer_size / 2.);
		for (unsigned int signal_i = 0; signal_i < buffer_size / 2; ++signal_i)
		{
			Y[signal_i] = (10 + PSD[signal_i]);
			//std::cout << PSD[signal_i] << std::endl;
			X[signal_i] = (double)signal_i;
		}
		*/
		
		/*
		double* PSD = recorder.getCumulativePSD(2); // PSD : Power Spectral Density
		std::vector<double> X; X.resize(buffer_size / 2.);
		std::vector<double> Y; Y.resize(buffer_size / 2.);
		for (unsigned int signal_i = 0; signal_i < buffer_size / 2; ++signal_i)
		{
			Y[signal_i] = (10 + PSD[signal_i]);
			X[signal_i] = (double)signal_i;
		}
		*/

		/*
		double* DEV = recorder.getDeviationDistribution(2);
		std::vector<double> X; X.resize(buffer_size/2.);
		std::vector<double> Y; Y.resize(buffer_size/2.);
		for (unsigned int signal_i = 0; signal_i < buffer_size/2; ++signal_i) 
		{
			Y[signal_i] = (10 + DEV[signal_i]);
			X[signal_i] = (double)signal_i;
		}
		*/

		/*
		std::vector<std::pair<size_t, double>> list_best_peaks = recorder.getPeaks(0.05, 20);
		
		for (std::pair<size_t, double> peak : list_best_peaks)
		{
			sf::CircleShape circle;
			circle.setFillColor(sf::Color::Red);
			circle.setOrigin(4, 4);
			circle.setRadius(4);
			circle.setPosition((double)(peak).first * 760. / (buffer_size/2) + 20, 400);
			window.draw(circle);
			
			//sf::Text text(std::string(peak.first), font);
			//text.setCharacterSize(10);
			//text.setStyle(sf::Text::Bold);
			//text.setPosition((double)(peak).first * 760. / (buffer_size / 2) + 20, 400);
			//window.draw(text);
			
		}
		*/
		
		std::vector<sf::Vertex> list_points = plot(Y, X, "title", 4000, 2);

		for (unsigned int line_i = 1; line_i < list_points.size(); ++line_i)
		{
			sf::Vertex line[] =
			{
				list_points[line_i - 1],
				list_points[line_i]
			};
			window.draw(line, 2, sf::Lines);
		}
		
		window.display();
		
		//===================================================//
		//delete[] PSD;
		//delete[] DEV;
	}
	recorder.stop();
	return 0;
}