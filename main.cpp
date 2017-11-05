#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <iostream>
#include <thread>
#include <chrono>
#include "plot.hpp"
//#include "fft.hpp"
#include "MusicRecorder.hpp"



int main()
{
	if (!sf::SoundBufferRecorder::isAvailable())
	{
		printf("error");
	}

	unsigned int height = 800+400;
	unsigned int length = 800;
	unsigned int margin = 20;

	sf::ContextSettings settings;
	settings.antialiasingLevel = 8;
	sf::RenderWindow window(sf::VideoMode(height, length), "FFT", sf::Style::Default, settings);
	sf::Font font;
	//font.loadFromFile("MonospaceTypewriter.ttf");

	double sample_rate = 20000.0;
	double sample_time = 4096.0 / 20000.0;

	MusicRecorder recorder(sample_rate, sample_time, 1.0 / 10.0);
	
	recorder.start((unsigned int)sample_rate);
	
	size_t past = 300;
	std::vector<std::vector<double>> list_Y; list_Y.resize(past);
	for (size_t i = 0; i < past; ++i)
	{
		std::vector<double> Y; Y.resize(168);
		list_Y[i] = Y;
	}
	size_t buffer_pos = 0;

	std::vector<double> X; X.resize(168);
	for (size_t i = 0; i < 168; ++i)
	{
		X[i] = (double)i;
	}

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


		//---------------------------------------------------//
		std::vector<double> freqs = recorder.getFreqs(2);

		for (size_t i = 0; i < 168; i++) 
		{
			list_Y[buffer_pos][i] = log(1.0 + freqs[i]/10.0);
		}
		
		for (size_t k = 1; k < past+1; ++k) 
		{
			std::vector<sf::Vertex> list_points = plot(list_Y[(buffer_pos + k) % past], X, "title", 168, 2.0, 2.1 * (1.0-((double)k / (double)past)));

			for (unsigned int line_i = 1; line_i < list_points.size(); ++line_i)
			{
				sf::Vertex line[] =
				{
					list_points[line_i - 1],
					list_points[line_i]
				};
				sf::Color color(127.0 ,
					            255.0 * k / past, 
					            255.0 * (1.0 - k / (3.0*past)), 
					            255.0 * (past + k) / (2.0*past));
				line[0].color = color;
				line[1].color = color;
				window.draw(line, 2, sf::Lines);
			}
		}

		buffer_pos = (buffer_pos+1) % past;
		/*
		std::vector<double> freqs = recorder.getFreqs(2);
		for (double freq : freqs)
			std::cout << freq << " - ";
		std::cout << std::endl;

		sf::Text note1(freq2note(recorder.getFreqs(2)[0], "fr") + freq2note(recorder.getFreqs(2)[1], "fr"), font);
		note1.setCharacterSize(150);
		note1.setStyle(sf::Text::Bold);
		window.draw(note1);
		*/

		//---------------------------------------------------//


		window.display();
		//===================================================//
	}
	recorder.stop();
	return 0;
}