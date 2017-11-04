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

	unsigned int height = 800;
	unsigned int length = 800;
	unsigned int margin = 20;

	sf::RenderWindow window(sf::VideoMode(height, length), "FFT");
	sf::Font font;
	font.loadFromFile("MonospaceTypewriter.ttf");

	double sample_rate = 20000.0;
	double sample_time = 4096.0 / 20000.0;

	MusicRecorder recorder(sample_rate, sample_time, 1.0 / 10.0);
	
	recorder.start((unsigned int)sample_rate);

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
		std::vector<double> X; X.resize(84);
		std::vector<double> Y; Y.resize(84);
		for (size_t i = 0; i < 84; i++) {
			Y[i] = freqs[i];
			X[i] = (double)i;
		}
			
		std::vector<sf::Vertex> list_points = plot(Y, X, "title", 84, 2);

		for (unsigned int line_i = 1; line_i < list_points.size(); ++line_i)
		{
			sf::Vertex line[] =
			{
				list_points[line_i - 1],
				list_points[line_i]
			};
			window.draw(line, 2, sf::Lines);
		}
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