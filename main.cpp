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

	unsigned int height = 800+400+800;
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
	
	size_t past = 500;
	std::vector<std::vector<double>> list_Y; list_Y.resize(past);
	for (size_t i = 0; i < past; ++i)
	{
		std::vector<double> Y; Y.resize(168);
		list_Y[i] = Y;
	}
	size_t buffer_pos = 0;

	std::vector<std::vector<double>> list_fY; list_fY.resize(past);
	for (size_t i = 0; i < past; ++i)
	{
		std::vector<double> fY; fY.resize(168);
		list_fY[i] = fY;
	}

	std::vector<std::vector<double>> list_smooth; list_smooth.resize(past);
	for (size_t i = 0; i < past; ++i)
	{
		std::vector<double> smooth; smooth.resize(168);
		list_smooth[i] = smooth;
	}

	std::vector<std::vector<double>> list_loc_max; list_loc_max.resize(past);
	for (size_t i = 0; i < past; ++i)
	{
		std::vector<double> loc_max; loc_max.resize(168);
		list_loc_max[i] = loc_max;
	}

	std::vector<double> X; X.resize(168);
	for (size_t i = 0; i < 168; ++i)
	{
		X[i] = (double)i;
	}

	while (window.isOpen())
	{
		auto tic = std::chrono::high_resolution_clock::now();
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
			list_Y[buffer_pos][i] = log(1.0+freqs[i])/10.0;
		}

		double sigma = 2.0;
		for (size_t i = 0; i < 168; i++)
		{
			list_smooth[buffer_pos][i] = 0.0;
		}
		for (int k = -5; k < 5 + 1; ++k)
		{
			double alpha = 1.0 / (sqrt(2.0 * PI) * sigma) * exp(-pow((double)k, 2.0) / (2.0*pow(sigma, 2.0)));
			for (size_t i = 0; i < 168; i++)
			{
				if (i + k >= 0 & i + k < 168)
				{
					list_smooth[buffer_pos][i] += list_Y[buffer_pos][i + k] * alpha;
				}		
			}
		}
		
		for (size_t i = 1; i < 167; i++)
		{
			list_loc_max[buffer_pos][i] = 0.0;
			if (list_Y[buffer_pos][i] > list_Y[buffer_pos][i+1] & list_Y[buffer_pos][i] > list_Y[buffer_pos][i - 1])
			{
				list_loc_max[buffer_pos][i] = 0.05*std::max(0.0, (list_Y[buffer_pos][i]/(0.001+list_smooth[buffer_pos][i])-1.0));
				/*
				if (list_Y[buffer_pos][i] < 0.5)
					list_loc_max[buffer_pos][i] = 0.0;
				else 
					list_loc_max[buffer_pos][i] = 0.02;
				*/
			}
			
		}
		for (size_t i = 0; i < 168; ++i)
		{
			/*
			list_fY[buffer_pos][i] = 0.0;
			for (size_t k = 0; k < 4; ++k) 
			{
				if (round(i + 12.0 * log((double)k) / log(2.0)) < 168)
				{
					list_fY[buffer_pos][i] += 1.0/(1.0+k) * list_loc_max[buffer_pos][round(i + 12.0 * log((double)k) / log(2.0))];
				}
			}
			list_fY[buffer_pos][i] = log(1.0 + list_fY[buffer_pos][i]/5.0)/2.0;
			*/
			//list_fY[buffer_pos][i] = 0.5*1.0/2.0 * (abs(list_Y[buffer_pos][i] - list_Y[(buffer_pos-1)%past][i]) + (list_Y[buffer_pos][i] - list_Y[(buffer_pos - 1) % past][i])) + 0.5*list_fY[(buffer_pos-1)%past][i];
		}

		for (size_t k = 1; k < past+1; ++k) 
		{
			double delta = 2.1 * (1.0 - ((double)k / (double)past));
			std::vector<sf::Vertex> list_points = plot(list_smooth[(buffer_pos + k) % past], X, "title", delta, 168, -delta, 2.0);

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

			list_points = plot(list_loc_max[(buffer_pos + k) % past], X, "title", delta + 4.0, 168, -delta, 2.0);

			for (unsigned int line_i = 1; line_i < list_points.size(); ++line_i)
			{
				sf::Vertex line[] =
				{
					list_points[line_i - 1],
					list_points[line_i]
				};
				sf::Color color(255.0 * (1.0 - k / (3.0*past)),
								255.0 * k / past,
								255.0 * 0.0,
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
		auto toc = std::chrono::high_resolution_clock::now();
		double elapsed_time = 1.e-6*std::chrono::duration_cast<std::chrono::nanoseconds>(toc - tic).count();
		std::this_thread::sleep_for((std::chrono::milliseconds(size_t(1000.0/30.0 - elapsed_time))));
		

	}
	recorder.stop();
	return 0;
}