#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <iostream>
#include <thread>
#include <chrono>
#include "fft3.hpp"
#include "plot.hpp"
#include "MyRecorder.hpp"

int main()
{
	if (!sf::SoundBufferRecorder::isAvailable())
	{
		printf("error");
	}

	sf::Time interval = sf::microseconds(100000);
	MyRecorder recorder(interval);

	unsigned int height = 800;
	unsigned int length = 800;
	unsigned int margin = 20;

	sf::RenderWindow window(sf::VideoMode(height, length), "title");
	//recorder.start();
	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}		

		recorder.start();
		std::this_thread::sleep_for(std::chrono::microseconds(20000));
		recorder.stop();

		//const sf::SoundBuffer& buffer = recorder.getBuffer();
		//const sf::Int16* samples = buffer.getSamples();
		//std::size_t count = buffer.getSampleCount();
		
		//recorder.start();
		//std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
		
		//std::cout << count << std::endl;

		/*CArray data;
		data.resize(count);
		for (unsigned int signal_i = 0; signal_i < count; ++signal_i) {
			Complex complex((double)(samples[signal_i]), 0.0);
			data[signal_i] = complex;
		}
		
		fft(data);
		*/
		CArray data = recorder.getFT();
		std::size_t count = recorder.getCount();

		std::vector<double> X;
		std::vector<double> Y;
		X.resize(count);
		Y.resize(count);
		
		for (unsigned int signal_i = 0; signal_i < count; ++signal_i) {
			X[signal_i] = log(1 + norm(data[signal_i]));
			Y[signal_i] = (double)signal_i;
		}
		
		std::vector<sf::Vertex> list_points = plot(X, Y, "title");
		
		window.clear();
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
		


		/*
		size_t chunk_size = 256;
		size_t sample_i = 0;
		while (sample_i+chunk_size < count) 
		{
			std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t3 - t1);

			std::chrono::microseconds time_span_ms = std::chrono::duration_cast<std::chrono::microseconds> (time_span);
			std::chrono::microseconds zero(50000);

			if ((recording_time - time_span_ms) < zero) { break; };

			CArray data_chunk;
			data_chunk.resize(chunk_size);
			for (unsigned int signal_i = 0; signal_i < chunk_size; ++signal_i) 
			{
				Complex complex((double)(samples[sample_i + signal_i]), 0.0);
				data_chunk[signal_i] = complex;
			}
			sample_i += chunk_size;
			fft(data_chunk);

			std::vector<double> X;
			std::vector<double> Y;
			X.resize(chunk_size);
			Y.resize(chunk_size);

			for (unsigned int signal_i = 0; signal_i < chunk_size; ++signal_i) {
				X[signal_i] = log(1 + norm(data_chunk[signal_i]));
				Y[signal_i] = (double)signal_i;
			}
			
			std::vector<sf::Vertex> list_points = plot(X, Y, "title");
			window.clear();
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
			// std::cout << data[0] << std::endl;
		}
		*/
		// std::cout << "ok" << std::endl;
		
		/*std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

		std::chrono::microseconds time_span_ms = std::chrono::duration_cast<std::chrono::microseconds> (time_span);
		
		std::cout << time_span_ms.count() << std::endl;
		std::this_thread::sleep_for(recording_time-time_span_ms);*/
		//recorder.stop();
	}
	
	return 0;
}