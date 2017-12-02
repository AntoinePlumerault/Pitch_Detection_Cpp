#ifndef AXIS_H
#define AXIS_H

#include <SFML/Graphics.hpp>
#include <vector>
#include <string>

#define PI 3.14159265359

class Axis {

public:
	Axis();
	Axis(sf::RenderWindow* window, std::vector<float> position, sf::Font font, float start = 0.0, float stop = 0.0, std::string label = "");
	float get_start() const;
	float get_stop() const;

	float get_length() const;
	std::vector<float> get_origin() const;

	void move(float x, float y);
	
	void ticks(std::vector<float> positions, std::vector<std::string> names, float tick_size = 2.0);
	void label(std::string label);
	void show() const;

private:
	sf::RenderWindow* _window;
	sf::Font _font;

	std::vector<float> _position;
	float _line_width;
	float _ticks_size;

	float _start;
	float _stop;
	
	std::vector<float> _ticks;
	std::vector<std::string> _ticks_label;
	
	std::string _label;
};

#endif
