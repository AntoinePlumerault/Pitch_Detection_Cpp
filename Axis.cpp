#include "Axis.hpp"

Axis::Axis() 
{

}

Axis::Axis(sf::RenderWindow* window, std::vector<float> position, sf::Font font, float start, float stop, std::string label) {
	_window = window;
	_position = position;
	_start = start;
	_stop = stop;
	_label = label;
	_font = font;
}

float Axis::get_start() const
{
	return _start;
}

float Axis::get_stop() const
{
	return _stop;
}

float Axis::get_length() const
{
	return sqrt(pow((_position[1] - _position[0]), 2) + pow((_position[3] - _position[2]), 2));
}

std::vector<float> Axis::get_origin() const
{
	std::vector<float> origin = {_position[0], _position[2]};
	return origin;
}

void Axis::move(float x, float y)
{
	_position[0] += x;
	_position[1] += x;
	_position[2] += y;
	_position[3] += y;

}

void Axis::ticks(std::vector<float> positions, std::vector<std::string> label, float tick_size) {
	_ticks = positions;
	_ticks_label = label;
	_ticks_size = tick_size;
}

void Axis::label(std::string label)
{
	_label = label;
}

void  Axis::show() const
{
	// Print Axis
	sf::Vertex axis[] =
	{
		sf::Vertex(sf::Vector2f(_position[0], _position[2])),
		sf::Vertex(sf::Vector2f(_position[1], _position[3]))
	};
	sf::Color color(255, 255, 255, 255);
	axis[0].color = color;
	axis[1].color = color;
	_window->draw(axis, 2, sf::Lines);

	float orientation = atan((_position[3] - _position[2]) / (_position[1] - _position[0]));
	if (_position[1] - _position[0] < 0) {
		orientation += (float) PI;
	}
	
	float length = get_length();
	float scale = _stop - _start;

	float delta_x = (float)-(_position[3] - _position[2]) / length * _ticks_size / 2.0;
	float delta_y = (float)-(_position[1] - _position[0]) / length * _ticks_size / 2.0;

	for (size_t i = 0; i < _ticks.size(); ++i) {

		// Print ticks
		float x = _position[1] * ((_ticks[i] - _start) / scale) + _position[0] * (1.0 - (_ticks[i] - _start) / scale);
		float y = _position[3] * ((_ticks[i] - _start) / scale) + _position[2] * (1.0 - (_ticks[i] - _start) / scale);
		
		sf::Vertex tick[] =
		{
			sf::Vertex(sf::Vector2f(x + delta_x, y + delta_y)),
			sf::Vertex(sf::Vector2f(x - delta_x, y - delta_y))
		};
		sf::Color color(255, 255, 255, 255);
		tick[0].color = color;
		tick[1].color = color;
		_window->draw(tick, 2.0, sf::Lines);

		// Print ticks labels
		sf::Text tick_label(_ticks_label[i], _font);
		tick_label.setPosition(x - delta_x, y - delta_y);
		tick_label.setCharacterSize(10);
		tick_label.setRotation(45);
		tick_label.setStyle(sf::Text::Regular);
		_window->draw(tick_label);
	}

	// Print axis label
	sf::Text label(_label, _font);
	label.setPosition((_position[1] + _position[0]) / 2.0 - 10*delta_x , (_position[3] + _position[2]) / 2.0 - 10*delta_y);
	label.setCharacterSize(20);
	label.setStyle(sf::Text::Regular);
	_window->draw(label);
}