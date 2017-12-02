#ifndef PLOT_H
#define PLOT_H

//#include <SFML/Graphics.hpp>
#include "Axis.hpp"

template<class T>
class Plot {

public:
	Plot(sf::RenderWindow* window, std::vector<T> X, std::vector<T> Y, sf::Font font);
	void axis(bool display);
	void xticks(std::vector<float> positions, std::vector<std::string> names, float tick_size = 2.0);
	void yticks(std::vector<float> positions, std::vector<std::string> names, float tick_size = 2.0);
	void xlabel(std::string label);
	void ylabel(std::string label);
	void move(float x, float y);
	void set_size(float delta_x, float delta_y);
	void show();

private:
	sf::RenderWindow* _window;
	std::vector<T> _X;
	std::vector<T> _Y;
	std::vector<float> _limits;
	Axis _xaxis;
	Axis _yaxis;
	bool _axis_display;
};


template<class T>
Plot<T>::Plot(sf::RenderWindow* window, std::vector<T> X, std::vector<T> Y, sf::Font font) {
	_window = window;
	_X = X;
	_Y = Y;
	float x_0 = 0.0;
	float y_0 = 0.0;
	float delta_x = 1000;
	float delta_y = 500;

	float max_X = *(std::max_element(std::begin(X), std::end(X)));
	float min_X = *(std::min_element(std::begin(X), std::end(X)));
	float max_Y = 0.5;// *(std::max_element(std::begin(Y), std::end(Y)));
	float min_Y = *(std::min_element(std::begin(Y), std::end(Y)));

	std::vector<float> xposition = { x_0 + 10, 
		                             x_0 + delta_x + 10,
		                             delta_y + y_0 + 10,
		                             delta_y + y_0 + 10 };

	Axis xaxis(window, xposition, font, min_X - 0.02*(max_X - min_X), max_X + 0.02*(max_X - min_X), "xlabel");
	_xaxis = xaxis;

	std::vector<float> yposition = { x_0 + 10, 
		                             x_0 + 10,
									 y_0 + 10, 
		                             y_0 + delta_y + 10 };
	Axis yaxis(window, yposition, font, min_Y - 0.02*(max_Y - min_Y), max_Y + 0.02*(max_Y - min_Y), "ylabel");
	_yaxis = yaxis;

	_axis_display = false;
}

template<class T>
void Plot<T>::move(float x, float y) {
	_xaxis.move(x, y);
	_yaxis.move(x, y);
}

template<class T>
void Plot<T>::set_size(float delta_x, float delta_y) {
	_delta_x = delta_x;
	_delta_y = delta_y;
}

template<class T>
void Plot<T>::axis(bool display)
{
	_axis_display = display;
}

template<class T>
void Plot<T>::xticks(std::vector<float> positions, std::vector<std::string> names, float tick_size) {
	_xaxis.ticks(positions, names, tick_size);
}

template<class T>
void Plot<T>::yticks(std::vector<float> positions, std::vector<std::string> names, float tick_size) {
	_yaxis.ticks(positions, names, tick_size);
}

template<class T>
void Plot<T>::xlabel(std::string label)
{
	_xaxis.label(label);
}

template<class T>
void Plot<T>::ylabel(std::string label)
{
	_yaxis.label(label);
}

template<class T>
void Plot<T>::show() {
	//std::cout << "aaa" << std::endl;
	float x_min = _xaxis.get_start();
	float y_min = _yaxis.get_start();
	float x_max = _xaxis.get_stop();
	float y_max = _yaxis.get_stop();
	//std::cout << " x_min " << x_min << " x_max " << x_max << " y_min " << y_min << " y_max " << y_max << std::endl;
	float x_range = x_max - x_min;
	float y_range = y_max - y_min;
	float x_0 = _xaxis.get_origin()[0];
	float y_0 = _xaxis.get_origin()[1];
	//std::cout << " x_range " << x_range << " y_range " << y_range << " x_0 " << x_0 << " y_0 " << y_0 << std::endl;
	
	float height = _yaxis.get_length();
	float length = _xaxis.get_length();
	//std::cout << height << "  ---  " << length << std::endl;

	std::vector<sf::Vertex> list_points;
	for (unsigned int i = 0; i < _X.size(); ++i)
	{
		float x = (std::max(x_min,std::min(x_max, (float)_X[i])) - x_min) / x_range * length + x_0;
		float y = -(std::max(y_min,std::min(y_max, (float)_Y[i])) - y_min) / y_range * height + y_0;
		//std::cout << x << "  aaa  " << y << std::endl;

		sf::Vertex point(sf::Vector2f(x, y));
		list_points.push_back(point);
	}

	for (unsigned int line_i = 1; line_i < list_points.size(); ++line_i)
	{
		sf::Vertex line[] =
		{
			list_points[line_i - 1],
			list_points[line_i]
		};
		sf::Color color(100,100,255,100);
		line[0].color = color;
		line[1].color = color;
		_window->draw(line, 2, sf::Lines);
	}

	if (_axis_display)
	{
		_xaxis.show();
		_yaxis.show();
	}
}
#endif