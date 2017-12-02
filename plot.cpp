/*#include "Plot.hpp"

template<class T>
Plot<T>::Plot(sf::RenderWindow* window, std::vector<T> X, std::vector<T> Y) {
	_window = window;
	_X = X;
	_Y = Y;
	_axis_display = false;
}

template<class T>
void Plot<T>::set_origin(float x, float y) {
	_x_0 = x;
	_y_0 = y;
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
void Plot<T>::show() {

	float x_min = xaxis().start();
	float y_min = yaxis().start();
	float x_max = xaxis().stop();
	float y_max = yaxis().stop();
	float x_range = x_max - x_min;
	float y_range = y_max - y_min;
	float x_0 = x_axis.origin()[0];
	float y_0 = x_axis.origin()[1];

	unsigned int height = _yaxis.get_length();
	unsigned int length = _xaxis.get_length();

	std::vector<sf::Vertex> list_points;
	for (unsigned int point_i = 0; point_i < X.size(); ++point_i)
	{
		float x = (min(x_min(max(x_max, X[i]))) - x_min) / x_range * length + x_0;
		float y = height - (min(y_min(max(y_max, Y[i]))) - y_min) / x_range * height - y_0;

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
		sf::Color color(255.0 * (1.0 - k / (3.0*past)),
			255.0 * k / past,
			255.0 * 0.0,
			255.0 * (past + k) / (2.0*past));
		line[0].color = color;
		line[1].color = color;
		window.draw(line, 2, sf::Lines);
	}

	if (_axis_display) 
	{
		_xaxis.show()
		_yaxis.show()
	}

}
*/