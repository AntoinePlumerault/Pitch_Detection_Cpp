#ifndef PLOT_OLD_H
#define PLOT_OLD_H

template<class T>
std::vector<sf::Vertex> plot(const std::vector<T> Y, const std::vector<T> X, const std::string title, const float x_min,const float x_range, const float y_min, const float y_range)
{
	T max_X = *(std::max_element(std::begin(X), std::end(X)));
	T min_X = *(std::min_element(std::begin(X), std::end(X)));
	T max_Y = *(std::max_element(std::begin(Y), std::end(Y)));
	T min_Y = *(std::min_element(std::begin(Y), std::end(Y)));
	
	T X_range = max_X - min_X;
	T Y_range = max_Y - min_Y;

	unsigned int height = 800;
	unsigned int length = 800;
	unsigned int margin = 20;

	float X_scale = (float)(length-2.0*margin) / x_range;
	float Y_scale = (float)(height-2.0*margin) / y_range;

	std::vector<sf::Vertex> list_points;
	for (unsigned int point_i = 0; point_i < X.size(); ++point_i)
	{
		float x = (float)(X[point_i] - x_min)*X_scale + margin + x_min*200;
		float y = (float)height - (float)(Y[point_i] - y_min)*Y_scale - margin;
		/*
		float x = (float)(X[point_i] - min_X)*X_scale + margin + delta*200;
		float y = (float)height - (float)(Y[point_i] - min_Y)*Y_scale - margin;
		*/
		sf::Vertex point (sf::Vector2f(x, y));
		list_points.push_back(point);
	}
	return list_points;
}

#endif