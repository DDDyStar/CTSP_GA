#include"CTSP_element.h"


Individual::Individual(int k) {
	this->fitness = 0;
	this->total_distance = 0;
	this->max_salesman_distance = 0;
	this->individual_cities;
	for (int i = 0; i < k; i++) {
		this->salesman_distance.push_back(0.);
	}
	
}

void Individual::print_cities() {
	for (int i = 0; i < this->individual_cities.size(); i++) {
		for (int j = 0; j < this->individual_cities[i].size(); j++) {
			cout << this->individual_cities[i][j] << " ";
		}
		cout << endl;
	}
}

void Individual::print_distance() {
	for (int i = 0; i < this->salesman_distance.size(); i++) {
		cout << this->salesman_distance[i] << " ";
	}
	cout << endl;
}

void Individual::assign_cities(vector<vector<int>> cities) {
	this->individual_cities = cities;
}




int dist(Point *a, Point *b)
{
	return int(sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2)) + 0.5);
}

