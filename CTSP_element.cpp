#include"CTSP_element.h"


Individual::Individual(int k) {
	this->salesman_num = k;
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

void Individual::update_by_dis() {
	long long int total = 0;
	long long int temp = 0;
	long long int max = 0;
	for (int i = 0; i < this->salesman_num; i++) {
		temp = this->salesman_distance[i];
		max = max > temp ? max : temp;
		total += temp;
	}
	this->total_distance = total;
	this->max_salesman_distance = max;
	this->fitness = 1.0 / (1.0 + total);
}



int dist(Point *a, Point *b)
{
	return int(sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2)) + 0.5);
}

