#pragma once
#include<iostream>
#include<vector>
#include<cmath>
#include<ctime>
#include<fstream>
using namespace std;



struct Point {
	int id;
	double x, y;
};

class Individual {
public:
	Individual(int k);
	double fitness;
	long long int total_distance;
	long long int max_salesman_distance;
	vector<vector<int>> individual_cities;
	vector<int> salesman_distance;
	void print_cities();
	void print_distance();
	void assign_cities(vector<vector<int>> cities);
};

int dist(Point *a, Point *b);

