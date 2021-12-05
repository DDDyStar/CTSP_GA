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
	Individual(int k, int id);
	double fitness;
	long long int total_distance;
	long long int max_salesman_distance;
	vector<vector<int>> individual_cities;
	vector<int> salesman_distance;
	int salesman_num;
	int ind_id = 0;
	void print_cities();
	void print_distance();
	void assign_cities(vector<vector<int>> cities);
	void update_by_dis(); // total_dis, max_dis, fit are updated according to salesman_dis
	void chromosome2seq(vector<int> &chromosome, vector<int> &color);
};

int dist(Point *a, Point *b);

