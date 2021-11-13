#include"utils_tool.h"

int getRandomInt(int rangeMax) {
	return round(rangeMax * rand() / (RAND_MAX + 1));
}

double get_random_double(int rangeMax) {
	return 1.0 * rangeMax * rand() / (RAND_MAX + 1);
}

void print_vector_int(vector<int> v) {
	for (int i = 0; i < v.size(); i++) {
		cout << v[i] << " ";
	}
	cout << endl;
}

