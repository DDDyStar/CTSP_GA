#include"utils_tool.h"

int getRandomInt(int rangeMax) {
	return round(rangeMax * rand() / (RAND_MAX + 1));
}

double get_random_double(int rangeMax) {
	return 1.0 * rangeMax * rand() / (RAND_MAX + 1);
}

