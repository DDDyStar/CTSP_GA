#include<iostream>
#include<string>
#include"CTSP_element.h"
#include"utils_tool.h"
#include<vector>
#include<algorithm>
using namespace std;

string DataName = "eil";
int CityNum = 21;
int SalesmanNum = 4;
int** CitiesDistance;
int** ColorMatrix;
vector<vector<int>> CityColor(CityNum, vector<int>());
vector<Point*> PointSet;
vector<Individual*> Population;
string datapath = "D:/cpp_ws/CTSP_GA/" + DataName + to_string(CityNum) + ".txt";
string colorpath = "D:/cpp_ws/CTSP_GA/" + DataName + to_string(CityNum) + "_" + to_string(SalesmanNum) + ".txt";

int PopulationNum = 30;
int InitAlgo = 1;
double SA_Probability = 0.3;

void print_dist_matrix() {
    for (int i = 0; i < CityNum; i++) {
        for (int j = 0; j < CityNum; j++) {
            cout << CitiesDistance[i][j] << " ";
        }
        cout << endl;
    }
}

void print_color_matrix() {
    for (int i = 0; i < SalesmanNum + 1; i++) {
        for (int j = 0; j < CityNum; j++) {
            cout << ColorMatrix[i][j] << " ";
        }
        cout << endl;
    }
}

void print_city_color(int n = 10) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < CityColor.at(i).size(); j++) {
            cout << CityColor.at(i).at(j) << " ";
        }
        cout << endl;
    }
}

Point* create_p(int id, double x, double y) {
    Point* p = new Point;
    p->id = id;
    p->x = x;
    p->y = y;
    return p;
}

void read_data() {
    CitiesDistance = new int* [CityNum];
    for (int i = 0; i < CityNum; i++) {
        CitiesDistance[i] = new int[CityNum];
    }

    PointSet.reserve(1000);
    PointSet.clear();

    ifstream cityData;
    cityData.open(datapath, ios::in);

    int id;
    double x, y;
    //for (int i = 1; i < CityNum; i++) {
    while (!cityData.eof()){
        cityData >> id;
        cityData >> x;
        cityData >> y;
        PointSet.push_back(create_p(id, x, y));
    }
    cityData.close();

    cout << PointSet.size() << endl;

    for (int i = 0; i < CityNum; i++) {
        for (int j = 0; j < CityNum; j++) {
            if (j < i) {
                CitiesDistance[i][j] = CitiesDistance[j][i];
            }
            else if (j == i) {
                CitiesDistance[i][j] = 0;
            }
            else {
                CitiesDistance[i][j] = dist(PointSet.at(i), PointSet.at(j));
            }
        }
    } 
}

void read_color() {
    ColorMatrix = new int* [SalesmanNum + 1];
    for (int i = 0; i < SalesmanNum+1; i++) {
        ColorMatrix[i] = new int[CityNum];
    }

    ifstream colorData;
    int temp_color;
    colorData.open(colorpath, ios::in);

    for (int i = 0; i < SalesmanNum+1; i++) {
        for (int j = 0; j < CityNum; j++) {
            colorData >> temp_color;
            ColorMatrix[i][j] = temp_color;
        }
    }
    
    for (int i = 0; i < CityNum; i++) {
        for (int j = 1; j < SalesmanNum + 1; j++) {
            if (ColorMatrix[j][i] == 1) {
                CityColor.at(i).push_back(j - 1);
            }
        }
    }
}

vector<vector<int>> initCitiesSeq() {
    vector<vector<int>> result;
    result.resize(SalesmanNum);
    for (int i = SalesmanNum; i < CityNum; i++) {
        int tempIndex = getRandomInt(CityColor[i].size());
        result[CityColor[i][tempIndex]].push_back(i);
    }
    for (int j = 0; j < SalesmanNum; j++) {
        random_shuffle(result[j].begin(), result[j].end());
        result[j].insert(result[j].begin(), j);
    }
    return result;
}

long long int getIndividualSalesDistance(vector<int> seq){
    if (seq.empty()) {
        return 0;
    }
    long long int result = 0;
    for (int i = 0; i < seq.size() - 1; i++) {
        result += CitiesDistance[seq[i]][seq[i + 1]];
    }
    
    result += CitiesDistance[seq[seq.size() - 1]][seq[0]];
    
    return result;
}

void updateIndividual(int populationIndex) {
    long long int total = 0;
    long long int temp = 0;
    long long int max = 0;
    for (int i = 0; i < SalesmanNum; i++) {
        temp = getIndividualSalesDistance(Population[populationIndex]->individual_cities[i]);
        Population[populationIndex]->salesman_distance[i] = temp;
        max = max > temp ? max : temp;
        total += temp;
    }
    Population[populationIndex]->total_distance = total;
    Population[populationIndex]->max_salesman_distance = max;
    Population[populationIndex]->fitness = 1.0 / (1.0 + total);
}

void updateIndividualByCost(Individual &ind) {

}

void greedy(vector<int> &seq) {
    int min, minIndex, tempIndex;
    for (int i = 0; i < seq.size()-1; i++) {
        for (int j = i+1; j < seq.size(); j++) {
            if (j == i+1) {
                min = CitiesDistance[seq[i]][seq[j]];
                minIndex = j;
            }
            else if (min > CitiesDistance[seq[i]][seq[j]]) {
                min = CitiesDistance[seq[i]][seq[j]];
                minIndex = j;
            }
        }
        tempIndex = seq[i + 1];
        seq[i + 1] = seq[minIndex];
        seq[minIndex] = tempIndex;
    }
}

void createPopulation() {
    for (int i = 0; i < PopulationNum; i++) {
        Population.push_back(new Individual(SalesmanNum));
    }
}

void clearPopulation() {
    for (int i = 0; i < PopulationNum; i++) {
        delete Population[i];
    }
    Population.clear();
}

void populationInit() {
    if (!Population.empty()) {
        clearPopulation();
    }
    createPopulation();
    for (int i = 0; i < Population.size(); i++) {
        Population[i]->individual_cities.clear();
        Population[i]->individual_cities = initCitiesSeq();
        if (InitAlgo == 1) {
            for (int j = 0; j < SalesmanNum; j++) {
                greedy(Population[i]->individual_cities[j]);
            }
        }
        updateIndividual(i);
    }
}

Individual getBestIndividual() {
    double max = 0;
    int maxIndex = 0;
    for (int i = 0; i < Population.size(); i++) {
        if (max < Population[i]->fitness) {
            max = Population[i]->fitness;
            maxIndex = i;
        }
    }
    return *(Population[maxIndex]);
}

int roulette() {
    vector<double> chooseTable;
    double sum = 0.0;
    for (int i = 0; i < Population.size(); i++) {
        chooseTable.push_back(Population[i]->fitness);
        sum = sum + Population[i]->fitness;
    }
    double randomCard = rand() / double(RAND_MAX);
    int indexResult = -1;
    double temp, sumTemp = 0.0;

    for (int i = 0; i < Population.size(); i++) {
        sumTemp += chooseTable[i];
        temp = sumTemp / sum;
        if (temp >= randomCard) {
            indexResult = 1;
            break;
        }
    }
    return indexResult;
}

void swap(vector<int>& seq, int index1, int index2) {
    vector<int> temp = seq;
    int length = index2 - index1;
    for (int i = index1; i <= index2; i++) {
        seq[i] = temp[index2-(i-index1)];
    }
}

void two_opt(vector<int>& seq) {
    if (seq.size() < 3) {
        return;
    }
    vector<int> R_part;
    int changeResult = 0;
    int t = getRandomInt(seq.size());
    R_part.insert(R_part.begin(), seq.begin(), seq.begin() + t);
    R_part.insert(R_part.begin(), seq.begin() + t, seq.end());
    seq = R_part;
    for (int i = 1; i < seq.size() - 2; i++) {
        for (int j = i + 1; j < seq.size() - 1; j++) {
            changeResult = -CitiesDistance[seq[i - 1]][seq[i]] - CitiesDistance[seq[j]][seq[j + 1]]
                + CitiesDistance[seq[i - 1]][seq[j]] + CitiesDistance[seq[i]][seq[j + 1]];
            if (changeResult < 0) {
                swap(seq, i, j);
                j = j + 1;
            }
        }
    }
}

void swap_cities(vector<int>& vec1, int index1, int index2) {
    int temp = vec1[index1];
    vec1[index1] = vec1[index2];
    vec1[index2] = temp;
}

void swap_cities(vector<int> &vec1, vector<int> &vec2, int index1, int index2) {
    int temp = vec1[index1];
    vec1[index1] = vec2[index2];
    vec2[index2] = temp;
}

void insert_opt(vector<int>& seq, string mode="hc") {
    // parameter: mode
    //      value: "hc" or "sa"
    if (seq.size() < 3) {
        return;
    }
    vector<int> R_part;
    int t = getRandomInt(seq.size());
    R_part.insert(R_part.begin(), seq.begin(), seq.begin() + t);
    R_part.insert(R_part.begin(), seq.begin() + t, seq.end());
    seq = R_part;
    int changeResult = 0.0;

    for (int i = 1; i < seq.size() - 3; i++) {
        for (int j = i + 2; j < seq.size() - 1; j++) {
            changeResult = -CitiesDistance[seq[i - 1]][seq[i]] - CitiesDistance[seq[i]][seq[i + 1]]
                - CitiesDistance[seq[j - 1]][seq[j]] - CitiesDistance[seq[j]][seq[j + 1]]
                + CitiesDistance[seq[i - 1]][seq[j]] + CitiesDistance[seq[j]][seq[i + 1]]
                + CitiesDistance[seq[j - 1]][seq[i]] + CitiesDistance[seq[i]][seq[j + 1]];
            if (changeResult < 0) {
                swap_cities(seq, i, j);
                j = j + 2;
            }
            /*else { // unfinished, SA mechanism
                if (mode == "sa") {
                    if (get_random_double(1) < SA_Probability) {

                    }
                }
            }*/
        }
    }
}

void mutation(Individual& individual) {
    // only change colors
    if (CityNum < 3) {
        return;
    }
    vector<int> seq, color, tx;
    // double chromosome seq, color, seq represents the tour
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < CityNum; j++) {
            for (int k = 0; k < SalesmanNum; k++) {
                if (j < individual.individual_cities[k].size()) {
                    seq.push_back(individual.individual_cities[k][j]);
                    color.push_back(k);
                }
            }
        }
        int locationL = getRandomInt(seq.size() / 2);
        int locationR = locationL + 1 + getRandomInt(seq.size() / 2);
        swap_cities(seq, locationL, locationR);

        individual.individual_cities.clear();
        individual.individual_cities.resize(SalesmanNum);

        int colorTemp1 = color[locationL];
        if (count(CityColor[seq[locationL]].begin(), CityColor[seq[locationL]].end(), colorTemp1) == 0) {
            colorTemp1 = CityColor[seq[locationL]][getRandomInt(CityColor[seq[locationL]].size())];
            color[locationL] = colorTemp1;
        }

        int colorTemp2 = color[locationR];
        if (count(CityColor[seq[locationR]].begin(), CityColor[seq[locationR]].end(), colorTemp2) == 0) {
            colorTemp2 = CityColor[seq[locationR]][getRandomInt(CityColor[seq[locationR]].size())];
            color[locationR] = colorTemp2;
        }

        for (int i = 0; i < seq.size(); i++) {
            individual.individual_cities[color[i]].push_back(i);
        }
        individual.individual_cities
    }
}


void test1() {
    Individual ind1(SalesmanNum);
    Population.push_back(&ind1);
    ind1.individual_cities = initCitiesSeq();
    updateIndividual(0);
    ind1.print_distance();
    for (int i = 0; i < Population[0]->individual_cities.size(); i++) {
        greedy(Population[0]->individual_cities[i]);
    }
    updateIndividual(0);
    ind1.print_distance();

}

void test2() {
    populationInit();
    for (int i = 0; i < Population.size(); i++) {
        cout << "-------" << endl;
        Population[i]->print_distance();
        Population[i]->print_cities();
    }
}

void test3() {
    vector<int> testSeq;
    for (int i = 0; i < 10; i++) {
        testSeq.push_back(i);
    }
    for (int i = 0; i < testSeq.size(); i++) {
        cout << testSeq[i] << " ";
    }
    cout << endl;
    swap(testSeq, 1, 2);
    for (int i = 0; i < testSeq.size(); i++) {
        cout << testSeq[i] << " ";
    }
    cout << endl;

    string s1 = "dyx";
    cout << (s1 == "dyx") << endl;
    for (int i = 0; i < 5; i++) {
        cout << get_random_double(1) << endl;
    }
}

void test4() {
    Individual ind1(SalesmanNum);
    int index = 0;
    vector<vector<int>> cities_seq;
    vector<int> temp_seq;
    for (int i = 0; i < SalesmanNum; i++) {
        temp_seq.clear();
        for (int j = 0; j < 5; j++) {
            temp_seq.push_back(index);
            index += 1;
        }
        if (i == 3) {
            temp_seq.push_back(index);
        }
        cities_seq.push_back(temp_seq);
    }
    ind1.assign_cities(cities_seq);
    ind1.print_cities();
    mutation(ind1);
}

void test_rand() {
    for (int i = 0; i < 10; i++) {
        cout << getRandomInt(2) << " ";
    }
    
}

int main() {
	read_data();
    read_color();
    test_rand();
	//print_dist_matrix();
    //print_color_matrix();
    //print_city_color(20);`
	return 0;
}