#include<iostream>
#include<string>
#include"CTSP_element.h"
#include"utils_tool.h"
#include<vector>
#include<algorithm>
#include<random>
using namespace std;

string DataName = "eil";
double ccp = 0.7;
double mp = 0.1;
int CityNum = 21;
int SalesmanNum = 4;
double MAX_time = 120.0;
int** CitiesDistance;
int** ColorMatrix;
vector<vector<int>> CityColor(CityNum, vector<int>());
vector<Point*> PointSet;
vector<Individual> Population;
string datapath = "D:/cpp_ws/CTSP_GA/" + DataName + to_string(CityNum) + ".txt";
string colorpath = "D:/cpp_ws/CTSP_GA/" + DataName + to_string(CityNum) + "_" + to_string(SalesmanNum) + ".txt";

int PopulationNum = 30;
int InitAlgo = 1;
double SA_Probability = 0.3;
Individual bestIndividual(SalesmanNum);
int localSearchChoice = 1;

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
        temp = getIndividualSalesDistance(Population[populationIndex].individual_cities[i]);
        Population[populationIndex].salesman_distance[i] = temp;
        max = max > temp ? max : temp;
        total += temp;
    }
    Population[populationIndex].total_distance = total;
    Population[populationIndex].max_salesman_distance = max;
    Population[populationIndex].fitness = 1.0 / (1.0 + total);
}

void updateIndividual(Individual &ind) {
    long long int total = 0;
    long long int temp = 0;
    long long int max = 0;
    for (int i = 0; i < SalesmanNum; i++) {
        temp = getIndividualSalesDistance(ind.individual_cities[i]);
        ind.salesman_distance[i] = temp;
        max = max > temp ? max : temp;
        total += temp;
    }
    ind.total_distance = total;
    ind.max_salesman_distance = max;
    ind.fitness = 1.0 / (1.0 + total);
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
        Population.push_back(Individual(SalesmanNum));
    }
}

//void clearPopulation() {
//    for (int i = 0; i < PopulationNum; i++) {
//        delete Population[i];
//    }
//    Population.clear();
//}

void populationInit() {
    if (!Population.empty()) {
        Population.clear();
    }
    createPopulation();
    for (int i = 0; i < Population.size(); i++) {
        Population[i].individual_cities.clear();
        Population[i].individual_cities = initCitiesSeq();
        if (InitAlgo == 1) {
            for (int j = 0; j < SalesmanNum; j++) {
                greedy(Population[i].individual_cities[j]);
            }
        }
        updateIndividual(i);
    }
}

Individual getBestIndividual() {
    double max = 0;
    int maxIndex = 0;
    for (int i = 0; i < Population.size(); i++) {
        if (max < Population[i].fitness) {
            max = Population[i].fitness;
            maxIndex = i;
        }
    }
    return Population[maxIndex];
}

int roulette(vector<Individual> population) {
    vector<double> chooseTable;
    double sum = 0.0;
    for (int i = 0; i < population.size(); i++) {
        chooseTable.push_back(population[i].fitness);
        sum = sum + population[i].fitness;
    }
    double randomCard = rand() / double(RAND_MAX);
    int indexResult = -1;
    double temp, sumTemp = 0.0;

    for (int i = 0; i < population.size(); i++) {
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

void swap(vector<int>& seq1, vector<int>& seq2, int index1, int index2) {
    int temp;
    for (int i = index1; i <= index2; i++) {
        temp = seq1[i];
        seq1[i] = seq2[i];
        seq2[i] = temp;
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
        
        for (int i = 0; i < SalesmanNum; i++) {
            individual.salesman_distance[i] = getIndividualSalesDistance(individual.individual_cities[i]);
        }
        individual.update_by_dis();
    }
}

Individual crossover(Individual& ind1, Individual& ind2) {

    vector<int> seq1, seq2, color1, color2;
    // get two chromesome
    for (int j = 0; j < CityNum; j++) {
        for (int k = 0; k < SalesmanNum; k++) {
            
            if (j < ind1.individual_cities[k].size()) {
                seq1.push_back(ind1.individual_cities[k][j]);
                color1.push_back(k);
            }

            if (j < ind2.individual_cities[k].size()) {
                seq2.push_back(ind2.individual_cities[k][j]);
                color2.push_back(k);
            }
        }
    }

    cout << "---seq1---" << endl;
    print_vector_int(seq1);
    print_vector_int(color1);
    cout << "---seq2---" << endl;
    print_vector_int(seq2);
    print_vector_int(color2);
    

    int left = getRandomInt(seq1.size() / 2);
    int right = left + 1 + getRandomInt(seq1.size() / 2);
    cout << "change position:" << left << " " << right << endl;
    vector<int> seg1, seg2, colorSeg1, colorSeg2;
    for (int i = left; i <= right; i++) {
        seg1.push_back(seq1[i]);
        seg2.push_back(seq2[i]);
        colorSeg1.push_back(color1[i]);
        colorSeg2.push_back(color2[i]);
    }
    // find the element appeared in the cross part of seq2
    vector<int> con1, con2;
    int num1 = 0, num2 = 0;
    for (int i = 0; i < left; i++) {
        num1 = count(seq2.begin() + left, seq2.begin() + right + 1, seq1[i]);
        num2 = count(seq1.begin() + left, seq1.begin() + right + 1, seq2[i]);
        if (num1 == 1) con1.push_back(i);
        if (num2 == 1) con2.push_back(i);
    }

    for (int i = right+1; i < CityNum; i++) {
        num1 = count(seq2.begin() + left, seq2.begin() + right + 1, seq1[i]);
        num2 = count(seq1.begin() + left, seq1.begin() + right + 1, seq2[i]);
        if (num1 == 1) con1.push_back(i);
        if (num2 == 1) con2.push_back(i);
    }

    // crossover
    swap(seq1, seq2, left, right);
    swap(color1, color2, left, right);
    
    
    // check
    int curCity, curColor;
    vector<int>::iterator corPos, segPos;
    cout << "checking..." << endl;

    for (int i = 0; i < con1.size(); i++) {
        curCity = seq1[con1[i]];
        curColor = color1[con1[i]];
        cout << curCity << endl;
        for (int j = left; j <= right; j++) {
            segPos = find(seg2.begin(), seg2.end(), seq1[con1[i]]);
            seq1[con1[i]] = seg1[segPos - seg2.begin()];
            color1[con1[i]] = colorSeg1[segPos - seg2.begin()];
            cout << "criterion" << endl;
            if (count(seg2.begin(), seg2.end(), seq1[con1[i]]) == 0) {
                corPos = find(seq2.begin(), seq2.begin() + left, seq1[con1[i]]);
                if (corPos != seq2.begin() + left) {
                    seq2[corPos - seq2.begin()] = curCity;
                    color2[corPos - seq2.begin()] = curColor;
                }
                else {
                    corPos = find(seq2.begin() + right + 1, seq2.end() + left, seq1[con1[i]]);
                    seq2[corPos - seq2.begin()] = curCity;
                    color2[corPos - seq2.begin()] = curColor;
                }
                break;
            }
        }
    }
    cout << "---seq1---" << endl;
    print_vector_int(seq1);
    print_vector_int(color1);
    cout << "---seq2---" << endl;
    print_vector_int(seq2);
    print_vector_int(color2);


    ind1.chromosome2seq(seq1, color1);
    updateIndividual(ind1);
    ind2.chromosome2seq(seq2, color2);
    updateIndividual(ind2);

    cout << ind1.fitness << " fit " << ind2.fitness << endl;

    if (ind1.fitness > ind2.fitness) {
        return ind1;
    }
    else {
        return ind2;
    }
}


void localSearch(Individual &ind) {
    for (int i = 0; i < ind.individual_cities.size(); i++) {
        switch (localSearchChoice)
        {
        case 1:
            insert_opt(ind.individual_cities[i]);
        default:
            break;
        }
    }
}
double GetTime()
{
    return (double)clock() / CLOCKS_PER_SEC;
}

int GA() {
    vector<Individual> nextGenerationPopulation;
    Individual ind = bestIndividual;
    int generationIndex = 1;
    double CPU_time, LastTime = GetTime();
    while (true) {
        nextGenerationPopulation.clear();
        localSearch(ind);
        updateIndividual(ind);
        nextGenerationPopulation.push_back(ind);
        while (nextGenerationPopulation.size() != Population.size()) {
            int A = roulette(Population);
            Individual suc = Population[A];
            if (rand() / RAND_MAX < ccp) {
                int B = roulette(Population);
                Individual rouletteSelectIndividualB = Population[B];
                while (A == B) {
                    B = roulette(Population);
                    Individual rouletteSelectIndividualB = Population[B];
                }
                suc = crossover(suc, rouletteSelectIndividualB);
            }
            if (rand() / RAND_MAX < mp) {
                mutation(suc);
            }
            nextGenerationPopulation.push_back(suc);
        }
        Population.clear();
        Population = nextGenerationPopulation;

        ind = getBestIndividual();

        if (ind.fitness > bestIndividual.fitness) {
            bestIndividual = ind;
        }

        CPU_time = fabs(GetTime() - LastTime);
        if (CPU_time >= MAX_time) { break; }
        else { generationIndex += 1; }
    }
    return generationIndex;
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

void test5() {
    vector<int> seq1, seq2;
    for (int i = 0; i < CityNum; i++) {
        seq1.push_back(i);
        seq2.push_back(i);
    }

    shuffle(seq2.begin(), seq2.end(), default_random_engine(10));

    int index = 0;
    vector<vector<int>> cities_seq1, cities_seq2;
    vector<int> temp_seq1, temp_seq2;
    for (int i = 0; i < SalesmanNum; i++) {
        temp_seq1.clear();
        temp_seq2.clear();
        for (int j = 0; j < 5; j++) {
            temp_seq1.push_back(seq1[index]);
            temp_seq2.push_back(seq2[index]);
            index += 1;
        }
        if (i == 3) {
            temp_seq1.push_back(seq1[index]);
            temp_seq2.push_back(seq2[index]);
        }
        cities_seq1.push_back(temp_seq1);
        cities_seq2.push_back(temp_seq2);
    }

    Individual ind1(SalesmanNum), ind2(SalesmanNum);
    ind1.assign_cities(cities_seq1);
    ind2.assign_cities(cities_seq2);

    cout << "cross over..." << endl;
    Individual new_ind = crossover(ind1, ind2);
}

void cross_test() {
    vector<int> seq1, seq2;
    for (int i = 0; i < 10; i++) {
        seq1.push_back(i);
        seq2.push_back(i);
    }
    shuffle(seq2.begin(), seq2.end(), default_random_engine(10));

    vector<vector<int>> cities1, cities2;
    int index = 0;
    for (int i = 0; i < SalesmanNum; i++) {
        for(int j=0; j<5; j++)
        seq1.push_back(i);
        seq2.push_back(i);
        index += 1;
    }

    int left = 3, right = 6;
    shuffle(seq2.begin(), seq2.end(), default_random_engine(10));
    Individual ind1(SalesmanNum), ind2(SalesmanNum);
}

void test_rand() {
    for (int i = 0; i < 10; i++) {
        cout << getRandomInt(2) << " ";
    }
    
}

int main() {
	read_data();
    read_color();
    test5();
	//print_dist_matrix();
    //print_color_matrix();
    //print_city_color(20);`
    
	return 0;
}