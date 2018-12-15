#include"util.h"
#include<random>
#include<iostream>
using namespace std;

double randomDouble(double min, double max)
{
	std::random_device rd;
	std::uniform_real_distribution<double> dist(min, max);
	std::default_random_engine generator(rd());

	double d = 0;//(max - min) * ((double)rand() / (double)RAND_MAX) + min;
	d = dist(generator);
	//cout << d << endl;
	return d;
}

int randomInt(int min, int max)
{
	std::random_device rd;
	std::uniform_int_distribution<int> dist(min, max);
	std::default_random_engine generator(rd());

	double d = 0;//(max - min) * ((double)rand() / (double)RAND_MAX) + min;
	d = dist(generator);
	//cout << d << endl;
	return d;
}
double eucDistance(Point a, Point b)
{
		return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y) + (a.z - b.z)*(a.z - b.z));
}