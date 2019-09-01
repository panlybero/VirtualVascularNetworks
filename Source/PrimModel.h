#pragma once
#include"Gene.h"
#include"Point.h"
#include<string>
#include<vector>

using namespace std;

int mainModel(Point* poiints, int n, int side, bool seeded = false , Gene* seed = nullptr);
Point* getExistingPoints(string path, int &v);
//Gene* generateRandomGene(const int size);
void geneToTree(Gene* gene, float** graph, int NofPoints, vector<int>** tree, int n);

void printTree(string path, vector<int>** tree, Point* points, int v, double fit);
Point* getDavidPoints(double side, double nearestRadius, int &n);
