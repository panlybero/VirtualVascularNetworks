#pragma once
#include"Gene.h"
#include<utility>
#include"Point.h"
#include"util.h"
#include<vector>
#include<iostream>
#include<future>
#include<string>
#include"FileHandler.h"

using namespace std;


class Network
{
public:
	Gene* gene;
	Gene* radGene;
	Point* points;
	Network(int nOfPoints, Point* pnts, bool seeded = false, Gene* seed = nullptr);
	Network(Network* other);
	~Network();
	void mutate(bool sure);
	void restrictedMutate();
	int getGeneLength();
	int getRadGeneLength();
	int getN();
	double getFitness(string mode);
	double getFit();
	double tpl;
	double mpl;
	double ztot;
	double cons;
	double spcfl;
	string optstr;
private:
	int n;
	int geneLength;
	int radGeneLength;
	double fit;

	

};

pair<Network*, Network*> cross(Network* p1, Network* p2);
//Gene* generateRandomGene(int size);
//Gene* generateRandomRadGene(int size);
double* getNetworkProbs(vector<Network*>* pool, string mode);
vector<int>**  geneToTree(Network* network);
void primMST(float** graph, int* parent, int n, bool bin = false);
double totalPathLength(vector<int>** tree, Point* points, int v);
double meanPathLength(vector<int>** tree, Point* points, int n);
double getNetworkImpedence(vector<int>** tree, Gene* radGene, Point* points, int NofPoints);
double getNetworkImpedence(vector<int>** tree, Network* network, int NofPoints, string mode);
void assignRadsToPoints(vector<int>** tree, Point* points, int v, Gene* radGene);
double totalPathLength(vector<int>** tree, Point* points, int v);
double spaceFilling(vector<int>** tree, Network* network, int NofPoints);


vector<Network*>* nextPopulation(vector<Network*>* pool, int generation, int totalgens, string mode);
vector<Network*>* nextSubPopulation(vector<Network*>* subpool, string selectmode, string mode);
vector<vector<Network*>*> dividePopulation(vector<Network*>* pool);
pair<int, int> chooseParentsFitProp(double* probs, int size);
pair<int, int> chooseParentsTournament(vector<Network*>* subpool, int size, string mode);
void printTree(string path, vector<int>** tree, Network* net, Point* points, int v,double fit);
double getConservation(vector<int>** tree, Network* network, int NofPoints, string mode);
double calcConstFlow(vector<int>** tree, Network* net, int curr, int par);
int countNiches(vector<Network*>* pool, int NofPoints);
int getLeaky();
void enforceAreaConservation(vector<int>** tree, Network* network, int NofPoints);

string makeTreeString(vector<int>** tree, Point* points, int v, Gene* radGene);
