#pragma once
#include<vector>
#include <utility> 
#include<iostream>
#include<random>
#include<queue>
#include<math.h>
#include<time.h>
#include"Gene.h"
#include"PrimModel.h"
#include"FileHandler.h"
#include"util.h"
#include"Network.h"


using namespace std;
class TreeGene
{
public:
	TreeGene(vector<int>** tree,int n);
	TreeGene(int n);
	TreeGene(TreeGene* other);
	~TreeGene();
	void randomizeTree();
	void dump();
	int size();
	void addNode(int i, int v);
	vector<int>* getNode(int i);
	int findintree(int v);
	vector<int>* findintree(vector<int>* v);
	bool isConnected();
	double fitness = 0;
	int generation = 0;
	vector<int>** getTree();
	Gene* getRadGene();
	void mutate();
	void mutateAnyPoint();
	void shuffleSubTree();
	double mutrate = 10;

private:
	vector<int>** m_tree;
	Gene* radGene;
	int m_n;
};

pair<TreeGene*, TreeGene*> makeChildren(TreeGene* p1, TreeGene* p2);


void printTree(string path, TreeGene* tree, Point* points, int v);
double getFitness(TreeGene* g, Point* points,  int NofPoints);
int gettheLeaks();

string stringTree(TreeGene * tree, Point * points, int v);