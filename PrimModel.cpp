// A C / C++ program for Prim's Minimum Spanning Tree (MST) algorithm. 
// The program is for adjacency matrix representation of the graph

#include <stdio.h>
#include<iostream>
#include <limits.h>
#include<random>

#include<ctime>
#include<vector>
#include"FileHandler.h"
#include"Gene.h"
#include<stack>
#include<direct.h>
#include <windows.h> 
#include"PrimModel.h"
#include"Network.h"
#include"util.h"
#include<math.h>
#include<queue>
#include<future>
#include<string>
#include<ctime>

using namespace std;// Number of vertices in the graph
#define V 50
#define SIDE 5
#define GEN 2000
#define TESTTAR 15
#define MUTRATE 15  // 0 means 0% chance,  100 means 100%
#define RUNS 5

#define NEWPOINTS true //deprecated
#define BIMODAL false
#define OPTFORRAD true
#define THREEDEE false


const int population = 175;
int mutations = 0;
Gene* best;
Gene* bestRad;
double bestFit=0;
vector<double> gens;

string mode = "pulsing";

struct Parents
{
public:
	Parents(): first(-1),second(-1)
	{

	}
	int first;
	int second;
};

// A utility function to find the vertex with minimum key value, from
// the set of vertices not yet included in MST


// Function to construct and print MST for a graph represented using adjacency
// matrix representation


void geneToTree(Gene* gene, float** graph, int NofPoints, vector<int>** tree, int n)
{
	/* Let us create the following graph
	2    3
	(0)--(1)--(2)
	|   / \   |
	6| 8/   \5 |7
	| /     \ |
	(3)-------(4)
	9          */


	const int geneLength = NofPoints*(NofPoints - 1) / 2;
	//double gene[geneLength];


	int index = NofPoints;
	//making tree 0-1 0-3, 1-2, 3-4


	float** graphprime = new float*[n];



	for (int i = 0; i < n; i++)
	{
		graphprime[i] = new float[n];
		for (int j = 0; j < n; j++)
		{
			graphprime[i][j] = graph[i][j];
		}

	}
	index = 0;
	for (int i = 0; i < n; i++) // Look at "Representing trees with genetic algorithms" 
	{
		for (int j = i + 1; j < n; j++)
		{

			graphprime[i][j] += 2 * gene->getPos(index);
			graphprime[j][i] += 2 * gene->getPos(index);
			index++;
		}
	}

	// Print the solution
	int* parent = new int[n];
	primMST(graphprime, parent,n, true);
	
	gene->dump();
	system("pause");

	//vector<int>* tree[V];
	for (int i = 0; i < n; i++)
	{
		tree[i] = new vector<int>();
	}

	for (int i = 1; i < n; i++)
	{
		tree[parent[i]]->push_back(i);
	}
	
	for (int i = 0; i < n; i++)
	{

		delete[] graphprime[i];

	}
	delete[] graphprime;
	delete[] parent;
}
/*
Gene* generateRandomGene(const int size)
{
	float* gene = new float[size];
	Gene* g = new Gene(gene, size);
	g->setMutrate(0.15);//randomDouble(0, 1));

	for (int i = 0; i < size; i++)
	{
		g->setPos(randomDouble(0, 2),i);

	}

	return g;
}
Gene* generateRandomRadGene(int size)
{
	float* f = new float[size];
	Gene* g = new Gene(f, size);
	g->setMutrate(0.15);//randomDouble(0, 1));
	if(OPTFORRAD)
	for (int i = 0; i < size; i++)
	{
		g->setPos(randomDouble(1, 5), i);
	}
	else 
	{
		for (int i = 0; i < size; i++)
		{
			g->setPos(3, i);
		}
	}

	return g;
}
*/
double networkModel(int index, int v, float** graph, Point* points, FileHandler* f, string folderPath, bool seeded , Gene* seed );

double getFitness(Gene* gene, Gene* radGene, Point* points, float** graph, int NofPoints);

double* getProbs(vector<Gene*>* pool, vector<Gene*>* radPool, Point* points, float** graph, int NofPoints);

Parents getParents(double* probs);

vector<Gene*> GenNewPop(vector<Gene*> pool, vector<Gene*> &radPool,Parents parents, int geneLength, Point* points, float** graph, int n);

pair<vector<Gene*>*, vector<Gene*>*> GenNewPopSplit(vector<Gene*>* pool, vector<Gene*>* radPool, int geneLength, Point* points, float** graph, int n);

void mutate(Gene* gene, Gene* radGene, int geneLength, int rGeneLength);

void biModMut(Gene* gene, Gene* radGene, int rGeneLength,int genelength);

double testFitness(double* gene, double target);

double* testProbs(vector<double*> pool);

void testmut(double* gene, int n);

float pathLength(vector<int> tree, Point* points);

//string makeTreeString(vector<int>** tree, Point* points, int v);


float pathLength(vector<int> path, Point* points)
{
	float pathL = 0.0f;
	for (int i = 0; i < path.size()-1; i++)
	{
		pathL += eucDistance(points[path.at(i)], points[path.at(i + 1)]);
	}

	return pathL;
}

double meanPathLength(vector<int>** tree, Point* points,int n)
{
	int nOfTips = 0;
	int prevPoint = 0;

	for (int i = 0; i < n; i++)
	{
		if (tree[i]->size() == 0)
		{
			nOfTips++;
		}
	}

	float currPathLength = 0.0f;
	float meanPathLength = 0.0f;
	int curr = 0;
	stack<int> s;
	s.push(0);
	vector<int> v;
	int ind = 0;

	while (!s.empty())
	{
		
		curr = s.top();
		s.pop();
		v.push_back(curr);
		for (int i = 0; i < tree[curr]->size(); i++)
		{
			s.push(tree[curr]->at(i));
			
		}

		if (tree[curr]->size() == 0)
		{

			meanPathLength += pathLength(v, points) / nOfTips;

			v.clear();
		}

	}
	return meanPathLength;
}


double eval(double* gene, int n)
{
	double tot = 0;
	for (int i = 0; i < n; i++)
	{
		tot += gene[i];
	}
	return tot;
}
//Code for this function adapted from Savage and Hunt(2016) with permission from the authors. 
Point* getDavidPoints(double side, double nearestRadius, int &n)
{
	//This function creates a circular dish and defines the service points using the method described in the Phys Rev paper. 

	int maxConsecutiveFailedInsertions = 10000;
	int consecutiveFailedInsertions(0);
	vector<Point> b;
	b.push_back(Point());
	b.back().x = side / 2.0;
	b.back().y = side / 2.0;

	while (consecutiveFailedInsertions < maxConsecutiveFailedInsertions) {
		double x(randomDouble(0, 1)*side), y(randomDouble(0, 1)*side);
		while (eucDistance(Point(side / 2.0, side / 2.0), Point(x, y)) > side / 2.0) {
			x = randomDouble(0, 1)*side;
			y = randomDouble(0, 1)*side;
		}
		bool validLocation(true);
		for (unsigned int i = 0; i < b.size(); i++) {
			if (eucDistance(Point(x, y), b[i]) < nearestRadius) {
				validLocation = false;
				i = b.size();
			}
		}
		if (validLocation) {
			b.push_back(Point());
			b.back().x = x;
			b.back().y = y;
			consecutiveFailedInsertions = 0;
		}
		else
			consecutiveFailedInsertions++;
	}
	Point* points = new Point[b.size()];
	for (int i = 0; i < b.size(); i++)
	{
		points[i] = b.at(i);
	}

	n = b.size();
	return points;
}

Point* getExistingPoints(string path, int &v)
{
	FileHandler fh;
	Point* points;
	if (THREEDEE)
	{
		points = fh.read3DPoints(path, v);
	}
	else
	{
		points = fh.readPoints(path, v);
	}
	 
	return points;
}
int talk = 0;
void printGen(vector<double> a, int n, FileHandler* f)
{
	
	talk++;
	cout << talk << endl;
	for (int i = 0; i <a.size(); i++)
	{
		f->print(to_string(a[i]));
	}
}

void runModel(int index, int v, float** graph, Point* points, FileHandler* f)
{
	int geneLength = v*(v - 1) / 2;
	int rGeneLength = v - 1;
	



	vector<Gene*>* radPool = new vector<Gene*>();;
	vector<Gene*>* pool = new vector<Gene*>();;
	
	best = new Gene(geneLength);
	bestRad = new Gene(rGeneLength);

	for (int i = 0; i < population; i++)
	{
		//cout << radPool.size();
		pool->push_back(generateRandomGene(geneLength));

		radPool->push_back(generateRandomRadGene(rGeneLength));
		
	}
	
	///////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////
	vector<int>** tree = new vector<int>*[v];


	double* probs = getProbs(pool,radPool, points, graph, v);


	Parents parents = getParents(probs);


	//cout << "hai" << endl;
	
	///////////////////////////////////////////////////////
	


	string allTrees = "";
	for (int i = 0; i < GEN; i++)
	{
		/*
		
		double* oldp = probs;
		delete[] oldp;

		probs = getProbs(pool, radPool,points, graph, v);

		parents = getParents(probs);
		
		pool = GenNewPop(pool, radPool, parents, geneLength, points, graph, v);
		
		//vector<Gene*> pool, vector<Gene*> &radPool,  int geneLength, Point* points, float** graph, int n
		*/

		vector<Gene*>* subpool1 = new vector<Gene*>();
		vector<Gene*>* subpool2 = new vector<Gene*>();
		vector<Gene*>* subpool3 = new vector<Gene*>();
		vector<Gene*>* subpool4 = new vector<Gene*>();


		vector<Gene*>* subradpool1 = new vector<Gene*>();
		vector<Gene*>* subradpool2 = new vector<Gene*>();
		vector<Gene*>* subradpool3 = new vector<Gene*>();
		vector<Gene*>* subradpool4 = new vector<Gene*>();

		//cout << "pool size at i " << i <<pool->size() << endl;
		//cout << "radpool size " << radPool->size() << endl;

		//radPool->at(7)->dump();
		//system("pause");
		for (int j = 0; j < pool->size(); j += 4)
		{
			//cout << j << endl;
			subpool1->push_back(pool->at(j));
			subradpool1->push_back(radPool->at(j));

			subpool2->push_back(pool->at(j+1));
			subradpool2->push_back(radPool->at(j + 1));

			subpool3->push_back(pool->at(j+2));
			subradpool3->push_back(radPool->at(j+2));

			subpool4->push_back(pool->at(j + 3));
			subradpool4->push_back(radPool->at(j + 3));


		}
		//cout << "after" << endl;
		//system("pause");




		future<pair<vector<Gene*>*, vector<Gene*>*>> t1 = async(GenNewPopSplit,subpool1,subradpool1,geneLength,points,graph,v);
		future<pair<vector<Gene*>*, vector<Gene*>*>> t2 = async(GenNewPopSplit, subpool2, subradpool2, geneLength, points, graph, v);
	
		future<pair<vector<Gene*>*, vector<Gene*>*>> t3 = async(GenNewPopSplit, subpool3, subradpool3, geneLength, points, graph, v);
		future<pair<vector<Gene*>*, vector<Gene*>*>> t4 = async(GenNewPopSplit, subpool4, subradpool4, geneLength, points, graph, v);


		pair<vector<Gene*>*, vector<Gene*>*> newpools1 = t1.get();
		pair<vector<Gene*>*, vector<Gene*>*> newpools2 = t2.get();

		pair<vector<Gene*>*, vector<Gene*>*> newpools3 = t3.get();
		pair<vector<Gene*>*, vector<Gene*>*> newpools4 = t4.get();


		subpool1 = newpools1.first;
		subradpool1 = newpools1.second;
		


		subpool2 = newpools2.first;
		subradpool2 = newpools2.second;

		subpool3 = newpools3.first;
		subradpool3 = newpools3.second;

		subpool4 = newpools4.first;
		subradpool4 = newpools4.second;


		pool->clear();
		radPool->clear();

		pool->insert(pool->end(), subpool1->begin(), subpool1->end());
		pool->insert(pool->end(), subpool2->begin(), subpool2->end());
		pool->insert(pool->end(), subpool3->begin(), subpool3->end());
		pool->insert(pool->end(), subpool4->begin(), subpool4->end());


		radPool->insert(radPool->end(), subradpool1->begin(), subradpool1->end());
		radPool->insert(radPool->end(), subradpool2->begin(), subradpool2->end());
		radPool->insert(radPool->end(), subradpool3->begin(), subradpool3->end());
		radPool->insert(radPool->end(), subradpool4->begin(), subradpool4->end());


		if (i > 0 && i % 100 == 0)
		{
			//system("CLS");
			cout << i <<" / "<< GEN << endl;

		}
		delete subpool1;
		delete subpool2;
		delete subradpool1;
		delete subradpool2;

	}
	
	double maxFitness = 0;
	double avgFitness = 0;
	int maxi = 0;
	double curr = 0;
	for (int i = 0; i < population; i++)
	{
		curr = getFitness(pool->at(i), radPool->at(i),points, graph, v);
		avgFitness += curr / population;
		if (curr>maxFitness)
		{
			maxFitness = curr;
			maxi = i;
		}
	}
	gens.push_back (maxFitness);
	//vector<int>** tree = new vector<int>*[v];
	
	//getFitness(pool.at(maxi), radPool.at(maxi), points, graph, v);
	radPool->at(maxi)->dump();
	geneToTree(pool->at(maxi),graph,v,tree,v);
	assignRadsToPoints(tree, points, v, radPool->at(maxi));
	//printTree("D:\\Generations\\Test\\final.txt", tree, points, v);


////////////////////////////////////////////////////////////////////////////////////
	
	printGen(gens, 0, f);
////////////////////////////////////////////////////////////////////////////////////
	for (int i = 0; i < population; i++)
	{
		delete pool->at(i);
		delete radPool->at(i);
	}
	FileHandler allT("D:\\Generations\\Test\\allT.txt");
	allT.print(allTrees);
	cout << gens.at(0) << "  " << gens.back() << endl;
}

int mainModel(Point* poiints, int n, int side, bool seeded, Gene* seed)
{
	int start_s = clock();
	
	//double graph[V][V];
	//Point points[V];

	srand(time(NULL));
	/////////////////////////////////////////////

	double nearestRadius = 1.0 / sqrt(acos(-1.0));
	int v = 8;
	Point* points;
	
	if (NEWPOINTS)
	{
		points = poiints;//getDavidPoints(side, nearestRadius, v);
		v = n;
	}
	else
	{

		points = getExistingPoints("best.txt", v);

	}

	

	float** graph = new float*[v];
	for (int i = 0; i < v; i++)
	{
		graph[i] = new float[v];
	}
	//Setup Graph

	for (int i = 0; i < v; i++)
	{
		for (int j = 0; j < v; j++)
		{
			graph[i][j] = 1.0f;
			if (i == j) graph[i][j] = 0.0f;
		}
	}
	//Setup Points

	/*
	Point* points = new Point[v];
	
	for (int i = 0; i < v; i++)
	{
		points[i].x = randomDouble(0, 10);
		points[i].y = randomDouble(0, 10);
	}
	*/

	

	

	//
	
	double* bests = new double[RUNS];

	for (int i = 0; i < RUNS; i++)
	{
		//runModel(0, v, graph, points, f);
		string foldername = to_string(side) + "_" + to_string(GEN) + "_" + mode + "_" + to_string(clock());
		string newFolder = "D:\\Generations\\Project22\\" + foldername;
		string newFile = newFolder + "\\file" + to_string(0) + ".txt";
		FileHandler* f = new FileHandler(newFile);
		bests[i] = networkModel(i, v, graph, points, f, newFolder,seeded, seed);//runModel(i, v, graph, points, f);
		points =  getDavidPoints(side, nearestRadius, v);
	}
	
	for (int i = 0; i < RUNS; i++)
	{
		cout << bests[i] << endl;
	}
	
	////////////////////////////////////////////
	

	/*
	geneToTree(pool.at(1), graph, V, tree);
	printTree("D:\\file4.txt", tree, points);
	cout << "4 " << totalPathLength(tree, points) << endl;
	
	*/


	//cout << "leaks = " << getLeaky() << endl;
	

	int stop_s = clock();
	cout << "time: " << ((stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000)/1000 <<" sec" << endl;

	cout << mutations << endl;
	system("pause");
	
	
	
	
	/*
	int geneLength = 15;
	vector<double*> pool;
	for (int i = 0; i < 10; i++)
	{
		pool.push_back(new double[15]);
		for (int j = 0; j < 15; j++)
		{
			pool.at(i)[j] = randomDouble(0, 10);
		}
	}
	double* probs;
	Parents parents;
	bool stop = false;
	for (int i = 0; !stop; i++)
	{
		probs = testProbs(pool);
		parents = getParents(probs);
		//cout << eval(pool.at(parents.first), geneLength)<<endl;
		if (eval(pool.at(parents.first), geneLength) <= 16 && eval(pool.at(parents.first), geneLength) >= 14)
		{
			stop = true;
		}
		pool = crossGenes(pool, parents, geneLength);
	}

	cout << eval(pool.at(parents.first), geneLength) << endl;
	*/
	return 0;
}

vector<Network*>* radiationEvent(vector<Network*>* pool)
{
	for (int i = 0; i < population; i++)
	{
		pool->at(i)->gene->setMutrate(1);
		pool->at(i)->radGene->setMutrate(1);
	}
	return pool;
}

vector<Network*>* stopRadiationEvent(vector<Network*>* pool)
{
	for (int i = 0; i < population; i++)
	{
		double genemut =(pool->at(i)->getGeneLength());
		genemut *= 0.001;
		double radgenemut =  (pool->at(i)->getRadGeneLength());
		radgenemut *= 0.001;
		pool->at(i)->gene->setMutrate(genemut);
		pool->at(i)->radGene->setMutrate(radgenemut);
	}
	return pool;
}

double networkModel(int index, int v, float** graph, Point* points, FileHandler* f, string folderPath, bool seeded = false, Gene* seed = nullptr)
{

	int geneLength = v*(v - 1) / 2;
	int rGeneLength = v - 1;
	int radiation_counter = 0;
	bool radiation_flag = false;
	double final_fit = 0;
	vector<Network*>* pool = new vector<Network*>();

	//best = new Network(v);
	if (seeded) 
	{
		for (int i = 0; i < population; i++)
		{
			//cout << radPool.size();
			pool->push_back(new Network(v, points, seeded, seed));
			if(i>0)
				pool->at(i)->mutate(false);


		}
	}
	else
	{
		for (int i = 0; i < population; i++)
		{
			//cout << radPool.size();
			pool->push_back(new Network(v, points));


		}
	}
	string allTrees = "";
	for (int i = 0; i < GEN; i++)
	{

		

		if (radiation_flag)
		{
			radiation_counter++;
		}
		if ((i+1) % 2000 == 20001)
		{
			cout << "Radiation Event: " << endl;
			radiation_flag = true;
			pool = radiationEvent(pool);
		}
		if (radiation_counter == 1000)
		{
			cout << "End of Radiation Event: " << endl;
			radiation_flag = false;
			radiation_counter = 0;
			pool = stopRadiationEvent(pool);
		}

		//if((i+1) % 1000 == 0) cout << "Niches: " << countNiches(pool, pool->at(0)->getN()) << endl;


		//if (i == GEN / 2) cout << "Switching to Fitness Propotional selection" << endl;
		pool = nextPopulation(pool, i, GEN, mode);
		double max = 0;
		int max_j = 0;
		

		if (i % 100 == 0)
		{
			double max = 0;
			double mut = 0;
			double mutrad = 0;
			int keep = 0;
			double avg_fit = 0;
			for (int j = 0; j < population; j++)
			{
				avg_fit += pool->at(j)->getFit() / population;

				if (pool->at(j)->getFit() > max)
				{
					max = pool->at(j)->getFit();
					keep = j;
					mut = pool->at(j)->gene->getMutrate();
					mutrad = pool->at(j)->radGene->getMutrate();
				}
			}
			gens.push_back(max);
			//cout<<index<<"/"<<RUNS<<" "<<i<<" / "<<GEN<< " bestfit " << max*100000.0 << " AvgFit "<<avg_fit* 100000.0 <<" Gene mutrate"<<mut<< " radGene mutrate " << mutrad <<endl;
			cout << index << "/" << RUNS << " " << i << " / " << GEN << " bestcost " << 1 / max << " AvgcCost " << 1 / avg_fit << " Gene mutrate" << mut << " radGene mutrate " << mutrad << endl;
			cout << pool->at(keep)->tpl << " + " << pool->at(keep)->mpl << " + " << pool->at(keep)->ztot << " + " << pool->at(keep)->cons << " + "<<pool->at(keep)->spcfl<< endl;

			//pool->at(keep)->radGene->dump();

			//cout << "Leaks " << getLeaky() << endl;
		}

		if (i == 1 || i == int(GEN/3) || i == int(2*GEN / 3) || i == GEN-1)
		{
			for (int j = 0; j < population; j++)
			{


				if (pool->at(j)->getFit() > max)
				{
					max = pool->at(j)->getFit();
					max_j = j;

				}
			}
			Network* best = pool->at(max_j);
			allTrees += makeTreeString(geneToTree(best), best->points, best->getN(),best->radGene);
			
		}
		/*
		if (i < GEN - 1)
		{
			cout << "del " << i << endl;
			system("pause");
			while (!oldPool->empty())
			{
				Network* tmp = oldPool->back();
				oldPool->pop_back();
				delete tmp;

			}
			cout << "deld " << i << endl;
			system("pause");
		}
		*/
	
	}
	pool->at(0)->getFitness(mode);
	string opt = pool->at(0)->optstr;
	int max_i = 0;
	double max_fit = 0;
	string newFolder = folderPath + "_" + opt;
	CreateDirectory(newFolder.c_str(), NULL);
	for (int i = 0; i < pool->size(); i++)
	{
		
		if (pool->at(i)->getFit() > max_fit)
		{
			max_fit = pool->at(i)->getFit();
			max_i = i;
		}
		
		Network* best = pool->at(i);
		
		vector<int>** tree = geneToTree(best);
		double tmpfit = best->getFitness(mode);
		
		//best->gene->dump();
		//best->radGene->dump();
		printTree(folderPath+"_"+ opt+ "\\best" + to_string(i) + ".txt", tree, best, best->points, best->getN(), tmpfit);

	}

	Network* best = pool->at(max_i);

	vector<int>** tree = geneToTree(best);
	best->getFitness(mode);
	opt = best->optstr;
	final_fit = best->getFitness(mode);
	printTree(newFolder+ "\\best" + ".txt", tree, best, best->points, best->getN(), final_fit);

	FileHandler allT("D:\\Generations\\Allts\\3ts.txt");
	allT.print(allTrees);
	//Network* best = pool->at(max_i);
	/*
	vector<int>** tree = geneToTree(best);
	best->getFitness(mode);
	best->gene->dump();
	//best->radGene->dump();
	printTree(folderPath+"\\best"+to_string(index)+".txt", tree, best, best->points, best->getN());
	//printGen(gens, 0, new FileHandler("D:\\Generations\\UpdatedModel\\gens.txt"));
	//system("pause");
	*/
	cout <<"size " <<pool->size() << endl;
	cout << "cleanup" << endl;
	while (!pool->empty())
	{
		Network* tmp = pool->back();
		pool->pop_back();
		delete tmp;
	}



	/*
	///////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////
	vector<int>** tree = new vector<int>*[v];


	double* probs = getNetworkProbs(pool);


	Parents parents = getParents(probs);


	//cout << "hai" << endl;

	///////////////////////////////////////////////////////



	string allTrees = "";
	for (int i = 0; i < GEN; i++)
	{
		/*

		double* oldp = probs;
		delete[] oldp;

		probs = getProbs(pool, radPool,points, graph, v);

		parents = getParents(probs);

		pool = GenNewPop(pool, radPool, parents, geneLength, points, graph, v);

		//vector<Gene*> pool, vector<Gene*> &radPool,  int geneLength, Point* points, float** graph, int n
		brek here

		vector<Gene*>* subpool1 = new vector<Gene*>();
		vector<Gene*>* subpool2 = new vector<Gene*>();
		vector<Gene*>* subpool3 = new vector<Gene*>();
		vector<Gene*>* subpool4 = new vector<Gene*>();


		vector<Gene*>* subradpool1 = new vector<Gene*>();
		vector<Gene*>* subradpool2 = new vector<Gene*>();
		vector<Gene*>* subradpool3 = new vector<Gene*>();
		vector<Gene*>* subradpool4 = new vector<Gene*>();

		//cout << "pool size at i " << i <<pool->size() << endl;
		//cout << "radpool size " << radPool->size() << endl;

		//radPool->at(7)->dump();
		//system("pause");
		for (int j = 0; j < pool->size(); j += 4)
		{
			//cout << j << endl;
			subpool1->push_back(pool->at(j));
			subradpool1->push_back(radPool->at(j));

			subpool2->push_back(pool->at(j + 1));
			subradpool2->push_back(radPool->at(j + 1));

			subpool3->push_back(pool->at(j + 2));
			subradpool3->push_back(radPool->at(j + 2));

			subpool4->push_back(pool->at(j + 3));
			subradpool4->push_back(radPool->at(j + 3));


		}
		//cout << "after" << endl;
		//system("pause");




		future<pair<vector<Gene*>*, vector<Gene*>*>> t1 = async(GenNewPopSplit, subpool1, subradpool1, geneLength, points, graph, v);
		future<pair<vector<Gene*>*, vector<Gene*>*>> t2 = async(GenNewPopSplit, subpool2, subradpool2, geneLength, points, graph, v);

		future<pair<vector<Gene*>*, vector<Gene*>*>> t3 = async(GenNewPopSplit, subpool3, subradpool3, geneLength, points, graph, v);
		future<pair<vector<Gene*>*, vector<Gene*>*>> t4 = async(GenNewPopSplit, subpool4, subradpool4, geneLength, points, graph, v);


		pair<vector<Gene*>*, vector<Gene*>*> newpools1 = t1.get();
		pair<vector<Gene*>*, vector<Gene*>*> newpools2 = t2.get();

		pair<vector<Gene*>*, vector<Gene*>*> newpools3 = t3.get();
		pair<vector<Gene*>*, vector<Gene*>*> newpools4 = t4.get();


		subpool1 = newpools1.first;
		subradpool1 = newpools1.second;



		subpool2 = newpools2.first;
		subradpool2 = newpools2.second;

		subpool3 = newpools3.first;
		subradpool3 = newpools3.second;

		subpool4 = newpools4.first;
		subradpool4 = newpools4.second;


		pool->clear();
		radPool->clear();

		pool->insert(pool->end(), subpool1->begin(), subpool1->end());
		pool->insert(pool->end(), subpool2->begin(), subpool2->end());
		pool->insert(pool->end(), subpool3->begin(), subpool3->end());
		pool->insert(pool->end(), subpool4->begin(), subpool4->end());


		radPool->insert(radPool->end(), subradpool1->begin(), subradpool1->end());
		radPool->insert(radPool->end(), subradpool2->begin(), subradpool2->end());
		radPool->insert(radPool->end(), subradpool3->begin(), subradpool3->end());
		radPool->insert(radPool->end(), subradpool4->begin(), subradpool4->end());


		if (i > 0 && i % 100 == 0)
		{
			//system("CLS");
			cout << i << " / " << GEN << endl;

		}
		delete subpool1;
		delete subpool2;
		delete subradpool1;
		delete subradpool2;

	}

	double maxFitness = 0;
	double avgFitness = 0;
	int maxi = 0;
	double curr = 0;
	for (int i = 0; i < population; i++)
	{
		curr = getFitness(pool->at(i), radPool->at(i), points, graph, v);
		avgFitness += curr / population;
		if (curr>maxFitness)
		{
			maxFitness = curr;
			maxi = i;
		}
	}
	gens.push_back(maxFitness);
	//vector<int>** tree = new vector<int>*[v];

	//getFitness(pool.at(maxi), radPool.at(maxi), points, graph, v);
	radPool->at(maxi)->dump();
	geneToTree(pool->at(maxi), graph, v, tree, v);
	assignRadsToPoints(tree, points, v, radPool->at(maxi));
	printTree("D:\\Generations\\Test\\final.txt", tree, points, v);


	////////////////////////////////////////////////////////////////////////////////////

	printGen(gens, 0, f);
	////////////////////////////////////////////////////////////////////////////////////
	for (int i = 0; i < population; i++)
	{
		delete pool->at(i);
		delete radPool->at(i);
	}
	FileHandler allT("D:\\Generations\\Test\\allT.txt");
	allT.print(allTrees);
	cout << gens.at(0) << "  " << gens.back() << endl;
	*/

return final_fit;
}

double* getProbs(vector<Gene*>* pool, vector<Gene*>* radPool,Point* points, float** graph, int NofPoints)
{
	double* fitPool = new double[pool->size()];
	double* probPool = new double[pool->size()];
	double totalFit = 0;
	double maxfit = 0;
	double maxfiti = 0;

	for (int i = 0; i < pool->size(); i++)
	{
		fitPool[i]=getFitness(pool->at(i), radPool->at(i),points, graph, NofPoints);
		totalFit += fitPool[i];
	}

	for (int i = 0; i < pool->size(); i++)
	{
		probPool[i] = fitPool[i] / totalFit;
		if (maxfit < fitPool[i])
		{
			maxfit = fitPool[i];
			maxfiti = i;
		}

	}
	gens.push_back(maxfit);
	delete[] fitPool;
	return probPool;
}

Parents getParents(double* probs)
{
	double value = randomDouble(0, 1);
	Parents parents;

	for (int i = 0; i<population; i++)
	{
		value -= probs[i];
		//cout << value << endl;
		if (value <= 0)
		{
			parents.first = i;
			break;
		}
	}
	value = randomDouble(0, 1);
	for (int i = 0; i<population; i++)
	{
		value -= probs[i];
		if (value <= 0)
		{
			parents.second = i;
			break;
		}
	}
	if (parents.first == -1) parents.first = population - 1;
	if (parents.second == -1) parents.second = population - 1;
	return parents;
}

pair<vector<Gene*>*, vector<Gene*>*> GenNewPopSplit(vector<Gene*>* pool, vector<Gene*>* radPool,  int geneLength, Point* points, float** graph, int n)
{

	int subpop = 25;

	vector<Gene*>* newpool = new vector<Gene*>();;
	vector<Gene*>* newradpool = new vector<Gene*>();;
	Point* temppoints = new Point[n];
	for (int i = 0; i < n; i++)
	{
		temppoints[i] = points[i];
	}


	//vector<double*> Oldpool = pool;

	//cout << parents.first << " " << parents.second << endl;
	double* probs = getProbs(pool, radPool, temppoints, graph, n);

	Parents parents = getParents(probs);
	Gene* p1 = pool->at(parents.first);
	Gene* p2 = pool->at(parents.second);


	Gene* p1rad = radPool->at(parents.first);
	Gene* p2rad = radPool->at(parents.second);

	//Elitism
	double maxFitness = 0;
	int maxi = 0;
	double curr = 0;

	for (int i = 0; i < subpop; i++)
	{
		curr = probs[i];

		//cout << 3 << endl;
		//system("pause");
		if (curr > maxFitness)
		{
			maxFitness = curr;
			maxi = i;
		}
	}
	//keep best genes



	//copyGene(best, pool.at(maxi));
	//copyGene(bestRad, radPool.at(maxi));


	Gene* a = new Gene(new float[geneLength], geneLength);
	//Gene* b = new Gene(new float[geneLength], geneLength);
	Gene* rada = new Gene(new float[n - 1], n - 1);
	//Gene* radb = new Gene(new float[n - 1], n - 1);

	copyGene(a, pool->at(maxi));
	//copyGene(b, pool->at(maxi));
	copyGene(rada, radPool->at(maxi));
	//copyGene(radb, radPool->at(maxi));

	newpool->push_back(a);
	//newpool->push_back(b);
	newradpool->push_back(rada);
	//newradpool->push_back(radb);


	// ONE POINT CROSS
	/*
	for (int i = 0; i < (population-1)/2; i++)
	{
	float* ch1 = new float[geneLength];
	float* ch2 = new float[geneLength];
	int point = rand() % geneLength;
	//cout << "hello "<< point << endl;
	for (int j = 0; j < geneLength; j++)
	{
	if (j<= point)
	{
	ch1[j] = p1->getPos(j);
	ch2[j] = p2->getPos(j);
	}
	else
	{
	ch1[j] = p2->getPos(j);
	ch2[j] = p1->getPos(j);
	}
	}


	Gene* chG1 = new Gene(ch1, geneLength);
	mutate(chG1, geneLength);
	Gene* chG2 = new Gene(ch2, geneLength);
	mutate(chG2, geneLength);


	newpool.push_back(chG1);
	newpool.push_back(chG2);
	//delete[] ch1;
	//delete[] ch2;
	}
	*/

	// TWO POINT CROSS

	for (int i = 0; i < (subpop-1)/2; i++)
	{
		float* ch1 = new float[geneLength];
		float* ch2 = new float[geneLength];
		float* ch1rad = new float[n-1];
		float* ch2rad = new float[n-1];
		int point1 = rand() % geneLength;
		int point2 = (rand() % (geneLength - point1)) + point1;
		int radPoint = rand() % (n - 1);

		//cout << "hello "<< point << endl;
		for (int j = 0; j < geneLength; j++)
		{
			if (j <= point1 || j >= point2)
			{
				ch1[j] = p1->getPos(j);
				ch2[j] = p2->getPos(j);
			}
			else
			{
				ch1[j] = p2->getPos(j);
				ch2[j] = p1->getPos(j);
			}
		}
		//One Point Cross for Rad
		for (int j = 0; j < n - 1; j++)
		{
			if (j <= radPoint)
			{
				ch1rad[j] = p1rad->getPos(j);
				ch2rad[j] = p2rad->getPos(j);

			}
			else
			{
				ch1rad[j] = p2rad->getPos(j);
				ch2rad[j] = p1rad->getPos(j);
			}
		}


		Gene* chG1 = new Gene(ch1, geneLength);
		Gene* chG2 = new Gene(ch2, geneLength);
		Gene* chG1rad = new Gene(ch1rad, n - 1);
		Gene* chG2rad = new Gene(ch2rad, n - 1);
		chG1->setMutrate(p1->getMutrate());
		chG2->setMutrate(p2->getMutrate());


		chG1rad->setMutrate(p1rad->getMutrate());
		chG2rad->setMutrate(p2rad->getMutrate());


		if (BIMODAL)
		{
			//	cout<< chG1->getMutrate() << endl;
			biModMut(chG1, chG1rad, n - 1, geneLength);
			biModMut(chG2, chG2rad, n - 1, geneLength);
		}
		else
		{
			mutate(chG1, chG1rad, geneLength, n - 1);
			mutate(chG2, chG2rad, geneLength, n - 1);

		}

		//


		newpool->push_back(chG1);
		newpool->push_back(chG2);

		newradpool->push_back(chG1rad);
		newradpool->push_back(chG2rad);

		
		//delete[] ch1;
		//delete[] ch2;

	}
	//int o = rand() % geneLength;

	while (!pool->empty())
	{
		Gene* temp = pool->back();
		Gene* radtemp = radPool->back();
		radPool->pop_back();
		pool->pop_back();
		delete temp;
		delete radtemp;
	}

	//newradpool.at(0)->dump();
	//delete radPool;
	pair<vector<Gene*>*, vector<Gene*>*> newpair = make_pair(newpool, newradpool);
	//radPool = newradpool;
	//cout << "radpool size = " << newradpool->size() << endl;
	//system("pause");
	delete[] temppoints;
	return newpair;
}

vector<Gene*> GenNewPop(vector<Gene*> pool, vector<Gene*> &radPool,Parents parents, int geneLength, Point* points, float** graph, int n)
{
	vector<Gene*> newpool;
	vector<Gene*> newradpool;
	Point* temppoints = new Point[n];
	for (int i = 0; i < n; i++)
	{
		temppoints[i] = points[i];
	}


	//vector<double*> Oldpool = pool;

	//cout << parents.first << " " << parents.second << endl;

	Gene* p1 = pool.at(parents.first);
	Gene* p2 = pool.at(parents.second);


	Gene* p1rad = radPool.at(parents.first);
	Gene* p2rad = radPool.at(parents.second);

	//Elitism
	double maxFitness = 0;
	int maxi = 0;
	double curr = 0;

	for (int i = 0; i < population; i++)
	{
		curr = getFitness(pool.at(i), radPool.at(i), points, graph, n);
		
		//cout << 3 << endl;
		//system("pause");
		if (curr > maxFitness)
		{
			maxFitness = curr;
			maxi = i;
		}
	}
	//keep best genes

	gens.push_back(maxFitness);
	

	copyGene(best, pool.at(maxi));
	copyGene(bestRad, radPool.at(maxi));


	Gene* a = new Gene(new float[geneLength], geneLength);
	Gene* b = new Gene(new float[geneLength], geneLength);
	Gene* rada = new Gene(new float[n-1], n-1);
	Gene* radb = new Gene(new float[n-1], n-1);

	copyGene(a, pool.at(maxi));
	copyGene(b, pool.at(maxi));
	copyGene(rada, radPool.at(maxi));
	copyGene(radb, radPool.at(maxi));

	newpool.push_back(a);
	newpool.push_back(b);
	newradpool.push_back(rada);
	newradpool.push_back(radb);


	// ONE POINT CROSS
	/*
	for (int i = 0; i < (population-1)/2; i++)
	{
		float* ch1 = new float[geneLength];
		float* ch2 = new float[geneLength];
		int point = rand() % geneLength;
		//cout << "hello "<< point << endl;
		for (int j = 0; j < geneLength; j++)
		{
			if (j<= point)
			{
				ch1[j] = p1->getPos(j);
				ch2[j] = p2->getPos(j);
			}
			else
			{
				ch1[j] = p2->getPos(j);
				ch2[j] = p1->getPos(j);
			}
		}


		Gene* chG1 = new Gene(ch1, geneLength);
		mutate(chG1, geneLength);
		Gene* chG2 = new Gene(ch2, geneLength);
		mutate(chG2, geneLength);


		newpool.push_back(chG1);
		newpool.push_back(chG2);
		//delete[] ch1;
		//delete[] ch2;
	}
		*/

	// TWO POINT CROSS
		
	for (int i = 0; i < (population - 1) / 2; i++)
	{
		float* ch1 = new float[geneLength];
		float* ch2 = new float[geneLength];
		float* ch1rad = new float[geneLength];
		float* ch2rad = new float[geneLength];
		int point1 = rand() % geneLength;
		int point2 = (rand() % (geneLength-point1))+point1;
		int radPoint = rand() % (n - 1);

		//cout << "hello "<< point << endl;
		for (int j = 0; j < geneLength; j++)
		{
			if (j <= point1 || j>=point2)
			{
				ch1[j] = p1->getPos(j);
				ch2[j] = p2->getPos(j);
			}	
			else
			{
				ch1[j] = p2->getPos(j);
				ch2[j] = p1->getPos(j);
			}
		}
		//One Point Cross for Rad
		for (int j = 0; j < n - 1; j++)
		{
			if (j <= radPoint)
			{
				ch1rad[j] = p1rad->getPos(j);
				ch2rad[j] = p2rad->getPos(j);
				
			}
			else
			{
				ch1rad[j] = p2rad->getPos(j);
				ch2rad[j] = p1rad->getPos(j);
			}
		}


		Gene* chG1 = new Gene(ch1, geneLength);
		Gene* chG2 = new Gene(ch2, geneLength);
		Gene* chG1rad = new Gene(ch1rad, n-1);
		Gene* chG2rad = new Gene(ch2rad, n-1);
		chG1->setMutrate(p1->getMutrate());
		chG2->setMutrate(p2->getMutrate());


		chG1rad->setMutrate(p1rad->getMutrate());
		chG2rad->setMutrate(p2rad->getMutrate());


		if (BIMODAL)
		{
		//	cout<< chG1->getMutrate() << endl;
			biModMut(chG1, chG1rad, n-1, geneLength);
			biModMut(chG2, chG2rad, n-1, geneLength);
		}
		else
		{
			mutate(chG1, chG1rad, geneLength, n-1);
			mutate(chG2, chG2rad, geneLength, n - 1);
			
		}
		
		//
		
		
		newpool.push_back(chG1);
		newpool.push_back(chG2);

		newradpool.push_back(chG1rad);
		newradpool.push_back(chG2rad);
		//delete[] ch1;
		//delete[] ch2;
	
	}
	int o = rand() % geneLength;

	while (!pool.empty())
	{
		Gene* temp = pool.back();
		Gene* radtemp = radPool.back();
		radPool.pop_back();
		pool.pop_back();
		delete temp;
		delete radtemp;
	}

	//newradpool.at(0)->dump();

	radPool = newradpool;
	
	
	return newpool;
}

void biModMut(Gene* gene, Gene* radGene, int rGeneLength, int genelength)
{
	int a = rand() % 100;
	int doit = randomDouble(0,1);
	float sigma = 1;
	double max = 0;

	for (int i = 0; i < genelength; i++)
	{
		if (max < gene->getPos(i)) max = gene->getPos(i);
	}
	
	normal_distribution<float> d0(0,sigma);
	normal_distribution<float> dmax(max, sigma);
	default_random_engine generator;
	//cout << radGene->getMutrate() << endl;
	int pos = rand() % genelength;
	if(doit< gene->getMutrate())
	if (a <= 50) // establish connection
	{
		gene->setPos(d0(generator), pos);
		//cout << gene->getPos(pos) << endl;
	}
	else // cut connection
	{
		gene->setPos(dmax(generator), pos);
	//	cout << gene->getPos(pos) << endl;
	}
	doit = randomDouble(0, 1);
	if (doit <  radGene->getMutrate())
	{
		pos = rand() % rGeneLength;
		radGene->setPos(randomDouble(1, 5), pos);
	}
	/*
	doit =  randomDouble(0, 1);
	if (doit < gene->getMutrate())
	{
		//gene->setMutrate(randomDouble(0, 1));
	}
	doit = randomDouble(0, 1);
	if (doit < radGene->getMutrate())
	{
		//radGene->setMutrate(randomDouble(0, 1));
	}
	*/

}

void mutate(Gene* gene, Gene* radGene, int geneLength, int rGeneLength)
{
	double rate = (double)MUTRATE/100;
	double ev = randomDouble(0, 1);
	double evr = randomDouble(0, 1);
	double pos = rand() % geneLength;
	double posr = rand() % rGeneLength;

	if (ev < gene->getMutrate())
	{
		gene->setPos(randomDouble(0, 2), pos);
	}

	if (evr < radGene->getMutrate() && OPTFORRAD)
	{
		radGene->setPos(randomDouble(1, 5), posr);
	}


}

double testFitness(double * gene, double target)
{
	double total = 0;
	for (int i = 0; i < 10; i++)
	{
		total += gene[i];
	}
	return abs(1/(total-target));
}

double* testProbs(vector<double*> pool)
{
	double totFit = 0;
	double* probs = new double[pool.size()];
	double* fitpool = new double[pool.size()];
	for (int i = 0; i < pool.size(); i++)
	{
		fitpool[i] = testFitness(pool[i], TESTTAR);
		totFit += testFitness(pool[i], TESTTAR);
	}
	for (int i = 0; i < pool.size(); i++)
	{
		probs[i] = fitpool[i] / totFit;
	}
	delete[] fitpool;
	return probs;
}

void testmut(double * gene, int n)
{
	for (int i = 0; i < n; i++)
	{
		//int a = rand() % 10;
		int b = rand() % 100;
		if (b <= n)
		{
			gene[b] = randomDouble(0,2);
		}

	}

}

void assignRadsToPoints(vector<int>** tree, Point* points, int v, Gene* radGene)
{
	points[0].rad = 100;

	for (int i = 1; i < v; i++)
	{
		for (int j = 0; j < tree[i]->size(); j++)
		{
			points[tree[i]->at(j)].rad = radGene->getPos(tree[i]->at(j) - 1);
		}
		if (tree[i]->size() == 1)
		{
			points[tree[i]->at(0)].rad = points[i].rad;
		}
	}
	/*

	for (int i = 0; i < v; i++)
	{
		for (int j = 0; j < tree[i]->size(); j++)
		{
			points[tree[i]->at(j)].rad = radGene->getPos(tree[i]->at(j) - 1);
		}
		
	}
	for (int i = 0; i < v; i++)
	{
		
		if (tree[i]->size() == 1)
		{
			points[tree[i]->at(0)].rad = points[i].rad;
		}
		
		
	}
	*/
	
}

string makeTreeString_dep(vector<int>** tree, Point* points, int v)
{
	string res = "";
	res+=to_string(v)+"\r\n";
	for (int i = 0; i < v; i++)
	{
		res+=to_string(i) + " " + to_string(points[i].x) + " " + to_string(points[i].y) + " " + to_string(points[i].rad)+"\r\n";
	}
	int connections = 0;
	for (int i = 0; i < v; i++)
	{
		for (int j = 0; j < tree[i]->size(); j++)
		{
			connections++;
		}
	}
	res+=(to_string(connections))+"\r\n";
	for (int i = 0; i < v; i++)
	{
		for (int j = 0; j < tree[i]->size(); j++)
		{
			res+=to_string(i) + "-" + to_string(tree[i]->at(j))+"\r\n";
		}
	}
	return res;
}

void printTree(string path,vector<int>** tree, Point* points, int v, double fit)
{
	FileHandler f(path);

	f.print(to_string(v));
	
	if (!THREEDEE)
		for (int i = 0; i < v; i++)
		{
			f.print(to_string(i) + " " + to_string(points[i].x) + " " + to_string(points[i].y) + " " + to_string(points[i].rad));
		}
	else {
		for (int i = 0; i < v; i++)
		{
			f.print(to_string(i) + " " + to_string(points[i].x) + " " + to_string(points[i].y) + " " + to_string(points[i].z) + " " + to_string(points[i].rad));
		}
	}
	int connections = 0;
	for (int i = 0; i < v; i++)
	{
		for (int j = 0; j < tree[i]->size(); j++)
		{
			connections++;
		}
	}
	f.print(to_string(connections));
	for (int i = 0; i < v; i++)
	{
		for (int j = 0; j < tree[i]->size(); j++)
		{
			f.print(to_string(i) + "-" + to_string(tree[i]->at(j)));
		}
	}
	f.print(to_string(fit));
	///////////////////////////////////////////////////////////


}

double getFitness(Gene* gene, Gene* radGene, Point* points, float** graph, int NofPoints)
{
	int C1 = 1;
	int C2 = 2;
	int C3 = 8;

	vector<int>** tree = new vector<int>*[NofPoints];
	//cout << 4 << endl;
	geneToTree(gene, graph, NofPoints, tree, NofPoints);
	//cout << 5 << endl;


	double meanpathlength = meanPathLength(tree, points, NofPoints);
	double totalpathlength = totalPathLength(tree, points, NofPoints);
	double totalImpedance = 0;
	if (OPTFORRAD) totalImpedance = getNetworkImpedence(tree, radGene, points, NofPoints);

	double cost = C1 * totalpathlength + C2 * meanpathlength + C3* totalImpedance;
	//cout << totalpathlength << "  " << meanpathlength << "  "<< cost << endl;
	for (int i = 0; i < NofPoints; i++)
	{
		delete tree[i];
	}
	delete[] tree;
	double fit = 1 / cost; // The larger the cost, the smaller the fitness
	return fit;
}
