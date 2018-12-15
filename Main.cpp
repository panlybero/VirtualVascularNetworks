#include"PrimModel.h"
#include"AugmentedGenome.h"
#include<iostream>
#include"TreeGene.h"
#include<time.h>
#include<ctime>

AugmentedGenome* treeToGenome(vector<int>** tree, Point* points, int size)
{
	AugmentedGenome* ag = new AugmentedGenome();
	vector<NodeGene*> ngenes;
	for (int i = 0; i < size; i++)
	{
		Point* p;
		vector<Point*>* out  = new vector<Point*>();
		NodeGene* ng;
		for (int j = 0; j < tree[i]->size(); j++)
		{
			out->push_back(&points[tree[i]->at(j)]);
		}

		if (i == 0) 
		{
			p = &points[0];
			ng = new NodeGene(p, &points[i], out);
		}
		else
		{
			for (int k = 0; k < size; k++)
			{
				for (int l = 0; l < tree[k]->size(); l++)
				{
					if (tree[k]->at(l) == i)
					{
						p = &points[i];
						ng = new NodeGene(p, &points[i], out);
					}
				}
			}
		}


	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < tree[i]->size(); j++)
		{

		}

	}

	return ag;
}
int population = 150;
int generations = 50000;
int rate = 20;
double side = 5;
std::random_device rd1{};
std::mt19937 gen1{ rd1() };
TreeGene* best;
TreeGene* first;
void runTreeGeneModel(Point* points, int n)
{
	string totalString = "";
	string allgens = "";
	totalString += to_string(n) + "\r\n";
	vector<TreeGene*> pool;
	vector<TreeGene*> newPool;
	//Generate Random Population
	uniform_real_distribution<float> probdist(0, 1);

	double initfit = 0;
	double finalmaxfit = 0;
	for (int i = 0; i < population; i++)
	{
		pool.push_back(new TreeGene(n));
		pool.back()->randomizeTree();
	}

	for (int i = 0; i < generations; i++)
	{
		
		if (i % 100 == 0)
		{
			system("CLS");
			cout << i << endl;
		}
		vector<double> probabilities;
		double totalfit = 0;
		double maxfit = -1;
		int maxfiti = 0;

		for (int j = 0; j < population; j++)
		{
			
			double currfit = getFitness(pool.at(j), points, n);
			
			totalfit += currfit;
			if (currfit > maxfit)
			{
				maxfit = currfit;
				maxfiti = j;
			}
			

		}

		if (i == 0) initfit = maxfit;
		for (int j = 0; j < population; j++)
		{
			probabilities.push_back((pool.at(j)->fitness) / totalfit);
		
		}
		//cout << "after probs" << endl;

		std::discrete_distribution<int> distribution(probabilities.begin(),probabilities.end()); 

		int ind1 = distribution(gen1);
		int ind2 = ind1;
		
		delete best;

		TreeGene* elite1 = new TreeGene(pool.at(maxfiti));
		TreeGene* elite2 = new TreeGene(pool.at(maxfiti));
		elite2->mutate();

		best = new TreeGene(elite1);
		allgens += to_string(best->fitness) + "\r\n";

		if (i % 100 == 0)
		{
			
			best->generation = i;
			totalString += stringTree(best, points, n) + "\r\n";
		}
		if (i == 0) first = new TreeGene(elite1);

		newPool.push_back(elite1);
		newPool.push_back(elite2);
		int count = population - 2;
		while (count >= 0)
		{
			ind1 = distribution(gen1);
			while (ind2 == ind1)
			{
				ind2 = (distribution(gen1));
			}
			TreeGene* p1 = pool.at(ind1);
			TreeGene* p2 = pool.at(ind2);
			pair<TreeGene*, TreeGene*> p = makeChildren(p1, p2);

			if (probdist(gen1) <= p.first->mutrate)	for (int q = 0; q < 1; q++)p.first->shuffleSubTree();
			if (probdist(gen1) <= p.second->mutrate) for (int q = 0; q < 1; q++)p.second->shuffleSubTree();
			newPool.push_back(p.first);
			newPool.push_back(p.second);
			count -= 2;

		}

		/*
		while (ind2 == ind1)
		{
			ind2 = (distribution(gen1));
		}

		TreeGene* p1 = pool.at(ind1);
		TreeGene* p2 = pool.at(ind2);

		

		for (int j = 0; j < (population-2) / 2; j++)
		{
			pair<TreeGene*, TreeGene*> p = makeChildren(p1, p2);
			
			if(probdist(gen1)<=p.first->mutrate)	for(int q=0;q<2;q++)p.first->mutateAnyPoint();
			if (probdist(gen1) <= p.second->mutrate) for (int q = 0; q<2; q++)p.second->mutateAnyPoint();
			newPool.push_back(p.first);
			newPool.push_back(p.second);
		}
		*/
		//cout << "out of mut" << endl;
		while (!pool.empty())
		{
			TreeGene* temp = pool.back();
			pool.pop_back();
			delete temp;
		}
		for (int k = 0; k < newPool.size(); k++)
		{
			pool.push_back(newPool.at(k));
	
		}
		while (!newPool.empty())
		{
			newPool.pop_back();
		}

		if (finalmaxfit < maxfit) finalmaxfit = maxfit;
		//cout << maxfit << endl;
		//cout << "end of gen" << endl;
	}

	cout << "Init fit = " << initfit << endl << "Final max fit = " << finalmaxfit << endl;

	//printTree("D:\\Generations\\first.txt", first, points, n);

	printTree("D:\\Generations\\final.txt", best, points, n);
	//FileHandler f("D:\\Generations\\All.txt");
	//f.print(totalString);
	//FileHandler f2("D:\\Generations\\gens.txt");
	//f2.print(allgens);
}


int main()
{
	int n = 0;

	Point* points = getDavidPoints(side, 1.0 / sqrt(acos(-1.0)), n);
	//Point* points = getExistingPoints("best.txt", n);
	/*
	srand(time(NULL));
	int start_s = clock();
	best = new TreeGene(10);
	
	runTreeGeneModel(points,n);

	cout << "leaks = " << gettheLeaks() << endl;

	int stop_s = clock();
	cout << "time: " << ((stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000) / 1000 << " sec" << endl;
	

	system("pause");
	
	
	cout << endl;
	*/
	int start_s = clock();

	mainModel(points, n, side);
	

	int stop_s = clock();
	cout << "time: " << ((stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000) / 1000 << " sec" << endl;
	

	/*
	AugmentedGenome ag;
	int v = 29;
	Point* points = getExistingPoints("D:\\file2.txt", v);

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

	Gene* g = generateRandomGene(v*(v - 1) / 2);

	vector<int>** vec = new vector<int>*[v];
	geneToTree(g, graph, v, vec, v);
	treeToGenome(vec, points, v);
	system("pause");


	ag.mutSwapConnection();

	for (int i = 0; i < ag.genes.size(); i++)
	{

		cout << "(" << ag.genes.at(i)->getUp()->getThis()->x<< "," << ag.genes.at(i)->getUp()->getThis()->y << ")"<< " --> " << "(" << ag.genes.at(i)->getDown()->getThis()->x << "," << ag.genes.at(i)->getDown()->getThis()->y << ")" << endl;
	}


	system("pause");
	*/
}