#include "Network.h"
#include <math.h> 
#include<algorithm>

#define THREEDEE false
int leaky = 0;
float heartSize = 100;

float myclamp(float v, float hi, float lo)
{
	return (v < lo) ? lo : (v > hi) ? hi : v;
}

Network::Network(int nOfPoints, Point* pnts)
{
	points = new Point[nOfPoints];
	
	for (int i = 0; i < nOfPoints; i++)
	{
		points[i] = Point(pnts[i]);
	}

	n = nOfPoints;
	geneLength = nOfPoints*(nOfPoints - 1) / 2;
	gene = generateRandomGene(geneLength);
	radGeneLength = n - 1;
	radGene = generateRandomRadGene(radGeneLength);
	
	gene->setMutrate(0.15);// myclamp(geneLength*0.001, 0.5, 0.1));
	radGene->setMutrate(0.15);//myclamp(radGeneLength*0.001, 0.5, 0.1));



	fit = 0;
	tpl = 0;
	mpl = 0;
	ztot = 0;
	cons = 0;
	leaky++;
}

Network::Network(Network * other) // Makes a copy (expensive...) 
{
	n = other->n;
	geneLength = other->geneLength;
	radGeneLength = other->getRadGeneLength();
	fit = other->fit;
	tpl = other->tpl;
	mpl = other->mpl;
	ztot = other->ztot;
	cons = other->cons;




	delete[] points;
	points = new Point[n];
	for (int i = 0; i < n; i++)
	{
		points[i] = Point(other->points[i]);
	}
	delete gene;
	gene = new Gene(geneLength);
	copyGene(gene, other->gene);
	delete radGene;
	radGene = new Gene(radGeneLength);
	copyGene(radGene, other->radGene);

	leaky++;
}

Network::~Network()
{
	delete[] points;
	delete gene;
	delete radGene;
	leaky--;
}

void Network::mutate()
{
	if (randomDouble(0, 1) < gene->getMutrate())
	{
		
		int mutpos = randomInt(0, geneLength-1);
		gene->setPos(randomDouble(0, heartSize), mutpos);
		//int mutperc = randomInt(1, 200);
		//double newval = max(gene->getPos(mutpos)*mutperc / 100, 1.0f);
		//gene->setPos(newval,mutpos);
		//double newmut = 1.0 / double(randomInt(1, 10));
		//gene->setMutrate(0.15);//newmut);//randomDouble(0, 1));
	}
	if (randomDouble(0, 1) < radGene->getMutrate())
	{

		int mutpos = randomInt(0, radGeneLength-1);

		int mutperc = randomInt(1, 200);
		//radGene->dump();
		//cout << mutpos << endl;
		radGene->setPos(randomDouble(1,heartSize), mutpos);
		//radGene->dump();
		//system("pause");

		//double newval = min(max((radGene->getPos(mutpos)*mutperc / 100), 1.0f),heartSize);
		//radGene->setPos(newval, mutpos);

		//double newmut = 1.0 / double(randomInt(1, 10));
		//radGene->setMutrate(0.15);//newmut);//randomDouble(0, 1));
	}
	if (randomDouble(0, 1) > 1.9)
	{
		double newmut = 1.0 / double(randomInt(1, 10));
		radGene->setMutrate(newmut);//randomDouble(0, 1));
		gene->setMutrate(newmut);
	}

}


void Network::restrictedMutate()
{
	vector<int>** tree = geneToTree(this);
	

	if (randomDouble(0, 1) < gene->getMutrate())
	{

		int mutpos = randomInt(0, geneLength - 1);
		gene->setPos(randomDouble(0, heartSize), mutpos);

	}
	if (randomDouble(0, 1) < radGene->getMutrate())
	{

		int mutpos = randomInt(0, radGeneLength - 1);
		for (int i = 0; i < getN(); i++)
		{
			for (int j = 0; j < tree[i]->size(); j++)
			{
				if (tree[i]->at(j) == mutpos)
				{
					radGene->setPos(randomDouble(1, (i==0)? heartSize:radGene->getPos(i-1)), mutpos);
				}
			}
		}
		

	}
	for (int i = 0; i < getN(); i++)
	{
		delete tree[i];
	}
	delete[] tree;
}
int Network::getGeneLength()
{
	return geneLength;
}

int Network::getRadGeneLength()
{
	return radGeneLength;
}

int Network::getN()
{
	return n;
}

double Network::getFitness(string mode) // Calculates and saves fitness in fit attribute 
{

	double costScale = 0.001;
	double C1 =  0.1;
	double C2 = 0;
	double C3 =  1;
	double C4 =  0.1;//1;
	vector<int>** tree = geneToTree(this);
	
	//enforceAreaConservation(tree, this, getN());
	
	double totalpathlength = 0;// totalPathLength(tree, points, getN());
	double meanpathlength = 0;// meanPathLength(tree, points, getN());
	double totalimpedence = 0;// getNetworkImpedence(tree, this, getN());
	double conservation = 0;// getAreaPreservation(tree, this, getN());
	//cout << C1*totalpathlength/(C3*totalimpedence) << endl;
	//system("pause");
	

	if (C1 != 0)  totalpathlength =  totalPathLength(tree, points, getN());
	if (C2 != 0) meanpathlength =  meanPathLength(tree, points, getN());
	if (C3 != 0) totalimpedence = getNetworkImpedence(tree, this, getN(), mode);
	if (C4 != 0) conservation = getConservation(tree, this, getN(),mode);
	 

	//if(randomDouble(0,1)>0.9999) cout << C1*totalpathlength << " + " << C2*meanpathlength << " + " << C3*totalimpedence << " + " << C4* conservation << endl;
	
	double cost = costScale * (C1*totalpathlength + C2*meanpathlength + C3*totalimpedence + C4* conservation);
	

	//cout << C1*totalpathlength / (C3*totalimpedence) << endl;
	//system("pause");
	fit = 1/cost;

	tpl = C1* totalpathlength;
	mpl = C2*meanpathlength;
	ztot = C3*totalimpedence;
	cons = C4* conservation;

	for (int i = 0; i < getN(); i++)
	{
		delete tree[i];
	}
	delete[] tree;

	return fit;
}

pair<Network*, Network*> cross(Network * p1, Network * p2) //Cross and produce 2 offspring per crossing
{

	Network* ch = new Network(p1);
	Network* ch2 = new Network(p2);
	double crossprob = 0.6;
	int pos = 0;
	//Gene
	if (randomDouble(0, 1) > crossprob)
	{


		pos = randomInt(0, ch->getGeneLength() - 1);
		int pos2 = randomInt(pos, ch->getGeneLength());
		for (int i = pos; i < ch->getGeneLength(); i++)
		{
			ch->gene->setPos(p2->gene->getPos(i), i);
			ch2->gene->setPos(p1->gene->getPos(i), i);
		}
	}
	if (randomDouble(0, 1) > crossprob)
	{
		//RadGene
		pos = randomInt(0, ch->getRadGeneLength());

		for (int i = pos; i < ch->getRadGeneLength(); i++)
		{
			ch->radGene->setPos(p2->radGene->getPos(i), i);
			ch2->radGene->setPos(p1->radGene->getPos(i), i);
		}
	}

	return make_pair(ch,ch2);
}

vector<int>** geneToTree(Network* network)
{
	/* Let us create the following graph
	2    3
	(0)--(1)--(2)
	|   / \   |
	6| 8/   \5 |7
	| /     \ |
	(3)-------(4)
	9          */
	int NofPoints = network->getN();


	const int geneLength = NofPoints*(NofPoints - 1) / 2;
	//double gene[geneLength];


	int index = NofPoints;
	//making tree 0-1 0-3, 1-2, 3-4


	float** graphprime = new float*[NofPoints];



	for (int i = 0; i < NofPoints; i++)
	{
		graphprime[i] = new float[NofPoints];
		for (int j = 0; j < NofPoints; j++)
		{
			graphprime[i][j] = 1;
			if (i == j) graphprime[i][j] = 0.0f;
			
		}

	}
	index = 0;
	for (int i = 0; i < NofPoints; i++) // Look at "Representing trees with genetic algorithms" 
	{
		for (int j = i + 1; j < NofPoints; j++)
		{

			graphprime[i][j] += 1 * network->gene->getPos(index);
			graphprime[j][i] += 1 * network->gene->getPos(index);
			index++;
		}
	}

	// Print the solution
	int* parent = new int[NofPoints];
	primMST(graphprime, parent, NofPoints);

	
	vector<int>** tree = new vector<int>*[NofPoints];

	//vector<int>* tree[V];
	for (int i = 0; i < NofPoints; i++)
	{
		tree[i] = new vector<int>();

	}

	for (int i = 1; i < NofPoints; i++)
	{
	
		tree[parent[i]]->push_back(i);
	}

	for (int i = 0; i < NofPoints; i++)
	{

		delete[] graphprime[i];

	}
	delete[] graphprime;
	delete[] parent;
	return tree;
}
int minKey(int key[], bool mstSet[], int n)
{
	// Initialize min value
	int min = INT_MAX, min_index;

	for (int v = 0; v < n; v++)
		if (mstSet[v] == false && key[v] < min)
			min = key[v], min_index = v;

	return min_index;
}


// A utility function to print the constructed MST stored in parent[]
void printMST(int parent[], int n, float** graph)
{
	printf("Edge   Weight\n");
	for (int i = 1; i < n; i++)
		printf("%d - %d    %f \n", parent[i], i, graph[i][parent[i]]);
}

void primMST(float** graph, int* parent, int n)
{
	//int parent[V]; // Array to store constructed MST
	int* key = new int[n];   // Key values used to pick minimum weight edge in cut
	bool* mstSet = new bool[n];  // To represent set of vertices not yet included in MST

								 // Initialize all keys as INFINITE
	for (int i = 0; i < n; i++)
		key[i] = INT_MAX, mstSet[i] = false;

	// Always include first 1st vertex in MST.
	key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
	parent[0] = -1; // First node is always root of MST 

					// The MST will have V vertices
	for (int count = 0; count < n - 1; count++)
	{
		// Pick the minimum key vertex from the set of vertices
		// not yet included in MST
		int u = minKey(key, mstSet, n);

		// Add the picked vertex to the MST Set
		mstSet[u] = true;

		// Update key value and parent index of the adjacent vertices of
		// the picked vertex. Consider only those vertices which are not yet
		// included in MST
		for (int v = 0; v < n; v++)
		{
			// graph[u][v] is non zero only for adjacent vertices of m
			// mstSet[v] is false for vertices not yet included in MST
			// Update the key only if graph[u][v] is smaller than key[v]
			bool a = graph[u][v];
			bool b = mstSet[v] == false;
			bool c = graph[u][v] < key[v];

			if (a && b && c)
				parent[v] = u, key[v] = graph[u][v];
		}
	}

	// print the constructed MST
	//printMST(parent, V, graph);

	delete[] key;
	delete[] mstSet;

}

double* getNetworkProbs(vector<Network*>* pool, string mode)
{
	double* fitness = new double[pool->size()]; // Essentially useless, may remove. Here only in case we need to look at just fitnesses. 
	double* probs = new double[pool->size()];
	double sum = 0;
	double curr = 0;
	double maxFit = 0;
	for (int i = 0; i < pool->size(); i++)
	{
		curr = pool->at(i)->getFitness(mode);
		sum += curr;
		fitness[i] = curr;

		if (maxFit < fitness[i])
		{
			maxFit = fitness[i];
		}

	}
	double meanFit = sum / pool->size();
	//cout << sum << endl;
	double var = 0;
	for (int i = 0; i < pool->size(); i++)
	{
		//cout << (fitness[i] - meanFit)*(fitness[i] - meanFit) << endl;
		var += (fitness[i] - meanFit)*(fitness[i] - meanFit);
		
	}
	//var /= pool->size();
	//cout << var << endl;
	//cout << s << endl;
	for (int i = 0; i < pool->size(); i++)
	{
		probs[i] = (fitness[i]/sum); 
		//cout << (probs[i]) << endl;

	}

	delete[] fitness;
	return probs;
}

double totalPathLength(vector<int>** tree, Point* points, int v)
{
	double t = 0;
	for (int i = 0; i < v; i++)
	{
		//cout << i << endl;
		for (int j = 0; j < tree[i]->size(); j++)
		{
			t += eucDistance(points[i], points[tree[i]->at(j)]);
		}
	}
	//cout << t << endl;
	return t;
}

double Network::getFit()//Warning, may be 0 or stale because getFitness has not been called yt
{
	return fit;
}

vector<Network*>* nextPopulation(vector<Network*>* pool, int generation, int totalgens, string mode)
{
	//break up using dividePopulation
	//proccess in parallel
	//join the new subpops
	//return that. 
	int population_size = pool->size();

	


	double max_fit = 0;
	int max_fit_i = 0;
	for (int i = 0; i < pool->size(); i++)
	{
		if (pool->at(i)->getFit() > max_fit)
		{
			max_fit = pool->at(i)->getFit();
			max_fit_i = i;
		}
	}
	Network* elite = new Network(pool->at(max_fit_i));
	//elite->gene->dump();
	//system("pause");
	vector<vector<Network*>*> divd = dividePopulation(pool);
	
	vector<future<vector<Network*>*>> nextdivd;
	vector<Network*>* newpool = new vector<Network*>();


	string selectmode = (generation > totalgens/2) ? "tournament" : "fitprop"; // first half tournament, then fitprop
	

	for (int i = 0; i < divd.size(); i++)
	{
		nextdivd.push_back(async(nextSubPopulation, divd.at(i), selectmode, mode));
	}
	/*
	while (!divd.empty())
	{
		vector<Network*>* tmp = divd.back();
		divd.pop_back();
		delete tmp;
	}
	*/

	divd.clear();

	for (int i = 0; i < nextdivd.size(); i++)
	{
		divd.push_back(nextdivd.at(i).get());

	}
	//make newpool through random ordering of subpool members
	while (newpool->size() < population_size)
	{
		int ch = randomInt(0, divd.size()-1);
		if (!divd.at(ch)->empty())
		{
			newpool->insert(newpool->begin(), divd.at(ch)->back());
			divd.at(ch)->pop_back();
		}
		else
		{
			divd.erase(divd.begin() + ch);
		}
		
	}
	while (!divd.empty())
	{
		vector<Network*>* tmp = divd.back();
		divd.pop_back();
		delete tmp;
	}
	Network* back = newpool->back();
	delete back;
	newpool->pop_back();
	newpool->push_back(elite);

	while (!pool->empty())
	{
		Network* tmp = pool->back();
		pool->pop_back();
		delete tmp;

	}
	delete pool;
	
	return newpool;
}

vector<Network*>* nextSubPopulation(vector<Network*>* subpool, string selectmode, string mode)
{
	//SubPop size is 25. 12 Pairs of parents will be chosen to produce 2 offspring each. Plus 1 for elitism = 25
	vector<Network*>* newpool = new vector<Network*>();
	selectmode = "tournament";
	
	//get probs
	//double* probs = getNetworkProbs(subpool);
	double max_parent_fit = 0;
	double max_parent_i = 0;
	double* probs = nullptr;
	if(selectmode == "fitprop")
	probs = getNetworkProbs(subpool, mode);

	for (int i = 0; i < 13; i++)
	{
		//choose parents

		pair<int, int> parents;
		if(selectmode == "fitprop") parents = chooseParentsFitProp(probs, subpool->size());

		if(selectmode == "tournament") parents = chooseParentsTournament(subpool, subpool->size(), mode);

		//max fit among parents
		if (max_parent_fit < subpool->at(parents.first)->getFit())
		{
			max_parent_fit = subpool->at(parents.first)->getFit();
			max_parent_i = parents.first;
		}  
		if (max_parent_fit < subpool->at(parents.second)->getFit())
		{
			max_parent_fit = subpool->at(parents.second)->getFit();
			max_parent_i = parents.second;
		}
		
		//cross
		pair<Network*, Network*> children = cross(subpool->at(parents.first), subpool->at(parents.second));

		//mutate
		children.first->restrictedMutate();
		children.second->restrictedMutate();

		//children.first->mutate();
		//children.second->mutate();
		//add
		newpool->push_back(children.first);

		if (i < 12)
		{
			newpool->push_back(children.second);
			
		}
		else
		{
			delete children.second;
		}

	}
	

	//clear old subpool;
	while (!subpool->empty())
	{
		Network* tmp = subpool->back();
		subpool->pop_back();
		delete tmp;
	}
	delete subpool;

	//return new subpop
	delete[] probs;
	return newpool;
}

pair<int, int> chooseParentsFitProp(double* probs, int size)
{
	double csum = 0;
	double ev = randomDouble(0, 1);
	

	int p1=0;
	int p2 = 0;
	for (int i = 0; i < size; i++)
	{
		csum += probs[i];
		if (ev < csum)
		{
			p1 = i;
			break;
		}
	}
	csum = 0;
	for (int i = 0; i < size; i++)
	{
		csum += probs[i];
		if (ev < csum)
		{
			p2 = i;
			break;
		}
	}
	
	return make_pair(p1, p2);
}

pair<int, int> chooseParentsTournament(vector<Network*>* subpool, int size, string mode)
{

	double k = 0.75;
	int p1 = 0;
	int p2 = 0;
	int c1 = randomInt(0, size - 1);
	int c2 = randomInt(0, size - 1);

	double fit_c1 = subpool->at(c1)->getFitness(mode);
	double fit_c2 = subpool->at(c2)->getFitness(mode);

	if (fit_c1 > fit_c2)
	{
		if (randomDouble(0, 1) < k)
		{
			p1 = c1;
		}
	}
	c1 = randomInt(0, size - 1);
	c2 = randomInt(0, size - 1);
	if (fit_c1 > fit_c2)
	{
		if (randomDouble(0, 1) < k)
		{
			p2 = c1;
		}
	}
	return make_pair(p1,p2);
}

vector<vector<Network*>*> dividePopulation(vector<Network*>* pool)
{
	//return a vector of subpopulations. 

	int population = pool->size(); 
	int subpopulation_size = 25;
	vector<vector<Network*>*> subpops;
	int nOfsub = population/subpopulation_size; //INTEGER DIVISION  	//std subpopulation size = 25 unless there is remainder;
	int lastSubSize = population - nOfsub * subpopulation_size;

	if (lastSubSize == 0) lastSubSize = subpopulation_size;

	for (int i = 0; i < nOfsub-1; i++)
	{
		subpops.push_back(new vector<Network*>);
		while (subpops.back()->size() < 25)
		{
			int choose = randomInt(0, pool->size()-1);//pool->size() changes in every iteration. 
			subpops.back()->push_back(pool->at(choose));
			vector<Network*>::iterator nth = pool->begin() + choose;
			pool->erase(nth);//removes from pool the one added to population
		}
		
	}
	subpops.push_back(new vector<Network*>);

	while (subpops.back()->size() < lastSubSize)
	{
		int choose = randomInt(0, pool->size()-1);//pool->size() changes in every iteration. 
		subpops.back()->push_back(pool->at(choose));
		vector<Network*>::iterator nth = pool->begin() + choose;
		pool->erase(nth);//removes from pool the one added to population
	}

	
	return subpops;
}

double getPulsingImp(Network* net, int curr)
{
	double imp = pow(net->radGene->getPos(curr - 1),2);
	return imp;

}

double getConstantImp(Network* net, pair<int,int> pc)
{
	double length = eucDistance(net->points[pc.first], net->points[pc.second]);
	
	//cout << pc.second << endl;

	double rad = net->radGene->getPos(pc.second - 1);

	if (rad == -1.0)
	{
		cout << "failed at " << pc.second << endl;
		cout << "length " << length << endl;
		system("pause");
	}
	double imp = length * pow(rad, -4);
	return imp;
}

double calcConstFlow(vector<int>** tree, Network* net, int curr, int par)
{
	double imp = 0;

	if (curr != 0)
		if (tree[curr]->size() == 0)
		{
			imp = getConstantImp(net, make_pair(par,curr));
			/*
			if (isinf(imp))
			{
			cout << "hi" << endl;
			}
			*/
			return  imp;
		}

	imp = (curr == 0) ? 100 : getConstantImp(net, make_pair(par,curr)); // come back here
	
	double temp = 0;

	if (tree[curr]->size() == 1)
	{
		imp += calcConstFlow(tree, net, tree[curr]->at(0),curr); //series
																 //imp/=

	}
	else
	{
		for (int i = 0; i < tree[curr]->size(); i++)
		{
			double k = calcConstFlow(tree, net, tree[curr]->at(i),curr);

			temp += 1 / k;

		}

		imp += 1 / temp;
	}
	/*
	if (isinf(imp))
	{
	for (int j = 0; j < net->getN(); j++)
	{
	cout << j;
	for (int l = 0; l < tree[j]->size(); l++)
	{
	cout << " " << tree[j]->at(l) << " ";
	}
	cout << endl;
	}
	net->radGene->dump();
	cout << curr << endl;
	cout << (net->radGene->getPos(curr - 1)*net->radGene->getPos(curr - 1)) << endl;
	system("pause");
	}
	*/
	return imp;
	/*
	double imp = 0;
	double tmp = 0;
	for (int i = 0; i < tree[curr]->size(); i++)
	{
		if (tree[i]->size() == 0)
		{
			tmp += 1 / getConstantImp(net, make_pair(curr, tree[curr]->at(i)));
		}
		else
		{
			tmp += 1 / calcConstFlow(tree, net, tree[curr]->at(i), curr);
		}
	}


	imp = (curr == 0) ? heartSize*heartSize : getConstantImp(net, make_pair(par,curr));
	double temp = 0;
	if (tree[curr]->size() == 1)
	{
		imp += calcConstFlow(tree, net, tree[curr]->at(0),curr); //series
																 //imp/=

	}
	else
	{
		for (int i = 0; i < tree[curr]->size(); i++)
		{
			double k = calcPulsingFlowImp(tree, net, tree[curr]->at(i));

			temp += 1 / k;

		}

		imp += 1 / temp;
	}
	return imp;
	*/
}

double calcPulsingFlowImp(vector<int>** tree, Network* net, int curr)
{

	double imp = 0;

	if(curr!=0)
	if (tree[curr]->size() == 0)
	{
		imp = getPulsingImp(net,curr);
		/*
		if (isinf(imp))
		{
			cout << "hi" << endl;
		}
		*/
		return  imp;
	}

	imp = (curr == 0) ? 0: getPulsingImp(net, curr);//heartSize*heartSize : getPulsingImp(net,curr);
	double temp = 0;
	if (tree[curr]->size() == 1)
	{
		net->radGene->setPos((curr ==0 )? heartSize: net->radGene->getPos(curr - 1), tree[curr]->at(0) - 1);
		imp += calcPulsingFlowImp(tree, net, tree[curr]->at(0)); //series
																	//imp/=
		
	}
	else
	{
		for (int i = 0; i < tree[curr]->size(); i++)
		{
			double k = calcPulsingFlowImp(tree, net, tree[curr]->at(i));

			temp += 1 / k;
			
		}
		
		imp += 1 / temp;
	}
	/*
	if (isinf(imp))
	{
		for (int j = 0; j < net->getN(); j++)
		{
			cout << j;
			for (int l = 0; l < tree[j]->size(); l++)
			{
				cout << " " << tree[j]->at(l) << " ";
			}
			cout << endl;
		}
		net->radGene->dump();
		cout << curr << endl;
		cout << (net->radGene->getPos(curr - 1)*net->radGene->getPos(curr - 1)) << endl;
		system("pause");
	}
	*/
	return imp;
}

double calcPulsingFlowImp(vector<int>** tree, Point* points, int curr)
{
	
	if (tree[curr]->size() == 0)
	{
		return  (points[curr].rad*points[curr].rad);
	}
	double imp = points[curr].rad * points[curr].rad;
	double temp = 0;
	if (tree[curr]->size() == 1)
	{
		imp += calcPulsingFlowImp(tree, points, tree[curr]->at(0)); //series
																	//imp/=
	}
	else
	{
		for (int i = 0; i < tree[curr]->size(); i++)
		{
			double k = calcPulsingFlowImp(tree, points, tree[curr]->at(i));;
			temp += 1 / k;
			
		}
		imp += 1 / temp;
	}
	
	return imp;
}

double localImpedanceMatching(vector<int>** tree, Network* network, int NofPoints, string mode)
{
	double impdiff = 0;
	double a_curr = 0;
	double a_children = 0;

	double deg = 0;
	if (mode == "pulsing") deg = 2;
	if (mode == "constant") deg = 3;

	for (int i = 0; i < NofPoints; i++)
	{
		a_curr = (i == 0) ? pow(heartSize,deg): pow(network->radGene->getPos(i - 1), deg);

		if (tree[i]->size() == 1)
		{
			network->radGene->setPos((i==0)? heartSize: network->radGene->getPos(i - 1), tree[i]->at(0) - 1);
		}

		for (int j = 0; j < tree[i]->size(); j++)
		{
			double tmp = pow(network->radGene->getPos(tree[i]->at(j) - 1), deg); //check if must multiply by length
			//if (mode == "constant") tmp *= eucDistance(network->points[i], network->points[tree[i]->at(j)]);
			a_children += tmp;//pow(network->radGene->getPos(tree[i]->at(j) - 1), deg);
		}
		
		impdiff += abs(a_curr - a_children)/a_curr;
		a_children = 0;
	}
	//if(randomDouble(0,1)>0.9999) cout << impdiff << endl;
	return impdiff;
}

double getNetworkImpedence(vector<int>** tree, Gene* radGene, Point* points, int NofPoints)
{
	//Distribute rad values from gene

	assignRadsToPoints(tree, points, NofPoints, radGene);

	//Get Impedence
	double imp = calcPulsingFlowImp(tree, points, 0);

	

	return imp;
}

double getConservation(vector<int>** tree, Network* network, int NofPoints, string mode)
{
	double apres = 0;
	if(mode == "pulsing") apres = localImpedanceMatching(tree, network, NofPoints,"pulsing");
	if(mode == "constant") apres = localImpedanceMatching(tree, network, NofPoints, "constant");
	
	return apres;
}

void enforceAreaConservation(vector<int>** tree, Network* network, int NofPoints)
{
	for (int i = 0; i < NofPoints; i++)
	{
		double parentrad = (i == 0)? heartSize : network->radGene->getPos(i-1);
		for (int j = 0; j < tree[i]->size(); j++)
		{
			double myrad = network->radGene->getPos(tree[i]->at(j) - 1);
			network->radGene->setPos(myclamp(myrad,parentrad,1), tree[i]->at(j) - 1);
		}
	}


	/*
	for (int i = 0; i < NofPoints; i++)
	{
		int daughters = tree[i]->size();
		double parentArea = (i == 0) ? pow(heartSize, 2) : pow(network->radGene->getPos(i - 1), 2);
		double darea = 0;
		//cout << daughters << endl;
		if (daughters == 1)
		{
			network->radGene->setPos(sqrt(parentArea), tree[i]->at(0) - 1);
			continue;
		}
		for (int j = 0; j < daughters - 1; j++)
		{
			darea+= pow(network->radGene->getPos(tree[i]->at(j) - 1), 2);
		}
		
		double lastDaughterArea = parentArea - darea;
		
		double lastDaughterRad = (lastDaughterArea<=1)? 1:sqrt(lastDaughterArea);
		if(daughters>1)
		network->radGene->setPos(lastDaughterRad, tree[i]->at(daughters - 1)-1);

	}
	*/
}

bool sameTree(vector<int>** tree1, vector<int>** tree2,int size)
{
	bool same = false;
	for (int i = 0; i < size; i++)
	{
		if (tree1[i]->size() != tree2[i]->size())
		{
			same = false;

			//cout << "different" << endl;
			return same;
		}
		for (int j = 0; j < tree1[i]->size(); j++)
		{

			if (tree1[i]->at(j) != tree2[i]->at(j))//maybe different order?
			{

				same = false;
				//cout << "different" << endl;
				return same;
				
			}
		}

	}
	//cout << "same" << endl;
	same = true;
	return same;
}

int countNiches(vector<Network*>* pool, int NofPoints)
{
	
	int n_niches = 0;
	int n = pool->size();
	bool* checked = new bool[n];
	bool madeNiche = false;

	for (int i = 0; i < n; i++)
	{
		checked[i] = false;
	}
	vector<int>** tree1 = nullptr; 
	vector<int>** tree2 = nullptr;
	vector<vector<int>**> niches;
	niches.reserve(pool->size());
	int treeSize = pool->at(0)->getN();
	for (int i = 0; i < n; i++)
	{
		tree1 = geneToTree(pool->at(i));
		if (i == 0)
		{
			niches.push_back(tree1);
		}else
		{
			bool found = false;
			for (int j = 0; j < niches.size(); j++)
			{
				if (!sameTree(niches.at(j), tree1, treeSize))
				{
					niches.push_back(tree1);
					break;
				}
			}
		}
			
	}
	/*
	cout << "checked" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << checked[i] << endl;
	}
	*/

	return niches.size();
}

double getNetworkImpedence(vector<int>** tree, Network* network, int NofPoints, string mode)
{
	//Distribute rad values from gene

	//assignRadsToPoints(tree, network->points, NofPoints, network->radGene);
	
	/*
	if (randomDouble(0, 1) > 0.999)
	{
		for (int i = 0; i < NofPoints; i++)
		{
			network->points[i].dump();
			cout << endl;
		}
		system("pause");
	}
	*/
	//Get Impedence
	double imp = 0;
	if(mode == "pulsing") imp = calcPulsingFlowImp(tree, network, 0);
	if (mode == "constant") imp = calcConstFlow(tree, network, 0, 0);
	//system("pause");

	if (isinf(imp))
	{
		cout << "a" << endl;
	}
	return imp;
}

int getLeaky()
{
	return leaky;
}

void printTree(string path, vector<int>** tree, Network* net, Point* points, int v)
{
	FileHandler f(path);

	f.print(to_string(v));

	if (!THREEDEE)
		for (int i = 0; i < v; i++)
		{
			double rad = (i == 0) ? heartSize : net->radGene->getPos(i-1);
			f.print(to_string(i) + " " + to_string(points[i].x) + " " + to_string(points[i].y) + " " + to_string(rad));
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
	///////////////////////////////////////////////////////////


}