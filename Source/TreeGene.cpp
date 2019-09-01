#include "TreeGene.h"
#define THREEDEE false

int theleaks = 0;
std::random_device rd{};
std::mt19937 gen{ rd() };

TreeGene::TreeGene(vector<int>** tree,int n)
{
	m_tree = tree;
	m_n = n;
	theleaks++;
}

TreeGene::TreeGene(int n)//Make empty tree
{
	m_tree = new vector<int>*[n];
	m_n = n;
	for (int i = 0; i < n; i++)
	{
		m_tree[i] = new vector<int>();
	}
	radGene = new Gene(n - 1);
	theleaks++;
}

TreeGene::TreeGene(TreeGene * other)
{
	theleaks++;
	m_n = other->size();
	fitness = other->fitness;
	m_tree = new vector<int>*[m_n];
	Gene* g = new Gene(m_n - 1);
	for (int i = 0; i < m_n; i++)
	{
		m_tree[i] = new vector<int>();
		vector<int>** tmp = other->getTree();
		for (int j = 0; j < tmp[i]->size(); j++)
		{
			m_tree[i]->push_back(tmp[i]->at(j));
		}
	}

	for (int i = 0; i < m_n - 1; i++)
	{
		g->setPos(other->radGene->getPos(i), i);
	}
	radGene = g;
}

TreeGene::~TreeGene()
{
	theleaks--;
	for (int i = 0; i < m_n; i++)
	{
		delete m_tree[i];
	}
	delete[] m_tree;
	delete radGene;
}

void TreeGene::randomizeTree()
{

	float** graph = new float*[m_n];
	for (int i = 0; i < m_n; i++)
	{
		graph[i] = new float[m_n];
	}
	//Setup Graph

	for (int i = 0; i < m_n; i++)
	{
		for (int j = 0; j < m_n; j++)
		{
			graph[i][j] = 1.0f;
			if (i == j) graph[i][j] = 0.0f;
		}
	}

	Gene* g = generateRandomGene(m_n*(m_n - 1) / 2);
	for (int i = 0; i < m_n; i++)
	{
		delete m_tree[i];
	}


	geneToTree(g, graph, m_n, m_tree, m_n);
	delete radGene;

	radGene = generateRandomGene(m_n - 1);


	for (int i = 0; i < m_n; i++)
	{
		delete[] graph[i];
	}
	delete[] graph;
	
}

void TreeGene::dump()
{
	for (int i = 0; i < m_n; i++)
	{
		cout << i << " ";
		for (int j = 0; j < m_tree[i]->size(); j++)
		{
			cout << m_tree[i]->at(j)<<" ";
		}
		cout << endl;
	}
	cout << endl;
}

int TreeGene::size()
{
	return m_n;
}

void TreeGene::addNode(int i, int v)
{
	if (i >= m_n) cout << "addNode out of range" << endl;
	m_tree[i]->push_back(v);
}

vector<int>* TreeGene::getNode(int i)
{
	return m_tree[i];
}

int TreeGene::findintree(int v) // return the parent of v
{

	for (int i = 0; i < m_n; i++)
	{
		
		if (find(m_tree[i]->begin(), m_tree[i]->end(), v) != m_tree[i]->end())
		{
			return i;
		}
		
	}
	 return -1;
}

vector<int>* TreeGene::findintree(vector<int>* v) //return the parents of all vertices in v
{

	queue<int> q;
	q.push(0);
	vector<int>* res = new vector<int>();
	if (v->size() == 0) return res;
	vector<int>::iterator it;
	for (int i = 0; i < v->size(); i++)
	{
		res->push_back(-1);
	}
	int found = 0;
	for (int i = 0; i < m_n; i++)
	{
		for (int j = 0; j < m_tree[i]->size(); j++)
		{
			it = find(v->begin(), v->end(), m_tree[i]->at(j));
			if ( it != v->end())
			{
				res->at(it-v->begin()) = i;
				found++;
			}
		}
	}
	/*
	while (!q.empty())
	{
		int curr = q.front();
		cout << curr << endl;
		q.pop();
		for (int j = 0; j < m_tree[curr]->size(); j++)
		{
			q.push(m_tree[curr]->at(j));

			for (int i = 0; i < v->size(); i++)
			{
				if (v->at(i) == m_tree[curr]->at(j))
				{
					res->at(i) = curr;
					found++;
				}
			}
		}

	}
	*/
	/*
	//cout << "v size = " << v->size() << " found = " << found;
	//cout << " v = ";
	for (int i = 0; i < v->size(); i++)
	{
		cout << v->at(i) << " ";
	}
	cout << endl << "res = ";
	for (int i = 0; i < res->size(); i++)
	{
		cout << res->at(i) << " ";
	}
	cout << endl;
	*/
	return res;
}

bool TreeGene::isConnected()
{
	queue<int> q;
	int counter = 0;
	q.push(0);
	counter++;

	while (!q.empty())
	{
		int curr = q.front();
		q.pop();

		for (int i = 0; i < m_tree[curr]->size(); i++)
		{
			q.push(m_tree[curr]->at(i));
			counter++;
		}
	}

	if (counter == m_n)
		return true;
	return false;
}

vector<int>** TreeGene::getTree()
{
	return m_tree;
}

Gene * TreeGene::getRadGene()
{
	return radGene;
}



void TreeGene::mutate()
{

	srand(time(NULL));
	uniform_int_distribution<> dist1(0, m_n - 2);
	uniform_int_distribution<> dist2(0, m_n - 1);
	for (int i = 0; i < 1; i++)
	{


		int radpos = dist1(gen);

		radGene->setPos(randomDouble(0, 5), radpos);

		mutrate = rand() % 100;

		int mutnode = dist2(gen);
		//cout << mutnode << endl;
		while (m_tree[mutnode]->size() != 0)
		{
			mutnode = dist2(gen);
		}


		int parent = findintree(mutnode);

		vector<int>* parentnode = getNode(parent);
		int keep = 0;
		for (int i = 0; i < parentnode->size(); i++)
		{
			if (parentnode->at(i) == mutnode)
			{
				keep = i;
				break;
			}
		}

		parentnode->erase(parentnode->begin() + keep);

		int newparent = dist2(gen);

		while (m_tree[newparent]->size() == 0 || newparent == mutnode)
		{
			newparent = dist2(gen);
		}
		//cout << "before temp" << endl;

		int tmp = rand() % m_tree[newparent]->size();
		//cout << "after temp" << endl;
		//cout << "before switch with" << endl;

		int switchwith = m_tree[newparent]->at(tmp);
		//cout << "after switch with" << endl;
		//cout << "before add" << endl;
		addNode(mutnode, switchwith);
		addNode(newparent, mutnode);
		//cout << "after add" << endl;
		//cout << "before erase" << endl;
		m_tree[newparent]->erase(m_tree[newparent]->begin() + tmp);
		//cout << "after erase" << endl;
	}

}

void TreeGene::mutateAnyPoint()
{
	//srand(time(NULL));
	uniform_int_distribution<> dist(0,m_n-1);
	uniform_real_distribution<float> dist2(0,1);
	mutrate = dist2(gen);

	int x = dist(gen);
	int count = 0;
	int notfound = false;
	while (m_tree[x]->size() < 2)
	{
		count++; 
		if (count > 1000) 
		{
			notfound = true; 
			break;
		}
		x = dist(gen);
	}
	if (notfound)
	{
		mutate();
		return;
	}
	//cout << "exited1" << endl;
	int y = rand() % m_tree[x]->size();

	vector<int>* rep1 = new vector<int>();
	rep1->push_back(m_tree[x]->at(y));
	vector<int>* cor1 = new vector<int>();
	//->push_back(dist(gen));//m_tree[x]->at(rand() % m_tree[x]->size()));
	cor1->push_back(m_tree[x]->at(dist(gen) % m_tree[x]->size()));
	//is rep1 in cor1 or cor 1 in rep1?
	bool placed = false;
	while(!placed)
	{
		count++; if (count > 1000) cout << "damn3";
		while (cor1->at(0) == rep1->at(0))
		{
			count++; if (count > 1000) cout << "damn2";
			//cor1->at(0) = dist(gen);//m_tree[x]->at(rand() % m_tree[x]->size());
			cor1->at(0) = m_tree[x]->at(dist(gen) % m_tree[x]->size());
			//cout << cor1->at(0) << " " << rep1->at(0) << endl;
		}

		//cout << "exited2" << endl;
		vector<int>* temp = getNode(rep1->at(0));
		vector<int>* temp2 = getNode(cor1->at(0));

		

		if (find(temp->begin(), temp->end(), cor1->at(0)) == temp->end() && find(temp2->begin(), temp2->end(), rep1->at(0)) == temp2->end()) // if you find cor at rep then replace cor with its parent
		{
			//dump();
			//cout << "moved " << m_tree[x]->at(y) << " to "<<cor1->at(0)<<endl;
			temp2->push_back(m_tree[x]->at(y));

			m_tree[x]->erase(m_tree[x]->begin()+y);
			//dump();
			placed = true;
		}
		else
		{
			//cor1->at(0) = dist(gen);//m_tree[x]->at(rand() % m_tree[x]->size());
			cor1->at(0) = m_tree[x]->at(dist(gen) % m_tree[x]->size());
		}
	}

	//cout << "exited3" << endl;
	mutate();
	if (false)//!isConnected())
	{
		cout << "disconnected" << endl;
		system("pause");
	}
	delete rep1;
	delete cor1;

}



pair<TreeGene*, TreeGene*> makeChildren(TreeGene* p1, TreeGene* p2)
{
	//cout << "in make children" << endl;
	


	int size = p1->size(); 
	
	normal_distribution<double> dist(log(size), 5);
	//uniform_int_distribution<> dist(0,size);

	TreeGene* c1 = new TreeGene(size);
	TreeGene* c2 = new TreeGene(size);
	
	c1->mutrate = p1->mutrate;
	c2->mutrate = p2->mutrate;

	bool* inc1 = new bool[size];
	bool* inc2 = new bool[size];
	inc1[0] = true;
	inc2[0] = true;
	for (int i = 1; i < size; i++)
	{
		inc1[i] = false;
		inc2[i] = false;
	}
	
	double doub= dist(gen);
	int cross = (int)doub;
	while (cross < 0 || cross >= size)
	{
		doub = dist(gen);
		cross = (int)doub;

	}


	//int cross = rand() % size;
	//cout << "cross at " << cross << endl;

	//cout << "before cross" << endl;
	//Before Cross
	for (int i = 0; i < cross; i++)
	{

		vector<int>* tmp1 = p1->getNode(i);
		vector<int>* tmp2 = p2->getNode(i);

		for (int j = 0; j < tmp1->size(); j++)
		{
			if (inc1[tmp1->at(j)] == false)
			{
				c1->addNode(i, tmp1->at(j));
				inc1[tmp1->at(j)] = true;
			}

		}

		for (int j = 0; j < tmp2->size(); j++)
		{
			if (inc2[tmp2->at(j)] == false)
			{
				c2->addNode(i, tmp2->at(j));
				inc2[tmp2->at(j)] = true;
			}
		}
	}
	//cout << "at cross" << endl;
	//Cross
	for (int l = cross; l < size; l++)
	{
		vector<int>* tmp1 = p2->getNode(l);
		vector<int>* tmp2 = p1->getNode(l);

		for (int j = 0; j < tmp2->size(); j++)
		{
			if (inc1[tmp2->at(j)] == false)
			{
				c1->addNode(l, tmp2->at(j));
				inc1[tmp2->at(j)] = true;
			}
		}
		for (int j = 0; j < tmp1->size(); j++)
		{
			if (inc2[tmp1->at(j)] == false)
			{
				c2->addNode(l, tmp1->at(j));
				inc2[tmp1->at(j)] = true;
			}
		}
	}
	//After Cross
	//cout << "after cross" << endl;
	/*
	for (int i = cross+1; i < size; i++)
	{
		tmp1 = p2->getNode(i);
		tmp2 = p1->getNode(i);

		for (int j = 0; j < tmp2->size(); j++)
		{
			if (inc1[tmp2->at(j)] == false)
			{
				//if(find(c1->getTree()[tmp2->at(j)]->begin(), c1->getTree()[tmp2->at(j)]->end(),i) == c1->getTree()[tmp2->at(j)]->end())
				c1->addNode(i, tmp2->at(j));


				inc1[tmp2->at(j)] = true;
			}
		}
		for (int j = 0; j < tmp1->size(); j++)
		{
			if (inc2[tmp1->at(j)] == false)
			{
				//if (find(c2->getTree()[tmp1->at(j)]->begin(), c2->getTree()[tmp1->at(j)]->end(), i) == c2->getTree()[tmp1->at(j)]->end())
				c2->addNode(i, tmp1->at(j));
				inc2[tmp1->at(j)] = true;
			}
		}
	}
	*/
	//cout << "before repair" << endl;
	//Repair
	vector<int>* rep1 = new vector<int>();
	vector<int>* rep2 = new vector<int>();

	for (int i = 0; i < size; i++)
	{
		if (inc1[i] == false)
		{
			rep1->push_back(i);
		}
		if (inc2[i] == false)
		{
			rep2->push_back(i);
		}
	}
	
	//cout << "before find" << endl;
	vector<int>* cor1 = p1->findintree(rep1);
	vector<int>* cor2 = p2->findintree(rep2);
	//cout << "after find" << endl;

	for (int i = 0; i < cor1->size(); i++)
	{
		
		vector<int>* temp = c1->getNode(rep1->at(i));
		if (find(temp->begin(), temp->end(), cor1->at(i)) != temp->end()) // if you find cor at rep then replace cor with its parent
		{
			cor1->at(i) = p2->findintree(cor1->at(i));
		}
	}
	//cout << "after first loop" << endl;
	for (int i = 0; i < cor2->size(); i++)
	{

		vector<int>* temp = c2->getNode(rep2->at(i));
		if (find(temp->begin(), temp->end(), cor2->at(i)) != temp->end()) // if you find cor at rep then replace cor with its parent
		{
			cor2->at(i) = p1->findintree(cor2->at(i));
		}
	}
	//cout << "after second loop" << endl;

	for (int i = 0; i < cor1->size(); i++)
	{
		//cout << cor1->at(i) << "  " << rep1->at(i) << endl;
		c1->addNode(cor1->at(i), rep1->at(i));
		
	}
	//cout << "inbetween" << endl;
//
	for (int i = 0; i < cor2->size(); i++)
	{
		//cout << cor2->at(i) << "  " << rep2->at(i) << endl;
		c2->addNode(cor2->at(i), rep2->at(i));

	}
	//cout << "after other loops" << endl;

	/*
	for (int i = 0; i < size; i++)
	{
		if (inc1[i] == false)
		{
			c1->addNode(p2->findintree(i), i);
		}
		if (inc2[i] == false)
		{
			c2->addNode(p1->findintree(i), i);
		}
	}
	*/
	pair<TreeGene*, TreeGene*> p = make_pair(c1, c2);
	

	//cout << "end make children" << endl;
	return p;
}

void TreeGene::shuffleSubTree()
{
	//cout << "to shuffle" << endl;
	uniform_int_distribution<> dist(0, m_n-1);
	int root = dist(gen);
	while(m_tree[root]->size()<1) root = dist(gen);
	int s = 0;
	vector<int> nodes;
	//selected subtree

	queue<int> q;
	q.push(root);
	while (!q.empty())
	{
		int c = q.front();
		q.pop();
		s++;
		nodes.push_back(c);
		for (int i = 0; i < m_tree[c]->size(); i++)
		{
			q.push(m_tree[c]->at(i));
		}
	}
	TreeGene* subtree = new TreeGene(s);
	subtree->randomizeTree();
	//subtree->dump();
	
	//dump();
	for (int i = 0; i < subtree->size(); i++)
	{
		m_tree[nodes[i]]->clear();
		for (int j = 0; j < subtree->getNode(i)->size(); j++)
		{
			m_tree[nodes[i]]->push_back(nodes[subtree->getNode(i)->at(j)]);
		}
	}

	//cout << "shuffled" << endl;
	//dump();
	delete subtree;
}

string stringTree( TreeGene * tree, Point * points, int v)
{
	string f = "";

	f+=(to_string(v))+"\r\n";
	f+=(to_string(tree->fitness))+"\r\n";
	f += to_string(tree->generation) + "\r\n";
	if (!THREEDEE)
		for (int i = 0; i < v; i++)
		{
			f+=(to_string(i) + " " + to_string(points[i].x) + " " + to_string(points[i].y) + " " + to_string(points[i].rad)) + "\r\n";
		}
	else {
		for (int i = 0; i < v; i++)
		{
			f+=(to_string(i) + " " + to_string(points[i].x) + " " + to_string(points[i].y) + " " + to_string(points[i].z) + " " + to_string(points[i].rad)) + "\r\n";
		}
	}
	int connections = 0;
	for (int i = 0; i < v; i++)
	{
		for (int j = 0; j < tree->getNode(i)->size(); j++)
		{
			connections++;
		}
	}
	f+=(to_string(connections)) + "\r\n";
	for (int i = 0; i < v; i++)
	{
		for (int j = 0; j < tree->getNode(i)->size(); j++)
		{
			f+=(to_string(i) + "-" + to_string(tree->getNode(i)->at(j))) + "\r\n";
		}
	}
	///////////////////////////////////////////////////////////
	return f;
}


void printTree(string path, TreeGene * tree, Point * points, int v)
{
	FileHandler f(path);

	f.print(to_string(v));
	f.print(to_string(tree->fitness));
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
		for (int j = 0; j < tree->getNode(i)->size(); j++)
		{
			connections++;
		}
	}
	f.print(to_string(connections));
	for (int i = 0; i < v; i++)
	{
		for (int j = 0; j < tree->getNode(i)->size(); j++)
		{
			f.print(to_string(i) + "-" + to_string(tree->getNode(i)->at(j)));
		}
	}
	///////////////////////////////////////////////////////////

}


double getFitness(TreeGene* g, Point * points, int NofPoints)
{
	double cost = 0;
	int C1 = 1;
	int C2 = 7;
	int C3 = 1;

	float totallength = totalPathLength(g->getTree(), points, NofPoints);
	float meanpathlength = meanPathLength(g->getTree(), points, NofPoints);
	float totalnetworkimpedence = getNetworkImpedence(g->getTree(), g->getRadGene(), points, NofPoints);
	
	cost = C1 * totallength + C2* meanpathlength + C3 * totalnetworkimpedence;
	
	g->fitness = 1/cost;

	return 1/cost;
}

int gettheLeaks()
{
	return theleaks;
}
