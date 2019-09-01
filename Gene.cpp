

#include"Gene.h"
#include<iostream>
#include"util.h"

using namespace std;
int leaks = 0;


Gene::Gene(float * g, int geneLength)
{
	data = g;
	m_geneLength = geneLength;
	leaks++;
	m_mutrate = 0.15;
}

Gene::Gene(int geneLength)
{
	data = new float[geneLength];
	m_geneLength = geneLength;
	leaks++;
	m_mutrate = 0.15;
}

Gene::~Gene()
{
	//cout << "deleting" << endl;
	leaks--;
	delete[] data;
}
int getLeaks()
{
	return leaks;
}
void copyGene(Gene * a, Gene * b)
{
	if (a->geneLength() != b->geneLength()) { cerr << "unequal gene lengths" << endl; return; }
	for (int i = 0; i < b->geneLength(); i++)
	{
		a->setPos(b->getPos(i), i);
	}
	a->setMutrate(b->getMutrate());
}

void Gene::setPos(float d, int pos)
{
	if (pos >= 0 && pos < m_geneLength)
	{
		data[pos] = d;
	}
	else {
		cerr << pos<<" invalid gene setPos" << endl;
		system("pause");
	}
}

float Gene::getPos(int pos)
{
	if (pos >= 0 && pos < m_geneLength)
	{
		return data[pos];
	}
	cerr << pos << " invalid gene getPos" << endl;
	cerr << "geneLength " << m_geneLength<<endl;
	system("pause");
	return -1.0;
}

int Gene::geneLength()
{
	return m_geneLength;
}

float Gene::getMutrate()
{
	return m_mutrate;
}

void Gene::setMutrate(float m)
{
	m_mutrate = m;
}

void Gene::dump()
{
	for (int i = 0; i < m_geneLength; i++)
	{
		cout << data[i] << " ";
	}
	cout << endl;
}

Gene* generateRandomRadGene(int size)
{
	float* f = new float[size];
	Gene* g = new Gene(f, size);
	g->setMutrate(0.15);//randomDouble(0, 1));


	int heartSize = 10;

	for (int i = 0; i < size; i++)
	{
		g->setPos(randomDouble(1, heartSize), i);
		//g->setPos(heartSize, i);
	}


	return g;
}
Gene* generateRandomGeneBinary(int size)
{
	float* f = new float[size];
	Gene* g = new Gene(f, size);
	g->setMutrate(0.15);//randomDouble(0, 1));

	int n = (1 + sqrt(1 + 8 * size)) / 2;

	int zers = n - 1;
	/*
	for (int i = 0; i < size; i++)
	{
		if (randomDouble(0, 1) > 0.5)
		{
			g->setPos(1, i);
		}
		else if(zers>=0)
		{
			g->setPos(0, i);
			zers--;
		}

	}
	*/
	for (int i = 0; i < size; i++)
	{
		g->setPos(1, i);
	}
	//cout << size << endl;
	while (zers >= 0)
	{
		int p = randomInt(0, size-1);
		while (g->getPos(p) == 0)
		{
			p = randomInt(0, size - 1);

		}
		g->setPos(0, p);
		zers--;
	}

	

	return g;
}
Gene* generateRandomGene(int size)
{
	float* gene = new float[size];
	Gene* g = new Gene(gene, size);
	g->setMutrate(0.15);//randomDouble(0, 1));

	for (int i = 0; i < size; i++)
	{
		g->setPos(randomDouble(0, 2), i);

	}

	return g;
}