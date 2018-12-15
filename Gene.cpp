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
	m_mutrate = 0.8;
}

Gene::Gene(int geneLength)
{
	data = new float[geneLength];
	m_geneLength = geneLength;
	leaks++;
	m_mutrate = 0.8;
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
int heartSize = 1000;
Gene* generateRandomRadGene(int size)
{
	float* f = new float[size];
	Gene* g = new Gene(f, size);
	g->setMutrate(0.15);//randomDouble(0, 1));

	for (int i = 0; i < size; i++)
	{
		g->setPos(randomDouble(1, heartSize), i);
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