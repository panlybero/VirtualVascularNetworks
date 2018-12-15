#include "AugmentedGenome.h"
#include<iostream>
int GlobalInov = 0;

NodeGene::NodeGene(Point * in, Point* here, vector<Point*>* out) : m_in(in),m_this(here),m_out(out)
{
	
}

NodeGene::~NodeGene()
{
	delete m_out;

}

void NodeGene::changeIn(Point * in)
{
	m_in = in;
}

Point * NodeGene::getIn()
{
	return m_in;
}

Point * NodeGene::getThis()
{
	return m_this;
}

vector<Point*>* NodeGene::getOut()
{
	return m_out;
}

ConnectGene::ConnectGene(NodeGene* up, NodeGene* down, bool active, int inovNumber) : m_up(up), m_down(down),m_active(active), m_inovNumber(inovNumber)
{

}

ConnectGene::~ConnectGene()
{
	delete m_up;
	delete m_down;

}

NodeGene * ConnectGene::getUp()
{
	return m_up;
}

NodeGene * ConnectGene::getDown()
{
	return m_down;
}

bool ConnectGene::isActive()
{
	return m_active;
}

int ConnectGene::getinovNumber()
{
	return m_inovNumber;
}

void ConnectGene::setUp(NodeGene * new_up) 
{
	m_up = new_up;
}

void ConnectGene::setDown(NodeGene * new_down)
{
	m_down = new_down;
}

void ConnectGene::setActive(bool active)
{
	m_active = active;
}

void ConnectGene::setInovNumber(int n)
{
	m_inovNumber = n;
}

AugmentedGenome::AugmentedGenome()
{

}

void AugmentedGenome::mutSwapConnection()
{
	double ran = randomDouble(0, 1);
	if (ran <= MUTRATE / 100)
	{
		
		double pos1 = 1;
		double pos2 = 1; //rand() % genes.size();

		if (genes.size() != 0)
		{
			pos1 = rand() % genes.size();
			pos2 = rand() % genes.size();
		}
		ConnectGene* cg1 = new ConnectGene(genes.at(pos1)->getUp(), genes.at(pos2)->getDown(), genes.at(pos1)->isActive(), genes.at(pos1)->getinovNumber());
		ConnectGene* cg2 = new ConnectGene(genes.at(pos2)->getUp(), genes.at(pos1)->getDown(), genes.at(pos2)->isActive(), genes.at(pos2)->getinovNumber());
		ConnectGene* d1 = genes.at(pos1);
		ConnectGene* d2 = genes.at(pos2);
		genes[pos1] = cg1;
		genes[pos2] = cg2;
		// d1;
		//delete d2;


	}
}
