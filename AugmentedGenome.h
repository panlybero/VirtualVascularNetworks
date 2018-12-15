#pragma once
#include<vector>
#include"Point.h"
#include"util.h"

#define MUTRATE 100

using namespace std;

class ConnectGene;
class NodeGene;

class AugmentedGenome
{

public:
	AugmentedGenome();
	vector<ConnectGene*> genes;
	void mutSwapConnection(); // must be done carefully to avoid cycles
	//void mutAddPoint();




};



class NodeGene
{
public: 
	NodeGene(Point* in, Point* here, vector<Point*>* out);
	~NodeGene();
	void changeIn(Point* in);
	Point* getIn();
	Point* getThis();
	vector<Point*>* getOut();

private:
	Point* m_in;
	Point* m_this;
	vector<Point*>* m_out;



};

class ConnectGene
{
public:
	ConnectGene(NodeGene* up , NodeGene* down, bool active, int inovNumber);
	~ConnectGene();
	NodeGene* getUp();
	NodeGene* getDown();
	bool isActive();
	int getinovNumber();

	void setUp(NodeGene* new_up);
	void setDown(NodeGene* new_down);
	void setActive(bool active);
	void setInovNumber(int n);

private:
	NodeGene* m_up;
	NodeGene* m_down;
	bool m_active;
	int m_inovNumber;
};
