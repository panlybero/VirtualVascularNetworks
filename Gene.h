#pragma once


class Gene
{
	
public:
	Gene(float* g, int geneLength);
	Gene(int geneLength);
	~Gene();
	
	

	void setPos(float d, int pos);
	float getPos(int pos);
	int geneLength();
	float getMutrate();
	void setMutrate(float m);
	void dump();

private:
	float* data;
	float m_mutrate;
	int m_geneLength;
};

void copyGene(Gene* a, Gene* b);
int getLeaks();

Gene* generateRandomGene(int size);
Gene* generateRandomRadGene(int size);