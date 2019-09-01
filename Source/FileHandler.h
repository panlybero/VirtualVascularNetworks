#ifndef FileHand
#define FileHand


#include<cstdio>
#include<string>
#include<fstream>
#include<vector>
#include"Point.h"
#include"Gene.h"
using namespace std;

class FileHandler
{
public:
	FileHandler(string path);
	FileHandler();
	~FileHandler();
	void print(string text);
	Gene * readGene(string path);
	Point* readPoints(string path, int &size);
	Point* read3DPoints(string path, int &size);

private:

	fstream m_file;
	bool m_test = false;

};

#endif // !FileHand
