#ifndef FileHand
#define FileHand


#include<cstdio>
#include<string>
#include<fstream>
#include"Point.h"
using namespace std;

class FileHandler
{
public:
	FileHandler(string path);
	FileHandler();
	~FileHandler();
	void print(string text);
	Point* readPoints(string path, int &size);
	Point* read3DPoints(string path, int &size);

private:

	fstream m_file;

};

#endif // !FileHand
