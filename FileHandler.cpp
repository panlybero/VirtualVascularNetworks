#include"FileHandler.h"
#include<iostream>

FileHandler::FileHandler(string path)
{

	m_file.open(path, std::fstream::binary | std::fstream::in | std::fstream::out | std::fstream::trunc);
	

	m_test = true;
}

FileHandler::FileHandler()
{
}

FileHandler::~FileHandler()
{
	if(m_test)
		m_file.close();

}

void FileHandler::print(string text)
{
	m_file << text << "\r\n";
}

Gene* FileHandler::readGene(string path)
{
	fstream f;
	f.open(path);
	char* s = new char[1000000];

	f.getline(s, 150);
	int n = stoi(s);

	f.getline(s, 1000000);

	//cout << s << endl;
	vector<int> vals;
	char * pch = NULL;
	char* next = NULL;
	char* pEnd;
	char seps[] = " ,\t\n";

	pch = strtok_s(s, seps, &next);
	vals.push_back((int)strtod(pch, &pEnd));

	while (pch != NULL)
	{
		//printf("%s\n", pch);
		pch = strtok_s(NULL, seps, &next);
		if (pch != NULL)
			vals.push_back((int)strtod(pch, &pEnd));

	}

	for (int i = 0; i < vals.size(); i++)
	{
		//cout << vals[i] << endl;
	}

	f.close();

	Gene* seed_gene = new Gene(n*(n-1)/2);
	for (int i = 0; i < vals.size(); i++) {
		seed_gene->setPos(vals[i], i);
	}

	delete[] s;
	return seed_gene;
}

Point * FileHandler::readPoints(string path, int &size)
{
	fstream f;
	f.open(path);
	char* s = new char[150];

	f.getline(s, 150);


	string st = s;

	int n = stoi(st);
	size = n;
	Point* points = new Point[n];
	vector<double> vals;
	for (int i = 0; i < n; i++)
	{
		Point p;
		f.getline(s, 100, '\n');
		
		st = s;
		char * pch = NULL;
		char* next = NULL;
		char* pEnd;
		char seps[] = " ,\t\n";
		
		pch = strtok_s(s, seps, &next);
		vals.push_back(strtod(pch, &pEnd));
	
		while (pch != NULL)
		{
			//printf("%s\n", pch);
			pch = strtok_s(NULL, seps,&next);
			if(pch!=NULL)
				vals.push_back(strtod(pch, &pEnd));
			
		}
		
		points[i].x = vals[1];
		points[i].y = vals[2];
		points[i].rad = 0;
		


		vals.clear();
	}
	f.close();
	
	return points;
}

Point * FileHandler::read3DPoints(string path, int & size)
{
	fstream f;
	f.open(path);
	char* s = new char[150];

	f.getline(s, 150);

	string st = s;


	int n = stoi(st);
	size = n;

	Point* points = new Point[size];
	// 1 234 456 999
	for (int i = 0; i < n; i++)
	{
		Point p;
		f.getline(s, 100, '\n');
		st = s;
		int ind = 0;
		int ind2 = 0;
		while (st[ind] != ' ') ind++;
		ind++;

		ind2 = ind;
		while (st[ind2] != ' ')	ind2++;

		char* pEnd;
		double x = strtod(st.substr(ind, ind2-1).c_str(), &pEnd);
		ind = ++ind2;


		while (st[ind2] != ' ')	ind2++;

		double y = strtod(st.substr(ind, ind2-1).c_str(), &pEnd);
		ind = ++ind2;


		while (st[ind2] != ' ')	ind2++;
		double z = strtod(st.substr(ind, ind2 - 1).c_str(), &pEnd);

		p.x = x;
		p.y = y;
		p.z = z;
		p.rad = 0;
		points[i] = p;
	}

	return points;
}
