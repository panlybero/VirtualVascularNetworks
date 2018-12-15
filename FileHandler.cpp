#include"FileHandler.h"
#include<iostream>
FileHandler::FileHandler(string path)
{

	m_file.open(path, std::fstream::binary | std::fstream::in | std::fstream::out | std::fstream::trunc);
	


}

FileHandler::FileHandler()
{
}

FileHandler::~FileHandler()
{
	m_file.close();
}

void FileHandler::print(string text)
{
	m_file << text << "\r\n";
}

Point * FileHandler::readPoints(string path, int &size)
{
	fstream f;
	f.open(path);
	char* s = new char[150];

	f.getline(s, 150);


	string st = s;

	cout << st << endl;
	system("pause");
	int n = stoi(st);
	size = n;

	Point* points = new Point[size];

	for (int i = 0; i < n; i++)
	{
		Point p;
		f.getline(s, 100, '\n');
		st = s;
		string g = st.substr(2, 8);
		char* pEnd;
		double d = strtod(g.c_str(),&pEnd);
		p.x = d;
		g = st.substr(11, 8);
		d = strtod(g.c_str(), &pEnd);
		p.y = d;
		points[i] = p;
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
