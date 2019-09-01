//#include "BESSEL.H"
#include <algorithm>
#include <cmath>
#include <complex>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <stack>
#include <sstream>
#include <streambuf>
#include <string>

#include "kiss.h"
#include "lodepng.h"

#include "drawTree.h"
#include "generalUtilities.h"
#include "terseStructs.h"
#include "tuneRadii.h"

using namespace std;

time_t startTime;

/*string niceTime(double t, bool roundSeconds = false){
	if(t >= 0.0){
		stringstream nT;
		int m(60);
		int h(60*m);
		int d(24*h);
		int w(7*d);
		int wt((int)floor(t/w));
		t -= w*wt;
		int dt((int)floor(t/d));
		t -= d*dt;
		int ht((int)floor(t/h));
		t -= h*ht;
		int mt((int)floor(t/m));
		t -= m*mt;
		if(wt > 0)
			nT << wt << "w ";
		if(wt > 0 || dt > 0)
			nT << dt << "d ";
		if(wt > 0 || dt > 0 || ht > 0)
			nT << ht << "h ";
		if(wt > 0 || dt > 0 || ht > 0 || mt > 0)
			nT << mt << "m ";
		if(roundSeconds)
			nT << ceil(t) << "s";
		else
			nT << t << "s";
		return nT.str();
	}else return "-(" + niceTime(-t) + ")";
}*/

string notIndeed(bool x){
	if(x)
		return "INDEED";
	return "NOT";
}

template<class T>
T min(const T &x, const T &y){
	if (x < y)
		return x;
	return y;
}
/*
template <class T>
string makeString(T thing){
	stringstream strstr;
	strstr << thing;
	return strstr.str();
}
*/
template <class T>
string makeVectorString(vector<T> thing){
	stringstream strstr;
	strstr << "(";
	if(thing.size() > 0){
		for(unsigned int i(0); i < thing.size() - 1; i++)
			strstr << thing[i] << ",";
		strstr << thing.back();
	}
	strstr << ")";
	return strstr.str();
}

template <class T>
string makeArrayString(T *thing, int size){
	stringstream strstr;
	strstr << "(";
	if(size > 0){
		for(int i(0); i < size - 1; i++)
			strstr << thing[i] << ",";
		strstr << thing[size - 1];
	}
	strstr << ")";
	return strstr.str();
}

/*
template <class T>
T maximum(vector<T> x){
	if(x.size() == 0)
		return T(-1000000);
	T m(x[0]);
	for(unsigned int i(1); i < x.size(); i++){
		if(m < x[i])
			m = x[i];
	}
	return m;
}

template <class T>
T minimum(vector<T> x){
	if(x.size() == 0)
		return T(1000000);
	T m(x[0]);
	for(unsigned int i(1); i < x.size(); i++){
		if(m > x[i])
			m = x[i];
	}
	return m;
}

template <class T>
double mean(vector<T> x){
	double sum(0.0);
	for(unsigned int i(0); i < x.size(); i++)
		sum += x[i];
	return sum/x.size();
}

template <class T>
double variance(vector<T> x){
	if(x.size() < 1)
		return -1.0;
	double sum(0.0);
	double sumSq(0.0);
	for(unsigned int i(0); i < x.size(); i++){
		sum += x[i];
		sumSq += x[i]*x[i];
	}
	sum /= x.size();
	return sumSq/x.size() - sum*sum;
}

template <class T>
double stDev(vector<T> x){
	if(x.size() < 2)
		return -1.0;
	double sum(0.0);
	double m(mean(x));
	for(unsigned int i(0); i < x.size(); i++){
		double term(x[i] - m);
		sum += term*term;
	}
	return sqrt(sum/(x.size() - 1));
}
*/

// I named it bifurcation, but it can really have more than two children; I'll name something Branch Point for contexts when higher branching order is likely
struct Bifurcation2D{
	//~Bifurcation2D(){ delete[] position; }
	Bifurcation2D(){ parentIndex = -1; position[0] = position[1] = 0.0; }
	//Bifurcation2D(Bifurcation2D &toBeCopied){
	//	parentIndex = toBeCopied.parentIndex;
	//	position[0] = toBeCopied.position[0];
	//	position[1] = toBeCopied.position[1];
	//	childIndex = toBeCopied.childIndex;
	//}
	Bifurcation2D(const Bifurcation2D &toBeCopied){
		parentIndex = toBeCopied.parentIndex;
		position[0] = toBeCopied.position[0];
		position[1] = toBeCopied.position[1];
		childIndex = toBeCopied.childIndex;
	}
	Bifurcation2D& operator=(const Bifurcation2D &toBeAssigned){
		parentIndex = toBeAssigned.parentIndex;
		position[0] = toBeAssigned.position[0];
		position[1] = toBeAssigned.position[1];
		childIndex = toBeAssigned.childIndex;
		return *this;
	}
	string tostring(){
		string str("pos = (" + makeString(position[0]) + ", " + makeString(position[1]) + "); par = " + makeString(parentIndex));
		if(childIndex.size() == 1)
			str += "; 1 child {";
		else
			str += "; " + makeString(childIndex.size()) + " children {";
		for(unsigned int c(0); c < childIndex.size(); c++)
			str += " " + makeString(childIndex[c]) + " ";
		return str + "}";
	}
	double position[2];
	int parentIndex;
	vector<int> childIndex;
};

// same as Bifurcation2D, but with a different name
struct BranchPoint2D{
	BranchPoint2D(){ parentIndex = -1; position[0] = position[1] = 0.0; }
	//Bifurcation2D(Bifurcation2D &toBeCopied){
	//	parentIndex = toBeCopied.parentIndex;
	//	position[0] = toBeCopied.position[0];
	//	position[1] = toBeCopied.position[1];
	//	childIndex = toBeCopied.childIndex;
	//}
	BranchPoint2D(const BranchPoint2D &toBeCopied){
		parentIndex = toBeCopied.parentIndex;
		position[0] = toBeCopied.position[0];
		position[1] = toBeCopied.position[1];
		childIndex = toBeCopied.childIndex;
	}
	BranchPoint2D& operator=(const BranchPoint2D &toBeAssigned){
		parentIndex = toBeAssigned.parentIndex;
		position[0] = toBeAssigned.position[0];
		position[1] = toBeAssigned.position[1];
		childIndex = toBeAssigned.childIndex;
		return *this;
	}
	string tostring(){
		string str("pos = (" + makeString(position[0]) + ", " + makeString(position[1]) + "); par = " + makeString(parentIndex));
		if(childIndex.size() == 1)
			str += "; 1 child {";
		else
			str += "; " + makeString(childIndex.size()) + " children {";
		for(unsigned int c(0); c < childIndex.size(); c++)
			str += " " + makeString(childIndex[c]) + " ";
		return str + "}";
	}
	double position[2];
	int parentIndex;
	vector<int> childIndex;
};

struct Bifurcation3D{
	Bifurcation3D(){ parentIndex = -1; position[0] = position[1] = position[2] = 0.0; }
	Bifurcation3D(const Bifurcation3D &toBeCopied){
		parentIndex = toBeCopied.parentIndex;
		position[0] = toBeCopied.position[0];
		position[1] = toBeCopied.position[1];
		position[2] = toBeCopied.position[2];
		childIndex = toBeCopied.childIndex;
	}
	Bifurcation3D& operator=(const Bifurcation3D &toBeAssigned){
		parentIndex = toBeAssigned.parentIndex;
		position[0] = toBeAssigned.position[0];
		position[1] = toBeAssigned.position[1];
		position[2] = toBeAssigned.position[2];
		childIndex = toBeAssigned.childIndex;
		return *this;
	}
	string tostring(){
		string str("pos = (" + makeString(position[0]) + ", " + makeString(position[1]) + ", " + makeString(position[2]) + "); par = " + makeString(parentIndex));
		if(childIndex.size() == 1)
			str += "; 1 child {";
		else
			str += "; " + makeString(childIndex.size()) + " children {";
		for(unsigned int c(0); c < childIndex.size(); c++)
			str += " " + makeString(childIndex[c]) + " ";
		return str + "}";
	}
	double position[3];
	int parentIndex;
	vector<int> childIndex;
};

struct BranchPoint3D{
	BranchPoint3D(){ parentIndex = -1; position[0] = position[1] = position[2] = 0.0; }
	BranchPoint3D(const BranchPoint3D &toBeCopied){
		parentIndex = toBeCopied.parentIndex;
		position[0] = toBeCopied.position[0];
		position[1] = toBeCopied.position[1];
		position[2] = toBeCopied.position[2];
		childIndex = toBeCopied.childIndex;
	}
	BranchPoint3D& operator=(const BranchPoint3D &toBeAssigned){
		parentIndex = toBeAssigned.parentIndex;
		position[0] = toBeAssigned.position[0];
		position[1] = toBeAssigned.position[1];
		position[2] = toBeAssigned.position[2];
		childIndex = toBeAssigned.childIndex;
		return *this;
	}
	string tostring(){
		string str("pos = (" + makeString(position[0]) + ", " + makeString(position[1]) + ", " + makeString(position[2]) + "); par = " + makeString(parentIndex));
		if(childIndex.size() == 1)
			str += "; 1 child {";
		else
			str += "; " + makeString(childIndex.size()) + " children {";
		for(unsigned int c(0); c < childIndex.size(); c++)
			str += " " + makeString(childIndex[c]) + " ";
		return str + "}";
	}
	double position[3];
	int parentIndex;
	vector<int> childIndex;
};

template <class T>
struct Pair{
	Pair(T a, T b){ x = a; y = b; }
	T x, y;
	bool sameOrder(const Pair &p){
		return p.x == x && p.y == y;
	}
	bool sameUnordered(const Pair &p){
		return (p.x == x && p.y == y) || (p.x == y && p.y == x);
	}
};

void testCopyConstructor(){
	Bifurcation2D a;
	a.position[0] = 1.0;
	a.position[0] = 2.0;
	a.parentIndex = 3;
	a.childIndex.push_back(4);
	a.childIndex.push_back(5);
	Bifurcation2D b(a);
	Bifurcation2D c = a;
	cout << "\n testCopyConstructor()" << endl;
	cout << "a = " << a.tostring() << endl;
	cout << "b = " << b.tostring() << endl;
	cout << "c = " << c.tostring() << endl;
	vector<Bifurcation2D> primary;
	primary.push_back(a);
	primary.push_back(b);
	primary.push_back(c);
	vector<Bifurcation2D> secondary(primary);
	vector<Bifurcation2D> tertiary = primary;
	cout << "\n primary:" << endl;
	for(unsigned int i(0); i < primary.size(); i++)
		cout << "i = " << i << " : \t" << primary[i].tostring() << endl;
	cout << "\n secondary:" << endl;
	for(unsigned int i(0); i < secondary.size(); i++)
		cout << "i = " << i << " : \t" << secondary[i].tostring() << endl;
	cout << "\n tertiary:" << endl;
	for(unsigned int i(0); i < tertiary.size(); i++)
		cout << "i = " << i << " : \t" << tertiary[i].tostring() << endl;
}

double separation2D(const Bifurcation2D &a, const Bifurcation2D &b){
	double dx(a.position[0] - b.position[0]), dy(a.position[1] - b.position[1]);
	return sqrt(dx*dx + dy*dy);
}

double separation2D(const BranchPoint2D &a, const BranchPoint2D &b){
	double dx(a.position[0] - b.position[0]), dy(a.position[1] - b.position[1]);
	return sqrt(dx*dx + dy*dy);
}

double separation2D(double x1, double y1, double x2, double y2){
	double dx(x1 - x2), dy(y1 - y2);
	return sqrt(dx*dx + dy*dy);
}

double separation3D(const Bifurcation3D &a, const Bifurcation3D &b){
	double dx(a.position[0] - b.position[0]), dy(a.position[1] - b.position[1]), dz(a.position[2] - b.position[2]);
	return sqrt(dx*dx + dy*dy + dz*dz);
}

double separation3D(const BranchPoint3D &a, const BranchPoint3D &b){
	double dx(a.position[0] - b.position[0]), dy(a.position[1] - b.position[1]), dz(a.position[2] - b.position[2]);
	return sqrt(dx*dx + dy*dy + dz*dz);
}

double separation3D(double x1, double y1, double z1, double x2, double y2, double z2){
	double dx(x1 - x2), dy(y1 - y2), dz(z1 - z2);
	return sqrt(dx*dx + dy*dy + dz*dz);
}

void addChild(vector<Bifurcation2D> &b, int parentIndex, double positionX, double positionY){
	b.push_back(Bifurcation2D());
	b.back().parentIndex = parentIndex;
	b.back().position[0] = positionX;
	b.back().position[1] = positionY;
	b[parentIndex].childIndex.push_back(b.size() - 1);
}

void addParent(vector<Bifurcation2D> &b, int child1Index, int child2Index, double positionX, double positionY){
	if(child1Index == child2Index)
		return;
	b.push_back(Bifurcation2D());
	b.back().position[0] = positionX;
	b.back().position[1] = positionY;
	b.back().childIndex.push_back(child1Index);
	b.back().childIndex.push_back(child2Index);
	b[child1Index].parentIndex = b[child2Index].parentIndex = b.size() - 1;
}

void addParent(vector<Bifurcation3D> &b, int child1Index, int child2Index, double positionX, double positionY, double positionZ){
	if(child1Index == child2Index)
		return;
	b.push_back(Bifurcation3D());
	b.back().position[0] = positionX;
	b.back().position[1] = positionY;
	b.back().position[2] = positionZ;
	b.back().childIndex.push_back(child1Index);
	b.back().childIndex.push_back(child2Index);
	b[child1Index].parentIndex = b[child2Index].parentIndex = b.size() - 1;
}

void addParent(vector<Bifurcation2D> &b, int child1Index, int child2Index){
	if(child1Index == child2Index)
		return;
	double positionX((b[child1Index].position[0] + b[child2Index].position[0])/2.0),
		positionY((b[child1Index].position[1] + b[child2Index].position[1])/2.0);
	b.push_back(Bifurcation2D());
	b.back().position[0] = positionX;
	b.back().position[1] = positionY;
	b.back().childIndex.push_back(child1Index);
	b.back().childIndex.push_back(child2Index);
	b[child1Index].parentIndex = b[child2Index].parentIndex = b.size() - 1;
}

void removeBifurcation(vector<Bifurcation2D> &b, int index){
	for(unsigned int c(0); c < b[index].childIndex.size(); c++){
		b[b[index].childIndex[c]].parentIndex = -1;
	}
	if(b[index].parentIndex > -1){
		for(unsigned int c(0); c < b[b[index].parentIndex].childIndex.size(); c++){
			if(b[b[index].parentIndex].childIndex[c] == index){
				b[b[index].parentIndex].childIndex.erase(b[b[index].parentIndex].childIndex.begin() + c);
				c = b[b[index].parentIndex].childIndex.size();
			}
		}
	}
	b.erase(b.begin() + index);
}

void insertBifurcation(vector<Bifurcation2D> &b, int parentIndex, int childIndex1, int childIndex2, double positionX, double positionY){
	b.push_back(Bifurcation2D());
	b.back().parentIndex = parentIndex;
	for(vector<int>::iterator it(b[parentIndex].childIndex.begin()); it != b[parentIndex].childIndex.end(); it++){
		if(*it == childIndex1){
			it = b[parentIndex].childIndex.erase(it);
			if(it == b[parentIndex].childIndex.end())
				break;
		}
	}
	for(vector<int>::iterator it(b[parentIndex].childIndex.begin()); it != b[parentIndex].childIndex.end(); it++){
		if(*it == childIndex2){
			it = b[parentIndex].childIndex.erase(it);
			if(it == b[parentIndex].childIndex.end())
				break;
			
		}
	}
	b[parentIndex].childIndex.push_back(b.size() - 1);
	b.back().childIndex.push_back(childIndex1);
	b.back().childIndex.push_back(childIndex2);
	b[childIndex1].parentIndex = b[childIndex2].parentIndex = b.size() - 1;
	b.back().position[0] = positionX;
	b.back().position[1] = positionY;
}

void insertBifurcation(vector<Bifurcation2D> &b, int parentIndex, int childIndex1, int childIndex2){
	double positionX((b[parentIndex].position[0] + b[childIndex1].position[0] + b[childIndex2].position[0])/3.0),
		positionY((b[parentIndex].position[1] + b[childIndex1].position[1] + b[childIndex2].position[1])/3.0);
	b.push_back(Bifurcation2D());
	b.back().parentIndex = parentIndex;
	for(vector<int>::iterator it(b[parentIndex].childIndex.begin()); it != b[parentIndex].childIndex.end(); it++){
		if(*it == childIndex1){
			it = b[parentIndex].childIndex.erase(it);
			if(it == b[parentIndex].childIndex.end())
				break;
			//it = b[parentIndex].childIndex.end();
		}
	}
	for(vector<int>::iterator it(b[parentIndex].childIndex.begin()); it != b[parentIndex].childIndex.end(); it++){
		if(*it == childIndex2){
			it = b[parentIndex].childIndex.erase(it);
			if(it == b[parentIndex].childIndex.end())
				break;
			//it = b[parentIndex].childIndex.end();
		}
	}
	b[parentIndex].childIndex.push_back(b.size() - 1);
	b.back().childIndex.push_back(childIndex1);
	b.back().childIndex.push_back(childIndex2);
	b[childIndex1].parentIndex = b[childIndex2].parentIndex = b.size() - 1;
	b.back().position[0] = positionX;
	b.back().position[1] = positionY;
}

double minX(const vector<Bifurcation2D> &b){
	double x(b[0].position[0]);
	for(unsigned int i(1); i < b.size(); i++){
		if(x > b[i].position[0])
			x = b[i].position[0];
	}
	return x;
}

double minX(const vector<BranchPoint2D> &b){
	double x(b[0].position[0]);
	for(unsigned int i(1); i < b.size(); i++){
		if(x > b[i].position[0])
			x = b[i].position[0];
	}
	return x;
}

double maxX(const vector<Bifurcation2D> &b){
	double x(b[0].position[0]);
	for(unsigned int i(1); i < b.size(); i++){
		if(x < b[i].position[0])
			x = b[i].position[0];
	}
	return x;
}

double maxX(const vector<BranchPoint2D> &b){
	double x(b[0].position[0]);
	for(unsigned int i(1); i < b.size(); i++){
		if(x < b[i].position[0])
			x = b[i].position[0];
	}
	return x;
}

double minY(const vector<Bifurcation2D> &b){
	double y(b[0].position[1]);
	for(unsigned int i(1); i < b.size(); i++){
		if(y > b[i].position[1])
			y = b[i].position[1];
	}
	return y;
}

double minY(const vector<BranchPoint2D> &b){
	double y(b[0].position[1]);
	for(unsigned int i(1); i < b.size(); i++){
		if(y > b[i].position[1])
			y = b[i].position[1];
	}
	return y;
}

double maxY(const vector<Bifurcation2D> &b){
	double y(b[0].position[1]);
	for(unsigned int i(1); i < b.size(); i++){
		if(y < b[i].position[1])
			y = b[i].position[1];
	}
	return y;
}

double maxY(const vector<BranchPoint2D> &b){
	double y(b[0].position[1]);
	for(unsigned int i(1); i < b.size(); i++){
		if(y < b[i].position[1])
			y = b[i].position[1];
	}
	return y;
}

double minX(const vector<Bifurcation3D> &b){
	double x(b[0].position[0]);
	for(unsigned int i(1); i < b.size(); i++){
		if(x > b[i].position[0])
			x = b[i].position[0];
	}
	return x;
}

double minX(const vector<BranchPoint3D> &b){
	double x(b[0].position[0]);
	for(unsigned int i(1); i < b.size(); i++){
		if(x > b[i].position[0])
			x = b[i].position[0];
	}
	return x;
}

double maxX(const vector<Bifurcation3D> &b){
	double x(b[0].position[0]);
	for(unsigned int i(1); i < b.size(); i++){
		if(x < b[i].position[0])
			x = b[i].position[0];
	}
	return x;
}

double maxX(const vector<BranchPoint3D> &b){
	double x(b[0].position[0]);
	for(unsigned int i(1); i < b.size(); i++){
		if(x < b[i].position[0])
			x = b[i].position[0];
	}
	return x;
}

double minY(const vector<Bifurcation3D> &b){
	double y(b[0].position[1]);
	for(unsigned int i(1); i < b.size(); i++){
		if(y > b[i].position[1])
			y = b[i].position[1];
	}
	return y;
}

double minY(const vector<BranchPoint3D> &b){
	double y(b[0].position[1]);
	for(unsigned int i(1); i < b.size(); i++){
		if(y > b[i].position[1])
			y = b[i].position[1];
	}
	return y;
}

double maxY(const vector<Bifurcation3D> &b){
	double y(b[0].position[1]);
	for(unsigned int i(1); i < b.size(); i++){
		if(y < b[i].position[1])
			y = b[i].position[1];
	}
	return y;
}

double maxY(const vector<BranchPoint3D> &b){
	double y(b[0].position[1]);
	for(unsigned int i(1); i < b.size(); i++){
		if(y < b[i].position[1])
			y = b[i].position[1];
	}
	return y;
}

double minZ(const vector<BranchPoint3D> &b){
	double z(b[0].position[2]);
	for(unsigned int i(1); i < b.size(); i++){
		if(z > b[i].position[2])
			z = b[i].position[2];
	}
	return z;
}

double maxZ(const vector<BranchPoint3D> &b){
	double z(b[0].position[2]);
	for(unsigned int i(1); i < b.size(); i++){
		if(z < b[i].position[2])
			z = b[i].position[2];
	}
	return z;
}

inline unsigned int imageIndex(int x, int y, int width, int height){
	if(x < 0 || y < 0 || x >= width || y >= height)
		return 0;
	return 4*(width*(height - y - 1) + x);
}

void drawLine(vector<unsigned char> &image, int width, int height, int r, int g, int b, double startX, double startY, double endX, double endY, double radius = -1.0){
	//cout << "\n drawLine(): from (" << startX << ", " << startY << ") to (" << endX << ", " << endY << ")" << endl;
	if(radius < 0.0)
		radius = double(width)/300.0;
	double separation(separation2D(startX, startY, endX, endY));
	for(double len(0.0); len < separation; len += 0.49){
		double cx(startX + (endX - startX)*len/separation),
			cy(startY + (endY - startY)*len/separation);
		for(int x((int)floor(cx - radius)); x <= (int)ceil(cx + radius); x++){
			for(int y((int)floor(cy - radius)); y <= (int)ceil(cy + radius); y++){
				if(separation2D(x, y, cx, cy) < radius){
					int index(imageIndex(x, y, width, height));
					image[index] = r;
					image[index + 1] = g;
					image[index + 2] = b;
					image[index + 3] = 255;
				}
			}
		}

		//cout << "\n\t drawLine(): coloring (" << x << ", " << y << ")" << endl;
		
	}
}

void drawLine(vector<unsigned char> &image, int width, int height, int rStart, int rEnd, int gStart, int gEnd, int bStart, int bEnd, double startX, double startY, double endX, double endY, double radius = -1.0){
	//cout << "\n drawLine(): from (" << startX << ", " << startY << ") to (" << endX << ", " << endY << ")" << endl;
	if(radius < 0.0)
		radius = double(width)/300.0;
	double separation(separation2D(startX, startY, endX, endY));
	for(double len(0.0); len < separation; len += 0.49){
		double cx(startX + (endX - startX)*len/separation),
			cy(startY + (endY - startY)*len/separation);
		for(int x((int)floor(cx - radius)); x <= (int)ceil(cx + radius); x++){
			for(int y((int)floor(cy - radius)); y <= (int)ceil(cy + radius); y++){
				if(separation2D(x, y, cx, cy) < radius){
					int index(imageIndex(x, y, width, height));
					double frac(len/separation);
					int r(rStart + (int)floor(frac*(rEnd - rStart))),
						g(gStart + (int)floor(frac*(gEnd - gStart))),
						b(bStart + (int)floor(frac*(bEnd - bStart)));
					image[index] = r;
					image[index + 1] = g;
					image[index + 2] = b;
					image[index + 3] = 255;
				}
			}
		}
	}
}

void drawLineNoOverwriteLessGreen(vector<unsigned char> &image, int width, int height, int rStart, int rEnd, int gStart, int gEnd, int bStart, int bEnd, double startX, double startY, double endX, double endY, double radius = -1.0){
	//cout << "\n drawLine(): from (" << startX << ", " << startY << ") to (" << endX << ", " << endY << ")" << endl;
	if(radius < 0.0)
		radius = double(width)/300.0;
	double separation(separation2D(startX, startY, endX, endY));
	for(double len(0.0); len < separation; len += 0.49){
		double cx(startX + (endX - startX)*len/separation),
			cy(startY + (endY - startY)*len/separation);
		for(int x((int)floor(cx - radius)); x <= (int)ceil(cx + radius); x++){
			for(int y((int)floor(cy - radius)); y <= (int)ceil(cy + radius); y++){
				if(separation2D(x, y, cx, cy) < radius){
					int index(imageIndex(x, y, width, height));
					double frac(len/separation);
					int r(rStart + (int)floor(frac*(rEnd - rStart))),
						g(gStart + (int)floor(frac*(gEnd - gStart))),
						b(bStart + (int)floor(frac*(bEnd - bStart)));
					if(image[index + 3] == 0 || image[index + 1] > g){
						image[index] = r;
						image[index + 1] = g;
						image[index + 2] = b;
						image[index + 3] = 255;
					}
				}
			}
		}
	}
}

bool containedIn(int element, const vector<int> &vec){
	for(unsigned int i(0); i < vec.size(); i++){
		if(vec[i] == element)
			return true;
	}
	return false;
}

bool containedIn(int element, vector<vector<int> > &vec){
	for(unsigned int i(0); i < vec.size(); i++){
		if(containedIn(element, vec[i]))
			return true;
	}
	return false;
}

unsigned int maxLevel(vector<Bifurcation2D> &b, int index, vector<int> &passedIndices){
	unsigned int mostLevels(0);
	passedIndices.push_back(index);
	for(unsigned int c(0); c < b[index].childIndex.size(); c++){
		unsigned int levels(0);
		if(!containedIn(b[index].childIndex[c], passedIndices))
			levels = maxLevel(b, b[index].childIndex[c], passedIndices);
		else{
			cout << "\nWarning maxLevel(): reticulation found!" << endl;
			//abort();
		}
		if(mostLevels < levels)
			mostLevels = levels;
	}
	return mostLevels + 1;
}

unsigned int maxLevel(const vector<BranchPoint2D> &b, int index, vector<int> &passedIndices){
	unsigned int mostLevels(0);
	passedIndices.push_back(index);
	for(unsigned int c(0); c < b[index].childIndex.size(); c++){
		unsigned int levels(0);
		if(!containedIn(b[index].childIndex[c], passedIndices))
			levels = maxLevel(b, b[index].childIndex[c], passedIndices);
		else{
			cout << "\nWarning maxLevel(): reticulation found!" << endl;
			//abort();
		}
		if(mostLevels < levels)
			mostLevels = levels;
	}
	return mostLevels + 1;
}

unsigned int maxLevel(vector<Bifurcation2D> &b, int index){
	vector<int> passedIndices;
	return maxLevel(b, index, passedIndices);
}

unsigned int maxLevel(const vector<BranchPoint2D> &b, int index){
	vector<int> passedIndices;
	return maxLevel(b, index, passedIndices);
}

unsigned int maxLevel(vector<Bifurcation3D> &b, int index, vector<int> &passedIndices){
	unsigned int mostLevels(0);
	passedIndices.push_back(index);
	for(unsigned int c(0); c < b[index].childIndex.size(); c++){
		unsigned int levels(0);
		if(!containedIn(b[index].childIndex[c], passedIndices))
			levels = maxLevel(b, b[index].childIndex[c], passedIndices);
		else{
			cout << "\nWarning maxLevel(): reticulation found!" << endl;
			//abort();
		}
		if(mostLevels < levels)
			mostLevels = levels;
	}
	return mostLevels + 1;
}

unsigned int maxLevel(const vector<BranchPoint3D> &b, int index, vector<int> &passedIndices){
	unsigned int mostLevels(0);
	passedIndices.push_back(index);
	for(unsigned int c(0); c < b[index].childIndex.size(); c++){
		unsigned int levels(0);
		if(!containedIn(b[index].childIndex[c], passedIndices))
			levels = maxLevel(b, b[index].childIndex[c], passedIndices);
		else{
			cout << "\nWarning maxLevel(): reticulation found!" << endl;
			//abort();
		}
		if(mostLevels < levels)
			mostLevels = levels;
	}
	return mostLevels + 1;
}

unsigned int maxLevel(vector<Bifurcation3D> &b, int index){
	vector<int> passedIndices;
	return maxLevel(b, index, passedIndices);
}

unsigned int maxLevel(const vector<BranchPoint3D> &b, int index){
	vector<int> passedIndices;
	return maxLevel(b, index, passedIndices);
}

void drawCircle(vector<unsigned char> &image, int width, int height, int r, int g, int b, double centerX, double centerY, double radius, double lineRadius = -1.0){
	//cout << "\n drawLine(): from (" << startX << ", " << startY << ") to (" << endX << ", " << endY << ")" << endl;
	if(lineRadius < 0.0)
		lineRadius = double(width)/300.0;
	double dtheta(1.0/(4.1*cos(-1.0)*radius));
	for(double angle(0.0); angle < 2.0*acos(-1.0) + dtheta; angle += dtheta){
		double cx(centerX + radius*cos(angle)),
			cy(centerY + radius*sin(angle));
		for(int x((int)floor(cx - lineRadius)); x <= (int)ceil(cx + lineRadius); x++){
			for(int y((int)floor(cy - lineRadius)); y <= (int)ceil(cy + lineRadius); y++){
				double sep(separation2D(x, y, centerX, centerY));
				if(radius - lineRadius < sep && sep < radius + lineRadius){
					int index(imageIndex(x, y, width, height));
						image[index] = r;
						image[index + 1] = g;
						image[index + 2] = b;
						image[index + 3] = 255;
				}
			}
		}
	}
}

void drawTriangle(vector<unsigned char> &image, int width, int height, int r, int g, int b, double centerX, double centerY, double radius, double lineRadius = -1.0){
	if(lineRadius < 0.0)
		lineRadius = double(width)/300.0;
	double upperLeftX((int)floor(centerX - radius)), upperLeftY((int)floor(centerY + radius)),
		upperRightX((int)floor(centerX + radius)), upperRightY((int)floor(centerY + radius)),
		lowerPointX((int)floor(centerX)), lowerPointY((int)floor(centerY - radius));
	//cout << "\n drawLine(): from (" << startX << ", " << startY << ") to (" << endX << ", " << endY << ")" << endl;
	drawLine(image, width, height, r, g, b, upperLeftX, upperLeftY, upperRightX, upperRightY, lineRadius);
	drawLine(image, width, height, r, g, b, upperLeftX, upperLeftY, lowerPointX, lowerPointY, lineRadius);
	drawLine(image, width, height, r, g, b, upperRightX, upperRightY, lowerPointX, lowerPointY, lineRadius);
}

void drawSquare(vector<unsigned char> &image, int width, int height, int r, int g, int b, double centerX, double centerY, double radius, double lineRadius = -1.0){
	if(lineRadius < 0.0)
		lineRadius = double(width)/300.0;
	double upperLeftX((int)floor(centerX - radius)), upperLeftY((int)floor(centerY + radius)),
		upperRightX((int)floor(centerX + radius)), upperRightY((int)floor(centerY + radius)),
		lowerLeftX((int)floor(centerX - radius)), lowerLeftY((int)floor(centerY - radius)),
		lowerRightX((int)floor(centerX + radius)), lowerRightY((int)floor(centerY - radius));
	//cout << "\n drawLine(): from (" << startX << ", " << startY << ") to (" << endX << ", " << endY << ")" << endl;
	drawLine(image, width, height, r, g, b, upperLeftX, upperLeftY, upperRightX, upperRightY, lineRadius);
	drawLine(image, width, height, r, g, b, upperLeftX, upperLeftY, lowerLeftX, lowerLeftY, lineRadius);
	drawLine(image, width, height, r, g, b, lowerRightX, lowerRightY, lowerLeftX, lowerLeftY, lineRadius);
	drawLine(image, width, height, r, g, b, lowerRightX, lowerRightY, upperRightX, upperRightY, lineRadius);
}

void drawBifurcationTree(string fn, vector<Bifurcation2D> &b, int width, int height, int border){
	vector<unsigned char> image;
	image.resize(width*height*4);
	for(int i(0); i < width*height*4; i++)
		image[i] = 0;
	double offsetX(minX(b)), offsetY(minY(b)), scaleX((width - 2*border)/(maxX(b) - offsetX)), scaleY((height - 2*border)/(maxY(b) - offsetY));
	if(scaleY > scaleX)
		scaleY = scaleX;
	if(scaleX > scaleY)
		scaleX = scaleY;
	//cout << "\ndrawBifurcationTree(): offset = (" << offsetX << ", " << offsetY << "), scale = (" << scaleX << ", " << scaleY << ")" << endl;
	int red(255), green(55), blue(155),
		redTip(0), greenTip(255), blueTip(155),
		redNode(55), greenNode(155), blueNode(255),
		redHeart(255), greenHeart(55), blueHeart(55);
	for(unsigned int i(0); i < b.size(); i++){
		for(unsigned int c(0); c < b[i].childIndex.size(); c++){
			double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border);
			double startX(centerX),//scaleX*(b[i].position[0] - offsetX) + border),
				startY(centerY),//scaleY*(b[i].position[1] - offsetY) + border),
				endX(scaleX*(b[b[i].childIndex[c]].position[0] - offsetX) + border),
				endY(scaleY*(b[b[i].childIndex[c]].position[1] - offsetY) + border);
			drawLine(image, width, height, red, green, blue, startX, startY, endX, endY);
		}
	}
	for(unsigned int i(0); i < b.size(); i++){
		double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border),
			standardRadius(border/2);
		if(i == 0)
			drawCircle(image, width, height, redHeart, greenHeart, blueHeart, centerX, centerY, 2*standardRadius);
		else if(b[i].childIndex.size() < 1)
			drawCircle(image, width, height, redTip, greenTip, blueTip, centerX, centerY, standardRadius);
		else
			drawCircle(image, width, height, redNode, greenNode, blueNode, centerX, centerY, standardRadius/8 + (int)floor(1.5*maxLevel(b, i)));
	}
	unsigned error = lodepng::encode(fn.c_str(), image, width, height);
	if(error)
		cout << "\nError snapshotPNG(): error " << error << " from lodepng" << endl;
}

void drawBentBifurcationTree(string fn, vector<Bifurcation2D> &b, double bendX, double bendY, int width, int height, int border){
	vector<unsigned char> image;
	image.resize(width*height*4);
	for(int i(0); i < width*height*4; i++)
		image[i] = 0;
	double offsetX(minX(b)), offsetY(minY(b)), scaleX((width - 2*border)/(maxX(b) - offsetX)), scaleY((height - 2*border)/(maxY(b) - offsetY));
	if(scaleY > scaleX)
		scaleY = scaleX;
	if(scaleX > scaleY)
		scaleX = scaleY;
	drawLine(image, width, height, 0, 0, 0, 0, height - 1, width - 1, height - 1, 0.5);
	drawLine(image, width, height, 0, 0, 0, width - 1, height - 1, width - 1, 0, 0.5);
	drawLine(image, width, height, 0, 0, 0, width - 1, 0, (int)floor(scaleX*(bendX - offsetX)), 0, 0.5);
	drawLine(image, width, height, 0, 0, 0, (int)floor(scaleX*(bendX - offsetX)), 0, (int)floor(scaleX*(bendX - offsetX)), (int)floor(scaleY*(bendY - offsetY)), 0.5);
	drawLine(image, width, height, 0, 0, 0, (int)floor(scaleX*(bendX - offsetX)), (int)floor(scaleY*(bendY - offsetY)), 0, (int)floor(scaleY*(bendY - offsetY)), 0.5);
	drawLine(image, width, height, 0, 0, 0, 0, (int)floor(scaleY*(bendY - offsetY)), 0, height - 1, 0.5);
	//cout << "\ndrawBifurcationTree(): offset = (" << offsetX << ", " << offsetY << "), scale = (" << scaleX << ", " << scaleY << ")" << endl;
	int red(255), green(55), blue(155),
		redTip(0), greenTip(255), blueTip(155),
		redNode(55), greenNode(155), blueNode(255),
		redHeart(255), greenHeart(55), blueHeart(55);
	for(unsigned int i(0); i < b.size(); i++){
		for(unsigned int c(0); c < b[i].childIndex.size(); c++){
			double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border);
			double startX(centerX),//scaleX*(b[i].position[0] - offsetX) + border),
				startY(centerY),//scaleY*(b[i].position[1] - offsetY) + border),
				endX(scaleX*(b[b[i].childIndex[c]].position[0] - offsetX) + border),
				endY(scaleY*(b[b[i].childIndex[c]].position[1] - offsetY) + border);
			drawLine(image, width, height, red, green, blue, startX, startY, endX, endY);
		}
	}
	for(unsigned int i(0); i < b.size(); i++){
		double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border),
			standardRadius(border/2);
		if(i == 0)
			drawCircle(image, width, height, redHeart, greenHeart, blueHeart, centerX, centerY, 2*standardRadius);
		else if(b[i].childIndex.size() < 1)
			drawCircle(image, width, height, redTip, greenTip, blueTip, centerX, centerY, standardRadius);
		else
			drawCircle(image, width, height, redNode, greenNode, blueNode, centerX, centerY, standardRadius/8 + (int)floor(1.5*maxLevel(b, i)));
	}
	unsigned error = lodepng::encode(fn.c_str(), image, width, height);
	if(error)
		cout << "\nError snapshotPNG(): error " << error << " from lodepng" << endl;
}

void drawBifurcationTree(string fn, vector<Bifurcation3D> &b, int width, int height, int border){
	vector<unsigned char> image;
	image.resize(width*height*4);
	for(int i(0); i < width*height*4; i++)
		image[i] = 0;
	double offsetX(minX(b)), offsetY(minY(b)), scaleX((width - 2*border)/(maxX(b) - offsetX)), scaleY((height - 2*border)/(maxY(b) - offsetY));
	if(scaleY > scaleX)
		scaleY = scaleX;
	if(scaleX > scaleY)
		scaleX = scaleY;
	//cout << "\ndrawBifurcationTree(): offset = (" << offsetX << ", " << offsetY << "), scale = (" << scaleX << ", " << scaleY << ")" << endl;
	int red(255), green(55), blue(155),
		redTip(0), greenTip(255), blueTip(155),
		redNode(55), greenNode(155), blueNode(255),
		redHeart(255), greenHeart(55), blueHeart(55);
	for(unsigned int i(0); i < b.size(); i++){
		for(unsigned int c(0); c < b[i].childIndex.size(); c++){
			double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border);
			double startX(centerX),//scaleX*(b[i].position[0] - offsetX) + border),
				startY(centerY),//scaleY*(b[i].position[1] - offsetY) + border),
				endX(scaleX*(b[b[i].childIndex[c]].position[0] - offsetX) + border),
				endY(scaleY*(b[b[i].childIndex[c]].position[1] - offsetY) + border);
			drawLine(image, width, height, red, green, blue, startX, startY, endX, endY);
		}
	}
	for(unsigned int i(0); i < b.size(); i++){
		double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border),
			standardRadius(border/2);
		if(i == 0)
			drawCircle(image, width, height, redHeart, greenHeart, blueHeart, centerX, centerY, 2*standardRadius);
		else if(b[i].childIndex.size() < 1)
			drawCircle(image, width, height, redTip, greenTip, blueTip, centerX, centerY, standardRadius);
		else
			drawCircle(image, width, height, redNode, greenNode, blueNode, centerX, centerY, standardRadius/8 + (int)floor(1.5*maxLevel(b, i)));
	}
	unsigned error = lodepng::encode(fn.c_str(), image, width, height);
	if(error)
		cout << "\nError snapshotPNG(): error " << error << " from lodepng" << endl;
}

double rFrac(double frac){
	if(frac > 1.0)
		frac = 1.0;
	if(frac < 0.2)
		return 1.0;
	if(frac < 0.4)
		return 2.0 - frac/0.2;//1.0 - (frac - 0.2)/0.2;
	if(frac < 0.8)
		return 0.0;
	return frac/0.2 - 4.0;// (frac - 0.8)/0.2
}

double gFrac(double frac){
	if(frac > 1.0)
		frac = 1.0;
	if(frac < 0.2)
		return frac/0.2;
	if(frac < 0.6)
		return 1.0;
	if(frac < 0.8)
		return 4.0 - frac/0.2;//1.0 - (frac - 0.6)/0.2
	return 0.0;
}

double bFrac(double frac){
	if(frac > 1.0)
		frac = 1.0;
	if(frac < 0.4)
		return 0.0;
	if(frac < 0.6)
		return frac/0.2 - 2.0;// (frac - 0.4)/0.2
	return 1.0;
}

// produces evenly spaces colors
void rainbow(int numColors, unsigned char *domR, unsigned char *domG, unsigned char *domB){
	for(int i(0); i < numColors; i++){
		domR[i] = (unsigned char)floor(255*rFrac(double(i)/double(numColors)));
		domG[i] = (unsigned char)floor(255*gFrac(double(i)/double(numColors)));
		domB[i] = (unsigned char)floor(255*bFrac(double(i)/double(numColors)));
	}
	//cout << "rainbow colors:" << endl;
	//for(int i(0); i < numColors; i++)
	//	cout << "\t" << (int)domR[i] << " " << (int)domG[i] << " " << (int)domB[i] << endl;
}

int trueCount(int size, bool *x){
	int count(0);
	for(int i(0); i < size; i++){
		if(x[i])
			count++;
	}
	return count;
}

void colorDomains(vector<unsigned char> &image, int width, int height, int border, const vector<BranchPoint2D> &b, double scaleX, double scaleY, double offsetX, double offsetY){
	unsigned long ii, jj, kk;
	kisseed(ii, jj, kk);
	kisset(11, 2222, 33333333, false);
	vector<BranchPoint2D> bTips;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].childIndex.size() == 0)
			bTips.push_back(b[i]);
	}
	double *cx = new double[bTips.size()];
	double *cy = new double[bTips.size()];
	for(unsigned int i(0); i < bTips.size(); i++){
		double centerX(scaleX*(bTips[i].position[0] - offsetX) + border),
			centerY(scaleY*(bTips[i].position[1] - offsetY) + border);
		cx[i] = centerX;
		cy[i] = centerY;
	}
	int **dominatedBy = new int*[width];
	for(int x(0); x < width; x++){
		dominatedBy[x] = new int[height];
		for(int y(0); y < height; y++){
			dominatedBy[x][y] = 0;
			double sep(separation2D(x, y, cx[0], cy[0]));
			for(unsigned int i(1); i < bTips.size(); i++){
				double sepTemp(separation2D(x, y, cx[i], cy[i]));
				if(sep > sepTemp){
					dominatedBy[x][y] = i;
					sep = sepTemp;
				}
			}
		}
	}
	bool **neighbors = new bool*[bTips.size()];
	for(unsigned int i(0); i < bTips.size(); i++){
		neighbors[i] = new bool[bTips.size()];
		for(unsigned int j(0); j < bTips.size(); j++)
			neighbors[i][j] = false;
	}
	for(int x(0); x < width; x++){
		for(int y(0); y < height; y++){
			if(x > 0){
				if(dominatedBy[x][y] != dominatedBy[x - 1][y])
					neighbors[dominatedBy[x][y]][dominatedBy[x - 1][y]] = neighbors[dominatedBy[x - 1][y]][dominatedBy[x][y]] = true;
			}
			if(x < width - 1){
				if(dominatedBy[x][y] != dominatedBy[x + 1][y])
					neighbors[dominatedBy[x][y]][dominatedBy[x + 1][y]] = neighbors[dominatedBy[x + 1][y]][dominatedBy[x][y]] = true;
			}

			if(y > 0){
				if(dominatedBy[x][y] != dominatedBy[x][y - 1])
					neighbors[dominatedBy[x][y]][dominatedBy[x][y - 1]] = neighbors[dominatedBy[x][y - 1]][dominatedBy[x][y]] = true;
			}
			if(y < height - 1){
				if(dominatedBy[x][y] != dominatedBy[x][y + 1])
					neighbors[dominatedBy[x][y]][dominatedBy[x][y + 1]] = neighbors[dominatedBy[x][y + 1]][dominatedBy[x][y]] = true;
			}
		}
	}
	int *degrees = new int[bTips.size()];
	for(unsigned int i(0); i < bTips.size(); i++)
		degrees[i] = 0;
	for(unsigned int i(0); i < bTips.size(); i++){
		for(unsigned int j(i + 1); j < bTips.size(); j++){
			if(neighbors[i][j]){
				degrees[i]++;
				degrees[j]++;
			}
		}
	}
	stack<int> colorOrder;
	int removeIndex(0);
	while(removeIndex > -1){
		removeIndex = -1;
		for(unsigned int i(0); i < bTips.size(); i++){
			if(degrees[i] > -1){
				if(removeIndex < 0)
					removeIndex = i;
				else if(degrees[removeIndex] > degrees[i])
					removeIndex = i;
			}
		}
		if(removeIndex > -1){
			colorOrder.push(removeIndex);
			degrees[removeIndex] = -1;
			for(unsigned int i(0); i < bTips.size(); i++){
				if(neighbors[removeIndex][i])
					degrees[i]--;
			}
		}
	}
	delete[] degrees;
	int numColors(4);
	int *domC = new int[bTips.size()];
	for(unsigned int i(0); i < bTips.size(); i++)
		domC[i] = -1;
	bool *colorTaken = new bool[numColors];
	while(!colorOrder.empty()){
		int i(colorOrder.top());
		colorOrder.pop();
		for(int c(0); c < numColors; c++)
			colorTaken[c] = false;
		for(unsigned int j(0); j < bTips.size(); j++){
			if(neighbors[i][j] && domC[j] > -1)
				colorTaken[domC[j]] = true;
		}
		vector<int> availableColors;
		for(int c(0); c < numColors; c++){
			if(!colorTaken[c])
				availableColors.push_back(c);
		}
		if(availableColors.size() < 1){
			cout << "\n Error colorDomains(): no available color!" << endl;
			domC[i] = abs(kiss()%numColors);
		}else
			domC[i] = minimum(availableColors);//;availableColors[abs(int(kiss()%availableColors.size()))];
	}
	delete[] colorTaken;
	/*
	int numColors(6);
	int *domC = new int[bTips.size()];
	for(unsigned int i(0); i < bTips.size(); i++){
		domC[i] = abs(kiss()%numColors);
	}
	int recolorCount(1), passCount(0), maxPasses(1000);
	bool *takenColors = new bool[numColors];
	while(recolorCount > 0 && passCount < maxPasses){
		recolorCount = 0;
		for(unsigned int i(0); i < bTips.size(); i++){
			for(unsigned int j(i + 1); j < bTips.size(); j++){
				if(neighbors[i][j] && domC[i] == domC[j]){
					for(int c(0); c < numColors; c++)
						takenColors[c] = false;
					for(unsigned int k(0); k < bTips.size(); k++){
						if(neighbors[i][k])
							takenColors[domC[k]] = true;
					}
					vector<int> availableColors;
					for(int c(0); c < numColors; c++){
						if(!takenColors[c])
							availableColors.push_back(c);
					}
					if(availableColors.size() < 1)
						domC[i] = abs(kiss()%numColors);
					else
						domC[i] = availableColors[kiss()%availableColors.size()];
					recolorCount++;
				}
			}
		}
		cout << "(" << recolorCount << ")" << endl;
		passCount++;
	}
	delete[] takenColors;
	*/
	/*
	int maxNeighbors(0);
	for(unsigned int i(0); i < bTips.size(); i++){
		int n(trueCount(bTips.size(), neighbors[i]));
		if(maxNeighbors < n)
			maxNeighbors = n;
	}
	for(int targetN(maxNeighbors); targetN > -1; targetN--){
		for(unsigned int i(0); i < bTips.size(); i++){
			if(trueCount(bTips.size(), neighbors[i]) == targetN){
				vector<int> takenColors;
				for(unsigned int j(0); j < bTips.size(); j++){
					if(neighbors[i][j] && domC[j] > -1)
						takenColors.push_back(domC[j]);
				}
				vector<int> possibleColors;
				for(unsigned int j(0); j < 4; j++){
					if(!containedIn(j, takenColors))
						possibleColors.push_back(j);
				}
				if(possibleColors.size() > 0)
					domC[i] = possibleColors[kiss()%possibleColors.size()];
				else{
					for(int j(4); j < maxNeighbors; j++){
						if(!containedIn(j, takenColors)){
							domC[i] = j;
							j = maxNeighbors;
						}
					}
				}
			}
		}
	}
	//for(unsigned int i(0); i < bTips.size(); i++){
	//	vector<int> takenColors;
	//	for(unsigned int j(0); j < i; j++){
	//		if(neighbors[i][j])
	//			takenColors.push_back(domC[j]);
	//	}
	//	domC[i] = 0;
	//	while(containedIn(domC[i], takenColors))
	//		domC[i]++;
	//}
	
	int numColors(0);
	for(unsigned int i(0); i < bTips.size(); i++){
		if(numColors < domC[i] + 1)
			numColors = domC[i] + 1;
	}
	int *colorCounts = new int[numColors];
	for(int c(0); c < numColors; c++)
		colorCounts[c] = 0;
	for(unsigned int i(0); i < bTips.size(); i++)
		colorCounts[domC[i]]++;
	cout << "initialColorCounts:" << endl;
	for(int c(0); c < numColors; c++)
		cout << "\t color " << c << " has " << colorCounts[c] << " instances" << endl;
	// balance colors
	for(unsigned int i(0); i < bTips.size(); i++){
		vector<int> neighborC;
		for(unsigned int j(0); j < bTips.size(); j++){
			if(j == i)
				continue;
			if(neighbors[i][j] && !containedIn(domC[j], neighborC))
				neighborC.push_back(domC[j]);
		}
		int minC(-1);
		for(int c(0); c < numColors; c++){
			if(!containedIn(c, neighborC)){
				minC = c;
				c = numColors;
			}
		}
		if(minC > -1){
			for(int c(0); c < numColors; c++){
				if(colorCounts[minC] > colorCounts[c] && !containedIn(c, neighborC))
					minC = c;
			}
			if(colorCounts[domC[i]] > colorCounts[domC[minC]] + 1){
				colorCounts[domC[i]]--;
				domC[i] = minC;
				colorCounts[domC[i]]++;
			}
		}
	}*/
	unsigned char *domR = new unsigned char[numColors];
	unsigned char *domG = new unsigned char[numColors];
	unsigned char *domB = new unsigned char[numColors];
	double *sat = new double[numColors]; // saturation
	for(int c(0); c < numColors; c++)
		sat[c] = 1.0 - (0.5*double(c)/double(numColors - 1));
	//rainbow(numColors, domR, domG, domB);
	domR[0] = 230; domG[0] = 97; domB[0] = 1;
	domR[3] = 253; domG[3] = 184; domB[3] = 99;
	domR[2] = 178; domG[2] = 171; domB[2] = 210;
	domR[1] = 94; domG[1] = 60; domB[1] = 153;
	for(int x(0*border); x < width - 0*border; x++){
		for(int y(0*border); y < height - 0*border; y++){
			double s(1.0);//s(sat[domC[dominatedBy[x][y]]]);
			image[imageIndex(x, y, width, height)] = (unsigned char)floor(s*domR[domC[dominatedBy[x][y]]]);
			image[imageIndex(x, y, width, height) + 1] = (unsigned char)floor(s*domG[domC[dominatedBy[x][y]]]);
			image[imageIndex(x, y, width, height) + 2] = (unsigned char)floor(s*domB[domC[dominatedBy[x][y]]]);
			image[imageIndex(x, y, width, height) + 3] = 100;//(unsigned char)floor(20 + 120*double(domC[dominatedBy[x][y]])/numColors);
		}
	}

	for(int x(0); x < width; x++)
		delete[] dominatedBy[x];
	delete[] dominatedBy;
	for(unsigned int i(0); i < bTips.size(); i++)
		delete[] neighbors[i];
	delete[] neighbors;
	delete[] cx;
	delete[] cy;
	delete[] domC;
	delete[] domR;
	delete[] domG;
	delete[] domB;
	//delete[] colorCounts;
	delete[] sat;
	kisset(ii, jj, kk, false);
}

void colorDomainsCirc(vector<unsigned char> &image, int width, int height, int border, const vector<BranchPoint2D> &b, double scaleX, double scaleY, double offsetX, double offsetY){
	unsigned long ii, jj, kk;
	kisseed(ii, jj, kk);
	kisset(11, 2222, 33333333, false);
	vector<BranchPoint2D> bTips;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].childIndex.size() == 0)
			bTips.push_back(b[i]);
	}
	double *cx = new double[bTips.size()];
	double *cy = new double[bTips.size()];
	for(unsigned int i(0); i < bTips.size(); i++){
		double centerX(scaleX*(bTips[i].position[0] - offsetX) + border),
			centerY(scaleY*(bTips[i].position[1] - offsetY) + border);
		cx[i] = centerX;
		cy[i] = centerY;
	}
	int **dominatedBy = new int*[width];
	for(int x(0); x < width; x++){
		dominatedBy[x] = new int[height];
		for(int y(0); y < height; y++){
			dominatedBy[x][y] = 0;
			double sep(separation2D(x, y, cx[0], cy[0]));
			for(unsigned int i(1); i < bTips.size(); i++){
				double sepTemp(separation2D(x, y, cx[i], cy[i]));
				if(sep > sepTemp){
					dominatedBy[x][y] = i;
					sep = sepTemp;
				}
			}
		}
	}
	bool **neighbors = new bool*[bTips.size()];
	for(unsigned int i(0); i < bTips.size(); i++){
		neighbors[i] = new bool[bTips.size()];
		for(unsigned int j(0); j < bTips.size(); j++)
			neighbors[i][j] = false;
	}
	for(int x(0); x < width; x++){
		for(int y(0); y < height; y++){
			if(x > 0){
				if(dominatedBy[x][y] != dominatedBy[x - 1][y])
					neighbors[dominatedBy[x][y]][dominatedBy[x - 1][y]] = neighbors[dominatedBy[x - 1][y]][dominatedBy[x][y]] = true;
			}
			if(x < width - 1){
				if(dominatedBy[x][y] != dominatedBy[x + 1][y])
					neighbors[dominatedBy[x][y]][dominatedBy[x + 1][y]] = neighbors[dominatedBy[x + 1][y]][dominatedBy[x][y]] = true;
			}

			if(y > 0){
				if(dominatedBy[x][y] != dominatedBy[x][y - 1])
					neighbors[dominatedBy[x][y]][dominatedBy[x][y - 1]] = neighbors[dominatedBy[x][y - 1]][dominatedBy[x][y]] = true;
			}
			if(y < height - 1){
				if(dominatedBy[x][y] != dominatedBy[x][y + 1])
					neighbors[dominatedBy[x][y]][dominatedBy[x][y + 1]] = neighbors[dominatedBy[x][y + 1]][dominatedBy[x][y]] = true;
			}
		}
	}
	int *degrees = new int[bTips.size()];
	for(unsigned int i(0); i < bTips.size(); i++)
		degrees[i] = 0;
	for(unsigned int i(0); i < bTips.size(); i++){
		for(unsigned int j(i + 1); j < bTips.size(); j++){
			if(neighbors[i][j]){
				degrees[i]++;
				degrees[j]++;
			}
		}
	}
	stack<int> colorOrder;
	int removeIndex(0);
	while(removeIndex > -1){
		removeIndex = -1;
		for(unsigned int i(0); i < bTips.size(); i++){
			if(degrees[i] > -1){
				if(removeIndex < 0)
					removeIndex = i;
				else if(degrees[removeIndex] > degrees[i])
					removeIndex = i;
			}
		}
		if(removeIndex > -1){
			colorOrder.push(removeIndex);
			degrees[removeIndex] = -1;
			for(unsigned int i(0); i < bTips.size(); i++){
				if(neighbors[removeIndex][i])
					degrees[i]--;
			}
		}
	}
	delete[] degrees;
	int numColors(4);
	int *domC = new int[bTips.size()];
	for(unsigned int i(0); i < bTips.size(); i++)
		domC[i] = -1;
	bool *colorTaken = new bool[numColors];
	while(!colorOrder.empty()){
		int i(colorOrder.top());
		colorOrder.pop();
		for(int c(0); c < numColors; c++)
			colorTaken[c] = false;
		for(unsigned int j(0); j < bTips.size(); j++){
			if(neighbors[i][j] && domC[j] > -1)
				colorTaken[domC[j]] = true;
		}
		vector<int> availableColors;
		for(int c(0); c < numColors; c++){
			if(!colorTaken[c])
				availableColors.push_back(c);
		}
		if(availableColors.size() < 1){
			cout << "\n Error colorDomains(): no available color!" << endl;
			domC[i] = abs(kiss()%numColors);
		}else
			domC[i] = minimum(availableColors);//;availableColors[abs(int(kiss()%availableColors.size()))];
	}
	delete[] colorTaken;
	unsigned char *domR = new unsigned char[numColors];
	unsigned char *domG = new unsigned char[numColors];
	unsigned char *domB = new unsigned char[numColors];
	double *sat = new double[numColors]; // saturation
	for(int c(0); c < numColors; c++)
		sat[c] = 1.0 - (0.5*double(c)/double(numColors - 1));
	//rainbow(numColors, domR, domG, domB);
	domR[0] = 230; domG[0] = 97; domB[0] = 1;
	domR[3] = 253; domG[3] = 184; domB[3] = 99;
	domR[2] = 178; domG[2] = 171; domB[2] = 210;
	domR[1] = 94; domG[1] = 60; domB[1] = 153;
	for(int x(0*border); x < width - 0*border; x++){
		for(int y(0*border); y < height - 0*border; y++){
			double s(1.0);//s(sat[domC[dominatedBy[x][y]]]);
			if(width*width/4 > (x - width/2)*(x - width/2) + (y - height/2)*(y - height/2)){
				image[imageIndex(x, y, width, height)] = (unsigned char)floor(s*domR[domC[dominatedBy[x][y]]]);
				image[imageIndex(x, y, width, height) + 1] = (unsigned char)floor(s*domG[domC[dominatedBy[x][y]]]);
				image[imageIndex(x, y, width, height) + 2] = (unsigned char)floor(s*domB[domC[dominatedBy[x][y]]]);
				image[imageIndex(x, y, width, height) + 3] = 100;//(unsigned char)floor(20 + 120*double(domC[dominatedBy[x][y]])/numColors);
			}
		}
	}

	for(int x(0); x < width; x++)
		delete[] dominatedBy[x];
	delete[] dominatedBy;
	for(unsigned int i(0); i < bTips.size(); i++)
		delete[] neighbors[i];
	delete[] neighbors;
	delete[] cx;
	delete[] cy;
	delete[] domC;
	delete[] domR;
	delete[] domG;
	delete[] domB;
	delete[] sat;
	kisset(ii, jj, kk, false);
}

void drawTree(string fn, const vector<BranchPoint2D> &b, int width, int height, int border){
	vector<unsigned char> image;
	image.resize(width*height*4);
	for(int i(0); i < width*height*4; i++)
		image[i] = 0;
	double offsetX(minX(b)), offsetY(minY(b)), scaleX((width - 2*border)/(maxX(b) - offsetX)), scaleY((height - 2*border)/(maxY(b) - offsetY));
	if(scaleY > scaleX)
		scaleY = scaleX;
	if(scaleX > scaleY)
		scaleX = scaleY;
	//cout << "\ndrawBifurcationTree(): offset = (" << offsetX << ", " << offsetY << "), scale = (" << scaleX << ", " << scaleY << ")" << endl;
	int red(111), green(55), blue(55), // for lines
		redTip(0), greenTip(111), blueTip(0),
		redNode(0), greenNode(0), blueNode(111),
		redHeart(111), greenHeart(0), blueHeart(0);
	colorDomains(image, width, height, border, b, scaleX, scaleY, offsetX, offsetY);
	for(unsigned int i(0); i < b.size(); i++){
		double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border);
		for(unsigned int c(0); c < b[i].childIndex.size(); c++){
			double startX(centerX),//scaleX*(b[i].position[0] - offsetX) + border),
				startY(centerY),//scaleY*(b[i].position[1] - offsetY) + border),
				endX(scaleX*(b[b[i].childIndex[c]].position[0] - offsetX) + border),
				endY(scaleY*(b[b[i].childIndex[c]].position[1] - offsetY) + border);
			drawLine(image, width, height, red, green, blue, startX, startY, endX, endY);
		}
	}
	for(unsigned int i(0); i < b.size(); i++){
		double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border),
			standardRadius(border/2);
		if(i == 0)
			drawTriangle(image, width, height, redHeart, greenHeart, blueHeart, centerX, centerY, 2*standardRadius);
		else if(b[i].childIndex.size() < 1)
			drawCircle(image, width, height, red, green, blue, centerX, centerY, double(width)*2.0/300.0, double(width)*2.0/300.0);
		//	drawSquare(image, width, height, redTip, greenTip, blueTip, centerX, centerY, standardRadius);
		//else
		//	drawCircle(image, width, height, redNode, greenNode, blueNode, centerX, centerY, standardRadius*4/8 + (int)floor(0*1.5*maxLevel(b, i)));
		
	}
	unsigned error = lodepng::encode(fn.c_str(), image, width, height);
	if(error)
		cout << "\nError snapshotPNG(): error " << error << " from lodepng" << endl;
}

void drawTreeCirc(string fn, const vector<BranchPoint2D> &b, int width, int height, int border){
	vector<unsigned char> image;
	image.resize(width*height*4);
	for(int i(0); i < width*height*4; i++)
		image[i] = 0;
	double offsetX(minX(b)), offsetY(minY(b)), scaleX((width - 2*border)/(maxX(b) - offsetX)), scaleY((height - 2*border)/(maxY(b) - offsetY));
	offsetX -= 0.05*width/scaleX;
	offsetY -= 0.05*height/scaleY;
	scaleX *= 0.9;
	scaleY *= 0.9;
	if(scaleY > scaleX)
		scaleY = scaleX;
	if(scaleX > scaleY)
		scaleX = scaleY;
	//cout << "\ndrawBifurcationTree(): offset = (" << offsetX << ", " << offsetY << "), scale = (" << scaleX << ", " << scaleY << ")" << endl;
	int red(111), green(55), blue(55), // for lines
		redTip(0), greenTip(111), blueTip(0),
		redNode(0), greenNode(0), blueNode(111),
		redHeart(111), greenHeart(0), blueHeart(0);
	colorDomainsCirc(image, width, height, border, b, scaleX, scaleY, offsetX, offsetY);
	for(unsigned int i(0); i < b.size(); i++){
		double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border);
		for(unsigned int c(0); c < b[i].childIndex.size(); c++){
			double startX(centerX),//scaleX*(b[i].position[0] - offsetX) + border),
				startY(centerY),//scaleY*(b[i].position[1] - offsetY) + border),
				endX(scaleX*(b[b[i].childIndex[c]].position[0] - offsetX) + border),
				endY(scaleY*(b[b[i].childIndex[c]].position[1] - offsetY) + border);
			drawLine(image, width, height, red, green, blue, startX, startY, endX, endY);
		}
	}
	for(unsigned int i(0); i < b.size(); i++){
		double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border),
			standardRadius(border/2);
		if(i == 0)
			drawTriangle(image, width, height, redHeart, greenHeart, blueHeart, centerX, centerY, 2*standardRadius);
		else if(b[i].childIndex.size() < 1)
			drawCircle(image, width, height, red, green, blue, centerX, centerY, double(width)*2.0/300.0, double(width)*2.0/300.0);
		//	drawSquare(image, width, height, redTip, greenTip, blueTip, centerX, centerY, standardRadius);
		//else
		//	drawCircle(image, width, height, redNode, greenNode, blueNode, centerX, centerY, standardRadius*4/8 + (int)floor(0*1.5*maxLevel(b, i)));
		
	}
	unsigned error = lodepng::encode(fn.c_str(), image, width, height);
	if(error)
		cout << "\nError snapshotPNG(): error " << error << " from lodepng" << endl;
}

void drawBifurcationTree(string fn, vector<Bifurcation2D> &b, vector<int> emphasis, int width, int height, int border){
	vector<unsigned char> image;
	image.resize(width*height*4);
	for(int i(0); i < width*height*4; i++)
		image[i] = 0;
	double offsetX(minX(b)), offsetY(minY(b)), scaleX((width - 2*border)/(maxX(b) - offsetX)), scaleY((height - 2*border)/(maxY(b) - offsetY));
	//cout << "\ndrawBifurcationTree(): offset = (" << offsetX << ", " << offsetY << "), scale = (" << scaleX << ", " << scaleY << ")" << endl;
	int red(255), green(55), blue(155),
		redTip(0), greenTip(255), blueTip(155),
		redNode(55), greenNode(155), blueNode(255),
		redHeart(255), greenHeart(55), blueHeart(55),
		redEmph(255), greenEmph(165), blueEmph(0);
	for(unsigned int i(0); i < b.size(); i++){
		double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border),
			standardRadius(border/2);
		if(containedIn(i, emphasis))
			drawCircle(image, width, height, redEmph, greenEmph, blueEmph, centerX, centerY, 3*standardRadius/2);
		else if(i == 0)
			drawCircle(image, width, height, redHeart, greenHeart, blueHeart, centerX, centerY, 2*standardRadius);
		else if(b[i].childIndex.size() < 1)
			drawCircle(image, width, height, redTip, greenTip, blueTip, centerX, centerY, standardRadius);
		else
			drawCircle(image, width, height, redNode, greenNode, blueNode, centerX, centerY, standardRadius/8 + (int)floor(1.5*maxLevel(b, i)));
		for(unsigned int c(0); c < b[i].childIndex.size(); c++){
			double startX(centerX),//scaleX*(b[i].position[0] - offsetX) + border),
				startY(centerY),//scaleY*(b[i].position[1] - offsetY) + border),
				endX(scaleX*(b[b[i].childIndex[c]].position[0] - offsetX) + border),
				endY(scaleY*(b[b[i].childIndex[c]].position[1] - offsetY) + border);
			drawLine(image, width, height, red, green, blue, startX, startY, endX, endY);
		}
	}
	unsigned error = lodepng::encode(fn.c_str(), image, width, height);
	if(error)
		cout << "\nError snapshotPNG(): error " << error << " from lodepng" << endl;
}

void drawTree(string fn, const vector<BranchPoint3D> &b, int width, int height, int border){
	vector<unsigned char> image;
	image.resize(width*height*4);
	for(int i(0); i < width*height*4; i++)
		image[i] = 0;
	double offsetX(minX(b)), offsetY(minY(b)), offsetZ(minZ(b)), largestZ(maxZ(b)),
		scaleX((width - 2*border)/(maxX(b) - offsetX)), scaleY((height - 2*border)/(maxY(b) - offsetY));
	if(scaleY > scaleX)
		scaleY = scaleX;
	if(scaleX > scaleY)
		scaleX = scaleY;
	//cout << "\ndrawBifurcationTree(): offset = (" << offsetX << ", " << offsetY << "), scale = (" << scaleX << ", " << scaleY << ")" << endl;
	int red(111), green(55), blue(55), // for lines; note that g will vary within lines
		redTip(0), greenTip(111), blueTip(0),
		redNode(0), greenNode(55), blueNode(111), // note that g will vary across nodes
		redHeart(111), greenHeart(0), blueHeart(0);
	//colorDomains(image, width, height, border, b, scaleX, scaleY, offsetX, offsetY);
	for(unsigned int i(0); i < b.size(); i++){
		double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border);
		for(unsigned int c(0); c < b[i].childIndex.size(); c++){
			double startX(centerX),//scaleX*(b[i].position[0] - offsetX) + border),
				startY(centerY),//scaleY*(b[i].position[1] - offsetY) + border),
				endX(scaleX*(b[b[i].childIndex[c]].position[0] - offsetX) + border),
				endY(scaleY*(b[b[i].childIndex[c]].position[1] - offsetY) + border);
			int greenStart((int)floor(green*(1.0 + (b[i].position[2] - offsetZ)/(largestZ - offsetZ)))),
				greenEnd((int)floor(green*(1.0 + (b[b[i].childIndex[c]].position[2] - offsetZ)/(largestZ - offsetZ))));
			drawLineNoOverwriteLessGreen(image, width, height, red, red, greenStart, greenEnd, blue, blue, startX, startY, endX, endY);
		}
	}
	for(unsigned int i(0); i < b.size(); i++){
		double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border),
			standardRadius(border/2);
		if(i == 0)
			drawTriangle(image, width, height, redHeart, greenHeart, blueHeart, centerX, centerY, 2*standardRadius);
		else if(b[i].childIndex.size() < 1)
			drawCircle(image, width, height, red, green, blue, centerX, centerY, double(width)*2.0/300.0, double(width)*2.0/300.0);

		//else if(b[i].childIndex.size() < 1)
		//	drawSquare(image, width, height, redTip, greenTip, blueTip, centerX, centerY, standardRadius);
		//else
		//	drawCircle(image, width, height, redNode, (int)floor(greenNode*(1.0 + (b[i].position[2] - offsetZ)/(largestZ - offsetZ))), blueNode, centerX, centerY, standardRadius*4/8 + (int)floor(0*1.5*maxLevel(b, i)));
		
	}
	unsigned error = lodepng::encode(fn.c_str(), image, width, height);
	if(error)
		cout << "\nError snapshotPNG(): error " << error << " from lodepng" << endl;
}

void drawBentTree(string fn, const vector<BranchPoint2D> &b, double bendX, double bendY, int width, int height, int border){
	vector<unsigned char> image;
	image.resize(width*height*4);
	for(int i(0); i < width*height*4; i++)
		image[i] = 0;
	double offsetX(minX(b)), offsetY(minY(b)), scaleX((width - 2*border)/(maxX(b) - offsetX)), scaleY((height - 2*border)/(maxY(b) - offsetY));
	if(scaleY > scaleX)
		scaleY = scaleX;
	if(scaleX > scaleY)
		scaleX = scaleY;
	//cout << "\ndrawBifurcationTree(): offset = (" << offsetX << ", " << offsetY << "), scale = (" << scaleX << ", " << scaleY << ")" << endl;
	int red(111), green(55), blue(55), // for lines
		redTip(0), greenTip(111), blueTip(0),
		redNode(0), greenNode(0), blueNode(111),
		redHeart(111), greenHeart(0), blueHeart(0);
	colorDomains(image, width, height, border, b, scaleX, scaleY, offsetX, offsetY);
	for(int x(0); x < scaleX*(bendX - offsetX); x++){
		for(int y(0); y < scaleY*(bendY - offsetY); y++){
			int i(imageIndex(x, y, width, height));
			image[i] = image[i + 1] = image[i + 2] = image[i + 3] = 0;
		}
	}
	drawLine(image, width, height, 0, 0, 0, 0, height - 1, width - 1, height - 1, 0.5);
	drawLine(image, width, height, 0, 0, 0, width - 1, height - 1, width - 1, 0, 0.5);
	drawLine(image, width, height, 0, 0, 0, width - 1, 0, (int)floor(scaleX*(bendX - offsetX)), 0, 0.5);
	drawLine(image, width, height, 0, 0, 0, (int)floor(scaleX*(bendX - offsetX)), 0, (int)floor(scaleX*(bendX - offsetX)), (int)floor(scaleY*(bendY - offsetY)), 0.5);
	drawLine(image, width, height, 0, 0, 0, (int)floor(scaleX*(bendX - offsetX)), (int)floor(scaleY*(bendY - offsetY)), 0, (int)floor(scaleY*(bendY - offsetY)), 0.5);
	drawLine(image, width, height, 0, 0, 0, 0, (int)floor(scaleY*(bendY - offsetY)), 0, height - 1, 0.5);
	for(unsigned int i(0); i < b.size(); i++){
		double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border);
		for(unsigned int c(0); c < b[i].childIndex.size(); c++){
			double startX(centerX),//scaleX*(b[i].position[0] - offsetX) + border),
				startY(centerY),//scaleY*(b[i].position[1] - offsetY) + border),
				endX(scaleX*(b[b[i].childIndex[c]].position[0] - offsetX) + border),
				endY(scaleY*(b[b[i].childIndex[c]].position[1] - offsetY) + border);
			drawLine(image, width, height, red, green, blue, startX, startY, endX, endY);
		}
	}
	for(unsigned int i(0); i < b.size(); i++){
		double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border),
			standardRadius(border/2);
		if(i == 0)
			drawTriangle(image, width, height, redHeart, greenHeart, blueHeart, centerX, centerY, 2*standardRadius);
		else if(b[i].childIndex.size() < 1)
			drawSquare(image, width, height, redTip, greenTip, blueTip, centerX, centerY, standardRadius);
		else
			drawCircle(image, width, height, redNode, greenNode, blueNode, centerX, centerY, standardRadius*4/8 + (int)floor(0*1.5*maxLevel(b, i)));
		
	}
	unsigned error = lodepng::encode(fn.c_str(), image, width, height);
	if(error)
		cout << "\nError snapshotPNG(): error " << error << " from lodepng" << endl;
}

void drawBifurcationTree(string fn, vector<Bifurcation3D> &b, vector<int> emphasis, int width, int height, int border){
	vector<unsigned char> image;
	image.resize(width*height*4);
	for(int i(0); i < width*height*4; i++)
		image[i] = 0;
	double offsetX(minX(b)), offsetY(minY(b)), scaleX((width - 2*border)/(maxX(b) - offsetX)), scaleY((height - 2*border)/(maxY(b) - offsetY));
	//cout << "\ndrawBifurcationTree(): offset = (" << offsetX << ", " << offsetY << "), scale = (" << scaleX << ", " << scaleY << ")" << endl;
	int red(255), green(55), blue(155),
		redTip(0), greenTip(255), blueTip(155),
		redNode(55), greenNode(155), blueNode(255),
		redHeart(255), greenHeart(55), blueHeart(55),
		redEmph(255), greenEmph(165), blueEmph(0);
	for(unsigned int i(0); i < b.size(); i++){
		double centerX(scaleX*(b[i].position[0] - offsetX) + border),
			centerY(scaleY*(b[i].position[1] - offsetY) + border),
			standardRadius(border/2);
		if(containedIn(i, emphasis))
			drawCircle(image, width, height, redEmph, greenEmph, blueEmph, centerX, centerY, 3*standardRadius/2);
		else if(i == 0)
			drawCircle(image, width, height, redHeart, greenHeart, blueHeart, centerX, centerY, 2*standardRadius);
		else if(b[i].childIndex.size() < 1)
			drawCircle(image, width, height, redTip, greenTip, blueTip, centerX, centerY, standardRadius);
		else
			drawCircle(image, width, height, redNode, greenNode, blueNode, centerX, centerY, standardRadius/8 + (int)floor(1.5*maxLevel(b, i)));
		for(unsigned int c(0); c < b[i].childIndex.size(); c++){
			double startX(centerX),//scaleX*(b[i].position[0] - offsetX) + border),
				startY(centerY),//scaleY*(b[i].position[1] - offsetY) + border),
				endX(scaleX*(b[b[i].childIndex[c]].position[0] - offsetX) + border),
				endY(scaleY*(b[b[i].childIndex[c]].position[1] - offsetY) + border);
			drawLine(image, width, height, red, green, blue, startX, startY, endX, endY);
		}
	}
	unsigned error = lodepng::encode(fn.c_str(), image, width, height);
	if(error)
		cout << "\nError snapshotPNG(): error " << error << " from lodepng" << endl;
}

string makeStringTree(vector<Bifurcation2D> b){
	stringstream strstr;
	strstr << endl << " i\t x\t y\t p\t c..." << endl;
	for(unsigned int i(0); i < b.size(); i++){
		strstr << i << "\t" << b[i].position[0] << "\t" << b[i].position[1] << "\t" << b[i].parentIndex;
		for(unsigned int c(0); c < b[i].childIndex.size(); c++)
			strstr << "\t" << b[i].childIndex[c];
		strstr << endl;
	}
	return strstr.str();
}

string makeStringTree(vector<BranchPoint2D> bc){
	stringstream strstr;
	strstr << endl << " i\t x\t y\t p\t c..." << endl;
	for(unsigned int i(0); i < bc.size(); i++){
		strstr << i << "\t" << bc[i].position[0] << "\t" << bc[i].position[1] << "\t" << bc[i].parentIndex;
		for(unsigned int c(0); c < bc[i].childIndex.size(); c++)
			strstr << "\t" << bc[i].childIndex[c];
		strstr << endl;
	}
	return strstr.str();
}

double cot(double x){ return 1.0/tan(x); }
double acot(double x){ return atan(1.0/x); }

double angleBetweenGrandchildren(int grandIndex, int child1Index, int child2Index, vector<Bifurcation2D> &b);
double angleBetweenGrandchildren3D(int grandIndex, int child1Index, int child2Index, vector<Bifurcation3D> &b);

double crossProduct2DSign(vector<Bifurcation2D> &b, int p, int c1, int c2){
	double u1(b[c1].position[0] - b[p].position[0]),
		u2(b[c1].position[1] - b[p].position[1]),
		v1(b[c2].position[0] - b[p].position[0]),
		v2(b[c2].position[1] - b[p].position[1]);
	double cProd(u1*v2 - u2*v1);
	if(cProd == 0.0)
		return 0.0;
	if(cProd < 0.0)
		return -1.0;
	return 1.0;
}

inline double areaHeron(double a, double b, double c){
	double s((a + b + c)/2.0);
	return sqrt(s*(s - a)*(s - b)*(s - c));
}

// Using the definition from http://faculty.evansville.edu/ck6/encyclopedia/ETC.html
inline double barycentricFermat(double a, double b, double c){
	double aSq(a*a), bSq(b*b), cSq(c*c);
	return aSq*aSq - 2.0*(bSq - cSq)*(bSq - cSq) + aSq*(bSq + cSq + 4.0*sqrt(3.0)*areaHeron(a, b, c));
}

vector<BranchPoint2D> consolidateDegenerates(const vector<Bifurcation2D> &bDegen, double thresh = 1.0e-6);

void settleBifurcation(int i, vector<Bifurcation2D> &b, double &x, double &y, double thresh = 1.0e-6, bool exclaim = false){
	if(b[i].childIndex.size() != 2){
		cout << "\n Error settleBifurcation()!: expected exactly 2 children" << endl;
		x = b[i].position[0];
		y = b[i].position[1];
		return;
	}
	if(b[i].parentIndex < 0){
		x = (b[b[i].childIndex[0]].position[0] + b[b[i].childIndex[1]].position[0])/2.0;
		y = (b[b[i].childIndex[0]].position[1] + b[b[i].childIndex[1]].position[1])/2.0;
		return;
	}
	double critAngle(2.0*acos(-1.0)/3.0);
	if(angleBetweenGrandchildren(b[i].parentIndex, b[i].childIndex[0], b[i].childIndex[1], b) >= critAngle){
		x = b[b[i].parentIndex].position[0];
		y = b[b[i].parentIndex].position[1];
		return;
	}else if(angleBetweenGrandchildren(b[i].childIndex[0], b[i].parentIndex, b[i].childIndex[1], b) >= critAngle){
		x = b[b[i].childIndex[0]].position[0];
		y = b[b[i].childIndex[0]].position[1];
		return;
	}else if(angleBetweenGrandchildren(b[i].childIndex[1], b[i].childIndex[0], b[i].parentIndex, b) >= critAngle){
		x = b[b[i].childIndex[1]].position[0];
		y = b[b[i].childIndex[1]].position[1];
		return;
	}
	double sideC(separation2D(b[b[i].parentIndex], b[b[i].childIndex[0]])),
		sideA(separation2D(b[b[i].childIndex[1]], b[b[i].childIndex[0]])),
		sideB(separation2D(b[b[i].parentIndex], b[b[i].childIndex[1]]));
	if(exclaim){
		cout << "\n settleBifurcation():\n\t " << i << ": sideA = " << sideA
			<< "\n\t sideB = " << sideB
			<< "\n\t sideC = " << sideC;
		if(sideA > sideB && sideA > sideC)
			cout << "\n\t\t A-B-C=" << sideA - sideB - sideC;
		else if(sideB > sideA && sideB > sideC)
			cout << "\n\t\t B-A-C=" << sideB - sideA - sideC;
		else
			cout << "\n\t\t C-A-B=" << sideC - sideA - sideB;
		cout << endl;
	}
	double sideThresh = thresh/10.0;
	if(sideA <= sideThresh && sideB > sideThresh && sideC > sideThresh){ // children overlap
		x = b[b[i].childIndex[0]].position[0];
		y = b[b[i].childIndex[0]].position[1];
		return;
	}
	if(sideB <= sideThresh || sideC <= sideThresh){ // at least one child is at the parent
		x = b[b[i].parentIndex].position[0];
		y = b[b[i].parentIndex].position[1];
		return;
	}
	if(abs(sideA - sideB - sideC) < sideThresh || abs(sideB - sideA - sideC) < sideThresh || abs(sideC - sideB - sideA) < sideThresh){ //zero area
		if(sideA >= sideB && sideA >= sideC){
			x = b[b[i].parentIndex].position[0];
			y = b[b[i].parentIndex].position[1];
			if(exclaim)
				cout << "\n sending " << i << " to parent" << endl;
		}else if(sideB >= sideA && sideB >= sideC){
			x = b[b[i].childIndex[0]].position[0];
			y = b[b[i].childIndex[0]].position[1];
			if(exclaim)
				cout << "\n sending " << i << " to child 0" << endl;
		}else{
			x = b[b[i].childIndex[1]].position[0];
			y = b[b[i].childIndex[1]].position[1];
			if(exclaim)
				cout << "\n sending " << i << " to child 1" << endl;
		}
		
		return;
	}
	double p(barycentricFermat(sideA, sideB, sideC)), q(barycentricFermat(sideB, sideC, sideA)), r(barycentricFermat(sideC, sideA, sideB));
	//if(p < 0.0 || q < 0.0 || r < 0.0)
	//	cout << "\n(" << p << ", " << q << ", " << r << ")\n\t{" << sideA << ", " << sideB << ", " << sideC << "}" << "\n\t\t[" << abs(sideA - sideB - sideC) << "]\n\t\t" << sideThresh;
	double sum(p + q + r);
	x = (p*b[b[i].parentIndex].position[0] + q*b[b[i].childIndex[0]].position[0] + r*b[b[i].childIndex[1]].position[0])/sum;
	y = (p*b[b[i].parentIndex].position[1] + q*b[b[i].childIndex[0]].position[1] + r*b[b[i].childIndex[1]].position[1])/sum;
	//cout << "\n\t [" << x << ", " << y << "]";
	//cout << "\nsettleBifurcation(): sides are " << sideA << "\t" << sideB << "\t" << sideC
	//	<< "\n\t p c1 c2 indices are " << b[i].parentIndex << "\t" << b[i].childIndex[0] << "\t" << b[i].childIndex[1]
	//	<< "\n\t barycentric coordinates are " << p << "\t" << q << "\t" << r
	//	<< "\n\t sum = " << sum
	//	<< "\n\t x = " << x << "\t y = " << y
	//	<< endl;

	/*
	int maxIts(1000), count(0);
	double bestx(b[b[i].parentIndex].position[0]/(1 + b[i].childIndex.size())),
		besty(b[b[i].parentIndex].position[0]/(1 + b[i].childIndex.size())),
		changeInBest(2.0*thresh), epsilon(thresh/2.0);
	for(unsigned int c(0); c < b[i].childIndex.size(); c++){
		bestx += b[b[i].childIndex[c]].position[0]/(1 + b[i].childIndex.size());
		besty += b[b[i].childIndex[c]].position[1]/(1 + b[i].childIndex.size());
	}
	//double bestx(b[i].position[0]), besty(b[i].position[1]), changeInBest(2.0*thresh), epsilon(0.001);
	//if(thresh > 1.0)
	//	epsilon = 2.0*thresh;

	while(changeInBest > thresh && count < maxIts){
		double separation(separation2D(bestx, besty, b[b[i].parentIndex].position[0], b[b[i].parentIndex].position[1]) + epsilon);//separation(separation2D(b[i], b[b[i].parentIndex]) + 0.001);
		double sumInverseSeparations(1.0/separation),
			xTemp(b[b[i].parentIndex].position[0]/separation),
			yTemp(b[b[i].parentIndex].position[1]/separation);
		for(unsigned int c(0); c < b[i].childIndex.size(); c++){
			separation = separation2D(bestx, besty, b[b[i].childIndex[c]].position[0], b[b[i].childIndex[c]].position[1]) + epsilon;//separation2D(b[i], b[b[i].childIndex[c]]) + 0.001;
			sumInverseSeparations += 1.0/separation;
			xTemp += b[b[i].childIndex[c]].position[0]/separation;
			yTemp += b[b[i].childIndex[c]].position[1]/separation;
		}
		xTemp /= sumInverseSeparations;
		yTemp /= sumInverseSeparations;
		double changeInBest(sqrt((xTemp - bestx)*(xTemp - bestx) + (yTemp - besty)*(yTemp - besty)));
		bestx = xTemp;
		besty = yTemp;
		count++;
	}
	x = bestx;
	y = besty;

	//if(crossProduct2DSign(b, b[i].parentIndex, b[i].childIndex[0], b[i].childIndex[1]) > 0.0){
	//	int temp(b[i].childIndex[0]);
	//	b[i].childIndex[0] = b[i].childIndex[1];
	//	b[i].childIndex[1] = temp;
	//}
	//double AB(separation2D(b[b[i].parentIndex], b[b[i].childIndex[0]])),
	//	AC(separation2D(b[b[i].parentIndex], b[b[i].childIndex[1]])),
	//	BC(separation2D(b[b[i].childIndex[0]], b[b[i].childIndex[1]])),
	//	phi(angleBetweenGrandchildren(b[i].parentIndex, b[i].childIndex[0], b[i].childIndex[1], b)),
	//	theta1(acos(-1.0)/3.0),
	//	theta2(theta1);
	//double theta(theta1 + theta2);
	//if(AB == 0.0)
	//	cout << "\nsettleBifurcations(): AB is zero!" << endl;
	//if(AC == 0.0)
	//	cout << "\nsettleBifurcations(): AC is zero!" << endl;
	//if(BC == 0.0)
	//	cout << "\nsettleBifurcations(): BC is zero!" << endl;
	//if(sin(theta1) == 0.0)
	//	cout << "\nsettleBifurcations(): sin(theta1) is zero!" << endl;
	//if(sin(theta - phi) == 0.0)
	//	cout << "\nsettleBifurcations(): sin(theta - phi) is zero!" << endl;
	//double gamma(acot((AB/AC)*(sin(theta2)/sin(theta1)/sin(theta - phi)) + cot(theta - phi)));
	//double chi((AC/AB)*(sin(theta2)/sin(theta1))*(sin(phi)*cot(gamma + phi - theta1) - cos(phi))); // L1/L2
	//if(chi == 0.0)
	//	cout << "\nsettleBifurcations(): chi is zero!" << endl;
	//double L1(chi*BC/sqrt(1.0 - 2.0*chi*cos(theta) + chi*chi));
	//double L2(L1/chi);//,
	//	//L0(L1*sin(theta1)/sin(theta1 - gamma));
	//double Cx(b[b[i].childIndex[1]].position[0]),
	//	Cy(b[b[i].childIndex[1]].position[1]),
	//	cosPhi2(((L2/BC) + (BC/L2) - (L1*chi/BC))/2.0),
	//	sinPhi2(L1*sin(theta)/BC);
	//double BCx(b[b[i].childIndex[0]].position[0] - Cx),
	//	BCy(b[b[i].childIndex[0]].position[1] - Cy);
	//x = Cx + (BCx*cosPhi2 - BCy*sinPhi2)*L2/BC;
	//y = Cy + (BCx*sinPhi2 + BCy*cosPhi2)*L2/BC;
	*/
}

void settleBifurcation(int i, vector<Bifurcation3D> &b, double &x, double &y, double &z, double thresh = 1.0e-6, bool exclaim = false){
	if(b[i].childIndex.size() != 2){
		cout << "\n Error settleBifurcation()!: expected exactly 2 children" << endl;
		x = b[i].position[0];
		y = b[i].position[1];
		z = b[i].position[2];
		return;
	}
	if(b[i].parentIndex < 0){
		x = (b[b[i].childIndex[0]].position[0] + b[b[i].childIndex[1]].position[0])/2.0;
		y = (b[b[i].childIndex[0]].position[1] + b[b[i].childIndex[1]].position[1])/2.0;
		z = (b[b[i].childIndex[0]].position[2] + b[b[i].childIndex[1]].position[2])/2.0;
		return;
	}
	double critAngle(2.0*acos(-1.0)/3.0);
	if(angleBetweenGrandchildren3D(b[i].parentIndex, b[i].childIndex[0], b[i].childIndex[1], b) >= critAngle){
		x = b[b[i].parentIndex].position[0];
		y = b[b[i].parentIndex].position[1];
		z = b[b[i].parentIndex].position[2];
		return;
	}else if(angleBetweenGrandchildren3D(b[i].childIndex[0], b[i].parentIndex, b[i].childIndex[1], b) >= critAngle){
		x = b[b[i].childIndex[0]].position[0];
		y = b[b[i].childIndex[0]].position[1];
		z = b[b[i].childIndex[0]].position[2];
		return;
	}else if(angleBetweenGrandchildren3D(b[i].childIndex[1], b[i].childIndex[0], b[i].parentIndex, b) >= critAngle){
		x = b[b[i].childIndex[1]].position[0];
		y = b[b[i].childIndex[1]].position[1];
		z = b[b[i].childIndex[1]].position[2];
		return;
	}
	double sideC(separation3D(b[b[i].parentIndex], b[b[i].childIndex[0]])),
		sideA(separation3D(b[b[i].childIndex[1]], b[b[i].childIndex[0]])),
		sideB(separation3D(b[b[i].parentIndex], b[b[i].childIndex[1]]));
	if(exclaim){
		cout << "\n settleBifurcation():\n\t " << i << ": sideA = " << sideA
			<< "\n\t sideB = " << sideB
			<< "\n\t sideC = " << sideC;
		if(sideA > sideB && sideA > sideC)
			cout << "\n\t\t A-B-C=" << sideA - sideB - sideC;
		else if(sideB > sideA && sideB > sideC)
			cout << "\n\t\t B-A-C=" << sideB - sideA - sideC;
		else
			cout << "\n\t\t C-A-B=" << sideC - sideA - sideB;
		cout << endl;
	}
	double sideThresh = thresh/10.0;
	if(sideA <= sideThresh && sideB > sideThresh && sideC > sideThresh){ // children overlap
		x = b[b[i].childIndex[0]].position[0];
		y = b[b[i].childIndex[0]].position[1];
		z = b[b[i].childIndex[0]].position[2];
		return;
	}
	if(sideB <= sideThresh || sideC <= sideThresh){ // at least one child is at the parent
		x = b[b[i].parentIndex].position[0];
		y = b[b[i].parentIndex].position[1];
		z = b[b[i].parentIndex].position[2];
		return;
	}
	if(abs(sideA - sideB - sideC) < sideThresh || abs(sideB - sideA - sideC) < sideThresh || abs(sideC - sideB - sideA) < sideThresh){ //zero area
		if(sideA >= sideB && sideA >= sideC){
			x = b[b[i].parentIndex].position[0];
			y = b[b[i].parentIndex].position[1];
			z = b[b[i].parentIndex].position[2];
			if(exclaim)
				cout << "\n sending " << i << " to parent" << endl;
		}else if(sideB >= sideA && sideB >= sideC){
			x = b[b[i].childIndex[0]].position[0];
			y = b[b[i].childIndex[0]].position[1];
			z = b[b[i].childIndex[0]].position[2];
			if(exclaim)
				cout << "\n sending " << i << " to child 0" << endl;
		}else{
			x = b[b[i].childIndex[1]].position[0];
			y = b[b[i].childIndex[1]].position[1];
			z = b[b[i].childIndex[1]].position[2];
			if(exclaim)
				cout << "\n sending " << i << " to child 1" << endl;
		}
		
		return;
	}
	double p(barycentricFermat(sideA, sideB, sideC)), q(barycentricFermat(sideB, sideC, sideA)), r(barycentricFermat(sideC, sideA, sideB));
	double sum(p + q + r);
	x = (p*b[b[i].parentIndex].position[0] + q*b[b[i].childIndex[0]].position[0] + r*b[b[i].childIndex[1]].position[0])/sum;
	y = (p*b[b[i].parentIndex].position[1] + q*b[b[i].childIndex[0]].position[1] + r*b[b[i].childIndex[1]].position[1])/sum;
	z = (p*b[b[i].parentIndex].position[2] + q*b[b[i].childIndex[0]].position[2] + r*b[b[i].childIndex[1]].position[2])/sum;
}

inline double det(	double r1c1, double r1c2,
					double r2c1, double r2c2	){
	return r1c1*r2c2 - r2c1*r1c2;
}

inline double det(	double r1c1, double r1c2, double r1c3,
					double r2c1, double r2c2, double r2c3,
					double r3c1, double r3c2, double r3c3	){
	return r1c1*det(r2c2, r2c3, r3c2, r3c3) - r1c2*det(r2c1, r2c3, r3c1, r3c3) + r1c3*det(r2c1, r2c2, r3c1, r3c2);
}

void circumscribingCircle(double x1, double y1, double x2, double y2, double x3, double y3, double &centerX, double &centerY, double &radius){
	double d1sq(x1*x1 + y1*y1), d2sq(x2*x2 + y2*y2), d3sq(x3*x3 + y3*y3);
	double a(det(x1, y1, 1.0, x2, y2, 1.0, x3, y3, 1.0)),
		bx(-det(d1sq, y1, 1.0, d2sq, y2, 1.0, d3sq, y3, 1.0)),
		by(det(d1sq, x1, 1.0, d2sq, x2, 1.0, d3sq, x3, 1.0)),
		c(det(d1sq, x1, y1, d2sq, x2, y2, d3sq, x3, y3));
	centerX = -0.5*bx/a;
	centerY = -0.5*by/a;
	radius = separation2D(centerX, centerY, x1, y1);//amiss: abs(0.5*sqrt(bx*bx + by*by - 4.0*a*c)/a);
	//cout << "\n circumscribingCircle(): radius is " << radius
	//	<< "\n point1 is " << separation2D(centerX, centerY, x1, y1) << " away from the center"
	//	<< "\n point2 is " << separation2D(centerX, centerY, x2, y2) << " away from the center"
	//	<< "\n point3 is " << separation2D(centerX, centerY, x3, y3) << " away from the center"
	//	<< endl;
}

double quadraticFormula(double a, double b, double c, bool takePositiveRoot = true){
	if(takePositiveRoot)
		return (-b + sqrt(b*b - 4.0*a*c))/(2.0*a);
	return (-b - sqrt(b*b - 4.0*a*c))/(2.0*a);
}

void circlesIntersectButNotNear(double x1, double y1, double r1, double x2, double y2, double r2, double xNot, double yNot, double &xIs, double &yIs){
	double deltaX(x2 - x1), deltaY(y2 - y1), deltaXSq(x2*x2 - x1*x1), deltaYSq(y2*y2 - y1*y1), deltaRSq(r2*r2 - r1*r1);
	//cout << "\n circlesIntersectButNotNear() 1:\n\t " << r1 << " at (" << x1 << ", " << y1
	//	<< "\n\t " << r2 << " at (" << x2 << ", " << y2
	//	<< "\n\t   NOT (" << xNot << ", " << yNot << ")"
	//	<< "\n\tdeltaX = " << deltaX << "\t\t deltaY = " << deltaY << endl;
	if(deltaX != 0.0){
		double ya(-deltaY/deltaX), yb(0.5*(deltaXSq + deltaYSq - deltaRSq)/deltaX);
		double a(1.0 + ya*ya), b(2.0*(ya*(yb - x1) - y1)), c((yb - x1)*(yb - x1) + y1*y1 - r1*r1);
		//cout << "\n circlesIntersectButNotNear() 2: ya = " << ya << ", yb = " << yb << " a = " << a << ", b = " << b << ", c = " << c << endl;
		double yAlpha(quadraticFormula(a, b, c)), yBeta(quadraticFormula(a, b, c, false));
		double xAlpha(ya*yAlpha + yb), xBeta(ya*yBeta + yb);
		//cout << "\n circlesIntersectButNotNear() 3: Alpha = (" << xAlpha << ", " << yAlpha << "), Beta = (" << xBeta << ", " << yBeta << ")" << endl;
		if(separation2D(xNot, yNot, xAlpha, yAlpha) > separation2D(xNot, yNot, xBeta, yBeta)){
			xIs = xAlpha;
			yIs = yAlpha;
		}else{
			xIs = xBeta;
			yIs = yBeta;
		}
	}else if(deltaY != 0.0){
		double xa(-deltaX/deltaY), xb(0.5*(deltaXSq + deltaYSq - deltaRSq)/deltaY);
		double a(1.0 + xa*xa), b(2.0*(xa*(xb - y1) - x1)), c((xb - y1)*(xb - y1) + x1*x1 - r1*r1);
		double xAlpha(quadraticFormula(a, b, c)), xBeta(quadraticFormula(a, b, c, false));
		double yAlpha(xa*xAlpha + xb), yBeta(xa*xBeta + xb);
		if(separation2D(xNot, yNot, xAlpha, yAlpha) > separation2D(xNot, yNot, xBeta, yBeta)){
			xIs = xAlpha;
			yIs = yAlpha;
		}else{
			xIs = xBeta;
			yIs = yBeta;
		}
	}else{
	// both should never be zero, but if they are, set as centers
		xIs = x1;
		yIs = y1;
	}
}

double angleFromSides(double oppositeSide, double adjacentSide1, double adjacentSide2){
	double a(oppositeSide), b(adjacentSide1), c(adjacentSide2);
	return abs(acos((b*b + c*c - a*a)/(2.0*b*c)));
}

void weightedFermat2D(double aw, double ax, double ay, double bw, double bx, double by, double cw, double cx, double cy, double &x, double &y){
	//cout << "\n weightedFermat2D(): " << aw << " for (" << ax << ", " << ay << ")"
	//	<< "\t\t\t\t" << bw << " for (" << bx << ", " << by << ")"
	//	<< "\t\t\t\t" << cw << " for (" << cx << ", " << cy << ")" << endl;
	if(aw >= bw + cw){
		x = ax;
		y = ay;
		return;
	}else if(bw >= aw + cw){
		x = bx;
		y = by;
		return;
	}else if(cw >= aw + bw){
		x = bx;
		y = by;
		return;
	}

	double P1P2(separation2D(ax, ay, bx, by));
	double P1P3(separation2D(ax, ay, cx, cy));
	double P2P3(separation2D(bx, by, cx, cy));

	double alpha1(angleFromSides(P2P3, P1P2, P1P3)), // angle between P1P2 and P1P3; opposite side is P2P3
		alpha2(angleFromSides(P1P3, P1P2, P2P3)), // angle between P1P2 and P2P3; opposite side is P1P3
		alpha3(angleFromSides(P1P2, P1P3, P2P3)), // angle between P1P3 and P2P3; opposite side is P1P2
		theta1(angleFromSides(aw, bw, cw)), // angle between lambda2 and lambda3
		theta2(angleFromSides(bw, aw, cw)), // angle between lambda1 and lambda3
		theta3(angleFromSides(cw, aw, bw)); // angle between lambda1 and lambda2

	if(alpha1 >= acos(-1.0) - theta1){
		x = ax;
		y = ay;
		return;
	}else if(alpha2 >= acos(-1.0) - theta2){
		x = bx;
		y = by;
		return;
	}else if(alpha3 >= acos(-1.0) - theta3){
		x = cx;
		y = cy;
		return;
	}

	double x1(0.0), y1(0.0), r1(0.0), x2(0.0), y2(0.0), r2(0.0), wx(0.0), wy(0.0);
	
	// P1, P2, and w is P3' that is the farther from P3
	double P2P3p(P1P2*aw/cw), P1P3p(P1P2*bw/cw);
	circlesIntersectButNotNear(ax, ay, P1P3p, bx, by, P2P3p, cx, cy, wx, wy);
	//cout << "weightedFermat2D(): (wx, wy)_1 = (" << wx << ", " << wy << ")" << endl;
	circumscribingCircle(ax, ay, bx, by, wx, wy, x1, y1, r1);

	// P1, P3, and w is P2' that is the farther from P2
	
	double P3P2p(P1P3*aw/bw), P1P2p(P1P3*cw/bw);
	circlesIntersectButNotNear(ax, ay, P1P2p, cx, cy, P3P2p, bx, by, wx, wy);
	//cout << "weightedFermat2D(): (wx, wy)_2 = (" << wx << ", " << wy << ")" << endl;
	circumscribingCircle(ax, ay, cx, cy, wx, wy, x2, y2, r2);

	circlesIntersectButNotNear(x1, y1, r1, x2, y2, r2, ax, ay, x, y);
}

void weightedFermat2DByVolume(double pRad, double px, double py, double aRad, double ax, double ay, double bRad, double bx, double by, double &x, double &y){
	double pw(acos(-1.0)*pRad*pRad), aw(acos(-1.0)*aRad*aRad), bw(acos(-1.0)*bRad*bRad);
	//cout << "\n weightedFermat2DByVolume(): pRad = " << pRad << endl;
	weightedFermat2D(pw, px, py, aw, ax, ay, bw, bx, by, x, y);
}

// "i" is the upstream Bifurcation2D point;
// "it" is the downstream branchPoint2D index
double largestNondegenerateDownstreamRadius(const vector<Bifurcation2D> &b, int i, terseSeg **bt, int it,
	const vector<BranchPoint2D> &bc, int *terseIndex, double consolidateThresh, double capRad){
	if(separation2D(b[i], b[b[i].parentIndex]) <= consolidateThresh){
		//look downstream -- consider all adjacent nondegenerate points
		if(b[i].childIndex.size() < 1)
			return capRad;
		double tempRad(bt[b[i].childIndex[0]]->rad);
		for(unsigned int c(1); c < b[i].childIndex.size(); c++){
			if(tempRad < bt[b[i].childIndex[c]]->rad)
				tempRad = bt[b[i].childIndex[c]]->rad;
		}
		return tempRad;
	}else // the child bifurcation is unique, the radius is already defined
		return bt[it]->rad;
}

// degenerate points are already consolidated; segment i has the downstream end as the Bifurcation2D of interest
// if degenerate with heart, takes largest child
double largestNondegenerateUpstreamRadius(const vector<Bifurcation2D> &b, int i, terseSeg **bt, int it,
	const vector<BranchPoint2D> &bc, int *terseIndex, double consolidateThresh, double capRad){
	if(bc[it].parentIndex == 0 && separation2D(bc[it], bc[bc[it].parentIndex]) <= consolidateThresh) // degenerate with heart
		return largestNondegenerateDownstreamRadius(b, 0, bt, it, bc, terseIndex, consolidateThresh, capRad);
	return bt[it]->rad;
}

int findCorrespondingTerseIndex(const vector<Bifurcation2D> &b, int i, const vector<BranchPoint2D> &bc, double consolidateThresh){
	int ic(0);
	double nearestSep(separation2D(b[i].position[0], b[i].position[1], bc[0].position[0], bc[0].position[1]));
	for(unsigned int j(1); j < bc.size(); j++){
		double sep(separation2D(b[i].position[0], b[i].position[1], bc[j].position[0], bc[j].position[1]));
		if(nearestSep + 0.0*consolidateThresh > sep){
			ic = j;
			nearestSep = sep;
		}
	}
	return ic;
}


// x is the initial guess
double MurrayExponent(double pRad, vector<double> cRad, double thresh = 1.0e-3, double x = 2.0){
	//cout << "\nMurrayExponent(): pRad = " << pRad << ", cRad = " << makeVectorString(cRad) << endl;
	double eps(2.0*thresh);
	while(abs(eps) > thresh){
		double f(pow(pRad, x)), df(x*pow(pRad, x - 1.0));
		for(unsigned int c(0); c < cRad.size(); c++){
			f -= pow(cRad[c], x);
			df -= x*pow(cRad[c], c - 1.0);
		}
		eps = f/df;
		x -= eps;
		//cout << "\n\t\t MurrayExponent(): eps = " << eps << ", x = " << x << endl;
	}
	//cout << "\n MurrayExponent(): x = " << x << endl;
	return x;
}

// returns index of element in b (NOT b[i].childIndex) that is not cNot (which is an index specifying an element of b, not b[i].childIndex)
int otherChildIndex(const vector<Bifurcation2D> &b, int i, int cNot){
	if(b[i].childIndex.size() < 2){
		cout << "\n Error otherChildIndex(): b[i = " << i << "] has ";
		if(b[i].childIndex.size() == 1)
			cout << "1 child." << endl;
		else
			cout << "no children." << endl;
		abort();
	}
	if(b[i].childIndex.size() > 2){
		cout << "\n Warning otherChildIndex(): b[i = " << i << "] has " << b[i].childIndex.size() << " children." << endl;
	}
	if(b[i].childIndex[0] == cNot)
		return b[i].childIndex[1];
	return b[i].childIndex[0];
}

/*double givenOrEffectiveBifurcationRadiusDown(int i, const vector<Bifurcation2D> &b, int *terseIndex, terseSeg **bt, const vector<BranchPoint2D> bc, double exponent, double heartRad, double capRad);

// assumes heart is b[0]
double givenOrEffectiveBifurcationRadiusUp(int i, const vector<Bifurcation2D> &b, int *terseIndex, terseSeg **bt, const vector<BranchPoint2D> bc, double exponent, double heartRad, double capRad){
	if(b[i].parentIndex < 1) // heart or above
		return heartRad;
	if(b[i].childIndex.size() < 2) // this should not happen going up, but I'll keep this here for anomalous cases
		return capRad;
	if(terseIndex[i] != terseIndex[b[i].parentIndex]) // non degenerate is easy
		return bt[terseIndex[i]]->rad;
	// segment is degenerate -- examine local context without backtracking down the tree
	double pRad(givenOrEffectiveBifurcationRadiusUp(b[i].parentIndex, b, terseIndex, bt, bc, exponent, heartRad, capRad));
	double cRad(givenOrEffectiveBifurcationRadiusDown(otherChildIndex(b, b[i].parentIndex, i), b, terseIndex, bt, bc, exponent, heartRad, capRad));
	return pow(pow(pRad, exponent) - pow(cRad, exponent), 1.0/exponent);
}*/

//assumes b[0] has exactly one child
terseSeg* getRoot(const vector<Bifurcation2D> &b, terseSeg** linList){
	return linList[b[0].childIndex[0]];
}


//assumes b[0] has exactly one child
terseSeg* getRoot(const vector<BranchPoint2D> &b, terseSeg** linList){
	return linList[b[0].childIndex[0]];
}

// assumes heart is b[0]
double givenOrEffectiveBifurcationRadiusDown(int i, const vector<Bifurcation2D> &b, int *terseIndex, terseSeg **bt, const vector<BranchPoint2D> &bc, double exponent, double heartRad, double capRad){
	if(b[i].parentIndex < 0) // heart or above; this should not happen going down, but I'll keep this here for anomalous cases
		return heartRad;
	if(b[i].childIndex.size() < 2) 
		return capRad;
	if(terseIndex[i] == 0)
		return heartRad;
	//cout << "\n givenOrEffectiveBifurcationRadiusDown(): i = " << i << endl;
	//cout << "\t\t bc.size() = " << bc.size() << endl;
	//cout << "\t\t terseIndex[i] = " << terseIndex[i] << endl;
	//cout << "\t\t b[i].parentIndex = " << b[i].parentIndex << endl;
	//cout << "\t\t terseIndex[b[i].parentIndex] = " << terseIndex[b[i].parentIndex] << endl;
	//cout << "\t bt is:" << makeStringTreeRadOnly(getRoot(bc, bt)) << endl;
 	if(terseIndex[i] != terseIndex[b[i].parentIndex])
		return bt[terseIndex[i]]->rad;
	//cout << "\n givenOrEffectiveBifurcationRadiusDown(): seg is degen" << endl;
	// segment is degenerate -- examine local context withouth backtracking up the tree
	double c0Rad(givenOrEffectiveBifurcationRadiusDown(b[i].childIndex[0], b, terseIndex, bt, bc, exponent, heartRad, capRad));
	double c1Rad(givenOrEffectiveBifurcationRadiusDown(b[i].childIndex[1], b, terseIndex, bt, bc, exponent, heartRad, capRad));
	//cout << "\n givenOrEffectiveBifurcationRadiusDown(): c0Rad = " << c0Rad
	//	<< ", c1Rad = " << c1Rad << ", exponent = " << exponent << endl;
	return pow(pow(c0Rad, exponent) + pow(c1Rad, exponent), 1.0/exponent);
}

double effectiveMurrayRadiusForBifurcation(int i, const vector<Bifurcation2D> &b, int *terseIndex, terseSeg **bt, const vector<BranchPoint2D> bc, double heartRad, double capRad){
	double pRad(heartRad);
	//cout << "\n effectiveMurrayRadiusForBifurcation(): i = " << i
	//	<< " \t terseIndex[i] = " << terseIndex[i]
	//	<< endl;
	if(terseIndex[i] > 0)
		pRad = bt[terseIndex[i]]->rad;
	vector<double> cRad;
	for(unsigned int c(0); c < bc[terseIndex[i]].childIndex.size(); c++)
		cRad.push_back(bt[bc[terseIndex[i]].childIndex[c]]->rad);
	//cout << "\n effectiveMurrayRadiusForBifurcation(): pRad = " << pRad << ", cRad = " << makeVectorString(cRad) << endl;
	double effectiveExponent(MurrayExponent(pRad, cRad));
	//pRad = givenOrEffectiveBifurcationRadiusUp(i, b, terseIndex, bt, bc, effectiveExponent, heartRad, capRad);
	return givenOrEffectiveBifurcationRadiusDown(i, b, terseIndex, bt, bc, effectiveExponent, heartRad, capRad);
}

double givenOrMurrayRadius(int i, const vector<Bifurcation2D> &b, int *terseIndex, terseSeg **bt, const vector<BranchPoint2D> bc, double heartRad, double capRad){
	if(terseIndex[i] != terseIndex[b[i].parentIndex]){
		//cout << "\n givenOrMurrayRadius(): terseIndex[i = " << i << "] = " << terseIndex[i]
		//	<< "\n\t with bc.size() = " << bc.size()
		//	<< "\n\t and terseIndex[b[i].parentIndex = " << b[i].parentIndex << "] = " << terseIndex[b[i].parentIndex]
		//	<< endl;
		//cout << "\n givenOrMurrayRadius(): bt is currently:\n";
		//for(unsigned int j(1); j < bc.size(); j++)
		//	cout << "\t\t" << j << "\t" << bt[j]->rad << endl;
		//cout << "\t\t i = " << i << "\t terseIndex[i] = " << terseIndex[i] << endl;
		if(terseIndex[i] > 0 && bc[terseIndex[i]].parentIndex > 0){
			//cout << "\n givenOrMurrayRadius(): bt[terseIndex[i]]->rad = " << bt[terseIndex[i]]->rad << endl;
			return bt[terseIndex[i]]->rad;
		}
		return heartRad;
	}
	//cout << "\n givenOrMurrayRadius(): finding effective radius..." << endl;
	return effectiveMurrayRadiusForBifurcation(i, b, terseIndex, bt, bc, heartRad, capRad);
}

// index of seg corresponds to index of the BranchPoint2D at the downstream end
void createSegAndAddChildren(const vector<BranchPoint2D> &b, terseSeg** linList, int index, double rad, double betaRad, double capRad){
	//cout << "\n createSegAndAddChildren(): creating new terseSeg at index " << index << endl;
	linList[index] = new terseSeg;
	linList[index]->len = separation2D(b[index], b[b[index].parentIndex]);
	if(b[index].childIndex.size() == 0)
		linList[index]->rad = capRad;
	else
		linList[index]->rad = rad;
	linList[index]->flow = heartFlowRate;
	if(b[index].parentIndex == 0)
		linList[index]->parent = NULL;
	else
		linList[index]->parent = linList[b[index].parentIndex];
	linList[index]->name = makeString(index);
	linList[index]->parentName = makeString(b[index].parentIndex);
	linList[index]->nchild = b[index].childIndex.size();
	linList[index]->child = new terseSeg*[b[index].childIndex.size()];
	for(unsigned int c(0); c < b[index].childIndex.size(); c++){
		createSegAndAddChildren(b, linList, b[index].childIndex[c], betaRad*rad, betaRad, capRad);
		linList[index]->child[c] = linList[b[index].childIndex[c]];
	}
}

//void deleteTerseSegTree(terseSeg **bt){
//	delete[] bt;
//}

terseSeg** createTerseSegTree(int maxSize, vector<BranchPoint2D> &b, double rootRad = 0.01, double betaRad = 0.95, double capRad = 0.005, double thresh = 1.0e-6, string initialDrawFn = "nowhere"){
	if(b[0].childIndex.size() != 1){
		cout << "\n Error createTerseSegTree(): b[0] must have exactly one child." << endl;
		return NULL;
	}
	terseSeg** linList = new terseSeg*[maxSize];
	createSegAndAddChildren(b, linList, b[0].childIndex[0], rootRad, betaRad, capRad);
	terseSeg* root(getRoot(b, linList));
	organizeBinary(root);
	updateBinary(root);
	//cout << "\n createTerseSegTree(): initial tree before optimizing is\n" << makeStringTreeRadOnly(root) << endl;
	int width(19), height(11), border(2),//int width(1920), height(1080), border(25),
		printEveryOptimization(0), maxOptimizationIts(10000);
	terseEnds *rootEnds = drawBinaryTreeByFlow(initialDrawFn, root, width, height, border);
	naivelyOptimizeWithConstantTips(rootEnds, thresh, maxOptimizationIts, capRad);
	//cout << "\n createTerseSegTree(): tree after optimizing is\n" << makeStringTreeRadOnly(root) << endl;
	return linList;
}

void recreateTerseSegTree(terseSeg **linList, const vector<BranchPoint2D> &b, double rootRad = 0.01, double betaRad = 0.95, double capRad = 0.005, double thresh = 1.0e-6, string initialDrawFn = "nowhere"){
	if(b[0].childIndex.size() != 1){
		cout << "\n Error recreateTerseSegTree(): b[0] must have exactly one child." << endl;
		return;
	}
	createSegAndAddChildren(b, linList, b[0].childIndex[0], rootRad, betaRad, capRad);
	terseSeg* root(getRoot(b, linList));
	organizeBinary(root);
	updateBinary(root);
	//cout << "\n createTerseSegTree(): initial tree is\n" << makeStringTree(root) << endl;
	int width(19), height(11), border(2),//int width(1920), height(1080), border(25),
		printEveryOptimization(0), maxOptimizationIts(10000);
	terseEnds *rootEnds = drawBinaryTreeByFlow(initialDrawFn, root, width, height, border);
	naivelyOptimizeWithConstantTips(rootEnds, thresh, maxOptimizationIts, capRad);
	//naivelyOptimizeWithConstantTipsAndRoot(rootEnds, thresh, maxOptimizationIts, printEveryOptimization, "nowhere", width, height, border);
}

// if degenerate, assume radius is given by largest of non-degenerate adjacent segments along same direction (down from children, up from parent)
double settleBifurcationWithRadius(int i, vector<Bifurcation2D> &b, terseSeg **bt, vector<BranchPoint2D> &bc,
	double consolidateThresh, int *terseIndex, double rootRad, double betaRad, double capRad, double radThresh){
	if(b[i].parentIndex < 0 || b[i].childIndex.size() < 1) // skip heart or tips
		return 0.0;
	//settle based on local info of bifurcation using weighted fermat
	int it(terseIndex[i]);
	double pRad(givenOrMurrayRadius(i, b, terseIndex, bt, bc, rootRad, capRad)),//largestNondegenerateUpstreamRadius(b, i, bt, it, bc, terseIndex, consolidateThresh, capRad)),
		c0Rad(givenOrMurrayRadius(b[i].childIndex[0], b, terseIndex, bt, bc, rootRad, capRad)),//largestNondegenerateDownstreamRadius(b, i, bt, bc[it].childIndex[0], bc, terseIndex, consolidateThresh, capRad)),
		c1Rad(givenOrMurrayRadius(b[i].childIndex[1], b, terseIndex, bt, bc, rootRad, capRad)),//largestNondegenerateDownstreamRadius(b, i, bt, bc[it].childIndex[1], bc, terseIndex, consolidateThresh, capRad)),
		preX(b[i].position[0]),
		preY(b[i].position[1]);
	//double nothing(givenOrMurrayRadius(i, b, terseIndex, bt, bc, rootRad, capRad));
	//cout << "\nsettleBifurcationWithRadius(): pRad = " << pRad
	//	<< " \t nothing = " << nothing
		//<< " found from method that returns " << givenOrMurrayRadius(i, b, terseIndex, bt, bc, rootRad, capRad)
	//	<< endl;

	//find weighted fermat point
	Bifurcation2D *par(&b[b[i].parentIndex]),
		*child0(&b[b[i].childIndex[0]]),
		*child1(&b[b[i].childIndex[1]]);
	weightedFermat2DByVolume(pRad, par->position[0], par->position[1],
		c0Rad, child0->position[0], child0->position[1],
		c1Rad, child1->position[0], child1->position[1],
		b[i].position[0], b[i].position[1]);

	bc = consolidateDegenerates(b, consolidateThresh);
	//deleteTerseSegTree(bt);
	recreateTerseSegTree(bt, bc, rootRad, betaRad, capRad, radThresh);
	for(unsigned int j(0); j < b.size(); j++)
		terseIndex[j] = findCorrespondingTerseIndex(b, j, bc, consolidateThresh);
	//cout << "\n settleBifurcationWithRadius(): updated bc to:\n" << makeStringTree(bc) << endl;
	//cout << "\n settleBifurcationWithRadius(): updated bt (bc.size() = " << bc.size() << ") to:"
	//	<< makeStringTreeRadOnly(getRoot(bc, bt)) << endl;
	//for(unsigned int j(1); j < bc.size(); j++)
	//	cout << j << "\t" << bt[j]->rad << endl;

	//cout << "\n settleBifurcationWithRadius(): i = " << i << "\t (" << preX << ", " << preY
	//	<< ") moved to (" << b[i].position[0] << ", " << b[i].position[1] << ")" << endl;
	return separation2D(b[i].position[0], b[i].position[1], preX, preY);
}

unsigned int globalCount(0);

string makeStringTreeOneLine(vector<Bifurcation2D> &b, int rootIndex = 0){
	if(b.size() < 1)
		return "";
	if(b[rootIndex].childIndex.size() == 0)
		return makeString(rootIndex);
	string stringTree("");
	for(unsigned int c(0); c < b[rootIndex].childIndex.size(); c++)
		stringTree += "(" + makeStringTreeOneLine(b, b[rootIndex].childIndex[c]) + ")";
	return stringTree;
}

string makeStringTreeOneLine(vector<BranchPoint2D> &b, int rootIndex = 0){
	if(b.size() < 1)
		return "";
	if(b[rootIndex].childIndex.size() == 0)
		return makeString(rootIndex);
	string stringTree("");
	for(unsigned int c(0); c < b[rootIndex].childIndex.size(); c++)
		stringTree += "(" + makeStringTreeOneLine(b, b[rootIndex].childIndex[c]) + ")";
	return stringTree;
}

string makeStringTreeOneLine(vector<Bifurcation3D> &b, int rootIndex = 0){
	if(b.size() < 1)
		return "";
	if(b[rootIndex].childIndex.size() == 0)
		return makeString(rootIndex);
	string stringTree("");
	for(unsigned int c(0); c < b[rootIndex].childIndex.size(); c++)
		stringTree += "(" + makeStringTreeOneLine(b, b[rootIndex].childIndex[c]) + ")";
	return stringTree;
}

string makeStringTreeOneLine(vector<BranchPoint3D> &b, int rootIndex = 0){
	if(b.size() < 1)
		return "";
	if(b[rootIndex].childIndex.size() == 0)
		return makeString(rootIndex);
	string stringTree("");
	for(unsigned int c(0); c < b[rootIndex].childIndex.size(); c++)
		stringTree += "(" + makeStringTreeOneLine(b, b[rootIndex].childIndex[c]) + ")";
	return stringTree;
}

Bifurcation2D RtildeTorres(const vector<Bifurcation2D> &b, vector<int> constantIndices, vector<double> w, int movingIndex){
	Bifurcation2D x;
	x.position[0] = x.position[1] = 0.0;
	for(unsigned int i(0); i < constantIndices.size(); i++){
		double sep(separation2D(b[movingIndex], b[constantIndices[i]]));
		if(sep > 0.0){
			x.position[0] += w[i]*(b[constantIndices[i]].position[0] - b[movingIndex].position[0])/sep;
			x.position[1] += w[i]*(b[constantIndices[i]].position[1] - b[movingIndex].position[1])/sep;
		}
	}
	return x;
}

Bifurcation3D RtildeTorres(const vector<Bifurcation3D> &b, vector<int> constantIndices, vector<double> w, int movingIndex){
	Bifurcation3D x;
	x.position[0] = x.position[1] = 0.0;
	for(unsigned int i(0); i < constantIndices.size(); i++){
		double sep(separation3D(b[movingIndex], b[constantIndices[i]]));
		if(sep > 0.0){
			x.position[0] += w[i]*(b[constantIndices[i]].position[0] - b[movingIndex].position[0])/sep;
			x.position[1] += w[i]*(b[constantIndices[i]].position[1] - b[movingIndex].position[1])/sep;
			x.position[2] += w[i]*(b[constantIndices[i]].position[2] - b[movingIndex].position[2])/sep;
		}
	}
	return x;
}

double betaTorres(const vector<Bifurcation2D> &b, vector<int> constantIndices, vector<double> w, int movingIndex){
	double gammaTorres(0.0);
	for(unsigned int i(0); i < constantIndices.size(); i++){
		Bifurcation2D Rtilde(RtildeTorres(b, constantIndices, w, movingIndex));
		double r(separation2D(Rtilde.position[0], Rtilde.position[1], 0.0, 0.0));
		if(separation2D(b[movingIndex], b[constantIndices[i]]) == 0.0 && r != 0.0)
			gammaTorres = w[i]/r;
	}
	if(gammaTorres > 1.0)
		return 1.0;
	return gammaTorres;
}

double betaTorres(const vector<Bifurcation3D> &b, vector<int> constantIndices, vector<double> w, int movingIndex){
	double gammaTorres(0.0);
	for(unsigned int i(0); i < constantIndices.size(); i++){
		Bifurcation3D Rtilde(RtildeTorres(b, constantIndices, w, movingIndex));
		double r(separation2D(Rtilde.position[0], Rtilde.position[1], 0.0, 0.0));
		if(separation3D(b[movingIndex], b[constantIndices[i]]) == 0.0 && r != 0.0)
			gammaTorres = w[i]/r;
	}
	if(gammaTorres > 1.0)
		return 1.0;
	return gammaTorres;
}

Bifurcation2D TtildeTorres(const vector<Bifurcation2D> &b, vector<int> constantIndices, vector<double> w, int movingIndex){
	//int constantIndicesIndex(-1);
	//for(unsigned int i(0); i < constantIndices.size(); i++){
	//	if(separation2D(b[movingIndex], b[constantIndices[i]]) == 0.0)
	//		constantIndicesIndex = i;
	//}
	double ATorres(0.0);
	for(unsigned int i(0); i < constantIndices.size(); i++){
		double sep(separation2D(b[movingIndex], b[constantIndices[i]]));
		if(sep > 0.0)
			ATorres += 0.5*w[i]/sep;
	}
	Bifurcation2D x;
	x.position[0] = x.position[1] = 0.0;
	for(unsigned int i(0); i < constantIndices.size(); i++){
		double sep(separation2D(b[movingIndex], b[constantIndices[i]]));
		if(sep > 0.0){
			x.position[0] += w[i]*b[constantIndices[i]].position[0]/sep;
			x.position[1] += w[i]*b[constantIndices[i]].position[1]/sep;
		}
	}
	x.position[0] /= 2.0*ATorres;
	x.position[1] /= 2.0*ATorres;
	return x;
}

Bifurcation3D TtildeTorres(const vector<Bifurcation3D> &b, vector<int> constantIndices, vector<double> w, int movingIndex){
	double ATorres(0.0);
	for(unsigned int i(0); i < constantIndices.size(); i++){
		double sep(separation3D(b[movingIndex], b[constantIndices[i]]));
		if(sep > 0.0)
			ATorres += 0.5*w[i]/sep;
	}
	Bifurcation3D x;
	x.position[0] = x.position[1] = 0.0;
	for(unsigned int i(0); i < constantIndices.size(); i++){
		double sep(separation3D(b[movingIndex], b[constantIndices[i]]));
		if(sep > 0.0){
			x.position[0] += w[i]*b[constantIndices[i]].position[0]/sep;
			x.position[1] += w[i]*b[constantIndices[i]].position[1]/sep;
			x.position[2] += w[i]*b[constantIndices[i]].position[2]/sep;
		}
	}
	x.position[0] /= 2.0*ATorres;
	x.position[1] /= 2.0*ATorres;
	x.position[2] /= 2.0*ATorres;
	return x;
}

Bifurcation2D TTorres(const vector<Bifurcation2D> &b, vector<int> constantIndices, vector<double> w, int movingIndex){
	Bifurcation2D Ttilde(TtildeTorres(b, constantIndices, w, movingIndex));
	double beta(betaTorres(b, constantIndices, w, movingIndex));
	Ttilde.position[0] *= 1.0 - beta;
	Ttilde.position[1] *= 1.0 - beta;
	Ttilde.position[0] += beta*b[movingIndex].position[0];
	Ttilde.position[1] += beta*b[movingIndex].position[1];
	return Ttilde;
}

Bifurcation3D TTorres(const vector<Bifurcation3D> &b, vector<int> constantIndices, vector<double> w, int movingIndex){
	Bifurcation3D Ttilde(TtildeTorres(b, constantIndices, w, movingIndex));
	double beta(betaTorres(b, constantIndices, w, movingIndex));
	Ttilde.position[0] *= 1.0 - beta;
	Ttilde.position[1] *= 1.0 - beta;
	Ttilde.position[2] *= 1.0 - beta;
	Ttilde.position[0] += beta*b[movingIndex].position[0];
	Ttilde.position[1] += beta*b[movingIndex].position[1];
	Ttilde.position[2] += beta*b[movingIndex].position[2];
	return Ttilde;
}

Bifurcation2D QTorres(const vector<Bifurcation2D> &b, vector<int> constantIndices, vector<double> w, int movingIndex){
	//for(unsigned int i(0); i < constantIndices.size(); i++){
	//	if(separation2D(b[movingIndex], b[constantIndices[i]]) == 0.0){
	//		Bifurcation2D T(TTorres(b, constantIndices, w, constantIndices[i]));
	//		return T;
	//	}
	//}
	return TTorres(b, constantIndices, w, movingIndex);
}

Bifurcation3D QTorres3D(const vector<Bifurcation3D> &b, vector<int> constantIndices, vector<double> w, int movingIndex){
	return TTorres(b, constantIndices, w, movingIndex);
}

void TorresWeiszfeld(vector<Bifurcation2D> &b, vector<int> constantIndices, vector<double> w, vector<int> movingIndices, double eps = 1.0e-6){
	double change(2.0*eps);
	int maxIts(1000), it(0);
	//cout << "\nTW from (" << b[movingIndices[0]].position[0] << ", " << b[movingIndices[0]].position[1] << ")" << endl;
	while(it < maxIts && eps < change){
		Bifurcation2D Qx(QTorres(b, constantIndices, w, movingIndices[0]));
		change = separation2D(b[movingIndices[0]], Qx);
		b[movingIndices[0]].position[0] = Qx.position[0];
		b[movingIndices[0]].position[1] = Qx.position[1];
	}
	//cout << "\t to (" << b[movingIndices[0]].position[0] << ", " << b[movingIndices[0]].position[1] << ")" << endl;
	if(it == maxIts)
		cout << "\n Error TorresWeiszfeld(): Did not converge! (Leaving the launder'd out to dry...)" << endl;
	else{
		for(unsigned int i(1); i < movingIndices.size(); i++){
			b[movingIndices[i]].position[0] = b[movingIndices[0]].position[0];
			b[movingIndices[i]].position[1] = b[movingIndices[0]].position[1];
		}
	}
}

void TorresWeiszfeld(vector<Bifurcation3D> &b, vector<int> constantIndices, vector<double> w, vector<int> movingIndices, double eps = 1.0e-6){
	double change(2.0*eps);
	int maxIts(1000), it(0);
	//cout << "\nTW from (" << b[movingIndices[0]].position[0] << ", " << b[movingIndices[0]].position[1] << ")" << endl;
	while(it < maxIts && eps < change){
		Bifurcation3D Qx(QTorres3D(b, constantIndices, w, movingIndices[0]));
		change = separation3D(b[movingIndices[0]], Qx);
		b[movingIndices[0]].position[0] = Qx.position[0];
		b[movingIndices[0]].position[1] = Qx.position[1];
		b[movingIndices[0]].position[2] = Qx.position[2];
	}
	//cout << "\t to (" << b[movingIndices[0]].position[0] << ", " << b[movingIndices[0]].position[1] << ")" << endl;
	if(it == maxIts)
		cout << "\n Error TorresWeiszfeld(): Did not converge! (Leaving the launder'd out to dry...)" << endl;
	else{
		for(unsigned int i(1); i < movingIndices.size(); i++){
			b[movingIndices[i]].position[0] = b[movingIndices[0]].position[0];
			b[movingIndices[i]].position[1] = b[movingIndices[0]].position[1];
			b[movingIndices[i]].position[2] = b[movingIndices[0]].position[2];
		}
	}
}

void findDegenerateGob(const vector<Bifurcation2D> &b, int seed, vector<int> &constantIndices, vector<int> &movingIndices, double thresh = 1.0e-6){
	constantIndices.clear();
	movingIndices.clear();
	movingIndices.push_back(seed);
	queue<int> toCheck;
	toCheck.push(seed);
	while(toCheck.size() > 0){
		int checking(toCheck.front());
		toCheck.pop();
		if(b[checking].parentIndex < 0)
			cout << "\n Error findDegenerateGob(): root has been labeled degenerate!" << endl;
		else{
			if(b[b[checking].parentIndex].parentIndex > -1 && separation2D(b[checking], b[b[checking].parentIndex]) <= thresh){
				if(!containedIn(b[checking].parentIndex, movingIndices)){
					movingIndices.push_back(b[checking].parentIndex);
					toCheck.push(b[checking].parentIndex);
				}
			} else
				constantIndices.push_back(b[checking].parentIndex);
		}
		if(b[checking].childIndex.size() < 1)
			cout << "\n Error findDegenerateGob(): tip has been labeled degenerate!" << endl;
		else{
			for(unsigned int c(0); c < b[checking].childIndex.size(); c++){
				if(b[b[checking].childIndex[c]].childIndex.size() > 0 && separation2D(b[checking], b[b[checking].childIndex[c]]) <= thresh){
					if(!containedIn(b[checking].childIndex[c], movingIndices)){
						movingIndices.push_back(b[checking].childIndex[c]);
						toCheck.push(b[checking].childIndex[c]);
					}
				} else
					constantIndices.push_back(b[checking].childIndex[c]);
			}
		}
	}
}

void findDegenerateGob(const vector<Bifurcation3D> &b, int seed, vector<int> &constantIndices, vector<int> &movingIndices, double thresh = 1.0e-6){
	constantIndices.clear();
	movingIndices.clear();
	movingIndices.push_back(seed);
	queue<int> toCheck;
	toCheck.push(seed);
	while(toCheck.size() > 0){
		int checking(toCheck.front());
		toCheck.pop();
		if(b[checking].parentIndex < 0)
			cout << "\n Error findDegenerateGob(): root has been labeled degenerate!" << endl;
		else{
			if(b[b[checking].parentIndex].parentIndex > -1 && separation3D(b[checking], b[b[checking].parentIndex]) <= thresh){
				if(!containedIn(b[checking].parentIndex, movingIndices)){
					movingIndices.push_back(b[checking].parentIndex);
					toCheck.push(b[checking].parentIndex);
				}
			} else
				constantIndices.push_back(b[checking].parentIndex);
		}
		if(b[checking].childIndex.size() < 1)
			cout << "\n Error findDegenerateGob(): tip has been labeled degenerate!" << endl;
		else{
			for(unsigned int c(0); c < b[checking].childIndex.size(); c++){
				if(b[b[checking].childIndex[c]].childIndex.size() > 0 && separation3D(b[checking], b[b[checking].childIndex[c]]) <= thresh){
					if(!containedIn(b[checking].childIndex[c], movingIndices)){
						movingIndices.push_back(b[checking].childIndex[c]);
						toCheck.push(b[checking].childIndex[c]);
					}
				} else
					constantIndices.push_back(b[checking].childIndex[c]);
			}
		}
	}
}

void findDegenerateGobs(const vector<Bifurcation2D> &b, vector<vector<int> > &boundaries, vector<vector<int> > &gobs, double thresh = 1.0e-6){
	bool *isGob = new bool[b.size()];
	for(unsigned int i(0); i < b.size(); i++)
		isGob[i] = false;
	for(unsigned int i(0); i < b.size(); i++){
		if(isGob[i])
			continue;
		if(b[i].parentIndex > -1 && b[i].childIndex.size() == 2){
			bool isInDegenerateGob(separation2D(b[i], b[b[i].parentIndex]) <= thresh && b[b[i].parentIndex].parentIndex > -1);
			for(unsigned int c(0); c < b[i].childIndex.size(); c++)
				isInDegenerateGob = isInDegenerateGob || (separation2D(b[i], b[b[i].childIndex[c]]) <= thresh && b[b[i].childIndex[c]].childIndex.size() > 0);
			if(isInDegenerateGob){
				vector<int> cInd, mInd;
				findDegenerateGob(b, i, cInd, mInd, thresh);
				gobs.push_back(mInd);
				boundaries.push_back(cInd);
				for(unsigned int j(0); j < mInd.size(); j++)
					isGob[mInd[j]] = true;
			}
		}
	}
	delete[] isGob;
}

void findDegenerateGobs(const vector<Bifurcation3D> &b, vector<vector<int> > &boundaries, vector<vector<int> > &gobs, double thresh = 1.0e-6){
	bool *isGob = new bool[b.size()];
	for(unsigned int i(0); i < b.size(); i++)
		isGob[i] = false;
	for(unsigned int i(0); i < b.size(); i++){
		if(isGob[i])
			continue;
		if(b[i].parentIndex > -1 && b[i].childIndex.size() == 2){
			bool isInDegenerateGob(separation3D(b[i], b[b[i].parentIndex]) <= thresh && b[b[i].parentIndex].parentIndex > -1);
			for(unsigned int c(0); c < b[i].childIndex.size(); c++)
				isInDegenerateGob = isInDegenerateGob || (separation3D(b[i], b[b[i].childIndex[c]]) <= thresh && b[b[i].childIndex[c]].childIndex.size() > 0);
			if(isInDegenerateGob){
				vector<int> cInd, mInd;
				findDegenerateGob(b, i, cInd, mInd, thresh);
				gobs.push_back(mInd);
				boundaries.push_back(cInd);
				for(unsigned int j(0); j < mInd.size(); j++)
					isGob[mInd[j]] = true;
			}
		}
	}
	delete[] isGob;
}

string recordString(const vector<Bifurcation2D> &b){
	string r("");
	for(unsigned int i(0); i < b.size(); i++){
		r += makeString(i) + "\t" + makeString(b[i].position[0]) + "\t" + makeString(b[i].position[1]);
		if(b[i].childIndex.size() > 0)
			r += "\t" + makeString(b[i].childIndex[0]);
		if(b[i].childIndex.size() > 1)
			r += "\t" + makeString(b[i].childIndex[1]);
		r += "\n";
	}
	return r;
}

string recordString(const vector<BranchPoint2D> &b){
	string r("");
	for(unsigned int i(0); i < b.size(); i++){
		r += makeString(i) + "\t" + makeString(b[i].position[0]) + "\t" + makeString(b[i].position[1]);
		for(unsigned int c(0); c < b[i].childIndex.size(); c++)
			r += "\t" + makeString(b[i].childIndex[c]);
		r += "\n";
	}
	return r;
}

vector<string> parseAtTabs(string ln){
	vector<string> parsedln;
	size_t pos(ln.find_first_of("\t"));
	while(pos != string::npos){
		parsedln.push_back(ln.substr(0, pos));
		ln = ln.substr(pos + 1);
		pos = ln.find_first_of("\t");
	}
	parsedln.push_back(ln);
	return parsedln;
}

// assumes bifurcations and 2D
vector<Bifurcation2D> readRecordedString(string s){
	stringstream strstr(s);
	string ln("");
	vector<Bifurcation2D> b;
	while(getline(strstr, ln)){
		vector<string> pln(parseAtTabs(ln)); // ignore pln[0]; it is the ordered index of b
		b.push_back(Bifurcation2D());
		b.back().position[0] = atof(pln[1].c_str());
		b.back().position[1] = atof(pln[2].c_str());
		for(unsigned int c(3); c < pln.size(); c++)
			b.back().childIndex.push_back(atoi(pln[c].c_str()));
	}
	return b;
}

// assumes branchpoints and 2D
vector<Bifurcation2D> readRecordedBranchPoint2DString(string s){
	stringstream strstr(s);
	string ln("");
	vector<Bifurcation2D> b;
	while(getline(strstr, ln)){
		vector<string> pln(parseAtTabs(ln)); // ignore pln[0]; it is the ordered index of b
		b.push_back(Bifurcation2D());
		b.back().position[0] = atof(pln[1].c_str());
		b.back().position[1] = atof(pln[2].c_str());
		for(unsigned int c(3); c < pln.size(); c++)
			b.back().childIndex.push_back(atoi(pln[c].c_str()));
	}
	return b;
}

double totalLength(vector<Bifurcation2D> &b){
	double len(0.0);
	for(unsigned int i(0); i < b.size(); i++){
		for(unsigned int c(0); c < b[i].childIndex.size(); c++)
			len += separation2D(b[i], b[b[i].childIndex[c]]);
	}
	return len;
}

double totalLength(vector<BranchPoint2D> &b){
	double len(0.0);
	for(unsigned int i(0); i < b.size(); i++){
		for(unsigned int c(0); c < b[i].childIndex.size(); c++)
			len += separation2D(b[i], b[b[i].childIndex[c]]);
	}
	return len;
}

double totalLength(vector<Bifurcation3D> &b){
	double len(0.0);
	for(unsigned int i(0); i < b.size(); i++){
		for(unsigned int c(0); c < b[i].childIndex.size(); c++)
			len += separation3D(b[i], b[b[i].childIndex[c]]);
	}
	return len;
}

double totalLength(vector<BranchPoint3D> &b){
	double len(0.0);
	for(unsigned int i(0); i < b.size(); i++){
		for(unsigned int c(0); c < b[i].childIndex.size(); c++)
			len += separation3D(b[i], b[b[i].childIndex[c]]);
	}
	return len;
}

double distanceFromHeart(vector<Bifurcation2D> &b, int index){
	double dist(0.0);
	while(index != 0){
		dist += separation2D(b[index], b[b[index].parentIndex]);
		index = b[index].parentIndex;
	}
	return dist;
}

double distanceFromHeart(vector<BranchPoint2D> &b, int index){
	double dist(0.0);
	while(index != 0){
		dist += separation2D(b[index], b[b[index].parentIndex]);
		index = b[index].parentIndex;
	}
	return dist;
}

double distanceFromHeart(vector<Bifurcation3D> &b, int index){
	double dist(0.0);
	while(index != 0){
		dist += separation3D(b[index], b[b[index].parentIndex]);
		index = b[index].parentIndex;
	}
	return dist;
}

double distanceFromHeart(vector<BranchPoint3D> &b, int index){
	double dist(0.0);
	while(index != 0){
		dist += separation3D(b[index], b[b[index].parentIndex]);
		index = b[index].parentIndex;
	}
	return dist;
}

// EXCLUDES heart (because distanceFromHeart is 0.0, and b.size() - 1 subtracts heart)
double averageDistanceFromHeart(vector<Bifurcation2D> &b){
	double distanceSum(0.0);
	for(unsigned int i(0); i < b.size(); i++)
		distanceSum += distanceFromHeart(b, i);
	if(b.size() > 1)
		return distanceSum/(b.size() - 1);
	return 0.0;
}

double maxDistanceFromHeart(vector<Bifurcation2D> &b){
	double maxDistance(-1.0);
	for(unsigned int i(0); i < b.size(); i++){
		double dist(distanceFromHeart(b, i));
		if(maxDistance < 0.0 || maxDistance < dist)
			maxDistance = dist;
	}
	return maxDistance;
}

// assumes heart is b[0];
bool connectedToHeart(vector<Bifurcation2D> &b, int index){
	while(index >= 0){
		if(index == 0)
			return true;
		index = b[index].parentIndex;
	}
	return false;
}

// assumes heart is b[0];
bool connectedToHeart(vector<BranchPoint2D> &b, int index){
	while(index >= 0){
		if(index == 0)
			return true;
		index = b[index].parentIndex;
	}
	return false;
}

// assumes heart is b[0];
bool connectedToHeart(vector<Bifurcation3D> &b, int index){
	while(index >= 0){
		if(index == 0)
			return true;
		index = b[index].parentIndex;
	}
	return false;
}

// assumes heart is b[0];
bool connectedToHeart(vector<BranchPoint3D> &b, int index){
	while(index >= 0){
		if(index == 0)
			return true;
		index = b[index].parentIndex;
	}
	return false;
}

// EXCLUDES heart (because distanceFromHeart is 0.0, and b.size() - 1 subtracts heart)
void maxAndAverageDistanceFromHeart(double &mDist, double &aDist, vector<Bifurcation2D> &b){
	mDist = -1.0;
	aDist = 0.0;
	int notConnected(0);
	for(unsigned int i(1); i < b.size(); i++){
		if(connectedToHeart(b, i)){
			double dist(distanceFromHeart(b, i));
			aDist += dist;
			if(mDist < 0.0 || mDist < dist)
				mDist = dist;
		}else
			notConnected++;
	}
	if(b.size() - notConnected > 1)
		aDist /= b.size() - notConnected - 1;
}

// EXCLUDES heart (because distanceFromHeart is 0.0, and b.size() - 1 subtracts heart)
void maxAndAverageDistanceFromHeart(double &mDist, double &aDist, vector<BranchPoint2D> &b){
	mDist = -1.0;
	aDist = 0.0;
	int notConnected(0);
	for(unsigned int i(1); i < b.size(); i++){
		if(connectedToHeart(b, i)){
			double dist(distanceFromHeart(b, i));
			aDist += dist;
			if(mDist < 0.0 || mDist < dist)
				mDist = dist;
		} else
			notConnected++;
	}
	if(b.size() - notConnected > 1)
		aDist /= b.size() - notConnected - 1;
}

// EXCLUDES heart (because distanceFromHeart is 0.0, and b.size() - 1 subtracts heart)
void maxAndAverageDistanceFromHeart(double &mDist, double &aDist, vector<Bifurcation3D> &b){
	mDist = -1.0;
	aDist = 0.0;
	int notConnected(0);
	for(unsigned int i(1); i < b.size(); i++){
		if(connectedToHeart(b, i)){
			double dist(distanceFromHeart(b, i));
			aDist += dist;
			if(mDist < 0.0 || mDist < dist)
				mDist = dist;
		}else
			notConnected++;
	}
	if(b.size() - notConnected > 1)
		aDist /= b.size() - notConnected - 1;
}

// EXCLUDES heart (because distanceFromHeart is 0.0, and b.size() - 1 subtracts heart)
void maxAndAverageDistanceFromHeart(double &mDist, double &aDist, vector<BranchPoint3D> &b){
	mDist = -1.0;
	aDist = 0.0;
	int notConnected(0);
	for(unsigned int i(1); i < b.size(); i++){
		if(connectedToHeart(b, i)){
			double dist(distanceFromHeart(b, i));
			aDist += dist;
			if(mDist < 0.0 || mDist < dist)
				mDist = dist;
		} else
			notConnected++;
	}
	if(b.size() - notConnected > 1)
		aDist /= b.size() - notConnected - 1;
}

double tripletMeasure(vector<Bifurcation2D> &b, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient){
	double totLen(-1.0), avePath(-1.0), maxPath(-1.0);
	if(totalLengthCoefficient != 0.0)
		totLen = totalLength(b);
	if(avePathLengthCoefficient != 0.0 || maxPathLengthCoefficient != 0.0)
		maxAndAverageDistanceFromHeart(maxPath, avePath, b);
	return totalLengthCoefficient*totLen + avePathLengthCoefficient*avePath + maxPathLengthCoefficient*maxPath;
}

double tripletMeasure(vector<Bifurcation3D> &b, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient){
	double totLen(-1.0), avePath(-1.0), maxPath(-1.0);
	if(totalLengthCoefficient != 0.0)
		totLen = totalLength(b);
	if(avePathLengthCoefficient != 0.0 || maxPathLengthCoefficient != 0.0)
		maxAndAverageDistanceFromHeart(maxPath, avePath, b);
	return totalLengthCoefficient*totLen + avePathLengthCoefficient*avePath + maxPathLengthCoefficient*maxPath;
}

double overrelaxationFactor(0.05);

void overrelax(vector<Bifurcation2D> &b, int i, double xPrev, double yPrev, double xAft, double yAft){
	//bool noOverrelax(false);
	//if(b[i].parentIndex > -1)
	//	noOverrelax = (xAft == b[b[i].parentIndex].position[0] && yAft == b[b[i].parentIndex].position[1]);
	//for(unsigned int c(0); c < b[i].childIndex.size(); c++)
	//	noOverrelax = noOverrelax || (xAft == b[b[i].childIndex[c]].position[0] && yAft == b[b[i].childIndex[c]].position[1]);
	//if(noOverrelax){
	//	b[i].position[0] = xAft;
	//	b[i].position[1] = yAft;
	//}else{
		b[i].position[0] = xAft + overrelaxationFactor*(xAft - xPrev);
		b[i].position[1] = yAft + overrelaxationFactor*(yAft - yPrev);
	//}
}

void overrelax(vector<Bifurcation3D> &b, int i, double xPrev, double yPrev, double zPrev, double xAft, double yAft, double zAft){
	b[i].position[0] = xAft + overrelaxationFactor*(xAft - xPrev);
	b[i].position[1] = yAft + overrelaxationFactor*(yAft - yPrev);
	b[i].position[2] = zAft + overrelaxationFactor*(zAft - zPrev);
}

void overrelax(vector<Bifurcation2D> &b, const vector<int> &gobs, double xPrev, double yPrev){
	if(gobs.size() > 0 && overrelaxationFactor != 0.0){
		double xAft(b[gobs[0]].position[0]), yAft(b[gobs[0]].position[1]);
		double xFinal(xAft + overrelaxationFactor*(xAft - xPrev));
		double yFinal(yAft + overrelaxationFactor*(yAft - yPrev));
		for(unsigned int g(0); g < gobs.size(); g++){
			b[gobs[g]].position[0] = xFinal;
			b[gobs[g]].position[1] = yFinal;
		}
	}
}

void overrelax(vector<Bifurcation3D> &b, const vector<int> &gobs, double xPrev, double yPrev, double zPrev){
	if(gobs.size() > 0 && overrelaxationFactor != 0.0){
		double xAft(b[gobs[0]].position[0]), yAft(b[gobs[0]].position[1]), zAft(b[gobs[0]].position[2]);
		double xFinal(xAft + overrelaxationFactor*(xAft - xPrev));
		double yFinal(yAft + overrelaxationFactor*(yAft - yPrev));
		double zFinal(zAft + overrelaxationFactor*(zAft - zPrev));
		for(unsigned int g(0); g < gobs.size(); g++){
			b[gobs[g]].position[0] = xFinal;
			b[gobs[g]].position[1] = yFinal;
			b[gobs[g]].position[2] = zFinal;
		}
	}
}

void settleBifurcations(vector<Bifurcation2D> &b, bool exclaim = false, double thresh = 1.0e-3, int maxIts = 10000){
	double change(thresh);
	double *prex = new double[b.size()];
	double *prey = new double[b.size()];
	for(int it(0); it < maxIts && change >= thresh; it++){
		//cout << "\n settleBifurcations(): b = " << makeStringTree(b);
		//if(exclaim)
		//	drawBifurcationTree("outputs/settleBifurcations_" + makeString(it) + ".png", b, 500, 500, 15);
		//change = 0.0;
		for(unsigned int i(0); i < b.size(); i++){
			prex[i] = b[i].position[0];
			prey[i] = b[i].position[1];
		}
		//for(unsigned int i(0); i < b.size(); i++){
		//	if(b[i].position[0] < 0.0 || b[i].position[0] > 4.0 || b[i].position[1] < 0.0 || b[i].position[1] > 4.0)
		//		cout << "\n START it " << it << " bi " << i << " at (" << b[i].position[0] << ", " << b[i].position[1] << ")" << endl;
		//}
		for(unsigned int i(0); i < b.size(); i++){
			if(b[i].parentIndex >= 0 && b[i].childIndex.size() == 2){
				if(exclaim){
					//cout << "\n settleBifurcations(): adjusting bifurcation " << i
					//	<< " with parent " << b[i].parentIndex
					//	<< " and children " << b[i].childIndex[0] << " and " << b[i].childIndex[1]
					//	<< endl;
					cout << "\nsettleBifurcations(): total length of b before is " << tripletMeasure(b, 1.0, 0.0, 0.0) << endl;
					drawTree("outputs/testSettling_" + makeString(it) + "_settle" + makeString(i) + "_before.png", consolidateDegenerates(b), 1000, 1000, 40);//drawBifurcationTree("outputs/testSettling_" + makeString(it) + "_settle" + makeString(i) + ".png", b, 500, 500, 15);
				}
				double x(0.0), y(0.0);
				double x1(b[i].position[0]), y1(b[i].position[1]);
				settleBifurcation(i, b, x, y, thresh);
				//if(x1 < 0.0 || x < 0.0 || y1 < 0.0 || y < 0.0 || x1 > 4.0 || x > 4.0 || y1 > 4.0 || y > 4.0){
				//	cout << "\n it " << it << " bi " << i << " F(" << globalCount << ") from (" << x1 << ", " << y1 << ")\t to (" << x << ", " << y << ")" << endl;
				//	
				//	drawBifurcationTree("outputs/F_" + makeString(globalCount) + ".png", b, 100, 100, 10);
				//	globalCount++;
				//}
				//change += separation2D(x, y, b[i].position[0], b[i].position[1]);
				//cout << "\n\tsubchange = " << change << " from (" << x << ", " << y << ") to (" << b[i].position[0] << ", " << b[i].position[1] << ")" << endl;
				overrelax(b, i, x1, y1, x, y);
				if(exclaim && separation2D(x, y, x1, y1) > 0.0){
					cout << "\nsettleBifurcations(): moved " << i << " by " << separation2D(x, y, x1, y1) << endl;
					cout << "\nsettleBifurcations(): total length of b is " << tripletMeasure(b, 1.0, 0.0, 0.0) << endl;
					drawTree("outputs/testSettling_" + makeString(it) + "_settle" + makeString(i) + "_consequently.png", consolidateDegenerates(b), 1000, 1000, 40);//drawBifurcationTree("outputs/testSettling_" + makeString(it) + "_settle" + makeString(i) + ".png", b, 500, 500, 15);
				}
				//b[i].position[0] = x;
				//b[i].position[1] = y;
			}
		}
		//for(unsigned int i(0); i < b.size(); i++){
		//	if(b[i].position[0] < 0.0 || b[i].position[0] > 4.0 || b[i].position[1] < 0.0 || b[i].position[1] > 4.0)
		//		cout << "\n MID it " << it << " bi " << i << " at (" << b[i].position[0] << ", " << b[i].position[1] << ")" << endl;
		//}
		//cout << "\nmidchange = " << change << endl;
		if(change < thresh || true){
			vector<vector<int> > gobs, boundaries;
			findDegenerateGobs(b, boundaries, gobs, thresh);
			for(unsigned int g(0); g < gobs.size(); g++){
				if(gobs[g].size() > 1){ // otherwise leave it to Fermat
					double x(b[gobs[g][0]].position[0]), y(b[gobs[g][0]].position[1]);
					vector<double> w(boundaries[g].size(), 1.0);
					TorresWeiszfeld(b, boundaries[g], w, gobs[g], thresh);
					overrelax(b, gobs[g], x, y);
					if(exclaim && separation2D(x, y, b[gobs[g][0]].position[0], b[gobs[g][0]].position[1]) > 0.0){
						cout << "\nsettleBifurcations(): moved gob of " << gobs[g].size() << " using " << gobs[g][0] << " by " << separation2D(x, y, b[gobs[g][0]].position[0], b[gobs[g][0]].position[1]) << endl;
						//drawBifurcationTree("outputs/testSettlingTW_" + makeString(it) + "_settle" + makeString(gobs[g][0]) + ".png", b, 500, 500, 15);
					}
					//double x2(b[gobs[g][0]].position[0]), y2(b[gobs[g][0]].position[1]);
					//if(x < 0.0 || x2 < 0.0 || y < 0.0 || y2 < 0.0 || x > 4.0 || x2 > 4.0 || y > 4.0 || y2 > 4.0)
					//	cout << "\n TW from (" << x << ", " << y << ")\t to (" << x2 << ", " << y2 << ")" << endl;
					//change += separation2D(x, y, b[gobs[g][0]].position[0], b[gobs[g][0]].position[1]);//gobs[g].size()*
				}
			}
		}
		//for(unsigned int i(0); i < b.size(); i++){
		//	if(b[i].position[0] < 0.0 || b[i].position[0] > 4.0 || b[i].position[1] < 0.0 || b[i].position[1] > 4.0)
		//		cout << "\n END it " << it << " bi " << i << " at (" << b[i].position[0] << ", " << b[i].position[1] << ")" << endl;
		//}
		//if(maxIts - it < 4){
		//	cout << "\nchange = " << change << endl;
		//	drawBifurcationTree("outputs/settleBifurcations_ExcessiveChange_" + makeString(globalCount + 0) + "_it" + makeString(it) + ".png", b, 1000, 1000, 50);
		//}
		if(exclaim)
			drawTree("outputs/settleBifurcations_" + makeString(it) + "_after.png", consolidateDegenerates(b), 1000, 1000, 40);//drawBifurcationTree("outputs/settleBifurcations_" + makeString(it) + "_after.png", b, 500, 500, 15);
		change = 0.0;
		for(unsigned int i(0); i < b.size(); i++)
			change += separation2D(prex[i], prey[i], b[i].position[0], b[i].position[1]);
		if(exclaim)
			cout << "  [sBit " << it << " c=" << change << "]  ";
	}

	if(change > thresh){
		globalCount++;
		cout << "\n settleBifurcations(): change = " << change << " >= thresh = " << thresh
			<< "\t globalCount = " << globalCount << "\n\t tree: " << makeStringTreeOneLine(b) << endl;
		//drawBifurcationTree("outputs/settleBifurcations_ExcessiveChange_" + makeString(globalCount) + ".png", b, 400, 300, 10);
	}
	delete[] prex;
	delete[] prey;
}

void settleBifurcations(vector<Bifurcation3D> &b, bool exclaim = false, double thresh = 1.0e-3, int maxIts = 10000){
	double change(thresh);
	double *prex = new double[b.size()];
	double *prey = new double[b.size()];
	double *prez = new double[b.size()];
	for(int it(0); it < maxIts && change >= thresh; it++){
		//cout << "\n settleBifurcations(): b = " << makeStringTree(b);
		//if(exclaim)
		//	drawBifurcationTree("outputs/settleBifurcations_" + makeString(it) + ".png", b, 500, 500, 15);
		//change = 0.0;
		for(unsigned int i(0); i < b.size(); i++){
			prex[i] = b[i].position[0];
			prey[i] = b[i].position[1];
			prez[i] = b[i].position[2];
		}
		//for(unsigned int i(0); i < b.size(); i++){
		//	if(b[i].position[0] < 0.0 || b[i].position[0] > 4.0 || b[i].position[1] < 0.0 || b[i].position[1] > 4.0)
		//		cout << "\n START it " << it << " bi " << i << " at (" << b[i].position[0] << ", " << b[i].position[1] << ")" << endl;
		//}
		for(unsigned int i(0); i < b.size(); i++){
			if(b[i].parentIndex >= 0 && b[i].childIndex.size() == 2){
				//if(exclaim){
					//cout << "\n settleBifurcations(): adjusting bifurcation " << i
					//	<< " with parent " << b[i].parentIndex
					//	<< " and children " << b[i].childIndex[0] << " and " << b[i].childIndex[1]
					//	<< endl;
					//drawBifurcationTree("outputs/testSettling_" + makeString(it) + "_settle" + makeString(i) + ".png", b, 500, 500, 15);
				//}
				double x(0.0), y(0.0), z(0.0);
				double x1(b[i].position[0]), y1(b[i].position[1]), z1(b[i].position[2]);
				settleBifurcation(i, b, x, y, z, thresh);
				if(exclaim && separation3D(x, y, z, x1, y1, z1) > 0.0){
					cout << "\nsettleBifurcations(): moved " << i << " by " << separation3D(x, y, z, x1, y1, z1) << endl;
					//drawBifurcationTree("outputs/testSettling_" + makeString(it) + "_settle" + makeString(i) + ".png", b, 500, 500, 15);
				}
				//if(x1 < 0.0 || x < 0.0 || y1 < 0.0 || y < 0.0 || x1 > 4.0 || x > 4.0 || y1 > 4.0 || y > 4.0){
				//	cout << "\n it " << it << " bi " << i << " F(" << globalCount << ") from (" << x1 << ", " << y1 << ")\t to (" << x << ", " << y << ")" << endl;
				//	
				//	drawBifurcationTree("outputs/F_" + makeString(globalCount) + ".png", b, 100, 100, 10);
				//	globalCount++;
				//}
				//change += separation2D(x, y, b[i].position[0], b[i].position[1]);
				//cout << "\n\tsubchange = " << change << " from (" << x << ", " << y << ") to (" << b[i].position[0] << ", " << b[i].position[1] << ")" << endl;
				overrelax(b, i, x1, y1, z1, x, y, z);
				//b[i].position[0] = x;
				//b[i].position[1] = y;
			}
		}
		//for(unsigned int i(0); i < b.size(); i++){
		//	if(b[i].position[0] < 0.0 || b[i].position[0] > 4.0 || b[i].position[1] < 0.0 || b[i].position[1] > 4.0)
		//		cout << "\n MID it " << it << " bi " << i << " at (" << b[i].position[0] << ", " << b[i].position[1] << ")" << endl;
		//}
		//cout << "\nmidchange = " << change << endl;
		if(change < thresh || true){
			vector<vector<int> > gobs, boundaries;
			findDegenerateGobs(b, boundaries, gobs, thresh);
			for(unsigned int g(0); g < gobs.size(); g++){
				if(gobs[g].size() > 1){ // otherwise leave it to Fermat
					double x(b[gobs[g][0]].position[0]), y(b[gobs[g][0]].position[1]), z(b[gobs[g][0]].position[2]);
					vector<double> w(boundaries[g].size(), 1.0);
					TorresWeiszfeld(b, boundaries[g], w, gobs[g], thresh);
					overrelax(b, gobs[g], x, y, z);
					if(exclaim && separation3D(x, y, z, b[gobs[g][0]].position[0], b[gobs[g][0]].position[1], b[gobs[g][0]].position[2]) > 0.0){
						cout << "\nsettleBifurcations(): moved gob of " << gobs[g].size() << " using " << gobs[g][0] << " by " << separation3D(x, y, z, b[gobs[g][0]].position[0], b[gobs[g][0]].position[1], b[gobs[g][0]].position[2]) << endl;
						//drawBifurcationTree("outputs/testSettlingTW_" + makeString(it) + "_settle" + makeString(gobs[g][0]) + ".png", b, 500, 500, 15);
					}
					//double x2(b[gobs[g][0]].position[0]), y2(b[gobs[g][0]].position[1]);
					//if(x < 0.0 || x2 < 0.0 || y < 0.0 || y2 < 0.0 || x > 4.0 || x2 > 4.0 || y > 4.0 || y2 > 4.0)
					//	cout << "\n TW from (" << x << ", " << y << ")\t to (" << x2 << ", " << y2 << ")" << endl;
					//change += separation2D(x, y, b[gobs[g][0]].position[0], b[gobs[g][0]].position[1]);//gobs[g].size()*
				}
			}
		}
		//for(unsigned int i(0); i < b.size(); i++){
		//	if(b[i].position[0] < 0.0 || b[i].position[0] > 4.0 || b[i].position[1] < 0.0 || b[i].position[1] > 4.0)
		//		cout << "\n END it " << it << " bi " << i << " at (" << b[i].position[0] << ", " << b[i].position[1] << ")" << endl;
		//}
		//if(maxIts - it < 4){
		//	cout << "\nchange = " << change << endl;
		//	drawBifurcationTree("outputs/settleBifurcations_ExcessiveChange_" + makeString(globalCount + 0) + "_it" + makeString(it) + ".png", b, 1000, 1000, 50);
		//}
		//if(exclaim)
		//	drawBifurcationTree("outputs/settleBifurcations_" + makeString(it) + "_after.png", b, 500, 500, 15);
		change = 0.0;
		for(unsigned int i(0); i < b.size(); i++)
			change += separation2D(prex[i], prey[i], b[i].position[0], b[i].position[1]);
		if(exclaim)
			cout << "  [sBit " << it << " c=" << change << "]  ";
	}

	if(change > thresh){
		globalCount++;
		cout << "\n settleBifurcations(): change = " << change << " >= thresh = " << thresh
			<< "\t globalCount = " << globalCount << "\n\t tree: " << makeStringTreeOneLine(b) << endl;
		//drawBifurcationTree("outputs/settleBifurcations_ExcessiveChange_" + makeString(globalCount) + ".png", b, 400, 300, 10);
	}
	delete[] prex;
	delete[] prey;
	delete[] prez;
}

//always modifies b
bool settleBifurcationsThresh(int checkEvery, double criticalMeasure, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient,
		vector<Bifurcation2D> &b, bool exclaim = false, double thresh = 1.0e-3, int maxIts = 10000){
	double change(thresh);
	double *prex = new double[b.size()];
	double *prey = new double[b.size()];
	double *preprex = new double[b.size()];
	double *preprey = new double[b.size()];// How much history will I need to end up keeping track of?
	for(unsigned int i(0); i < b.size(); i++){
		prex[i] = b[i].position[0];
		prey[i] = b[i].position[1];
	}
	for(int it(0); it < maxIts && change >= thresh; it++){
		for(unsigned int i(0); i < b.size(); i++){
			preprex[i] = prex[i];
			preprey[i] = prey[i];
			prex[i] = b[i].position[0];
			prey[i] = b[i].position[1];
		}
		for(unsigned int i(0); i < b.size(); i++){
			if(b[i].parentIndex >= 0 && b[i].childIndex.size() == 2){
				double x(0.0), y(0.0);
				settleBifurcation(i, b, x, y, thresh, false && exclaim);
				overrelax(b, i, prex[i], prey[i], x, y);
				if(it < 50 && exclaim && separation2D(x, y, prex[i], prey[i]) > 0.0){
					cout << "\nsettleBifurcationsThresh(): moved " << i << " by " << (1.0 + overrelaxationFactor)*separation2D(x, y, prex[i], prey[i])
						<< " from (" << prex[i] << ", " << prey[i] << ") to (" << b[i].position[0] << ", " << b[i].position[1] << ")" << endl;
					drawBifurcationTree("outputs/testSettling_" + makeString(it) + "_settle" + makeString(i) + ".png", b, 500, 500, 15);
				}
				
			}
		}
		if(change < thresh || true){
			vector<vector<int> > gobs, boundaries;
			findDegenerateGobs(b, boundaries, gobs, thresh);
			for(unsigned int g(0); g < gobs.size(); g++){
				if(gobs[g].size() > 1){ // otherwise leave it to Fermat
					//double x(b[gobs[g][0]].position[0]), y(b[gobs[g][0]].position[1]);
					vector<double> w(boundaries[g].size(), 1.0);
					TorresWeiszfeld(b, boundaries[g], w, gobs[g], thresh);
					overrelax(b, gobs[g], prex[gobs[g][0]], prey[gobs[g][0]]);
					if(it < 50 && exclaim && separation2D(prex[gobs[g][0]], prey[gobs[g][0]], b[gobs[g][0]].position[0], b[gobs[g][0]].position[1]) > 0.0){
						cout << "\nsettleBifurcationsThresh(): moved gob of " << gobs[g].size() << " consisting of " << makeVectorString(gobs[g]) << " using " << gobs[g][0] << " by " << (1.0 + overrelaxationFactor)*separation2D(prex[gobs[g][0]], prey[gobs[g][0]], b[gobs[g][0]].position[0], b[gobs[g][0]].position[1]) << endl;
						drawBifurcationTree("outputs/testSettlingTW_" + makeString(it) + "_settle" + makeString(gobs[g][0]) + ".png", b, 500, 500, 15);
					}
				}
			}
		}
		change = 0.0;
		for(unsigned int i(0); i < b.size(); i++){
				change += separation2D(prex[i], prey[i], b[i].position[0], b[i].position[1]);
				if(it == 49 && i == 18 && exclaim){
					cout << "\n " << i << "(" << prex[i] << ", " << prey[i] << ") to ("
						<< b[i].position[0] << ", " << b[i].position[1] << ")\n\t change now " << change << endl;
				}
		}
		double preChange(0.0);
		for(unsigned int i(0); i < b.size(); i++)
			preChange += separation2D(preprex[i], preprey[i], b[i].position[0], b[i].position[1]);
		if(change > preChange)
			change = preChange;
		if(exclaim)
			cout << "  [sBTit " << it << " c=" << change << "]  ";
		if(it > 3 && (it + 1)%checkEvery == 0){
			double meas(tripletMeasure(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
			if(change < 0.05*(meas - criticalMeasure)){
				delete[] prex;
				delete[] prey;
				delete[] preprex;
				delete[] preprey;
				return false;
			}
			checkEvery *= 2;
		}
	}
	

	if(change > thresh){
		globalCount++;
		cout << "\n settleBifurcationsThresh(): change = " << change << " >= thresh = " << thresh
			<< "\t globalCount = " << globalCount << "\n\t tree: " << makeStringTreeOneLine(b) << endl;
		drawBifurcationTree("outputs/settleBifurcations_ExcessiveChange_" + makeString(globalCount) + ".png", b, 400, 300, 10);
	}
	delete[] prex;
	delete[] prey;
	delete[] preprex;
	delete[] preprey;
	return true;
}

//always modifies b
bool settleBifurcationsThresh(int checkEvery, double criticalMeasure, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient,
		vector<Bifurcation3D> &b, bool exclaim = false, double thresh = 1.0e-3, int maxIts = 10000){
	double change(thresh);
	double *prex = new double[b.size()];
	double *prey = new double[b.size()];
	double *prez = new double[b.size()];
	double *preprex = new double[b.size()];
	double *preprey = new double[b.size()];// How much history will I need to end up keeping track of?
	double *preprez = new double[b.size()];
	for(unsigned int i(0); i < b.size(); i++){
		prex[i] = b[i].position[0];
		prey[i] = b[i].position[1];
		prez[i] = b[i].position[2];
	}
	for(int it(0); it < maxIts && change >= thresh; it++){
		for(unsigned int i(0); i < b.size(); i++){
			preprex[i] = prex[i];
			preprey[i] = prey[i];
			preprez[i] = prez[i];
			prex[i] = b[i].position[0];
			prey[i] = b[i].position[1];
			prez[i] = b[i].position[2];
		}
		for(unsigned int i(0); i < b.size(); i++){
			if(b[i].parentIndex >= 0 && b[i].childIndex.size() == 2){
				double x(0.0), y(0.0), z(0.0);
				settleBifurcation(i, b, x, y, z, thresh, false && exclaim);
				overrelax(b, i, prex[i], prey[i], prez[i], x, y, z);
				if(it < 50 && exclaim && separation3D(x, y, z, prex[i], prey[i], prez[i]) > 0.0){
					cout << "\nsettleBifurcationsThresh(): moved " << i << " by " << (1.0 + overrelaxationFactor)*separation3D(x, y, z, prex[i], prey[i], prez[i])
						<< " from (" << prex[i] << ", " << prey[i] << ", " << prez[i] << ") to (" << b[i].position[0] << ", " << b[i].position[1] << ", " << b[i].position[2] << ")" << endl;
					//drawBifurcationTree("outputs/testSettling_" + makeString(it) + "_settle" + makeString(i) + ".png", b, 500, 500, 15);
				}
				
			}
		}
		if(change < thresh || true){
			vector<vector<int> > gobs, boundaries;
			findDegenerateGobs(b, boundaries, gobs, thresh);
			for(unsigned int g(0); g < gobs.size(); g++){
				if(gobs[g].size() > 1){ // otherwise leave it to Fermat
					//double x(b[gobs[g][0]].position[0]), y(b[gobs[g][0]].position[1]);
					vector<double> w(boundaries[g].size(), 1.0);
					TorresWeiszfeld(b, boundaries[g], w, gobs[g], thresh);
					overrelax(b, gobs[g], prex[gobs[g][0]], prey[gobs[g][0]], prez[gobs[g][0]]);
					if(it < 50 && exclaim && separation3D(prex[gobs[g][0]], prey[gobs[g][0]], prez[gobs[g][0]], b[gobs[g][0]].position[0], b[gobs[g][0]].position[1], b[gobs[g][0]].position[2]) > 0.0){
						cout << "\nsettleBifurcationsThresh(): moved gob of " << gobs[g].size() << " consisting of " << makeVectorString(gobs[g]) << " using " << gobs[g][0] << " by " << (1.0 + overrelaxationFactor)*separation3D(prex[gobs[g][0]], prey[gobs[g][0]], prez[gobs[g][0]], b[gobs[g][0]].position[0], b[gobs[g][0]].position[1], b[gobs[g][0]].position[2]) << endl;
						//drawBifurcationTree("outputs/testSettlingTW_" + makeString(it) + "_settle" + makeString(gobs[g][0]) + ".png", b, 500, 500, 15);
					}
				}
			}
		}
		change = 0.0;
		for(unsigned int i(0); i < b.size(); i++){
				change += separation3D(prex[i], prey[i], prez[i], b[i].position[0], b[i].position[1], b[i].position[2]);
				if(it == 49 && i == 18 && exclaim){
					cout << "\n " << i << "(" << prex[i] << ", " << prey[i] << ", " << prez[i] << ") to ("
						<< b[i].position[0] << ", " << b[i].position[1] << ", " << b[i].position[2] << ")\n\t change now " << change << endl;
				}
		}
		double preChange(0.0);
		for(unsigned int i(0); i < b.size(); i++)
			preChange += separation3D(preprex[i], preprey[i], preprez[i], b[i].position[0], b[i].position[1], b[i].position[2]);
		if(change > preChange)
			change = preChange;
		if(exclaim)
			cout << "  [sBTit " << it << " c=" << change << "]  ";
		if(it > 3 && (it + 1)%checkEvery == 0){
			double meas(tripletMeasure(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
			if(change < 0.05*(meas - criticalMeasure)){
				delete[] prex;
				delete[] prey;
				delete[] prez;
				delete[] preprex;
				delete[] preprey;
				delete[] preprez;
				return false;
			}
			checkEvery *= 2;
		}
	}

	if(change > thresh){
		globalCount++;
		cout << "\n settleBifurcationsThresh(): change = " << change << " >= thresh = " << thresh
			<< "\t globalCount = " << globalCount << "\n\t tree: " << makeStringTreeOneLine(b) << endl;
		drawBifurcationTree("outputs/settleBifurcations_ExcessiveChange_" + makeString(globalCount) + ".png", b, 400, 300, 10);
	}
	delete[] prex;
	delete[] prey;
	delete[] preprex;
	delete[] preprey;
	return true;
}

void settleBifurcationsOLD(vector<Bifurcation2D> &b, bool exclaim = false, double thresh = 1.0e-6, int maxIts = 10000){
	double change(thresh);
	for(int it(0); it < maxIts && change >= thresh; it++){
		//cout << "\n settleBifurcations(): b = " << makeStringTree(b);
		if(exclaim)
			drawBifurcationTree("outputs/settleBifurcations_" + makeString(it) + ".png", b, 400, 300, 10);
		change = 0.0;
		for(unsigned int i(0); i < b.size(); i++){
			if(b[i].parentIndex >= 0 && b[i].childIndex.size() == 2){
				if(exclaim){
					cout << "\n settleBifurcations(): adjusting bifurcation " << i 
						<< " with parent " << b[i].parentIndex
						<< " and children " << b[i].childIndex[0] << " and " << b[i].childIndex[1]
						<< endl;
					drawBifurcationTree("outputs/testSettling_" + makeString(it) + "_settle" + makeString(i) + ".png", b, 400, 300, 10);
				}
				
				// separation(separation2D(b[i], b[b[i].parentIndex]) + 0.001);
				//double sumInverseSeparations(1.0/separation),
				//	x(b[b[i].parentIndex].position[0]/separation),
				//	y(b[b[i].parentIndex].position[1]/separation);
				//for(unsigned int c(0); c < b[i].childIndex.size(); c++){
				//	separation = separation2D(b[i], b[b[i].childIndex[c]]) + 0.001;
				//	sumInverseSeparations += 1.0/separation;
				//	x += b[b[i].childIndex[c]].position[0]/separation;
				//	y += b[b[i].childIndex[c]].position[1]/separation;
				//}
				//x /= sumInverseSeparations;
				//y /= sumInverseSeparations;

				double x(0.0), y(0.0);
				settleBifurcation(i, b, x, y, thresh);

				//cout << "\n settleBifurcations(): i = " << i << "; (x, y) = (" << x << ", " << y << ")" << endl;
				
				change += separation2D(x, y, b[i].position[0], b[i].position[1]);
				b[i].position[0] = x;
				b[i].position[1] = y;
			}
		}
		if(exclaim)
			drawBifurcationTree("outputs/settleBifurcations_" + makeString(it) + "_after.png", b, 400, 300, 10);
	}
	if(change > thresh){
		globalCount++;
		cout << "\n settleBifurcations(): change = " << change << " >= thresh = " << thresh
			<< "\t globalCount = " << globalCount << "\n\t tree: " << makeStringTreeOneLine(b) << endl;
		drawBifurcationTree("outputs/settleBifurcations_ExcessiveChange_" + makeString(globalCount) + ".png", b, 400, 300, 10);
	}
}

double settleBifurcationSubtree(int i, vector<Bifurcation2D> &b, bool exclaim = false, double thresh = 1.0e-6){
	double change(0.0);
	if(b[i].parentIndex >= 0 && b[i].childIndex.size() == 2){
		//double separation(separation2D(b[i], b[b[i].parentIndex]) + 0.001);
		//double sumInverseSeparations(1.0/separation),
		//	x(b[b[i].parentIndex].position[0]/separation),
		//	y(b[b[i].parentIndex].position[1]/separation);
		//for(unsigned int c(0); c < b[i].childIndex.size(); c++){
		//	separation = separation2D(b[i], b[b[i].childIndex[c]]) + 0.001;
		//	sumInverseSeparations += 1.0/separation;
		//	x += b[b[i].childIndex[c]].position[0]/separation;
		//	y += b[b[i].childIndex[c]].position[1]/separation;
		//}
		//x /= sumInverseSeparations;
		//y /= sumInverseSeparations;

		double x(0.0), y(0.0);
		settleBifurcation(i, b, x, y, thresh);

		change = separation2D(x, y, b[i].position[0], b[i].position[1]);
		b[i].position[0] = x;
		b[i].position[1] = y;
	}
	if(b[i].childIndex.size() > 0){
		for(unsigned int c(0); c < b[i].childIndex.size(); c++)
			change += settleBifurcationSubtree(b[i].childIndex[c], b, exclaim);
	}
	return change;
}

double settleBifurcationSubtree(int i, vector<Bifurcation3D> &b, bool exclaim = false, double thresh = 1.0e-6){
	double change(0.0);
	if(b[i].parentIndex >= 0 && b[i].childIndex.size() == 2){

		double x(0.0), y(0.0);
		settleBifurcation(i, b, x, y, thresh);

		change = separation2D(x, y, b[i].position[0], b[i].position[1]);
		b[i].position[0] = x;
		b[i].position[1] = y;
	}
	if(b[i].childIndex.size() > 0){
		for(unsigned int c(0); c < b[i].childIndex.size(); c++)
			change += settleBifurcationSubtree(b[i].childIndex[c], b, exclaim);
	}
	return change;
}

void settleBifurcationsFromSource(int sourceIndex, vector<Bifurcation2D> &b, bool exclaim = false, double thresh = 1.0e-6, int maxIts = 1000){
	double change(thresh);
	for(int it(0); it < maxIts && change >= thresh; it++){
		change = settleBifurcationSubtree(sourceIndex, b);
		if(exclaim)
			drawBifurcationTree("outputs/testSettling_" + makeString(sourceIndex) + "_" + makeString(it) + ".png", b, 400, 300, 10);
	}
	if(change > thresh)
		cout << "\n settleBifurcationsFromSource(): change = " << change << " >= thresh = " << thresh << endl;
}

void settleBifurcationsFromSource(int sourceIndex, vector<Bifurcation3D> &b, double thresh = 1.0e-6, int maxIts = 2000){
	double change(thresh);
	for(int it(0); it < maxIts && change >= thresh; it++){
		change = settleBifurcationSubtree(sourceIndex, b);
	}
	if(change > thresh)
		cout << "\n settleBifurcationsFromSource(): change = " << change << " >= thresh = " << thresh << endl;
}

void testSettle(){
	int width(80), height(60), border(5);
	vector<Bifurcation2D> b;
	b.push_back(Bifurcation2D());
	addChild(b, 0, 1.0, 1.0);
	addChild(b, 0, 1.0, -2.0);
	insertBifurcation(b, 0, 1, 2, 1.0, 2.0);
	//cout << "\n bifurcation tree: " << makeStringTree(b) << endl;
	//drawBifurcationTree("outputs/testSettleInitial.png", b, width, height, border);
	//settleBifurcations(b, false);
	//drawBifurcationTree("outputs/testSettleFinal.png", b, width, height, border);

	int count(0);
	cout << "\n\n count\t x1\t y1\t x2\t y2" << endl;
	for(double x1(-2); x1 < 3.0; x1 += 4.0){
		for(double y1(-1.0); y1 < 2.0; y1 += 2.0){
			for(double x2(-1.0); x2 < 2.0; x2 += 2.0){
				for(double y2(-2.0); y2 < 3.0; y2 += 4.0){
					b[1].position[0] = x1;
					b[1].position[1] = y1;
					b[2].position[0] = x2;
					b[2].position[1] = y2;
					cout << count << "\t" << x1 << "\t" << y1 << "\t" << x2 << "\t" << y2 << endl;
					settleBifurcations(b);
					drawBifurcationTree("outputs/testSettleFinal" + makeString(count) + ".png", b, width, height, border);
					count++;
				}
			}
		}
	}

	double thetaInc(0.1);
	int thetaCounter(0);
	cout << "\n label\t angle\t note heart = (" << b[0].position[0] << ", " << b[0].position[1] << ")" << endl;
	for(double theta(thetaInc); theta < acos(-1.0); theta += thetaInc){
		b[2].position[0] = b[1].position[0]*cos(theta) - b[1].position[1]*sin(theta);
		b[2].position[1] = b[1].position[0]*sin(theta) + b[1].position[1]*cos(theta);
		settleBifurcations(b);
		drawBifurcationTree("outputs/testSettleAngle" + makeString(thetaCounter) + ".png", b, width, height, border);
		cout << thetaCounter << "\t" << theta << endl;
		thetaCounter++;
	}
}

void segregateHalves(map<double, vector<int> > &m, vector<int> &a, vector<int> &b){
	unsigned int size(0);
	for(map<double, vector<int> >::iterator it(m.begin()); it != m.end(); it++)
		size += it->second.size();
	for(map<double, vector<int> >::iterator it(m.begin()); it != m.end(); it++){
		for(unsigned int j(0); j < it->second.size(); j++){
			if(a.size() < size/2)
				a.push_back(it->second[j]);
			else
				b.push_back(it->second[j]);
		}
	}
}

void symmetricallySegregateChildIndex(vector<Bifurcation2D> &b, double settleThresh, int i = 0){
	//return;
	if(b[i].childIndex.size() < 3) // nothing needs to be done; tips and penultimate bifurcations addressed elsewhere
		return;
	double tipCMx(0.0), tipCMy(0.0);
	for(unsigned int c(0); c < b[i].childIndex.size(); c++){
		tipCMx += b[b[i].childIndex[c]].position[0];
		tipCMy += b[b[i].childIndex[c]].position[1];
	}
	tipCMx /= b[i].childIndex.size();
	tipCMy /= b[i].childIndex.size();
	double slope(0.0), inter(0.0);
	vector<int> childIndex0, childIndex1; // segregated groups of child 0 and 1
	if(tipCMx == b[i].position[0]){ // [choose or is] vertical line: segregate by x value; in future, segregate to make less eccentric
		map<double, vector<int> > xPos;
		for(unsigned int c(0); c < b[i].childIndex.size(); c++)
			xPos[b[b[i].childIndex[c]].position[0]].push_back(b[i].childIndex[c]);
		segregateHalves(xPos, childIndex0, childIndex1);
	}else{ // sort by intercept
		double slope((tipCMy - b[i].position[1])/(tipCMx - b[i].position[0]));
		map<double, vector<int> > inter;
		for(unsigned int c(0); c < b[i].childIndex.size(); c++){
			double intercept(b[b[i].childIndex[c]].position[1] - (slope*b[b[i].childIndex[c]].position[0]));
			inter[intercept].push_back(b[i].childIndex[c]);
		}
		segregateHalves(inter, childIndex0, childIndex1);
	}
	b[i].childIndex.clear();

	//cout << "\n symmetricallySegregateChildIndex(): two groups for " << i << " are"
	//	<< "\n\t\t\t" << makeVectorString(childIndex0)
	//	<< "\n\t\t\t" << makeVectorString(childIndex1)
	//	<< endl;

	if(childIndex0.size() > 2){
		b.push_back(Bifurcation2D());
		b.back().parentIndex = i;
		b.back().childIndex = childIndex0;
		vector<double> w(childIndex0.size() + 1, 1.0);
		vector<int> movInd(1, b.size() - 1), constInd(childIndex0);
		constInd.push_back(i);
		b[i].childIndex.push_back(movInd[0]);
		TorresWeiszfeld(b, constInd, w, movInd, settleThresh);
		symmetricallySegregateChildIndex(b, settleThresh, b.size() - 1);
	}else if(childIndex0.size() == 2){
		b.push_back(Bifurcation2D());
		b.back().parentIndex = i;
		b[i].childIndex.push_back(b.size() - 1);
		b.back().childIndex = childIndex0;
		b[childIndex0[0]].parentIndex = b[childIndex0[1]].parentIndex = b.size() - 1;
		settleBifurcation(b.size() - 1, b, b.back().position[0], b.back().position[1], settleThresh);
	}else if(childIndex0.size() == 1){
		b[i].childIndex.push_back(childIndex0[0]);
		b[childIndex0[0]].parentIndex = i;
	}else{
		cout << "\n Error symmetricallySegregateChildIndex(): childIndex0 is empty." << endl;
		abort();
	}
	
	if(childIndex1.size() > 2){
		b.push_back(Bifurcation2D());
		b.back().parentIndex = i;
		b.back().childIndex = childIndex1;
		vector<double> w(childIndex1.size() + 1, 1.0);
		vector<int> movInd(1, b.size() - 1), constInd(childIndex1);
		constInd.push_back(i);
		b[i].childIndex.push_back(movInd[0]);
		TorresWeiszfeld(b, constInd, w, movInd, settleThresh);
		symmetricallySegregateChildIndex(b, settleThresh, b.size() - 1);
	}else if(childIndex1.size() == 2){
		b.push_back(Bifurcation2D());
		b.back().parentIndex = i;
		b[i].childIndex.push_back(b.size() - 1);
		b.back().childIndex = childIndex1;
		b[childIndex1[0]].parentIndex = b[childIndex1[1]].parentIndex = b.size() - 1;
		settleBifurcation(b.size() - 1, b, b.back().position[0], b.back().position[1], settleThresh);
	}else if(childIndex1.size() == 1){
		b[i].childIndex.push_back(childIndex1[0]);
		b[childIndex1[0]].parentIndex = i;
	}else{
		cout << "\n Error symmetricallySegregateChildIndex(): childIndex1 is empty." << endl;
		abort();
	}
}

void introduceBifurcationsSymmetricallyFromHeart(vector<Bifurcation2D> &b, double settleThresh){
	b.push_back(Bifurcation2D());
	b.back().parentIndex = 0;
	b[0].childIndex.push_back(b.size() - 1);
	for(unsigned int i(1); i < b.size() - 1; i++)
		b.back().childIndex.push_back(i);
	symmetricallySegregateChildIndex(b, settleThresh, b.size() - 1);
	settleBifurcations(b, false, settleThresh);
}

void introduceBifurcationsFromHeart(vector<Bifurcation2D> &b, int width = 300, int height = 300, int border = 12){
	unsigned int nchild(b[0].childIndex.size());
	cout << "\n\n nch\t totLength\t aveFromHeart\t maxFromheart" << endl;
	while(nchild > 1){
		/*
		double closestSeparation(-1.0);
		int closest1(-1), closest2(-1);
		int leastLevels(-1);
		for(int i(0); i < (int)nchild; i++){
			for(int j(i + 1); j < (int)nchild; j++){
				double separation(separation2D(b[b[0].childIndex[i]], b[b[0].childIndex[j]]));
				int levels(maxLevel(b, i) + maxLevel(b, j));
				if(closestSeparation < 0.0 || (closestSeparation > separation && leastLevels > levels)){
					closestSeparation = separation;
					closest1 = i;
					closest2 = j;
					leastLevels = levels;
				}
			}
		}
		insertBifurcation(b, 0, b[0].childIndex[closest1], b[0].childIndex[closest2]);
		settleBifurcations(b, false);
		drawBifurcationTree("outputs/testInsertBifurcations_" + makeString(nchild) + "heartChildren.png", b, width, height, border);
		nchild = b[0].childIndex.size();
		*/
		double *angles = new double[nchild];
		for(unsigned int c(0); c < nchild; c++)
			angles[c] = atan2(b[b[0].childIndex[c]].position[1] - b[0].position[1], b[b[0].childIndex[c]].position[0] - b[0].position[0]);
		int *orderedAngles = new int[nchild];
		for(unsigned int c(0); c < nchild; c++)
			orderedAngles[c] = c;
		for(unsigned int c(0); c < nchild; c++){
			for(unsigned int i(0); i < nchild - 1; i++){
				if(angles[orderedAngles[i]] > angles[orderedAngles[i + 1]]){
					int temp(orderedAngles[i]);
					orderedAngles[i] = orderedAngles[i + 1];
					orderedAngles[i + 1] = temp;
				}
			}
		}
		double closestAdjacency(-1.0), closestAngle(-1.0);
		int firstJoiner(-1);
		for(unsigned int i(0); i < nchild - 1; i++){
			double separation(separation2D(b[b[0].childIndex[orderedAngles[i]]], b[b[0].childIndex[orderedAngles[i + 1]]])), dtheta(angles[orderedAngles[i + 1]] - angles[orderedAngles[i]]);
			separation += 2.0e0*abs(1.0*maxLevel(b, b[0].childIndex[orderedAngles[i]]) - maxLevel(b, b[0].childIndex[orderedAngles[i + 1]]));
			//cout << "\n separation starting at " << i << " = " << separation << endl;
			if(firstJoiner < 0 || (closestAdjacency > separation) ){//&& closestAngle > dtheta)){
				firstJoiner = i;
				closestAdjacency = separation;
				closestAngle = dtheta;
			}
		}
		//cout << "\n CHOSE i = " << firstJoiner;
		//double separation(separation2D(b[b[0].childIndex[nchild - 1]], b[b[0].childIndex[0]])), dtheta(angles[orderedAngles[nchild - 1]] - angles[orderedAngles[0]]);
		//dtheta = 2*acos(-1.0) - dtheta;
		//if(firstJoiner < 0 || (closestAdjacency > separation && closestAngle > dtheta)){
		//	firstJoiner = nchild - 1;
		//	closestAdjacency = separation;
		//	closestAngle = dtheta;
		//}
		unsigned int ip0(firstJoiner), ip1(firstJoiner + 1);
		if(ip1 >= nchild)
			ip1 -= nchild;
		insertBifurcation(b, 0, b[0].childIndex[orderedAngles[ip0]], b[0].childIndex[orderedAngles[ip1]]);

		/*double *sumAngleSeparations = new double[nchild];
		double bestPairing(-1.0);
		for(unsigned int c(0); c < nchild; c++){
			sumAngleSeparations[c] = 0.0;
			for(unsigned int i(0); i < nchild - 1; i += 2){
				unsigned int ip0(c + i);
				if(ip0 >= nchild)
					ip0 -= nchild;
				unsigned int ip1(ip0 + 1);
				if(ip1 >= nchild)
					ip1 -= nchild;
				double dangle(angles[orderedAngles[ip1]] - angles[orderedAngles[ip0]]);
				if(dangle < 0.0)
					dangle = 2.0*acos(-1.0) - dangle;
				sumAngleSeparations[c] += dangle;
			}
			if(c == 0 || bestPairing > sumAngleSeparations[c])
				bestPairing = sumAngleSeparations[c];
		}
		int startPairing(-1);
		for(unsigned int c(0); c < nchild; c++){
			if(sumAngleSeparations[c] == bestPairing)
				startPairing = c;
		}
		cout << "\n orderedAngles:" << endl;
		for(unsigned int c(0); c < nchild; c++)
			cout << " child " << orderedAngles[c] << " at (" << b[b[0].childIndex[orderedAngles[c]]].position[0] << ", " << b[b[0].childIndex[orderedAngles[c]]].position[1] << ") has angle " << angles[orderedAngles[c]] << endl;
		cout << "\n sumAngleSeparations:" << endl;
		for(unsigned int c(0); c < nchild; c++)
			cout << c << "\t" << sumAngleSeparations[c] << endl;
		cout << "\n bestPairing = " << bestPairing << " starting with child " << startPairing << endl;
		cout << "\n total path length before iteration = " << totalLength(b) << endl;
		int *childIndex = new int[nchild];
		for(unsigned int c(0); c < nchild; c++)
			childIndex[c] = b[0].childIndex[c];
		for(unsigned int i(0); i < nchild - 1; i += 2){
			unsigned int ip0(startPairing + i);
			if(ip0 >= nchild)
				ip0 -= nchild;
			unsigned int ip1(ip0 + 1);
			if(ip1 >= nchild)
				ip1 -= nchild;
			insertBifurcation(b, 0, childIndex[orderedAngles[ip0]], childIndex[orderedAngles[ip1]]);
		}
		*/
		settleBifurcations(b, false);
		//cout << "\n total path length after iteration = " << totalLength(b) << endl;
		double maxDist(0.0), aveDist(0.0);
		for(unsigned int i(0); i < b.size(); i++){
			double dist(distanceFromHeart(b, i));
			if(maxDist < dist)
				maxDist = dist;
			aveDist += dist/b.size();
		}
		//cout << "\n average distance from heart = " << aveDist
		//	<< "\n maximum distance from heart = " << maxDist << endl;
		cout << nchild << "\t" << totalLength(b) << "\t" << aveDist << "\t" << maxDist << endl;
		drawBifurcationTree("outputs/testInsertBifurcations_" + makeString(nchild) + "heartChildren.png", b, width, height, border);
		nchild = b[0].childIndex.size();
		delete[] angles;
		delete[] orderedAngles;
		//delete[] sumAngleSeparations;
		//delete[] childIndex;
		
	}
}

int sourceIndex(vector<Bifurcation2D> &b, int index){
	while(b[index].parentIndex >= 0)
		index = b[index].parentIndex;
	return index;
}

//assumes b[0] is the heart
//COULD BE FASTER by only updating all nodes with the subset of nodes that are affected by the previous addition!!!!
void introduceBifurcationsFromSourcesSLOW(vector<Bifurcation2D> &b, bool print = false, int width = 300, int height = 300, int border = 12){
	bool disjoint(true);
	vector<int> sources;
	for(unsigned int i(1); i < b.size(); i++){
		if(b[i].parentIndex < 0)
			sources.push_back(i);
	}
	while(sources.size() > 0){
		//vector<int> sources;
		//for(unsigned int i(1); i < b.size(); i++){
		//	int source(sourceIndex(b, i));
		//	if(!containedIn(source, sources))
		//		sources.push_back(source);
		//}
		//cout << "\n sources: ";
		//for(unsigned int i(0); i < sources.size(); i++)
		//	cout << "\t" << sources[i];
		//cout << endl;

		double closestSeparation(-1.0);
		int join1(-1), join2(-1);
		for(unsigned int i(0); i < sources.size(); i++){
			for(unsigned int j(i + 1); j < sources.size(); j++){
				double sep(separation2D(b[sources[i]], b[sources[j]]));
				//sep += 10.0e0*abs(1.0*maxLevel(b, sources[i]) - maxLevel(b, sources[j]));
				if(closestSeparation < 0.0 ||
						(closestSeparation > sep) ){
					closestSeparation = sep;
					join1 = sources[i];
					join2 = sources[j];
				}
			}
		}
		//check if -1.0; if so, connect to heart, otherwise connect closest
		if(closestSeparation < 0.0){
			b[0].childIndex.push_back(sources[0]);
			b[sources[0]].parentIndex = 0;
			sources.clear();
		}else{
			addParent(b, join1, join2, 0.5*(b[join1].position[0] + b[join2].position[0]), 0.5*(b[join1].position[1] + b[join2].position[1]));
			//add the new bifurcation to sources; remove join1 and join2 from sources
			sources.push_back(b.size() - 1);
			for(unsigned int i(0); i < sources.size(); i++){
				if(sources[i] == join1 || sources[i] == join2){
					sources.erase(sources.begin() + i);
					i--;
				}
			}
		}
		settleBifurcations(b, false);
		if(print){
			cout << "\n sources.size() = " << sources.size()
				<< "\n\t join1 = " << join1 << "\t join2 = " << join2
				<< "\n closestSeparation = " << closestSeparation
				<< endl;
			drawBifurcationTree("outputs/introBiFromSources_" + makeString(sources.size()) + "heartChildren.png", b, width, height, border);
		}

		//update disjoint
		//disjoint = false;
		//for(unsigned int i(0); i < b.size(); i++){
		//	if(sourceIndex(b, i) != 0){
		//		disjoint = true;
		//		//cout << "\n bifurcation " << i << " is not connected to the heart; it terminates at " << sourceIndex(b, i) << endl;
		//	}
		//}
	}
}

// assumes BIfurcations
bool equivalentTrees(const vector<Bifurcation2D> &b1, const vector<Bifurcation2D> &b2, unsigned int index1 = 0, unsigned int index2 = 0){
	//cout << "\nequivalentTrees():\n\t\t\t" << makeStringTreeOneLine(b1) << "\n\t\t\t" << makeStringTreeOneLine(b2) << endl;
	if(b1.size() != b2.size())
		return false;
	if(index1 >= b1.size() && index2 >= b2.size())
		return false;
	if(index1 >= b1.size() || index2 >= b2.size())
		return false;
	if(b1[index1].childIndex.size() != b2[index2].childIndex.size())
		return false;
	if(b1[index1].childIndex.size() == 0)
		return index1 == index2;
	if(b1[index1].childIndex.size() == 1)
		return equivalentTrees(b1, b2, b1[index1].childIndex[0], b2[index2].childIndex[0]);
	for(unsigned int c(0); c < b1[index1].childIndex.size(); c++){
		unsigned int otherC((c + 1)%2);
		for(unsigned int cc(0); cc < b2[index2].childIndex.size(); cc++){
			unsigned int otherCC((cc + 1)%2);
			if(equivalentTrees(b1, b2, b1[index1].childIndex[c], b2[index2].childIndex[cc])
				&& equivalentTrees(b1, b2, b1[index1].childIndex[otherC], b2[index2].childIndex[otherCC]))
				return true;
		}
	}
	return false;
}

bool equivalentTrees(const vector<vector<Bifurcation2D> > &v, const vector<Bifurcation2D> &b){
	for(unsigned int i(0); i < v.size(); i++){
		if(equivalentTrees(v[i], b))
			return true;
	}
	return false;
}

void addTree(vector<vector<Bifurcation2D> > &v, vector<double> &ms, vector<Bifurcation2D> &b, double m){
	for(unsigned int i(0); i < ms.size(); i++){
		if(ms[i] > m){
			ms.insert(ms.begin() + i, m);
			v.insert(v.begin() + i, b);
			return;
		}
	}
	ms.push_back(m);
	v.push_back(b);
}

// assumes b1 and b2 have been properly consolidated
bool equivalentTrees(const vector<BranchPoint2D> &b1, const vector<BranchPoint2D> &b2, double thresh = 1.0e-6, unsigned int index1 = 0, unsigned int index2 = 0){
	//cout << "\nequivalentTrees():\n\t\t\t" << makeStringTreeOneLine(b1) << "\n\t\t\t" << makeStringTreeOneLine(b2) << endl;
	if(b1.size() != b2.size())
		return false;
	if(index1 >= b1.size() || index2 >= b2.size())
		return false;
	//cout << "\n testing childIndex sizes . . ." << endl;
	if(b1[index1].childIndex.size() != b2[index2].childIndex.size())
		return false;
	//cout << "\t same childIndex sizes for index1 = " << index1 << " and index2 = " << index2 << " of " << b1[index1].childIndex.size() << endl;
	if(b1[index1].childIndex.size() == 0)
		return separation2D(b1[index1], b2[index2]) < thresh;
	if(b1[index1].childIndex.size() == 1)
		return equivalentTrees(b1, b2, thresh, b1[index1].childIndex[0], b2[index2].childIndex[0]);
	bool *foundMatchingBranch1 = new bool[b1[index1].childIndex.size()];
	bool *foundMatchingBranch2 = new bool[b1[index1].childIndex.size()];
	for(unsigned int c(0); c < b1[index1].childIndex.size(); c++)
		foundMatchingBranch1[c] = foundMatchingBranch2[c]  = false;
	for(unsigned int c(0); c < b1[index1].childIndex.size(); c++){
		for(unsigned int cc(0); cc < b2[index2].childIndex.size() && !foundMatchingBranch1[c]; cc++){
			if(foundMatchingBranch2[cc])
				continue;
			if(equivalentTrees(b1, b2, thresh, b1[index1].childIndex[c], b2[index2].childIndex[cc]))
				foundMatchingBranch1[c] = foundMatchingBranch2[cc] = true;
		}
	}
	bool allMatched(true);
	for(unsigned int c(0); c < b1[index1].childIndex.size(); c++){
		//cout << "\n index1 " << index1 << " child " << c << " is " << notIndeed(foundMatchingBranch1[c]) << " matched." << endl;
		allMatched = allMatched && foundMatchingBranch1[c];
	}
	delete[] foundMatchingBranch1;
	delete[] foundMatchingBranch2;
	return allMatched;
}

int orderIndexOfFollowingNode(int firstNode, vector<int> &indexOfFollowingNode, vector<double> &distToNearest){
	bool changed(true);
	while(changed){
		
		changed = false;
		int curNode(firstNode), previousNode(-1), followingNode(indexOfFollowingNode[curNode]);
		while(followingNode > -1){

			//cout << "\n\t\t firstNode = " << firstNode
			//	<< "\n\t\t previousNode = " << previousNode << ";\t curNode = " << curNode << ";\t followingNode = " << followingNode
			//	<< "\n\t\t\t\t dist{cur->" << indexOfNearest[curNode] << "} = " << distToNearest[curNode]
			//	<< ";\t dist{fol->" << indexOfNearest[followingNode] << "} = " << distToNearest[followingNode]
			//	<< endl;
			//	if(followingNode > -1)
			//		cout << "\t\t\t\t indexOfFollowingNode[fol=" << followingNode << "] = " << indexOfFollowingNode[followingNode] << endl;
			//check to make sure not out of order; correct if out of order
			if(distToNearest[curNode] > distToNearest[followingNode]){
				if(curNode == firstNode)
					firstNode = followingNode;
				if(previousNode > -1)
					indexOfFollowingNode[previousNode] = followingNode;
				//int temp = indexOfFollowingNode[curNode];
				indexOfFollowingNode[curNode] = indexOfFollowingNode[followingNode];
				indexOfFollowingNode[followingNode] = curNode;
				previousNode = followingNode;
				//do not change curNode
				changed = true;
			}else{
				previousNode = curNode;
				curNode = indexOfFollowingNode[curNode];
			}
			followingNode = indexOfFollowingNode[curNode];
		}
	}
	return firstNode;
}

// assumes heart is 0
void introduceBifurcationsFromSourcesRandomly(vector<Bifurcation2D> &b){
	vector<int> indicesToConnect;
	for(unsigned int i(1); i < b.size(); i++)
		indicesToConnect.push_back(i);
	while(indicesToConnect.size() > 1){
		int childIndexIndex(kiss()%indicesToConnect.size());
		int childIndex1(indicesToConnect[childIndexIndex]);
		indicesToConnect.erase(indicesToConnect.begin() + childIndexIndex);
		childIndexIndex =  kiss()%indicesToConnect.size();
		int childIndex2(indicesToConnect[childIndexIndex]);
		indicesToConnect.erase(indicesToConnect.begin() + childIndexIndex);
		addParent(b, childIndex1, childIndex2);
		indicesToConnect.push_back(b.size() - 1);
		//cout << "\nintroduceBifurcationsFromSourcesRandomly(): " << childIndex1 << " and " << childIndex2 << " are now children of " << b.size() - 1 << endl;
	}
	b[0].childIndex.push_back(indicesToConnect[0]);
	b[indicesToConnect[0]].parentIndex = 0;
	settleBifurcations(b);
}

// assumes heart is 0
void introduceBifurcationsFromSourcesRandomlyNoSettle(vector<Bifurcation2D> &b){
	vector<int> indicesToConnect;
	for(unsigned int i(1); i < b.size(); i++)
		indicesToConnect.push_back(i);
	while(indicesToConnect.size() > 1){
		int childIndexIndex(kiss()%indicesToConnect.size());
		int childIndex1(indicesToConnect[childIndexIndex]);
		indicesToConnect.erase(indicesToConnect.begin() + childIndexIndex);
		childIndexIndex =  kiss()%indicesToConnect.size();
		int childIndex2(indicesToConnect[childIndexIndex]);
		indicesToConnect.erase(indicesToConnect.begin() + childIndexIndex);
		addParent(b, childIndex1, childIndex2);
		indicesToConnect.push_back(b.size() - 1);
		//cout << "\nintroduceBifurcationsFromSourcesRandomly(): " << childIndex1 << " and " << childIndex2 << " are now children of " << b.size() - 1 << endl;
	}
	b[0].childIndex.push_back(indicesToConnect[0]);
	b[indicesToConnect[0]].parentIndex = 0;
	//settleBifurcations(b);
}

void introduceBifurcationsFromSources(vector<Bifurcation2D> &b, bool print = false, int width = 300, int height = 300, int border = 12){
	if(print)
		cout << "\n\t Creating list of nearest neighbors . . ." << endl;
	//create list of nearest neighbors
	vector<double> distToNearest(b.size(), -1.0);
	vector<int> indexOfNearest(b.size(), -1), indexOfFollowingNode;
	int firstNode(1);
	for(unsigned int i(0); i < b.size(); i++){
		indexOfFollowingNode.push_back(i + 1);
		if(i > 0){
			for(unsigned int j(i + 1); j < b.size(); j++){
				double dist(separation2D(b[i], b[j]));
				if(dist < distToNearest[i] || distToNearest[i] < 0){
					distToNearest[i] = dist;
					indexOfNearest[i] = j;
				}
				if(dist < distToNearest[j] || distToNearest[j] < 0){
					distToNearest[j] = dist;
					indexOfNearest[j] = i;
				}
			}
		}
	}
	indexOfFollowingNode.back() = -1;
	if(print){
		cout << "\n indexOfFollowingNode: ";
		for(int curNode(firstNode); curNode > -1; curNode = indexOfFollowingNode[curNode])
			cout << "\n\t" << curNode << "\t (" << indexOfNearest[curNode] << " dist=" << distToNearest[curNode] << ")";
		cout << endl;

		cout << "\n\t\t . . . List of nearest neighbors created." << endl;
		cout << "\n indexOfFollowingNode: ";
		for(int curNode(firstNode); curNode > -1; curNode = indexOfFollowingNode[curNode])
			cout << "\n\t" << curNode << "\t (" << indexOfNearest[curNode] << " dist=" << distToNearest[curNode] << ")";
		cout << endl;
		cout << "\n\t Ordering list of nearest neighbors . . ." << endl;
	}

	//order list of nearest neighbors

	firstNode = orderIndexOfFollowingNode(firstNode, indexOfFollowingNode, distToNearest);

	if(print){
		cout << "\n\t\t . . . List of nearest neighbors ordered." << endl;
		cout << "\n indexOfFollowingNode: ";
		for(int curNode(firstNode); curNode > -1; curNode = indexOfFollowingNode[curNode])
			cout << "\n\t" << curNode << "\t (" << indexOfNearest[curNode] << " dist=" << distToNearest[curNode] << ")";
		cout << endl;
		cout << "\n\t Introducing bifurcations . . ." << endl;
	}

	//introduce bifurcations until there is only one bifurcation in the linked list;
	//	well, until there are only two left, at which point the algorithm is trivial
	while(indexOfFollowingNode[firstNode] > -1){
		if(print){
			cout << "\n adding Bifurcations indexOfFollowingNode: ";
			for(int curNode(firstNode); curNode > -1; curNode = indexOfFollowingNode[curNode])
				cout << "\n\t" << curNode << "\t (" << indexOfNearest[curNode] << " dist=" << distToNearest[curNode] << ")";
			cout << endl;
		}
		int join1(firstNode), join2(indexOfNearest[firstNode]);
		addParent(b, join1, join2, 0.5*(b[join1].position[0] + b[join2].position[0]), 0.5*(b[join1].position[1] + b[join2].position[1]));
		settleBifurcationsFromSource(b.size() - 1, b, false);
		distToNearest.push_back(-1.0);
		indexOfNearest.push_back(-1);
		if(print){
			cout << "\n\t\t\t b.size() = " << b.size()
				<< ";\t join1 = " << join1
				<< ";\t join2 = " << join2
				<< endl;
		}
		int curNode(firstNode);
		//modify list of nodes (indexOfFollowingNode)
		while(curNode > -1){
			int followingNode(indexOfFollowingNode[curNode]);
			if(firstNode == join1 || firstNode == join2)
				firstNode = followingNode;
			else if(followingNode == join1 || followingNode == join2) // remove followingNode from indexOfFollowingNode list
				indexOfFollowingNode[curNode] = indexOfFollowingNode[followingNode];
			else{//check nearest distances to/from the new node
				double dist(separation2D(b[curNode], b.back()));
				if(dist < distToNearest.back() || distToNearest.back() < 0.0){
					distToNearest.back() = dist;
					indexOfNearest.back() = curNode;
				}
				if(dist < distToNearest[curNode]){
					if(print){
						cout << "\n new node is now nearest to " << curNode << " since "
							<< dist << " < " << distToNearest[curNode] << endl;
					}
					distToNearest[curNode] = dist;
					indexOfNearest[curNode] = indexOfNearest.size() - 1;
				}
			}
			curNode = indexOfFollowingNode[curNode];
		}
		
		//if(indexOfFollowingNode[firstNode] > -1){
			//prehaps unnecessary
			//double dist(separation2D(b[firstNode], b.back()));
			//if(dist < distToNearest.back() || distToNearest.back() < 0.0){
			//	distToNearest.back() = dist;
			//	indexOfNearest.back() = firstNode;
			//}
			//if(dist < distToNearest[firstNode]){
			//	distToNearest[firstNode] = dist;
			//	indexOfNearest[firstNode] = indexOfNearest.size() - 1;
			//}

		indexOfFollowingNode.push_back(firstNode);
		firstNode = indexOfFollowingNode.size() - 1;
		firstNode = orderIndexOfFollowingNode(firstNode, indexOfFollowingNode, distToNearest);
			//make sure the removed nodes were not the closest node to a remaining node
			curNode = firstNode;
			while(curNode > -1){
				if(indexOfNearest[curNode] == join1 || indexOfNearest[curNode] == join2){
					distToNearest[curNode] = -1.0;
					int cur(firstNode);
					while(cur > -1){
						if(cur == curNode){
							cur = indexOfFollowingNode[cur];
							continue;
						}
						double dist(separation2D(b[curNode], b[cur]));
						if(dist < distToNearest[curNode] || distToNearest[curNode] < 0.0){
							distToNearest[curNode] = dist;
							indexOfNearest[curNode] = cur;
						}
						cur = indexOfFollowingNode[cur];
					}
				}
				curNode = indexOfFollowingNode[curNode];
			}

			firstNode = orderIndexOfFollowingNode(firstNode, indexOfFollowingNode, distToNearest);
			//perhaps unnecessary
			//distToNearest[indexOfNearest.back()] = distToNearest.back();
			//indexOfNearest[indexOfNearest.back()] = indexOfNearest.size() - 1;
			//indexOfFollowingNode.push_back(-2);
			//if(distToNearest.back() <= distToNearest[firstNode]){
			//	indexOfFollowingNode.back() = firstNode;
			//	firstNode = indexOfFollowingNode.size() - 1;
			//}else for(int curNode(firstNode); indexOfFollowingNode.back() < -1; curNode = indexOfFollowingNode[curNode]){
			//	if(indexOfFollowingNode[curNode] < 0){
			//		indexOfFollowingNode.back() = -1;
			//		indexOfFollowingNode[curNode] = indexOfFollowingNode.size() - 1;
			//	}else if(distToNearest.back() < distToNearest[indexOfFollowingNode[curNode]]){
			//		indexOfFollowingNode.back() = indexOfFollowingNode[curNode];
			//		indexOfFollowingNode[curNode] = indexOfFollowingNode.size() - 1;
			//	}
			//}
		//}
		if(print){
			int count(1);
			for(int curNode(firstNode); indexOfFollowingNode[curNode] > -1; curNode = indexOfFollowingNode[curNode])
				count++;
			cout << "\n indexOfFollowingNode has " << count << " elements." << endl;
			drawBifurcationTree("outputs/introBiFromSourcesFAST_" + makeString(count) + "heartChildren.png", b, width, height, border);
		}
		firstNode = orderIndexOfFollowingNode(firstNode, indexOfFollowingNode, distToNearest);
	}
	//drawBifurcationTree("outputs/introduceBifurcationsFromSources1.png", b, width, height, border);
	//int join1(firstNode), join2(b.size() - 1);
	//cout << "\n introduceBifurcationsFromSources(): join1 = " << join1 << ", join2 = " << join2 << endl;
	//addParent(b, join1, join2, 0.5*(b[join1].position[0] + b[join2].position[0]), 0.5*(b[join1].position[1] + b[join2].position[1]));
	//drawBifurcationTree("outputs/introduceBifurcationsFromSources2.png", b, width, height, border);
	
	//connect remaining source to heart and resettle all bifurcations
	b[0].childIndex.push_back(b.size() - 1);
	b.back().parentIndex = 0;
	//drawBifurcationTree("outputs/introduceBifurcationsFromSources3.png", b, width, height, border);
	//cout << "\n introduceBifurcationsFromSources(): invoking settleBifurcations one last time . . ." << endl;
	settleBifurcations(b, false);
	

	if(print)
		cout << "\n\t\t . . . Bifurcations added." << endl;

	//OLD and SLOW below
	//bool disjoint(true);
	//while(disjoint){
	//	vector<int> sources;
	//	for(unsigned int i(1); i < b.size(); i++){
	//		int source(sourceIndex(b, i));
	//		if(!containedIn(source, sources))
	//			sources.push_back(source);
	//	}
	//	//cout << "\n sources: ";
	//	//for(unsigned int i(0); i < sources.size(); i++)
	//	//	cout << "\t" << sources[i];
	//	//cout << endl;
	//
	//	double closestSeparation(-1.0);
	//	int join1(-1), join2(-1);
	//	for(unsigned int i(0); i < sources.size(); i++){
	//		for(unsigned int j(i + 1); j < sources.size(); j++){
	//			double sep(separation2D(b[sources[i]], b[sources[j]]));
	//			//sep += 10.0e0*abs(1.0*maxLevel(b, sources[i]) - maxLevel(b, sources[j]));
	//			if(closestSeparation < 0.0 ||
	//					(closestSeparation > sep) ){
	//				closestSeparation = sep;
	//				join1 = sources[i];
	//				join2 = sources[j];
	//			}
	//		}
	//	}
	//	//check if -1; if so, connect to heart, otherwise connect closest
	//	if(closestSeparation < 0.0){
	//		b[0].childIndex.push_back(sources[0]);
	//		b[sources[0]].parentIndex = 0;
	//	}else
	//		addParent(b, join1, join2, 0.5*(b[join1].position[0] + b[join2].position[0]), 0.5*(b[join1].position[1] + b[join2].position[1]));
	//	settleBifurcations(b, false);
	//	if(print)
	//		drawBifurcationTree("outputs/introBiFromSources_" + makeString(sources.size()) + "heartChildren.png", b, width, height, border);
	//
	//	//update disjoint
	//	disjoint = false;
	//	for(unsigned int i(0); i < b.size(); i++){
	//		if(sourceIndex(b, i) != 0){
	//			disjoint = true;
	//			//cout << "\n bifurcation " << i << " is not connected to the heart; it terminates at " << sourceIndex(b, i) << endl;
	//		}
	//	}
	//}
}

void introduceBifurcationsFromSources(vector<Bifurcation3D> &b){
	//create list of nearest neighbors
	vector<double> distToNearest(b.size(), -1.0);
	vector<int> indexOfNearest(b.size(), -1), indexOfFollowingNode;
	int firstNode(1);
	for(unsigned int i(0); i < b.size(); i++){
		indexOfFollowingNode.push_back(i + 1);
		if(i > 0){
			for(unsigned int j(i + 1); j < b.size(); j++){
				double dist(separation3D(b[i], b[j]));
				if(dist < distToNearest[i] || distToNearest[i] < 0){
					distToNearest[i] = dist;
					indexOfNearest[i] = j;
				}
				if(dist < distToNearest[j] || distToNearest[j] < 0){
					distToNearest[j] = dist;
					indexOfNearest[j] = i;
				}
			}
		}
	}
	indexOfFollowingNode.back() = -1;

	//order list of nearest neighbors

	firstNode = orderIndexOfFollowingNode(firstNode, indexOfFollowingNode, distToNearest);

	//introduce bifurcations until there is only one bifurcation in the linked list;
	//	well, until there are only two left, at which point the algorithm is trivial
	while(indexOfFollowingNode[firstNode] > -1){
		int join1(firstNode), join2(indexOfNearest[firstNode]);
		addParent(b, join1, join2, 0.5*(b[join1].position[0] + b[join2].position[0]), 0.5*(b[join1].position[1] + b[join2].position[1]), 0.5*(b[join1].position[2] + b[join2].position[2]));
		settleBifurcationsFromSource(b.size() - 1, b);
		distToNearest.push_back(-1.0);
		indexOfNearest.push_back(-1);
		int curNode(firstNode);
		//modify list of nodes (indexOfFollowingNode)
		while(curNode > -1){
			int followingNode(indexOfFollowingNode[curNode]);
			if(firstNode == join1 || firstNode == join2)
				firstNode = followingNode;
			else if(followingNode == join1 || followingNode == join2) // remove followingNode from indexOfFollowingNode list
				indexOfFollowingNode[curNode] = indexOfFollowingNode[followingNode];
			else{//check nearest distances to/from the new node
				double dist(separation3D(b[curNode], b.back()));
				if(dist < distToNearest.back() || distToNearest.back() < 0.0){
					distToNearest.back() = dist;
					indexOfNearest.back() = curNode;
				}
				if(dist < distToNearest[curNode]){
					distToNearest[curNode] = dist;
					indexOfNearest[curNode] = indexOfNearest.size() - 1;
				}
			}
			curNode = indexOfFollowingNode[curNode];
		}


		indexOfFollowingNode.push_back(firstNode);
		firstNode = indexOfFollowingNode.size() - 1;
		firstNode = orderIndexOfFollowingNode(firstNode, indexOfFollowingNode, distToNearest);
			//make sure the removed nodes were not the closest node to a remaining node
			curNode = firstNode;
			while(curNode > -1){
				if(indexOfNearest[curNode] == join1 || indexOfNearest[curNode] == join2){
					distToNearest[curNode] = -1.0;
					int cur(firstNode);
					while(cur > -1){
						if(cur == curNode){
							cur = indexOfFollowingNode[cur];
							continue;
						}
						double dist(separation3D(b[curNode], b[cur]));
						if(dist < distToNearest[curNode] || distToNearest[curNode] < 0.0){
							distToNearest[curNode] = dist;
							indexOfNearest[curNode] = cur;
						}
						cur = indexOfFollowingNode[cur];
					}
				}
				curNode = indexOfFollowingNode[curNode];
			}

			firstNode = orderIndexOfFollowingNode(firstNode, indexOfFollowingNode, distToNearest);
		firstNode = orderIndexOfFollowingNode(firstNode, indexOfFollowingNode, distToNearest);
	}
	b[0].childIndex.push_back(b.size() - 1);
	b.back().parentIndex = 0;
	settleBifurcations(b, false);
}

void testInsertBifurcations(){
	//kisset(22989, 23205, 2126);
	int NX(30), NY(NX), heartX(1*NX/2), heartY(0*NY/2);
	double noiseFactor(0.3);
	int width(NX*100), height(NY*100), border(width/25);
	if(height < width)
		border = height/25;
	vector<Bifurcation2D> b;
	b.push_back(Bifurcation2D());
	b.back().position[0] = 1.0*heartX;
	b.back().position[1] = 1.0*heartY;
	for(int x(0); x < NX; x++){
		for(int y(0); y < NY; y++){
			if(x != heartX || y != heartY){
				double positionX(1.0*x + noiseFactor*(rkiss() - 0.5)), positionY(1.0*y + noiseFactor*(rkiss() - 0.5));
				b.push_back(Bifurcation2D());
				b.back().position[0] = positionX;
				b.back().position[1] = positionY;
				//addChild(b, 0, positionX, positionY);
			}
		}
	}
	//cout << "\n bifurcation tree: " << makeStringTree(b) << endl;
	drawBifurcationTree("outputs/testInsertBifurcationsInitial.png", b, width, height, border);
	vector<Bifurcation2D> B(b);
	introduceBifurcationsFromHeart(b);
	time_t startOld(time(NULL));
	introduceBifurcationsFromSourcesSLOW(b, false);
	time_t endOld(time(NULL));
	time_t startNew(time(NULL));
	introduceBifurcationsFromSources(B, false);//used to be B
	time_t endNew(time(NULL));
	cout << "\n\n Old algorithm to introduce bifurcations took " << niceTime(difftime(endOld, startOld))
		<< "\n New algorithm to introduce bifurcations took " << niceTime(difftime(endNew, startNew)) << endl;
	drawBifurcationTree("outputs/testInsertBifurcationsFinal.png", b, width, height, border);
	drawBifurcationTree("outputs/testInsertBifurcationsFinalFAST.png", B, width, height, border); // used to be B
}

vector<Bifurcation2D> initializeGridWithNoise(int NX, int NY, int heartX, int heartY, double noiseFactor = 0.0){
	vector<Bifurcation2D> b;
	b.push_back(Bifurcation2D());
	b.back().position[0] = 1.0*heartX;
	b.back().position[1] = 1.0*heartY;
	for(int x(0); x < NX; x++){
		for(int y(0); y < NY; y++){
			if(x != heartX || y != heartY){
				double positionX(1.0*x + noiseFactor*(rkiss() - 0.5)), positionY(1.0*y + noiseFactor*(rkiss() - 0.5));
				b.push_back(Bifurcation2D());
				b.back().position[0] = positionX;
				b.back().position[1] = positionY;
			}
		}
	}
	return b;
}

double* getBinBottoms(double minX, double maxX, int numBins){
	double xInc((maxX - minX)/(numBins - 1));
	double *binBottoms = new double[numBins];
	if(xInc <= 0.0)
		return binBottoms;
	for(int i(0); i < numBins; i++)
		binBottoms[i] = minX + xInc*i;
	return binBottoms;
}

//assumes bins are even so that xInc = binBottoms[1] - binBOtoms[0]
double* probabilityDistribution(vector<double> &x, double *binBottoms, int numBins){
	double *p = new double[numBins];
	if(x.size() < 1)
		return p;
	double xInc(binBottoms[1] - binBottoms[0]);
	if(xInc == 0.0)
		return p;
	int *binCounts = new int[numBins];
	for(int i(0); i < numBins; i++)
		binCounts[i] = 0;
	for(unsigned int i(0); i < x.size(); i++){
		int bin(0);
		while(bin < numBins - 1 && binBottoms[bin + 1] < x[i])
			bin++;
		binCounts[bin]++;
	}
	for(int bin(0); bin < numBins; bin++){
		p[bin] = double(binCounts[bin])/x.size()/xInc;
	}
	delete[] binCounts;
	return p;
}

map<double, double> probabilityDistribution(vector<double> &x, int numBins){
	if(numBins < 2)
		numBins = 2;
	map<double, double> m;
	if(x.size() < 1)
		return m;
	double minX(x[0]), maxX(x[0]);
	for(unsigned int i(1); i < x.size(); i++){
		if(minX > x[i])
			minX = x[i];
		if(maxX < x[i])
			maxX = x[i];
	}
	double xInc((maxX - minX)/(numBins - 1));
	if(xInc == 0.0)
		return m;
	double *binBottoms = new double[numBins];
	int *binCounts = new int[numBins];
	for(int i(0); i < numBins; i++){
		binBottoms[i] = minX + xInc*i;
		//cout << "\n binBottoms[" << i << "] = " << binBottoms[i] << endl;
		binCounts[i] = 0;
	}
	for(unsigned int i(0); i < x.size(); i++){
		int bin(0);
		while(bin < numBins - 1 && binBottoms[bin + 1] < x[i])
			bin++;
		//cout << "\n  adding " << x[i] << " to bin " << bin << endl;
		binCounts[bin]++;
	}
	for(int bin(0); bin < numBins; bin++){
		//cout << "\n bin = " << bin << "\t\t\tbinCount = " << binCounts[bin] << endl;
		m[minX + (bin + 0.5)*xInc] = double(binCounts[bin])/x.size()/xInc;
	}
	delete[] binBottoms;
	delete[] binCounts;
	return m;
}

void printMap(map<double, double> &m, string name = ""){
	if(name.compare("") == 0)
		name = "x";
	cout << "\n" << name << "\tP(" << name << ")" << endl;
	for(map<double, double>::iterator it(m.begin()); it != m.end(); it++)
		cout << it->first << "\t" << it->second << endl;
}

void printMap(ofstream &outFile, map<double, double> &m){
	for(map<double, double>::iterator it(m.begin()); it != m.end(); it++)
		outFile << it->first << "\t" << it->second << endl;
	outFile << endl;
}

void testLengthAnalysis(){
	kisset(786, 28325, 12198);
	int NX(36), NY(36), heartx(1*NX/2), hearty(0*NY/2);
	int numBins(60);
	double thresh(0.01);
	double noiseFactor(0.3);
	int width(NX*100), height(NY*100), border(width/25);
	if(height < width)
		border = height/25;
	vector<Bifurcation2D> b = initializeGridWithNoise(NX, NY, heartx, hearty, noiseFactor);
	drawBifurcationTree("outputs/testLengthAnalysisInitial.png", b, width, height, border);
	introduceBifurcationsFromSources(b, true);
	drawBifurcationTree("outputs/testLengthAnalysisFinal.png", b, width, height, border);
	vector<double> lambdaLs;
	int ignoreCount(0);
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].childIndex.size() == 2){
			double len1(separation2D(b[i], b[b[i].childIndex[0]])), len2(separation2D(b[i], b[b[i].childIndex[1]]));
			if(len1 < thresh || len2 < thresh){
				cout << "Ignoring children of b[" << i << "] with len1 = " << len1 << " and len2 = " << len2 << endl;
				ignoreCount++;
				continue;
			}
			if(len1 > len2)
				lambdaLs.push_back(len2/len1);
			else
				lambdaLs.push_back(len1/len2);
		}
	}
	cout << "\nIgnored " << 100.0*double(ignoreCount)/b.size() << "% of bifurcations." << endl;
	cout << "\nlambdaLs:" << endl;
	for(unsigned int i(0); i < lambdaLs.size(); i++)
		cout << lambdaLs[i] << endl;
	map<double, double> lambdaLDensity = probabilityDistribution(lambdaLs, numBins);
	map<double, double>::iterator it(lambdaLDensity.begin());
	it++;
	double checkInc(it->first - lambdaLDensity.begin()->first), checkSum(0.0);
	for(map<double, double>::iterator it(lambdaLDensity.begin()); it != lambdaLDensity.end(); it++)
		checkSum += it->second*checkInc;
	cout << "\n\n Note checkSum = " << checkSum << "  (with checkInc = " << checkInc << ")" << endl;
	printMap(lambdaLDensity, "len");
	lambdaLDensity = probabilityDistribution(lambdaLs, 3*numBins/2);
	cout << "\n\n with " << 3*numBins/2 << " bins:" << endl;
	printMap(lambdaLDensity, "len");
	lambdaLDensity = probabilityDistribution(lambdaLs, 4*numBins/2);
	cout << "\n\n with " << 4*numBins/2 << " bins:" << endl;
	printMap(lambdaLDensity, "len");
	lambdaLDensity = probabilityDistribution(lambdaLs, 6*numBins/2);
	cout << "\n\n with " << 6*numBins/2 << " bins:" << endl;
	printMap(lambdaLDensity, "len");
}

void multipleLengths(){
	int startNX(5), incNX(5), startNY(5), incNY(5), numN(10), numGrids(100), heartxFact(1), heartyFact(0);
	string uniquer("NX");
	int numBins(20);
	double thresh(0.1);
	double noiseFactor(0.3);
	cout << "\n startNX = " << startNX << "\n incNX = " << incNX
		<< "\n startNY = " << startNY << "\n incNY = " << incNY
		<< "\n numN = " << numN << "\n numGrids = " << numGrids
		<< "\n heartxFact = " << heartxFact << "\n heartyFact = " << heartyFact
		<< "\n uniquer = " << uniquer << "\n numBins = " << numBins
		<< "\n thresh = " << thresh << "\n noiseFactor = " << noiseFactor
		<< endl;
	string coarseFn("outputs/" + uniquer + "_Ncoarse.dat");
	ofstream coarse(coarseFn.c_str());
	string fineFn("outputs/" + uniquer + "_Nfine.dat");
	ofstream fine(fineFn.c_str());
	string ultrafineFn("outputs/" + uniquer + "_Nultrafine.dat");
	ofstream ultrafine(ultrafineFn.c_str());
	for(int N(0); N < numN; N++){
		int NX(startNX + N*incNX), NY(startNY + N*incNY);
		vector<double> lambdaLs;
		int ignoreCount(0);
		vector<Bifurcation2D> b;
		for(int g(0); g < numGrids; g++){
			b = initializeGridWithNoise(NX, NY, heartxFact*NX/2, heartyFact*NY/2, noiseFactor);
			introduceBifurcationsFromSources(b);
			for(unsigned int i(0); i < b.size(); i++){
				if(b[i].childIndex.size() == 2){
					double len1(separation2D(b[i], b[b[i].childIndex[0]])), len2(separation2D(b[i], b[b[i].childIndex[1]]));
					if(len1 < thresh || len2 < thresh){
						ignoreCount++;
						continue;
					}
					if(len1 > len2)
						lambdaLs.push_back(len2/len1);
					else
						lambdaLs.push_back(len1/len2);
				}
			}
			cout << "*";
		}
		cout << "\nIgnored " << 100.0*double(ignoreCount)/(b.size()*numGrids) << "% of bifurcations for grid " << NX << "x" << NY << "." << endl;
		map<double, double> lambdaLDensity = probabilityDistribution(lambdaLs, numBins);
		printMap(coarse, lambdaLDensity);
		lambdaLDensity = probabilityDistribution(lambdaLs, 3*numBins/2);
		printMap(fine, lambdaLDensity);
		lambdaLDensity = probabilityDistribution(lambdaLs, 4*numBins/2);
		printMap(ultrafine, lambdaLDensity);
	}
	coarse.close();
	fine.close();
	ultrafine.close();
}

vector<double> validBranchingLengthRatios(vector<Bifurcation2D> &b, int &ignoreCount, double thresh = 0.01){
	vector<double> lambdaLs;
	ignoreCount = 0;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].childIndex.size() == 2){
			double len1(separation2D(b[i], b[b[i].childIndex[0]])), len2(separation2D(b[i], b[b[i].childIndex[1]]));
			if(len1 < thresh || len2 < thresh){
				ignoreCount++;
				continue;
			}
			if(len1 > len2)
				lambdaLs.push_back(len2/len1);
			else
				lambdaLs.push_back(len1/len2);
		}
	}
	return lambdaLs;
}

vector<double> validBranchingLengthRatios(vector<BranchPoint2D> &b, int &ignoreCount, double thresh = 0.01){
	vector<double> lambdaLs;
	ignoreCount = 0;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].childIndex.size() > 0){
			for(unsigned int c(0); c < b[i].childIndex.size() - 1; c++){
				double len1(separation2D(b[i], b[b[i].childIndex[c]]));
				for(unsigned int cc(c + 1); cc < b[i].childIndex.size(); cc++){
					double len2(separation2D(b[i], b[b[i].childIndex[cc]]));
					if(len1 >= thresh && len2 >= thresh){
						if(len1 > len2)
							lambdaLs.push_back(len2/len1);
						else
							lambdaLs.push_back(len1/len2);
					}else
						ignoreCount++;
				}
			}
		}
	}
	return lambdaLs;
}

vector<double> validBranchingLengthRatiosWithConstant(vector<BranchPoint2D> &b, int &ignoreCount, double capConst, double thresh = 0.01){
	vector<double> lambdaLs;
	ignoreCount = 0;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].childIndex.size() > 0){
			for(unsigned int c(0); c < b[i].childIndex.size() - 1; c++){
				double len1(separation2D(b[i], b[b[i].childIndex[c]]));
				if(b[b[i].childIndex[c]].childIndex.size() < 1)
					len1 += capConst;
				for(unsigned int cc(c + 1); cc < b[i].childIndex.size(); cc++){
					double len2(separation2D(b[i], b[b[i].childIndex[cc]]));
					if(b[b[i].childIndex[cc]].childIndex.size() < 1)
						len2 += capConst;
					if(len1 >= thresh && len2 >= thresh){
						if(len1 > len2)
							lambdaLs.push_back(len2/len1);
						else
							lambdaLs.push_back(len1/len2);
					}else
						ignoreCount++;
				}
			}
		}
	}
	return lambdaLs;
}

vector<double> validBranchingLengthRatiosIgnoreTips(vector<BranchPoint2D> &b, int &ignoreCount, double thresh = 0.01){
	vector<double> lambdaLs;
	ignoreCount = 0;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].childIndex.size() > 0){
			for(unsigned int c(0); c < b[i].childIndex.size() - 1; c++){
				double len1(separation2D(b[i], b[b[i].childIndex[c]]));
				for(unsigned int cc(c + 1); cc < b[i].childIndex.size(); cc++){
					double len2(separation2D(b[i], b[b[i].childIndex[cc]]));
					if(b[b[i].childIndex[c]].childIndex.size() < 1 || b[b[i].childIndex[cc]].childIndex.size() < 1){
						ignoreCount++;
						continue;
					}
					if(len1 >= thresh && len2 >= thresh){
						if(len1 > len2)
							lambdaLs.push_back(len2/len1);
						else
							lambdaLs.push_back(len1/len2);
					}else
						ignoreCount++;
				}
			}
		}
	}
	return lambdaLs;
}

vector<double> validParentChildLengthRatios(vector<Bifurcation2D> &b, int &ignoreCount, double thresh = 0.01){
	vector<double> gammas;
	ignoreCount = 0;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].parentIndex >= 0 && b[i].childIndex.size() == 2){
			double lenp(separation2D(b[i], b[b[i].parentIndex])),
				len1(separation2D(b[i], b[b[i].childIndex[0]])),
				len2(separation2D(b[i], b[b[i].childIndex[1]]));
			if(lenp < thresh){
				ignoreCount += 2;
				continue;
			}
			if(len1 < thresh)
				ignoreCount++;
			else
				gammas.push_back(len1/lenp);
			if(len2 < thresh)
				ignoreCount++;
			else
				gammas.push_back(len2/lenp);
		}
	}
	return gammas;
}

vector<double> validParentChildLengthRatios(vector<BranchPoint2D> &b, int &ignoreCount, double thresh = 0.01){
	vector<double> gammas;
	ignoreCount = 0;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].parentIndex > -1){
			for(unsigned int c(0); c < b[i].childIndex.size(); c++){
				double lenp(separation2D(b[i], b[b[i].parentIndex])),
					lenc(separation2D(b[i], b[b[i].childIndex[c]]));
				if(lenp >= thresh && lenc >= thresh)
					gammas.push_back(lenc/lenp);
				else
					ignoreCount++;
			}
		}
	}
	return gammas;
}

vector<double> allParentChildLengthRatios(vector<BranchPoint2D> &b){
	vector<double> gammas;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].parentIndex > -1){
			for(unsigned int c(0); c < b[i].childIndex.size(); c++){
				double lenp(separation2D(b[i], b[b[i].parentIndex])),
					lenc(separation2D(b[i], b[b[i].childIndex[c]]));
				if(lenp > 0.0)
					gammas.push_back(lenc/lenp);
			}
		}
	}
	return gammas;
}

int segsToHeart(const vector<BranchPoint2D> &b, int i){
	int segCount(0);
	while(b[i].parentIndex > -1){
		segCount++;
		i = b[i].parentIndex;
	}
	return segCount;
}

void tipSetAdd(const vector<BranchPoint2D> &b, int i, vector<int> &tSet){
	if(b[i].childIndex.size() < 1)
		tSet.push_back(i);
	else for(unsigned int c(0); c < b[i].childIndex.size(); c++)
		tipSetAdd(b, b[i].childIndex[c], tSet);
}

vector<int> tipSet(const vector<BranchPoint2D> &b, int i){
	vector<int> tSet;
	tipSetAdd(b, i, tSet);
	return tSet;
}

void validParentChildLengthRatiosTipsLevelsToHeartLengths(const vector<BranchPoint2D> &b, int &ignoreCount, vector<double> &gammas, vector<int> &tipCounts, vector<int> &levelsToHeart, vector<double> &lengths, double thresh = 0.01){
	gammas.clear();
	tipCounts.clear();
	levelsToHeart.clear();
	ignoreCount = 0;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].parentIndex > -1){
			for(unsigned int c(0); c < b[i].childIndex.size(); c++){
				double lenp(separation2D(b[i], b[b[i].parentIndex])),
					lenc(separation2D(b[i], b[b[i].childIndex[c]]));
				if(lenp >= thresh && lenc >= thresh){
					gammas.push_back(lenc/lenp);
					vector<int> tSet(tipSet(b, i));
					tipCounts.push_back(tSet.size());
					levelsToHeart.push_back(segsToHeart(b, i));
					lengths.push_back(lenp);
				}else
					ignoreCount++;
			}
		}
	}
}

int allNumTips(const vector<Bifurcation2D> &b, int startIndex, int *nTips){
	if(b[startIndex].childIndex.size() < 1)
		nTips[startIndex] = 1;
	else{
		nTips[startIndex] = 0;
		for(unsigned int c(0); c < b[startIndex].childIndex.size(); c++)
			nTips[startIndex] += allNumTips(b, b[startIndex].childIndex[c], nTips);
	}
	return nTips[startIndex];
}

int* numTips(const vector<Bifurcation2D> &b){
	int *nTips = new int[b.size()];
	for(unsigned int i(0); i < b.size(); i++)
		nTips[i] = 0;
	allNumTips(b, 0, nTips);
	return nTips;
}

int allNumTips(const vector<BranchPoint2D> &b, int startIndex, int *nTips){
	if(b[startIndex].childIndex.size() < 1)
		nTips[startIndex] = 1;
	else{
		nTips[startIndex] = 0;
		for(unsigned int c(0); c < b[startIndex].childIndex.size(); c++)
			nTips[startIndex] += allNumTips(b, b[startIndex].childIndex[c], nTips);
	}
	return nTips[startIndex];
}

int* numTips(const vector<BranchPoint2D> &b){
	int *nTips = new int[b.size()];
	for(unsigned int i(0); i < b.size(); i++)
		nTips[i] = 0;
	allNumTips(b, 0, nTips);
	return nTips;
}

int numTips(const vector<BranchPoint2D> &b, int i){
	if(b[i].childIndex.size() < 1)
		return 1;
	int count(0);
	for(unsigned int c(0); c < b[i].childIndex.size(); c++)
		count += numTips(b, b[i].childIndex[c]);
	return count;
}

vector<double> validTipCountRatios(vector<BranchPoint2D> &b, int &ignoreCount, double thresh = 0.01){
	vector<double> tipRatios;
	ignoreCount = 0;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].childIndex.size() > 0){
			for(unsigned int c(0); c < b[i].childIndex.size() - 1; c++){
				int tips1(numTips(b, b[i].childIndex[c]));
				double len1(separation2D(b[i], b[b[i].childIndex[c]]));
				for(unsigned int cc(c + 1); cc < b[i].childIndex.size(); cc++){
					int tips2(numTips(b, b[i].childIndex[cc]));
					double len2(separation2D(b[i], b[b[i].childIndex[cc]]));
					if(len1 >= thresh || len2 >= thresh){ // to be consistent with lambdaLs
						if(tips1 < tips2)
							tipRatios.push_back(double(tips1)/double(tips2));
						else
							tipRatios.push_back(double(tips2)/double(tips1));
					}else
						ignoreCount++;
				}
			}
		}
	}
	return tipRatios;
}

vector<double> validBranchingLengthRatios(vector<Bifurcation3D> &b, int &ignoreCount, double thresh = 0.01){
	vector<double> lambdaLs;
	ignoreCount = 0;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].childIndex.size() == 2){
			double len1(separation3D(b[i], b[b[i].childIndex[0]])), len2(separation3D(b[i], b[b[i].childIndex[1]]));
			if(len1 < thresh || len2 < thresh){
				ignoreCount++;
				continue;
			}
			if(len1 > len2)
				lambdaLs.push_back(len2/len1);
			else
				lambdaLs.push_back(len1/len2);
		}
	}
	return lambdaLs;
}

vector<double> validBranchingLengthRatios(vector<BranchPoint3D> &b, int &ignoreCount, double thresh = 0.01){
	vector<double> lambdaLs;
	ignoreCount = 0;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].childIndex.size() > 0){
			for(unsigned int c(0); c < b[i].childIndex.size() - 1; c++){
				double len1(separation3D(b[i], b[b[i].childIndex[c]]));
				for(unsigned int cc(c + 1); cc < b[i].childIndex.size(); cc++){
					double len2(separation3D(b[i], b[b[i].childIndex[cc]]));
					if(len1 >= thresh || len2 >= thresh){
						if(len1 > len2)
							lambdaLs.push_back(len2/len1);
						else
							lambdaLs.push_back(len1/len2);
					}else
						ignoreCount++;
				}
			}
		}
	}
	return lambdaLs;
}

vector<double> validParentChildLengthRatios(vector<Bifurcation3D> &b, int &ignoreCount, double thresh = 0.01){
	vector<double> gammas;
	ignoreCount = 0;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].parentIndex >= 0 && b[i].childIndex.size() == 2){
			double lenp(separation3D(b[i], b[b[i].parentIndex])),
				len1(separation3D(b[i], b[b[i].childIndex[0]])),
				len2(separation3D(b[i], b[b[i].childIndex[1]]));
			if(len1 < thresh)
				ignoreCount++;
			else
				gammas.push_back(len1/lenp);
			if(len2 < thresh)
				ignoreCount++;
			else
				gammas.push_back(len2/lenp);
		}
	}
	return gammas;
}

vector<double> validParentChildLengthRatios(vector<BranchPoint3D> &b, int &ignoreCount, double thresh = 0.01){
	vector<double> gammas;
	ignoreCount = 0;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].parentIndex > -1){
			for(unsigned int c(0); c < b[i].childIndex.size(); c++){
				double lenp(separation3D(b[i], b[b[i].parentIndex])),
					lenc(separation3D(b[i], b[b[i].childIndex[c]]));
				if(lenp >= thresh && lenc >= thresh)
					gammas.push_back(lenc/lenp);
				else
					ignoreCount++;
			}
		}
	}
	return gammas;
}

int numTips(const vector<BranchPoint3D> &b, int i){
	if(b[i].childIndex.size() < 1)
		return 1;
	int count(0);
	for(unsigned int c(0); c < b[i].childIndex.size(); c++)
		count += numTips(b, b[i].childIndex[c]);
	return count;
}

vector<double> validTipCountRatios(vector<BranchPoint3D> &b, int &ignoreCount, double thresh = 0.01){
	vector<double> tipRatios;
	ignoreCount = 0;
	for(unsigned int i(0); i < b.size(); i++){
		if(b[i].childIndex.size() > 0){
			for(unsigned int c(0); c < b[i].childIndex.size() - 1; c++){
				int tips1(numTips(b, b[i].childIndex[c]));
				double len1(separation3D(b[i], b[b[i].childIndex[c]]));
				for(unsigned int cc(c + 1); cc < b[i].childIndex.size(); cc++){
					int tips2(numTips(b, b[i].childIndex[cc]));
					double len2(separation3D(b[i], b[b[i].childIndex[cc]]));
					if(len1 >= thresh || len2 >= thresh){ // to be consistent with lambdaLs
						if(tips1 < tips2)
							tipRatios.push_back(double(tips1)/double(tips2));
						else
							tipRatios.push_back(double(tips2)/double(tips1));
					}else
						ignoreCount++;
				}
			}
		}
	}
	return tipRatios;
}

void replaceChildWith(vector<int> &childIndex, int oldChildIndex, int newChildIndex){
	for(unsigned int i(0); i < childIndex.size(); i++){
		if(childIndex[i] == oldChildIndex){
			childIndex[i] = newChildIndex;
			return;
		}
	}
	cout << "\n Error replaceChildWith(): Unable to find oldChildIndex = " << oldChildIndex
		<< " to replace with newChildIndex = " << newChildIndex
		<< endl;
}

void heirarchicalSwaps(vector<Bifurcation2D> &b, string uniquerPre = "", string uniquerPost = "", bool print = false,
		bool requireShorterLength = true, bool requireShorterAveDist = true, bool requireShorterMaxDist = true){
	bool changed(true);
	double currentLength(totalLength(b)), currentAveDist(-1.0), currentMaxDist(-1.0);
	maxAndAverageDistanceFromHeart(currentMaxDist, currentAveDist, b);
	double initLength(currentLength), initAveDist(currentAveDist), initMaxDist(currentMaxDist);
	if(print){
		cout << "\n heirarchicalSwaps(): initial currentLength = " << currentLength
			<< "\t currentAveDist = " << currentAveDist << "\t currentMaxDist = " << currentMaxDist;
	}
	vector<int> swappedAtLeastOnce;

	int *numSwaps = new int[b.size()];
	for(unsigned int i(0); i < b.size(); i++)
		numSwaps[i] = 0;

	cout << endl;
	while(changed){
		changed = false;
		bool *itSwaps = new bool[b.size()];
		for(unsigned int i(0); i < b.size(); i++)
			itSwaps[i] = false;
		for(unsigned int i(0); i < b.size(); i++){
			int pIndex(b[i].parentIndex);
			if(pIndex > -1 && b[i].childIndex.size() > 0){
				int rIndex(b[b[i].parentIndex].parentIndex);
				if(rIndex > -1){
					if(print)
						cout << "\n heirarchicalSwaps(): Bifurcation " << i << " is a candidate for swapping." << endl;
					int bestS(-1);
					double shortestLength(currentLength), bestAveDist(currentAveDist), bestMaxDist(currentMaxDist);
					for(unsigned int sEnum(0); sEnum < b[rIndex].childIndex.size(); sEnum++){
						if(b[rIndex].childIndex[sEnum] == pIndex)
							continue;
						vector<Bifurcation2D> bTemp(b);
						int sIndex(b[rIndex].childIndex[sEnum]);
						Bifurcation2D *r( &(bTemp[rIndex]) ), *p( &(bTemp[pIndex]) ), *t( &(bTemp[i]) ), *s( &(bTemp[sIndex]) );
						// t_par(p->r)
						t->parentIndex = rIndex;
						//r_child(-s, +t)
						//s_par(r->p)
						//p_child(-t, +s)
						replaceChildWith(r->childIndex, sIndex, i);
						s->parentIndex = pIndex;
						replaceChildWith(p->childIndex, i, sIndex);
						settleBifurcations(bTemp);
						//evaluate
						double tLen(totalLength(bTemp)), aDist(-1.0), mDist(-1.0);
						maxAndAverageDistanceFromHeart(mDist, aDist, bTemp);
						int betterCount(0), lastBad(-1);
						if(print){
							cout << "\n HeirarchicalSwaps(): {tLen/aDist/mDist} \t best:"
								<< "\t" << shortestLength << " \t" << bestAveDist << " \t" << bestMaxDist
								<< "\n                                         \t this:"
								<< "\t" << tLen << " \t" << aDist << " \t" << mDist
								<< endl;
							cout << "                                                   ";
							if(tLen < shortestLength){ betterCount++; cout << "\t better Length     "; }
							else{ lastBad = 0; cout << "\t                "; }
							if(aDist < bestAveDist){ betterCount++; cout << "\t better aveDist     "; }
							else{ lastBad = 1; cout << "\t                 "; }
							if(mDist <= bestMaxDist){ betterCount++; cout << "\t better maxDist     "; }
							else{ lastBad = 2; cout << "\t                 "; }
							cout << endl;
						}
						bool goodSwap(true);
						if((requireShorterLength && shortestLength < tLen)
								|| (requireShorterAveDist && bestAveDist < aDist)
								|| (requireShorterMaxDist && bestMaxDist < mDist)
								|| (!requireShorterLength && !requireShorterAveDist && !requireShorterMaxDist))
							goodSwap = false;
						//if(shortestLength > tLen && bestAveDist > aDist && bestMaxDist >= mDist){
						//cout << "\n goodSwap = " << goodSwap << endl;
						if(goodSwap){
							shortestLength = tLen;
							bestAveDist = aDist;
							bestMaxDist = mDist;
							bestS = sIndex;
							if(print)
								cout << "\t\t ** Beneficial Swap Identified **" << endl;
						}else if(print && betterCount > 1){
							string names[] = {"total_length", "average_distance", "max_distance"};
							string fn("outputs/" + uniquerPre + "heirarchicalSwaps_REJECT(" + names[lastBad] + ")_" + makeString(i) + "-" + makeString(bestS) + uniquerPost + ".png");
							vector<int> emphasis;
							emphasis.push_back(i);
							emphasis.push_back(sIndex);
							drawBifurcationTree(fn, bTemp, emphasis, 300, 300, 10);
						}
						//return s configuration -- unneed(u)d since a new bTemp will be created each time
						//replaceChildWith(r->childIndex, i, sIndex);
						//s->parentIndex = rIndex;
						//replaceChildWith(p->childIndex, sIndex, i);
					}
					//check if better configuration possible
					if(bestS > -1){
						if(print)
							cout << "\n heirarchicalSwaps(): swapping " << i << " with " << bestS << endl;
						b[i].parentIndex = rIndex;
						replaceChildWith(b[rIndex].childIndex, bestS, i);
						b[bestS].parentIndex = pIndex;
						replaceChildWith(b[pIndex].childIndex, i, bestS);
						settleBifurcations(b);
						currentLength = totalLength(b);
						maxAndAverageDistanceFromHeart(currentMaxDist, currentAveDist, b);
						vector<int> emphasis;
						emphasis.push_back(i);
						emphasis.push_back(bestS);
						string fn("outputs/" + uniquerPre + "heirarchicalSwaps_" + makeString(i) + "-" + makeString(bestS) + uniquerPost + ".png");
						if(print)
							drawBifurcationTree(fn, b, emphasis, 300, 300, 10);
						if(!containedIn(i, swappedAtLeastOnce))
							swappedAtLeastOnce.push_back(i);
						if(!containedIn(bestS, swappedAtLeastOnce))
							swappedAtLeastOnce.push_back(bestS);
						changed = true;
						numSwaps[i]++;
						numSwaps[bestS]++;
						itSwaps[i] = true;
						itSwaps[bestS] = true;
						//i--;
					}
				}
			}
		}
		for(unsigned int i(0); i < b.size(); i++){
			if(itSwaps[i])
				cout << "*";
			else
				cout << ".";
		}
		cout << endl;
		delete[] itSwaps;
	}
	if(b.size() > 0)
		cout << "\n" << numSwaps[0];
	for(unsigned int i(1); i < b.size(); i++)
		cout << "." << numSwaps[i];
	cout << endl;
	delete[] numSwaps;
	cout << "\n heirarchicalSwaps(): final currentLength = " << currentLength
		<< "\t currentAveDist = " << currentAveDist << "\t currentMaxDist = " << currentMaxDist
		<< "\n                                          (" << currentLength - initLength
		<< ")\t                (" << currentAveDist - initAveDist
													<< ")\t            (" << currentMaxDist - initMaxDist << ")"
		<< endl;
	drawBifurcationTree("outputs/" + uniquerPre + "heirarchicalSwaps_AllSwaps" + uniquerPost + ".png", b, swappedAtLeastOnce, 300, 300, 10);
}

void randomHeirarchicalSwaps(vector<Bifurcation2D> &b, string uniquerPre = "", string uniquerPost = "", bool print = false,
	bool requireShorterLength = true, bool requireShorterAveDist = true, bool requireShorterMaxDist = true){
	bool changed(true);
	double currentLength(totalLength(b)), currentAveDist(-1.0), currentMaxDist(-1.0);
	maxAndAverageDistanceFromHeart(currentMaxDist, currentAveDist, b);
	double initLength(currentLength), initAveDist(currentAveDist), initMaxDist(currentMaxDist);
	if(print){
		cout << "\n randomHeirarchicalSwaps(): initial currentLength = " << currentLength
			<< "\t currentAveDist = " << currentAveDist << "\t currentMaxDist = " << currentMaxDist;
	}
	vector<int> swappedAtLeastOnce;

	int *numSwaps = new int[b.size()];
	for(unsigned int i(0); i < b.size(); i++)
		numSwaps[i] = 0;

	cout << endl;
	while(changed){
		changed = false;
		bool *itSwaps = new bool[b.size()];
		for(unsigned int i(0); i < b.size(); i++)
			itSwaps[i] = false;
		vector<int> indices;
		for(unsigned int i(0); i < b.size(); i++)
			indices.push_back(i);
		while(indices.size() > 0){
			int randomIndexIndex(kiss()%indices.size());
			int i(indices[randomIndexIndex]);
			//cout << " " << i << " ";
			indices.erase(indices.begin() + randomIndexIndex);
			int pIndex(b[i].parentIndex);
			if(pIndex > -1 && b[i].childIndex.size() > 0){
				int rIndex(b[b[i].parentIndex].parentIndex);
				if(rIndex > -1){
					if(print)
						cout << "\n randomHeirarchicalSwaps(): Bifurcation " << i << " is a candidate for swapping." << endl;
					int bestS(-1);
					double shortestLength(currentLength), bestAveDist(currentAveDist), bestMaxDist(currentMaxDist);
					for(unsigned int sEnum(0); sEnum < b[rIndex].childIndex.size(); sEnum++){
						if(b[rIndex].childIndex[sEnum] == pIndex)
							continue;
						vector<Bifurcation2D> bTemp(b);
						int sIndex(b[rIndex].childIndex[sEnum]);
						Bifurcation2D *r(&(bTemp[rIndex])), *p(&(bTemp[pIndex])), *t(&(bTemp[i])), *s(&(bTemp[sIndex]));
						// t_par(p->r)
						t->parentIndex = rIndex;
						//r_child(-s, +t)
						//s_par(r->p)
						//p_child(-t, +s)
						replaceChildWith(r->childIndex, sIndex, i);
						s->parentIndex = pIndex;
						replaceChildWith(p->childIndex, i, sIndex);
						settleBifurcations(bTemp);
						//evaluate
						double tLen(totalLength(bTemp)), aDist(-1.0), mDist(-1.0);
						maxAndAverageDistanceFromHeart(mDist, aDist, bTemp);
						int betterCount(0), lastBad(-1);
						if(print){
							cout << "\n randomHeirarchicalSwaps(): {tLen/aDist/mDist} \t best:"
								<< "\t" << shortestLength << " \t" << bestAveDist << " \t" << bestMaxDist
								<< "\n                                         \t this:"
								<< "\t" << tLen << " \t" << aDist << " \t" << mDist
								<< endl;
							cout << "                                                   ";
							if(tLen < shortestLength){ betterCount++; cout << "\t better Length     "; } else{ lastBad = 0; cout << "\t                "; }
							if(aDist < bestAveDist){ betterCount++; cout << "\t better aveDist     "; } else{ lastBad = 1; cout << "\t                 "; }
							if(mDist <= bestMaxDist){ betterCount++; cout << "\t better maxDist     "; } else{ lastBad = 2; cout << "\t                 "; }
							cout << endl;
						}
						bool goodSwap(true);
						if((requireShorterLength && shortestLength < tLen)
							|| (requireShorterAveDist && bestAveDist < aDist)
							|| (requireShorterMaxDist && bestMaxDist < mDist)
							|| (!requireShorterLength && !requireShorterAveDist && !requireShorterMaxDist))
							goodSwap = false;
						if(goodSwap){
							shortestLength = tLen;
							bestAveDist = aDist;
							bestMaxDist = mDist;
							bestS = sIndex;
							if(print)
								cout << "\t\t ** Beneficial Swap Identified **" << endl;
						} else if(print && betterCount > 1){
							string names[] = {"total_length", "average_distance", "max_distance"};
							string fn("outputs/" + uniquerPre + "randomheirarchicalSwaps_REJECT(" + names[lastBad] + ")_" + makeString(i) + "-" + makeString(bestS) + uniquerPost + ".png");
							vector<int> emphasis;
							emphasis.push_back(i);
							emphasis.push_back(sIndex);
							drawBifurcationTree(fn, bTemp, emphasis, 300, 300, 10);
						}
					}
					//check if better configuration possible
					if(bestS > -1){
						if(print)
							cout << "\n randomHeirarchicalSwaps(): swapping " << i << " with " << bestS << endl;
						b[i].parentIndex = rIndex;
						replaceChildWith(b[rIndex].childIndex, bestS, i);
						b[bestS].parentIndex = pIndex;
						replaceChildWith(b[pIndex].childIndex, i, bestS);
						settleBifurcations(b);
						currentLength = totalLength(b);
						maxAndAverageDistanceFromHeart(currentMaxDist, currentAveDist, b);
						vector<int> emphasis;
						emphasis.push_back(i);
						emphasis.push_back(bestS);
						string fn("outputs/" + uniquerPre + "randomHeirarchicalSwaps_" + makeString(i) + "-" + makeString(bestS) + uniquerPost + ".png");
						if(print)
							drawBifurcationTree(fn, b, emphasis, 300, 300, 10);
						if(!containedIn(i, swappedAtLeastOnce))
							swappedAtLeastOnce.push_back(i);
						if(!containedIn(bestS, swappedAtLeastOnce))
							swappedAtLeastOnce.push_back(bestS);
						changed = true;
						numSwaps[i]++;
						numSwaps[bestS]++;
						itSwaps[i] = true;
						itSwaps[bestS] = true;
						//i--;
					}
				}
			}
		}
		for(unsigned int i(0); i < b.size(); i++){
			if(itSwaps[i])
				cout << "*";
			else
				cout << ".";
		}
		cout << endl;
		delete[] itSwaps;
	}
	if(b.size() > 0)
		cout << "\n" << numSwaps[0];
	for(unsigned int i(1); i < b.size(); i++)
		cout << "." << numSwaps[i];
	cout << endl;
	delete[] numSwaps;
	cout << "\n heirarchicalSwaps(): final currentLength = " << currentLength
		<< "\t currentAveDist = " << currentAveDist << "\t currentMaxDist = " << currentMaxDist
		<< "\n                                          (" << currentLength - initLength
		<< ")\t                (" << currentAveDist - initAveDist
		<< ")\t            (" << currentMaxDist - initMaxDist << ")"
		<< endl;
	drawBifurcationTree("outputs/" + uniquerPre + "heirarchicalSwaps_AllSwaps" + uniquerPost + ".png", b, swappedAtLeastOnce, 300, 300, 10);
}

void randomHeirarchicalSwapsForPartialNetwork(vector<Bifurcation2D> &b, string uniquerPre = "", string uniquerPost = "", bool print = false,
	bool requireShorterLength = true, bool requireShorterAveDist = true, bool requireShorterMaxDist = true){
	bool changed(true);
	double currentLength(totalLength(b)), currentAveDist(-1.0), currentMaxDist(-1.0);
	maxAndAverageDistanceFromHeart(currentMaxDist, currentAveDist, b);
	double initLength(currentLength), initAveDist(currentAveDist), initMaxDist(currentMaxDist);
	if(print){
		cout << "\n randomHeirarchicalSwaps(): initial currentLength = " << currentLength
			<< "\t currentAveDist = " << currentAveDist << "\t currentMaxDist = " << currentMaxDist;
	}
	vector<int> swappedAtLeastOnce;

	int *numSwaps = new int[b.size()];
	for(unsigned int i(0); i < b.size(); i++)
		numSwaps[i] = 0;

	cout << endl;
	while(changed){
		changed = false;
		bool *itSwaps = new bool[b.size()];
		for(unsigned int i(0); i < b.size(); i++)
			itSwaps[i] = false;
		vector<int> indices;
		for(unsigned int i(0); i < b.size(); i++){
			if(connectedToHeart(b, i) && b[i].parentIndex > -1 && b[i].childIndex.size() > 0)
				indices.push_back(i);
		}
		//cout << "\nindices:" << endl;
		//for(unsigned int i(0); i < indices.size(); i++)
		//	cout << indices[i] << "(parentIndex = " << b[indices[i]].parentIndex << ")" << endl;
		while(indices.size() > 0){
			int randomIndexIndex(kiss()%indices.size());
			int i(indices[randomIndexIndex]);
			//cout << " " << i << " ";
			indices.erase(indices.begin() + randomIndexIndex);
			int pIndex(b[i].parentIndex);
			if(pIndex > -1 && b[i].childIndex.size() > 0){
				int rIndex(b[b[i].parentIndex].parentIndex);
				if(rIndex > -1){
					if(print)
						cout << "\n randomHeirarchicalSwaps(): Bifurcation " << i << " is a candidate for swapping." << endl;
					int bestS(-1);
					double shortestLength(currentLength), bestAveDist(currentAveDist), bestMaxDist(currentMaxDist);
					for(unsigned int sEnum(0); sEnum < b[rIndex].childIndex.size(); sEnum++){
						if(b[rIndex].childIndex[sEnum] == pIndex)
							continue;
						vector<Bifurcation2D> bTemp(b);
						int sIndex(b[rIndex].childIndex[sEnum]);
						Bifurcation2D *r(&(bTemp[rIndex])), *p(&(bTemp[pIndex])), *t(&(bTemp[i])), *s(&(bTemp[sIndex]));
						// t_par(p->r)
						t->parentIndex = rIndex;
						//r_child(-s, +t)
						//s_par(r->p)
						//p_child(-t, +s)
						replaceChildWith(r->childIndex, sIndex, i);
						s->parentIndex = pIndex;
						replaceChildWith(p->childIndex, i, sIndex);
						settleBifurcations(bTemp);
						//evaluate
						double tLen(totalLength(bTemp)), aDist(-1.0), mDist(-1.0);
						maxAndAverageDistanceFromHeart(mDist, aDist, bTemp);
						int betterCount(0), lastBad(-1);
						if(print){
							cout << "\n randomHeirarchicalSwaps(): {tLen/aDist/mDist} \t best:"
								<< "\t" << shortestLength << " \t" << bestAveDist << " \t" << bestMaxDist
								<< "\n                                         \t this:"
								<< "\t" << tLen << " \t" << aDist << " \t" << mDist
								<< endl;
							cout << "                                                   ";
							if(tLen < shortestLength){ betterCount++; cout << "\t better Length     "; } else{ lastBad = 0; cout << "\t                "; }
							if(aDist < bestAveDist){ betterCount++; cout << "\t better aveDist     "; } else{ lastBad = 1; cout << "\t                 "; }
							if(mDist <= bestMaxDist){ betterCount++; cout << "\t better maxDist     "; } else{ lastBad = 2; cout << "\t                 "; }
							cout << endl;
						}
						bool goodSwap(true);
						if((requireShorterLength && shortestLength <= tLen)
							|| (requireShorterAveDist && bestAveDist <= aDist)
							|| (requireShorterMaxDist && bestMaxDist < mDist)
							|| (!requireShorterLength && !requireShorterAveDist && !requireShorterMaxDist))
							goodSwap = false;
						if(goodSwap){
							shortestLength = tLen;
							bestAveDist = aDist;
							bestMaxDist = mDist;
							bestS = sIndex;
							if(print)
								cout << "\t\t ** Beneficial Swap Identified **" << endl;
						} else if(print && betterCount > 1){
							string names[] = {"total_length", "average_distance", "max_distance"};
							string fn("outputs/" + uniquerPre + "randomheirarchicalSwaps_REJECT(" + names[lastBad] + ")_" + makeString(i) + "-" + makeString(bestS) + uniquerPost + ".png");
							vector<int> emphasis;
							emphasis.push_back(i);
							emphasis.push_back(sIndex);
							drawBifurcationTree(fn, bTemp, emphasis, 300, 300, 10);
						}
					}
					//check if better configuration possible
					if(bestS > -1){
						if(print)
							cout << "\n randomHeirarchicalSwaps(): swapping " << i << " with " << bestS << endl;
						b[i].parentIndex = rIndex;
						replaceChildWith(b[rIndex].childIndex, bestS, i);
						b[bestS].parentIndex = pIndex;
						replaceChildWith(b[pIndex].childIndex, i, bestS);
						settleBifurcations(b);
						currentLength = totalLength(b);
						maxAndAverageDistanceFromHeart(currentMaxDist, currentAveDist, b);
						vector<int> emphasis;
						emphasis.push_back(i);
						emphasis.push_back(bestS);
						string fn("outputs/" + uniquerPre + "randomHeirarchicalSwaps_" + makeString(i) + "-" + makeString(bestS) + uniquerPost + ".png");
						if(print)
							drawBifurcationTree(fn, b, emphasis, 300, 300, 10);
						if(!containedIn(i, swappedAtLeastOnce))
							swappedAtLeastOnce.push_back(i);
						if(!containedIn(bestS, swappedAtLeastOnce))
							swappedAtLeastOnce.push_back(bestS);
						changed = true;
						numSwaps[i]++;
						numSwaps[bestS]++;
						itSwaps[i] = true;
						itSwaps[bestS] = true;
						//i--;
					}
				}
			}
		}
		for(unsigned int i(0); i < b.size(); i++){
			if(itSwaps[i])
				cout << "*";
			else
				cout << ".";
		}
		cout << endl;
		delete[] itSwaps;
	}
	if(b.size() > 0)
		cout << "\n" << numSwaps[0];
	for(unsigned int i(1); i < b.size(); i++)
		cout << "." << numSwaps[i];
	cout << endl;
	delete[] numSwaps;
	cout << "\n heirarchicalSwaps(): final currentLength = " << currentLength
		<< "\t currentAveDist = " << currentAveDist << "\t currentMaxDist = " << currentMaxDist
		<< "\n                                          (" << currentLength - initLength
		<< ")\t                (" << currentAveDist - initAveDist
		<< ")\t            (" << currentMaxDist - initMaxDist << ")"
		<< endl;
	drawBifurcationTree("outputs/" + uniquerPre + "heirarchicalSwaps_AllSwaps" + uniquerPost + ".png", b, swappedAtLeastOnce, 300, 300, 10);
}

void testHeirarchicalSwaps(){
	kisset(24611, 19046, 7706);
	int NX(10), NY(NX), heartx(1*NX/2), hearty(0*NY/2);
	//double thresh(0.01);
	double noiseFactor(0.3);
	int width(NX*100), height(NY*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300;
	border = 10;
	cout << "\n\n NX = " << NX << "\n NY = " << NY << "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n noiseFactor = " << noiseFactor << "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<Bifurcation2D> b = initializeGridWithNoise(NX, NY, heartx, hearty, noiseFactor);
	drawBifurcationTree("outputs/testHeirarchicalSwapsInitial.png", b, width, height, border);
	introduceBifurcationsFromSources(b);
	drawBifurcationTree("outputs/testHeirarchicalSwapsTreeBuilt.png", b, width, height, border);
	heirarchicalSwaps(b, "", "", true);
	drawBifurcationTree("outputs/testHeirarchicalSwapsAfterSwapping.png", b, width, height, border);
}

void testIntroduceBifurcationsSpeedup(){
	kisset(24611, 19046, 7706);
	int startN(3), endN(20), incN(1), numN(10);
	//double thresh(0.01);
	double noiseFactor(0.3);
	//int width(NX*100), height(NY*100), border(min(width, height)/25);
	//width = height = 300;
	//border = 10;
	cout << "\n\n startN = " << startN << "\n endN = " << endN << "\n incN = " << incN << "\n numN = " << numN
		//<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n noiseFactor = " << noiseFactor
		//<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	cout << "\n\n N\t <T>_new\t <T>_old" << endl;
	for(int N(startN); N <= endN; N += incN){
		int NX(N), NY(N), heartx(1*NX/2), hearty(0*NY/2);
		cout << N;
		kisset(24611, 19046, 7706, false);
		time_t st(time(NULL));
		for(int n(0); n < numN; n++){
			vector<Bifurcation2D> b = initializeGridWithNoise(NX, NY, heartx, hearty, noiseFactor);
			introduceBifurcationsFromSources(b);
		}
		cout << "\t" << difftime(time(NULL), st)/numN;
		kisset(24611, 19046, 7706, false);
		st = time(NULL);
		for(int n(0); n < numN; n++){
			vector<Bifurcation2D> b = initializeGridWithNoise(NX, NY, heartx, hearty, noiseFactor);
			introduceBifurcationsFromSourcesSLOW(b);
		}
		cout << "\t" << difftime(time(NULL), st)/numN;
		cout << endl;
	}
}

void heirarchicalSwapAnalysis(){
	kisset(24611, 19046, 7706);
	int numBins(16), NX(20), NY(NX), heartx(1*NX/2), hearty(0*NY/2);
	double thresh(0.01), noiseFactor(0.3);
	int width(NX*100), height(NY*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300;
	border = 10;

	cout << "\n heirarchicalSwapAnalysis()\n numBins = " << numBins << "\n NX = " << NX << "\n NY = NY"
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n thresh = " << thresh << "\n noiseFactor = " << noiseFactor
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;

	vector<Bifurcation2D> b = initializeGridWithNoise(NX, NY, heartx, hearty, noiseFactor);
	drawBifurcationTree("outputs/testHeirarchicalSwapAnalysisInitial.png", b, width, height, border);
	introduceBifurcationsFromSources(b, false);
	int ignoredTree(0);
	vector<double> lambdaLsTree = validBranchingLengthRatios(b, ignoredTree, thresh);
	drawBifurcationTree("outputs/testHeirarchicalSwapAnalysisTree.png", b, width, height, border);
	heirarchicalSwaps(b, "", "", false);
	int ignoredSwap(0);
	vector<double> lambdaLsSwap = validBranchingLengthRatios(b, ignoredSwap, thresh);
	drawBifurcationTree("outputs/testHeirarchicalSwapAnalysisSwapped.png", b, width, height, border);
	
	cout << "\n Tree ignored " << 100.0*ignoredTree/b.size() << "% of segments"
		<< "\n Swap ignored " << 100.0*ignoredSwap/b.size() << "% of segments"
		<< endl;

	cout << "\nnumBins = " << numBins;
	map<double, double> treeDensity = probabilityDistribution(lambdaLsTree, numBins);
	printMap(treeDensity, "tree");
	map<double, double> swapDensity = probabilityDistribution(lambdaLsSwap, numBins);
	printMap(treeDensity, "swap");

	numBins = 3*numBins/2;
	cout << "\nnumBins = " << numBins;
	treeDensity = probabilityDistribution(lambdaLsTree, numBins);
	printMap(treeDensity, "tree");
	swapDensity = probabilityDistribution(lambdaLsSwap, numBins);
	printMap(treeDensity, "swap");

	numBins = 3*numBins/2;
	cout << "\nnumBins = " << numBins;
	treeDensity = probabilityDistribution(lambdaLsTree, numBins);
	printMap(treeDensity, "tree");
	swapDensity = probabilityDistribution(lambdaLsSwap, numBins);
	printMap(treeDensity, "swap");
}

//if root, zero levels between; if root is parent, one level between; etc . . .
unsigned int levelsBetweenRoot(int index, vector<Bifurcation2D> &b){
	unsigned int levelsBetween(0);
	while(b[index].parentIndex > -1){
		index = b[index].parentIndex;
		levelsBetween++;
	}
	return levelsBetween;
}

//root is level zero (0)
unsigned int* levelsFromHeart(vector<Bifurcation2D> &b){
	unsigned int *levels = new unsigned int[b.size()];
	for(unsigned int i(0); i < b.size(); i++)
		levels[i] = levelsBetweenRoot(i, b);
	return levels;
}

double dotProduct(int a1, int a2, int b1, int b2, vector<Bifurcation2D> &b){
	double x1(b[a2].position[0] - b[a1].position[0]), y1(b[a2].position[1] - b[a1].position[1]),
		x2(b[b2].position[0] - b[b1].position[0]), y2(b[b2].position[1] - b[b1].position[1]);
	return x1*x2 + y1*y2;
}

double dotProduct(int a1, int a2, int b1, int b2, vector<Bifurcation3D> &b){
	double x1(b[a2].position[0] - b[a1].position[0]), y1(b[a2].position[1] - b[a1].position[1]), z1(b[a2].position[2] - b[a1].position[2]),
		x2(b[b2].position[0] - b[b1].position[0]), y2(b[b2].position[1] - b[b1].position[1]), z2(b[b2].position[2] - b[b1].position[2]);
	return x1*x2 + y1*y2 + z1*z2;
}

double angleBetween(int parentIndex, int childIndex, int grandchildIndex, vector<Bifurcation2D> &b){
	double lenPar(separation2D(b[parentIndex], b[childIndex])), lenChild(separation2D(b[childIndex], b[grandchildIndex]));
	return acos(dotProduct(parentIndex, childIndex, childIndex, grandchildIndex, b)/lenPar/lenChild);
}

double angleBetweenGrandchildren(int grandIndex, int child1Index, int child2Index, vector<Bifurcation2D> &b){
	double lenPar(separation2D(b[grandIndex], b[child1Index])), lenChild(separation2D(b[grandIndex], b[child2Index]));
	return acos(dotProduct(grandIndex, child1Index, grandIndex, child2Index, b)/lenPar/lenChild);
}

double angleBetweenGrandchildren3D(int grandIndex, int child1Index, int child2Index, vector<Bifurcation3D> &b){
	double lenPar(separation3D(b[grandIndex], b[child1Index])), lenChild(separation3D(b[grandIndex], b[child2Index]));
	return acos(dotProduct(grandIndex, child1Index, grandIndex, child2Index, b)/lenPar/lenChild);
}

unsigned int numTips(int index, vector<Bifurcation2D> &b){
	unsigned int tips(0);
	if(b[index].childIndex.size() < 1)
		return 1;
	for(unsigned int c(0); c < b[index].childIndex.size(); c++)
		tips += numTips(b[index].childIndex[c], b);
	return tips;
}

void generalAnalysis(){
	kisset(24611, 19046, 7706);
	int numBins(25), numGrid(3), NX(10), NY(NX), heartx(1*NX/2), hearty(0*NY/2);
	string uniquer("gaTEST" + makeString(NX));
	double thresh(0.001), noiseFactor(0.3);
	int width(NX*100), height(NY*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300;
	border = 10;

	cout << "\n heirarchicalSwapAnalysis()\n numBins = " << numBins << "\n NX = " << NX << "\n NY = " << NY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n thresh = " << thresh << "\n noiseFactor = " << noiseFactor
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;

	vector<Bifurcation2D> b = initializeGridWithNoise(NX, NY, heartx, hearty, noiseFactor);
	//drawBifurcationTree("outputs/" + uniquer + "generalAnalysisGrid.png", b, width, height, border);
	introduceBifurcationsFromSources(b);
	drawBifurcationTree("outputs/" + uniquer + "generalAnalysisInitialBuild.png", b, width, height, border);
	for(int n(0); n < numGrid; n++){
		vector<Bifurcation2D> B(b);
		randomHeirarchicalSwaps(B);
		double aveDist(0.0), maxDist(0.0);
		maxAndAverageDistanceFromHeart(maxDist, aveDist, B);
		cout << "\n n = " << n << "\n total length = " << totalLength(B)
			<< "\n average distance = " << aveDist
			<< "\n max distance = " << maxDist << endl;
		drawBifurcationTree("outputs/" + uniquer + "generalAnalysisSwapped" + makeString(n) + ".png", B, width, height, border);
		int ignored(0);
		vector<double> lambdaLs = validBranchingLengthRatios(B, ignored, thresh);
		cout << "\n ignored " << 100.0*ignored/B.size() << "% of segments for lambdaLs" << endl;

		cout << "\nnumBins = " << numBins;
		map<double, double> lambdaLsDensity = probabilityDistribution(lambdaLs, numBins);
		printMap(lambdaLsDensity, "lambdaLs");
		vector<double> gammas;
		vector<double> angles;
		vector<unsigned int> childTips;
		unsigned int *levels = levelsFromHeart(B);
		unsigned int maxLvl(maxLevel(B, 0));
		vector<double> *gammasLevel = new vector<double>[maxLvl - 2];
		vector<double> *anglesLevel = new vector<double>[maxLvl - 2];
		for(unsigned int i(0); i < B.size(); i++){
			for(unsigned int c(0); c < B[i].childIndex.size(); c++){
				int childIndex(B[i].childIndex[c]);
				double parentLength(separation2D(B[i], B[childIndex]));
				for(unsigned int cc(0); cc < B[childIndex].childIndex.size(); cc++){
					int grandchildIndex(B[childIndex].childIndex[cc]);
					double childLength(separation2D(B[childIndex], B[grandchildIndex]));
					if(parentLength >= thresh && childLength >= thresh){
						double gamma(childLength/parentLength);
						double angle(angleBetween(i, childIndex, grandchildIndex, B));
						unsigned int tips(numTips(childIndex, B));
						gammas.push_back(gamma);
						angles.push_back(angle);
						childTips.push_back(tips);
						gammasLevel[levels[i]].push_back(gamma);
						anglesLevel[levels[i]].push_back(angle);
					}
				}
			}
		}
		// relation of angle with level and with gamma; ratio of child angles with ratio of number of tips of children
		cout << "\n<gamma> = " << mean(gammas)
			<< "\n stDev[angle] = " << stDev(angles)
			<< endl;
		map<double, double> gammasDensity = probabilityDistribution(gammas, numBins);
		printMap(gammasDensity, "gammas");
		cout << "\n lvl\t <gamma(level)> \t\t #[gammas(lvl)] \t\t stDev[angle(lvl)]" << endl;
		string gammasLevelFn("outputs/" + uniquer + "_" + makeString(n) + "_" + "levelGamma.dat");
		ofstream outGammasLevel(gammasLevelFn.c_str());
		string gammasAngleFn("outputs/" + uniquer + "_" + makeString(n) + "_" +  +"gammaAngle.dat");
		ofstream outGammasAngle(gammasAngleFn.c_str());
		for(unsigned int lvl(0); lvl < maxLvl - 2; lvl++){
			cout << lvl << "\t" << mean(gammasLevel[lvl]) << "\t\t" << gammasLevel[lvl].size()
				<< "\t\t" << stDev(anglesLevel[lvl]) << endl;
			for(unsigned int i(0); i < gammasLevel[lvl].size(); i++){
				outGammasLevel << lvl << "\t" << gammasLevel[lvl][i] << endl;
				outGammasAngle << gammasLevel[lvl][i] << "\t" << anglesLevel[lvl][i] << endl;
			}
			outGammasLevel << endl;
			outGammasAngle << endl;
		}
		outGammasLevel.close();
		outGammasAngle.close();
		string tipsAngleFn("outputs/" + uniquer + "_" + makeString(n) + "_" +  +"tipsAngle.dat");
		ofstream outTipsAngle(tipsAngleFn.c_str());
		for(unsigned int i(0); i < childTips.size(); i++)
			outTipsAngle << childTips[i] << "\t" << angles[i] << endl;
		outTipsAngle.close();
		delete[] levels;
	}
}

void testSwapConditions(){
	kisset(24611, 19046, 7706);
	int numGrid(5), NX(5), NY(NX), heartx(1*NX/2), hearty(0*NY/2);
	string uniquer("tsc_noTipSwap" + makeString(NX));
	//double thresh(0.01);
	double noiseFactor(0.3);
	int width(NX*100), height(NY*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300;
	border = 10;
	cout << "\n\n numGrid = " << numGrid << "\n NX = " << NX << "\n NY = " << NY << "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n noiseFactor = " << noiseFactor << "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	double *initialLength = new double[numGrid];
	double *initialAveDist = new double[numGrid];
	double *initialMaxDist = new double[numGrid];
	double **finalLength = new double*[numGrid];
	double **finalAveDist = new double*[numGrid];
	double **finalMaxDist = new double*[numGrid];
	for(int n(0); n < numGrid; n++){
		finalLength[n] = new double[7];
		finalAveDist[n] = new double[7];
		finalMaxDist[n] = new double[7];
		vector<Bifurcation2D> b = initializeGridWithNoise(NX, NY, heartx, hearty, noiseFactor);
		drawBifurcationTree("outputs/" + uniquer + "_" + makeString(n) + "_testSwapConditionsInitial.png", b, width, height, border);
		introduceBifurcationsFromSources(b);
		initialLength[n] = totalLength(b);
		maxAndAverageDistanceFromHeart(initialMaxDist[n], initialAveDist[n], b);
		drawBifurcationTree("outputs/" + uniquer + "_" + makeString(n) + "_testSwapConditionsTreeBuilt.png", b, width, height, border);
		for(int i(1); i < 8; i++){
			vector<Bifurcation2D> B(b);
			string combo(makeString(i/4) + makeString((i/2)%2) + makeString(i%2));
			heirarchicalSwaps(B, uniquer + "_" + makeString(n) + "_", "_" + combo, false, (i/4) > 0, (i/2)%2 > 0, i%2 > 0);
			finalLength[n][i - 1] = totalLength(B);
			maxAndAverageDistanceFromHeart(finalMaxDist[n][i - 1], finalAveDist[n][i - 1], B);
			drawBifurcationTree("outputs/" + uniquer + "_" + makeString(n) + "_testSwapConditionsAfterSwapping" + combo + ".png", B, width, height, border);
		}
	}

	cout << "\n\n total length (4's {leftmost} digit)\n init";
	for(int i(1); i < 8; i++)
		cout << "\t  " << makeString(i/4) + makeString((i/2)%2) + makeString(i%2);
	cout << endl;
	for(int n(0); n < numGrid; n++){
		cout << initialLength[n];
		for(int i(1); i < 8; i++)
			cout << "\t" << finalLength[n][i - 1];
		cout << endl;
	}

	cout << "\n\n average distance from heart (2's {middle} digit)\n init";
	for(int i(1); i < 8; i++)
		cout << "\t  " << makeString(i/4) + makeString((i/2)%2) + makeString(i%2);
	cout << endl;
	for(int n(0); n < numGrid; n++){
		cout << initialAveDist[n];
		for(int i(1); i < 8; i++)
			cout << "\t" << finalAveDist[n][i - 1];
		cout << endl;
	}

	cout << "\n\n maximum distance from heart (1's {rightmost} digit)\n init";
	for(int i(1); i < 8; i++)
		cout << "\t  " << makeString(i/4) + makeString((i/2)%2) + makeString(i%2);
	cout << endl;
	for(int n(0); n < numGrid; n++){
		cout << initialMaxDist[n];
		for(int i(1); i < 8; i++)
			cout << "\t" << finalMaxDist[n][i - 1];
		cout << endl;
	}

	for(int n(0); n < numGrid; n++){
		delete[] finalLength[n];
		delete[] finalAveDist[n];
		delete[] finalMaxDist[n];
	}
	delete[] initialLength;
	delete[] initialAveDist;
	delete[] initialMaxDist;
	delete[] finalLength;
	delete[] finalAveDist;
	delete[] finalMaxDist;
}

void swappingSequence(){
	kisset(24611, 19046, 7706);
	int numGrid(500), NX(8), NY(NX), heartx(1*NX/2), hearty(0*NY/2);
	string uniquer("ss_" + makeString(NX));
	//double thresh(0.01);
	double noiseFactor(0.3);
	int width(NX*100), height(NY*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300;
	border = 10;
	cout << "\n\n numGrid = " << numGrid << "\n NX = " << NX << "\n NY = " << NY << "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n noiseFactor = " << noiseFactor << "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<Bifurcation2D> b = initializeGridWithNoise(NX, NY, heartx, hearty, noiseFactor);
	introduceBifurcationsFromSources(b);
	drawBifurcationTree("outputs/" + uniquer + "_swappingSequenceTreeBuilt.png", b, width, height, border);
	double *lengths = new double[numGrid];
	double *aveDists = new double[numGrid];
	double *maxDists = new double[numGrid];
	
	vector<Bifurcation2D> *B = new vector<Bifurcation2D>[numGrid];
	for(int n(0); n < numGrid; n++){
		B[n] = b;
		randomHeirarchicalSwaps(B[n], uniquer + "_" + makeString(n) + "_");
		lengths[n] = totalLength(B[n]);
		maxAndAverageDistanceFromHeart(maxDists[n], aveDists[n], B[n]);
	}

	int *orderedGrids = new int[numGrid];
	for(int n(0); n < numGrid; n++)
		orderedGrids[n] = n;
	bool changed(true);
	while(changed){
		changed = false;
		for(int n(0); n < numGrid - 1; n++){
			int i(orderedGrids[n]);
			int j(orderedGrids[n + 1]);
			double iMetric(lengths[i] + aveDists[i] + maxDists[i]),
				jMetric(lengths[j] + aveDists[j] + maxDists[j]);
			if(iMetric > jMetric){
				orderedGrids[n] = j;
				orderedGrids[n + 1] = i;
				changed = true;
			}
		}
	}

	double initAveDist(0.0), initMaxDist(0.0);
	maxAndAverageDistanceFromHeart(initMaxDist, initAveDist, b);
	cout << "\n\n #\t len\t aveDist\t maxDist\n"
		<< "i \t" << totalLength(b) << "\t" << initAveDist << "\t" << initMaxDist
		<< "\n -- \t -------------- \t -------------- \t -------------- "
		<< endl;
	for(int n(0); n < numGrid; n++){
		cout << lengths[orderedGrids[n]] << "\t" << aveDists[orderedGrids[n]] << "\t" << maxDists[orderedGrids[n]] << endl;
		drawBifurcationTree("outputs/" + uniquer + "_rankedSwaps" + makeString(n) + ".png", B[n], width, height, border);
	}

	delete[] orderedGrids;
	delete[] B;
	delete[] lengths;
	delete[] aveDists;
	delete[] maxDists;
}

void exploreSwapsAndChoose(vector<Bifurcation2D> &b, int numSwapGrids = 1, bool print = false){
	double initialLength(totalLength(b)), initialAveDist(0.0), initialMaxDist(0.0);
	maxAndAverageDistanceFromHeart(initialMaxDist, initialAveDist, b);
	vector<Bifurcation2D> bestB(b);
	bool foundBetterSwap(false);
	double bestLength(initialLength), bestAveDist(initialAveDist), bestMaxDist(initialMaxDist);
	for(int sg(0); sg < numSwapGrids; sg++){
		vector<Bifurcation2D> B(b);
		randomHeirarchicalSwapsForPartialNetwork(B);
		double length(totalLength(B)), aveDist(0.0), maxDist(0.0);
		maxAndAverageDistanceFromHeart(maxDist, aveDist, B);
		if(bestLength > length && bestAveDist > aveDist && bestMaxDist >= maxDist){
			bestLength = length;
			bestAveDist = aveDist;
			bestMaxDist = maxDist;
			bestB = B;
			foundBetterSwap = true;
		}
	}
	if(foundBetterSwap)
		b = bestB;
	if(print){
		cout << "\n exploreSwapsAndChoose():" << endl;
		if(foundBetterSwap)
			cout << "\t Found a better configuration through swapping.";
		else
			cout << "\t No better configuration found through swapping.";
		cout << "                                        \t       length      \t      aveDist      \t      maxDist"
			<<  "\n                            initial  \t " << initialLength << "\t" << initialAveDist << "\t" << initialMaxDist
			<<  "\n                            best     \t " << bestLength << "\t" << bestAveDist << "\t" << bestMaxDist << endl;
	}
}

// assumes b[0] is heart and that there are no connections present yet
void introduceBifurcationsGrowingFromConnected(vector<Bifurcation2D> &b, int numSwapGrids = 1){
	vector<int> notConnected;
	int *nearestConnection = new int[b.size()];
	double *nearestConnectionSep = new double[b.size()];
	nearestConnection[0] = 0;
	nearestConnectionSep[0] = 0.0;
	for(unsigned int i(1); i < b.size(); i++){
		notConnected.push_back(i);
		nearestConnectionSep[i] = separation2D(b[0], b[i]);
		nearestConnection[i] = 0;
	}
	int notConnectedIndex(0);
	for(unsigned int i(1); i < notConnected.size(); i++){
		if(nearestConnectionSep[notConnected[notConnectedIndex]] > nearestConnectionSep[notConnected[i]]){
			notConnectedIndex = i;
		}
	}
	b[0].childIndex.push_back(notConnected[notConnectedIndex]);
	b[notConnected[notConnectedIndex]].parentIndex = 0;
	settleBifurcations(b);
	for(unsigned int i(0); i < notConnected.size(); i++){
		double sep(separation2D(b[notConnected[notConnectedIndex]], b[notConnected[i]]));
		if(nearestConnectionSep[notConnected[i]] > sep){
			nearestConnectionSep[notConnected[i]] = sep;
			nearestConnection[notConnected[i]] = notConnected[notConnectedIndex];
		}
	}
	notConnected.erase(notConnected.begin() + notConnectedIndex);
	while(notConnected.size() > 0){
		drawBifurcationTree("outputs/introduceBifurcationsGrowingFromHeart" + makeString(notConnected.size()) + ".png", b, 300, 300, 10);
		int notConnectedIndex(0);
		for(unsigned int i(1); i < notConnected.size(); i++){
			if(nearestConnectionSep[notConnected[notConnectedIndex]] > nearestConnectionSep[notConnected[i]]){
				notConnectedIndex = i;
			}
		}
		int chosenC(-1);
		if(b[nearestConnection[notConnected[notConnectedIndex]]].parentIndex > -1){
			//cout << " adding bifurcation between nearest and parent . . ." << endl;
			insertBifurcation(b, b[nearestConnection[notConnected[notConnectedIndex]]].parentIndex, nearestConnection[notConnected[notConnectedIndex]], notConnected[notConnectedIndex]);
		}else if(b[nearestConnection[notConnected[notConnectedIndex]]].childIndex.size() > 0){
			double nearestChildSep(separation2D(b[b[nearestConnection[notConnected[notConnectedIndex]]].childIndex[0]], b[notConnected[notConnectedIndex]]));
			chosenC = 0;
			for(unsigned int c(1); c < b[nearestConnection[notConnectedIndex]].childIndex.size(); c++){
				double sep(separation2D(b[b[nearestConnection[notConnectedIndex]].childIndex[c]], b[notConnected[notConnectedIndex]]));
				if(nearestChildSep > sep){
					nearestChildSep = sep;
					chosenC = c;
				}
			}
			//cout << "\n adding bifurcation betwen nearest and child . . ." << endl;
			insertBifurcation(b, nearestConnection[notConnected[notConnectedIndex]], b[nearestConnection[notConnected[notConnectedIndex]]].childIndex[chosenC], notConnected[notConnectedIndex]);
		}else{
			cout << "\n Error: nearestConnection that is notConnected (notConnectedIndex) has no family." << endl;
			abort();
		}
			
		settleBifurcations(b);
		for(unsigned int i(0); i < notConnected.size(); i++){
			double sep(separation2D(b[notConnected[notConnectedIndex]], b[notConnected[i]]));
			if(nearestConnectionSep[notConnected[i]] > sep){
				nearestConnectionSep[notConnected[i]] = sep;
				nearestConnection[notConnected[i]] = notConnected[notConnectedIndex];
			}
		}
		exploreSwapsAndChoose(b, numSwapGrids, true);
		notConnected.erase(notConnected.begin() + notConnectedIndex);
	}
}

double tripletMeasure(vector<BranchPoint2D> &b, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient){
	double totLen(-1.0), avePath(-1.0), maxPath(-1.0);
	if(totalLengthCoefficient != 0.0)
		totLen = totalLength(b);
	if(avePathLengthCoefficient != 0.0 || maxPathLengthCoefficient != 0.0)
		maxAndAverageDistanceFromHeart(maxPath, avePath, b);
	return totalLengthCoefficient*totLen + avePathLengthCoefficient*avePath + maxPathLengthCoefficient*maxPath;
}

double tripletMeasure(vector<BranchPoint3D> &b, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient){
	double totLen(-1.0), avePath(-1.0), maxPath(-1.0);
	if(totalLengthCoefficient != 0.0)
		totLen = totalLength(b);
	if(avePathLengthCoefficient != 0.0 || maxPathLengthCoefficient != 0.0)
		maxAndAverageDistanceFromHeart(maxPath, avePath, b);
	return totalLengthCoefficient*totLen + avePathLengthCoefficient*avePath + maxPathLengthCoefficient*maxPath;
}

// assumes b[0] is heart and that there are no connections present yet
void introduceBifurcationsGrowingFromHeart(vector<Bifurcation2D> &b, int numSwapGrids = 1){
	vector<int> notConnected;
	int *nearestConnection = new int[b.size()];
	double *nearestConnectionSep = new double[b.size()];
	nearestConnection[0] = 0;
	nearestConnectionSep[0] = 0.0;
	for(unsigned int i(1); i < b.size(); i++){
		notConnected.push_back(i);
		nearestConnectionSep[i] = separation2D(b[0], b[i]);
		nearestConnection[i] = 0;
	}
	int notConnectedIndex(0);
	for(unsigned int i(1); i < notConnected.size(); i++){
		if(nearestConnectionSep[notConnected[notConnectedIndex]] > nearestConnectionSep[notConnected[i]]){
			notConnectedIndex = i;
		}
	}
	b[0].childIndex.push_back(notConnected[notConnectedIndex]);
	b[notConnected[notConnectedIndex]].parentIndex = 0;
	settleBifurcations(b);
	for(unsigned int i(0); i < notConnected.size(); i++){
		double sep(separation2D(b[notConnected[notConnectedIndex]], b[notConnected[i]]));
		if(nearestConnectionSep[notConnected[i]] > sep){
			nearestConnectionSep[notConnected[i]] = sep;
			nearestConnection[notConnected[i]] = notConnected[notConnectedIndex];
		}
	}
	notConnected.erase(notConnected.begin() + notConnectedIndex);
	while(notConnected.size() > 0){
		drawBifurcationTree("outputs/introduceBifurcationsGrowingFromHeart" + makeString(notConnected.size()) + ".png", b, 300, 300, 10);
		int notConnectedIndex(0);
		for(unsigned int i(1); i < notConnected.size(); i++){
			double nearestSepFromHeart(separation2D(b[0], b[notConnected[notConnectedIndex]])),
				sepFromHeart(separation2D(b[0], b[notConnected[i]]));
			if(nearestSepFromHeart > sepFromHeart){
				notConnectedIndex = i;
			}
		}
		int chosenC(-1);
		if(b[nearestConnection[notConnected[notConnectedIndex]]].parentIndex > -1){
			//cout << " adding bifurcation between nearest and parent . . ." << endl;
			insertBifurcation(b, b[nearestConnection[notConnected[notConnectedIndex]]].parentIndex, nearestConnection[notConnected[notConnectedIndex]], notConnected[notConnectedIndex]);
		} else if(b[nearestConnection[notConnected[notConnectedIndex]]].childIndex.size() > 0){
			double nearestChildSep(separation2D(b[b[nearestConnection[notConnected[notConnectedIndex]]].childIndex[0]], b[notConnected[notConnectedIndex]]));
			chosenC = 0;
			for(unsigned int c(1); c < b[nearestConnection[notConnectedIndex]].childIndex.size(); c++){
				double sep(separation2D(b[b[nearestConnection[notConnectedIndex]].childIndex[c]], b[notConnected[notConnectedIndex]]));
				if(nearestChildSep > sep){
					nearestChildSep = sep;
					chosenC = c;
				}
			}
			//cout << "\n adding bifurcation betwen nearest and child . . ." << endl;
			insertBifurcation(b, nearestConnection[notConnected[notConnectedIndex]], b[nearestConnection[notConnected[notConnectedIndex]]].childIndex[chosenC], notConnected[notConnectedIndex]);
		} else{
			cout << "\n Error: nearestConnection that is notConnected (notConnectedIndex) has no family." << endl;
			abort();
		}

		settleBifurcations(b);
		for(unsigned int i(0); i < notConnected.size(); i++){
			double sep(separation2D(b[notConnected[notConnectedIndex]], b[notConnected[i]]));
			if(nearestConnectionSep[notConnected[i]] > sep){
				nearestConnectionSep[notConnected[i]] = sep;
				nearestConnection[notConnected[i]] = notConnected[notConnectedIndex];
			}
		}
		exploreSwapsAndChoose(b, numSwapGrids, true);
		notConnected.erase(notConnected.begin() + notConnectedIndex);
	}
}

void gradualGrowthFromConnectedWithSwaps(){
	kisset(24611, 19046, 7706);
	int numSwapGrids(10), NX(10), NY(NX), heartx(1*NX/2), hearty(0*NY/2);
	string uniquer("ggfcws" + makeString(NX));
	double thresh(0.001), noiseFactor(0.3);
	int width(NX*100), height(NY*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300;
	border = 10;

	cout << "\n heirarchicalSwapAnalysis()\n uniquer = " << uniquer << "\n NX = " << NX << "\n NY = " << NY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n thresh = " << thresh << "\n noiseFactor = " << noiseFactor
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;

	vector<Bifurcation2D> b = initializeGridWithNoise(NX, NY, heartx, hearty, noiseFactor);
	drawBifurcationTree("outputs/" + uniquer + "gradualGrowthWithSwapsGrid.png", b, width, height, border);
	introduceBifurcationsGrowingFromConnected(b, numSwapGrids);
	drawBifurcationTree("outputs/" + uniquer + "gradualGrowthWithSwapsTreeBuilt.png", b, width, height, border);
}

void gradualGrowthFromHeartWithSwaps(){
	kisset(24611, 19046, 7706);
	int numSwapGrids(10), NX(10), NY(NX), heartx(1*NX/2), hearty(0*NY/2);
	string uniquer("ggfhws" + makeString(NX));
	double thresh(0.001), noiseFactor(0.3);
	int width(NX*100), height(NY*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300;
	border = 10;

	cout << "\n heirarchicalSwapAnalysis()\n uniquer = " << uniquer << "\n NX = " << NX << "\n NY = " << NY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n thresh = " << thresh << "\n noiseFactor = " << noiseFactor
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;

	vector<Bifurcation2D> b = initializeGridWithNoise(NX, NY, heartx, hearty, noiseFactor);
	drawBifurcationTree("outputs/" + uniquer + "gradualGrowthWithSwapsGrid.png", b, width, height, border);
	introduceBifurcationsGrowingFromHeart(b, numSwapGrids);
	drawBifurcationTree("outputs/" + uniquer + "gradualGrowthWithSwapsTreeBuilt.png", b, width, height, border);
}

unsigned long myFact(unsigned long x){
	if(x < 2)
		return 1;
	return x*myFact(x - 1);
}

unsigned long myCeilSlash(double x){
	if(x < 0.0){
		cout << "\nError myCeilSlash(): x = " << x << " < 0.0" << endl;
		return 0;
	}
	if(x == ceil(x))
		return (unsigned long)abs(x) + 1;
	return (unsigned long)ceil(abs(x));
}

unsigned long myFloorSlash(double x){
	if(x < 0.0){
		cout << "\nError myFloorSlash(): x = " << x << " < 0.0" << endl;
		return 0;
	}
	if(x == floor(x) && x >= 1.0)
		return (unsigned long)abs(x) - 1;
	return (unsigned long)floor(abs(x));
}

double myLog(double base, double arg){
	return log(arg)/log(base);
}

unsigned long myBinomialCoefficientFACT(unsigned long n, unsigned long k){
	if(k > n){
		cout << "\nError myBinomialCoefficientFACT(): n = " << n << " but k = " << k << endl;
		return 0;
	}
	return (myFact(n)/myFact(k))/myFact(n - k);
}

unsigned long myBinomialCoefficientBAD(unsigned int n, unsigned int k){
	if(k > n){
		cout << "\nError myBinomialCoefficient(): n = " << n << " but k = " << k << endl;
		return 0;
	}
	vector<unsigned int> x;
	for(unsigned int i(k + 1); i <= n; i++)
		x.push_back(i);
	for(unsigned int i(k); i > 1; i--){
		for(unsigned int j(0); j < x.size(); j++){
			if(x[j]%i == 0){
				x[j] /= i;
				j = x.size();
			}
		}
	}
	unsigned long prod(1);
	for(unsigned int j(0); j < x.size(); j++)
		prod *= x[j];
	return prod;
}

// FROM https://gist.github.com/jeetsukumaran, specifically https://gist.github.com/jeetsukumaran/5392166
/**
* Calculates the binomial coefficient, $\choose{n, k}$, i.e., the number of
* distinct sets of $k$ elements that can be sampled with replacement from a
* population of $n$ elements.
*
* @tparam T
* Numeric type. Defaults to unsigned long.
* @param n
* Population size.
* @param k
* Number of elements to sample without replacement.
*
* @return
* The binomial coefficient, $\choose{n, k}$.
*
* @note
* Modified from: http://etceterology.com/fast-binomial-coefficients
*/
//template <class T>// T = unsigned long
unsigned long binomial_coefficient(unsigned long n, unsigned long k) {
	unsigned long i;
	unsigned long b;
	if(0 == k || n == k) {
		return 1;
	}
	if(k > n) {
		return 0;
	}
	if(k > (n - k)) {
		k = n - k;
	}
	if(1 == k) {
		return n;
	}
	b = 1;
	for(i = 1; i <= k; ++i) {
		b *= (n - (k - i));
		//doh impossible! if(b < 0) return -1; /* Overflow */
		b /= i;
	}
	return b;
}

unsigned int N_cap(double r, double r_cap, double beta, double lambda_R){
	unsigned int imin((unsigned int)ceil(myLog(beta*lambda_R, r_cap/r))),
	//unsigned int imin(myCeilSlash(myLog(beta*lambda_R, r_cap/r))),
		//imax(myCeilSlash(myLog(beta, r_cap/r)));
		imax((unsigned int)ceil(myLog(beta, r_cap/r)));
	//cout << "\n N_cap(): imin = " << imin << ", imax = " << imax << endl;
	unsigned int sum(0);
	for(unsigned int i(imin); i < imax + 1; i++){
		unsigned int jmin((unsigned int)ceil(myLog(lambda_R, pow(beta, -(int)i)*r_cap/r))),
			//jmax(floor(1 + myLog(lambda_R, pow(beta, 1 - (int)i)*r_cap/r)));
			jmax(myFloorSlash(1 + myLog(lambda_R, pow(beta, 1 - (int)i)*r_cap/r)));
		//if(jmin < 0) // <- won't happen, but I'll leave it in case jmin crosses the road (why?  to get to . . .)
		//	jmin = 0;
		if(jmin > jmax){
			//cout << "\n N_cap(): NOTE THAT jmin = " << jmin << " > " << jmax << " = jmax; jmin -> 0" << endl;
			jmin = 0;
		}
		if(jmax > i)
			jmax = i;
		//cout << "\n\t N_cap(): jmin = " << jmin << ", jmax = " << jmax << endl;
		for(unsigned int j(jmin); j < jmax + 1; j++){
			if(i > 0){
				if(j > 0 && j <= i)
					sum += binomial_coefficient(i - 1, j - 1);
				if(r_cap/r < pow(beta, double(i - 1))*pow(lambda_R, double(j)) && j < i)
					sum += binomial_coefficient(i - 1, j);
			}
		}
	}
	return sum;
}

unsigned int countTips(double r, double r_cap, double beta, double lambda_R){
	if(r <= r_cap)
		return 1;
	return countTips(r*beta, r_cap, beta, lambda_R) + countTips(r*beta*lambda_R, r_cap, beta, lambda_R);
}

void asymTheory(){
	double r(1.0), r_capSmall(1.0e-1), r_capLarge(5*r_capSmall), r_capInc(0.5*r_capSmall),
		betaSmall(0.3), betaLarge(0.9), betaInc(0.2),
		lambda_RSmall(0.3), lambda_RLarge(0.9), lambda_RInc(0.2);
	cout << "\n r = " << r
		<< "\n r_cap in [" << r_capSmall << ", " << r_capLarge << "] by " << r_capInc
		<< "\n beta in [" << betaSmall << ", " << betaLarge << "] by " << betaInc
		<< "\n lambda_R in [" << lambda_RSmall << ", " << lambda_RLarge << "] by " << lambda_RInc
		<< endl;
	cout.precision(4);
	cout << "\n\n r\t r_cap\t beta\t lambda_R\t N_cap(Theory)\t N_cap(Count)" << endl;
	for(double r_cap(r_capLarge); r_cap > r_capSmall - r_capInc/2.0; r_cap -= r_capInc){
		for(double beta(betaSmall); beta < betaLarge + betaInc/2.0; beta += betaInc){
			for(double lambda_R(lambda_RSmall); lambda_R < lambda_RLarge + lambda_RInc/2.0; lambda_R += lambda_RInc){
				unsigned int N_capTheory(N_cap(r, r_cap, beta, lambda_R)), N_capCount(countTips(r, r_cap, beta, lambda_R));
				if(N_capTheory != N_capCount)
					cout << "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" << endl;
				cout << r << "\t" << r_cap << "\t" << beta << "\t" << lambda_R << "\t"
					<< N_capTheory << "\t" << N_capCount << endl;
				if(N_capTheory != N_capCount)
					cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
			}
		}
	}
}

vector<Bifurcation2D> initializeEvenRandomSpacing(int maxConsecutiveFailedInsertions, double maxX, double maxY, double nearestRadius, double heartX, double heartY){
	int consecutiveFailedInsertions(0);
	vector<Bifurcation2D> b;
	b.push_back(Bifurcation2D());
	b.back().position[0] = heartX;
	b.back().position[1] = heartY;

	while(consecutiveFailedInsertions < maxConsecutiveFailedInsertions){
		double x(r2kiss()*maxX), y(r2kiss()*maxY);
		bool validLocation(true);
		for(unsigned int i(0); i < b.size(); i++){
			if(separation2D(x, y, b[i].position[0], b[i].position[1]) < nearestRadius){
				validLocation = false;
				i = b.size();
			}
		}
		if(validLocation){
			b.push_back(Bifurcation2D());
			b.back().position[0] = x;
			b.back().position[1] = y;
			consecutiveFailedInsertions = 0;
		}else
			consecutiveFailedInsertions++;
	}

	return b;
}

vector<Bifurcation3D> initializeEvenRandomSpacing(int maxConsecutiveFailedInsertions, double maxX, double maxY, double maxZ, double nearestRadius, double heartX, double heartY, double heartZ){
	int consecutiveFailedInsertions(0);
	vector<Bifurcation3D> b;
	b.push_back(Bifurcation3D());
	b.back().position[0] = heartX;
	b.back().position[1] = heartY;
	b.back().position[2] = heartZ;

	while(consecutiveFailedInsertions < maxConsecutiveFailedInsertions){
		double x(r2kiss()*maxX), y(r2kiss()*maxY), z(rkiss()*maxZ);
		bool validLocation(true);
		for(unsigned int i(0); i < b.size(); i++){
			if(separation3D(x, y, z, b[i].position[0], b[i].position[1], b[i].position[2]) < nearestRadius){
				validLocation = false;
				i = b.size();
			}
		}
		if(validLocation){
			b.push_back(Bifurcation3D());
			b.back().position[0] = x;
			b.back().position[1] = y;
			b.back().position[2] = z;
			consecutiveFailedInsertions = 0;
		}else
			consecutiveFailedInsertions++;
	}
	return b;
}

vector<Bifurcation3D> initializeEvenRandomSpacingSphere(int maxConsecutiveFailedInsertions, double sphereRadius, double nearestRadius){
	int consecutiveFailedInsertions(0);
	vector<Bifurcation3D> b;
	b.push_back(Bifurcation3D());
	b.back().position[0] = 0.0;
	b.back().position[1] = 0.0;
	b.back().position[2] = 0.0;

	while(consecutiveFailedInsertions < maxConsecutiveFailedInsertions){
		// more efficient to use spherical coordinates, but cartesian is fast enough
		double x(2.0*(r2kiss() - 0.5)*sphereRadius), y(2.0*(r2kiss() - 0.5)*sphereRadius), z(2.0*(r2kiss() - 0.5)*sphereRadius);
		while(separation3D(x, y, z, 0.0, 0.0, 0.0) > sphereRadius){
			x = 2.0*(r2kiss() - 0.5)*sphereRadius;
			y = 2.0*(r2kiss() - 0.5)*sphereRadius;
			z = 2.0*(r2kiss() - 0.5)*sphereRadius;
		}
		bool validLocation(true);
		for(unsigned int i(0); i < b.size(); i++){
			if(separation3D(x, y, z, b[i].position[0], b[i].position[1], b[i].position[2]) < nearestRadius){
				validLocation = false;
				i = b.size();
			}
		}
		if(validLocation){
			b.push_back(Bifurcation3D());
			b.back().position[0] = x;
			b.back().position[1] = y;
			b.back().position[2] = z;
			consecutiveFailedInsertions = 0;
		}else
			consecutiveFailedInsertions++;
	}
	return b;
}

void evenRandomSpacingTest(){

	kisset(24611, 19046, 7706);
	int maxFailedInsertions(10000);
	double sizeX(10.0), sizeY(10.0), nearestRadius(1.0);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	string uniquer("evenRandomSpacing" + makeString(sizeX));
	int width((int)floor(sizeX)*10), height((int)floor(sizeY)*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300;
	border = 10;

	cout << "\n evenRandomSpacingTest()\n uniquer = " << uniquer
		<< "\n sizeX = " << sizeX << "\n sizeY = " << sizeY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< endl;

	vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty));
	drawBifurcationTree("outputs/" + uniquer + "_initial.png", b, width, height, border);
	introduceBifurcationsFromSources(b);
	drawBifurcationTree("outputs/" + uniquer + "_built.png", b, width, height, border);
}

// box goes to (maxX, maxY), bend inner at (minX, minY)
vector<Bifurcation2D> initializeEvenRandomSpacingSingleBend(int maxConsecutiveFailedInsertions, double minX, double minY, double maxX, double maxY, double nearestRadius, double heartX, double heartY){
	int consecutiveFailedInsertions(0);
	vector<Bifurcation2D> b;
	b.push_back(Bifurcation2D());
	b.back().position[0] = heartX;
	b.back().position[1] = heartY;

	//maxConsecutiveFailedInsertions *= (int)ceil(maxX*maxY/((maxX*maxY) - (minX*minY)));
	while(consecutiveFailedInsertions < maxConsecutiveFailedInsertions){
		double x(r2kiss()*maxX), y(r2kiss()*maxY);
		while(x < minX && y < minY){
			x = r2kiss()*maxX;
			y = r2kiss()*maxY;
		}
		bool validLocation(true);
		for(unsigned int i(0); i < b.size(); i++){
			if(separation2D(x, y, b[i].position[0], b[i].position[1]) < nearestRadius){
				validLocation = false;
				i = b.size();
			}
		}
		if(validLocation){
			b.push_back(Bifurcation2D());
			b.back().position[0] = x;
			b.back().position[1] = y;
			consecutiveFailedInsertions = 0;
		}else
			consecutiveFailedInsertions++;
	}

	return b;
}

vector<Bifurcation2D> initializeEvenRandomSpacingCircle(int maxConsecutiveFailedInsertions, double side, double nearestRadius){
	int consecutiveFailedInsertions(0);
	vector<Bifurcation2D> b;
	b.push_back(Bifurcation2D());
	b.back().position[0] = side/2.0;
	b.back().position[1] = side/2.0;
	
	while(consecutiveFailedInsertions < maxConsecutiveFailedInsertions){
		double x(r2kiss()*side), y(r2kiss()*side);
		while(separation2D(side/2.0, side/2.0, x, y) > side/2.0){
			x = r2kiss()*side;
			y = r2kiss()*side;
		}
		bool validLocation(true);
		for(unsigned int i(0); i < b.size(); i++){
			if(separation2D(x, y, b[i].position[0], b[i].position[1]) < nearestRadius){
				validLocation = false;
				i = b.size();
			}
		}
		if(validLocation){
			b.push_back(Bifurcation2D());
			b.back().position[0] = x;
			b.back().position[1] = y;
			consecutiveFailedInsertions = 0;
		}else
			consecutiveFailedInsertions++;
	}

	return b;
}

vector<BranchPoint2D> consolidateDegenerates(const vector<Bifurcation2D> &b, double thresh);

void initializeEvenRandomSpacingSingleBendTest(){
	int numNetworks(1), maxConsecutiveFailedInsertions(10000), width(1000), height(width), border(5);
	double minX(10.0), minY(minX), maxX(2.0*minX), maxY(2.0*minY), nearestRadius(0.5/sqrt(acos(-1.0)));
	double consolidateThresh(0.1), settleThresh(0.001), heartx(nearestRadius), hearty((minY + maxY)/2.0);
	for(int n(0); n < numNetworks; n++){
		vector<Bifurcation2D> b(initializeEvenRandomSpacingSingleBend(maxConsecutiveFailedInsertions, minX, minY, maxX, maxY, nearestRadius, heartx, hearty));
		drawBentBifurcationTree("outputs/bendTest_" + makeString(n) + ".png", b, minX, minY, width, height, border);
		introduceBifurcationsFromSources(b);
		settleBifurcations(b, false, settleThresh);
		drawBentBifurcationTree("outputs/bendTest_" + makeString(n) + "_built.png", b, minX, minY, width, height, border);
		vector<BranchPoint2D> br(consolidateDegenerates(b, consolidateThresh));
		drawBentTree("outputs/bendTest_" + makeString(n) + "_consol.png", br, minX, minY, width, height, border);
	}
}

void copyBifurcation2DTreeFromTo(vector<Bifurcation2D> &copyFrom, vector<Bifurcation2D> &copyTo){
	if(copyFrom.size() != copyTo.size())
		cout << "\n Error copyBifurcation2DTreeFromTo(): implemented such that copyFrom.size() = " << copyFrom.size() << " is the same as copyTo.size() = " << copyTo.size() << endl;
	else{
		for(unsigned int i(0); i < copyFrom.size(); i++)
			copyTo[i] = copyFrom[i];
	}
}

//b[0] is the heart
vector<Bifurcation2D> simulatedAnnealing(vector<Bifurcation2D> b, int annealIts, double temperature, int &totalSwaps,
		double totalLengthCoefficient = 1.0, double avePathLengthCoefficient = 1.0, double maxPathLengthCoefficient = 1.0){

	double currentLength(totalLength(b)), currentAveDist(-1.0), currentMaxDist(-1.0);
	maxAndAverageDistanceFromHeart(currentMaxDist, currentAveDist, b);
	//double initLength(currentLength), initAveDist(currentAveDist), initMaxDist(currentMaxDist);

	//vector<int> swappedAtLeastOnce;

	//int *numSwaps = new int[b.size()];
	//for(unsigned int i(0); i < b.size(); i++)
	//	numSwaps[i] = 0;

	//cout << endl;
	vector<int> candidateIndices;
	for(unsigned int i(1); i < b.size(); i++){
		if(b[i].childIndex.size() > 0)
			candidateIndices.push_back(i);
	}
	totalSwaps = 0;
	for(int annealIt(0); annealIt < annealIts; annealIt++){
		vector<int> leftToVisit(candidateIndices);
		while(leftToVisit.size() > 0){
			int randomIndexIndex(kiss()%leftToVisit.size());
			int i(leftToVisit[randomIndexIndex]);
			leftToVisit.erase(leftToVisit.begin() + randomIndexIndex);
			int pIndex(b[i].parentIndex);
			if(pIndex > -1 && b[i].childIndex.size() > 0){
				int rIndex(b[b[i].parentIndex].parentIndex);
				if(rIndex > -1){
					unsigned int sEnum(kiss()%b[rIndex].childIndex.size());
					
					if(b[rIndex].childIndex[sEnum] == pIndex)
						continue;
					vector<Bifurcation2D> bTemp(b);
					int sIndex(b[rIndex].childIndex[sEnum]);
					Bifurcation2D *r(&(bTemp[rIndex])), *p(&(bTemp[pIndex])), *t(&(bTemp[i])), *s(&(bTemp[sIndex]));
					// t_par(p->r)
					t->parentIndex = rIndex;
					//r_child(-s, +t)
					//s_par(r->p)
					//p_child(-t, +s)
					replaceChildWith(r->childIndex, sIndex, i);
					s->parentIndex = pIndex;
					replaceChildWith(p->childIndex, i, sIndex);
					settleBifurcations(bTemp);
					
					//evaluate

					double tempLength(totalLength(bTemp)), tempAveDist(-1.0), tempMaxDist(-1.0);
					maxAndAverageDistanceFromHeart(tempMaxDist, tempAveDist, bTemp);
					double costChange(totalLengthCoefficient*(tempLength - currentLength)
							+ avePathLengthCoefficient*(tempAveDist - currentAveDist)
							+ maxPathLengthCoefficient*(tempMaxDist - currentMaxDist));
					bool acceptSwap(false);
					if(temperature <= 0.0)
						acceptSwap = costChange < 0.0;
					else
						acceptSwap = rkiss() < exp(-costChange/temperature);
					if(acceptSwap){
						//cout << "\n\t\t simulatedAnnealing(): accepted swap with costChange = " << costChange
						//	<< " at temperature = " << temperature << endl;
						copyBifurcation2DTreeFromTo(bTemp, b);
						currentLength = tempLength;
						currentAveDist = tempAveDist;
						currentMaxDist = tempMaxDist;
						totalSwaps++;
					}
				}
			}
		}
	}
	return b;
}

// assumes heart is bDegen[0]; will set heart as b[0]
vector<BranchPoint2D> consolidateDegenerates(const vector<Bifurcation2D> &bDegen, double thresh){
	vector<vector<int> > boundaries, gobs;
	findDegenerateGobs(bDegen, boundaries, gobs, thresh);
	int *consolidatedLabels = new int[bDegen.size()];
	for(unsigned int i(0); i < bDegen.size(); i++)
		consolidatedLabels[i] = -1;
	int unclaimedIndex(0);
	//cout << "\n\n bDegen\t b" << endl;
	for(unsigned int i(0); i < bDegen.size(); i++){
		if(consolidatedLabels[i] < 0){
			consolidatedLabels[i] = unclaimedIndex;
			for(unsigned int j(0); j < gobs.size(); j++){
				if(containedIn(i, gobs[j])){
					for(unsigned int k(0); k < gobs[j].size(); k++)
						consolidatedLabels[gobs[j][k]] = unclaimedIndex;
				}
			}
			unclaimedIndex++;
		}
		//cout << i << "\t" << consolidatedLabels[i] << endl;
	}
	vector<BranchPoint2D> b(unclaimedIndex, BranchPoint2D());
	for(unsigned int i(0); i < bDegen.size(); i++){
		b[consolidatedLabels[i]].position[0] = bDegen[i].position[0];
		b[consolidatedLabels[i]].position[1] = bDegen[i].position[1];
		for(unsigned int c(0); c < bDegen[i].childIndex.size(); c++){
			int cIndex(consolidatedLabels[bDegen[i].childIndex[c]]);
			if(consolidatedLabels[i] != cIndex){
				//cout << "\nAdding " << cIndex << " (formerly " << bDegen[i].childIndex[c] << ") to " << consolidatedLabels[i] << " (formerly " << i << ")" << endl;
				b[cIndex].parentIndex = consolidatedLabels[i];
				b[consolidatedLabels[i]].childIndex.push_back(cIndex);
			}
		}
	}
	//cout << "\n without degenerates:" << endl;
	//for(unsigned int i(0); i < b.size(); i++)
	//	cout << i << "\t" << b[i].tostring() << endl;
	delete[] consolidatedLabels;
	return b;
}

// assumes heart is bDegen[0]; will set heart as b[0]
vector<BranchPoint3D> consolidateDegenerates(const vector<Bifurcation3D> &bDegen, double thresh = 1.0e-6){
	vector<vector<int> > boundaries, gobs;
	findDegenerateGobs(bDegen, boundaries, gobs, thresh);
	int *consolidatedLabels = new int[bDegen.size()];
	for(unsigned int i(0); i < bDegen.size(); i++)
		consolidatedLabels[i] = -1;
	int unclaimedIndex(0);
	//cout << "\n\n bDegen\t b" << endl;
	for(unsigned int i(0); i < bDegen.size(); i++){
		if(consolidatedLabels[i] < 0){
			consolidatedLabels[i] = unclaimedIndex;
			for(unsigned int j(0); j < gobs.size(); j++){
				if(containedIn(i, gobs[j])){
					for(unsigned int k(0); k < gobs[j].size(); k++)
						consolidatedLabels[gobs[j][k]] = unclaimedIndex;
				}
			}
			unclaimedIndex++;
		}
		//cout << i << "\t" << consolidatedLabels[i] << endl;
	}
	vector<BranchPoint3D> b(unclaimedIndex, BranchPoint3D());
	for(unsigned int i(0); i < bDegen.size(); i++){
		b[consolidatedLabels[i]].position[0] = bDegen[i].position[0];
		b[consolidatedLabels[i]].position[1] = bDegen[i].position[1];
		b[consolidatedLabels[i]].position[2] = bDegen[i].position[2];
		for(unsigned int c(0); c < bDegen[i].childIndex.size(); c++){
			int cIndex(consolidatedLabels[bDegen[i].childIndex[c]]);
			if(consolidatedLabels[i] != cIndex){
				//cout << "\nAdding " << cIndex << " (formerly " << bDegen[i].childIndex[c] << ") to " << consolidatedLabels[i] << " (formerly " << i << ")" << endl;
				b[cIndex].parentIndex = consolidatedLabels[i];
				b[consolidatedLabels[i]].childIndex.push_back(cIndex);
			}
		}
	}
	//cout << "\n without degenerates:" << endl;
	//for(unsigned int i(0); i < b.size(); i++)
	//	cout << i << "\t" << b[i].tostring() << endl;
	delete[] consolidatedLabels;
	return b;
}

//b[0] is the heart
vector<Bifurcation2D> simulatedConsolidatedAnnealing(vector<Bifurcation2D> b, int annealIts, double temperature, int &totalSwaps,
	double settleThresh, double consolidateThresh,
	double totalLengthCoefficient = 1.0, double avePathLengthCoefficient = 1.0, double maxPathLengthCoefficient = 1.0){

	double currentLength(totalLength(b)), currentAveDist(-1.0), currentMaxDist(-1.0);
	maxAndAverageDistanceFromHeart(currentMaxDist, currentAveDist, b);
	//double initLength(currentLength), initAveDist(currentAveDist), initMaxDist(currentMaxDist);

	//vector<int> swappedAtLeastOnce;

	//int *numSwaps = new int[b.size()];
	//for(unsigned int i(0); i < b.size(); i++)
	//	numSwaps[i] = 0;

	//cout << endl;

	vector<Bifurcation2D> bBest(b);
	vector<BranchPoint2D> bcBest(consolidateDegenerates(b));
	double measureBest(tripletMeasure(bcBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));

	vector<int> candidateIndices;
	for(unsigned int i(1); i < b.size(); i++){
		if(b[i].childIndex.size() > 1)
			candidateIndices.push_back(i);
	}
	totalSwaps = 0;
	for(int annealIt(0); annealIt < annealIts; annealIt++){
		vector<int> leftToVisit(candidateIndices);
		while(leftToVisit.size() > 0){
			int randomIndexIndex(kiss()%leftToVisit.size());
			int i(leftToVisit[randomIndexIndex]);
			leftToVisit.erase(leftToVisit.begin() + randomIndexIndex);
			int pIndex(b[i].parentIndex);
			if(pIndex > -1 && b[i].childIndex.size() > 0){
				int rIndex(b[b[i].parentIndex].parentIndex);
				if(rIndex > -1){
					unsigned int sEnum(kiss()%b[rIndex].childIndex.size());

					if(b[rIndex].childIndex[sEnum] == pIndex)
						continue;
					vector<Bifurcation2D> bTemp(b);
					int sIndex(b[rIndex].childIndex[sEnum]);
					Bifurcation2D *r(&(bTemp[rIndex])), *p(&(bTemp[pIndex])), *t(&(bTemp[i])), *s(&(bTemp[sIndex]));
					// t_par(p->r)
					t->parentIndex = rIndex;
					//r_child(-s, +t)
					//s_par(r->p)
					//p_child(-t, +s)
					replaceChildWith(r->childIndex, sIndex, i);
					s->parentIndex = pIndex;
					replaceChildWith(p->childIndex, i, sIndex);
					settleBifurcations(bTemp, false, settleThresh);

					//evaluate

					vector<BranchPoint2D> bTempc(consolidateDegenerates(bTemp, consolidateThresh));
					double tempLength(totalLength(bTempc)), tempAveDist(-1.0), tempMaxDist(-1.0);
					maxAndAverageDistanceFromHeart(tempMaxDist, tempAveDist, bTempc);
					double costChange(totalLengthCoefficient*(tempLength - currentLength)
						+ avePathLengthCoefficient*(tempAveDist - currentAveDist)
						+ maxPathLengthCoefficient*(tempMaxDist - currentMaxDist));
					bool acceptSwap(false);
					if(temperature <= 0.0)
						acceptSwap = costChange < 0.0;
					else
						acceptSwap = rkiss() < exp(-costChange/temperature);
					if(acceptSwap){
						double tempMeasure(totalLengthCoefficient*tempLength + avePathLengthCoefficient*tempAveDist + maxPathLengthCoefficient*tempMaxDist < measureBest);
						if(measureBest > tempMeasure){
							measureBest = tempMeasure;
							bBest = b;
						}
						//cout << "\n\t\t simulatedAnnealing(): accepted swap with costChange = " << costChange
						//	<< " at temperature = " << temperature << endl;
						copyBifurcation2DTreeFromTo(bTemp, b);
						currentLength = tempLength;
						currentAveDist = tempAveDist;
						currentMaxDist = tempMaxDist;
						totalSwaps++;
					}
				}
			}
		}
	}
	return bBest;
}

vector<Bifurcation2D> simulatedCooling(vector<Bifurcation2D> b, double initialTemperature, int numTemperatureIncrements,
		int numTopKeep,
		int numIndependentAnnealingsPerTemperature, int annealingIterations,
		double totalLengthCoefficient = 1.0, double avePathLengthCoefficient = 1.0, double maxPathLengthCoefficient = 1.0,
		bool printReport = false){
	if(initialTemperature <= 0.0)
		numTemperatureIncrements = 1;
	else if(numTemperatureIncrements < 1)
		numTemperatureIncrements = 2;
	//vector<Bifurcation2D> bestB(b);
	//double bestLength(totalLength(bestB)), bestAveDist(-1.0), bestMaxDist(-1.0);
	//maxAndAverageDistanceFromHeart(bestMaxDist, bestAveDist, bestB);
	vector<Bifurcation2D>* bestBs = new vector<Bifurcation2D>[numTopKeep];
	double *bestLengths = new double[numTopKeep];
	double *bestAveDists = new double[numTopKeep];
	double *bestMaxDists = new double[numTopKeep];
	double *bestMeasures = new double[numTopKeep];
	bestBs[0] = b;
	bestLengths[0] = totalLength(bestBs[0]);
	bestAveDists[0] = bestMaxDists[0] = -1.0;
	maxAndAverageDistanceFromHeart(bestMaxDists[0], bestAveDists[0], bestBs[0]);
	bestMeasures[0] = totalLengthCoefficient*bestLengths[0] + avePathLengthCoefficient*bestAveDists[0] + maxPathLengthCoefficient*bestMaxDists[0];
	for(int k(1); k < numTopKeep; k++){
		bestBs[k] = bestBs[0];
		bestLengths[k] = bestLengths[0];
		bestAveDists[k] = bestAveDists[0];
		bestMaxDists[k] = bestMaxDists[0];
		bestMeasures[k] = bestMeasures[0];
	}
	
	for(int i(0); i < numTemperatureIncrements; i++){
		double temperature(initialTemperature - i*(initialTemperature/(numTemperatureIncrements - 1)));
		if(i == numTemperatureIncrements - 1)
			temperature = 0.0;
		vector<int> totSwap;

		vector<Bifurcation2D> *bestBsTemp = new vector<Bifurcation2D>[numTopKeep];
		double *bestTempLengths = new double[numTopKeep];
		double *bestTempAveDists = new double[numTopKeep];
		double *bestTempMaxDists = new double[numTopKeep];
		double *bestTempMeasures = new double[numTopKeep];
		for(int k(0); k < numTopKeep; k++){
			bestBsTemp[k] = bestBs[k];
			bestTempLengths[k] = bestLengths[k];
			bestTempAveDists[k] = bestAveDists[k];
			bestTempMaxDists[k] = bestMaxDists[k];
			bestTempMeasures[k] = bestMeasures[k];
		}
		for(int k(0); k < numTopKeep; k++){
			//cout << "\n simulatedCooling(): temperature = " << temperature << endl;
			for(int j(0); j < numIndependentAnnealingsPerTemperature; j++){
				int totalSwaps(-1);
				vector<Bifurcation2D> bAnneal = simulatedAnnealing(bestBs[k], annealingIterations, temperature, totalSwaps, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
				//cout << "\n\t simulatedCooling(): annealing " << j + 1 << " of " << numIndependentAnnealingsPerTemperature
				//	<< "; performed " << totalSwaps << " total swaps." << endl;
				totSwap.push_back(totalSwaps);
				double annealLength(totalLength(bAnneal)), annealAveDist(-1.0), annealMaxDist(-1.0);
				maxAndAverageDistanceFromHeart(annealMaxDist, annealAveDist, bAnneal);
				double annealMeasure(totalLengthCoefficient*annealLength + avePathLengthCoefficient*annealAveDist + maxPathLengthCoefficient*annealMaxDist);
				
				if(annealMeasure < bestTempMeasures[numTopKeep - 1]){
					bool isRedundant(false);
					for(int m(0); m < numTopKeep && !isRedundant; m++)
						isRedundant = !isRedundant && equivalentTrees(bAnneal, bestBsTemp[m]);
					if(isRedundant)
						continue;
					int betterThan(numTopKeep - 1);
					for(int m(numTopKeep - 2); m > -1; m--){
						if(annealMeasure < bestTempMeasures[m]){
							betterThan = m;
						}
					}
					for(int m(numTopKeep - 1); m > betterThan; m--){
						bestBsTemp[m] = bestBsTemp[m - 1];
						bestTempLengths[m] = bestTempLengths[m - 1];
						bestTempAveDists[m] = bestTempAveDists[m - 1];
						bestTempMaxDists[m] = bestTempMaxDists[m - 1];
						bestTempMeasures[m] = bestTempMeasures[m - 1];
					}
					bestBsTemp[betterThan] = bAnneal;
					bestTempLengths[betterThan] = annealLength;
					bestTempAveDists[betterThan] = annealAveDist;
					bestTempMaxDists[betterThan] = annealMaxDist;
					bestTempMeasures[betterThan] = annealMeasure;
				}
			}
		}
		for(int k(0); k < numTopKeep; k++){
			bestBs[k] = bestBsTemp[k];
			bestLengths[k] = bestTempLengths[k];
			bestAveDists[k] = bestTempAveDists[k];
			bestMaxDists[k] = bestTempMaxDists[k];
			bestMeasures[k] = bestTempMeasures[k];
		}
		delete[] bestBsTemp;
		delete[] bestTempLengths;
		delete[] bestTempAveDists;
		delete[] bestTempMaxDists;
		if(printReport)
			cout << "\n\t simulatedCooling(): [temperature = " << temperature << "] <totSwaps> = " << mean(totSwap) << " +/- " << stDev(totSwap) << endl;
	}
	vector<Bifurcation2D> bestB(bestBs[0]);
	delete[] bestBs;
	delete[] bestLengths;
	delete[] bestAveDists;
	delete[] bestMaxDists;
	return bestB;
}

// b is a guess, fully built
vector<BranchPoint2D> simulatedConsolidatedCooling(vector<Bifurcation2D> b, double initialTemperature, int numTemperatureIncrements,
	int numTopKeep, int numIndependentAnnealingsPerTemperature, int annealingIterations,
	double settleThresh, double consolidateThresh,
	double totalLengthCoefficient = 1.0, double avePathLengthCoefficient = 1.0, double maxPathLengthCoefficient = 1.0,
	bool printReport = false){
	if(initialTemperature <= 0.0)
		numTemperatureIncrements = 1;
	else if(numTemperatureIncrements < 1)
		numTemperatureIncrements = 2;
	//vector<Bifurcation2D> bestB(b);
	//double bestLength(totalLength(bestB)), bestAveDist(-1.0), bestMaxDist(-1.0);
	//maxAndAverageDistanceFromHeart(bestMaxDist, bestAveDist, bestB);
	vector<Bifurcation2D>* bestBs = new vector<Bifurcation2D>[numTopKeep];
	double *bestLengths = new double[numTopKeep];
	double *bestAveDists = new double[numTopKeep];
	double *bestMaxDists = new double[numTopKeep];
	double *bestMeasures = new double[numTopKeep];
	bestBs[0] = b;
	bestLengths[0] = totalLength(bestBs[0]);
	bestAveDists[0] = bestMaxDists[0] = -1.0;
	maxAndAverageDistanceFromHeart(bestMaxDists[0], bestAveDists[0], bestBs[0]);
	bestMeasures[0] = totalLengthCoefficient*bestLengths[0] + avePathLengthCoefficient*bestAveDists[0] + maxPathLengthCoefficient*bestMaxDists[0];
	for(int k(1); k < numTopKeep; k++){
		bestBs[k] = bestBs[0];
		bestLengths[k] = bestLengths[0];
		bestAveDists[k] = bestAveDists[0];
		bestMaxDists[k] = bestMaxDists[0];
		bestMeasures[k] = bestMeasures[0];
	}

	for(int i(0); i < numTemperatureIncrements; i++){
		double temperature(initialTemperature - i*(initialTemperature/(numTemperatureIncrements - 1)));
		if(i == numTemperatureIncrements - 1)
			temperature = 0.0;
		vector<int> totSwap;

		vector<Bifurcation2D> *bestBsTemp = new vector<Bifurcation2D>[numTopKeep];
		double *bestTempLengths = new double[numTopKeep];
		double *bestTempAveDists = new double[numTopKeep];
		double *bestTempMaxDists = new double[numTopKeep];
		double *bestTempMeasures = new double[numTopKeep];
		for(int k(0); k < numTopKeep; k++){
			bestBsTemp[k] = bestBs[k];
			bestTempLengths[k] = bestLengths[k];
			bestTempAveDists[k] = bestAveDists[k];
			bestTempMaxDists[k] = bestMaxDists[k];
			bestTempMeasures[k] = bestMeasures[k];
		}
		for(int k(0); k < numTopKeep; k++){
			for(int j(0); j < numIndependentAnnealingsPerTemperature; j++){
				int totalSwaps(-1);
				vector<Bifurcation2D> bAnneal = simulatedConsolidatedAnnealing(bestBs[k], annealingIterations, temperature, totalSwaps, settleThresh, consolidateThresh, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
				//cout << "\n\t simulatedCooling(): annealing " << j + 1 << " of " << numIndependentAnnealingsPerTemperature
				//	<< "; performed " << totalSwaps << " total swaps." << endl;
				totSwap.push_back(totalSwaps);
				double annealLength(totalLength(bAnneal)), annealAveDist(-1.0), annealMaxDist(-1.0);
				maxAndAverageDistanceFromHeart(annealMaxDist, annealAveDist, bAnneal);
				double annealMeasure(totalLengthCoefficient*annealLength + avePathLengthCoefficient*annealAveDist + maxPathLengthCoefficient*annealMaxDist);

				if(annealMeasure < bestTempMeasures[numTopKeep - 1]){
					bool isRedundant(false);
					for(int m(0); m < numTopKeep && !isRedundant; m++)
						isRedundant = !isRedundant && equivalentTrees(bAnneal, bestBsTemp[m]);
					if(isRedundant)
						continue;
					int betterThan(numTopKeep - 1);
					for(int m(numTopKeep - 2); m > -1; m--){
						if(annealMeasure < bestTempMeasures[m]){
							betterThan = m;
						}
					}
					for(int m(numTopKeep - 1); m > betterThan; m--){
						bestBsTemp[m] = bestBsTemp[m - 1];
						bestTempLengths[m] = bestTempLengths[m - 1];
						bestTempAveDists[m] = bestTempAveDists[m - 1];
						bestTempMaxDists[m] = bestTempMaxDists[m - 1];
						bestTempMeasures[m] = bestTempMeasures[m - 1];
					}
					bestBsTemp[betterThan] = bAnneal;
					bestTempLengths[betterThan] = annealLength;
					bestTempAveDists[betterThan] = annealAveDist;
					bestTempMaxDists[betterThan] = annealMaxDist;
					bestTempMeasures[betterThan] = annealMeasure;
				}
			}
		}
		for(int k(0); k < numTopKeep; k++){
			bestBs[k] = bestBsTemp[k];
			bestLengths[k] = bestTempLengths[k];
			bestAveDists[k] = bestTempAveDists[k];
			bestMaxDists[k] = bestTempMaxDists[k];
			bestMeasures[k] = bestTempMeasures[k];
		}
		delete[] bestBsTemp;
		delete[] bestTempLengths;
		delete[] bestTempAveDists;
		delete[] bestTempMaxDists;
		if(printReport)
			cout << "\n\t simulatedCooling(): [temperature = " << temperature << "] <totSwaps> = " << mean(totSwap) << " +/- " << stDev(totSwap) << endl;
	}
	vector<BranchPoint2D> bestB(consolidateDegenerates(bestBs[0]));
	delete[] bestBs;
	delete[] bestLengths;
	delete[] bestAveDists;
	delete[] bestMaxDists;
	return bestB;
}

void annealingTest(){
	kisset(24611, 19046, 7706);
	int maxFailedInsertions(100000), numTempIncs(3), numTopKeep(15), numIndAnnealings(50), annealIts(50);
	double temp(0.5), sizeX(7.0), sizeY(7.0), nearestRadius(1.0);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	string uniquer("evenRandomSpacingForAnnealing_" + makeString(sizeX));
	int width((int)floor(sizeX)*10), height((int)floor(sizeY)*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300/2;
	border = 10/2;

	cout << "\n annealingTest()\n uniquer = " << uniquer
		<< "\n sizeX = " << sizeX << "\n sizeY = " << sizeY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< "\n temp = " << temp << "\n numTempIncs = " << numTempIncs
		<< "\n numIndAnnealings = " << numIndAnnealings << "\n annealIts = " << annealIts
		<< endl;

	vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty));
	introduceBifurcationsFromSources(b);
	drawBifurcationTree("outputs/" + uniquer + "_initial.png", b, width, height, border);

	int totalSwaps(-1);
	double initLength(totalLength(b)), initAveDist(-1.0), initMaxDist(-1.0);
	maxAndAverageDistanceFromHeart(initMaxDist, initAveDist, b);
	vector<Bifurcation2D> bAnnealed = simulatedAnnealing(b, annealIts, temp, totalSwaps);
	drawBifurcationTree("outputs/" + uniquer + "_annealed.png", bAnnealed, width, height, border);
	vector<Bifurcation2D> bCooled = simulatedCooling(b, temp, numTempIncs, numTopKeep, numIndAnnealings, annealIts);
	drawBifurcationTree("outputs/" + uniquer + "_cooled.png", bCooled, width, height, border);
	cout << "\n\t\t cooling 1 . . ." << endl;
	vector<Bifurcation2D> bCooled1 = simulatedCooling(b, temp, numTempIncs, numTopKeep, numIndAnnealings, annealIts);
	drawBifurcationTree("outputs/" + uniquer + "_cooled1.png", bCooled1, width, height, border);
	cout << "\n\t\t cooling 2 . . ." << endl;
	vector<Bifurcation2D> bCooled2 = simulatedCooling(b, temp, numTempIncs, numTopKeep, numIndAnnealings, annealIts);
	drawBifurcationTree("outputs/" + uniquer + "_cooled2.png", bCooled2, width, height, border);
	cout << "\n\t\t cooling 3 . . ." << endl;
	vector<Bifurcation2D> bCooled3 = simulatedCooling(b, temp, numTempIncs, numTopKeep, numIndAnnealings, annealIts);
	drawBifurcationTree("outputs/" + uniquer + "_cooled3.png", bCooled3, width, height, border);

	double checkLength(totalLength(b)), checkAveDist(-1.0), checkMaxDist(-1.0);
	maxAndAverageDistanceFromHeart(checkMaxDist, checkAveDist, b);
	double annealLength(totalLength(bAnnealed)), annealAveDist(-1.0), annealMaxDist(-1.0);
	maxAndAverageDistanceFromHeart(annealMaxDist, annealAveDist, bAnnealed);
	double coolLength(totalLength(bCooled)), coolAveDist(-1.0), coolMaxDist(-1.0);
	maxAndAverageDistanceFromHeart(coolMaxDist, coolAveDist, bCooled);
	double cool1Length(totalLength(bCooled1)), cool1AveDist(-1.0), cool1MaxDist(-1.0);
	maxAndAverageDistanceFromHeart(cool1MaxDist, cool1AveDist, bCooled1);
	double cool2Length(totalLength(bCooled2)), cool2AveDist(-1.0), cool2MaxDist(-1.0);
	maxAndAverageDistanceFromHeart(cool2MaxDist, cool2AveDist, bCooled2);
	double cool3Length(totalLength(bCooled3)), cool3AveDist(-1.0), cool3MaxDist(-1.0);
	maxAndAverageDistanceFromHeart(cool3MaxDist, cool3AveDist, bCooled3);
	cout << "\n\n totalSwaps = " << totalSwaps << "\n\n version\t totalLength\t aveDist\t maxDist\t\t evaluation" << endl;
	cout << " init\t" << initLength << "\t" << initAveDist << "\t" << initMaxDist
		<< "\t\t" << initLength + initAveDist + initMaxDist
		<< "\n check\t" << checkLength << "\t" << checkAveDist << "\t" << checkMaxDist
		<< "\t\t" << checkLength + checkAveDist + checkMaxDist
		<< "\n anneal\t" << annealLength << "\t" << annealAveDist << "\t" << annealMaxDist
		<< "\t\t" << annealLength + annealAveDist + annealMaxDist
		<< "\n cool\t" << coolLength << "\t" << coolAveDist << "\t" << coolMaxDist
		<< "\t\t" << coolLength + coolAveDist + coolMaxDist
		<< "\n cool1\t" << cool1Length << "\t" << cool1AveDist << "\t" << cool1MaxDist
		<< "\t\t" << cool1Length + cool1AveDist + cool1MaxDist
		<< "\n cool2\t" << cool2Length << "\t" << cool2AveDist << "\t" << cool2MaxDist
		<< "\t\t" << cool2Length + cool2AveDist + cool2MaxDist
		<< "\n cool3\t" << cool3Length << "\t" << cool3AveDist << "\t" << cool3MaxDist
		<< "\t\t" << cool3Length + cool3AveDist + cool3MaxDist
		<< endl;
}

// b[0] is the heart
// not quite perfect: double-counts symmetric cases
void exhaustiveHierarchyConfigurationSearch(vector<Bifurcation2D> &b, int numBest, vector<Bifurcation2D> *bestBs, double *bestBMeasures,
	bool report = true, double totalLengthCoefficient = 1.0, double avePathLengthCoefficient = 1.0, double maxPathLengthCoefficient = 1.0, bool topCall = true,
	double settleThresh = 1.0e-3){
	//cout << "\n exhaustiveHierarchyConfigurationSearch()" << endl;
	int maxIts(200000);
	vector<unsigned int> roots;
	for(unsigned int i(1); i < b.size(); i++){
		if(b[i].parentIndex < 0)
			roots.push_back(i);
	}
	if(roots.size() == 1){
		b[0].childIndex.push_back(roots[0]);
		b[roots[0]].parentIndex = 0;
		vector<Bifurcation2D> bTemp(b);
		for(int i(0); i < numBest; i++){
			if(equivalentTrees(bestBs[i], b)){ // no need to explore
				b[0].childIndex.clear();
				b[roots[0]].parentIndex = -1;
				return;
			}
		}
		settleBifurcations(bTemp, false, settleThresh, maxIts);
		double totLen(totalLength(bTemp)), avePath(-1.0), maxPath(-1.0);
		maxAndAverageDistanceFromHeart(maxPath, avePath, bTemp);
		double tempMeasure(totalLengthCoefficient*totLen + avePathLengthCoefficient*avePath + maxPathLengthCoefficient*maxPath);
		for(int i(0); i < numBest; i++){
			if(bestBMeasures[i] < 0.0){
				bestBs[i] = bTemp;
				bestBMeasures[i] = tempMeasure;
				i = numBest;
			}else if(tempMeasure < bestBMeasures[i]){
				for(int j(numBest - 1); j > i; j--){
					if(bestBMeasures[j - 1] >= 0.0){
						bestBMeasures[j] = bestBMeasures[j - 1];
						bestBs[j] = bestBs[j - 1];
					}
				}
				bestBMeasures[i] = tempMeasure;
				bestBs[i] = bTemp;
				i = numBest;
			}
		}
		b[0].childIndex.clear();
		b[roots[0]].parentIndex = -1;
	}
	int totalSteps(0), stepCount(0);
	for(unsigned int i(0); i < roots.size(); i++)
		totalSteps += roots.size() - i - 1;
	time_t local_start_t(time(NULL));
	for(unsigned int i(0); i < roots.size(); i++){
		for(unsigned int j(i + 1); j < roots.size(); j++){
			addParent(b, roots[i], roots[j]);
			exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, report, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, false, settleThresh);
			removeBifurcation(b, b.size() - 1);
			stepCount++;
			if(topCall && report){
				time_t time_now(time(NULL));
				cout << endl << "exhaustiveHierarchyConfigurationSearch(): completed " << stepCount << " of " << totalSteps
					<< "\n\t\t" << niceTime(difftime(time_now, local_start_t)) << " passed"
					<< "\n\t\t" << niceTime(difftime(time_now, local_start_t)*double(totalSteps - stepCount)/double(stepCount)) << " to go"
					<< endl;
			}
		}
	}
}

//int tempGlobal(0);

void exhaustiveConsolidatedHierarchyConfigurationSearch(vector<Bifurcation2D> &b, int numBest, vector<BranchPoint2D> *bestBs, double *bestBMeasures,
	bool report = true, double totalLengthCoefficient = 1.0, double avePathLengthCoefficient = 1.0, double maxPathLengthCoefficient = 1.0,
	double settleThresh = 1.0e-3, double consolidateThresh = 1.0e-3, bool topCall = true){
	//cout << "\n exhaustiveHierarchyConfigurationSearch()" << endl;
	if(topCall){
		for(int i(0); i < numBest; i++){
			bestBMeasures[i] = -1.0;
			bestBs[i] = vector<BranchPoint2D>();
		}
	}
	int maxIts(200000);
	vector<unsigned int> roots;
	for(unsigned int i(1); i < b.size(); i++){
		if(b[i].parentIndex < 0)
			roots.push_back(i);
	}
	if(roots.size() == 1){
		b[0].childIndex.push_back(roots[0]);
		b[roots[0]].parentIndex = 0;
		vector<Bifurcation2D> bTemp(b);
		//for(int i(0); i < numBest; i++){
		//	if(equivalentTrees(bestBs[i], b)){ // no need to explore
		//		b[0].childIndex.clear();
		//		b[roots[0]].parentIndex = -1;
		//		return;
		//	}
		//}
		//cout << "\n searching " << makeStringTreeOneLine(bTemp) << endl;
		settleBifurcations(bTemp, false, settleThresh, maxIts);
		vector<BranchPoint2D> br(consolidateDegenerates(bTemp, consolidateThresh));
		//cout << "\n biTree = " << makeStringTreeOneLine(bTemp)
		//	<< "\n\t brTree = " << makeStringTreeOneLine(br) << endl;
		double totLen(totalLength(br)), avePath(-1.0), maxPath(-1.0);
		maxAndAverageDistanceFromHeart(maxPath, avePath, br);
		double tempMeasure(totalLengthCoefficient*totLen + avePathLengthCoefficient*avePath + maxPathLengthCoefficient*maxPath);
		bool foundEquivalent(false);
		for(int i(0); i < numBest && !foundEquivalent; i++){
			if(equivalentTrees(bestBs[i], br)){
				foundEquivalent = true;
				if(bestBMeasures[i] > tempMeasure){
					bestBMeasures[i] = tempMeasure;
					bestBs[i] = br;
				}
				if(i > 0){
					if(bestBMeasures[i] < bestBMeasures[i - 1]){
						for(int j(i); j > 1 && bestBMeasures[j] < bestBMeasures[j - 1]; j--){
							double meas(bestBMeasures[j]);
							bestBMeasures[j] = bestBMeasures[j - 1];
							bestBMeasures[j - 1] = meas;
							vector<BranchPoint2D> v(bestBs[j]);
							bestBs[j] = bestBs[j - 1];
							bestBs[j - 1] = v;
						}
					}
				}
				/*if(i > 0){
					if(bestBMeasures[i - 1] > tempMeasure){
						cout << "\n exhaustiveConsolidatedHierarchyConfigurationSearch(): Note that config is the same as "
							<< i << " but config's measure is " << tempMeasure << " while config " << i - 1 << " has measure "
							<< bestBMeasures[i - 1] << endl;
					}
				}
				if(i < numBest - 1){
					if(bestBMeasures[i + 1] >= 0.0 && bestBMeasures[i + 1] < tempMeasure){
						cout << "\n exhaustiveConsolidatedHierarchyConfigurationSearch(): Note that config is the same as "
							<< i << " but config's measure is " << tempMeasure << " while config " << i + 1 << " has measure "
							<< bestBMeasures[i + 1] << endl;
					}
				}
				drawTree("outputs/exhCons_dup" + makeString(tempGlobal) + "_keep.png", bestBs[i], 300, 300, 12);
				drawTree("outputs/exhCons_dup" + makeString(tempGlobal) + "_ignore.png", br, 300, 300, 12);
				tempGlobal++;*/
			}
		}
		if(!foundEquivalent){
			for(int i(0); i < numBest; i++){
				if(bestBMeasures[i] < 0.0){
					bestBs[i] = br;
					bestBMeasures[i] = tempMeasure;
					i = numBest;
				}else if(tempMeasure < bestBMeasures[i]){
					for(int j(numBest - 1); j > i; j--){
						if(bestBMeasures[j - 1] >= 0.0){
							bestBMeasures[j] = bestBMeasures[j - 1];
							bestBs[j] = bestBs[j - 1];
						}
					}
					bestBMeasures[i] = tempMeasure;
					bestBs[i] = br;
					i = numBest;
				}
			}
		}
		b[0].childIndex.clear();
		b[roots[0]].parentIndex = -1;
	}
	int totalSteps(0), stepCount(0);
	for(unsigned int i(0); i < roots.size(); i++)
		totalSteps += roots.size() - i - 1;
	time_t local_start_t(time(NULL));
	for(unsigned int i(0); i < roots.size(); i++){
		for(unsigned int j(i + 1); j < roots.size(); j++){
			addParent(b, roots[i], roots[j]);
			exhaustiveConsolidatedHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, report, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, false);
			removeBifurcation(b, b.size() - 1);
			stepCount++;
			if(topCall && report){
				time_t time_now(time(NULL));
				cout << endl << "exhaustiveConsolidatedHierarchyConfigurationSearch(): completed " << stepCount << " of " << totalSteps
					<< "\n\t\t" << niceTime(difftime(time_now, local_start_t)) << " passed"
					<< "\n\t\t" << niceTime(difftime(time_now, local_start_t)*double(totalSteps - stepCount)/double(stepCount)) << " to go"
					<< endl;
			}
		}
	}
}

unsigned long myDoubleFact(int n){
	if(n < -1){
		cout << "\nError myDoubleFact(): expecting an argument no less than -1, but n = " << n << endl;
		return 0;
	}
	if(n == -1 || n == 0 || n == 1)
		return 1;
	return n*myDoubleFact(n - 2);
}

unsigned long numUniqueTrees(int n){
	return myDoubleFact(2*n - 3);
}

void exhaustiveSearchTest(){
	kisset(24611, 19046, 7706);
	int maxFailedInsertions(100000), numBest(100);
	double sizeX(2.0), sizeY(2.0), nearestRadius(1.0);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	string uniquer("exhaustiveSearchTest_" + makeString(sizeX));
	int width((int)floor(sizeX)*10), height((int)floor(sizeY)*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300/2;
	border = 10/2;

	cout << "\n exhaustiveSearchTest()\n uniquer = " << uniquer
		<< "\n sizeX = " << sizeX << "\n sizeY = " << sizeY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n numBest = " << numBest
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< endl;
	vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty));
	drawBifurcationTree("outputs/" + uniquer + "_initial.png", b, width, height, border);
	vector<Bifurcation2D> *bestBs = new vector<Bifurcation2D>[numBest];
	double *bestBMeasures = new double[numBest];
	for(int i(0); i < numBest; i++)
		bestBMeasures[i] = -1.0;
	exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures);
	cout << "\n\n bestBs\t measure\t normalized_measure" << endl;
	for(int i(0); i < numBest; i++){
		if(bestBMeasures[i] >= 0.0){
			cout << makeStringTreeOneLine(bestBs[i]) << "\t" << bestBMeasures[i]
			<< "\t" << bestBMeasures[i]/bestBMeasures[0] << endl;
			drawBifurcationTree("outputs/" + uniquer + "_ranked" + makeString(i) + ".png", bestBs[i], width, height, border);
		}
	}

	delete[] bestBs;
	delete[] bestBMeasures;
}

void exhaustiveTrends(){
	kisset(24611, 19046, 7706);
	int maxFailedInsertions(100000), numBest(2704), numNetworks(1000), maxFailedNetworks(10000), targetServiceVolumes(8);
	double nearestRadius(1.0);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	string uniquer("exhaustiveTrends" + makeString(sizeX) + "_2704_");
	int width((int)floor(sizeX)*10), height((int)floor(sizeY)*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300/2;
	border = 10/2;

	cout << "\n exhaustiveTrends()\n uniquer = " << uniquer
		<< "\n numNetworks = " << numNetworks << "\n maxFailedNetworks = " << maxFailedNetworks
		<< "\n targetServiceVolumes = " << targetServiceVolumes
		<< "\n sizeX = " << sizeX << "\n sizeY = " << sizeY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n numBest = " << numBest
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< endl << endl;

	int networkCount(0), failedNetworks(0);
	double **bestBMeasures = new double*[numNetworks];
	vector<Bifurcation2D> **bestBs = new vector<Bifurcation2D>*[numNetworks];
	cout << endl;
	while(networkCount < numNetworks && failedNetworks < maxFailedNetworks){
		vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty));
		if(b.size() == targetServiceVolumes + 1){
			cout << "{" << failedNetworks << "}";
			cout.flush();
			//drawBifurcationTree("outputs/" + uniquer + "_initial.png", b, width, height, border);
			bestBs[networkCount] = new vector<Bifurcation2D>[numBest];
			bestBMeasures[networkCount] = new double[numBest];
			for(int i(0); i < numBest; i++)
				bestBMeasures[networkCount][i] = -1.0;
			exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs[networkCount], bestBMeasures[networkCount], false);
			failedNetworks = 0;
			networkCount++;
		}else{
			failedNetworks++;
			//cout << "failed: b.size() = " << b.size() << " != " << targetServiceVolumes + 1 << endl;
		}
	}
	cout << "\n\n rank\t ave_measure\t +/-\t ave_norm_measure\t +/-" << endl;
	for(int i(0); i < numBest; i++){
		vector<double> measures, normMeasures;
		for(int n(0); n < networkCount; n++){
			if(bestBMeasures[n][i] >= 0.0){
				measures.push_back(bestBMeasures[n][i]);
				normMeasures.push_back(bestBMeasures[n][i]/bestBMeasures[n][0]);
			}
		}
		if(measures.size() > 0){
			cout << i + 1 << "\t" << mean(measures) << "\t" << stDev(measures)
				<< "\t" << mean(normMeasures) << "\t" << stDev(normMeasures) << endl;
		}else
			cout << i + 1 << "\t X\t X\t X\t X" << endl;
	}
	cout << "\n\n measures..." << endl;
	for(int i(0); i < numBest; i++){
		cout << i + 1;
		for(int n(0); n < networkCount; n++){
			if(bestBMeasures[n][i] >= 0.0)
				cout << "\t" << makeStringTreeOneLine(bestBs[n][i]) << "\t" << bestBMeasures[n][i];
			else
				cout << "\t X\t X";
		}
		cout << endl;
	}
	cout << "\n\n normalized_measures..." << endl;
	for(int i(0); i < numBest; i++){
		cout << i + 1;
		for(int n(0); n < networkCount; n++){
			if(bestBMeasures[n][i] >= 0.0)
				cout << "\t" << makeStringTreeOneLine(bestBs[n][i]) << "\t" << bestBMeasures[n][i]/bestBMeasures[n][0];
			else
				cout << "\tX \tX";
		}
		cout << endl;
	}
	for(int n(0); n < networkCount; n++){
		delete[] bestBs[n];
		delete[] bestBMeasures[n];
	}

	delete[] bestBs;
	delete[] bestBMeasures;
}

void exhaustiveTrendsMem(){
	kisset(24611, 19046, 7706);
	// 1 2 3  4  5   6    7      8      9       10
	// 1 1 3 15 105 945 10395 135135 2027025 34459425
	int maxFailedInsertions(100000), numBest(15), numNetworks(1), maxFailedNetworks(10000), targetServiceVolumes(4);
	double nearestRadius(1.0);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	double totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0);
	string uniquer("exhaustiveTrendsMem_" + makeString(targetServiceVolumes) + "_" + makeString(totalLengthCoefficient) + "TL_" + makeString(avePathLengthCoefficient) + "AP_" + makeString(maxPathLengthCoefficient) + "MP");
	int width(120), height(120), border(10);

	cout << "\n exhaustiveTrendsMem()\n uniquer = " << uniquer
		<< "\n numNetworks = " << numNetworks << "\n maxFailedNetworks = " << maxFailedNetworks
		<< "\n targetServiceVolumes = " << targetServiceVolumes
		<< "\n sizeX = " << sizeX << "\n sizeY = " << sizeY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n numBest = " << numBest
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< endl << endl;

	int networkCount(0), failedNetworks(0);
	double **bestBMeasures = new double*[numNetworks];
	cout << endl;
	while(networkCount < numNetworks && failedNetworks < maxFailedNetworks){
		vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty));
		if(b.size() == targetServiceVolumes + 1){
			cout << "{" << failedNetworks << "}";
			cout.flush();
			//drawBifurcationTree("outputs/" + uniquer + "_initial.png", b, width, height, border);
			vector<Bifurcation2D> *bestBs = new vector<Bifurcation2D>[numBest];
			bestBMeasures[networkCount] = new double[numBest];
			for(int i(0); i < numBest; i++)
				bestBMeasures[networkCount][i] = -1.0;
			exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures[networkCount], false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			if(networkCount == 0){
				cout << "\n\n configurations for first network:" << endl;
				for(int j(0); j < numBest; j++){
					drawBifurcationTree(uniquer + "_" + makeString(j) + ".png", bestBs[j], width, height, border);
					cout << j + 1 << "\t" << makeStringTreeOneLine(bestBs[j]) << endl;
				}
				cout << endl;
			}
			delete[] bestBs;
			failedNetworks = 0;
			networkCount++;
		} else{
			failedNetworks++;
			//cout << "failed: b.size() = " << b.size() << " != " << targetServiceVolumes + 1 << endl;
		}
	}
	cout << "\n\n rank\t ave_measure\t +/-\t ave_norm_measure\t +/-" << endl;
	for(int i(0); i < numBest; i++){
		vector<double> measures, normMeasures;
		for(int n(0); n < networkCount; n++){
			if(bestBMeasures[n][i] >= 0.0){
				measures.push_back(bestBMeasures[n][i]);
				normMeasures.push_back(bestBMeasures[n][i]/bestBMeasures[n][0]);
			}
		}
		if(measures.size() > 0){
			cout << i + 1 << "\t" << mean(measures) << "\t" << stDev(measures)
				<< "\t" << mean(normMeasures) << "\t" << stDev(normMeasures) << endl;
		} else
			cout << i + 1 << "\t X\t X\t X\t X" << endl;
	}
	cout << "\n\n rank\t measures..." << endl;
	for(int i(0); i < numBest; i++){
		vector<double> measures;
		cout << i + 1;
		for(int n(0); n < networkCount; n++){
			if(bestBMeasures[n][i] >= 0.0)
				cout << "\t" << bestBMeasures[n][i];
		}
		cout << endl;
	}
	cout << "\n\n rank\t norm_measures..." << endl;
	for(int i(0); i < numBest; i++){
		vector<double> measures;
		cout << i + 1;
		for(int n(0); n < networkCount; n++){
			if(bestBMeasures[n][i] >= 0.0)
				cout << "\t" << bestBMeasures[n][i]/bestBMeasures[n][0];
		}
		cout << endl;
	}
	for(int n(0); n < networkCount; n++)
		delete[] bestBMeasures[n];
	delete[] bestBMeasures;
}

void exhaustiveProblemSettle(){
	kisset(24611, 19046, 7706);
	int maxFailedInsertions(100000), numBest(1), numNetworks(52), maxFailedNetworks(10000), targetServiceVolumes(7);
	double nearestRadius(1.0);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	string uniquer("exhaustiveProblemSettle" + makeString(sizeX));
	int width((int)floor(sizeX)*10), height((int)floor(sizeY)*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300/2;
	border = 10/2;

	cout << "\n exhaustiveTrends()\n uniquer = " << uniquer
		<< "\n numNetworks = " << numNetworks << "\n maxFailedNetworks = " << maxFailedNetworks
		<< "\n targetServiceVolumes = " << targetServiceVolumes
		<< "\n sizeX = " << sizeX << "\n sizeY = " << sizeY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n numBest = " << numBest
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< endl << endl;


	
	int networkCount(0), failedNetworks(0);
	cout << endl;
	while(networkCount < numNetworks && failedNetworks < maxFailedNetworks){
		vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty));
		if(b.size() == targetServiceVolumes + 1){
			cout << "{" << failedNetworks << "}";
			cout.flush();
			if(networkCount == 50){
				addParent(b, 3, 6); // makes 8
				addParent(b, 2, 8); // makes 9
				addParent(b, 1, 7); // makes 10
				addParent(b, 9, 10); // makes 11
				addParent(b, 5, 11); // makes 12
				addParent(b, 4, 12); // makes 13;
				b[13].parentIndex = 0;
				b[0].childIndex.push_back(13);
				cout << "\n made tree: " << makeStringTreeOneLine(b) << endl;
				settleBifurcations(b, true, 1.0e-6, 10);
				//drawBifurcationTree("outputs/" + uniquer + "_initial.png", b, width, height, border);
				//double *bestBMeasures = new double[numBest];
				//vector<Bifurcation2D> *bestBs = new vector<Bifurcation2D>[numBest];
				//for(int i(0); i < numBest; i++)
				//	bestBMeasures[i] = -1.0;
				//exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, false);
				//delete[] bestBMeasures;
				//delete[] bestBs;
			}
			failedNetworks = 0;
			networkCount++;
		}else{
			failedNetworks++;
			//cout << "failed: b.size() = " << b.size() << " != " << targetServiceVolumes + 1 << endl;
		}
	}
	
}

int siblingIndex(vector<Bifurcation2D> &b, int i){
	if(b[i].parentIndex < 0)
		return -1;
	if(b[b[i].parentIndex].childIndex.size() != 2)
		return -1;
	if(b[b[i].parentIndex].childIndex[0] == i)
		return b[b[i].parentIndex].childIndex[1];
	return b[b[i].parentIndex].childIndex[0];
}

bool sibling_has_two_children(vector<Bifurcation2D> &b, int i){
	int sib(siblingIndex(b, i));
	if(sib < 0)
		return false;
	return b[sib].childIndex.size() == 2;
}

// assumes that a valid value for i is given and that i > 0 (i.e. i is not the root)
void swapWithSiblingChild(vector<Bifurcation2D> &b, int i, bool withFirstChild = true){
	int sib(siblingIndex(b, i));
	int cSwap(1);
	if(withFirstChild)
		cSwap = 0;
	int cIndex(b[sib].childIndex[cSwap]);
	for(unsigned int c(0); c < b[b[i].parentIndex].childIndex.size(); c++){
		if(b[b[i].parentIndex].childIndex[c] == i)
			b[b[i].parentIndex].childIndex[c] = cIndex;
	}
	for(unsigned int c(0); c < b[b[cIndex].parentIndex].childIndex.size(); c++){
		if(b[b[cIndex].parentIndex].childIndex[c] == cIndex)
			b[b[cIndex].parentIndex].childIndex[c] = i;
	}
	int tempP(b[i].parentIndex);
	b[i].parentIndex = b[cIndex].parentIndex;
	b[cIndex].parentIndex = tempP;
}

vector<int> findShortestDescendingPathToZero(vector<int> *descendingPaths, vector<int> *neighborRanks, int start){
	int bestNeighbor(-1);
	for(unsigned int i(0); i < neighborRanks[start].size(); i++){
		if(descendingPaths[neighborRanks[start][i]].size() > 0){
			if(bestNeighbor < 0)
				bestNeighbor = i;
			else if(descendingPaths[neighborRanks[start][bestNeighbor]].size() > descendingPaths[neighborRanks[start][i]].size())
				bestNeighbor = i;
		}
	}
	vector<int> path;
	if(bestNeighbor >= 0)
		path.insert(path.end(), descendingPaths[neighborRanks[start][bestNeighbor]].begin(), descendingPaths[neighborRanks[start][bestNeighbor]].end());
	if(bestNeighbor >= 0 || start == 0)
		path.push_back(start);
	return path;

	/*vector<int> path;
	path.push_back(start);
	if(start == 0)
		return path;
	vector<int> shortestPath;
	for(unsigned int i(0); i < neighborRanks[start].size(); i++){
		if(neighborRanks[start][i] >= start)
			continue;
		vector<int> trialPath(findShortestDescendingPathToZero(neighborRanks, neighborRanks[start][i]));
		//cout << "\n start = " << start << ", looking at neighbor at " << neighborRanks[start][i] << ", shortestPath.size() = " << shortestPath.size() << ", trialPath.size() = " << trialPath.size() << endl;
		if(shortestPath.size() < 1 || (trialPath.size() > 0 && trialPath.size() < shortestPath.size()))
			shortestPath = trialPath;
	}
	if(shortestPath.size() < 1)
		return shortestPath;
	path.insert(path.end(), shortestPath.begin(), shortestPath.end());
	return path;*/
}

void rankedExchangeConnectivity(){
	//kisset(24611, 19046, 7706);
	// 1 2 3  4  5   6    7      8      9       10
	// 1 1 3 15 105 945 10395 135135 2027025 34459425
	int maxFailedInsertions(100000), numBest(10395), maxFailedNetworks(10000), targetServiceVolumes(7);
	double nearestRadius(1.0);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	double totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0);
	string uniquer("outputs/exhaustiveTrendsMem_" + makeString(targetServiceVolumes) + "_" + makeString(totalLengthCoefficient) + "TL_" + makeString(avePathLengthCoefficient) + "AP_" + makeString(maxPathLengthCoefficient) + "MP");
	int width(120), height(120), border(10);

	cout << "\n exhaustiveTrendsMem()\n uniquer = " << uniquer
		<< "\n maxFailedNetworks = " << maxFailedNetworks
		<< "\n targetServiceVolumes = " << targetServiceVolumes
		<< "\n sizeX = " << sizeX << "\n sizeY = " << sizeY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n numBest = " << numBest
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< endl << endl;

	int failedNetworks(0);
	cout << endl;
	while(failedNetworks < maxFailedNetworks){
		vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty));
		if(b.size() == targetServiceVolumes + 1){
			cout << "{" << failedNetworks << "}";
			cout.flush();
			//drawBifurcationTree("outputs/" + uniquer + "_initial.png", b, width, height, border);
			vector<Bifurcation2D> *bestBs = new vector<Bifurcation2D>[numBest];
			double *bestBMeasures = new double[numBest];
			for(int i(0); i < numBest; i++)
				bestBMeasures[i] = -1.0;
			exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, true, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			cout << "\n Finding neighbor ranks..." << endl;
			vector<int> *neighborRanks = new vector<int>[numBest]; // the best is rank 0, the second best is rank 1, etc . . .
			for(int j(0); j < numBest; j++){
				for(unsigned int k(1); k < bestBs[j].size(); k++){
					if(sibling_has_two_children(bestBs[j], k)){
						vector<Bifurcation2D> tempB(bestBs[j]);
						swapWithSiblingChild(tempB, k);
						for(int m(j + 1); m < numBest; m++){
							if(equivalentTrees(tempB, bestBs[m])){
								neighborRanks[j].push_back(m);
								neighborRanks[m].push_back(j);
								m = numBest;
							}
						}
						tempB = bestBs[j];
						swapWithSiblingChild(tempB, k, false);
						for(int m(j + 1); m < numBest; m++){
							if(equivalentTrees(tempB, bestBs[m])){
								neighborRanks[j].push_back(m);
								neighborRanks[m].push_back(j);
								m = numBest;
							}
						}
					}
				}
				sort(neighborRanks[j].begin(), neighborRanks[j].end());
			}
			cout << "\n\t Found neighbor ranks." << endl;
			cout << "\n Finding best paths..." << endl;
			vector<int> *pathToBest = new vector<int>[numBest];
			pathToBest[0].push_back(0);
			queue<int> toCheck;
			toCheck.push(0);
			while(toCheck.size() > 0){
				int j(toCheck.front());
				toCheck.pop();
				for(unsigned int k(0); k < neighborRanks[j].size(); k++){
					int n(neighborRanks[j][k]);
					if(pathToBest[n].size() < 1){
						pathToBest[n] = pathToBest[j];
						pathToBest[n].push_back(n);
						toCheck.push(n);
					}
				}
			}
			cout << "\n\t Found best paths." << endl;
			cout << "\n Finding descending paths..." << endl;
			vector<int> *descendingPaths = new vector<int>[numBest];
			for(int j(0); j < numBest; j++){
				if(j%100 == 0)
					cout << "(d." << j << ")";
				descendingPaths[j] = findShortestDescendingPathToZero(descendingPaths, neighborRanks, j);
			}
			cout << "\n\t Found descending paths." << endl;
			int descendingCount(0), noDescendingCount(0);
			cout << "\n\n rank\t resc_rank\t tree\t measure\t norm_measure\t "
				"num_path_to_best_edges\t path_to_best\t "
				"num_descending_path_edges\t descending_path\t "
				"min_neighbor_rank\t max_neighbor_rank\t "
				"ave_neighbor_rank\t std_dev_neighbor_rank\t "
				"num_neighbors\t neighbor_ranks..." << endl;
			for(int j(0); j < numBest; j++){
				if(j < 100)
					drawBifurcationTree(uniquer + "_" + makeString(j) + ".png", bestBs[j], width, height, border);
				cout << j + 1 << "\t" << double(j + 1)/double(numBest + 1) << "\t" << makeStringTreeOneLine(bestBs[j]) << "\t" << bestBMeasures[j] << "\t" << bestBMeasures[j]/bestBMeasures[0]
					<< "\t" << pathToBest[j].size() - 1 << "\t" << makeVectorString(pathToBest[j])
					<< "\t" << (int)descendingPaths[j].size() - 1 << "\t" << makeVectorString(descendingPaths[j])
					<< "\t" << minimum(neighborRanks[j]) << "\t" << maximum(neighborRanks[j])
					<< "\t" << mean(neighborRanks[j]) << "\t" << stDev(neighborRanks[j])
					<< "\t" << neighborRanks[j].size() << "\t" << makeVectorString(neighborRanks[j]) << endl;
				if(descendingPaths[j].size() > 0)
					descendingCount++;
				else
					noDescendingCount++;
			}
			cout << endl << endl;
			cout << descendingCount << " of " << descendingCount + noDescendingCount
				<< " can descend to 0 (" << 100.0*double(descendingCount)/double(descendingCount + noDescendingCount) << "%)" << endl;
			cout << noDescendingCount << " of " << descendingCount + noDescendingCount
				<< " cannot descend to 0 (" << 100.0*double(noDescendingCount)/double(descendingCount + noDescendingCount) << "%)" << endl;
			delete[] bestBs;
			delete[] bestBMeasures;
			delete[] neighborRanks;
			delete[] pathToBest;
			delete[] descendingPaths;
			failedNetworks = maxFailedNetworks;
		} else
			failedNetworks++;
	}
}

void createSegAndAddChildren(vector<Bifurcation2D> &b, terseSeg** linList, int index, double rad, double betaRad, double capRad){
	linList[index] = new terseSeg;
	linList[index]->len = separation2D(b[index], b[b[index].parentIndex]);
	if(b[index].childIndex.size() == 0)
		linList[index]->rad = capRad;
	else
		linList[index]->rad = rad;
	linList[index]->flow = heartFlowRate;
	if(b[index].parentIndex == 0)
		linList[index]->parent = NULL;
	else
		linList[index]->parent = linList[b[index].parentIndex];
	linList[index]->name = makeString(index);
	linList[index]->nchild = b[index].childIndex.size();
	linList[index]->child = new terseSeg*[b[index].childIndex.size()];
	for(unsigned int c(0); c < b[index].childIndex.size(); c++){
		createSegAndAddChildren(b, linList, b[index].childIndex[c], betaRad*rad, betaRad, capRad);
		linList[index]->child[c] = linList[b[index].childIndex[c]];
	}
}

terseSeg** createTerseSegTree(vector<Bifurcation2D> &b, double rootRad = 0.01, double betaRad = 0.95, double capRad = 0.005, double thresh = 1.0e-6, string initialDrawFn = "createTerseSegTree.png"){
	if(b[0].childIndex.size() != 1){
		cout << "\n Error createTerseSegTree(): b[0] must have exactly one child." << endl;
		return NULL;
	}
	terseSeg** linList = new terseSeg*[b.size()];
	createSegAndAddChildren(b, linList, b[0].childIndex[0], rootRad, betaRad, capRad);
	terseSeg* root(getRoot(b, linList));
	organizeBinary(root);
	updateBinary(root);
	//cout << "\n createTerseSegTree(): initial tree is\n" << makeStringTree(root) << endl;
	int width(1920), height(1080), border(25),
		printEveryOptimization(0), maxOptimizationIts(10000);
	terseEnds *rootEnds = drawBinaryTreeByFlow(initialDrawFn, root, width, height, border);
	naivelyOptimizeWithConstantTipsAndRoot(rootEnds, thresh, maxOptimizationIts, printEveryOptimization, "nowhere", width, height, border);
	return linList;
}

void updateTerseSegs(vector<Bifurcation2D> &b, terseSeg** linList, double thresh = 1.0e-6, string initialDrawFn = "updateTerseSegs.png"){
	for(unsigned int i(1); i < b.size(); i++)
		linList[i]->len = separation2D(b[i], b[b[i].parentIndex]);
	int width(1920), height(1080), border(25),
		printEveryOptimization(0), maxOptimizationIts(10000);
	terseEnds *rootEnds = drawBinaryTreeByFlow(initialDrawFn, getRoot(b, linList), width, height, border);
	naivelyOptimizeWithConstantTipsAndRoot(rootEnds, thresh, maxOptimizationIts, printEveryOptimization, "nowhere", width, height, border);
}

void testTerseSeg(){
	kisset(123, 4567, 89);
	int maxFailedInsertions(10000), targetServiceVolumes(6), width(1920), height(width), border(25);
	double nearestRadius(1.0);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	vector<Bifurcation2D> b(initializeEvenRandomSpacingCircle(maxFailedInsertions, sizeX, nearestRadius));
	introduceBifurcationsFromSources(b);
	vector<BranchPoint2D> bc(consolidateDegenerates(b));
	terseSeg** terseB = createTerseSegTree(b.size(), bc);

	drawBifurcationTree("outputs/testTerseSeg_Bifurcation2DTree.png", b, width, height, border);
	drawTree("outputs/testTerseSeg_BranchPoint2DTree.png", bc, width, height, border);
	terseEnds *rootEnds = drawBinaryTreeByFlow("outputs/testTerseSeg_terseSegTree.png", getRoot(bc, terseB), width, height, border);

	cout << "\n\n Bifurcation2D tree:\n" << makeStringTree(b) << endl;
	cout << "\n\n BranchPoint2D tree:\n" << makeStringTree(bc) << endl;
	cout << "\n terseSeg tree:\n" << makeStringTree(getRoot(bc, terseB)) << endl;
	cout << "\n terseEnds tree:\n" << makeStringTree(rootEnds) << endl;
}

// assumes values are strictly ascending
double binLowerBound(int size, int i, double *values){
	if(size < 2 || i >= size)// cannot handle a single point (or no point)
		return -1.0;
	if(i == 0)
		return values[0] - ((values[1] - values[0])/2.0);
	return (values[i - 1] + values[i])/2.0;
}

// assumes values are strictly ascending
double binUpperBound(int size, int i, double *values){
	if(size < 2 || i >= size)// cannot handle a single point (or no point)
		return -1.0;
	if(i == size - 1)
		return values[size - 1] + ((values[size - 1] - values[size - 2])/2.0);
	return (values[i] + values[i + 1])/2.0;
}

double binWidth(int size, int i, double *values){
	return binUpperBound(size, i, values) - binLowerBound(size, i, values);
}

// assumes a normalized distribution and values are strictly ascending
double mean(int size, double *values, double *density){
	if(size < 1)
		return 0.0;
	if(size < 2)
		return values[0];
	double m(0.0);
	for(int i(0); i < size; i++)
		m += values[i]*density[i]*binWidth(size, i, values);
	return m;
}

// assumes a normalized distribution and values are strictly ascending
double variance(int size, double *values, double *density){
	if(size < 2)
		return 0.0;
	double v(0.0);
	for(int i(0); i < size; i++)
		v += values[i]*values[i]*density[i]*binWidth(size, i, values);
	double m(mean(size, values, density));
	return v - (m*m);
}

// assumes values are strictly ascending
double* normalizedDensity(int size, double *values, double *density){
	double *normd = new double[size];
	if(size < 2) // cannot handle single point (or no point)
		return normd;
	double sum(0.0);
	for(int i(0); i < size; i++)
		sum += density[i]*binWidth(size, i, values);
	if(sum == 0.0)
		return normd;
	for(int i(0); i < size; i++)
		normd[i] = density[i]/sum;
	return normd;
}

// assumes a normalized distribution and values are strictly ascending
map<double, double> zeroMeanAndUnitVariance(int size, double *values, double *density){
	map<double, double> rescaled;
	double m(mean(size, values, density)), v(variance(size, values, density));
	if(v <= 0.0)
		return rescaled;
	for(int i(0); i < size; i++)
		rescaled[(values[i] - m)/sqrt(v)] = density[i]*sqrt(v);
	return rescaled;
}

void testRescaling(){
	int numPoints(10);
	double *x = new double[numPoints];
	double *y = new double[numPoints];
	for(int i(0); i < numPoints; i++){
		x[i] = 0.1*(i - (numPoints/3));
		y[i] = 1 + numPoints/2 - abs(numPoints/2 - i) + rkiss();
	}
	double *ny = normalizedDensity(numPoints, x, y);
	cout << "\n\n i\t x\t y\t ny" << endl;
	for(int i(0); i < numPoints; i++)
		cout << i << "\t" << x[i] << "\t" << y[i] << "\t" << ny[i] << endl;
	double m(mean(numPoints, x, ny));
	cout << "\n mean x-y = " << m << endl;
	double v(variance(numPoints, x, ny));
	cout << " variance x-y = " << v << endl;
	map<double, double> rxy = zeroMeanAndUnitVariance(numPoints, x, ny);
	double *rx = new double[numPoints];
	double *ry = new double[numPoints];
	int index(0);
	cout << "\n\n i\t rx\t ry" << endl;
	for(map<double, double>::iterator it(rxy.begin()); it != rxy.end(); it++){
		rx[index] = it->first;
		ry[index] = it->second;
		cout << index << "\t" << rx[index] << "\t" << ry[index] << endl;
		index++;
	}
	double rm(mean(numPoints, rx, ry));
	cout << "\n rescaled mean = " << rm << endl;
	double rv(variance(numPoints, rx, ry));
	cout << " rescaled variance = " << rv << endl;
	delete[] x;
	delete[] y;
	delete[] ny;
	delete[] rx; 
	delete[] ry;
}

// assumes x is strictly ascending
double* localDerivatives(int size, double *x, double *y){
	double *derivative = new double[size];
	if(size < 2) // cannot handle a single point (or no point)
		return derivative;
	derivative[0] = (y[1] - y[0])/(x[1] - x[0]);
	for(int i(1); i < size - 1; i++)
		derivative[i] = (y[i + 1] - y[i - 1])/(x[i + 1] - x[i - 1]);
	derivative[size - 1] = (y[size - 1] - y[size - 2])/(x[size - 1] - x[size - 2]);
	return derivative;
}

double linearInterpolation(double xLow, double x0, double xHigh, double yLow, double yHigh){
	if(x0 < xLow || x0 > xHigh)
		return 0.0;
	if(x0 == xHigh)
		return yHigh;
	double frac((x0 - xLow)/(xHigh - xLow));
	return yLow + frac*(yHigh - yLow);
}

void rescaledDerivatives(){
	kisset(15895, 26460, 31928);
	// kisset(24611, 19046, 7706);
	// 1 2 3  4  5   6    7      8      9       10
	// 1 1 3 15 105 945 10395 135135 2027025 34459425
	int maxFailedInsertions(100000), numBest(10395), numNetworks(300), maxFailedNetworks(10000), targetServiceVolumes(7);
	double nearestRadius(1.0);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	double totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0);
	string uniquer("rescaledDerivatives" + makeString(sizeX));
	int width((int)floor(sizeX)*10), height((int)floor(sizeY)*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300/2;
	border = 10/2;

	cout << "\n rescaledDerivatives()\n uniquer = " << uniquer
		<< "\n numNetworks = " << numNetworks << "\n maxFailedNetworks = " << maxFailedNetworks
		<< "\n targetServiceVolumes = " << targetServiceVolumes
		<< "\n sizeX = " << sizeX << "\n sizeY = " << sizeY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n numBest = " << numBest
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< endl << endl;

	int networkCount(0), failedNetworks(0);
	double **bestBMeasures = new double*[numNetworks];
	vector<Bifurcation2D> **bestBs = new vector<Bifurcation2D>*[numNetworks];
	cout << endl;
	while(networkCount < numNetworks && failedNetworks < maxFailedNetworks){
		vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty));
		if(b.size() == targetServiceVolumes + 1){
			cout << "{" << failedNetworks << "}";
			cout.flush();
			//drawBifurcationTree("outputs/" + uniquer + "_initial.png", b, width, height, border);
			bestBs[networkCount] = new vector<Bifurcation2D>[numBest];
			bestBMeasures[networkCount] = new double[numBest];
			for(int i(0); i < numBest; i++)
				bestBMeasures[networkCount][i] = -1.0;
			exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs[networkCount], bestBMeasures[networkCount], false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			failedNetworks = 0;
			networkCount++;
		}else{
			failedNetworks++;
			//cout << "failed: b.size() = " << b.size() << " != " << targetServiceVolumes + 1 << endl;
		}
	}
	double *ranks = new double[numBest];
	for(int i(0); i < numBest; i++)
		ranks[i] = (i + 1.0)/numBest;
	map<double, double> *rescaleds = new map<double, double>[numNetworks];
	double minRescaledMeasure(0.0), maxRescaledMeasure(0.0);
	for(int n(0); n < numNetworks; n++){
		double *derivative = localDerivatives(numBest, ranks, bestBMeasures[n]);
		double *normD = normalizedDensity(numBest, bestBMeasures[n], derivative);
		delete[] derivative;
		rescaleds[n] = zeroMeanAndUnitVariance(numBest, bestBMeasures[n], normD);
		delete[] normD;
		if(n == 0){
			minRescaledMeasure = rescaleds[0].begin()->first;
			maxRescaledMeasure = rescaleds[0].rbegin()->first;
		}else{
			if(minRescaledMeasure > rescaleds[n].begin()->first)
				minRescaledMeasure = rescaleds[n].begin()->first;
			if(maxRescaledMeasure < rescaleds[n].rbegin()->first)
				maxRescaledMeasure = rescaleds[n].rbegin()->first;
		}
	}
	for(int i(0); i < numBest; i++)
		ranks[i] = minRescaledMeasure + (i*(maxRescaledMeasure - minRescaledMeasure)/(numBest - 1.0));
	map<double, double>::iterator *its = new map<double, double>::iterator[numNetworks];
	for(int n(0); n < numNetworks; n++)
		its[n] = rescaleds[n].begin();
	cout << "\n\n rescaled_rank\t ave_drr/dm\t +/-" << endl;
	for(int i(0); i < numBest; i++){
		vector<double> drrdm;
		for(int n(0); n < numNetworks; n++){
			while(next(its[n])->first < ranks[i] && next(its[n]) != rescaleds[n].end())
				its[n]++;
			double interp(0.0);
			if(its[n]->first <= ranks[i] && ranks[i] <= next(its[n])->first)
				interp = linearInterpolation(its[n]->first, ranks[i], next(its[n])->first, its[n]->second, next(its[n])->second);
			drrdm.push_back(interp);
		}
		cout << ranks[i] << "\t" << mean(drrdm) << "\t" << stDev(drrdm) << endl;
	}
	delete[] its;
	cout << "\n\n rank\t ave_measure\t +/-\t ave_norm_measure\t +/-" << endl;
	for(int i(0); i < numBest; i++){
		vector<double> measures, normMeasures;
		for(int n(0); n < networkCount; n++){
			if(bestBMeasures[n][i] >= 0.0){
				measures.push_back(bestBMeasures[n][i]);
				normMeasures.push_back(bestBMeasures[n][i]/bestBMeasures[n][0]);
			}
		}
		if(measures.size() > 0){
			cout << i + 1 << "\t" << mean(measures) << "\t" << stDev(measures)
				<< "\t" << mean(normMeasures) << "\t" << stDev(normMeasures) << endl;
		} else
			cout << i + 1 << "\t X\t X\t X\t X" << endl;
	}
	cout << "\n\n measures..." << endl;
	for(int i(0); i < numBest; i++){
		cout << i + 1;
		for(int n(0); n < networkCount; n++){
			if(bestBMeasures[n][i] >= 0.0)
				cout << "\t" << makeStringTreeOneLine(bestBs[n][i]) << "\t" << bestBMeasures[n][i];
			else
				cout << "\t X\t X";
		}
		cout << endl;
	}
	cout << "\n\n normalized_measures..." << endl;
	for(int i(0); i < numBest; i++){
		cout << i + 1;
		for(int n(0); n < networkCount; n++){
			if(bestBMeasures[n][i] >= 0.0)
				cout << "\t" << makeStringTreeOneLine(bestBs[n][i]) << "\t" << bestBMeasures[n][i]/bestBMeasures[n][0];
			else
				cout << "\tX \tX";
		}
		cout << endl;
	}
	for(int n(0); n < numNetworks; n++){
		delete[] bestBs[n];
		delete[] bestBMeasures[n];
	}

	delete[] bestBs;
	delete[] bestBMeasures;
	delete[] rescaleds;
	delete[] ranks;
}

double nthMapKey(map<double, double> m, int n){
	map<double, double>::iterator it(m.begin());
	for(int i(0); i < n; i++)
		it++;
	return it->first;
}

void rescaledDerivativesAndAnalysis(){
	kisset(15895, 26460, 31928);
	//kisset(24611, 19046, 7706);
	// 1 2 3  4  5   6    7      8      9       10
	// 1 1 3 15 105 945 10395 135135 2027025 34459425
	int maxFailedInsertions(100000), numBest(945), numNetworks(100), maxFailedNetworks(10000), targetServiceVolumes(6);
	double nearestRadius(1.0);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0), lengthThresh(1.0e-6);
	double totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0);
	string uniquer("outputs/reDerAA_" + makeString(targetServiceVolumes) + "_n" + makeString(numNetworks) + "_");
	int width((int)floor(sizeX)*10), height((int)floor(sizeY)*100), border(width/25), numBins(100), numToDraw(5);
	if(numToDraw > numNetworks)
		numToDraw = numNetworks;
	if(numToDraw > numBest)
		numToDraw = numBest;
	if(numBins < 5)
		numBins = 5;
	if(height < width)
		border = height/25;
	width = height = 300/2;
	border = 10/2;

	cout << "\n rescaledDerivatives()\n uniquer = " << uniquer
		<< "\n numNetworks = " << numNetworks << "\n maxFailedNetworks = " << maxFailedNetworks
		<< "\n targetServiceVolumes = " << targetServiceVolumes
		<< "\n lengthThresh = " << lengthThresh
		<< "\n sizeX = " << sizeX << "\n sizeY = " << sizeY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n numBest = " << numBest << "\n numBins = " << numBins
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< endl << endl;

	vector<double> **validLambdaLs = new vector<double>*[numNetworks];
	vector<double> **validGammas = new vector<double>*[numNetworks];
	int **ignoredCount = new int*[numNetworks];

	int networkCount(0), failedNetworks(0);
	double **bestBMeasures = new double*[numNetworks];
	vector<Bifurcation2D> **bestBs = new vector<Bifurcation2D>*[numNetworks];
	cout << endl;
	while(networkCount < numNetworks && failedNetworks < maxFailedNetworks){
		vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty));
		if(b.size() == targetServiceVolumes + 1){
			cout << "{" << failedNetworks << "}";
			cout.flush();
			bestBs[networkCount] = new vector<Bifurcation2D>[numBest];
			bestBMeasures[networkCount] = new double[numBest];
			for(int i(0); i < numBest; i++)
				bestBMeasures[networkCount][i] = -1.0;
			exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs[networkCount], bestBMeasures[networkCount], false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			if(networkCount < numToDraw){
				drawBifurcationTree(uniquer + makeString(networkCount) + "_initial.png", b, width, height, border);
				for(int draw(0); draw < numToDraw; draw++)
					drawBifurcationTree(uniquer + makeString(networkCount) + "_rank" + makeString(draw) + ".png", bestBs[networkCount][draw], width, height, border);

			}
			validLambdaLs[networkCount] = new vector<double>[numBest];
			validGammas[networkCount] = new vector<double>[numBest];
			ignoredCount[networkCount] = new int[numBest];
			for(int i(0); i < numBest; i++){
				validLambdaLs[networkCount][i] = validBranchingLengthRatios(bestBs[networkCount][i], ignoredCount[networkCount][i], lengthThresh);
				validGammas[networkCount][i] = validParentChildLengthRatios(bestBs[networkCount][i], ignoredCount[networkCount][i], lengthThresh);
			}
			failedNetworks = 0;
			networkCount++;
		}else{
			failedNetworks++;
			//cout << "failed: b.size() = " << b.size() << " != " << targetServiceVolumes + 1 << endl;
		}
	}

	string longFn(uniquer + "lambdaLbyRank.dat");
	ofstream lambdaLByRankAveragedOverRealizationsOut(longFn.c_str());
	cout << "\n lambda_L distribution by ranks averaged over realizations\n [lambda_L\t P(lambda_L)] . . ." << endl;
	double **lambda_LByRanks = new double*[numBest];
	double lambda_LMin(validLambdaLs[0][0][0]), lambda_LMax(validLambdaLs[0][0][0]);
	for(int i(0); i < numBest; i++){
		for(int n(0); n < numNetworks; n++){
			double tempMin(minimum(validLambdaLs[n][i])), tempMax(maximum(validLambdaLs[n][i]));
			if(lambda_LMin > tempMin)
				lambda_LMin = tempMin;
			if(lambda_LMax < tempMax)
				lambda_LMax = tempMax;
		}
	}
	double *lambdaLBinBottoms = getBinBottoms(lambda_LMin, lambda_LMax, numBins);
	for(int i(0); i < numBest; i++){
		vector<double> lambdaLTemp;
		for(int n(0); n < numNetworks; n++)
			lambdaLTemp.insert(lambdaLTemp.end(), validLambdaLs[n][i].begin(), validLambdaLs[n][i].end());
		if(i == 0)
			cout << "\n The first distribution is composed of " << lambdaLTemp.size() << " points." << endl;
		lambda_LByRanks[i] = probabilityDistribution(lambdaLTemp, lambdaLBinBottoms, numBins);
	}
	lambdaLByRankAveragedOverRealizationsOut << "lambda_L";
	for(int i(0); i < numBest; i++)
		lambdaLByRankAveragedOverRealizationsOut << "\t P(l_L_" << i << ")";
	lambdaLByRankAveragedOverRealizationsOut << endl;
	double lambdaLBinInc(lambdaLBinBottoms[1] - lambdaLBinBottoms[0]);
	for(int b(0); b < numBins; b++){
		lambdaLByRankAveragedOverRealizationsOut << lambdaLBinBottoms[b] + lambdaLBinInc/2.0;
		for(int i(0); i < numBest; i++){
			lambdaLByRankAveragedOverRealizationsOut << "\t" << lambda_LByRanks[i][b];
		}
		lambdaLByRankAveragedOverRealizationsOut << endl;
	}
	lambdaLByRankAveragedOverRealizationsOut.close();
	cout << "\n lambda_L distribution over realizations for all ranks\n lambda_L\t P(lambda_L)" << endl;
	for(int n(0); n < numNetworks; n++){

	}
	cout << "\n gamma distribution over ranks for many realizations\n lambda_L\t P(gamma)" << endl;
	for(int i(0); i < numBest; i++){

	}
	cout << "\n gamma distribution over realizations for all ranks\n lambda_L\t P(gamma)" << endl;
	for(int n(0); n < numNetworks; n++){

	}
	cout << "\n gamma distribution over ranks for many realizations\n lambda_L\t P(gamma)" << endl;
	for(int i(0); i < numBest; i++){

	}
	string ignoreCountFn(uniquer + "ignoreCount.dat");
	ofstream ignoreCountOut(ignoreCountFn.c_str());
	ignoreCountOut << "ignoredCounts:\n\t mean";
	cout << "\n\n rank\t ave_ignoredCounts" << endl;
	for(int n(0); n < numNetworks; n++)
		ignoreCountOut << "\t n" << n;
	ignoreCountOut << endl << "mean";
	for(int n(0); n < numNetworks; n++){
		vector<int> ignoredCountTemp;
		for(int i(0); i < numBest; i++)
			ignoredCountTemp.push_back(ignoredCount[n][i]);
		ignoreCountOut << "\t" << mean(ignoredCountTemp);
	}
	ignoreCountOut << endl;
	for(int i(0); i < numBest; i++){
		ignoreCountOut << i;
		cout << i;
		vector<int> ignoredCountTemp;
		for(int n(0); n < numNetworks; n++)
			ignoredCountTemp.push_back(ignoredCount[n][i]);
		ignoreCountOut << "\t" << mean(ignoredCountTemp);
		cout << "\t" << mean(ignoredCountTemp) << endl;
		for(int n(0); n < numNetworks; n++)
			ignoreCountOut << "\t" << ignoredCount[n][i];
		ignoreCountOut << endl;
	}
	

	double *ranks = new double[numBest];
	for(int i(0); i < numBest; i++)
		ranks[i] = (i + 1.0)/numBest;
	map<double, double> *rescaleds = new map<double, double>[numNetworks];
	double minRescaledMeasure(0.0), maxRescaledMeasure(0.0);
	for(int n(0); n < numNetworks; n++){
		double *derivative = localDerivatives(numBest, ranks, bestBMeasures[n]);
		double *normD = normalizedDensity(numBest, bestBMeasures[n], derivative);
		delete[] derivative;
		rescaleds[n] = zeroMeanAndUnitVariance(numBest, bestBMeasures[n], normD);
		delete[] normD;
		if(n == 0){
			minRescaledMeasure = rescaleds[0].begin()->first;
			maxRescaledMeasure = rescaleds[0].rbegin()->first;
		}else{
			if(minRescaledMeasure > rescaleds[n].begin()->first)
				minRescaledMeasure = rescaleds[n].begin()->first;
			if(maxRescaledMeasure < rescaleds[n].rbegin()->first)
				maxRescaledMeasure = rescaleds[n].rbegin()->first;
		}
	}
	ignoreCountOut.close();
	for(int i(0); i < numBest; i++)
		ranks[i] = minRescaledMeasure + (i*(maxRescaledMeasure - minRescaledMeasure)/(numBest - 1.0));
	map<double, double>::iterator *its = new map<double, double>::iterator[numNetworks];
	for(int n(0); n < numNetworks; n++)
		its[n] = rescaleds[n].begin();
	string drrdmFn(uniquer + "drrdm.dat");
	ofstream drrdmOut(drrdmFn.c_str());
	drrdmOut << "\n\n rescaled_rank\t ave_drr/dm\t +/-" << endl;
	for(int i(0); i < numBest; i++){
		vector<double> drrdm;
		for(int n(0); n < numNetworks; n++){
			while(next(its[n])->first < ranks[i] && next(its[n]) != rescaleds[n].end())
				its[n]++;
			double interp(0.0);
			if(its[n]->first <= ranks[i] && ranks[i] <= next(its[n])->first)
				interp = linearInterpolation(its[n]->first, ranks[i], next(its[n])->first, its[n]->second, next(its[n])->second);
			drrdm.push_back(interp);
		}
		drrdmOut << ranks[i] << "\t" << mean(drrdm) << "\t" << stDev(drrdm) << endl;
	}
	drrdmOut.close();
	delete[] its;
	string aveMeasuresFn(uniquer + "aveMeasuresOut.dat");
	ofstream aveMeasuresOut(aveMeasuresFn.c_str());
	aveMeasuresOut << "\n\n rank\t ave_measure\t +/-\t ave_norm_measure\t +/-" << endl;
	for(int i(0); i < numBest; i++){
		vector<double> measures, normMeasures;
		for(int n(0); n < networkCount; n++){
			if(bestBMeasures[n][i] >= 0.0){
				measures.push_back(bestBMeasures[n][i]);
				normMeasures.push_back(bestBMeasures[n][i]/bestBMeasures[n][0]);
			}
		}
		if(measures.size() > 0){
			aveMeasuresOut << i + 1 << "\t" << mean(measures) << "\t" << stDev(measures)
				<< "\t" << mean(normMeasures) << "\t" << stDev(normMeasures) << endl;
		}else
			aveMeasuresOut << i + 1 << "\t X\t X\t X\t X" << endl;
	}
	aveMeasuresOut.close();
	string measuresFn(uniquer + "measures.dat");
	ofstream measuresOut(measuresFn.c_str());
	measuresOut << "\n\n measures..." << endl;
	for(int i(0); i < numBest; i++){
		measuresOut << i + 1;
		for(int n(0); n < networkCount; n++){
			if(bestBMeasures[n][i] >= 0.0)
				measuresOut << "\t" << makeStringTreeOneLine(bestBs[n][i]) << "\t" << bestBMeasures[n][i];
			else
				measuresOut << "\t X\t X";
		}
		measuresOut << endl;
	}
	measuresOut.close();
	string normMeasuresFn(uniquer + "normMeasures.dat");
	ofstream normMeasuresOut(normMeasuresFn.c_str());
	normMeasuresOut << "\n\n normalized_measures..." << endl;
	for(int i(0); i < numBest; i++){
		normMeasuresOut << i + 1;
		for(int n(0); n < networkCount; n++){
			if(bestBMeasures[n][i] >= 0.0)
				normMeasuresOut << "\t" << makeStringTreeOneLine(bestBs[n][i]) << "\t" << bestBMeasures[n][i]/bestBMeasures[n][0];
			else
				normMeasuresOut << "\tX \tX";
		}
		normMeasuresOut << endl;
	}
	normMeasuresOut.close();
	for(int n(0); n < numNetworks; n++){
		delete[] bestBs[n];
		delete[] bestBMeasures[n];
		delete[] validLambdaLs[n];
		delete[] validGammas[n];
		delete[] ignoredCount[n];
	}

	delete[] bestBs;
	delete[] bestBMeasures;
	delete[] validLambdaLs;
	delete[] validGammas;
	delete[] ignoredCount;
	delete[] rescaleds;
	delete[] ranks;
}

void capillariesPerArea(){
	double sideStart(1.4), sideInc(0.02), sideEnd(4.0), nearestRadius(1.0);
	int numNetworks(1000), maxFailedInsertions(100000);
	cout << "\n\n sideStart = " << sideStart << "\n sideInc = " << sideInc << "\n sideEnd = " << sideEnd
		<< "\n nearestRadius = " << nearestRadius << "\n numNetworks = " << numNetworks
		<< "\n maxFailedInsertions = " << maxFailedInsertions
		<< endl;

	cout << "\n\n side\t <numTips>\t +/-\t area\t tips/area" << endl;
	for(double side(sideStart); side < sideEnd + sideInc/2.0; side += sideInc){
		vector<unsigned int> numTips;
		double heartx(side/2.0), hearty(nearestRadius/2.0);
		for(int n(0); n < numNetworks; n++){
			vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, side, side, nearestRadius, heartx, hearty));
			numTips.push_back(b.size());
		}
		double m(mean(numTips));
		cout << side << "\t" << m << "\t" << stDev(numTips) << "\t" << side*side << "\t" << m/side/side << endl;
	}
}

void capillariesPerVolume(){
	double sideStart(2.0), sideInc(0.5), sideEnd(10.0), nearestRadius(1.0);
	int numNetworks(100), maxFailedInsertions(100000);
	cout << "\n\n sideStart = " << sideStart << "\n sideInc = " << sideInc << "\n sideEnd = " << sideEnd
		<< "\n nearestRadius = " << nearestRadius << "\n numNetworks = " << numNetworks
		<< "\n maxFailedInsertions = " << maxFailedInsertions
		<< endl;

	cout << "\n\n side\t <numTips>\t +/-\t vol\t tips/vol" << endl;
	for(double side(sideStart); side < sideEnd + sideInc/2.0; side += sideInc){
		vector<unsigned int> numTips;
		double heartx(side/2.0), hearty(nearestRadius/2.0), heartz(side/2.0);
		for(int n(0); n < numNetworks; n++){
			vector<Bifurcation3D> b(initializeEvenRandomSpacing(maxFailedInsertions, side, side, side, nearestRadius, heartx, hearty, heartz));
			numTips.push_back(b.size());
		}
		double m(mean(numTips));
		cout << side << "\t" << m << "\t" << stDev(numTips) << "\t" << side*side*side << "\t" << m/side/side/side << endl;
	}
}

// assumes the heart is 0
vector<Bifurcation2D> purelyDescendingSearch(vector<vector<Bifurcation2D> > &history, unsigned int numToContinue, vector<Bifurcation2D> &b, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient,
	bool print = false, string uniquer = "", int width = 300, int height = 300, int border = 12){
	for(unsigned int i(0); i < history.size(); i++){
		if(equivalentTrees(b, history[i]))
			return b;
	}
	history.push_back(b);
	double initialMeasure(tripletMeasure(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	double bestMeasure(initialMeasure);
	vector<Bifurcation2D> bestB(b);
	vector<vector<Bifurcation2D> > helpfulSwaps;
	for(unsigned int i(1); i < b.size(); i++){ // skip the heart at 0
		if(sibling_has_two_children(b, i)){
			vector<Bifurcation2D> bCopy(b);
			swapWithSiblingChild(bCopy, i);
			settleBifurcations(bCopy);
			double tempMeasure(tripletMeasure(bCopy, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
			//cout << "/ts.tM=" << tempMeasure << "/";
			if(tempMeasure < initialMeasure)
				helpfulSwaps.push_back(bCopy);
			bCopy = b;
			swapWithSiblingChild(bCopy, i, false);
			settleBifurcations(bCopy);
			tempMeasure = tripletMeasure(bCopy, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			if(tempMeasure < initialMeasure)
				helpfulSwaps.push_back(bCopy);
		}
	}
	while(helpfulSwaps.size() > numToContinue)
		helpfulSwaps.erase(helpfulSwaps.begin() + (kiss()%helpfulSwaps.size()));
	for(unsigned int s(0); s < helpfulSwaps.size(); s++){
		vector<Bifurcation2D> tempB = purelyDescendingSearch(history, numToContinue, helpfulSwaps[s], totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
		double tempMeasure = tripletMeasure(tempB, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
		if(bestMeasure > tempMeasure){
			if(print){
				cout << "\n purelyDescendingSearch(): accepting swap that improves bestMeasure from "
					<< bestMeasure << " to " << tempMeasure << "\n\t tree: " << makeStringTreeOneLine(tempB) << endl;
			}
			bestMeasure = tempMeasure;
			bestB = tempB;
		}
	}
	return bestB;
}

// assumes the heart is 0
vector<BranchPoint2D> purelyDescendingConsolidatedSearch(vector<vector<Bifurcation2D> > &history, unsigned int numToContinue, vector<Bifurcation2D> &b, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient,
	double settleThresh = 1.0e-6, double consolidateThresh = 0.01, bool print = false, string uniquer = "", int width = 300, int height = 300, int border = 12){
	for(unsigned int i(0); i < history.size(); i++){
		if(equivalentTrees(b, history[i]))
			return consolidateDegenerates(b);
	}
	history.push_back(b);
	vector<BranchPoint2D> bc(consolidateDegenerates(b, consolidateThresh));
	double initialMeasure(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	double bestMeasure(initialMeasure);
	vector<BranchPoint2D> bestB(bc);
	vector<vector<Bifurcation2D> > helpfulSwaps;
	for(unsigned int i(1); i < b.size(); i++){ // skip the heart at 0
		if(sibling_has_two_children(b, i)){
			vector<Bifurcation2D> bCopy(b);
			swapWithSiblingChild(bCopy, i);
			settleBifurcations(bCopy, false, settleThresh);
			vector<BranchPoint2D> bc(consolidateDegenerates(bCopy, consolidateThresh));
			double tempMeasure(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
			//cout << "/ts.tM=" << tempMeasure << "/";
			if(tempMeasure < initialMeasure)
				helpfulSwaps.push_back(bCopy);
			bCopy = b;
			swapWithSiblingChild(bCopy, i, false);
			settleBifurcations(bCopy, false, settleThresh);
			bc = consolidateDegenerates(bCopy, consolidateThresh);
			tempMeasure = tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			if(tempMeasure < initialMeasure)
				helpfulSwaps.push_back(bCopy);
		}
	}
	while(helpfulSwaps.size() > numToContinue)
		helpfulSwaps.erase(helpfulSwaps.begin() + (kiss()%helpfulSwaps.size()));
	for(unsigned int s(0); s < helpfulSwaps.size(); s++){
		vector<BranchPoint2D> tempB = purelyDescendingConsolidatedSearch(history, numToContinue, helpfulSwaps[s], totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, print, uniquer, width, height, border);
		double tempMeasure = tripletMeasure(tempB, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
		if(bestMeasure > tempMeasure){
			if(print){
				cout << "\n purelyDescendingSearch(): accepting swap that improves bestMeasure from "
					<< bestMeasure << " to " << tempMeasure << "\n\t tree: " << makeStringTreeOneLine(tempB) << endl;
			}
			bestMeasure = tempMeasure;
			bestB = tempB;
		}
	}
	return bestB;
}

// assumes the heart is 0
vector<Bifurcation2D> purelyDescendingSearchOLD2(vector<Bifurcation2D> &b, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient,
		bool print = false, string uniquer = "", int width = 300, int height = 300, int border = 12){
	double initialMeasure(tripletMeasure(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	double bestMeasure(initialMeasure);
	vector<Bifurcation2D> bestB(b);
	for(unsigned int i(1); i < b.size(); i++){ // skip the heart at 0
		if(sibling_has_two_children(b, i)){
			vector<Bifurcation2D> bCopy(b);
			swapWithSiblingChild(bCopy, i);
			settleBifurcations(bCopy);
			double tempMeasure(tripletMeasure(bCopy, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
			//cout << "/ts.tM=" << tempMeasure << "/";
			if(tempMeasure < initialMeasure){
				//cout << "\nav" + makeStringTreeOneLine(bCopy);  cout.flush();
				vector<Bifurcation2D> tempB = purelyDescendingSearchOLD2(bCopy, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
				//cout << "a^"; cout.flush();
				tempMeasure = tripletMeasure(tempB, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
				if(bestMeasure > tempMeasure){
					if(print){
						cout << "\n purelyDescendingSearch(): accepting swap of " << i << " that improves bestMeasure from "
							<< bestMeasure << " to " << tempMeasure << "\n\t tree: " << makeStringTreeOneLine(tempB) << endl;
					}
					bestMeasure = tempMeasure;
					bestB = tempB;
				}
			}
			bCopy = b;
			swapWithSiblingChild(bCopy, i, false);
			settleBifurcations(bCopy);
			tempMeasure = tripletMeasure(bCopy, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			if(tempMeasure < initialMeasure){
				//cout << "\nbv" + makeStringTreeOneLine(bCopy);  cout.flush();
				vector<Bifurcation2D> tempB = purelyDescendingSearchOLD2(bCopy, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
				//cout << "b^"; cout.flush();
				tempMeasure = tripletMeasure(tempB, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
				if(bestMeasure > tempMeasure){
					if(print){
						cout << "\n purelyDescendingSearch(): accepting swap of " << i << " that improves bestMeasure from "
							<< bestMeasure << " to " << tempMeasure << "\n\t tree: " << makeStringTreeOneLine(tempB) << endl;
					}
					bestMeasure = tempMeasure;
					bestB = tempB;
				}
			}
		}
	}
	return bestB;
}

// assumes the heart is 0
vector<Bifurcation2D> purelyDescendingSearchOLD(vector<vector<Bifurcation2D> > &history, vector<Bifurcation2D> &b, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient,
	bool print = false, string uniquer = "", int width = 300, int height = 300, int border = 12){
	double initialMeasure(tripletMeasure(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	double bestMeasure(initialMeasure);
	vector<Bifurcation2D> bestB(b);
	for(unsigned int i(1); i < b.size(); i++){ // skip the heart at 0
		if(sibling_has_two_children(b, i)){
			vector<Bifurcation2D> bCopy(b);
			swapWithSiblingChild(bCopy, i);
			settleBifurcations(bCopy);
			double tempMeasure(tripletMeasure(bCopy, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
			//cout << "/ts.tM=" << tempMeasure << "/";
			if(tempMeasure < initialMeasure){
				//cout << "\nav" + makeStringTreeOneLine(bCopy);  cout.flush();
				if(history.size() > 0){
					for(unsigned int i(0); i < history.size() - 1; i++){
						if(equivalentTrees(history[i], bCopy))
							continue;
					}
				}
				history.push_back(bCopy);
				vector<Bifurcation2D> tempB = purelyDescendingSearchOLD(history, bCopy, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
				//cout << "a^"; cout.flush();
				tempMeasure = tripletMeasure(tempB, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
				if(bestMeasure > tempMeasure){
					if(print){
						cout << "\n purelyDescendingSearch(): accepting swap of " << i << " that improves bestMeasure from "
							<< bestMeasure << " to " << tempMeasure << "\n\t tree: " << makeStringTreeOneLine(tempB) << endl;
					}
					bestMeasure = tempMeasure;
					bestB = tempB;
				}
			}
			bCopy = b;
			swapWithSiblingChild(bCopy, i, false);
			settleBifurcations(bCopy);
			tempMeasure = tripletMeasure(bCopy, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			if(tempMeasure < initialMeasure){
				//cout << "\nbv" + makeStringTreeOneLine(bCopy);  cout.flush();
				if(history.size() > 0){
					for(unsigned int i(0); i < history.size() - 1; i++){
						if(equivalentTrees(history[i], bCopy))
							continue;
					}
				}
				history.push_back(bCopy);
				vector<Bifurcation2D> tempB = purelyDescendingSearchOLD(history, bCopy, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
				//cout << "b^"; cout.flush();
				tempMeasure = tripletMeasure(tempB, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
				if(bestMeasure > tempMeasure){
					if(print){
						cout << "\n purelyDescendingSearch(): accepting swap of " << i << " that improves bestMeasure from "
							<< bestMeasure << " to " << tempMeasure << "\n\t tree: " << makeStringTreeOneLine(tempB) << endl;
					}
					bestMeasure = tempMeasure;
					bestB = tempB;
				}
			}
		}
	}
	return bestB;
}

vector<Bifurcation2D> optimizeByPurelyDescendingSearchOLD(unsigned int numRandomSeeds, vector<Bifurcation2D> &bOrig, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient,
		unsigned int numToContinue, bool print = false, string uniquer = "", int width = 300, int height = 300, int border = 12){
	vector<Bifurcation2D> bNearest(bOrig);
	introduceBifurcationsFromSources(bNearest); // alternatively, choose a totally random configuration
	vector<vector<Bifurcation2D> > historyNearest;
	vector<Bifurcation2D> bBest(purelyDescendingSearch(historyNearest, numToContinue, bNearest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, print, uniquer, width, height, border));
	double measureBest(tripletMeasure(bNearest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	if(print)
		drawBifurcationTree(uniquer + "_bBest_from_bNearest.png", bBest, width, height, border);
	vector<int> indicesToSwap;
	for(unsigned int i(1); i < bOrig.size(); i++){
		indicesToSwap.push_back(i);
	}
	while(indicesToSwap.size() > numRandomSeeds)
		indicesToSwap.erase(indicesToSwap.begin() + (kiss()%indicesToSwap.size()));
	for(unsigned int i(1); i < bOrig.size(); i++){//skip the heart
		if(sibling_has_two_children(bNearest, i)){
			vector<Bifurcation2D> bTemp(bNearest);
			swapWithSiblingChild(bTemp, i);
			vector<vector<Bifurcation2D> > historyRand;
			bTemp = purelyDescendingSearch(historyRand, numToContinue, bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, print, uniquer, width, height, border);
			double tempMeasure(tripletMeasure(bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
			if(measureBest > tempMeasure){
				measureBest = tempMeasure;
				bBest = bTemp;
			}
			bTemp = bNearest;
			swapWithSiblingChild(bTemp, i, false);
			bTemp = purelyDescendingSearch(historyRand, numToContinue, bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, print, uniquer, width, height, border);
			tempMeasure = tripletMeasure(bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			if(measureBest > tempMeasure){
				measureBest = tempMeasure;
				bBest = bTemp;
			}
		}
	}
	/*vector<Bifurcation2D> *bSeeds = new vector<Bifurcation2D>[numRandomSeeds];
	int maxConsecutiveFailures(1000), failureCount(0), successCount(0);
	while(successCount < numRandomSeeds && failureCount < maxConsecutiveFailures){
		vector<Bifurcation2D> bRand(bOrig);
		introduceBifurcationsFromSourcesRandomly(bRand);
		bool isRedundant(false);
		for(int c(0); c < successCount && !isRedundant; c++)
			isRedundant = !isRedundant && equivalentTrees(bRand, bSeeds[c]);
		if(isRedundant)
			failureCount++;
		else{
			if(print){
				cout << "{s" << failureCount << "}" << makeStringTreeOneLine(bRand);
				drawBifurcationTree(uniquer + "_rs" + makeString(successCount) + ".png", bRand, width, height, border);
			}
			failureCount = 0;
			bSeeds[successCount] = bRand;
			vector<vector<Bifurcation2D> > historyRand;
			vector<Bifurcation2D> bRandOpt(purelyDescendingSearch(numToContinue, bRand, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, print, uniquer, width, height, border));
			double measureTemp(tripletMeasure(bRandOpt, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
			if(measureBest > measureTemp){
				if(print)
					cout << "[improving measure from " << measureBest << " to " << measureTemp << " in random]";
				measureBest = measureTemp;
				bBest = bRandOpt;
			}
			successCount++;
		}
	}
	delete[] bSeeds;*/
	return bBest;
}

bool isContainedInHistory(const vector<Bifurcation2D> &b, const vector<vector<Bifurcation2D> > &history){
	for(unsigned int i(0); i < history.size(); i++){
		if(equivalentTrees(b, history[i]))
			return true;
	}
	return false;
}

vector<Bifurcation2D> optimizeByPurelyDescendingSearch(unsigned int numRandomSeeds, vector<Bifurcation2D> &bOrig, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient,
	unsigned int numToContinue, bool print = false, string uniquer = "", int width = 300, int height = 300, int border = 12){
	vector<Bifurcation2D> bNearest(bOrig);
	introduceBifurcationsFromSources(bNearest); // alternatively, choose a totally random configuration
	vector<vector<Bifurcation2D> > history;
	vector<Bifurcation2D> bBest(purelyDescendingSearch(history, numToContinue, bNearest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, print, uniquer, width, height, border));
	double measureBest(tripletMeasure(bNearest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	if(print)
		drawBifurcationTree(uniquer + "_bBest_from_bNearest.png", bBest, width, height, border);
	vector<int> indicesToSwap;
	double measureOrig(tripletMeasure(bOrig, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	for(unsigned int s(0); s < numRandomSeeds; s++){
		bool isRedundant(true);
		vector<Bifurcation2D> bSeed;
		for(int i(0); i < 10 && isRedundant; i++){
			bSeed = bOrig;
			introduceBifurcationsFromSourcesRandomly(bSeed);
			isRedundant = !isContainedInHistory(bSeed, history);
		}
		if(!isRedundant){
			settleBifurcations(bSeed);
			bSeed = purelyDescendingSearch(history, numToContinue, bSeed, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, print, uniquer, width, height, border);
			double measureSeed(tripletMeasure(bSeed, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
			if(measureBest > measureSeed){
				measureBest = measureSeed;
				bBest = bSeed;
			}
		}
	}
	/*for(unsigned int i(1); i < bOrig.size(); i++){//skip the heart
		if(sibling_has_two_children(bNearest, i)){
			vector<Bifurcation2D> bTemp(bNearest);
			swapWithSiblingChild(bTemp, i);
			settleBifurcations(bTemp);
			if(tripletMeasure(bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) > measureOrig)
				indicesToSwap.push_back(i);
		}
	}
	while(indicesToSwap.size() > numRandomSeeds)
		indicesToSwap.erase(indicesToSwap.begin() + (kiss()%indicesToSwap.size()));
	for(unsigned int s(0); s < indicesToSwap.size(); s++){
		int i(indicesToSwap[s]);
		if(sibling_has_two_children(bNearest, i)){
			vector<Bifurcation2D> bTemp(bNearest);
			swapWithSiblingChild(bTemp, i);
			settleBifurcations(bTemp);
			vector<vector<Bifurcation2D> > historyRand;
			bTemp = purelyDescendingSearch(historyRand, numToContinue, bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, print, uniquer, width, height, border);
			double tempMeasure(tripletMeasure(bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
			if(measureBest > tempMeasure){
				measureBest = tempMeasure;
				bBest = bTemp;
			}
			bTemp = bNearest;
			swapWithSiblingChild(bTemp, i, false);
			settleBifurcations(bTemp);
			bTemp = purelyDescendingSearch(historyRand, numToContinue, bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, print, uniquer, width, height, border);
			tempMeasure = tripletMeasure(bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			if(measureBest > tempMeasure){
				measureBest = tempMeasure;
				bBest = bTemp;
			}
		}
	}*/
	return bBest;
}

vector<BranchPoint2D> optimizeByPurelyDescendingConsolidatedSearch(unsigned int numRandomSeeds, vector<Bifurcation2D> &bOrig, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient,
	unsigned int numToContinue, double settleThresh, double consolidateThresh, bool print = false, string uniquer = "", int width = 300, int height = 300, int border = 12){
	vector<Bifurcation2D> bNearest(bOrig);
	introduceBifurcationsFromSources(bNearest); // alternatively, choose a totally random configuration
	vector<vector<Bifurcation2D> > history;
	vector<BranchPoint2D> bBest(purelyDescendingConsolidatedSearch(history, numToContinue, bNearest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, print, uniquer, width, height, border));
	double measureBest(tripletMeasure(bBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	if(print)
		drawTree(uniquer + "_bBest_from_bNearest.png", bBest, width, height, border);
	vector<int> indicesToSwap;
	double measureOrig(tripletMeasure(bOrig, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	for(unsigned int s(0); s < numRandomSeeds; s++){
		bool isRedundant(true);
		vector<Bifurcation2D> bSeed;
		for(int i(0); i < 10 && isRedundant; i++){
			bSeed = bOrig;
			introduceBifurcationsFromSourcesRandomly(bSeed);
			isRedundant = !isContainedInHistory(bSeed, history);
		}
		if(!isRedundant){
			settleBifurcations(bSeed);
			vector<BranchPoint2D> bSeedc = purelyDescendingConsolidatedSearch(history, numToContinue, bSeed, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, print, uniquer, width, height, border);
			double measureSeed(tripletMeasure(bSeedc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
			if(measureBest > measureSeed){
				//cout << "optimizeByPurelyDescendingConsolidatedSearch(): accepting " << measureSeed << " to replace " << measureBest << endl;
				measureBest = measureSeed;
				bBest = bSeedc;
			}
		}
	}
	/*for(unsigned int i(1); i < bOrig.size(); i++){//skip the heart
	if(sibling_has_two_children(bNearest, i)){
	vector<Bifurcation2D> bTemp(bNearest);
	swapWithSiblingChild(bTemp, i);
	settleBifurcations(bTemp);
	if(tripletMeasure(bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) > measureOrig)
	indicesToSwap.push_back(i);
	}
	}
	while(indicesToSwap.size() > numRandomSeeds)
	indicesToSwap.erase(indicesToSwap.begin() + (kiss()%indicesToSwap.size()));
	for(unsigned int s(0); s < indicesToSwap.size(); s++){
	int i(indicesToSwap[s]);
	if(sibling_has_two_children(bNearest, i)){
	vector<Bifurcation2D> bTemp(bNearest);
	swapWithSiblingChild(bTemp, i);
	settleBifurcations(bTemp);
	vector<vector<Bifurcation2D> > historyRand;
	bTemp = purelyDescendingSearch(historyRand, numToContinue, bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, print, uniquer, width, height, border);
	double tempMeasure(tripletMeasure(bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	if(measureBest > tempMeasure){
	measureBest = tempMeasure;
	bBest = bTemp;
	}
	bTemp = bNearest;
	swapWithSiblingChild(bTemp, i, false);
	settleBifurcations(bTemp);
	bTemp = purelyDescendingSearch(historyRand, numToContinue, bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, print, uniquer, width, height, border);
	tempMeasure = tripletMeasure(bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	if(measureBest > tempMeasure){
	measureBest = tempMeasure;
	bBest = bTemp;
	}
	}
	}*/
	return bBest;
}

void exhaustiveGreed(){
	//kisset(30677, 16433, 20738);
	//kisset(123960, 29603, 31767);
	// 1 2 3  4  5   6    7      8      9       10
	// 1 1 3 15 105 945 10395 135135 2027025 34459425
	int maxFailedInsertions(100000), numBest(100), maxFailedNetworks(10000), targetServiceVolumes(8), numRandSeeds(2),
		tempIncs(2), numTopKeep(2), numIndAnnealing(2), annealingIts(2), numToContinue(3);
	double nearestRadius(1.0), initTemp(1.0);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0), lengthThresh(1.0e-6);
	double totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0);
	string uniquer("outputs/exhaustiveGreed_" + makeString(targetServiceVolumes) + "_");
	int width(300), height(300), border(12);//, numBins(100), numToDraw(5);

	cout << "\n rescaledDerivatives()\n uniquer = " << uniquer
		<< "\n maxFailedNetworks = " << maxFailedNetworks << "\n numRandSeeds = " << numRandSeeds
		<< "\n targetServiceVolumes = " << targetServiceVolumes
		<< "\n lengthThresh = " << lengthThresh
		<< "\n sizeX = " << sizeX << "\n sizeY = " << sizeY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n numBest = " << numBest //<< "\n numBins = " << numBins
		<< "\n initTemp = " << initTemp << "\n tempIncs = " << tempIncs
		<< "\n numTopKeep = " << numTopKeep << "\n numIndAnnealing = " << numIndAnnealing
		<< "\n annealingIts = " << annealingIts
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< endl;

	int failedNetworks(0);
	while(failedNetworks < maxFailedNetworks){
		vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty));
		if(b.size() == targetServiceVolumes + 1){
			cout << "\n{c" << failedNetworks << "}" << endl;
			time_t coolStart = time(NULL);
			vector<Bifurcation2D> seedCool(b);
			introduceBifurcationsFromSources(seedCool);
			settleBifurcations(seedCool);
			vector<Bifurcation2D> bCool(simulatedCooling(seedCool, initTemp, tempIncs, numTopKeep, numIndAnnealing, annealingIts, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
			drawBifurcationTree(uniquer + "_cool.png", bCool, width, height, border);
			cout << "\n Cool search found optimal configuration in " <<  niceTime(1.0*(time(NULL) - coolStart)) << ":\n\t" + makeStringTreeOneLine(bCool) + "\t\t" + makeString(tripletMeasure(bCool, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)) << endl;
			time_t greedyStart(time(NULL));
			vector<Bifurcation2D> bGreedy(optimizeByPurelyDescendingSearch(numRandSeeds, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, numToContinue, true, uniquer + "pDS"));
			drawBifurcationTree(uniquer + "_greedy.png", bGreedy, width, height, border);
			cout << "\n Greedy search found optimal configuration in " <<  niceTime(1.0*(time(NULL) - greedyStart)) << ":\n\t" + makeStringTreeOneLine(bGreedy) + "\t\t" + makeString(tripletMeasure(bGreedy, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)) << endl;
			if(targetServiceVolumes < 9){
				time_t startExhaustive(time(NULL));
				vector<Bifurcation2D> *bestBs = new vector<Bifurcation2D>[numBest];
				double *bestBMeasures = new double[numBest];
				exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
				cout << "\n Exhaustive search took " << niceTime(1.0*(time(NULL) - startExhaustive)) << " to complete." << endl;
				drawBifurcationTree(uniquer + "_exhaustive.png", bestBs[0], width, height, border);
				cout << "\n Exhaustive search found optimal configuration:\n\t" + makeStringTreeOneLine(bestBs[0]) + "\t\t" + makeString(tripletMeasure(bestBs[0], totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)) << endl;
				cout << "\n Exhaustive search optimal configuration is ";
				if(equivalentTrees(bestBs[0], bCool))
					cout << "INDEED";
				else
					cout << "NOT";
				cout << " the same as the cool search optimal tree." << endl;
				cout << "\n Exhaustive search optimal configuration is ";
				if(equivalentTrees(bestBs[0], bGreedy))
					cout << "INDEED";
				else
					cout << "NOT";
				cout << " the same as the greedy search optimal tree." << endl;
				int matchConfigCool(-1);
				for(int i(0); i < numBest; i++)
					if(equivalentTrees(bestBs[i], bCool))
						matchConfigCool = i;
				if(matchConfigCool > -1)
					cout << "\n bCool matches the configuration with rank " << matchConfigCool + 1 << " out of " << numBest << " (1 is best)." << endl;
				else
					cout << "\n No configuration matching bCool was found in the best " << numBest << " configurations." << endl;
				int matchConfigGreedy(-1);
				for(int i(0); i < numBest; i++)
					if(equivalentTrees(bestBs[i], bGreedy))
						matchConfigGreedy = i;
				if(matchConfigGreedy > -1)
					cout << "\n bGreedy matches the configuration with rank " << matchConfigGreedy + 1 << " out of " << numBest << " (1 is best)." << endl;
				else
					cout << "\n No configuration matching bGreedy was found in the best " << numBest << " configurations." << endl;

				delete[] bestBs;
				delete[] bestBMeasures;
			}else
				cout << "\n No exhaustive hierarchy search performed." << endl;
			double measureGreedy(tripletMeasure(bGreedy, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
			double measureCool(tripletMeasure(bCool, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
			if(equivalentTrees(bGreedy, bCool))
				cout << "\nbGreedy and bCool are THE SAME." << endl;
			else if(measureGreedy < measureCool)
				cout << "\nmeasureGreedy < measureCool\n\t" << measureGreedy << " < " << measureCool << endl;
			else
				cout << "\nmeasureGreedy > measureCool\n\t" << measureGreedy << ">" << measureCool << endl;
			failedNetworks = maxFailedNetworks;
		}else
			failedNetworks++;
	}
			
}

vector<Bifurcation2D> initializeEvenRandomSpacingServiceVolumes(int maxFailedInsertions, double sizeX, double sizeY, double nearestRadius, double heartx, double hearty, int targetServiceVolumes, int maxFailedNetworks){
	for(int i(0); i < maxFailedNetworks; i++){
		vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty));
		if(b.size() - 1 == targetServiceVolumes)
			return b;
	}
	return vector<Bifurcation2D>();
}

void greedySuccessRate(){
	vector<int> targetServiceVolumes, numToContinue, numRandSeeds; // all quantities should be ascending
	for(int i(0); i < 3; i++) targetServiceVolumes.push_back(i + 4);
	for(int i(0); i < 3; i++) numToContinue.push_back(i + 1);
	for(int i(0); i < 2; i++) numRandSeeds.push_back(i);
	int numNetworks(10);
	int maxFailedInsertions(100000), numBest(200), maxFailedNetworks(10000);
	double nearestRadius(1.0);
	double totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0);
	string uniquer("outputs/greedySuccessRate4-7someRand_c010_");
	int width(300), height(300), border(12);

	string outFileFn(uniquer + "summary.txt");
	ofstream outFile(outFileFn.c_str());
	outFile << "\n rescaledDerivatives()\n uniquer = " << uniquer
		<< "\n maxFailedNetworks = " << maxFailedNetworks
		<< "\n numNetworks = " << numNetworks
		<< "\n numRandSeeds = " << makeVectorString(numRandSeeds)
		<< "\n numToContinue = " << makeVectorString(numToContinue)
		<< "\n targetServiceVolumes = " << makeVectorString(targetServiceVolumes)
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n numBest = " << numBest
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< endl;

	int ****greedyRanks = new int***[targetServiceVolumes.size()];
	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		greedyRanks[tsvi] = new int**[numNetworks];
		double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes[tsvi])*acos(-1.0))/2.0), sizeY(sizeX);
		double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
		cout << "\n targetServiceVolumes = " << targetServiceVolumes[tsvi]
			<< " \t sizeX = sizeY = " << sizeX
			<< " \t heartx = " << heartx << " \t hearty = " << hearty << endl;
		for(int n(0); n < numNetworks; n++){
			greedyRanks[tsvi][n] = new int*[numToContinue.size()];
			vector<Bifurcation2D> b(initializeEvenRandomSpacingServiceVolumes(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty, targetServiceVolumes[tsvi], maxFailedNetworks));
			vector<Bifurcation2D> *bestBs = new vector<Bifurcation2D>[numBest];
			double *bestBMeasures = new double[numBest];
			exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
				greedyRanks[tsvi][n][ntci] = new int[numRandSeeds.size()];
				cout << "\t numToContinue = " << numToContinue[ntci] << " \t (network " << n + 1 << " of " << numNetworks << " for " << targetServiceVolumes[tsvi] << " capillaries.)" << endl;
				for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++){
					cout << "\t\t numRandSeeds = " << numRandSeeds[nrsi] << endl;
					vector<Bifurcation2D> bGreedy(optimizeByPurelyDescendingSearch(numRandSeeds[nrsi], b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, numToContinue[ntci]));
					double greedyMeasure(tripletMeasure(bGreedy, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
					if(targetServiceVolumes[tsvi] < 9){
						time_t startExhaustive(time(NULL));
						
						//cout << "\n Exhaustive search took " << niceTime(1.0*(time(NULL) - startExhaustive)) << " to complete." << endl;
						//cout << "\n Exhaustive search optimal configuration is ";
						//if(equivalentTrees(bestBs[0], bGreedy))
						//	cout << "INDEED";
						//else
						//	cout << "NOT";
						//cout << " the same as the greedy search optimal tree." << endl;
						int matchConfigGreedy(-1);
						for(int i(0); i < numBest && matchConfigGreedy < 0; i++){
							if(equivalentTrees(bestBs[i], bGreedy))
								matchConfigGreedy = i;
						}
						greedyRanks[tsvi][n][ntci][nrsi] = matchConfigGreedy + 1;
						if(matchConfigGreedy > -1){
							cout << "\t\t\trank " << matchConfigGreedy + 1 << " \t exM[" << matchConfigGreedy << "]:" << bestBMeasures[matchConfigGreedy]
								<< " \t gM:" << greedyMeasure << " \t exM[0]:" << bestBMeasures[0] << endl;
						}else
							cout << "\t\t\trank -X-" << endl;
					}else
						cout << "\t\t\t (No exhaustive hierarchy)" << endl;
				}
			}
			delete[] bestBs;
			delete[] bestBMeasures;
		}
	}

	outFile << "\n\n v-targetServiceVolumes/numToContinue(mostRand)->";
	for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++)
		outFile << "\t" << numToContinue[ntci] << "\t +/- \t %-first";
	outFile << endl;
	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		outFile << targetServiceVolumes[tsvi];
		for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
			vector<int> tempRanks;
			int firstsCount(0);
			for(int n(0); n < numNetworks; n++){
				int gRank(greedyRanks[tsvi][n][ntci][numRandSeeds.size() - 1]);
				if(gRank > -1){
					tempRanks.push_back(gRank);
					if(gRank == 1) firstsCount++;
				}
			}
			outFile << "\t" << mean(tempRanks) << "\t" << stDev(tempRanks) << "\t" << 100.0*double(firstsCount)/double(numNetworks) << "%";
		}
		outFile << endl;
	}

	outFile << "\n\n v-targetServiceVolumes/numRandSeeds(mostToCont)->";
	for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++)
		outFile << "\t" << numRandSeeds[nrsi] << "\t +/- \t %-first";
	outFile << endl;
	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		outFile << targetServiceVolumes[tsvi];
		for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++){
			vector<int> tempRanks;
			int firstsCount(0);
			for(int n(0); n < numNetworks; n++){
				int gRank(greedyRanks[tsvi][n][numToContinue.size() - 1][nrsi]);
				tempRanks.push_back(gRank);
				if(gRank == 1) firstsCount++;
			}
			outFile << "\t" << mean(tempRanks) << "\t" << stDev(tempRanks) << "\t" << 100.0*double(firstsCount)/double(numNetworks) << "%";
		}
		outFile << endl;
	}

	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		outFile << "\n\n targetServiceVolumes = " << targetServiceVolumes[tsvi] << endl
			<< "v-numToContinue/numRandSeeds->";
		for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++)
			outFile << "\t" << numRandSeeds[nrsi] << "\t +/- \t %-first";
		outFile << endl;
		for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
			outFile << numToContinue[ntci];
			for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++){
				vector<int> tempRanks;
				int firstsCount(0);
				for(int n(0); n < numNetworks; n++){
					int gRank(greedyRanks[tsvi][n][ntci][nrsi]);
					tempRanks.push_back(gRank);
					if(gRank == 1) firstsCount++;
				}
				outFile << "\t" << mean(tempRanks) << "\t" << stDev(tempRanks) << "\t" << 100.0*double(firstsCount)/double(numNetworks) << "%";
			}
			outFile << endl;
		}
	}


	outFile.close();
	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		for(int n(0); n < numNetworks; n++){
			for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++)
				delete[] greedyRanks[tsvi][n][ntci];
			delete[] greedyRanks[tsvi][n];
		} 
		delete[] greedyRanks[tsvi];
	}
	delete[] greedyRanks;
}

bool subtreeOfContains(const vector<Bifurcation2D> &b, int x, int y){
	//queue<int> toCheck;
	//toCheck.push(x);
	//while(toCheck.size() > 0){
	//	int z(toCheck.front());
	//	toCheck.pop();
	//	if(z == y)
	//		return true;
	//	for(unsigned int c(0); c < b[z].childIndex.size(); c++)
	//		toCheck.push(b[z].childIndex[c]);
	//}
	if(x == y) return true;
	int i(y);
	while(b[i].parentIndex > -1){
		if(b[i].parentIndex == x)
			return true;
		i = b[i].parentIndex;
	}
	return false;
}

bool subtreeOfContains(const vector<Bifurcation3D> &b, int x, int y){
	//queue<int> toCheck;
	//toCheck.push(x);
	//while(toCheck.size() > 0){
	//	int z(toCheck.front());
	//	toCheck.pop();
	//	if(z == y)
	//		return true;
	//	for(unsigned int c(0); c < b[z].childIndex.size(); c++)
	//		toCheck.push(b[z].childIndex[c]);
	//}
	if(x == y) return true;
	int i(y);
	while(b[i].parentIndex > -1){
		if(b[i].parentIndex == x)
			return true;
		i = b[i].parentIndex;
	}
	return false;
}

bool intersectingSubtrees(const vector<Bifurcation2D> &b, int x, int y){
	if(subtreeOfContains(b, x, y)) return true;
	if(subtreeOfContains(b, y, x)) return true;
	return false;
}

bool intersectingSubtrees(const vector<Bifurcation3D> &b, int x, int y){
	if(subtreeOfContains(b, x, y)) return true;
	if(subtreeOfContains(b, y, x)) return true;
	return false;
}

void removeFromVector(int i, vector<int> &v){
	for(unsigned int x(0); x < v.size(); x++){
		if(v[x] == i){
			v.erase(v.begin() + x);
			return;
		}
	}
}

bool isValidSingleRegraft(vector<Bifurcation2D> &b, int p, int c, bool keepFirstChild = true){
	// don't regraft heart     don't regraft tips/heart       don't parent heart      regrafting will do nothing
	if(b[p].parentIndex < 0 || b[p].childIndex.size() != 2 || b[c].parentIndex < 0 || containedIn(c, b[p].childIndex)
			|| p >= (int)b.size() || c >= (int)b.size() || p == c){
			// bad index        bad index                  nonsense
		return false;
	}
	//cout << "\n isValidSingleRegraft(): not nonsense" << endl;
	// no reticulations:
	if(keepFirstChild){
		//cout << "\n isValidSingleRegraft(): checking child 0" << endl;
		if(subtreeOfContains(b, b[p].childIndex[0], c))//intersectingSubtrees(b, c, b[p].childIndex[0]))
			return false;
	}else{
		//cout << "\n isValidSingleRegraft(): checking child 1" << endl;
		if(subtreeOfContains(b, b[p].childIndex[1], c))//intersectingSubtrees(b, c, b[p].childIndex[1]))
			return false;
	}
	return true;
}

bool isValidSingleRegraft(vector<Bifurcation3D> &b, int p, int c, bool keepFirstChild = true){
	// don't regraft heart     don't regraft tips/heart       don't parent heart      regrafting will do nothing
	if(b[p].parentIndex < 0 || b[p].childIndex.size() != 2 || b[c].parentIndex < 0 || containedIn(c, b[p].childIndex)
			|| p >= (int)b.size() || c >= (int)b.size() || p == c){
			// bad index        bad index                  nonsense
		return false;
	}
	// no reticulations:
	if(keepFirstChild){
		if(subtreeOfContains(b, b[p].childIndex[0], c))//intersectingSubtrees(b, c, b[p].childIndex[0]))
			return false;
	}else{
		if(subtreeOfContains(b, b[p].childIndex[1], c))//intersectingSubtrees(b, c, b[p].childIndex[1]))
			return false;
	}
	return true;
}

//moves p to be the parent of c; returns false if regrafting is impossible
//if keepFirstChild is true, then p takes child[0] with it
bool regraftSingle(vector<Bifurcation2D> &b, int p, int c, bool keepFirstChild = true){
	if(!isValidSingleRegraft(b, p, c, keepFirstChild))
		return false;
	if(keepFirstChild){
		b[b[p].parentIndex].childIndex.push_back(b[p].childIndex[1]);
		b[b[p].childIndex[1]].parentIndex = b[p].parentIndex;
		b[p].childIndex.erase(b[p].childIndex.begin() + 1);
	}else{
		b[b[p].parentIndex].childIndex.push_back(b[p].childIndex[0]);
		b[b[p].childIndex[0]].parentIndex = b[p].parentIndex;
		b[p].childIndex.erase(b[p].childIndex.begin());
	}
	removeFromVector(p, b[b[p].parentIndex].childIndex);
	removeFromVector(c, b[b[c].parentIndex].childIndex);
	b[b[c].parentIndex].childIndex.push_back(p);
	b[p].parentIndex = b[c].parentIndex;
	b[c].parentIndex = p;
	b[p].childIndex.push_back(c);
	return true;
}

bool regraftSingle(vector<Bifurcation3D> &b, int p, int c, bool keepFirstChild = true){
	if(!isValidSingleRegraft(b, p, c, keepFirstChild))
		return false;
	if(keepFirstChild){
		b[b[p].parentIndex].childIndex.push_back(b[p].childIndex[1]);
		b[b[p].childIndex[1]].parentIndex = b[p].parentIndex;
		b[p].childIndex.erase(b[p].childIndex.begin() + 1);
	}else{
		b[b[p].parentIndex].childIndex.push_back(b[p].childIndex[0]);
		b[b[p].childIndex[0]].parentIndex = b[p].parentIndex;
		b[p].childIndex.erase(b[p].childIndex.begin());
	}
	removeFromVector(p, b[b[p].parentIndex].childIndex);
	removeFromVector(c, b[b[c].parentIndex].childIndex);
	b[b[c].parentIndex].childIndex.push_back(p);
	b[p].parentIndex = b[c].parentIndex;
	b[c].parentIndex = p;
	b[p].childIndex.push_back(c);
	return true;
}

//swaps x and y
void regraft(vector<Bifurcation2D> &b, int x, int y){
	if(b[x].parentIndex < 0 || b[y].parentIndex < 0 || b[x].parentIndex == b[y].parentIndex)
		return;
	for(unsigned int c(0); c < b[b[x].parentIndex].childIndex.size(); c++){
		if(b[b[x].parentIndex].childIndex[c] == x){
			b[b[x].parentIndex].childIndex[c] = y;
			c = b[b[x].parentIndex].childIndex.size();
		}
	}
	for(unsigned int c(0); c < b[b[y].parentIndex].childIndex.size(); c++){
		if(b[b[y].parentIndex].childIndex[c] == y){
			b[b[y].parentIndex].childIndex[c] = x;
			c = b[b[y].parentIndex].childIndex.size();
		}
	}
	int foster(b[x].parentIndex);
	b[x].parentIndex = b[y].parentIndex;
	b[y].parentIndex = foster;
}

//swaps x and y
void regraft(vector<Bifurcation3D> &b, int x, int y){
	if(b[x].parentIndex < 0 || b[y].parentIndex < 0 || b[x].parentIndex == b[y].parentIndex)
		return;
	for(unsigned int c(0); c < b[b[x].parentIndex].childIndex.size(); c++){
		if(b[b[x].parentIndex].childIndex[c] == x){
			b[b[x].parentIndex].childIndex[c] = y;
			c = b[b[x].parentIndex].childIndex.size();
		}
	}
	for(unsigned int c(0); c < b[b[y].parentIndex].childIndex.size(); c++){
		if(b[b[y].parentIndex].childIndex[c] == y){
			b[b[y].parentIndex].childIndex[c] = x;
			c = b[b[y].parentIndex].childIndex.size();
		}
	}
	int foster(b[x].parentIndex);
	b[x].parentIndex = b[y].parentIndex;
	b[y].parentIndex = foster;
}

// assumes the heart is 0
vector<BranchPoint2D> optimizeByGlobalConsolidatedGrafting(const vector<Bifurcation2D> &bOrig, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient, double settleThresh, double consolidateThresh, int maxPasses = 10000, int settleIts = 10000){
	vector<Bifurcation2D> b(bOrig);
	introduceBifurcationsFromSources(b);
	settleBifurcations(b, false, settleThresh, settleIts);
	vector<BranchPoint2D> bBest(consolidateDegenerates(b, consolidateThresh));
	double mBest(tripletMeasure(bBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	int swapCount(1), passCount(0);
	vector<int> swappable;
	for(unsigned int i(1); i < b.size(); i++) // not the heart
		swappable.push_back(i);
	vector<Pair<int> > lastPass(1, Pair<int>(-1, -1));
	while(lastPass.size() > 0 && passCount < maxPasses){
		lastPass.clear();
		//cout << "\n passCount = " << passCount << endl;
		swapCount = 0;
		vector<int> toCheck(swappable);
		while(toCheck.size() > 0){
			int checkingIndex(kiss()%toCheck.size());
			int checking(toCheck[checkingIndex]);
			for(unsigned int p(0); p < swappable.size(); p++){
				if(swappable[p] == checking || b[swappable[p]].parentIndex == b[checking].parentIndex)
					continue;
				if(intersectingSubtrees(b, swappable[p], checking))
					continue;
				//cout << "(ps" << passCount << "c" << checking << "p" << p << ") "; cout.flush();
				vector<Bifurcation2D> bCopy(b);
				regraft(bCopy, checking, swappable[p]);
				bool reportSettle(false && passCount == 2 && checking == 16 && p == 8);
				if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy, reportSettle, settleThresh, settleIts)){
					vector<BranchPoint2D> bc(consolidateDegenerates(bCopy, consolidateThresh));
					double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
					if(mBest > m){
						//cout << "\n " << passCount << " accepting regraft of " << checking << " with " << swappable[p]
						//	<< "\n\t to get from measure " << mBest << " to measure " << m << endl;
						swapCount++;
						mBest = m;
						bBest = bc;
						b = bCopy;
						lastPass.push_back(Pair<int>(checking, swappable[p]));
					}//else
						//regraft(b, checking, swappable[p]);
				}
			}
			toCheck.erase(toCheck.begin() + checkingIndex);
		}
		passCount++;
		for(unsigned int i(0); i < lastPass.size(); i++){
			for(unsigned int j(i + 1); j < lastPass.size(); j++){
				if(lastPass[i].sameUnordered(lastPass[j])){
					lastPass.erase(lastPass.begin() + j);
					lastPass.erase(lastPass.begin() + i);
					i--;
					j = lastPass.size();
				}
			}
		}
	}
	/*swapCount = 1;
	while(swapCount > 0 && passCount < maxPasses){
		swapCount = 0;
		vector<int> toCheck(swappable);
		while(toCheck.size() > 0){
			int checkingIndex(kiss()%toCheck.size());
			int checking(toCheck[checkingIndex]);
			for(int lev(1); lev < 10; lev++){ // not the heart
				int c(checking);
				for(int ell(0); ell < lev; ell++){
					if(b[c].parentIndex > -1)
						c = b[c].parentIndex;
				}
				if(c > -1 && b[c].parentIndex == checking || intersectingSubtrees(b, c, checking))
					continue;
				//cout << "(ps" << passCount << "c" << checking << "p" << p << ") "; cout.flush();
				vector<Bifurcation2D> bCopy(b), bCopy2(b);
				regraftSingle(bCopy, checking, c);
				regraftSingle(bCopy2, checking, c, false);
				bool reportSettle(false && passCount == 2 && checking == 16 && c == 8);
				if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy, reportSettle, settleThresh, settleIts)){
					vector<BranchPoint2D> bc(consolidateDegenerates(bCopy, consolidateThresh));
					double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
					if(mBest > m){
						swapCount++;
						mBest = m;
						bBest = bc;
						b = bCopy;
					}
				}
				if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy2, reportSettle, settleThresh, settleIts)){
					vector<BranchPoint2D> bc(consolidateDegenerates(bCopy2, consolidateThresh));
					double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
					if(mBest > m){
						swapCount++;
						mBest = m;
						bBest = bc;
						b = bCopy2;
					}
				}
			}
			toCheck.erase(toCheck.begin() + checkingIndex);
		}
		passCount++;
	}*/
	if(passCount == maxPasses){
		cout << "\n optimizeByGlobalConsolidatedGrafting(): after " << passCount << " passes, still performed " << swapCount << " swaps." << endl;
	}
	return bBest;
}

int globalTemp(0);

// assumes the heart is 0
vector<BranchPoint2D> optimizeByGlobalConsolidatedGraftingNEW(const vector<Bifurcation2D> &bOrig, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient, double settleThresh, double consolidateThresh, double sepThresh, int maxPasses = 10000, int settleIts = 10000){
	vector<Bifurcation2D> b(bOrig);
	introduceBifurcationsFromSources(b);
	settleBifurcations(b, false, settleThresh, settleIts);
	vector<BranchPoint2D> bBest(consolidateDegenerates(b, consolidateThresh));
	double mBest(tripletMeasure(bBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	int swapCount(1), passCount(0);
	vector<int> swappable;
	for(unsigned int i(1); i < b.size(); i++) // not the heart
		swappable.push_back(i);
	vector<Pair<int> > lastPass(1, Pair<int>(-1, -1));
	while(lastPass.size() > 0 && passCount < maxPasses){
		lastPass.clear();
		//cout << "\n passCount = " << passCount << endl;
		swapCount = 0;
		vector<int> toCheck(swappable);
		while(toCheck.size() > 0){
			int checkingIndex(kiss()%toCheck.size());
			int checking(toCheck[checkingIndex]);
			for(unsigned int p(0); p < swappable.size(); p++){
				if(swappable[p] == checking || b[swappable[p]].parentIndex == b[checking].parentIndex)
					continue;
				if(intersectingSubtrees(b, swappable[p], checking))
					continue;
				//cout << "(ps" << passCount << "c" << checking << "p" << p << ") "; cout.flush();
				vector<Bifurcation2D> bCopy(b);
				regraft(bCopy, checking, swappable[p]);
				bool reportSettle(false && passCount == 2 && checking == 16 && p == 8);
				if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy, reportSettle, settleThresh, settleIts)){
					vector<BranchPoint2D> bc(consolidateDegenerates(bCopy, consolidateThresh));
					double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
					if(mBest > m){
						//cout << "\n " << passCount << " accepting regraft of " << checking << " with " << swappable[p]
						//	<< "\n\t to get from measure " << mBest << " to measure " << m << endl;
						swapCount++;
						mBest = m;
						bBest = bc;
						b = bCopy;
						lastPass.push_back(Pair<int>(checking, swappable[p]));
					}
				}
			}
			toCheck.erase(toCheck.begin() + checkingIndex);
		}
		passCount++;
		for(unsigned int i(0); i < lastPass.size(); i++){
			for(unsigned int j(i + 1); j < lastPass.size(); j++){
				if(lastPass[i].sameUnordered(lastPass[j])){
					lastPass.erase(lastPass.begin() + j);
					lastPass.erase(lastPass.begin() + i);
					i--;
					j = lastPass.size();
				}
			}
		}
	}
	swapCount = 1;
	double **seps = new double*[b.size()];
	for(unsigned int i(0); i < b.size(); i++)
		seps[i] = new double[b.size()];
	while(swapCount > 0 && passCount < maxPasses){
		for(unsigned int i(1); i < swappable.size(); i++){ // skipping the heart
			for(unsigned int j(i + 1); j < b.size(); j++)
				seps[i][j] = seps[j][i] = separation2D(b[i], b[j]);
		}
		swapCount = 0;
		vector<int> toCheck(swappable);
		while(toCheck.size() > 0){
			int checkingIndex(kiss()%toCheck.size());
			int checking(toCheck[checkingIndex]);
			toCheck.erase(toCheck.begin() + checkingIndex);
			vector<int> toCheckWith(toCheck);
			while(toCheckWith.size() > 0){ // not the heart
				int cIndex(kiss()%toCheckWith.size());
				int c(toCheckWith[cIndex]);
				toCheckWith.erase(toCheckWith.begin() + cIndex);
				if(seps[checking][c] < sepThresh){
					//cout << "(ps" << passCount << "ch" << checking << "c" << c << ") "; cout.flush();
					vector<Bifurcation2D> bCopy(b), bCopy2(b);
					bool reportSettle(false && passCount == 2 && checking == 15 && c == 16);
					if(regraftSingle(bCopy, checking, c)){
						//drawBifurcationTree("outputs/oBGCGNew_1pre" + makeString(checking) + "k" + makeString(c) + "_" + makeString(globalTemp) + ".png", bCopy, 300, 300, 10);
						if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy, reportSettle, settleThresh, settleIts)){
							vector<BranchPoint2D> bc(consolidateDegenerates(bCopy, consolidateThresh));
							double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
							//if(globalTemp < 100){
							//	string sss("");
							//	if(mBest > m)
							//		sss = "*";
							//	drawTree("outputs/oBGCGNew" + makeString(checking) + "k" + makeString(c) + "_" + makeString(globalTemp) + sss + ".png", bc, 300, 300, 10);
							//	globalTemp++;
							//}
							if(mBest > m){
								swapCount++;
								mBest = m;
								bBest = bc;
								b = bCopy;
							}
						}
					}
					if(regraftSingle(bCopy2, checking, c, false)){
						//drawBifurcationTree("outputs/oBGCGNew_2pre" + makeString(checking) + "k" + makeString(c) + "_" + makeString(globalTemp) + ".png", bCopy2, 300, 300, 10);
						if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy2, reportSettle, settleThresh, settleIts)){
							vector<BranchPoint2D> bc(consolidateDegenerates(bCopy2, consolidateThresh));
							double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
							if(mBest > m){
								swapCount++;
								mBest = m;
								bBest = bc;
								b = bCopy2;
							}
						}
					}
				}
			}
		}
		passCount++;
	}

	delete[] seps;
	if(passCount == maxPasses){
		cout << "\n optimizeByGlobalConsolidatedGrafting(): after " << passCount << " passes, still performed " << swapCount << " swaps." << endl;
	}
	return bBest;
}

double minimumTipRatio(const vector<Bifurcation2D> &b, int *numT, int startIndex = 0){
	double minTipRatio(-1.0);
	//cout << "minimumTipRatio(): " << b[startIndex].childIndex.size() << " children" << endl;
	if(b[startIndex].childIndex.size() < 1)
		return -1.0;
	for(unsigned int c(0); c < b[startIndex].childIndex.size(); c++){
		double ratioUnderChild(minimumTipRatio(b, numT, b[startIndex].childIndex[c]));
		if(ratioUnderChild >= 0.0 && (minTipRatio < 0.0 || minTipRatio > ratioUnderChild))
			minTipRatio = ratioUnderChild;
		for(unsigned int cc(c + 1); cc < b[startIndex].childIndex.size(); cc++){
			double tipRatio(double(numT[b[startIndex].childIndex[c]])/double(numT[b[startIndex].childIndex[cc]]));
			if(tipRatio > 1.0)
				tipRatio = 1.0/tipRatio;
			if(minTipRatio < 0.0 || minTipRatio > tipRatio)
				minTipRatio = tipRatio;
		}
	}
	return minTipRatio;
}

double minimumTipRatio(const vector<Bifurcation2D> &b){
	int *numT = numTips(b);
	double mTR(minimumTipRatio(b, numT));
	delete[] numT;
	return mTR;
}

double adjustedMinimumTipRatio(const vector<Bifurcation2D> &b, int *numT, int startIndex = 0){
	double minTipRatio(-1.0);
	//cout << "minimumTipRatio(): " << b[startIndex].childIndex.size() << " children" << endl;
	if(b[startIndex].childIndex.size() < 1)
		return -1.0;
	for(unsigned int c(0); c < b[startIndex].childIndex.size(); c++){
		double ratioUnderChild(adjustedMinimumTipRatio(b, numT, b[startIndex].childIndex[c]));
		if(ratioUnderChild >= 0.0 && (minTipRatio < 0.0 || minTipRatio > ratioUnderChild))
			minTipRatio = ratioUnderChild;
		for(unsigned int cc(c + 1); cc < b[startIndex].childIndex.size(); cc++){
			double tipRatio(1.0);
			if(numT[b[startIndex].childIndex[c]] < numT[b[startIndex].childIndex[cc]])
				tipRatio = double(numT[b[startIndex].childIndex[c]] + 1)/double(numT[b[startIndex].childIndex[cc]]);
			else if(numT[b[startIndex].childIndex[c]] > numT[b[startIndex].childIndex[cc]])
				tipRatio = double(numT[b[startIndex].childIndex[cc]] + 1)/double(numT[b[startIndex].childIndex[c]]);
			if(tipRatio > 1.0) // should never happen
				tipRatio = 1.0/tipRatio;
			if(minTipRatio < 0.0 || minTipRatio > tipRatio)
				minTipRatio = tipRatio;
		}
	}
	return minTipRatio;
}

// similar to minimumTipRatio, but adds 1 to the smaller branch if asymmetric
double adjustedMinimumTipRatio(const vector<Bifurcation2D> &b){
	int *numT = numTips(b);
	double mTR(adjustedMinimumTipRatio(b, numT));
	delete[] numT;
	return mTR;
}

double minimumTipRatio(const vector<BranchPoint2D> &b, int *numT, int startIndex = 0){
	double minTipRatio(-1.0);
	//cout << "minimumTipRatio(): " << b[startIndex].childIndex.size() << " children" << endl;
	if(b[startIndex].childIndex.size() < 1)
		return -1.0;
	for(unsigned int c(0); c < b[startIndex].childIndex.size(); c++){
		double ratioUnderChild(minimumTipRatio(b, numT, b[startIndex].childIndex[c]));
		if(ratioUnderChild >= 0.0 && (minTipRatio < 0.0 || minTipRatio > ratioUnderChild))
			minTipRatio = ratioUnderChild;
		for(unsigned int cc(c + 1); cc < b[startIndex].childIndex.size(); cc++){
			double tipRatio(double(numT[b[startIndex].childIndex[c]])/double(numT[b[startIndex].childIndex[cc]]));
			if(tipRatio > 1.0)
				tipRatio = 1.0/tipRatio;
			if(minTipRatio < 0.0 || minTipRatio > tipRatio)
				minTipRatio = tipRatio;
		}
	}
	return minTipRatio;
}

double minimumTipRatio(const vector<BranchPoint2D> &b){
	int *numT = numTips(b);
	return minimumTipRatio(b, numT);
}

// assumes the heart is 0
vector<BranchPoint2D> optimizeByGlobalConsolidatedGraftingTipRatioThresh(const vector<Bifurcation2D> &bOrig, vector<Bifurcation2D> &bBi, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient, double settleThresh, double consolidateThresh, double sepThresh, double tipRatioThresh, int maxPasses = 10000, int settleIts = 10000){
	vector<Bifurcation2D> b(bOrig);
	//introduceBifurcationsFromSources(b);
	introduceBifurcationsSymmetricallyFromHeart(b, settleThresh);
	settleBifurcations(b, false, settleThresh, settleIts);
	bBi = b;
	vector<BranchPoint2D> bBest(consolidateDegenerates(b, consolidateThresh));
	double mBest(tripletMeasure(bBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	int swapCount(1), passCount(0);
	vector<int> swappable;
	for(unsigned int i(1); i < b.size(); i++) // not the heart
		swappable.push_back(i);
	vector<Pair<int> > lastPass(1, Pair<int>(-1, -1));
	while(lastPass.size() > 0 && passCount < maxPasses){
		lastPass.clear();
		//cout << "\n passCount = " << passCount << endl;
		swapCount = 0;
		vector<int> toCheck(swappable);
		while(toCheck.size() > 0){
			int checkingIndex(kiss()%toCheck.size());
			int checking(toCheck[checkingIndex]);
			for(unsigned int p(0); p < swappable.size(); p++){
				if(swappable[p] == checking || b[swappable[p]].parentIndex == b[checking].parentIndex)
					continue;
				if(intersectingSubtrees(b, swappable[p], checking))
					continue;
				//cout << "(ps" << passCount << "c" << checking << "p" << p << ") "; cout.flush();
				vector<Bifurcation2D> bCopy(b);
				regraft(bCopy, checking, swappable[p]);
				if(minimumTipRatio(bCopy) < tipRatioThresh)
					continue;
				bool reportSettle(false && passCount == 2 && checking == 16 && p == 8);
				if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy, reportSettle, settleThresh, settleIts)){
					vector<BranchPoint2D> bc(consolidateDegenerates(bCopy, consolidateThresh));
					double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
					if(mBest > m){
						//cout << "\n " << passCount << " accepting regraft of " << checking << " with " << swappable[p]
						//	<< "\n\t to get from measure " << mBest << " to measure " << m << endl;
						swapCount++;
						mBest = m;
						bBest = bc;
						bBi = b;
						b = bCopy;
						lastPass.push_back(Pair<int>(checking, swappable[p]));
					}
				}
			}
			toCheck.erase(toCheck.begin() + checkingIndex);
		}
		passCount++;
		for(unsigned int i(0); i < lastPass.size(); i++){
			for(unsigned int j(i + 1); j < lastPass.size(); j++){
				if(lastPass[i].sameUnordered(lastPass[j])){
					lastPass.erase(lastPass.begin() + j);
					lastPass.erase(lastPass.begin() + i);
					i--;
					j = lastPass.size();
				}
			}
		}
	}
	swapCount = 1;
	double **seps = new double*[b.size()];
	for(unsigned int i(0); i < b.size(); i++)
		seps[i] = new double[b.size()];
	while(swapCount > 0 && passCount < maxPasses){
		for(unsigned int i(1); i < swappable.size(); i++){ // skipping the heart
			for(unsigned int j(i + 1); j < b.size(); j++)
				seps[i][j] = seps[j][i] = separation2D(b[i], b[j]);
		}
		swapCount = 0;
		vector<int> toCheck(swappable);
		while(toCheck.size() > 0){
			int checkingIndex(kiss()%toCheck.size());
			int checking(toCheck[checkingIndex]);
			toCheck.erase(toCheck.begin() + checkingIndex);
			vector<int> toCheckWith(toCheck);
			while(toCheckWith.size() > 0){ // not the heart
				int cIndex(kiss()%toCheckWith.size());
				int c(toCheckWith[cIndex]);
				toCheckWith.erase(toCheckWith.begin() + cIndex);
				if(seps[checking][c] < sepThresh){
					//cout << "(ps" << passCount << "ch" << checking << "c" << c << ") "; cout.flush();
					vector<Bifurcation2D> bCopy(b), bCopy2(b);
					bool reportSettle(false && passCount == 2 && checking == 15 && c == 16);
					if(regraftSingle(bCopy, checking, c)){
						if(minimumTipRatio(bCopy) < tipRatioThresh)
							continue;
						//drawBifurcationTree("outputs/oBGCGNew_1pre" + makeString(checking) + "k" + makeString(c) + "_" + makeString(globalTemp) + ".png", bCopy, 300, 300, 10);
						if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy, reportSettle, settleThresh, settleIts)){
							vector<BranchPoint2D> bc(consolidateDegenerates(bCopy, consolidateThresh));
							double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
							//if(globalTemp < 100){
							//	string sss("");
							//	if(mBest > m)
							//		sss = "*";
							//	drawTree("outputs/oBGCGNew" + makeString(checking) + "k" + makeString(c) + "_" + makeString(globalTemp) + sss + ".png", bc, 300, 300, 10);
							//	globalTemp++;
							//}
							if(mBest > m){
								swapCount++;
								mBest = m;
								bBest = bc;
								bBi = b;
								b = bCopy;
							}
						}
					}
					if(regraftSingle(bCopy2, checking, c, false)){
						if(minimumTipRatio(bCopy) < tipRatioThresh)
							continue;
						//drawBifurcationTree("outputs/oBGCGNew_2pre" + makeString(checking) + "k" + makeString(c) + "_" + makeString(globalTemp) + ".png", bCopy2, 300, 300, 10);
						if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy2, reportSettle, settleThresh, settleIts)){
							vector<BranchPoint2D> bc(consolidateDegenerates(bCopy2, consolidateThresh));
							double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
							if(mBest > m){
								swapCount++;
								mBest = m;
								bBest = bc;
								bBi = b;
								b = bCopy2;
							}
						}
					}
				}
			}
		}
		passCount++;
	}

	delete[] seps;
	if(passCount == maxPasses){
		cout << "\n optimizeByGlobalConsolidatedGrafting(): after " << passCount << " passes, still performed " << swapCount << " swaps." << endl;
	}
	return bBest;
}

// assumes the heart is 0
vector<BranchPoint2D> optimizeByGlobalConsolidatedGraftingAndTracking(vector<vector<Bifurcation2D> > &bestBs, vector<double> &bestBMeasures, const vector<Bifurcation2D> &bOrig, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient, double settleThresh, double consolidateThresh, double sepThresh, int maxPasses = 10000, int settleIts = 10000){
	vector<Bifurcation2D> b(bOrig);
	introduceBifurcationsFromSources(b);
	settleBifurcations(b, false, settleThresh, settleIts);
	vector<BranchPoint2D> bBest(consolidateDegenerates(b, consolidateThresh));
	double mBest(tripletMeasure(bBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	int swapCount(1), passCount(0);
	vector<int> swappable;
	for(unsigned int i(1); i < b.size(); i++) // not the heart
		swappable.push_back(i);
	vector<Pair<int> > lastPass(1, Pair<int>(-1, -1));
	while(lastPass.size() > 0 && passCount < maxPasses){
		lastPass.clear();
		//cout << "\n passCount = " << passCount << endl;
		swapCount = 0;
		vector<int> toCheck(swappable);
		while(toCheck.size() > 0){
			int checkingIndex(kiss()%toCheck.size());
			int checking(toCheck[checkingIndex]);
			for(unsigned int p(0); p < swappable.size(); p++){
				if(swappable[p] == checking || b[swappable[p]].parentIndex == b[checking].parentIndex)
					continue;
				if(intersectingSubtrees(b, swappable[p], checking))
					continue;
				//cout << "(ps" << passCount << "c" << checking << "p" << p << ") "; cout.flush();
				vector<Bifurcation2D> bCopy(b);
				regraft(bCopy, checking, swappable[p]);
				bool reportSettle(false && passCount == 2 && checking == 16 && p == 8);
				if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy, reportSettle, settleThresh, settleIts)){
					vector<BranchPoint2D> bc(consolidateDegenerates(bCopy, consolidateThresh));
					double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));

					if(!equivalentTrees(bestBs, bCopy))
						addTree(bestBs, bestBMeasures, bCopy, m);

					if(mBest > m){
						//cout << "\n " << passCount << " accepting regraft of " << checking << " with " << swappable[p]
						//	<< "\n\t to get from measure " << mBest << " to measure " << m << endl;
						swapCount++;
						mBest = m;
						bBest = bc;
						b = bCopy;
						lastPass.push_back(Pair<int>(checking, swappable[p]));
					}
				}
			}
			toCheck.erase(toCheck.begin() + checkingIndex);
		}
		passCount++;
		for(unsigned int i(0); i < lastPass.size(); i++){
			for(unsigned int j(i + 1); j < lastPass.size(); j++){
				if(lastPass[i].sameUnordered(lastPass[j])){
					lastPass.erase(lastPass.begin() + j);
					lastPass.erase(lastPass.begin() + i);
					i--;
					j = lastPass.size();
				}
			}
		}
	}
	swapCount = 1;
	double **seps = new double*[b.size()];
	for(unsigned int i(0); i < b.size(); i++)
		seps[i] = new double[b.size()];
	while(swapCount > 0 && passCount < maxPasses){
		for(unsigned int i(1); i < swappable.size(); i++){ // skipping the heart
			for(unsigned int j(i + 1); j < b.size(); j++)
				seps[i][j] = seps[j][i] = separation2D(b[i], b[j]);
		}
		swapCount = 0;
		vector<int> toCheck(swappable);
		while(toCheck.size() > 0){
			int checkingIndex(kiss()%toCheck.size());
			int checking(toCheck[checkingIndex]);
			toCheck.erase(toCheck.begin() + checkingIndex);
			vector<int> toCheckWith(toCheck);
			while(toCheckWith.size() > 0){ // not the heart
				int cIndex(kiss()%toCheckWith.size());
				int c(toCheckWith[cIndex]);
				toCheckWith.erase(toCheckWith.begin() + cIndex);
				if(seps[checking][c] < sepThresh){
					//cout << "(ps" << passCount << "ch" << checking << "c" << c << ") "; cout.flush();
					vector<Bifurcation2D> bCopy(b), bCopy2(b);
					bool reportSettle(false && passCount == 2 && checking == 15 && c == 16);
					if(regraftSingle(bCopy, checking, c)){
						//drawBifurcationTree("outputs/oBGCGNew_1pre" + makeString(checking) + "k" + makeString(c) + "_" + makeString(globalTemp) + ".png", bCopy, 300, 300, 10);
						if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy, reportSettle, settleThresh, settleIts)){
							vector<BranchPoint2D> bc(consolidateDegenerates(bCopy, consolidateThresh));
							double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
							if(!equivalentTrees(bestBs, bCopy))
								addTree(bestBs, bestBMeasures, bCopy, m);
							//if(globalTemp < 100){
							//	string sss("");
							//	if(mBest > m)
							//		sss = "*";
							//	drawTree("outputs/oBGCGNew" + makeString(checking) + "k" + makeString(c) + "_" + makeString(globalTemp) + sss + ".png", bc, 300, 300, 10);
							//	globalTemp++;
							//}
							if(mBest > m){
								swapCount++;
								mBest = m;
								bBest = bc;
								b = bCopy;
							}
						}
					}
					if(regraftSingle(bCopy2, checking, c, false)){
						//drawBifurcationTree("outputs/oBGCGNew_2pre" + makeString(checking) + "k" + makeString(c) + "_" + makeString(globalTemp) + ".png", bCopy2, 300, 300, 10);
						if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy2, reportSettle, settleThresh, settleIts)){
							vector<BranchPoint2D> bc(consolidateDegenerates(bCopy2, consolidateThresh));
							double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
							if(mBest > m){
								swapCount++;
								mBest = m;
								bBest = bc;
								b = bCopy2;
							}
						}
					}
				}
			}
		}
		passCount++;
	}

	delete[] seps;
	if(passCount == maxPasses){
		cout << "\n optimizeByGlobalConsolidatedGrafting(): after " << passCount << " passes, still performed " << swapCount << " swaps." << endl;
	}
	return bBest;
}

// assumes the heart is 0
vector<BranchPoint3D> optimizeByGlobalConsolidatedGraftingNEW(const vector<Bifurcation3D> &bOrig, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient, double settleThresh, double consolidateThresh, double sepThresh, int maxPasses = 10000, int settleIts = 10000){
	vector<Bifurcation3D> b(bOrig);
	introduceBifurcationsFromSources(b);
	settleBifurcations(b, false, settleThresh, settleIts);
	vector<BranchPoint3D> bBest(consolidateDegenerates(b, consolidateThresh));
	double mBest(tripletMeasure(bBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	int swapCount(1), passCount(0);
	vector<int> swappable;
	for(unsigned int i(1); i < b.size(); i++) // not the heart
		swappable.push_back(i);
	vector<Pair<int> > lastPass(1, Pair<int>(-1, -1));
	while(lastPass.size() > 0 && passCount < maxPasses){
		lastPass.clear();
		//cout << "\n passCount = " << passCount << endl;
		swapCount = 0;
		vector<int> toCheck(swappable);
		while(toCheck.size() > 0){
			int checkingIndex(kiss()%toCheck.size());
			int checking(toCheck[checkingIndex]);
			for(unsigned int p(0); p < swappable.size(); p++){
				if(swappable[p] == checking || b[swappable[p]].parentIndex == b[checking].parentIndex)
					continue;
				if(intersectingSubtrees(b, swappable[p], checking))
					continue;
				//cout << "(ps" << passCount << "c" << checking << "p" << p << ") "; cout.flush();
				vector<Bifurcation3D> bCopy(b);
				regraft(bCopy, checking, swappable[p]);
				bool reportSettle(false && passCount == 2 && checking == 16 && p == 8);
				if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy, reportSettle, settleThresh, settleIts)){
					vector<BranchPoint3D> bc(consolidateDegenerates(bCopy, consolidateThresh));
					double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
					if(mBest > m){
						//cout << "\n " << passCount << " accepting regraft of " << checking << " with " << swappable[p]
						//	<< "\n\t to get from measure " << mBest << " to measure " << m << endl;
						swapCount++;
						mBest = m;
						bBest = bc;
						b = bCopy;
						lastPass.push_back(Pair<int>(checking, swappable[p]));
					}
				}
			}
			toCheck.erase(toCheck.begin() + checkingIndex);
		}
		passCount++;
		for(unsigned int i(0); i < lastPass.size(); i++){
			for(unsigned int j(i + 1); j < lastPass.size(); j++){
				if(lastPass[i].sameUnordered(lastPass[j])){
					lastPass.erase(lastPass.begin() + j);
					lastPass.erase(lastPass.begin() + i);
					i--;
					j = lastPass.size();
				}
			}
		}
	}
	swapCount = 1;
	double **seps = new double*[b.size()];
	for(unsigned int i(0); i < b.size(); i++)
		seps[i] = new double[b.size()];
	while(swapCount > 0 && passCount < maxPasses){
		for(unsigned int i(1); i < swappable.size(); i++){ // skipping the heart
			for(unsigned int j(i + 1); j < b.size(); j++)
				seps[i][j] = seps[j][i] = separation3D(b[i], b[j]);
		}
		swapCount = 0;
		vector<int> toCheck(swappable);
		while(toCheck.size() > 0){
			int checkingIndex(kiss()%toCheck.size());
			int checking(toCheck[checkingIndex]);
			toCheck.erase(toCheck.begin() + checkingIndex);
			vector<int> toCheckWith(toCheck);
			while(toCheckWith.size() > 0){ // not the heart
				int cIndex(kiss()%toCheckWith.size());
				int c(toCheckWith[cIndex]);
				toCheckWith.erase(toCheckWith.begin() + cIndex);
				if(seps[checking][c] < sepThresh){
					//cout << "(ps" << passCount << "ch" << checking << "c" << c << ") "; cout.flush();
					vector<Bifurcation3D> bCopy(b), bCopy2(b);
					bool reportSettle(false && passCount == 2 && checking == 15 && c == 16);
					if(regraftSingle(bCopy, checking, c)){
						//drawBifurcationTree("outputs/oBGCGNew_1pre" + makeString(checking) + "k" + makeString(c) + "_" + makeString(globalTemp) + ".png", bCopy, 300, 300, 10);
						if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy, reportSettle, settleThresh, settleIts)){
							vector<BranchPoint3D> bc(consolidateDegenerates(bCopy, consolidateThresh));
							double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
							//if(globalTemp < 100){
							//	string sss("");
							//	if(mBest > m)
							//		sss = "*";
							//	drawTree("outputs/oBGCGNew" + makeString(checking) + "k" + makeString(c) + "_" + makeString(globalTemp) + sss + ".png", bc, 300, 300, 10);
							//	globalTemp++;
							//}
							if(mBest > m){
								swapCount++;
								mBest = m;
								bBest = bc;
								b = bCopy;
							}
						}
					}
					if(regraftSingle(bCopy2, checking, c, false)){
						//drawBifurcationTree("outputs/oBGCGNew_2pre" + makeString(checking) + "k" + makeString(c) + "_" + makeString(globalTemp) + ".png", bCopy2, 300, 300, 10);
						if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy2, reportSettle, settleThresh, settleIts)){
							vector<BranchPoint3D> bc(consolidateDegenerates(bCopy2, consolidateThresh));
							double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
							if(mBest > m){
								swapCount++;
								mBest = m;
								bBest = bc;
								b = bCopy2;
							}
						}
					}
				}
			}
		}
		passCount++;
	}

	delete[] seps;
	if(passCount == maxPasses){
		cout << "\n optimizeByGlobalConsolidatedGrafting(): after " << passCount << " passes, still performed " << swapCount << " swaps." << endl;
	}
	return bBest;
}

// assumes the heart is 0
vector<BranchPoint3D> optimizeByGlobalConsolidatedGrafting(const vector<Bifurcation3D> &bOrig, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient, double settleThresh, double consolidateThresh, int maxPasses = 10000, int settleIts = 10000){
	vector<Bifurcation3D> b(bOrig);
	introduceBifurcationsFromSources(b);
	settleBifurcations(b, false, settleThresh, settleIts);
	vector<BranchPoint3D> bBest(consolidateDegenerates(b, consolidateThresh));
	double mBest(tripletMeasure(bBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	int swapCount(1), passCount(0);
	vector<int> swappable;
	for(unsigned int i(1); i < b.size(); i++) // not the heart
		swappable.push_back(i);
	vector<Pair<int> > lastPass(1, Pair<int>(-1, -1));
	while(lastPass.size() > 0 && passCount < maxPasses){
		lastPass.clear();
		//cout << "\n passCount = " << passCount << endl;
		swapCount = 0;
		vector<int> toCheck(swappable);
		while(toCheck.size() > 0){
			int checkingIndex(kiss()%toCheck.size());
			int checking(toCheck[checkingIndex]);
			for(unsigned int p(0); p < swappable.size(); p++){
				if(swappable[p] == checking || b[swappable[p]].parentIndex == b[checking].parentIndex)
					continue;
				if(intersectingSubtrees(b, swappable[p], checking))
					continue;
				//cout << "(ps" << passCount << "c" << checking << "p" << p << ") "; cout.flush();
				vector<Bifurcation3D> bCopy(b);
				regraft(bCopy, checking, swappable[p]);
				bool reportSettle(false && passCount == 2 && checking == 16 && p == 8);
				if(settleBifurcationsThresh(2, mBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, bCopy, reportSettle, settleThresh, settleIts)){
					vector<BranchPoint3D> bc(consolidateDegenerates(bCopy, consolidateThresh));
					double m(tripletMeasure(bc, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
					if(mBest > m){
						//cout << "\n " << passCount << " accepting regraft of " << checking << " with " << swappable[p]
						//	<< "\n\t to get from measure " << mBest << " to measure " << m << endl;
						swapCount++;
						mBest = m;
						bBest = bc;
						b = bCopy;
						lastPass.push_back(Pair<int>(checking, swappable[p]));
					}//else
						//regraft(b, checking, swappable[p]);
				}
			}
			toCheck.erase(toCheck.begin() + checkingIndex);
		}
		passCount++;
		for(unsigned int i(0); i < lastPass.size(); i++){
			for(unsigned int j(i + 1); j < lastPass.size(); j++){
				if(lastPass[i].sameUnordered(lastPass[j])){
					lastPass.erase(lastPass.begin() + j);
					lastPass.erase(lastPass.begin() + i);
					i--;
					j = lastPass.size();
				}
			}
		}
	}
	if(passCount == maxPasses){
		cout << "\n optimizeByGlobalConsolidatedGrafting(): after " << passCount << " passes, still performed " << swapCount << " swaps." << endl;
	}
	return bBest;
}

vector<BranchPoint2D> optimizeByGlobalConsolidatedGrafting(int numOpt, const vector<Bifurcation2D> &bOrig, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient, double settleThresh, double consolidateThresh, int maxPasses = 10000, int settleIts = 10000){
	vector<BranchPoint2D> bBest(optimizeByGlobalConsolidatedGrafting(bOrig, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, maxPasses, settleIts));
	double measureBest(tripletMeasure(bBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	for(int i(1); i < numOpt; i++){
		//cout << "oBGCG(" << i << ") "; cout.flush();
		//kissprint();
		vector<BranchPoint2D> bTemp(optimizeByGlobalConsolidatedGrafting(bOrig, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, maxPasses, settleIts));
		double measureTemp(tripletMeasure(bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
		if(measureBest > measureTemp){
			measureBest = measureTemp;
			bBest = bTemp;
		}
	}
	return bBest;
}

vector<BranchPoint2D> optimizeByGlobalConsolidatedGraftingWithSingle(int numOpt, const vector<Bifurcation2D> &bOrig, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient, double settleThresh, double consolidateThresh, double sepThresh, int maxPasses = 10000, int settleIts = 10000){
	vector<BranchPoint2D> bBest(optimizeByGlobalConsolidatedGraftingNEW(bOrig, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, sepThresh, maxPasses, settleIts));
	double measureBest(tripletMeasure(bBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	for(int i(1); i < numOpt; i++){
		//cout << "oBGCG(" << i << ") "; cout.flush();
		//kissprint();
		vector<BranchPoint2D> bTemp(optimizeByGlobalConsolidatedGraftingNEW(bOrig, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, sepThresh, maxPasses, settleIts));
		double measureTemp(tripletMeasure(bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
		if(measureBest > measureTemp){
			measureBest = measureTemp;
			bBest = bTemp;
		}
	}
	return bBest;
}

vector<BranchPoint2D> optimizeByGlobalConsolidatedGraftingWithSingleAndTipRatioThresh(int numOpt, const vector<Bifurcation2D> &bOrig, vector<Bifurcation2D> &bBi, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient, double settleThresh, double consolidateThresh, double sepThresh, double tipRatioThresh, int maxPasses = 10000, int settleIts = 10000){
	vector<BranchPoint2D> bBest(optimizeByGlobalConsolidatedGraftingTipRatioThresh(bOrig, bBi, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, sepThresh, tipRatioThresh, maxPasses, settleIts));
	double measureBest(tripletMeasure(bBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	for(int i(1); i < numOpt; i++){
		//cout << "oBGCG(" << i << ") "; cout.flush();
		//kissprint();
		vector<Bifurcation2D> bBiTemp;
		vector<BranchPoint2D> bTemp(optimizeByGlobalConsolidatedGraftingTipRatioThresh(bOrig, bBiTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, sepThresh, tipRatioThresh, maxPasses, settleIts));
		double measureTemp(tripletMeasure(bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
		if(measureBest > measureTemp){
			measureBest = measureTemp;
			bBest = bTemp;
			bBi = bBiTemp;
		}
	}
	return bBest;
}

vector<BranchPoint2D> optimizeByGlobalConsolidatedGraftingWithSingleAndTracking(vector<vector<Bifurcation2D> > &bestBs, vector<double> &bestBMeasures, int numOpt, const vector<Bifurcation2D> &bOrig, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient, double settleThresh, double consolidateThresh, double sepThresh, int maxPasses = 10000, int settleIts = 10000){
	vector<BranchPoint2D> bBest(optimizeByGlobalConsolidatedGraftingAndTracking(bestBs, bestBMeasures, bOrig, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, sepThresh, maxPasses, settleIts));
	double measureBest(tripletMeasure(bBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	for(int i(1); i < numOpt; i++){
		//cout << "oBGCG(" << i << ") "; cout.flush();
		//kissprint();
		vector<BranchPoint2D> bTemp(optimizeByGlobalConsolidatedGraftingAndTracking(bestBs, bestBMeasures, bOrig, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, sepThresh, maxPasses, settleIts));
		double measureTemp(tripletMeasure(bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
		if(measureBest > measureTemp){
			measureBest = measureTemp;
			bBest = bTemp;
		}
	}
	return bBest;
}

vector<BranchPoint3D> optimizeByGlobalConsolidatedGraftingWithSingle(int numOpt, const vector<Bifurcation3D> &bOrig, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient, double settleThresh, double consolidateThresh, double sepThresh, int maxPasses = 10000, int settleIts = 10000){
	vector<BranchPoint3D> bBest(optimizeByGlobalConsolidatedGraftingNEW(bOrig, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, sepThresh, maxPasses, settleIts));
	double measureBest(tripletMeasure(bBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	for(int i(1); i < numOpt; i++){
		//cout << "oBGCG(" << i << ") "; cout.flush();
		//kissprint();
		vector<BranchPoint3D> bTemp(optimizeByGlobalConsolidatedGraftingNEW(bOrig, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, sepThresh, maxPasses, settleIts));
		double measureTemp(tripletMeasure(bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
		if(measureBest > measureTemp){
			measureBest = measureTemp;
			bBest = bTemp;
		}
	}
	return bBest;
}

vector<BranchPoint3D> optimizeByGlobalConsolidatedGrafting(int numOpt, const vector<Bifurcation3D> &bOrig, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient, double settleThresh, double consolidateThresh, int maxPasses = 10000, int settleIts = 10000){
	vector<BranchPoint3D> bBest(optimizeByGlobalConsolidatedGrafting(bOrig, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, maxPasses, settleIts));
	double measureBest(tripletMeasure(bBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	for(int i(1); i < numOpt; i++){
		//cout << "oBGCG(" << i << ") "; cout.flush();
		//kissprint();
		vector<BranchPoint3D> bTemp(optimizeByGlobalConsolidatedGrafting(bOrig, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, maxPasses, settleIts));
		double measureTemp(tripletMeasure(bTemp, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
		if(measureBest > measureTemp){
			measureBest = measureTemp;
			bBest = bTemp;
		}
	}
	return bBest;
}

void greedyConsolidatedSuccessRate(){
	kisset(6577, 880, 4353);
	vector<int> targetServiceVolumes, numToContinue, numRandSeeds; // all quantities should be ascending
	for(int i(0); i < 2; i++) targetServiceVolumes.push_back(i + 4);
	for(int i(0); i < 2; i++) numToContinue.push_back(2*i + 1);
	for(int i(0); i < 2; i++) numRandSeeds.push_back(3*i);
	bool doDraw(false);
	int numNetworks(10), numTempIncs(3), annealingIts(20), numTopKeepMult(2), numIndAnnMult(5);
	int maxFailedInsertions(100000), numBest(200), maxFailedNetworks(10000);
	double nearestRadius(1.0), settleThresh(1.0e-3), consolidateThresh(0.01), initT(1.0);
	double totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0);
	string uniquer("outputs/greedySuccessRate4-6_");
	int width(300), height(300), border(12);

	string outFileFn(uniquer + "summary.txt");
	ofstream outFile(outFileFn.c_str());
	outFile << "\n rescaledDerivatives()\n uniquer = " << uniquer
		<< "\n maxFailedNetworks = " << maxFailedNetworks
		<< "\n numNetworks = " << numNetworks
		<< "\n numRandSeeds = " << makeVectorString(numRandSeeds)
		<< "\n numToContinue = " << makeVectorString(numToContinue)
		<< "\n targetServiceVolumes = " << makeVectorString(targetServiceVolumes)
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n numBest = " << numBest
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< endl;

	int ****greedyRanks = new int***[targetServiceVolumes.size()];
	int **graftRanks = new int*[targetServiceVolumes.size()];
	int ****coolRanks = new int***[targetServiceVolumes.size()];
	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		greedyRanks[tsvi] = new int**[numNetworks];
		graftRanks[tsvi] = new int[numNetworks];
		coolRanks[tsvi] = new int**[numNetworks];
		double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes[tsvi])*acos(-1.0))/2.0), sizeY(sizeX);
		double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
		cout << "\n targetServiceVolumes = " << targetServiceVolumes[tsvi]
			<< " \t sizeX = sizeY = " << sizeX
			<< " \t heartx = " << heartx << " \t hearty = " << hearty << endl;
		for(int n(0); n < numNetworks; n++){
			greedyRanks[tsvi][n] = new int*[numToContinue.size()];
			coolRanks[tsvi][n] = new int*[numToContinue.size()];
			vector<Bifurcation2D> b(initializeEvenRandomSpacingServiceVolumes(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty, targetServiceVolumes[tsvi], maxFailedNetworks));
			vector<BranchPoint2D> *bestBs = new vector<BranchPoint2D>[numBest];
			double *bestBMeasures = new double[numBest];
			exhaustiveConsolidatedHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			vector<BranchPoint2D> bGraft(optimizeByGlobalConsolidatedGrafting(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh));
			double graftMeasure(tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
			if(graftMeasure < 0.999*bestBMeasures[0])
				cout << "\nGraft is better than exhaustive!" << endl;
			int matchConfigGraft(-1);
			for(int i(0); i < numBest && matchConfigGraft < 0; i++){
				if(equivalentTrees(bestBs[i], bGraft))
					matchConfigGraft = i;
			}
			graftRanks[tsvi][n] = matchConfigGraft + 1;
			cout << " graftRank = " << graftRanks[tsvi][n] << " \t graftMeasure = " << graftMeasure << endl;
			if(doDraw){
				drawTree(uniquer + "_n" + makeString(n) + "_exhaustive0.png", bestBs[0], width, height, border);
				drawTree(uniquer + "_n" + makeString(n) + "_graft_er" + makeString(graftRanks[tsvi][n]) + ".png", bGraft, width, height, border);
			}
			for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
				greedyRanks[tsvi][n][ntci] = new int[numRandSeeds.size()];
				coolRanks[tsvi][n][ntci] = new int[numRandSeeds.size()];
				cout << "\t numToContinue = " << numToContinue[ntci] << " \t (network " << n + 1 << " of " << numNetworks << " for " << targetServiceVolumes[tsvi] << " capillaries.)" << endl;
				for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++){
					cout << "\t\t numRandSeeds = " << numRandSeeds[nrsi] << endl;
					vector<BranchPoint2D> bGreedy(optimizeByPurelyDescendingConsolidatedSearch(numRandSeeds[nrsi], b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, numToContinue[ntci], settleThresh, consolidateThresh));
					double greedyMeasure(tripletMeasure(bGreedy, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
					if(greedyMeasure < 0.999*bestBMeasures[0])
						cout << "\nGreedy is better than exhaustive!" << endl;
					vector<Bifurcation2D> bCool(b);
					introduceBifurcationsFromSources(bCool);
					settleBifurcations(bCool, false, settleThresh);
					vector<BranchPoint2D> bCooled = simulatedConsolidatedCooling(bCool, initT, numTempIncs, numTopKeepMult*numToContinue[ntci], numIndAnnMult*numRandSeeds[nrsi], annealingIts, settleThresh, consolidateThresh, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
					double coolMeasure(tripletMeasure(bCooled, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
					if(coolMeasure < 0.999*bestBMeasures[0])
						cout << "\nCooling is better than exhaustive!" << endl;
					if(targetServiceVolumes[tsvi] < 7){
						time_t startExhaustive(time(NULL));

						//cout << "\n Exhaustive search took " << niceTime(1.0*(time(NULL) - startExhaustive)) << " to complete." << endl;
						//cout << "\n Exhaustive search optimal configuration is ";
						//if(equivalentTrees(bestBs[0], bGreedy))
						//	cout << "INDEED";
						//else
						//	cout << "NOT";
						//cout << " the same as the greedy search optimal tree." << endl;
						int matchConfigGreedy(-1);
						for(int i(0); i < numBest && matchConfigGreedy < 0; i++){
							if(equivalentTrees(bestBs[i], bGreedy))
								matchConfigGreedy = i;
						}
						greedyRanks[tsvi][n][ntci][nrsi] = matchConfigGreedy + 1;
						if(matchConfigGreedy > -1){
							cout << "\t\t\trank " << matchConfigGreedy + 1 << " \t exM[" << matchConfigGreedy << "]:" << bestBMeasures[matchConfigGreedy]
								<< " \t gM:" << greedyMeasure << " \t exM[0]:" << bestBMeasures[0] << endl;
						}else
							cout << "\t\t\trank -X-" << endl;
						int matchConfigCool(-1);
						for(int i(0); i < numBest && matchConfigCool < 0; i++){
							if(equivalentTrees(bestBs[i], bGreedy))
								matchConfigCool = i;
						}
						coolRanks[tsvi][n][ntci][nrsi] = matchConfigCool + 1;
						if(matchConfigGreedy > -1){
							cout << "\t\t\trank " << matchConfigCool + 1 << " \t exM[" << matchConfigCool << "]:" << bestBMeasures[matchConfigCool]
								<< " \t cM:" << coolMeasure << " \t exM[0]:" << bestBMeasures[0] << endl;
						}else
							cout << "\t\t\trank -X-" << endl;
					}else
						cout << "\t\t\t (No exhaustive hierarchy)" << endl;
					if(doDraw){
						drawTree(uniquer + "_n" + makeString(n) + "_greedy_ntc" + makeString(numToContinue[ntci]) + "_rs" + makeString(numRandSeeds[nrsi]) + "_er" + makeString(graftRanks[tsvi][n]) + ".png", bGreedy, width, height, border);
						drawTree(uniquer + "_n" + makeString(n) + "_cooled_ntk" + makeString(numTopKeepMult*numToContinue[ntci]) + "_nia" + makeString(numIndAnnMult*numRandSeeds[nrsi]) + "_er" + makeString(graftRanks[tsvi][n]) + ".png", bCooled, width, height, border);
					}
				}
			}
			delete[] bestBs;
			delete[] bestBMeasures;
		}
	}

	outFile << "\n\n v-targetServiceVolumes/numToContinue(mostRand)->";
	for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++)
		outFile << "\t" << numToContinue[ntci] << "\t +/- \t %-first\t %-unkn";
	outFile << endl;
	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		outFile << targetServiceVolumes[tsvi];
		for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
			vector<int> tempRanks;
			int firstsCount(0), unknownCount(0);
			for(int n(0); n < numNetworks; n++){
				int gRank(greedyRanks[tsvi][n][ntci][numRandSeeds.size() - 1]);
				if(gRank > 0){
					tempRanks.push_back(gRank);
					if(gRank == 1) firstsCount++;
				} else
					unknownCount++;
			}
			outFile << "\t" << mean(tempRanks) << "\t" << stDev(tempRanks) << "\t" << 100.0*double(firstsCount)/double(numNetworks) << "%" << "\t" << 100.0*double(unknownCount)/double(numNetworks) << "%";
		}
		outFile << endl;
	}

	outFile << "\n\n v-targetServiceVolumes/numRandSeeds(mostToCont)->";
	for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++)
		outFile << "\t" << numRandSeeds[nrsi] << "\t +/- \t %-first\t %-unkn";
	outFile << endl;
	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		outFile << targetServiceVolumes[tsvi];
		for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++){
			vector<int> tempRanks;
			int firstsCount(0);
			int unknownCount(0);
			for(int n(0); n < numNetworks; n++){
				int gRank(greedyRanks[tsvi][n][numToContinue.size() - 1][nrsi]);
				if(gRank > 0)
					tempRanks.push_back(gRank);
				if(gRank == 1) firstsCount++;
				if(gRank == 0) unknownCount++;
			}
			outFile << "\t" << mean(tempRanks) << "\t" << stDev(tempRanks) << "\t" << 100.0*double(firstsCount)/double(numNetworks) << "%" << "\t" << 100.0*double(unknownCount)/double(numNetworks) << "%";
		}
		outFile << endl;
	}

	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		outFile << "\n\n targetServiceVolumes = " << targetServiceVolumes[tsvi] << endl
			<< "v-numToContinue/numRandSeeds->";
		for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++)
			outFile << "\t" << numRandSeeds[nrsi] << "\t +/- \t %-first\t %-unkn";
		outFile << endl;
		for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
			outFile << numToContinue[ntci];
			for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++){
				vector<int> tempRanks;
				int firstsCount(0), unknownCount(0);
				for(int n(0); n < numNetworks; n++){
					int gRank(greedyRanks[tsvi][n][ntci][nrsi]);
					if(gRank > 0)
					tempRanks.push_back(gRank);
					if(gRank == 1) firstsCount++;
					if(gRank == 0) unknownCount++;
				}
				outFile << "\t" << mean(tempRanks) << "\t" << stDev(tempRanks) << "\t" << 100.0*double(firstsCount)/double(numNetworks) << "%" << "\t" << 100.0*double(unknownCount)/double(numNetworks) << "%";
			}
			outFile << endl;
		}
	}
	
	outFile << "\n\n v-targetServiceVolumes/numTopKeep(mostIndAnneal)->";
	for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++)
		outFile << "\t" << numTopKeepMult*numToContinue[ntci] << "\t +/- \t %-first\t %-unkn";
	outFile << endl;
	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		outFile << targetServiceVolumes[tsvi];
		for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
			vector<int> tempRanks;
			int firstsCount(0), unknownCount(0);
			for(int n(0); n < numNetworks; n++){
				int cRank(coolRanks[tsvi][n][ntci][numRandSeeds.size() - 1]);
				if(cRank > 0){
					tempRanks.push_back(cRank);
					if(cRank == 1) firstsCount++;
				}else
					unknownCount++;
			}
			outFile << "\t" << mean(tempRanks) << "\t" << stDev(tempRanks) << "\t" << 100.0*double(firstsCount)/double(numNetworks) << "%" << "\t" << 100.0*double(unknownCount)/double(numNetworks) << "%";
		}
		outFile << endl;
	}

	outFile << "\n\n v-targetServiceVolumes/indAnneal(mostTopKeep)->";
	for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++)
		outFile << "\t" << numIndAnnMult*numRandSeeds[nrsi] << "\t +/- \t %-first\t %-unkn";
	outFile << endl;
	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		outFile << targetServiceVolumes[tsvi];
		for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++){
			vector<int> tempRanks;
			int firstsCount(0);
			int unknownCount(0);
			for(int n(0); n < numNetworks; n++){
				int cRank(coolRanks[tsvi][n][numToContinue.size() - 1][nrsi]);
				if(cRank > 0)
					tempRanks.push_back(cRank);
				if(cRank == 1) firstsCount++;
				if(cRank == 0) unknownCount++;
			}
			outFile << "\t" << mean(tempRanks) << "\t" << stDev(tempRanks) << "\t" << 100.0*double(firstsCount)/double(numNetworks) << "%" << "\t" << 100.0*double(unknownCount)/double(numNetworks) << "%";
		}
		outFile << endl;
	}

	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		outFile << "\n\n targetServiceVolumes = " << targetServiceVolumes[tsvi] << endl
			<< "v-numTopKeep/numIndAnneal->";
		for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++)
			outFile << "\t" << numIndAnnMult*numRandSeeds[nrsi] << "\t +/- \t %-first\t %-unkn";
		outFile << endl;
		for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
			outFile << numTopKeepMult*numToContinue[ntci];
			for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++){
				vector<int> tempRanks;
				int firstsCount(0), unknownCount(0);
				for(int n(0); n < numNetworks; n++){
					int cRank(coolRanks[tsvi][n][ntci][nrsi]);
					if(cRank > 0)
						tempRanks.push_back(cRank);
					if(cRank == 1) firstsCount++;
					if(cRank == 0) unknownCount++;
				}
				outFile << "\t" << mean(tempRanks) << "\t" << stDev(tempRanks) << "\t" << 100.0*double(firstsCount)/double(numNetworks) << "%" << "\t" << 100.0*double(unknownCount)/double(numNetworks) << "%";
			}
			outFile << endl;
		}
	}

	outFile << "\n\n targetServiceVolumes\t <graftRank>\t +/-\t %-first\t %-unkn"<< endl;
	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		int firstsCount(0), unknownCount(0);
		vector<int> tempRanks;
		for(int n(0); n < numNetworks; n++){
			if(graftRanks[tsvi][n] > 0)
				tempRanks.push_back(graftRanks[tsvi][n]);
			if(graftRanks[tsvi][n] == 1) firstsCount++;
			if(graftRanks[tsvi][n] == 0) unknownCount++;
		}
		outFile << targetServiceVolumes[tsvi] << "\t" << mean(tempRanks) << "\t" << stDev(tempRanks) << "\t" << 100.0*double(firstsCount)/double(numNetworks) << "%" << "\t" << 100.0*double(unknownCount)/double(numNetworks) << "%" << endl;
	}


	outFile.close();
	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		for(int n(0); n < numNetworks; n++){
			for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
				delete[] greedyRanks[tsvi][n][ntci];
				delete[] coolRanks[tsvi][n][ntci];
			}
			delete[] greedyRanks[tsvi][n];
			delete[] coolRanks[tsvi][n];
		}
		delete[] greedyRanks[tsvi];
		delete[] coolRanks[tsvi];
		delete[] graftRanks[tsvi];
	}
	delete[] greedyRanks;
	delete[] coolRanks;
	delete[] graftRanks;
}

void greedyConsolidatedCompetition(){
	kisset(6577, 880, 4353);
	vector<int> targetServiceVolumes, numToContinue, numRandSeeds; // all quantities should be ascending
	for(int i(0); i < 3; i++) targetServiceVolumes.push_back(5*i + 10);
	for(int i(1); i < 2; i++) numToContinue.push_back(2*i + 1);
	for(int i(1); i < 2; i++) numRandSeeds.push_back(2*i);
	bool doDraw(true);
	int numNetworks(3), numTempIncs(2), annealingIts(10), numTopKeepMult(1), numIndAnnMult(2), graftPasses(100);
	int maxFailedInsertions(100000), maxFailedNetworks(10000);
	double nearestRadius(1.0), settleThresh(1.0e-3), consolidateThresh(0.1), initT(1.0);
	double totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0);
	string uniquer("outputs/greedyConsolidatedCompetition_");
	int width(300), height(300), border(12);

	string outFileFn(uniquer + "summary.txt");
	ofstream outFile(outFileFn.c_str());
	outFile << "\n rescaledDerivatives()\n uniquer = " << uniquer
		<< "\n maxFailedNetworks = " << maxFailedNetworks
		<< "\n numNetworks = " << numNetworks
		<< "\n numRandSeeds = " << makeVectorString(numRandSeeds)
		<< "\n numToContinue = " << makeVectorString(numToContinue)
		<< "\n targetServiceVolumes = " << makeVectorString(targetServiceVolumes)
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< "\n numTempIncs = " << numTempIncs << "\n annealingIts = " << annealingIts
		<< "\n numTopKeepMult = " << numTopKeepMult << "\n numIndAnnMult = " << numIndAnnMult
		<< "\n graftPasses = " << graftPasses
		<< "\n settleThresh = " << settleThresh << "\n consolidateThresh = " << consolidateThresh
		<< endl;

	double **guessMeasures = new double*[targetServiceVolumes.size()];
	double **graftMeasures = new double*[targetServiceVolumes.size()];
	double ****greedyMeasures = new double***[targetServiceVolumes.size()];
	double ****coolMeasures = new double***[targetServiceVolumes.size()];
	vector<BranchPoint2D> **guessTrees = new vector<BranchPoint2D>*[targetServiceVolumes.size()];
	vector<BranchPoint2D> **graftTrees = new vector<BranchPoint2D>*[targetServiceVolumes.size()];
	vector<BranchPoint2D> ****greedyTrees = new vector<BranchPoint2D>***[targetServiceVolumes.size()];
	vector<BranchPoint2D> ****coolTrees = new vector<BranchPoint2D>***[targetServiceVolumes.size()];
	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes[tsvi])*acos(-1.0))/2.0), sizeY(sizeX);
		double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
		cout << "\n targetServiceVolumes = " << targetServiceVolumes[tsvi]
			<< " \t sizeX = sizeY = " << sizeX
			<< " \t heartx = " << heartx << " \t hearty = " << hearty << endl;
		guessMeasures[tsvi] = new double[numNetworks];
		graftMeasures[tsvi] = new double[numNetworks];
		greedyMeasures[tsvi] = new double**[numNetworks];
		coolMeasures[tsvi] = new double**[numNetworks];
		guessTrees[tsvi] = new vector<BranchPoint2D>[numNetworks];
		graftTrees[tsvi] = new vector<BranchPoint2D>[numNetworks];
		greedyTrees[tsvi] = new vector<BranchPoint2D>**[numNetworks];
		coolTrees[tsvi] = new vector<BranchPoint2D>**[numNetworks];
		for(int n(0); n < numNetworks; n++){
			vector<Bifurcation2D> b(initializeEvenRandomSpacingServiceVolumes(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty, targetServiceVolumes[tsvi], maxFailedNetworks));
			vector<Bifurcation2D> bCopy(b);
			introduceBifurcationsFromSources(bCopy);
			guessTrees[tsvi][n] = consolidateDegenerates(bCopy, consolidateThresh);
			guessMeasures[tsvi][n] = tripletMeasure(guessTrees[tsvi][n], totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			cout << "\n\t\t\t\t\t\t Guess measure = " << guessMeasures[tsvi][n] << endl;
			if(doDraw){
				drawTree(uniquer + makeString(targetServiceVolumes[tsvi]) + "_n" + makeString(n) + "_guess.png", guessTrees[tsvi][n], width, height, border);
			}
			time_t tempStart(time(NULL));
			graftTrees[tsvi][n] = optimizeByGlobalConsolidatedGrafting(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, graftPasses);
			graftMeasures[tsvi][n] = tripletMeasure(graftTrees[tsvi][n], totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			cout << "\n\t\t Graft took " << niceTime(difftime(time(NULL), tempStart)) << "\t\t Graft measure = " << graftMeasures[tsvi][n] << endl;
			if(doDraw){
				drawTree(uniquer + makeString(targetServiceVolumes[tsvi]) + "_n" + makeString(n) + "_graft.png", graftTrees[tsvi][n], width, height, border);
			}
			greedyMeasures[tsvi][n] = new double*[numToContinue.size()];
			coolMeasures[tsvi][n] = new double*[numToContinue.size()];
			greedyTrees[tsvi][n] = new vector<BranchPoint2D>*[numToContinue.size()];
			coolTrees[tsvi][n] = new vector<BranchPoint2D>*[numToContinue.size()];
			for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
				cout << "\t numToContinue = " << numToContinue[ntci] << " \t (network " << n + 1 << " of " << numNetworks << " for " << targetServiceVolumes[tsvi] << " capillaries.)" << endl;
				greedyMeasures[tsvi][n][ntci] = new double[numRandSeeds.size()];
				coolMeasures[tsvi][n][ntci] = new double[numRandSeeds.size()];
				greedyTrees[tsvi][n][ntci] = new vector<BranchPoint2D>[numRandSeeds.size()];
				coolTrees[tsvi][n][ntci] = new vector<BranchPoint2D>[numRandSeeds.size()];
				for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++){
					cout << "\t\t numRandSeeds = " << numRandSeeds[nrsi] << endl;
					tempStart = time(NULL);
					greedyTrees[tsvi][n][ntci][nrsi] = optimizeByPurelyDescendingConsolidatedSearch(numRandSeeds[nrsi], b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, numToContinue[ntci], settleThresh, consolidateThresh);
					greedyMeasures[tsvi][n][ntci][nrsi] = tripletMeasure(greedyTrees[tsvi][n][ntci][nrsi], totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
					cout << "\t\t Greedy took " << niceTime(difftime(time(NULL), tempStart)) << "\t\t Greedy measure = " << greedyMeasures[tsvi][n][ntci][nrsi] << endl;
					vector<Bifurcation2D> bCool(b);
					introduceBifurcationsFromSources(bCool);
					settleBifurcations(bCool, false, settleThresh);
					tempStart = time(NULL);
					coolTrees[tsvi][n][ntci][nrsi] = simulatedConsolidatedCooling(bCool, initT, numTempIncs, numTopKeepMult*numToContinue[ntci], numIndAnnMult*numRandSeeds[nrsi], annealingIts, settleThresh, consolidateThresh, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
					coolMeasures[tsvi][n][ntci][nrsi] = tripletMeasure(coolTrees[tsvi][n][ntci][nrsi], totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
					cout << "\t\t Cool took " << niceTime(difftime(time(NULL), tempStart)) << "\t\t Cool measure = " << coolMeasures[tsvi][n][ntci][nrsi] << endl;
					if(doDraw){
						drawTree(uniquer + makeString(targetServiceVolumes[tsvi]) + "_n" + makeString(n) + "_greedy_ntc" + makeString(numToContinue[ntci]) + "_rs" + makeString(numRandSeeds[nrsi]) + ".png", greedyTrees[tsvi][n][ntci][nrsi], width, height, border);
						drawTree(uniquer + makeString(targetServiceVolumes[tsvi]) + "_n" + makeString(n) + "_cooled_ntk" + makeString(numTopKeepMult*numToContinue[ntci]) + "_nia" + makeString(numIndAnnMult*numRandSeeds[nrsi]) + ".png", coolTrees[tsvi][n][ntci][nrsi], width, height, border);
					}
				}
			}
		}
	}

	outFile << "\n\n tsv\t n\t guess\t graft";
	for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
		for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++)
			outFile << "\t greed_c" << numToContinue[ntci] << "_r" << numRandSeeds[nrsi];
	}
	for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
		for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++)
			outFile << "\t cool_c" << numTopKeepMult*numToContinue[ntci] << "_r" << numIndAnnMult*numRandSeeds[nrsi];
	}
	outFile << endl;
	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		for(int n(0); n < numNetworks; n++){
			outFile << targetServiceVolumes[tsvi] << "\t" << n << "\t" << guessMeasures[tsvi][n] << "\t" << graftMeasures[tsvi][n];
			for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
				for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++)
					outFile << "\t" << greedyMeasures[tsvi][n][ntci][nrsi];
			}
			for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
				for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++)
					outFile << "\t" << coolMeasures[tsvi][n][ntci][nrsi];
			}
			outFile << endl;
		}
	}

	outFile << "\n\n tsv\t n\t guess\t graft";
	for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
		for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++)
			outFile << "\t greed_c" << numToContinue[ntci] << "_r" << numRandSeeds[nrsi];
	}
	for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
		for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++)
			outFile << "\t cool_c" << numTopKeepMult*numToContinue[ntci] << "_r" << numIndAnnMult*numRandSeeds[nrsi];
	}
	outFile << endl;
	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		for(int n(0); n < numNetworks; n++){
			outFile << targetServiceVolumes[tsvi] << "\t" << n;
			double minMeasure(guessMeasures[tsvi][n]);
			vector<BranchPoint2D> minTree(guessTrees[tsvi][n]);
			if(minMeasure > graftMeasures[tsvi][n]){
				minMeasure = graftMeasures[tsvi][n];
				minTree = graftTrees[tsvi][n];
			}
			for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
				for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++){
					if(minMeasure > greedyMeasures[tsvi][n][ntci][nrsi]){
						minMeasure = greedyMeasures[tsvi][n][ntci][nrsi];
						minTree = greedyTrees[tsvi][n][ntci][nrsi];
					}
					if(minMeasure > coolMeasures[tsvi][n][ntci][nrsi]){
						minMeasure = coolMeasures[tsvi][n][ntci][nrsi];
						minTree = coolTrees[tsvi][n][ntci][nrsi];
					}
				}
			}
			if(guessMeasures[tsvi][n] == minMeasure)
				outFile << "\t*";
			else if(equivalentTrees(guessTrees[tsvi][n], minTree)){
				if(guessMeasures[tsvi][n] < 1.001*minMeasure)
					outFile << "\to";
				else
					outFile << "\t|";
			}else
				outFile << "\tx";
			if(graftMeasures[tsvi][n] == minMeasure)
				outFile << "\t*";
			else if(equivalentTrees(graftTrees[tsvi][n], minTree)){
				if(graftMeasures[tsvi][n] < 1.001*minMeasure)
					outFile << "\to";
				else
					outFile << "\t|";
			}else
				outFile << "\tx";
			for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
				for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++){
					if(greedyMeasures[tsvi][n][ntci][nrsi] == minMeasure)
						outFile << "\t*";
					else if(equivalentTrees(greedyTrees[tsvi][n][ntci][nrsi], minTree)){
						if(greedyMeasures[tsvi][n][ntci][nrsi] < 1.001*minMeasure)
							outFile << "\to";
						else
							outFile << "\t|";
					}else
						outFile << "\tx";
				}
			}
			for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
				for(unsigned int nrsi(0); nrsi < numRandSeeds.size(); nrsi++){
					if(coolMeasures[tsvi][n][ntci][nrsi] == minMeasure)
						outFile << "\t*";
					else if(equivalentTrees(coolTrees[tsvi][n][ntci][nrsi], minTree)){
						if(coolMeasures[tsvi][n][ntci][nrsi] < 1.001*minMeasure)
							outFile << "\to";
						else
							outFile << "\t|";
					}else
						outFile << "\tx";
				}
			}
			outFile << endl;
		}
	}
	outFile.close();

	for(unsigned int tsvi(0); tsvi < targetServiceVolumes.size(); tsvi++){
		for(int n(0); n < numNetworks; n++){
			for(unsigned int ntci(0); ntci < numToContinue.size(); ntci++){
				delete[] greedyMeasures[tsvi][n][ntci];
				delete[] coolMeasures[tsvi][n][ntci];
				delete[] greedyTrees[tsvi][n][ntci];
				delete[] coolTrees[tsvi][n][ntci];
			}
			delete[] greedyMeasures[tsvi][n];
			delete[] coolMeasures[tsvi][n];
			delete[] greedyTrees[tsvi][n];
			delete[] coolTrees[tsvi][n];
		}
		delete[] greedyMeasures[tsvi];
		delete[] coolMeasures[tsvi];
		delete[] graftMeasures[tsvi];
		delete[] greedyTrees[tsvi];
		delete[] coolTrees[tsvi];
		delete[] graftTrees[tsvi];
	}
	delete[] greedyMeasures;
	delete[] coolMeasures;
	delete[] graftMeasures;
	delete[] greedyTrees;
	delete[] coolTrees;
	delete[] graftTrees;
}

string makeStringThree(double a, double b, double c){
	return makeString(ceil(a)) + makeString(ceil(b)) + makeString(ceil(c));
}

void exhaustiveOptimal(){
	kisset(29530, 13947, 20995);
	kisset(10394, 26468, 26630);
	string uniquer = "outputs/exhaustiveOptimal";
	int targetServiceVolumes(4), numBest(15), maxFailedInsertions(100000), maxFailedNetworks(1000),
		width(300), height(300), border(12);
	double totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0),
		nearestRadius(1.0);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	vector<Bifurcation2D> b(initializeEvenRandomSpacingServiceVolumes(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty, targetServiceVolumes, maxFailedNetworks));
	vector<Bifurcation2D> *bestBs = new vector<Bifurcation2D>[numBest];
	double *bestBMeasures = new double[numBest];
	totalLengthCoefficient = 1.0;
	avePathLengthCoefficient = 0.0;
	maxPathLengthCoefficient = 0.0;
	exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	for(int i(0); i < numBest; i++){
		drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_r" + makeString(i + 1) + ".png", bestBs[i], width, height, border);
		drawTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_r" + makeString(i + 1) + "_cons.png", consolidateDegenerates(bestBs[i], 0.01), width, height, border);
		cout << "configuration " << makeStringTreeOneLine(bestBs[i]) << " has total length " << tripletMeasure(bestBs[i], 1.0, 0.0, 0.0) << endl;
	}
	cout << "\n best tree: " << makeStringTreeOneLine(bestBs[0])
		<< "\nc" << makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n measure 100 = " << tripletMeasure(bestBs[0], 1.0, 0.0, 0.0)
		<< "\n measure 010 = " << tripletMeasure(bestBs[0], 0.0, 1.0, 0.0)
		<< "\n measure 001 = " << tripletMeasure(bestBs[0], 0.0, 0.0, 1.0)
		<< endl;
	cout << "\n second best tree: " << makeStringTreeOneLine(bestBs[1])
		<< "\nc" << makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n measure 100 = " << tripletMeasure(bestBs[1], 1.0, 0.0, 0.0)
		<< "\n measure 010 = " << tripletMeasure(bestBs[1], 0.0, 1.0, 0.0)
		<< "\n measure 001 = " << tripletMeasure(bestBs[1], 0.0, 0.0, 1.0)
		<< endl;
	cout << "\n third best tree: " << makeStringTreeOneLine(bestBs[2])
		<< "\nc" << makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n measure 100 = " << tripletMeasure(bestBs[2], 1.0, 0.0, 0.0)
		<< "\n measure 010 = " << tripletMeasure(bestBs[2], 0.0, 1.0, 0.0)
		<< "\n measure 001 = " << tripletMeasure(bestBs[2], 0.0, 0.0, 1.0)
		<< endl;
	totalLengthCoefficient = 0.0;
	avePathLengthCoefficient = 1.0;
	maxPathLengthCoefficient = 0.0;
	exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	for(int i(0); i < numBest; i++)
		drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_r" + makeString(i + 1) + ".png", bestBs[i], width, height, border);
	cout << "\n best tree: " << makeStringTreeOneLine(bestBs[0])
		<< "\nc" << makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n measure 100 = " << tripletMeasure(bestBs[0], 1.0, 0.0, 0.0)
		<< "\n measure 010 = " << tripletMeasure(bestBs[0], 0.0, 1.0, 0.0)
		<< "\n measure 001 = " << tripletMeasure(bestBs[0], 0.0, 0.0, 1.0)
		<< endl;
	totalLengthCoefficient = 0.0;
	avePathLengthCoefficient = 0.0;
	maxPathLengthCoefficient = 1.0;
	exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	for(int i(0); i < numBest; i++)
		drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_r" + makeString(i + 1) + ".png", bestBs[i], width, height, border);
	cout << "\n best tree: " << makeStringTreeOneLine(bestBs[0])
		<< "\nc" << makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n measure 100 = " << tripletMeasure(bestBs[0], 1.0, 0.0, 0.0)
		<< "\n measure 010 = " << tripletMeasure(bestBs[0], 0.0, 1.0, 0.0)
		<< "\n measure 001 = " << tripletMeasure(bestBs[0], 0.0, 0.0, 1.0)
		<< endl;
	totalLengthCoefficient = 1.0;
	avePathLengthCoefficient = 1.0;
	maxPathLengthCoefficient = 0.0;
	exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	for(int i(0); i < numBest; i++)
		drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_r" + makeString(i + 1) + ".png", bestBs[i], width, height, border);
	cout << "\n best tree: " << makeStringTreeOneLine(bestBs[0])
		<< "\nc" << makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n measure 100 = " << tripletMeasure(bestBs[0], 1.0, 0.0, 0.0)
		<< "\n measure 010 = " << tripletMeasure(bestBs[0], 0.0, 1.0, 0.0)
		<< "\n measure 001 = " << tripletMeasure(bestBs[0], 0.0, 0.0, 1.0)
		<< endl;
	totalLengthCoefficient = 1.0;
	avePathLengthCoefficient = 0.0;
	maxPathLengthCoefficient = 1.0;
	exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	for(int i(0); i < numBest; i++)
		drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_r" + makeString(i + 1) + ".png", bestBs[i], width, height, border);
	cout << "\n best tree: " << makeStringTreeOneLine(bestBs[0])
		<< "\nc" << makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n measure 100 = " << tripletMeasure(bestBs[0], 1.0, 0.0, 0.0)
		<< "\n measure 010 = " << tripletMeasure(bestBs[0], 0.0, 1.0, 0.0)
		<< "\n measure 001 = " << tripletMeasure(bestBs[0], 0.0, 0.0, 1.0)
		<< endl;
	totalLengthCoefficient = 0.0;
	avePathLengthCoefficient = 1.0;
	maxPathLengthCoefficient = 1.0;
	exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	for(int i(0); i < numBest; i++)
		drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_r" + makeString(i + 1) + ".png", bestBs[i], width, height, border);
	cout << "\n best tree: " << makeStringTreeOneLine(bestBs[0])
		<< "\nc" << makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n measure 100 = " << tripletMeasure(bestBs[0], 1.0, 0.0, 0.0)
		<< "\n measure 010 = " << tripletMeasure(bestBs[0], 0.0, 1.0, 0.0)
		<< "\n measure 001 = " << tripletMeasure(bestBs[0], 0.0, 0.0, 1.0)
		<< endl;
	totalLengthCoefficient = 1.0;
	avePathLengthCoefficient = 1.0;
	maxPathLengthCoefficient = 1.0;
	exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	for(int i(0); i < numBest; i++)
		drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_r" + makeString(i + 1) + ".png", bestBs[i], width, height, border);
	cout << "\n best tree: " << makeStringTreeOneLine(bestBs[0])
		<< "\nc" << makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n measure 100 = " << tripletMeasure(bestBs[0], 1.0, 0.0, 0.0)
		<< "\n measure 010 = " << tripletMeasure(bestBs[0], 0.0, 1.0, 0.0)
		<< "\n measure 001 = " << tripletMeasure(bestBs[0], 0.0, 0.0, 1.0)
		<< endl;

	delete[] bestBs;
	delete[] bestBMeasures;
}

// assumes the heart is 0
vector<Bifurcation2D> optimizeByTipSwapping(const vector<Bifurcation2D> &bOrig, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient){
	vector<Bifurcation2D> b(bOrig);
	introduceBifurcationsFromSources(b);
	settleBifurcations(b);
	vector<Bifurcation2D> bBest(b);
	double mBest(tripletMeasure(bBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	int swapCount(1), passCount(0);
	vector<int> swappable;
	for(unsigned int i(1); i < bOrig.size(); i++) // only tips, like bOrig less the heart
		swappable.push_back(i);
	while(swapCount > 0 && passCount < 10000){
		swapCount = 0;
		vector<int> toCheck(swappable);
		while(toCheck.size() > 0){
			int checkingIndex(kiss()%toCheck.size());
			int checking(toCheck[checkingIndex]);
			for(unsigned int p(0); p < swappable.size(); p++){
				if(swappable[p] == checking || b[swappable[p]].parentIndex == b[checking].parentIndex)
					continue;
				regraft(b, checking, swappable[p]);
				settleBifurcations(b);
				double m(tripletMeasure(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
				if(mBest > m){
					swapCount++;
					mBest = m;
					bBest = b;
				}else
					regraft(b, checking, swappable[p]);
			}
			toCheck.erase(toCheck.begin() + checkingIndex);
		}
		passCount++;
	}
	return bBest;
}

// assumes the heart is 0
vector<Bifurcation2D> optimizeByGlobalGrafting(const vector<Bifurcation2D> &bOrig, double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient){
	vector<Bifurcation2D> b(bOrig);
	introduceBifurcationsFromSources(b);
	settleBifurcations(b);
	vector<Bifurcation2D> bBest(b);
	double mBest(tripletMeasure(bBest, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	int swapCount(1), passCount(0);
	vector<int> swappable;
	for(unsigned int i(1); i < b.size(); i++) // not the heart
		swappable.push_back(i);
	while(swapCount > 0 && passCount < 10000){
		swapCount = 0;
		vector<int> toCheck(swappable);
		while(toCheck.size() > 0){
			int checkingIndex(kiss()%toCheck.size());
			int checking(toCheck[checkingIndex]);
			for(unsigned int p(0); p < swappable.size(); p++){
				if(swappable[p] == checking || b[swappable[p]].parentIndex == b[checking].parentIndex)
					continue;
				if(intersectingSubtrees(b, swappable[p], checking))
					continue;
				regraft(b, checking, swappable[p]);
				settleBifurcations(b);
				double m(tripletMeasure(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
				if(mBest > m){
					swapCount++;
					mBest = m;
					bBest = b;
				}else
					regraft(b, checking, swappable[p]);
			}
			toCheck.erase(toCheck.begin() + checkingIndex);
		}
		passCount++;
	}
	return bBest;
}

int findRank(const vector<Bifurcation2D> &b, const vector<Bifurcation2D> *ranks, int rankSize){
	for(int i(0); i < rankSize; i++){
		if(equivalentTrees(b, ranks[i]))
			return i + 1;
	}
	return -1;
}

void tipSwapping(){
	kisset(22442, 28410, 29790);
	string uniquer = "outputs/swapping_";
	int targetServiceVolumes(6), numBest(945), maxFailedInsertions(100000), maxFailedNetworks(1000),
		width(300), height(300), border(12);
	double totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0),
		nearestRadius(1.0);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	vector<Bifurcation2D> b(initializeEvenRandomSpacingServiceVolumes(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty, targetServiceVolumes, maxFailedNetworks));
	vector<Bifurcation2D> bGuess(b);
	introduceBifurcationsFromSources(bGuess);
	settleBifurcations(bGuess);
	drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_guess.png", bGuess, width, height, border);

	vector<Bifurcation2D> *bestBs0 = new vector<Bifurcation2D>[numBest];
	double *bestBMeasures0 = new double[numBest];
	totalLengthCoefficient = 1.0;
	avePathLengthCoefficient = 0.0;
	maxPathLengthCoefficient = 0.0;
	time_t start_temp(time(NULL));
	exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs0, bestBMeasures0, false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	cout << "\n exhaustive took " << niceTime(difftime(time(NULL), start_temp)) << endl;
	start_temp = time(NULL);
	vector<Bifurcation2D> bTips(optimizeByTipSwapping(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	cout << "\n tips took " << niceTime(difftime(time(NULL), start_temp)) << endl;
	start_temp = time(NULL);
	vector<Bifurcation2D> bGlob(optimizeByGlobalGrafting(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
	cout << "\n glob took " << niceTime(difftime(time(NULL), start_temp)) << endl;
	drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_tip.png", bTips, width, height, border);
	drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_glob.png", bGlob, width, height, border);
	int guessRank(findRank(bGuess, bestBs0, numBest));
	int tipsRank(findRank(bTips, bestBs0, numBest));
	int globRank(findRank(bGlob, bestBs0, numBest));
	if(tipsRank > -1){
		drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_tipExhEqu.png", bestBs0[tipsRank - 1], width, height, border);
		cout << "\n equivalent exhaustive configuration to tips has measure " << bestBMeasures0[tipsRank - 1] << endl;
	}
	if(globRank > -1){
		drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_globExhEqu.png", bestBs0[globRank - 1], width, height, border);
		cout << "\n equivalent exhaustive configuration to glob has measure " << bestBMeasures0[globRank - 1] << endl;
	}
	cout << "\n\nc" << makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< " starting with rank " << guessRank << " and measure " << tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n\t tips rank " << tipsRank << " and measure " << tripletMeasure(bTips, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n\t glob rank " << globRank << " and measure " << tripletMeasure(bGlob, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n\t best exhaustive configuration (rank 1) has measure " << bestBMeasures0[0]
		<< endl;
	delete[] bestBs0;
	delete[] bestBMeasures0;
	
	vector<Bifurcation2D> *bestBs1 = new vector<Bifurcation2D>[numBest];
	double *bestBMeasures1 = new double[numBest];
	totalLengthCoefficient = 0.0;
	avePathLengthCoefficient = 1.0;
	maxPathLengthCoefficient = 0.0;
	start_temp = time(NULL);
	exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs1, bestBMeasures1, false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	cout << "\n exhaustive took " << niceTime(difftime(time(NULL), start_temp)) << endl;
	start_temp = time(NULL);
	bTips = optimizeByTipSwapping(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	cout << "\n tips took " << niceTime(difftime(time(NULL), start_temp)) << endl;
	start_temp = time(NULL);
	bGlob = optimizeByGlobalGrafting(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	cout << "\n glob took " << niceTime(difftime(time(NULL), start_temp)) << endl;
	drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_tip.png", bTips, width, height, border);
	drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_glob.png", bGlob, width, height, border);
	guessRank = findRank(bGuess, bestBs1, numBest);
	tipsRank = findRank(bTips, bestBs1, numBest);
	globRank = findRank(bGlob, bestBs1, numBest);
	if(tipsRank > -1){
		drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_tipExhEqu.png", bestBs0[tipsRank - 1], width, height, border);
		cout << "\n equivalent exhaustive configuration to tips has measure " << bestBMeasures1[tipsRank - 1] << endl;
	}
	if(globRank > -1){
		drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_globExhEqu.png", bestBs0[globRank - 1], width, height, border);
		cout << "\n equivalent exhaustive configuration to glob has measure " << bestBMeasures1[globRank - 1] << endl;
	}
	cout << "\n\nc" << makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< " starting with rank " << guessRank << " and measure " << tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n\t tips rank " << tipsRank << " and measure " << tripletMeasure(bTips, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n\t glob rank " << globRank << " and measure " << tripletMeasure(bGlob, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n\t best exhaustive configuration (rank 1) has measure " << bestBMeasures1[0]
		<< endl;
	delete[] bestBs1;
	delete[] bestBMeasures1;

	vector<Bifurcation2D> *bestBs2 = new vector<Bifurcation2D>[numBest];
	double *bestBMeasures2 = new double[numBest];
	totalLengthCoefficient = 1.0;
	avePathLengthCoefficient = 1.0;
	maxPathLengthCoefficient = 0.0;
	start_temp = time(NULL);
	exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs2, bestBMeasures2, false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	cout << "\n exhaustive took " << niceTime(difftime(time(NULL), start_temp)) << endl;
	start_temp = time(NULL);
	bTips = optimizeByTipSwapping(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	cout << "\n tips took " << niceTime(difftime(time(NULL), start_temp)) << endl;
	start_temp = time(NULL);
	bGlob = optimizeByGlobalGrafting(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	cout << "\n glob took " << niceTime(difftime(time(NULL), start_temp)) << endl;
	drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_tip.png", bTips, width, height, border);
	drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_glob.png", bGlob, width, height, border);
	guessRank = findRank(bGuess, bestBs2, numBest);
	tipsRank = findRank(bTips, bestBs2, numBest);
	globRank = findRank(bGlob, bestBs2, numBest);
	if(tipsRank > -1){
		drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_tipExhEqu.png", bestBs0[tipsRank - 1], width, height, border);
		cout << "\n equivalent exhaustive configuration to tips has measure " << bestBMeasures2[tipsRank - 1] << endl;
	}
	if(globRank > -1){
		drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_c" + makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) + "_globExhEqu.png", bestBs0[globRank - 1], width, height, border);
		cout << "\n equivalent exhaustive configuration to glob has measure " << bestBMeasures2[globRank - 1] << endl;
	}
	cout << "\n\nc" << makeStringThree(totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< " starting with rank " << guessRank << " and measure " << tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n\t tips rank " << tipsRank << " and measure " << tripletMeasure(bTips, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n\t glob rank " << globRank << " and measure " << tripletMeasure(bGlob, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient)
		<< "\n\t best exhaustive configuration (rank 1) has measure " << bestBMeasures2[0]
		<< endl;
	delete[] bestBs2;
	delete[] bestBMeasures2;
}

void testTorres(){
	string uniquer("outputs/testTorres_");
	int width(300), height(300), border(12);
	vector<Bifurcation2D> b;
	b.push_back(Bifurcation2D());
	b.push_back(Bifurcation2D());
	b.push_back(Bifurcation2D());
	b.back().position[0] = 1.0;
	b.back().position[1] = 0.0;
	b.push_back(Bifurcation2D());
	b.back().position[0] = 0.0;
	b.back().position[1] = 1.0;
	b[0].childIndex.push_back(1);
	b[1].parentIndex = 0;
	b[1].childIndex.push_back(2);
	b[1].childIndex.push_back(3);
	b[2].parentIndex = b[3].parentIndex = 1;
	vector<int> constInd;
	constInd.push_back(0);
	constInd.push_back(2);
	constInd.push_back(3);
	vector<int> movingInd;
	movingInd.push_back(1);
	vector<double> w(3, 1.0);
	vector<Bifurcation2D> bCopy(b);
	vector<Bifurcation2D> bCopy2(b);
	drawBifurcationTree(uniquer + "4initial.png", b, width, height, border);
	TorresWeiszfeld(b, constInd, w, movingInd);
	drawBifurcationTree(uniquer + "4Torres.png", b, width, height, border);
	settleBifurcationsOLD(bCopy);
	drawBifurcationTree(uniquer + "4Fermat.png", bCopy, width, height, border);
	settleBifurcations(bCopy2);
	drawBifurcationTree(uniquer + "4Settle.png", bCopy2, width, height, border);

	b.clear(); bCopy.clear(); bCopy2.clear();  constInd.clear(); movingInd.clear(); w.clear();
	b.push_back(Bifurcation2D());
	b.back().position[0] = 0.5;
	b.back().position[1] = 0.0;
	b.push_back(Bifurcation2D());
	b.push_back(Bifurcation2D());
	b.back().position[0] = 1.0;
	b.back().position[1] = 1.0;
	b.push_back(Bifurcation2D());
	b.back().position[0] = 1.0;
	b.back().position[1] = 0.0;
	b.push_back(Bifurcation2D());
	b.back().position[0] = 0.0;
	b.back().position[1] = 1.0;
	addParent(b, 1, 2, 1.0, 0.5);
	addParent(b, 3, 4, 0.5, 1.0);
	addParent(b, 5, 6, 0.0, 0.5);
	b[7].parentIndex = 0;
	b[0].childIndex.push_back(7);
	bCopy = b;
	drawBifurcationTree(uniquer + "8initial.png", b, width, height, border);
	settleBifurcationsOLD(bCopy);
	drawBifurcationTree(uniquer + "8Fermat.png", bCopy, width, height, border);
	vector<int> allMovingInd;
	for(unsigned int i(0); i < bCopy.size(); i++){
		if(containedIn(i, allMovingInd))
			continue;
		if(bCopy[i].parentIndex > -1 && bCopy[i].childIndex.size() == 2){
			bool isInDegenerateGob(separation2D(bCopy[i], bCopy[bCopy[i].parentIndex]) <= 1.0e-6 && bCopy[bCopy[i].parentIndex].parentIndex > -1);
			for(unsigned int c(0); c < bCopy[i].childIndex.size(); c++)
				isInDegenerateGob = isInDegenerateGob || (separation2D(bCopy[i], bCopy[bCopy[i].childIndex[c]]) <= 1.0e-6 && bCopy[bCopy[i].childIndex[c]].childIndex.size() > 0);
			if(isInDegenerateGob){
				vector<int> cInd, mInd;
				findDegenerateGob(bCopy, i, cInd, mInd);
				cout << "\n identified degenerate gob consisting of " << makeVectorString(mInd)
					<< " constrained by " << makeVectorString(cInd) << endl;
				allMovingInd.insert(allMovingInd.end(), mInd.begin(), mInd.end());
			}
		}
	}
	for(int i(0); i < 5; i++)
		constInd.push_back(i);
	for(int i(5); i < 8; i++)
		movingInd.push_back(i);
	cout << "\n Should have found moving indices to be " << makeVectorString(movingInd) << " and constrained by " << makeVectorString(constInd) << endl;
	w = vector<double>(5, 1.0);
	bCopy2 = b;
	TorresWeiszfeld(b, constInd, w, movingInd);
	drawBifurcationTree(uniquer + "8Torres.png", b, width, height, border);
	TorresWeiszfeld(bCopy, constInd, w, movingInd);
	drawBifurcationTree(uniquer + "8TorresAfterFermat.png", bCopy, width, height, border);
	settleBifurcations(bCopy2);
	drawBifurcationTree(uniquer + "8Settle.png", bCopy2, width, height, border);
	
}

void testConsolidation(){
	kisset(123, 45678, 9012); // with crossing for 5 service volumes, with near-consolidation for 7 service volumes
	string uniquer = "outputs/consol_";
	int targetServiceVolumes(5), numBest(5), maxFailedInsertions(100000), maxFailedNetworks(1000),
		width(300), height(300), border(12);
	double totalLengthCoefficient(0.0),
		avePathLengthCoefficient(1.0),
		maxPathLengthCoefficient(0.0),
		nearestRadius(1.0),
		settleThresh(1.0e-3);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	vector<Bifurcation2D> b(initializeEvenRandomSpacingServiceVolumes(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty, targetServiceVolumes, maxFailedNetworks));
	//introduceBifurcationsFromSources(b);
	//settleBifurcations(b);
	//vector<Bifurcation2D> bOpt = optimizeByGlobalGrafting(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	vector<Bifurcation2D> *bestBs = new vector<Bifurcation2D>[numBest];
	double *bestBMeasures = new double[numBest];
	exhaustiveHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, true, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, true, settleThresh);
	vector<Bifurcation2D> bOpt = bestBs[0];
	drawBifurcationTree(uniquer + makeString(targetServiceVolumes) + "_Bi.png", bOpt, width, height, border);
	vector<BranchPoint2D> br(consolidateDegenerates(bOpt, settleThresh));
	drawTree(uniquer + makeString(targetServiceVolumes) + "_Br.png", br, width, height, border);
	vector<BranchPoint2D> brGenerous(consolidateDegenerates(bOpt, 0.1));
	drawTree(uniquer + makeString(targetServiceVolumes) + "_BrGenerous.png", brGenerous, width, height, border);
	cout << "\n Bifurcation2D tree: \t" << makeStringTreeOneLine(bOpt)
		<< "\n BranchPoint2D tree: \t" << makeStringTreeOneLine(br)
		<< "\n Generous BP2D tree: \t" << makeStringTreeOneLine(brGenerous)
		<< "\n note that BP2D trees are " << notIndeed(equivalentTrees(br, brGenerous, settleThresh)) << " the same"
		<< endl;
	cout << "\n br:" << endl;
	for(unsigned int i(0); i < br.size(); i++)
		cout << i << "\t" << br[i].tostring() << endl;
	cout << "\n brGenerous:" << endl;
	for(unsigned int i(0); i < brGenerous.size(); i++)
		cout << i << "\t" << brGenerous[i].tostring() << endl;
	delete[] bestBs;
	delete[] bestBMeasures;
}

void testExhaustiveConsolidated(){
	kisset(10394, 26468, 26630);
	string uniquer = "outputs/consol_c001_";
	int targetServiceVolumes(4), numBest(15), numDraw(1),
		maxFailedInsertions(100000), maxFailedNetworks(1000),
		width(1000), height(1000), border(50);
	double totalLengthCoefficient(0.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(1.0),
		nearestRadius(1.0),
		settleThresh(1.0e-3),
		consolidateThresh(0.01);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	vector<Bifurcation2D> b(initializeEvenRandomSpacingServiceVolumes(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty, targetServiceVolumes, maxFailedNetworks));
	//introduceBifurcationsFromSources(b);
	//settleBifurcations(b);
	//vector<Bifurcation2D> bOpt = optimizeByGlobalGrafting(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	vector<BranchPoint2D> *bestBs = new vector<BranchPoint2D>[numBest];
	double *bestBMeasures = new double[numBest];
	exhaustiveConsolidatedHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, true, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh);
	cout << "\n\n rank-1\t measure\t mi/m0" << endl;
	for(int i(0); i < numBest; i++){
		if(bestBMeasures[i] < 0.0)
			i = numBest;
		else{
			if(i < numDraw)
				drawTree(uniquer + makeString(targetServiceVolumes) + "_" + makeString(i) + ".png", bestBs[i], width, height, border);
			cout << i << "\t" << bestBMeasures[i] << "\t" << bestBMeasures[i]/bestBMeasures[0] << endl;
		}
	}
}

bool oneContainsTheOther(const vector<int> &s1, const vector<int> &s2){
	if(s1.size() < s2.size()){
		for(unsigned int i(0); i < s1.size(); i++){
			if(!containedIn(s1[i], s2))
				return false;
		}
		return true;
	}
	for(unsigned int i(0); i < s2.size(); i++){
		if(!containedIn(s2[i], s1))
			return false;
	}
	return true;
}

bool firstContainsSecond(const vector<int> &s1, const vector<int> &s2){
	for(unsigned int i(0); i < s2.size(); i++){
		if(!containedIn(s2[i], s1))
			return false;
	}
	return true;
}

double subtreeSimilarity(const vector<BranchPoint2D> &b1, const vector<BranchPoint2D> &b2){
	int maxCount1(0), maxCount2(0);
	vector<vector<int> > subtrees1, subtrees2;
	for(unsigned int i(0); i < b1.size(); i++){
		if(b1[i].childIndex.size() > 1) // exclude heart and tips, include higher order branch points
			subtrees1.push_back(tipSet(b1, i));
	}
	for(unsigned int i(0); i < b2.size(); i++){
		if(b2[i].childIndex.size() > 1) // exclude heart and tips, include higher order branch points
			subtrees2.push_back(tipSet(b2, i));
	}
	for(unsigned int i(0); i < subtrees1.size(); i++){
		for(unsigned int j(0); j < subtrees1.size(); j++){
			if(firstContainsSecond(subtrees1[i], subtrees1[j]))
				maxCount1++;
			if(firstContainsSecond(subtrees1[j], subtrees1[i]))
				maxCount1++;
		}
	}
	for(unsigned int i(0); i < subtrees2.size(); i++){
		for(unsigned int j(0); j < subtrees2.size(); j++){
			if(firstContainsSecond(subtrees2[i], subtrees2[j]))
				maxCount2++;
			if(firstContainsSecond(subtrees2[j], subtrees2[i]))
				maxCount2++;
		}
	}
	int maxCount(maxCount1);
	if(maxCount < maxCount2)
		maxCount = maxCount2;
	int count(0);
	for(unsigned int i(0); i < subtrees1.size(); i++){
		for(unsigned int j(0); j < subtrees2.size(); j++){
			if(firstContainsSecond(subtrees1[i], subtrees2[j]))
				count++;
			if(firstContainsSecond(subtrees2[j], subtrees1[i]))
				count++;
		}
	}
	return double(count)/double(maxCount);
}

void exhaustiveConsolidatedMeasures(){
	kisset(10394, 26468, 26630);
	string uniquer = "exhSim_";
	int numNetworks(1000), targetServiceVolumes(7), numBest(10395),
		maxFailedInsertions(100000), maxFailedNetworks(10000),
		width(300), height(300), border(width/40), numDraw(0),
		numAveBins(numBest);
	double totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0),
		nearestRadius(1.0),
		settleThresh(1.0e-3),
		consolidateThresh(0.01);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	vector<BranchPoint2D> **bestBs = new vector<BranchPoint2D>*[numNetworks];
	double **bestBMeasures = new double*[numNetworks];
	int previousProgress(0), numUpdates(10);
	for(int n(0); n < numNetworks; n++){
		vector<Bifurcation2D> b(initializeEvenRandomSpacingServiceVolumes(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty, targetServiceVolumes, maxFailedNetworks));
		bestBs[n] = new vector<BranchPoint2D>[numBest];
		bestBMeasures[n] = new double[numBest];
		exhaustiveConsolidatedHierarchyConfigurationSearch(b, numBest, bestBs[n], bestBMeasures[n], false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh);
		if((n + 1)*numUpdates/numNetworks > previousProgress){
			previousProgress++;
			double elapsedTime(difftime(time(NULL), startTime));
			cout << "\n Completed network " << n + 1 << " of " << numNetworks << " in " << niceTime(elapsedTime)
				<< "\n\t Time remaining: " << niceTime((elapsedTime/double(n + 1))*(numNetworks - n - 1)) << endl;
		}
	}
	int *maxRank = new int[numNetworks];
	for(int n(0); n < numNetworks; n++){
		for(int i(0); i < numBest; i++){
			if(bestBMeasures[n][i] >= 0.0)
				maxRank[n] = i + 1;
		}
	}

	double **sims = new double*[numNetworks];
	for(int n(0); n < numNetworks; n++){
		sims[n] = new double[numBest];
		for(int i(0); i < numBest; i++){
			if(bestBMeasures[n][i] >= 0.0)
				sims[n][i] = subtreeSimilarity(bestBs[n][i], bestBs[n][0]);
			else
				sims[n][i] = -1.0;
		}
	}

	cout << "\n\n RANKED BY MEASURE\n v-rr(0,1]\t interp_meas/meas[0] \t simToBest" << endl;
	cout << "0.0\t1.0\t1.0" << endl;
	for(int i(0); i < numAveBins; i++){
		double binCenter((i + 0.5)/numAveBins);
		double sum(0.0), simSum(0.0);
		for(int n(0); n < numNetworks; n++){
			int binBelow((int)floor(binCenter*maxRank[n]));
			//cout << " [ n" << n << " binBelow=" << binBelow << " maxRank[n]=" << maxRank[n] << " ] ";
			if(binBelow > maxRank[n] - 1)
				cout << "\n Bad bin!" << endl;
			else if(binBelow > maxRank[n] - 2){
				sum += bestBMeasures[n][binBelow]/bestBMeasures[n][0];
				simSum += sims[n][binBelow];
			}else{
				sum += linearInterpolation(double(binBelow)/maxRank[n], binCenter, double(binBelow + 1)/maxRank[n], bestBMeasures[n][binBelow], bestBMeasures[n][binBelow + 1])/bestBMeasures[n][0];
				double simBelow(sims[n][binBelow]), simNext(sims[n][binBelow + 1]);
				simSum += linearInterpolation(double(binBelow)/maxRank[n], binCenter, double(binBelow + 1)/maxRank[n], simBelow, simNext);
			}
		}
		cout << binCenter << "\t" << sum/numNetworks << "\t" << simSum/numNetworks << endl;
	}

	// re-order for decreasing similarity
	for(int n(0); n < numNetworks; n++){
		bool changed(true);
		while(changed){
			changed = false;
			for(int i(0); i < numBest - 1; i++){
				if(bestBMeasures[n][i] < 0.0)
					i = numBest;
				else{
					if(sims[n][i] < sims[n][i + 1]){
						changed = true;
						swap(bestBs[n][i], bestBs[n][i + 1]);
						swap(bestBMeasures[n][i], bestBMeasures[n][i + 1]);
						swap(sims[n][i], sims[n][i + 1]);
					}
				}
			}
		}
	}

	cout << "\n\n RANKED BY SIMILARITY\n v-rr(0,1]\t interp_meas/meas[0] \t simToBest" << endl;
	cout << "0.0\t1.0\t1.0" << endl;
	for(int i(0); i < numAveBins; i++){
		double binCenter((i + 0.5)/numAveBins);
		double sum(0.0), simSum(0.0);
		for(int n(0); n < numNetworks; n++){
			int binBelow((int)floor(binCenter*maxRank[n]));
			//cout << " [ n" << n << " binBelow=" << binBelow << " maxRank[n]=" << maxRank[n] << " ] ";
			if(binBelow > maxRank[n] - 1)
				cout << "\n Bad bin!" << endl;
			else if(binBelow > maxRank[n] - 2){
				sum += bestBMeasures[n][binBelow]/bestBMeasures[n][0];
				simSum += sims[n][binBelow];
			}else{
				sum += linearInterpolation(double(binBelow)/maxRank[n], binCenter, double(binBelow + 1)/maxRank[n], bestBMeasures[n][binBelow], bestBMeasures[n][binBelow + 1])/bestBMeasures[n][0];
				double simBelow(sims[n][binBelow]), simNext(sims[n][binBelow + 1]);
				simSum += linearInterpolation(double(binBelow)/maxRank[n], binCenter, double(binBelow + 1)/maxRank[n], simBelow, simNext);
			}
		}
		cout << binCenter << "\t" << sum/numNetworks << "\t" << simSum/numNetworks << endl;
	}

	vector<double> uniqueCount;
	for(int n(0); n < numNetworks; n++)
		uniqueCount.push_back(double(maxRank[n]));
	map<double, double> maxRankDist(probabilityDistribution(uniqueCount, maximum(numNetworks, maxRank) - minimum(numNetworks, maxRank) + 1));
	printMap(maxRankDist, "max_rank");


	/*
	cout << "\n\n v-rank/meausures->";
	for(int n(0); n < numNetworks; n++)
		cout << "\tn" << n;
	cout << endl;
	for(int i(0); i < numBest; i++){
		cout << i + 1;
		bool noValidNetworks(true);
		for(int n(0); n < numNetworks; n++){
			if(bestBMeasures[n][i] < 0.0)
				cout << "\t -";
			else{
				noValidNetworks = false;
				if(i < numDraw || maxRank[n] - i - 1 < numDraw ||
					abs(maxRank[n]/2 - i) < 2)
					drawTree(uniquer + makeString(targetServiceVolumes) + " _n" + makeString(n) + "_r" + makeString(i + 1) + ".png", bestBs[n][i], width, height, border);
				cout << "\t" << bestBMeasures[n][i];
			}
		}
		cout << endl;
		if(noValidNetworks)
			i = numBest;
	}
	cout << "\n\n";
	for(int n(0); n < numNetworks; n++)
		cout << "\trr" << n << "\trm" << n;
	cout << endl;
	for(int i(0); i < numBest; i++){
		bool noValidNetworks(true);
		for(int n(0); n < numNetworks; n++){
			cout << "\t" << double(i + 1)/maxRank[n];

			if(bestBMeasures[n][i] < 0.0)
				cout << "\t -";
			else{
				noValidNetworks = false;
				cout << "\t" << bestBMeasures[n][i]/bestBMeasures[n][0];
			}
		}
		cout << endl;
		if(noValidNetworks)
			i = numBest;
	}
	*/

	for(int n(0); n < numNetworks; n++){
		delete[] bestBs[n];
		delete[] bestBMeasures[n];
		delete[] sims[n];
	}
	delete[] bestBs;
	delete[] bestBMeasures;
	delete[] sims;
	delete[] maxRank;
}

void graftScaling(){
	//kisset (7446, 13195, 32375);
	//kisset (118, 23596, 28568);
	//kisset (31742, 27991, 19257);
	kisset (23777, 9548, 22058);
	string uniquer("graftScaling_");
	bool doDraw(true);
	double areaStart(10.0), areaInc(10.0), areaEnd(100),
		settleThresh(0.001), consolidateThresh(0.01),
		nearestRadius(1.0/sqrt(acos(-1.0))),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(9.0),
		maxPathLengthCoefficient(0.0);
	int numNetworks(10), numOpt(5), graftPasses(100), settleIts(10000),
		maxFailedInsertions(10000), width(1000), height(width), border(width/20);
	cout << "\n uniquer = " << uniquer
		<< "\n areaStart = " << areaStart << "\n areaInc = " << areaInc << "\n areaEnd = " << areaEnd
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<vector<double> > times;
	vector<vector<unsigned int> > terminalServiceVolumes;
	for(double area(areaStart); area < areaEnd + areaInc/2.0; area += areaInc){
		double side(sqrt(area/2.0)), heartx(side/2.0), hearty(nearestRadius/2.0);//side(sqrt(area))
		times.push_back(vector<double>());
		terminalServiceVolumes.push_back(vector<unsigned int>());
		for(int n(0); n < numNetworks; n++){
			vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, side, side, nearestRadius, heartx, hearty));
			//vector<Bifurcation2D> bCopy(b);
			//introduceBifurcationsFromSources(bCopy);
			//vector<BranchPoint2D> bGuess = consolidateDegenerates(bCopy, consolidateThresh);
			//double guessMeasure = tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			//cout << "\n\t\t\t\t\t\t Guess measure = " << guessMeasure << endl;
			//if(doDraw)
			//	drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_guess.png", bGuess, width, height, border);
			time_t tempStart(time(NULL));
			vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGrafting(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, graftPasses, settleIts);
			double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			times.back().push_back(difftime(time(NULL), tempStart));
			terminalServiceVolumes.back().push_back(b.size() - 1);
			cout << "\n\t\t Graft took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = " << graftMeasure << endl;
			if(doDraw)
				drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_graft.png", bGraft, width, height, border);
		}
	}
	cout << "\n\n area\t <TSV>\t +/-\t <time>\t +/-" << endl;
	for(unsigned int a(0); a < times.size(); a++){
		cout << areaStart + a*areaInc
			<< "\t" << mean(terminalServiceVolumes[a]) << "\t" << stDev(terminalServiceVolumes[a])
			<< "\t" << mean(times[a]) << "\t" << stDev(times[a])
			<< endl;
	}
}

void printReport(string name, int numBins, vector<vector<map<double, double> > > &lambdaL, vector<map<double, double> > &areaLambdaL,
	map<double, double> &allLambdaL, double areaStart, double areaInc, int numNetworks){
	cout << "\n\n" << name << "\t all_" << lambdaL.size() << "_areas/vols";
	for(unsigned int a(0); a < lambdaL.size(); a++)
		cout << "\t" << name << "\t area/vol=" << areaStart + a*areaInc;
	cout << endl;
	for(int i(0); i < numBins; i++){
		double k(nthMapKey(allLambdaL, i));
		cout << k << "\t" << allLambdaL[k];
		for(unsigned int a(0); a < lambdaL.size(); a++){
			k = nthMapKey(areaLambdaL[a], i);
			cout << "\t" << k << "\t" << areaLambdaL[a][k];
		}
		cout << endl;
	}

	for(unsigned int a(0); a < lambdaL.size(); a++){
		cout << "\n\n\t area/vol = " << areaStart + a*areaInc << endl;
		cout << name << "\t all_" << numNetworks*lambdaL.size() << "_networks\t" << name << "\t all_" << numNetworks << "_networks";
		for(int n(0); n < numNetworks; n++)
			cout << "\t" << name << "\t n" << n;
		for(int i(0); i < numBins; i++){
			double k(nthMapKey(allLambdaL, i));
			cout << k << "\t" << allLambdaL[k];
			k = nthMapKey(areaLambdaL[a], i);
			cout << "\t" << k << "\t" << areaLambdaL[a][k];
			for(int n(0); n < numNetworks; n++){
				k = nthMapKey(lambdaL[a][n], i);
				cout << "\t" << k << "\t" << lambdaL[a][n][k];
			}
			cout << endl;
		}
	}
}

void lambdaLGammaOpt(){
	kisset (7446, 13195, 32375);
	string uniquer("lambdaLGammaOptOneNetSquare_c100_");// "OneFour" for oneFour, "Square" for square
	bool doDraw(true);
	double areaStart(50.0), areaInc(areaStart), areaEnd(areaStart),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0),
		localRad(2.0*nearestRadius);
	int numNetworks(100), numOpt(10), graftPasses(100), settleIts(10000), numBins(50),
		maxFailedInsertions(10000), width(1000), height(1000), border(10);// width 250 for oneFour
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		<< "\n note that area is really halved"
		<< "\n areaStart = " << areaStart << "\n areaInc = " << areaInc << "\n areaEnd = " << areaEnd
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n localRad = " << localRad
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses << "\n numBins = " << numBins
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<vector<double> > times, measures;
	vector<vector<vector<double> > >lambdaLs, gammas, tipRatios;
	vector<vector<unsigned int> > terminalServiceVolumes;
	vector<vector<int> > ignoreCounts, ignoreCountsGamma;
	for(double area(areaStart); area < areaEnd + areaInc/2.0; area += areaInc){
		double side(sqrt(area/2.0)), heartx(0.5*side), hearty(nearestRadius/2.0);//heartx 0.25*side for oneFour
		times.push_back(vector<double>());
		measures.push_back(vector<double>());
		lambdaLs.push_back(vector<vector<double> >());
		gammas.push_back(vector<vector<double> >());
		tipRatios.push_back(vector<vector<double> >());
		terminalServiceVolumes.push_back(vector<unsigned int>());
		ignoreCounts.push_back(vector<int>());
		ignoreCountsGamma.push_back(vector<int>());
		for(int n(0); n < numNetworks; n++){ // maxFailedInsertions, 0.5*side, 2.0*side... for oneFour
			vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, 1.0*side, 1.0*side, nearestRadius, heartx, hearty));
			vector<Bifurcation2D> bCopy(b);
			introduceBifurcationsFromSources(bCopy);
			vector<BranchPoint2D> bGuess = consolidateDegenerates(bCopy, consolidateThresh);
			//double guessMeasure = tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			//cout << "\n\t\t\t\t\t\t Guess measure = " << guessMeasure << endl;
			if(doDraw)
				drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_guess.png", bGuess, width, height, border);
			time_t tempStart(time(NULL));
			vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingle(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, graftPasses, settleIts);
			double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			times.back().push_back(difftime(time(NULL), tempStart));
			measures.back().push_back(graftMeasure);
			terminalServiceVolumes.back().push_back(b.size() - 1);
			cout << "\n\t\t Graft took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = " << graftMeasure << endl;
			if(doDraw)
				drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_graft.png", bGraft, width, height, border);
			ignoreCounts.back().push_back(0);
			ignoreCountsGamma.back().push_back(0);
			lambdaLs.back().push_back(validBranchingLengthRatios(bGraft, ignoreCounts.back().back(), lengthThresh));
			gammas.back().push_back(validParentChildLengthRatios(bGraft, ignoreCountsGamma.back().back(), lengthThresh));
			tipRatios.back().push_back(validTipCountRatios(bGraft, ignoreCounts.back().back(), lengthThresh)); // ignoreCounts should be the same as for lambdaL
		}
	}
	cout << "\n\n area\t <TSV>\t +/-\t fracIgn_rat\t +/-\t fracIgn_gam\t +/-\t measure\t +/-\t <time>\t +/-" << endl;
	for(unsigned int a(0); a < times.size(); a++){
		vector<double> ignoreFracs, ignoreFracsGamma;
		for(int n(0); n < numNetworks; n++){
			ignoreFracs.push_back(double(ignoreCounts[a][n])/double(terminalServiceVolumes[a][n] - 1));
			ignoreFracsGamma.push_back(double(ignoreCountsGamma[a][n])/double(2*terminalServiceVolumes[a][n] - 2));
		}
		cout << areaStart + a*areaInc
			<< "\t" << mean(terminalServiceVolumes[a]) << "\t" << stDev(terminalServiceVolumes[a])
			<< "\t" << mean(ignoreFracs) << "\t" << stDev(ignoreFracs)
			<< "\t" << mean(ignoreFracsGamma) << "\t" << stDev(ignoreFracsGamma)
			<< "\t" << mean(measures[a]) << "\t" << stDev(measures[a])
			<< "\t" << mean(times[a]) << "\t" << stDev(times[a])
			<< endl;
	}

	vector<vector<map<double, double> > > lambdaL, gamma, tipRatio; // single networks
	vector<map<double, double> > areaLambdaL, areaGamma, areaTipRatio; // within area, over networks
	vector<double> allLambdaLs, allGammas, allTipRatios; // over all areas
	for(unsigned int a(0); a < lambdaLs.size(); a++){
		lambdaL.push_back(vector<map<double, double> >());
		gamma.push_back(vector<map<double, double> >());
		tipRatio.push_back(vector<map<double, double> >());
		vector<double> areaLambdaLs, areaGammas, areaTipRatios; // within area, over networks
		for(int n(0); n < numNetworks; n++){
			areaLambdaLs.insert(areaLambdaLs.end(), lambdaLs[a][n].begin(), lambdaLs[a][n].end());
			areaGammas.insert(areaGammas.end(), gammas[a][n].begin(), gammas[a][n].end());
			areaTipRatios.insert(areaTipRatios.end(), tipRatios[a][n].begin(), tipRatios[a][n].end());
			lambdaL.back().push_back(probabilityDistribution(lambdaLs[a][n], numBins));
			gamma.back().push_back(probabilityDistribution(gammas[a][n], numBins));
			tipRatio.back().push_back(probabilityDistribution(tipRatios[a][n], numBins));
		}
		areaLambdaL.push_back(probabilityDistribution(areaLambdaLs, numBins));
		areaGamma.push_back(probabilityDistribution(areaGammas, numBins));
		areaTipRatio.push_back(probabilityDistribution(areaTipRatios, numBins));
		allLambdaLs.insert(allLambdaLs.end(), areaLambdaLs.begin(), areaLambdaLs.end());
		allGammas.insert(allGammas.end(), areaGammas.begin(), areaGammas.end());
		allTipRatios.insert(allTipRatios.end(), areaTipRatios.begin(), areaTipRatios.end());
	}
	map<double, double> allLambdaL(probabilityDistribution(allLambdaLs, numBins)),
		allGamma(probabilityDistribution(allGammas, numBins)),
		allTipRatio(probabilityDistribution(allTipRatios, numBins));

	printReport("lambdaL", numBins, lambdaL, areaLambdaL, allLambdaL, areaStart, areaInc, numNetworks);
	printReport("gamma", numBins, gamma, areaGamma, allGamma, areaStart, areaInc, numNetworks);
	printReport("tipRatio", numBins, tipRatio, areaTipRatio, allTipRatio, areaStart, areaInc, numNetworks);

	
}

void lambdaLGammaOptTips(){
	kisset (7446, 13195, 32375);
	string uniquer("lambdaLGammaOptOneNetSquare_c100_");// "OneFour" for oneFour, "Square" for square
	bool doDraw(true);
	double areaStart(101.0), areaInc(areaStart), areaEnd(areaStart),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0),
		localRad(2.0*nearestRadius),
		tipRatioThresh(0.01);
	int numNetworks(1), numOpt(1), graftPasses(100), settleIts(10000), numBins(50),
		maxFailedInsertions(10000), width(1000), height(1000), border(10);// width 250 for oneFour
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		<< "\n note that area is really halved"
		<< "\n areaStart = " << areaStart << "\n areaInc = " << areaInc << "\n areaEnd = " << areaEnd
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n localRad = " << localRad << "\n tipRatioThresh = " << tipRatioThresh
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses << "\n numBins = " << numBins
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<vector<double> > times, measures;
	vector<vector<vector<double> > >lambdaLs, gammas, tipRatios;
	vector<vector<unsigned int> > terminalServiceVolumes;
	vector<vector<int> > ignoreCounts, ignoreCountsGamma;
	for(double area(areaStart); area < areaEnd + areaInc/2.0; area += areaInc){
		double side(sqrt(area/2.0)), heartx(0.5*side), hearty(nearestRadius/2.0);//heartx 0.25*side for oneFour
		times.push_back(vector<double>());
		measures.push_back(vector<double>());
		lambdaLs.push_back(vector<vector<double> >());
		gammas.push_back(vector<vector<double> >());
		tipRatios.push_back(vector<vector<double> >());
		terminalServiceVolumes.push_back(vector<unsigned int>());
		ignoreCounts.push_back(vector<int>());
		ignoreCountsGamma.push_back(vector<int>());
		for(int n(0); n < numNetworks; n++){ // maxFailedInsertions, 0.5*side, 2.0*side... for oneFour
			vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, 1.0*side, 1.0*side, nearestRadius, heartx, hearty));
			vector<Bifurcation2D> bCopy(b);
			introduceBifurcationsFromSources(bCopy);
			vector<BranchPoint2D> bGuess = consolidateDegenerates(bCopy, consolidateThresh);
			cout << "\n bCopy has minimumTipRatio = " << minimumTipRatio(bCopy) 
				<< "\n bGuess has minimumTipRatio = " << minimumTipRatio(bGuess) << endl;
			//double guessMeasure = tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			//cout << "\n\t\t\t\t\t\t Guess measure = " << guessMeasure << endl;
			if(doDraw)
				drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_guess.png", bGuess, width, height, border);
			time_t tempStart(time(NULL));
			vector<Bifurcation2D> bBi;
			vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingleAndTipRatioThresh(numOpt, b, bBi, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, tipRatioThresh, graftPasses, settleIts);
			double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			times.back().push_back(difftime(time(NULL), tempStart));
			measures.back().push_back(graftMeasure);
			terminalServiceVolumes.back().push_back(b.size() - 1);
			cout << "\n\t\t Graft took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = "
				<< graftMeasure << "\n\t\t minimumTipRatio(bBi) = " << minimumTipRatio(bBi)
				<< "\n\t\t minimumTipRatio(bGraft) = " << minimumTipRatio(bGraft) << endl;
			if(doDraw)
				drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_graft.png", bGraft, width, height, border);
			ignoreCounts.back().push_back(0);
			ignoreCountsGamma.back().push_back(0);
			lambdaLs.back().push_back(validBranchingLengthRatios(bGraft, ignoreCounts.back().back(), lengthThresh));
			gammas.back().push_back(validParentChildLengthRatios(bGraft, ignoreCountsGamma.back().back(), lengthThresh));
			tipRatios.back().push_back(validTipCountRatios(bGraft, ignoreCounts.back().back(), lengthThresh)); // ignoreCounts should be the same as for lambdaL
		}
	}
	cout << "\n\n area\t <TSV>\t +/-\t fracIgn_rat\t +/-\t fracIgn_gam\t +/-\t measure\t +/-\t <time>\t +/-" << endl;
	for(unsigned int a(0); a < times.size(); a++){
		vector<double> ignoreFracs, ignoreFracsGamma;
		for(int n(0); n < numNetworks; n++){
			ignoreFracs.push_back(double(ignoreCounts[a][n])/double(terminalServiceVolumes[a][n] - 1));
			ignoreFracsGamma.push_back(double(ignoreCountsGamma[a][n])/double(2*terminalServiceVolumes[a][n] - 2));
		}
		cout << areaStart + a*areaInc
			<< "\t" << mean(terminalServiceVolumes[a]) << "\t" << stDev(terminalServiceVolumes[a])
			<< "\t" << mean(ignoreFracs) << "\t" << stDev(ignoreFracs)
			<< "\t" << mean(ignoreFracsGamma) << "\t" << stDev(ignoreFracsGamma)
			<< "\t" << mean(measures[a]) << "\t" << stDev(measures[a])
			<< "\t" << mean(times[a]) << "\t" << stDev(times[a])
			<< endl;
	}

	vector<vector<map<double, double> > > lambdaL, gamma, tipRatio; // single networks
	vector<map<double, double> > areaLambdaL, areaGamma, areaTipRatio; // within area, over networks
	vector<double> allLambdaLs, allGammas, allTipRatios; // over all areas
	for(unsigned int a(0); a < lambdaLs.size(); a++){
		lambdaL.push_back(vector<map<double, double> >());
		gamma.push_back(vector<map<double, double> >());
		tipRatio.push_back(vector<map<double, double> >());
		vector<double> areaLambdaLs, areaGammas, areaTipRatios; // within area, over networks
		for(int n(0); n < numNetworks; n++){
			areaLambdaLs.insert(areaLambdaLs.end(), lambdaLs[a][n].begin(), lambdaLs[a][n].end());
			areaGammas.insert(areaGammas.end(), gammas[a][n].begin(), gammas[a][n].end());
			areaTipRatios.insert(areaTipRatios.end(), tipRatios[a][n].begin(), tipRatios[a][n].end());
			lambdaL.back().push_back(probabilityDistribution(lambdaLs[a][n], numBins));
			gamma.back().push_back(probabilityDistribution(gammas[a][n], numBins));
			tipRatio.back().push_back(probabilityDistribution(tipRatios[a][n], numBins));
		}
		areaLambdaL.push_back(probabilityDistribution(areaLambdaLs, numBins));
		areaGamma.push_back(probabilityDistribution(areaGammas, numBins));
		areaTipRatio.push_back(probabilityDistribution(areaTipRatios, numBins));
		allLambdaLs.insert(allLambdaLs.end(), areaLambdaLs.begin(), areaLambdaLs.end());
		allGammas.insert(allGammas.end(), areaGammas.begin(), areaGammas.end());
		allTipRatios.insert(allTipRatios.end(), areaTipRatios.begin(), areaTipRatios.end());
	}
	map<double, double> allLambdaL(probabilityDistribution(allLambdaLs, numBins)),
		allGamma(probabilityDistribution(allGammas, numBins)),
		allTipRatio(probabilityDistribution(allTipRatios, numBins));

	printReport("lambdaL", numBins, lambdaL, areaLambdaL, allLambdaL, areaStart, areaInc, numNetworks);
	printReport("gamma", numBins, gamma, areaGamma, allGamma, areaStart, areaInc, numNetworks);
	printReport("tipRatio", numBins, tipRatio, areaTipRatio, allTipRatio, areaStart, areaInc, numNetworks);

	
}

void lambdaLGammaOptNoSearch(){
	kisset (7446, 13195, 32375);
	string uniquer("lambdaLGammaOptNoSearchOneFour_");// OneFour for oneFour
	bool doDraw(true);
	double areaStart(400.0), areaInc(400.0), areaEnd(4000.0),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),
		localRad(2.0*nearestRadius);
	int numNetworks(100), numOpt(10), graftPasses(100), settleIts(10000), numBins(50),
		maxFailedInsertions(10000), width(250), height(1000), border(height/40); // 250 for OneFour
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		<< "\n note that area is really halved"
		<< "\n areaStart = " << areaStart << "\n areaInc = " << areaInc << "\n areaEnd = " << areaEnd
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n localRad = " << localRad
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses << "\n numBins = " << numBins
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<vector<double> > times, measures100, measures010, measures001;
	vector<vector<vector<double> > >lambdaLs, gammas, tipRatios;
	vector<vector<unsigned int> > terminalServiceVolumes;
	vector<vector<int> > ignoreCounts, ignoreCountsGamma;
	for(double area(areaStart); area < areaEnd + areaInc/2.0; area += areaInc){
		double side(sqrt(area/2.0)), heartx(side/4.0), hearty(nearestRadius/2.0);//side(sqrt(area)) // side/4.0 for oneFour
		times.push_back(vector<double>());
		measures100.push_back(vector<double>());
		measures010.push_back(vector<double>());
		measures001.push_back(vector<double>());
		lambdaLs.push_back(vector<vector<double> >());
		gammas.push_back(vector<vector<double> >());
		tipRatios.push_back(vector<vector<double> >());
		terminalServiceVolumes.push_back(vector<unsigned int>());
		ignoreCounts.push_back(vector<int>());
		ignoreCountsGamma.push_back(vector<int>());
		for(int n(0); n < numNetworks; n++){
			time_t tempStart(time(NULL)); // 0.5*side, 2.0*side for oneFour
			vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, 0.5*side, 2.0*side, nearestRadius, heartx, hearty));//initializeEvenRandomSpacingCircle(maxFailedInsertions, side, nearestRadius));//
			introduceBifurcationsFromSources(b);
			//vector<Bifurcation2D> bCopy(b);
			//introduceBifurcationsFromSources(bCopy);
			//vector<BranchPoint2D> bGuess = consolidateDegenerates(bCopy, consolidateThresh);
			//double guessMeasure = tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			//cout << "\n\t\t\t\t\t\t Guess measure = " << guessMeasure << endl;
			//if(doDraw)
			//	drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_guess.png", bGuess, width, height, border);
			//vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingle(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, graftPasses, settleIts);
			vector<BranchPoint2D> bCons(consolidateDegenerates(b));
			times.back().push_back(difftime(time(NULL), tempStart));
			measures100.back().push_back(tripletMeasure(bCons, 1.0, 0.0, 0.0));
			measures010.back().push_back(tripletMeasure(bCons, 0.0, 1.0, 0.0));
			measures001.back().push_back(tripletMeasure(bCons, 0.0, 0.0, 1.0));
			terminalServiceVolumes.back().push_back(b.size() - 1);
			cout << "\n\t\t Graft took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs" << endl;
			if(doDraw)
				drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_graft.png", bCons, width, height, border);
			ignoreCounts.back().push_back(0);
			ignoreCountsGamma.back().push_back(0);
			lambdaLs.back().push_back(validBranchingLengthRatios(bCons, ignoreCounts.back().back(), lengthThresh));
			gammas.back().push_back(validParentChildLengthRatios(bCons, ignoreCountsGamma.back().back(), lengthThresh));
			tipRatios.back().push_back(validTipCountRatios(bCons, ignoreCounts.back().back(), lengthThresh)); // ignoreCounts should be the same as for lambdaL
		}
	}
	cout << "\n\n area\t <TSV>\t +/-\t fracIgn_rat\t +/-\t fracIgn_gam\t +/-\t measure100\t +/-\t measure010\t +/-\t measure001\t +/-\t <time>\t +/-" << endl;
	for(unsigned int a(0); a < times.size(); a++){
		vector<double> ignoreFracs, ignoreFracsGamma;
		for(int n(0); n < numNetworks; n++){
			ignoreFracs.push_back(double(ignoreCounts[a][n])/double(terminalServiceVolumes[a][n] - 1));
			ignoreFracsGamma.push_back(double(ignoreCountsGamma[a][n])/double(2*terminalServiceVolumes[a][n] - 2));
		}
		cout << areaStart + a*areaInc
			<< "\t" << mean(terminalServiceVolumes[a]) << "\t" << stDev(terminalServiceVolumes[a])
			<< "\t" << mean(ignoreFracs) << "\t" << stDev(ignoreFracs)
			<< "\t" << mean(ignoreFracsGamma) << "\t" << stDev(ignoreFracsGamma)
			<< "\t" << mean(measures100[a]) << "\t" << stDev(measures100[a])
			<< "\t" << mean(measures010[a]) << "\t" << stDev(measures010[a])
			<< "\t" << mean(measures001[a]) << "\t" << stDev(measures001[a])
			<< "\t" << mean(times[a]) << "\t" << stDev(times[a])
			<< endl;
	}

	vector<vector<map<double, double> > > lambdaL, gamma, tipRatio; // single networks
	vector<map<double, double> > areaLambdaL, areaGamma, areaTipRatio; // within area, over networks
	vector<double> allLambdaLs, allGammas, allTipRatios; // over all areas
	for(unsigned int a(0); a < lambdaLs.size(); a++){
		lambdaL.push_back(vector<map<double, double> >());
		gamma.push_back(vector<map<double, double> >());
		tipRatio.push_back(vector<map<double, double> >());
		vector<double> areaLambdaLs, areaGammas, areaTipRatios; // within area, over networks
		for(int n(0); n < numNetworks; n++){
			areaLambdaLs.insert(areaLambdaLs.end(), lambdaLs[a][n].begin(), lambdaLs[a][n].end());
			areaGammas.insert(areaGammas.end(), gammas[a][n].begin(), gammas[a][n].end());
			areaTipRatios.insert(areaTipRatios.end(), tipRatios[a][n].begin(), tipRatios[a][n].end());
			lambdaL.back().push_back(probabilityDistribution(lambdaLs[a][n], numBins));
			gamma.back().push_back(probabilityDistribution(gammas[a][n], numBins));
			tipRatio.back().push_back(probabilityDistribution(tipRatios[a][n], numBins));
		}
		areaLambdaL.push_back(probabilityDistribution(areaLambdaLs, numBins));
		areaGamma.push_back(probabilityDistribution(areaGammas, numBins));
		areaTipRatio.push_back(probabilityDistribution(areaTipRatios, numBins));
		allLambdaLs.insert(allLambdaLs.end(), areaLambdaLs.begin(), areaLambdaLs.end());
		allGammas.insert(allGammas.end(), areaGammas.begin(), areaGammas.end());
		allTipRatios.insert(allTipRatios.end(), areaTipRatios.begin(), areaTipRatios.end());
	}
	map<double, double> allLambdaL(probabilityDistribution(allLambdaLs, numBins)),
		allGamma(probabilityDistribution(allGammas, numBins)),
		allTipRatio(probabilityDistribution(allTipRatios, numBins));

	printReport("lambdaL", numBins, lambdaL, areaLambdaL, allLambdaL, areaStart, areaInc, numNetworks);
	printReport("gamma", numBins, gamma, areaGamma, allGamma, areaStart, areaInc, numNetworks);
	printReport("tipRatio", numBins, tipRatio, areaTipRatio, allTipRatio, areaStart, areaInc, numNetworks);

	
}

void lambdaL3D(){
	kisset (7446, 13195, 32375);
	string uniquer("lambdaL3D_c195_");
	bool doDraw(true);
	double sideStart(2.0), sideInc(1.0), sideEnd(8.0),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(9.0),
		maxPathLengthCoefficient(5.0),
		localRad(1.0*nearestRadius);
	int numNetworks(100), numOpt(5), graftPasses(200), settleIts(10000), numBins(27),
		maxFailedInsertions(10000), width(1000), height(1000), border(width/40);
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		<< "\n note that area is really halved"
		<< "\n sideStart = " << sideStart << "\n sideInc = " << sideInc << "\n sideEnd = " << sideEnd
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n localRad = " << localRad
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses << "\n numBins = " << numBins
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<vector<double> > times, measures;
	vector<vector<vector<double> > >lambdaLs, gammas, tipRatios;
	vector<vector<unsigned int> > terminalServiceVolumes;
	vector<vector<int> > ignoreCounts, ignoreCountsGamma;
	for(double side(sideStart); side < sideEnd + sideInc/2.0; side += sideInc){
		double heartx(side/2.0), hearty(side/2.0), heartz(side/2.0);//side(sqrt(area))
		times.push_back(vector<double>());
		measures.push_back(vector<double>());
		lambdaLs.push_back(vector<vector<double> >());
		gammas.push_back(vector<vector<double> >());
		tipRatios.push_back(vector<vector<double> >());
		terminalServiceVolumes.push_back(vector<unsigned int>());
		ignoreCounts.push_back(vector<int>());
		ignoreCountsGamma.push_back(vector<int>());
		for(int n(0); n < numNetworks; n++){
			vector<Bifurcation3D> b(initializeEvenRandomSpacing(maxFailedInsertions, side, side, side, nearestRadius, heartx, hearty, heartz));
			//vector<Bifurcation2D> bCopy(b);
			//introduceBifurcationsFromSources(bCopy);
			//vector<BranchPoint2D> bGuess = consolidateDegenerates(bCopy, consolidateThresh);
			//double guessMeasure = tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			//cout << "\n\t\t\t\t\t\t Guess measure = " << guessMeasure << endl;
			//if(doDraw)
			//	drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_guess.png", bGuess, width, height, border);
			time_t tempStart(time(NULL));
			vector<BranchPoint3D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingle(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, graftPasses, settleIts);
			double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			times.back().push_back(difftime(time(NULL), tempStart));
			measures.back().push_back(graftMeasure);
			terminalServiceVolumes.back().push_back(b.size() - 1);
			cout << "\n\t\t Graft took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = " << graftMeasure << endl;
			if(doDraw)
				drawTree(uniquer + makeString((int)ceil(side)) + "_n" + makeString(n) + "_graft.png", bGraft, width, height, border);
			ignoreCounts.back().push_back(0);
			ignoreCountsGamma.back().push_back(0);
			lambdaLs.back().push_back(validBranchingLengthRatios(bGraft, ignoreCounts.back().back(), lengthThresh));
			gammas.back().push_back(validParentChildLengthRatios(bGraft, ignoreCountsGamma.back().back(), lengthThresh));
			tipRatios.back().push_back(validTipCountRatios(bGraft, ignoreCounts.back().back(), lengthThresh)); // ignoreCounts should be the same as for lambdaL
		}
	}
	cout << "\n\n area\t <TSV>\t +/-\t fracIgn_rat\t +/-\t fracIgn_gam\t +/-\t measure\t +/-\t <time>\t +/-" << endl;
	for(unsigned int a(0); a < times.size(); a++){
		vector<double> ignoreFracs, ignoreFracsGamma;
		for(int n(0); n < numNetworks; n++){
			ignoreFracs.push_back(double(ignoreCounts[a][n])/double(terminalServiceVolumes[a][n] - 1));
			ignoreFracsGamma.push_back(double(ignoreCountsGamma[a][n])/double(2*terminalServiceVolumes[a][n] - 2));
		}
		cout << sideStart + a*sideInc
			<< "\t" << mean(terminalServiceVolumes[a]) << "\t" << stDev(terminalServiceVolumes[a])
			<< "\t" << mean(ignoreFracs) << "\t" << stDev(ignoreFracs)
			<< "\t" << mean(ignoreFracsGamma) << "\t" << stDev(ignoreFracsGamma)
			<< "\t" << mean(measures[a]) << "\t" << stDev(measures[a])
			<< "\t" << mean(times[a]) << "\t" << stDev(times[a])
			<< endl;
	}

	vector<vector<map<double, double> > > lambdaL, gamma, tipRatio; // single networks
	vector<map<double, double> > areaLambdaL, areaGamma, areaTipRatio; // within area, over networks
	vector<double> allLambdaLs, allGammas, allTipRatios; // over all areas
	for(unsigned int a(0); a < lambdaLs.size(); a++){
		lambdaL.push_back(vector<map<double, double> >());
		gamma.push_back(vector<map<double, double> >());
		tipRatio.push_back(vector<map<double, double> >());
		vector<double> areaLambdaLs, areaGammas, areaTipRatios; // within area, over networks
		for(int n(0); n < numNetworks; n++){
			areaLambdaLs.insert(areaLambdaLs.end(), lambdaLs[a][n].begin(), lambdaLs[a][n].end());
			areaGammas.insert(areaGammas.end(), gammas[a][n].begin(), gammas[a][n].end());
			areaTipRatios.insert(areaTipRatios.end(), tipRatios[a][n].begin(), tipRatios[a][n].end());
			lambdaL.back().push_back(probabilityDistribution(lambdaLs[a][n], numBins));
			gamma.back().push_back(probabilityDistribution(gammas[a][n], numBins));
			tipRatio.back().push_back(probabilityDistribution(tipRatios[a][n], numBins));
		}
		areaLambdaL.push_back(probabilityDistribution(areaLambdaLs, numBins));
		areaGamma.push_back(probabilityDistribution(areaGammas, numBins));
		areaTipRatio.push_back(probabilityDistribution(areaTipRatios, numBins));
		allLambdaLs.insert(allLambdaLs.end(), areaLambdaLs.begin(), areaLambdaLs.end());
		allGammas.insert(allGammas.end(), areaGammas.begin(), areaGammas.end());
		allTipRatios.insert(allTipRatios.end(), areaTipRatios.begin(), areaTipRatios.end());
	}
	map<double, double> allLambdaL(probabilityDistribution(allLambdaLs, numBins)),
		allGamma(probabilityDistribution(allGammas, numBins)),
		allTipRatio(probabilityDistribution(allTipRatios, numBins));

	printReport("lambdaL", numBins, lambdaL, areaLambdaL, allLambdaL, sideStart, sideInc, numNetworks);
	printReport("gamma", numBins, gamma, areaGamma, allGamma, sideStart, sideInc, numNetworks);
	printReport("tipRatio", numBins, tipRatio, areaTipRatio, allTipRatio, sideStart, sideInc, numNetworks);

	
}

void lambdaLSphere(){
	kisset (7446, 13195, 32375);
	string uniquer("outputs/lambdaLSphere_c195_");
	bool doDraw(true);
	double sideStart(3.0), sideInc(1.0), sideEnd(4.0),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(9.0),
		maxPathLengthCoefficient(5.0),
		localRad(1.0*nearestRadius);
	int numNetworks(2), numOpt(5), graftPasses(200), settleIts(10000), numBins(27),
		maxFailedInsertions(10000), width(1000), height(1000), border(width/40);
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		<< "\n sideStart = " << sideStart << "\n sideInc = " << sideInc << "\n sideEnd = " << sideEnd
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n localRad = " << localRad
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses << "\n numBins = " << numBins
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<vector<double> > times, measures;
	vector<vector<vector<double> > >lambdaLs, gammas, tipRatios;
	vector<vector<unsigned int> > terminalServiceVolumes;
	vector<vector<int> > ignoreCounts, ignoreCountsGamma;
	for(double side(sideStart); side < sideEnd + sideInc/2.0; side += sideInc){
		double heartx(side/2.0), hearty(side/2.0), heartz(side/2.0);//side(sqrt(area))
		times.push_back(vector<double>());
		measures.push_back(vector<double>());
		lambdaLs.push_back(vector<vector<double> >());
		gammas.push_back(vector<vector<double> >());
		tipRatios.push_back(vector<vector<double> >());
		terminalServiceVolumes.push_back(vector<unsigned int>());
		ignoreCounts.push_back(vector<int>());
		ignoreCountsGamma.push_back(vector<int>());
		for(int n(0); n < numNetworks; n++){
			vector<Bifurcation3D> b(initializeEvenRandomSpacingSphere(maxFailedInsertions, side/2.0, nearestRadius));
			//vector<Bifurcation2D> bCopy(b);
			//introduceBifurcationsFromSources(bCopy);
			//vector<BranchPoint2D> bGuess = consolidateDegenerates(bCopy, consolidateThresh);
			//double guessMeasure = tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			//cout << "\n\t\t\t\t\t\t Guess measure = " << guessMeasure << endl;
			//if(doDraw)
			//	drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_guess.png", bGuess, width, height, border);
			time_t tempStart(time(NULL));
			vector<BranchPoint3D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingle(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, graftPasses, settleIts);
			double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			times.back().push_back(difftime(time(NULL), tempStart));
			measures.back().push_back(graftMeasure);
			terminalServiceVolumes.back().push_back(b.size() - 1);
			cout << "\n\t\t Graft took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = " << graftMeasure << endl;
			if(doDraw)
				drawTree(uniquer + makeString((int)ceil(side)) + "_n" + makeString(n) + "_graft.png", bGraft, width, height, border);
			ignoreCounts.back().push_back(0);
			ignoreCountsGamma.back().push_back(0);
			lambdaLs.back().push_back(validBranchingLengthRatios(bGraft, ignoreCounts.back().back(), lengthThresh));
			gammas.back().push_back(validParentChildLengthRatios(bGraft, ignoreCountsGamma.back().back(), lengthThresh));
			tipRatios.back().push_back(validTipCountRatios(bGraft, ignoreCounts.back().back(), lengthThresh)); // ignoreCounts should be the same as for lambdaL
		}
	}
	cout << "\n\n area\t <TSV>\t +/-\t fracIgn_rat\t +/-\t fracIgn_gam\t +/-\t measure\t +/-\t <time>\t +/-" << endl;
	for(unsigned int a(0); a < times.size(); a++){
		vector<double> ignoreFracs, ignoreFracsGamma;
		for(int n(0); n < numNetworks; n++){
			ignoreFracs.push_back(double(ignoreCounts[a][n])/double(terminalServiceVolumes[a][n] - 1));
			ignoreFracsGamma.push_back(double(ignoreCountsGamma[a][n])/double(2*terminalServiceVolumes[a][n] - 2));
		}
		cout << sideStart + a*sideInc
			<< "\t" << mean(terminalServiceVolumes[a]) << "\t" << stDev(terminalServiceVolumes[a])
			<< "\t" << mean(ignoreFracs) << "\t" << stDev(ignoreFracs)
			<< "\t" << mean(ignoreFracsGamma) << "\t" << stDev(ignoreFracsGamma)
			<< "\t" << mean(measures[a]) << "\t" << stDev(measures[a])
			<< "\t" << mean(times[a]) << "\t" << stDev(times[a])
			<< endl;
	}

	vector<vector<map<double, double> > > lambdaL, gamma, tipRatio; // single networks
	vector<map<double, double> > areaLambdaL, areaGamma, areaTipRatio; // within area, over networks
	vector<double> allLambdaLs, allGammas, allTipRatios; // over all areas
	for(unsigned int a(0); a < lambdaLs.size(); a++){
		lambdaL.push_back(vector<map<double, double> >());
		gamma.push_back(vector<map<double, double> >());
		tipRatio.push_back(vector<map<double, double> >());
		vector<double> areaLambdaLs, areaGammas, areaTipRatios; // within area, over networks
		for(int n(0); n < numNetworks; n++){
			areaLambdaLs.insert(areaLambdaLs.end(), lambdaLs[a][n].begin(), lambdaLs[a][n].end());
			areaGammas.insert(areaGammas.end(), gammas[a][n].begin(), gammas[a][n].end());
			areaTipRatios.insert(areaTipRatios.end(), tipRatios[a][n].begin(), tipRatios[a][n].end());
			lambdaL.back().push_back(probabilityDistribution(lambdaLs[a][n], numBins));
			gamma.back().push_back(probabilityDistribution(gammas[a][n], numBins));
			tipRatio.back().push_back(probabilityDistribution(tipRatios[a][n], numBins));
		}
		areaLambdaL.push_back(probabilityDistribution(areaLambdaLs, numBins));
		areaGamma.push_back(probabilityDistribution(areaGammas, numBins));
		areaTipRatio.push_back(probabilityDistribution(areaTipRatios, numBins));
		allLambdaLs.insert(allLambdaLs.end(), areaLambdaLs.begin(), areaLambdaLs.end());
		allGammas.insert(allGammas.end(), areaGammas.begin(), areaGammas.end());
		allTipRatios.insert(allTipRatios.end(), areaTipRatios.begin(), areaTipRatios.end());
	}
	map<double, double> allLambdaL(probabilityDistribution(allLambdaLs, numBins)),
		allGamma(probabilityDistribution(allGammas, numBins)),
		allTipRatio(probabilityDistribution(allTipRatios, numBins));

	printReport("lambdaL", numBins, lambdaL, areaLambdaL, allLambdaL, sideStart, sideInc, numNetworks);
	printReport("gamma", numBins, gamma, areaGamma, allGamma, sideStart, sideInc, numNetworks);
	printReport("tipRatio", numBins, tipRatio, areaTipRatio, allTipRatio, sideStart, sideInc, numNetworks);
}

void allSearchedConfig(){
	kisset (7446, 13195, 32375);
	string uniquer("allSearchedConfigLen5-20_c195_");
	bool doDraw(true);
	double areaStart(5.0), areaInc(5.0), areaEnd(20.0),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(9.0),
		maxPathLengthCoefficient(5.0),
		localRad(2.0*nearestRadius);
	int numNetworks(20), numOpt(5), graftPasses(200), settleIts(10000), numBins(25), numAveBins(1000),
		maxFailedInsertions(10000), width(1000), height(width), border(width/40);
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		<< "\n note that area is NOT halved"
		<< "\n areaStart = " << areaStart << "\n areaInc = " << areaInc << "\n areaEnd = " << areaEnd
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n localRad = " << localRad
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses << "\n numBins = " << numBins
		<< "\n numAveBins = "
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<vector<double> > times, measures, sims;
	vector<vector<vector<double> > >lambdaLs, gammas, tipRatios, lengths;
	vector<vector<vector<int> > > tipCounts, levelsToHeart;
	vector<vector<unsigned int> > terminalServiceVolumes;
	vector<vector<int> > ignoreCounts, ignoreCountsGamma;
	vector<vector<vector<Bifurcation2D> >* > bestBs;
	vector<vector<double>* > bestBMeasures;
	vector<vector<double>* > bestBSims;
	for(double area(areaStart); area < areaEnd + areaInc/2.0; area += areaInc){
		double side(sqrt(area)), heartx(side/2.0), hearty(nearestRadius/2.0);//side(sqrt(area))
		times.push_back(vector<double>());
		measures.push_back(vector<double>());
		lambdaLs.push_back(vector<vector<double> >());
		gammas.push_back(vector<vector<double> >());
		tipCounts.push_back(vector<vector<int> >());
		levelsToHeart.push_back(vector<vector<int> >());
		lengths.push_back(vector<vector<double> >());
		tipRatios.push_back(vector<vector<double> >());
		terminalServiceVolumes.push_back(vector<unsigned int>());
		ignoreCounts.push_back(vector<int>());
		ignoreCountsGamma.push_back(vector<int>());
		bestBs.push_back(new vector<vector<Bifurcation2D> >[numNetworks]);
		bestBMeasures.push_back(new vector<double>[numNetworks]);
		bestBSims.push_back(new vector<double>[numNetworks]);
		for(int n(0); n < numNetworks; n++){
			vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, side, side, nearestRadius, heartx, hearty));
			//vector<Bifurcation2D> bCopy(b);
			//introduceBifurcationsFromSources(bCopy);
			//vector<BranchPoint2D> bGuess = consolidateDegenerates(bCopy, consolidateThresh);
			//double guessMeasure = tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			//cout << "\n\t\t\t\t\t\t Guess measure = " << guessMeasure << endl;
			//if(doDraw)
			//	drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_guess.png", bGuess, width, height, border);
			time_t tempStart(time(NULL));
			vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingleAndTracking(bestBs.back()[n], bestBMeasures.back()[n], numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, graftPasses, settleIts);
			for(unsigned int i(0); i < bestBs.back()[n].size(); i++){
				if(bestBMeasures.back()[n][i] < 0.0)
					bestBSims.back()[n].push_back(-1.0);
				else
					bestBSims.back()[n].push_back(subtreeSimilarity(consolidateDegenerates(bestBs.back()[n][0]), consolidateDegenerates(bestBs.back()[n][i])));
			}
			double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			times.back().push_back(difftime(time(NULL), tempStart));
			measures.back().push_back(graftMeasure);
			terminalServiceVolumes.back().push_back(b.size() - 1);
			cout << "\n\t\t Graft took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = " << graftMeasure << endl;
			if(doDraw)
				drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_graft.png", bGraft, width, height, border);
			ignoreCounts.back().push_back(0);
			ignoreCountsGamma.back().push_back(0);
			lambdaLs.back().push_back(validBranchingLengthRatios(bGraft, ignoreCounts.back().back(), lengthThresh));
			//gammas.back().push_back(validParentChildLengthRatios(bGraft, ignoreCountsGamma.back().back(), lengthThresh));
			gammas.back().push_back(vector<double>());
			tipCounts.back().push_back(vector<int>());
			levelsToHeart.back().push_back(vector<int>());
			lengths.back().push_back(vector<double>());
			validParentChildLengthRatiosTipsLevelsToHeartLengths(bGraft, ignoreCountsGamma.back().back(), gammas.back().back(), tipCounts.back().back(), levelsToHeart.back().back(), lengths.back().back(), lengthThresh);
			tipRatios.back().push_back(validTipCountRatios(bGraft, ignoreCounts.back().back(), lengthThresh)); // ignoreCounts should be the same as for lambdaL
		}
	}

	cout << "\n\n area\t <TSV>\t +/-\t fracIgn_rat\t +/-\t fracIgn_gam\t +/-\t measure\t +/-\t <time>\t +/-" << endl;
	for(unsigned int a(0); a < times.size(); a++){
		vector<double> ignoreFracs, ignoreFracsGamma;
		for(int n(0); n < numNetworks; n++){
			ignoreFracs.push_back(double(ignoreCounts[a][n])/double(terminalServiceVolumes[a][n] - 1));
			ignoreFracsGamma.push_back(double(ignoreCountsGamma[a][n])/double(2*terminalServiceVolumes[a][n] - 2));
		}
		cout << areaStart + a*areaInc
			<< "\t" << mean(terminalServiceVolumes[a]) << "\t" << stDev(terminalServiceVolumes[a])
			<< "\t" << mean(ignoreFracs) << "\t" << stDev(ignoreFracs)
			<< "\t" << mean(ignoreFracsGamma) << "\t" << stDev(ignoreFracsGamma)
			<< "\t" << mean(measures[a]) << "\t" << stDev(measures[a])
			<< "\t" << mean(times[a]) << "\t" << stDev(times[a])
			<< endl;
	}


	vector<int*> maxRank;
	vector<map<double, double> > searchedConfigDists;
	for(unsigned int a(0); a < bestBs.size(); a++){
		maxRank.push_back(new int[numNetworks]);
		vector<double> areaRanks;
		for(int n(0); n < numNetworks; n++){
			for(unsigned int i(0); i < bestBMeasures[a][n].size(); i++){
				if(bestBMeasures[a][n][i] >= 0.0)
					maxRank[a][n] = i + 1;
			}
			areaRanks.push_back(maxRank[a][n]);
		}
		cout << "\n\n area = " << areaStart + a*areaInc << endl;
		searchedConfigDists.push_back(probabilityDistribution(areaRanks, 2 + numNetworks/3));
		printMap(searchedConfigDists.back(), "numUnique");
	}
	cout << "\n\n v-rr(0,1]";
	for(unsigned int a(0); a < bestBs.size(); a++)
		cout << "\t {interp_meas/meas[0]}_" << areaStart + a*areaInc << "\t {interp_sim}_" << areaStart + a*areaInc;
	cout << "\n0.0";
	for(unsigned int a(0); a < bestBs.size(); a++)
		cout << "\t1.0\t1.0";
	cout << endl;
	for(int i(0); i < numAveBins; i++){
		double binCenter((i + 0.5)/numAveBins);
		cout << binCenter;
		for(unsigned int a(0); a < bestBs.size(); a++){
			double sum(0.0);
			for(int n(0); n < numNetworks; n++){
				int binBelow((int)floor(binCenter*maxRank[a][n]));
				//cout << " [ n" << n << " binBelow=" << binBelow << " maxRank[n]=" << maxRank[n] << " ] ";
				if(binBelow > maxRank[a][n] - 1)
					cout << "\n Bad bin!" << endl;
				else if(binBelow > maxRank[a][n] - 2)
					sum += bestBMeasures[a][n][binBelow]/bestBMeasures[a][n][0];
				else
					sum += linearInterpolation(double(binBelow)/maxRank[a][n], binCenter, double(binBelow + 1)/maxRank[a][n], bestBMeasures[a][n][binBelow], bestBMeasures[a][n][binBelow + 1])/bestBMeasures[a][n][0];
			}
			cout << "\t" << sum/numNetworks;
			double simSum(0.0);
			for(int n(0); n < numNetworks; n++){
				int binBelow((int)floor(binCenter*maxRank[a][n]));
				//cout << " [ n" << n << " binBelow=" << binBelow << " maxRank[n]=" << maxRank[n] << " ] ";
				if(binBelow > maxRank[a][n] - 1)
					cout << "\n Bad bin! (sim)" << endl;
				else if(binBelow > maxRank[a][n] - 2)
					simSum += bestBSims[a][n][binBelow];
				else
					simSum += linearInterpolation(double(binBelow)/maxRank[a][n], binCenter, double(binBelow + 1)/maxRank[a][n], bestBSims[a][n][binBelow], bestBSims[a][n][binBelow + 1]);
			}
			cout << "\t" << simSum/numNetworks;
		}
		cout << endl;
	}

	vector<vector<map<double, double> > > lambdaL, gamma, tipRatio; // single networks
	vector<map<double, double> > areaLambdaL, areaGamma, areaTipRatio; // within area, over networks
	vector<double> allLambdaLs, allGammas, allTipRatios; // over all areas
	for(unsigned int a(0); a < lambdaLs.size(); a++){
		lambdaL.push_back(vector<map<double, double> >());
		gamma.push_back(vector<map<double, double> >());
		tipRatio.push_back(vector<map<double, double> >());
		vector<double> areaLambdaLs, areaGammas, areaTipRatios; // within area, over networks
		for(int n(0); n < numNetworks; n++){
			areaLambdaLs.insert(areaLambdaLs.end(), lambdaLs[a][n].begin(), lambdaLs[a][n].end());
			areaGammas.insert(areaGammas.end(), gammas[a][n].begin(), gammas[a][n].end());
			areaTipRatios.insert(areaTipRatios.end(), tipRatios[a][n].begin(), tipRatios[a][n].end());
			lambdaL.back().push_back(probabilityDistribution(lambdaLs[a][n], numBins));
			gamma.back().push_back(probabilityDistribution(gammas[a][n], numBins));
			tipRatio.back().push_back(probabilityDistribution(tipRatios[a][n], numBins));
		}
		areaLambdaL.push_back(probabilityDistribution(areaLambdaLs, numBins));
		areaGamma.push_back(probabilityDistribution(areaGammas, numBins));
		areaTipRatio.push_back(probabilityDistribution(areaTipRatios, numBins));
		allLambdaLs.insert(allLambdaLs.end(), areaLambdaLs.begin(), areaLambdaLs.end());
		allGammas.insert(allGammas.end(), areaGammas.begin(), areaGammas.end());
		allTipRatios.insert(allTipRatios.end(), areaTipRatios.begin(), areaTipRatios.end());
	}
	map<double, double> allLambdaL(probabilityDistribution(allLambdaLs, numBins)),
		allGamma(probabilityDistribution(allGammas, numBins)),
		allTipRatio(probabilityDistribution(allTipRatios, numBins));

	cout << "\n\n Scaling data for gamma" << endl;
	cout << "\n areas:";
	for(unsigned int a(0); a < bestBs.size(); a++)
		cout << " " << areaStart + a*areaInc;
	cout << endl;
	cout << "\n gamma\t N_tips\t levelsToHeart\t length" << endl;
	for(unsigned int a(0); a < bestBs.size(); a++){
		for(int n(0); n < numNetworks; n++){
			for(unsigned int i(0); i < gammas[a][n].size(); i++)
				cout << gammas[a][n][i] << "\t" << tipCounts[a][n][i] << "\t" << levelsToHeart[a][n][i] << "\t" << lengths[a][n][i] << endl;
		}
		cout << endl;
	}

	printReport("lambdaL", numBins, lambdaL, areaLambdaL, allLambdaL, areaStart, areaInc, numNetworks);
	printReport("gamma", numBins, gamma, areaGamma, allGamma, areaStart, areaInc, numNetworks);
	printReport("tipRatio", numBins, tipRatio, areaTipRatio, allTipRatio, areaStart, areaInc, numNetworks);

	for(unsigned int a(0); a < bestBs.size(); a++){
		delete[] bestBMeasures[a];
		delete[] bestBs[a];
		delete[] bestBSims[a];
	}
}

void allSearchedConfigPathAndRandom(){
	kisset (7446, 13195, 32375);
	string uniquer("allSearchedConfigPathAndRandomOneCont100_c100_");
	bool doDraw(true);
	double areaStart(10.0), areaInc(areaStart), areaEnd(areaStart),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0),
		localRad(2.0*nearestRadius);
	int numNetworks(1), numOpt(100), numRand(0), numToContinue(1), graftPasses(200), settleIts(10000), numBins(25), numAveBins(1000),
		maxFailedInsertions(10000), width(1000), height(width), border(width/40);
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		<< "\n note that area is INDEED halved"
		<< "\n areaStart = " << areaStart << "\n areaInc = " << areaInc << "\n areaEnd = " << areaEnd
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n localRad = " << localRad << "\n numRand = " << numRand
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses << "\n numBins = " << numBins
		<< "\n numAveBins = " << numAveBins << "\n numToContinue = " << numToContinue
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<vector<double> > times;//, measures, sims;
	//vector<vector<vector<double> > >lambdaLs, gammas, tipRatios, lengths;
	//vector<vector<vector<int> > > tipCounts, levelsToHeart;
	vector<vector<unsigned int> > terminalServiceVolumes;
	//vector<vector<int> > ignoreCounts, ignoreCountsGamma;
	vector<vector<vector<Bifurcation2D> >** > bestBs;
	vector<vector<double>** > bestBMeasures; // bestBMeasures(area, network, opt.or.rand, config)
	vector<vector<double>** > bestBSims, bestBSimsAllRuns;
	vector<int*> bestOpt;
	for(double area(areaStart); area < areaEnd + areaInc/2.0; area += areaInc){
		double side(sqrt(area/2.0)), heartx(side/2.0), hearty(nearestRadius/2.0);//side(sqrt(area))
		times.push_back(vector<double>());
		//measures.push_back(vector<double>());
		terminalServiceVolumes.push_back(vector<unsigned int>());
		bestBs.push_back(new vector<vector<Bifurcation2D> >*[numNetworks]);
		bestBMeasures.push_back(new vector<double>*[numNetworks]);
		bestBSims.push_back(new vector<double>*[numNetworks]);
		bestBSimsAllRuns.push_back(new vector<double>*[numNetworks]);
		bestOpt.push_back(new int[numNetworks]);
		for(int n(0); n < numNetworks; n++){
			vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, side, side, nearestRadius, heartx, hearty));
			terminalServiceVolumes.back().push_back(b.size() - 1);
			vector<Bifurcation2D> bCopy(b);
			introduceBifurcationsFromSources(bCopy);
			vector<BranchPoint2D> bGuess = consolidateDegenerates(bCopy, consolidateThresh);
			double guessMeasure = tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			cout << "\n\t\t\t\t\t\t Guess measure = " << guessMeasure << endl;
			if(doDraw)
				drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_guess.png", bGuess, width, height, border);
			bestBs.back()[n] = new vector<vector<Bifurcation2D> >[numOpt + numRand];
			bestBMeasures.back()[n] = new vector<double>[numOpt + numRand];
			bestBSims.back()[n] = new vector<double>[numOpt + numRand];
			bestBSimsAllRuns.back()[n] = new vector<double>[numOpt + numRand];
			for(int opt(0); opt < numOpt; opt++){
				time_t tempStart(time(NULL));
				vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingleAndTracking(bestBs.back()[n][opt], bestBMeasures.back()[n][opt], 1, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, graftPasses, settleIts);
				for(unsigned int i(0); i < bestBs.back()[n][opt].size(); i++){
					if(bestBMeasures.back()[n][opt][i] < 0.0)
						bestBSims.back()[n][opt].push_back(-1.0);
					else
						bestBSims.back()[n][opt].push_back(subtreeSimilarity(consolidateDegenerates(bestBs.back()[n][opt][0]), consolidateDegenerates(bestBs.back()[n][opt][i])));
				}
				double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
				times.back().push_back(difftime(time(NULL), tempStart));
				//measures.back().push_back(graftMeasure);
				cout << "\n\t\t Graft took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = " << graftMeasure << endl;
				if(doDraw)
					drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_opt" + makeString(opt) + "_graft.png", bGraft, width, height, border);
			}
			for(int opt(numOpt); opt < numOpt + numRand; opt++){
				time_t tempStart(time(NULL));
				vector<Bifurcation2D> bCopy(b);
				introduceBifurcationsFromSourcesRandomly(bCopy);
				purelyDescendingConsolidatedSearch(bestBs.back()[n][opt], numToContinue, bCopy, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh);
				for(unsigned int i(0); i < bestBs.back()[n][opt].size(); i++){
					vector<BranchPoint2D> bestBci(consolidateDegenerates(bestBs.back()[n][opt][i], consolidateThresh));
					bestBMeasures.back()[n][opt].push_back(tripletMeasure(bestBci, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient));
				}
				bool changed(true);
				while(changed){
					changed = false;
					for(unsigned int i(0); i < bestBs.back()[n][opt].size() - 1; i++){
						if(bestBMeasures.back()[n][opt][i] > bestBMeasures.back()[n][opt][i + 1]){
							swap(bestBMeasures.back()[n][opt][i], bestBMeasures.back()[n][opt][i + 1]);
							swap(bestBs.back()[n][opt][i], bestBs.back()[n][opt][i + 1]);
							changed = true;
						}
					}
				}
				for(unsigned int i(0); i < bestBs.back()[n][opt].size(); i++)
					bestBSims.back()[n][opt].push_back(subtreeSimilarity(consolidateDegenerates(bestBs.back()[n][opt][0]), consolidateDegenerates(bestBs.back()[n][opt][i])));
				times.back().push_back(difftime(time(NULL), tempStart));
				cout << "\n\t\t Purely Descending took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs\t\t measure = " << bestBMeasures.back()[n][opt][0] << endl;
				if(doDraw)
					drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_rand" + makeString(opt) + "_descending.png", consolidateDegenerates(bestBs.back()[n][opt][0]), width, height, border);
			}
			bestOpt.back()[n] = 0;
			for(int opt(1); opt < numOpt + numRand; opt++){
				if(bestBMeasures.back()[n][bestOpt.back()[n]][0] > bestBMeasures.back()[n][opt][0])
					bestOpt.back()[n] = opt;
			}
			for(int opt(0); opt < numOpt + numRand; opt++){
				for(unsigned int i(0); i < bestBs.back()[n][opt].size(); i++)
					bestBSimsAllRuns.back()[n][opt].push_back(subtreeSimilarity(consolidateDegenerates(bestBs.back()[n][bestOpt.back()[n]][0]), consolidateDegenerates(bestBs.back()[n][opt][i])));
			}
		}
	}

	// Analysis...
	for(unsigned int a(0); a < bestBs.size(); a++){
		double area(areaStart + a*areaInc);
		for(int n(0); n < numNetworks; n++){
			for(int opt(0); opt < numOpt + numRand; opt++){
				cout << "\n\n area = " << area
					<< "\t\t TSVs = " << terminalServiceVolumes[a][n]
					<< "\n\t n = " << n
					<< "\n\t opt = " << opt;
				if(opt >= numOpt)
					cout << " (rand)" << endl;
				else cout << endl;
				cout << "\n rr\t meas\t meas/meas0(within)\t meas/meas0(overall)\t sim(within)\t sim(overall)" << endl;
				for(unsigned int i(0); i < bestBs[a][n][opt].size(); i++){
					cout << double(i)/double(bestBs[a][n][opt].size() - 1)
						<< "\t" << bestBMeasures[a][n][opt][i]
						<< "\t" << bestBMeasures[a][n][opt][i]/bestBMeasures[a][n][opt][0]
						<< "\t" << bestBMeasures[a][n][opt][i]/bestBMeasures[a][n][bestOpt[a][n]][0]
						<< "\t" << bestBSims[a][n][opt][i]
						<< "\t" << bestBSimsAllRuns[a][n][opt][i]
						<< endl;
				}
			}
		}
	}

	// ADAPTED FROM EXHAUSTIVE ANALYSIS...  // bestBMeasures(area, network, opt.or.rand, config)



	for(unsigned int a(0); a < bestBs.size(); a++){
		for(int n(0); n < numNetworks; n++){
			cout << "\n\n area = " << areaStart + a*areaInc << ", network " << n
				<< " RANKED BY MEASURE\n v-rr(0,1]\t interp_meas/meas[0]_overall \t simToBest_overall" << endl;
			double bestNetworkBMeasure(bestBMeasures[a][n][0][0]);
			for(int opt(1); opt < numOpt + numRand; opt++){
				if(bestNetworkBMeasure > bestBMeasures[a][n][opt][0])
					bestNetworkBMeasure = bestBMeasures[a][n][opt][0];
			}
			//cout << "0.0\t1.0\t1.0" << endl;
			for(int i(0); i < numAveBins; i++){
				double binCenter((i + 0.5)/numAveBins);
				double sum(0.0), simSum(0.0);
				for(int opt(0); opt < numOpt + numRand; opt++){
					int maxRank = bestBs[a][n][opt].size();
					int binBelow((int)floor(binCenter*maxRank));
					//cout << " [ n" << n << " binBelow=" << binBelow << " maxRank[n]=" << maxRank[n] << " ] ";
					if(binBelow > maxRank - 1)
						cout << "\n Bad bin!" << endl;
					else if(binBelow > maxRank - 2){
						sum += bestBMeasures[a][n][opt][binBelow]/bestNetworkBMeasure;
						simSum += bestBSimsAllRuns[a][n][opt][binBelow];
					}else{
						sum += linearInterpolation(double(binBelow)/maxRank, binCenter, double(binBelow + 1)/maxRank, bestBMeasures[a][n][opt][binBelow], bestBMeasures[a][n][opt][binBelow + 1])/bestNetworkBMeasure;
						double simBelow(bestBSimsAllRuns[a][n][opt][binBelow]), simNext(bestBSimsAllRuns[a][n][opt][binBelow + 1]);
						simSum += linearInterpolation(double(binBelow)/maxRank, binCenter, double(binBelow + 1)/maxRank, simBelow, simNext);
					}
				}
				cout << binCenter << "\t" << sum/(numOpt + numRand) << "\t" << simSum/(numOpt + numRand) << endl;
			}

			vector<double> uniqueCount;
			unsigned int minConfig(bestBs[a][n][0].size()), maxConfig(minConfig);
			for(int opt(0); opt < numOpt + numRand; opt++){
				uniqueCount.push_back(double(bestBs[a][n][opt].size()));
				if(minConfig > bestBs[a][n][opt].size())
					minConfig = bestBs[a][n][opt].size();
				if(maxConfig < bestBs[a][n][opt].size())
					maxConfig = bestBs[a][n][opt].size();
			}
			unsigned int uniqueBins(numOpt + numRand);
			if(uniqueBins > maxConfig - minConfig)
				uniqueBins = maxConfig - minConfig;
			map<double, double> maxRankDist(probabilityDistribution(uniqueCount, uniqueBins));
			printMap(maxRankDist, "max_rank");

		}
	}

	/*
	// re-order for decreasing similarity
	for(int n(0); n < numNetworks; n++){
		bool changed(true);
		while(changed){
			changed = false;
			for(int i(0); i < numBest - 1; i++){
				if(bestBMeasures[n][i] < 0.0)
					i = numBest;
				else{
					if(sims[n][i] < sims[n][i + 1]){
						changed = true;
						swap(bestBs[n][i], bestBs[n][i + 1]);
						swap(bestBMeasures[n][i], bestBMeasures[n][i + 1]);
						swap(sims[n][i], sims[n][i + 1]);
					}
				}
			}
		}
	}

	cout << "\n\n RANKED BY SIMILARITY\n v-rr(0,1]\t interp_meas/meas[0] \t simToBest" << endl;
	cout << "0.0\t1.0\t1.0" << endl;
	for(int i(0); i < numAveBins; i++){
		double binCenter((i + 0.5)/numAveBins);
		double sum(0.0), simSum(0.0);
		for(int n(0); n < numNetworks; n++){
			int binBelow((int)floor(binCenter*maxRank[n]));
			//cout << " [ n" << n << " binBelow=" << binBelow << " maxRank[n]=" << maxRank[n] << " ] ";
			if(binBelow > maxRank[n] - 1)
				cout << "\n Bad bin!" << endl;
			else if(binBelow > maxRank[n] - 2){
				sum += bestBMeasures[n][binBelow]/bestBMeasures[n][0];
				simSum += sims[n][binBelow];
			}else{
				sum += linearInterpolation(double(binBelow)/maxRank[n], binCenter, double(binBelow + 1)/maxRank[n], bestBMeasures[n][binBelow], bestBMeasures[n][binBelow + 1])/bestBMeasures[n][0];
				double simBelow(sims[n][binBelow]), simNext(sims[n][binBelow + 1]);
				simSum += linearInterpolation(double(binBelow)/maxRank[n], binCenter, double(binBelow + 1)/maxRank[n], simBelow, simNext);
			}
		}
		cout << binCenter << "\t" << sum/numNetworks << "\t" << simSum/numNetworks << endl;
	}
	*/

	// Nothing has been cleaned up yet...
}

void scalingDataForGammaNoSearch(){
	kisset (7446, 13195, 32375);
	string uniquer("outputs/scalingDataForGammaNoSearch200-1000_");
	bool doDraw(true);
	double areaStart(200.0), areaInc(areaStart), areaEnd(1000.0),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(9.0), // optimization not performed
		maxPathLengthCoefficient(5.0),
		localRad(2.0*nearestRadius);
	int numNetworks(10), numOpt(5), graftPasses(200), settleIts(10000), numBins(25), numAveBins(1000),
		maxFailedInsertions(10000), width(1000), height(width), border(width/40);
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		<< "\n note that area is NOT halved"
		<< "\n areaStart = " << areaStart << "\n areaInc = " << areaInc << "\n areaEnd = " << areaEnd
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n localRad = " << localRad
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses << "\n numBins = " << numBins
		<< "\n numAveBins = "
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<vector<double> > times, measures, sims;
	vector<vector<vector<double> > >lambdaLs, gammas, tipRatios, lengths;
	vector<vector<vector<int> > > tipCounts, levelsToHeart;
	vector<vector<unsigned int> > terminalServiceVolumes;
	vector<vector<int> > ignoreCounts, ignoreCountsGamma;
	vector<vector<vector<Bifurcation2D> >* > bestBs;
	vector<vector<double>* > bestBMeasures;
	vector<vector<double>* > bestBSims;
	for(double area(areaStart); area < areaEnd + areaInc/2.0; area += areaInc){
		double side(sqrt(area)), heartx(side/2.0), hearty(nearestRadius/2.0);//side(sqrt(area))
		times.push_back(vector<double>());
		measures.push_back(vector<double>());
		lambdaLs.push_back(vector<vector<double> >());
		gammas.push_back(vector<vector<double> >());
		tipCounts.push_back(vector<vector<int> >());
		levelsToHeart.push_back(vector<vector<int> >());
		lengths.push_back(vector<vector<double> >());
		tipRatios.push_back(vector<vector<double> >());
		terminalServiceVolumes.push_back(vector<unsigned int>());
		ignoreCounts.push_back(vector<int>());
		ignoreCountsGamma.push_back(vector<int>());
		bestBs.push_back(new vector<vector<Bifurcation2D> >[numNetworks]);
		bestBMeasures.push_back(new vector<double>[numNetworks]);
		bestBSims.push_back(new vector<double>[numNetworks]);
		for(int n(0); n < numNetworks; n++){
			vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, side, side, nearestRadius, heartx, hearty));
			//vector<Bifurcation2D> bCopy(b);
			//introduceBifurcationsFromSources(bCopy);
			//vector<BranchPoint2D> bGuess = consolidateDegenerates(bCopy, consolidateThresh);
			//double guessMeasure = tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			//cout << "\n\t\t\t\t\t\t Guess measure = " << guessMeasure << endl;
			//if(doDraw)
			//	drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_guess.png", bGuess, width, height, border);
			
			//time_t tempStart(time(NULL));
			//vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingleAndTracking(bestBs.back()[n], bestBMeasures.back()[n], numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, graftPasses, settleIts);
			//for(unsigned int i(0); i < bestBs.back()[n].size(); i++){
			//	if(bestBMeasures.back()[n][i] < 0.0)
			//		bestBSims.back()[n].push_back(-1.0);
			//	else
			//		bestBSims.back()[n].push_back(subtreeSimilarity(consolidateDegenerates(bestBs.back()[n][0]), consolidateDegenerates(bestBs.back()[n][i])));
			//}
			//double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			//times.back().push_back(difftime(time(NULL), tempStart));
			//measures.back().push_back(graftMeasure);

			vector<Bifurcation2D> bCopy(b);
			introduceBifurcationsFromSources(bCopy);
			vector<BranchPoint2D> bNoSearch(consolidateDegenerates(bCopy));
			terminalServiceVolumes.back().push_back(b.size() - 1);
			//cout << "\n\t\t Graft took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = " << graftMeasure << endl;
			if(doDraw)
				drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_noSearch.png", bNoSearch, width, height, border);
			ignoreCounts.back().push_back(0);
			ignoreCountsGamma.back().push_back(0);
			lambdaLs.back().push_back(validBranchingLengthRatios(bNoSearch, ignoreCounts.back().back(), lengthThresh));
			//gammas.back().push_back(validParentChildLengthRatios(bGraft, ignoreCountsGamma.back().back(), lengthThresh));
			gammas.back().push_back(vector<double>());
			tipCounts.back().push_back(vector<int>());
			levelsToHeart.back().push_back(vector<int>());
			lengths.back().push_back(vector<double>());
			validParentChildLengthRatiosTipsLevelsToHeartLengths(bNoSearch, ignoreCountsGamma.back().back(), gammas.back().back(), tipCounts.back().back(), levelsToHeart.back().back(), lengths.back().back(), lengthThresh);
			tipRatios.back().push_back(validTipCountRatios(bNoSearch, ignoreCounts.back().back(), lengthThresh)); // ignoreCounts should be the same as for lambdaL
		}
	}

	/*cout << "\n\n area\t <TSV>\t +/-\t fracIgn_rat\t +/-\t fracIgn_gam\t +/-\t measure\t +/-\t <time>\t +/-" << endl;
	for(unsigned int a(0); a < times.size(); a++){
		vector<double> ignoreFracs, ignoreFracsGamma;
		for(int n(0); n < numNetworks; n++){
			ignoreFracs.push_back(double(ignoreCounts[a][n])/double(terminalServiceVolumes[a][n] - 1));
			ignoreFracsGamma.push_back(double(ignoreCountsGamma[a][n])/double(2*terminalServiceVolumes[a][n] - 2));
		}
		cout << areaStart + a*areaInc
			<< "\t" << mean(terminalServiceVolumes[a]) << "\t" << stDev(terminalServiceVolumes[a])
			<< "\t" << mean(ignoreFracs) << "\t" << stDev(ignoreFracs)
			<< "\t" << mean(ignoreFracsGamma) << "\t" << stDev(ignoreFracsGamma)
			<< "\t" << mean(measures[a]) << "\t" << stDev(measures[a])
			<< "\t" << mean(times[a]) << "\t" << stDev(times[a])
			<< endl;
	}


	vector<int*> maxRank;
	vector<map<double, double> > searchedConfigDists;
	for(unsigned int a(0); a < bestBs.size(); a++){
		maxRank.push_back(new int[numNetworks]);
		vector<double> areaRanks;
		for(int n(0); n < numNetworks; n++){
			for(unsigned int i(0); i < bestBMeasures[a][n].size(); i++){
				if(bestBMeasures[a][n][i] >= 0.0)
					maxRank[a][n] = i + 1;
			}
			areaRanks.push_back(maxRank[a][n]);
		}
		cout << "\n\n area = " << areaStart + a*areaInc << endl;
		searchedConfigDists.push_back(probabilityDistribution(areaRanks, 2 + numNetworks/3));
		printMap(searchedConfigDists.back(), "numUnique");
	}
	cout << "\n\n v-rr(0,1]";
	for(unsigned int a(0); a < bestBs.size(); a++)
		cout << "\t {interp_meas/meas[0]}_" << areaStart + a*areaInc << "\t {interp_sim}_" << areaStart + a*areaInc;
	cout << "\n0.0";
	for(unsigned int a(0); a < bestBs.size(); a++)
		cout << "\t1.0\t1.0";
	cout << endl;
	for(int i(0); i < numAveBins; i++){
		double binCenter((i + 0.5)/numAveBins);
		cout << binCenter;
		for(unsigned int a(0); a < bestBs.size(); a++){
			double sum(0.0);
			for(int n(0); n < numNetworks; n++){
				int binBelow((int)floor(binCenter*maxRank[a][n]));
				//cout << " [ n" << n << " binBelow=" << binBelow << " maxRank[n]=" << maxRank[n] << " ] ";
				if(binBelow > maxRank[a][n] - 1)
					cout << "\n Bad bin!" << endl;
				else if(binBelow > maxRank[a][n] - 2)
					sum += bestBMeasures[a][n][binBelow]/bestBMeasures[a][n][0];
				else
					sum += linearInterpolation(double(binBelow)/maxRank[a][n], binCenter, double(binBelow + 1)/maxRank[a][n], bestBMeasures[a][n][binBelow], bestBMeasures[a][n][binBelow + 1])/bestBMeasures[a][n][0];
			}
			cout << "\t" << sum/numNetworks;
			double simSum(0.0);
			for(int n(0); n < numNetworks; n++){
				int binBelow((int)floor(binCenter*maxRank[a][n]));
				//cout << " [ n" << n << " binBelow=" << binBelow << " maxRank[n]=" << maxRank[n] << " ] ";
				if(binBelow > maxRank[a][n] - 1)
					cout << "\n Bad bin! (sim)" << endl;
				else if(binBelow > maxRank[a][n] - 2)
					simSum += bestBSims[a][n][binBelow];
				else
					simSum += linearInterpolation(double(binBelow)/maxRank[a][n], binCenter, double(binBelow + 1)/maxRank[a][n], bestBSims[a][n][binBelow], bestBSims[a][n][binBelow + 1]);
			}
			cout << "\t" << simSum/numNetworks;
		}
		cout << endl;
	}*/

	vector<vector<map<double, double> > > lambdaL, gamma, tipRatio; // single networks
	vector<map<double, double> > areaLambdaL, areaGamma, areaTipRatio; // within area, over networks
	vector<double> allLambdaLs, allGammas, allTipRatios; // over all areas
	for(unsigned int a(0); a < lambdaLs.size(); a++){
		lambdaL.push_back(vector<map<double, double> >());
		gamma.push_back(vector<map<double, double> >());
		tipRatio.push_back(vector<map<double, double> >());
		vector<double> areaLambdaLs, areaGammas, areaTipRatios; // within area, over networks
		for(int n(0); n < numNetworks; n++){
			areaLambdaLs.insert(areaLambdaLs.end(), lambdaLs[a][n].begin(), lambdaLs[a][n].end());
			areaGammas.insert(areaGammas.end(), gammas[a][n].begin(), gammas[a][n].end());
			areaTipRatios.insert(areaTipRatios.end(), tipRatios[a][n].begin(), tipRatios[a][n].end());
			lambdaL.back().push_back(probabilityDistribution(lambdaLs[a][n], numBins));
			gamma.back().push_back(probabilityDistribution(gammas[a][n], numBins));
			tipRatio.back().push_back(probabilityDistribution(tipRatios[a][n], numBins));
		}
		areaLambdaL.push_back(probabilityDistribution(areaLambdaLs, numBins));
		areaGamma.push_back(probabilityDistribution(areaGammas, numBins));
		areaTipRatio.push_back(probabilityDistribution(areaTipRatios, numBins));
		allLambdaLs.insert(allLambdaLs.end(), areaLambdaLs.begin(), areaLambdaLs.end());
		allGammas.insert(allGammas.end(), areaGammas.begin(), areaGammas.end());
		allTipRatios.insert(allTipRatios.end(), areaTipRatios.begin(), areaTipRatios.end());
	}
	map<double, double> allLambdaL(probabilityDistribution(allLambdaLs, numBins)),
		allGamma(probabilityDistribution(allGammas, numBins)),
		allTipRatio(probabilityDistribution(allTipRatios, numBins));

	cout << "\n\n Scaling data for gamma" << endl;
	cout << "\n areas:";
	for(unsigned int a(0); a < bestBs.size(); a++)
		cout << " " << areaStart + a*areaInc;
	cout << endl;
	cout << "\n gamma\t N_tips\t levelsToHeart\t length" << endl;
	for(unsigned int a(0); a < bestBs.size(); a++){
		for(int n(0); n < numNetworks; n++){
			for(unsigned int i(0); i < gammas[a][n].size(); i++)
				cout << gammas[a][n][i] << "\t" << tipCounts[a][n][i] << "\t" << levelsToHeart[a][n][i] << "\t" << lengths[a][n][i] << endl;
		}
		cout << endl;
	}

	/*printReport("lambdaL", numBins, lambdaL, areaLambdaL, allLambdaL, areaStart, areaInc, numNetworks);
	printReport("gamma", numBins, gamma, areaGamma, allGamma, areaStart, areaInc, numNetworks);
	printReport("tipRatio", numBins, tipRatio, areaTipRatio, allTipRatio, areaStart, areaInc, numNetworks);

	for(unsigned int a(0); a < bestBs.size(); a++){
		delete[] bestBMeasures[a];
		delete[] bestBs[a];
		delete[] bestBSims[a];
	}
	*/

	// nothing has been cleaned up yet . . .
}

void optConsistency(){
	kisset (7446, 13195, 32375);
	string uniquer("outputs/optConsistency_c195_");
	bool doDraw(true);
	double area(10.0),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(9.0),
		maxPathLengthCoefficient(5.0);
	int numNetworks(1), numOpt(10), numRepeats(5), graftPasses(100), settleIts(10000),
		maxFailedInsertions(10000), width(500), height(width), border(width/20);
	cout << "\n uniquer = " << uniquer
		<< "\n area = " << area
		<< "\n settleThresh = " << settleThresh << "\n consolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;

	vector<vector<double> > times, measures;
	vector<int>terminalServiceVolumes;
	double side(sqrt(area/2.0)), heartx(side/2.0), hearty(nearestRadius/2.0);//side(sqrt(area))
	for(int n(0); n < numNetworks; n++){
		times.push_back(vector<double>());
		measures.push_back(vector<double>());
		vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, side, side, nearestRadius, heartx, hearty));
		terminalServiceVolumes.push_back(b.size());
		for(int r(0); r < numRepeats; r++){
			time_t tempStart(time(NULL));
			vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGrafting(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, graftPasses, settleIts);
			double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			times.back().push_back(difftime(time(NULL), tempStart));
			measures.back().push_back(graftMeasure);
			cout << "\n\t\t Graft took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = " << graftMeasure << endl;
			if(doDraw)
				drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_r" + makeString(r) + "_graft.png", bGraft, width, height, border);
		}
	}

	cout << "\n\n n\t r\t meaure" << endl;
	for(int n(0); n < numNetworks; n++){
		for(int r(0); r < numRepeats; r++)
			cout << n << "\t" << r << "\t" << measures[n][r] << endl;
	}

	cout << "\n\n <TSV> =\t" << mean(terminalServiceVolumes) << "\t +/t\t" << stDev(terminalServiceVolumes)
		<< "\n\nn\t <meas>\t +/-\t time\t +/-" << endl;
	for(int n(0); n < numNetworks; n++){
		cout << n << "\t" << mean(measures[n]) << "\t" << stDev(measures[n])
			<< "\t" << mean(times[n]) << "\t" << stDev(times[n])
			<< endl;
	}
}

void testOverrelaxation(){
	//kisset (7446, 13195, 32375);
	kisset (9267, 2654, 4259);
	string uniquer("testOverrelaxation_c195_");
	bool doDraw(true);
	double area(50.0),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(9.0),
		maxPathLengthCoefficient(5.0),
		relaxStart(0.0), relaxInc(0.03), relaxEnd(0.16);
	int numNetworks(10), numOpt(5), graftPasses(100), settleIts(1000),
		maxFailedInsertions(10000), width(500), height(width), border(width/20);
	cout << "\n uniquer = " << uniquer
		<< "\n area = " << area
		<< "\n settleThresh = " << settleThresh << "\n consolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;

	vector<vector<double> > times, measures;
	vector<int> terminalServiceVolumes;
	double side(sqrt(area/2.0)), heartx(side/2.0), hearty(nearestRadius/2.0);//side(sqrt(area))
	for(int n(0); n < numNetworks; n++){
		times.push_back(vector<double>());
		measures.push_back(vector<double>());
		unsigned long ii, jj, kk;
		kisseed(ii, jj, kk);
		kisset(ii, jj, kk);
		vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, side, side, nearestRadius, heartx, hearty));
		terminalServiceVolumes.push_back(b.size());
		for(double relax(relaxStart); relax < relaxEnd + relaxInc/2.0; relax += relaxInc){
			time_t tempStart(time(NULL));
			overrelaxationFactor = relax;
			//kisset(1, 222, 333333);
			//kisset(543816315, 1769402960, 1382217032); // takes a while
			//kisset(3096135201, 1012496761, 1822534974); // doesn't settle
			//kisset(3427749460, 1192877312, 138985079); // exceeds passes
			vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGrafting(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, graftPasses, settleIts);
			double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			times.back().push_back(difftime(time(NULL), tempStart));
			measures.back().push_back(graftMeasure);
			cout << "\n\t\t n = " << n << "\t relax = " << relax << " \tGraft took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = " << graftMeasure << endl;
			if(doDraw)
				drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_r" + makeString(relax) + "_graft.png", bGraft, width, height, border);
		}
	}

	cout << "\n\n n\t <time>\t +/-" << endl;
	int relaxIndex(0);
	for(double relax(relaxStart); relax < relaxEnd + relaxInc/2.0; relax += relaxInc){
		vector<double> relaxTimes;
		for(int n(0); n < numNetworks; n++)
			relaxTimes.push_back(times[n][relaxIndex]);
		cout << relax << "\t" << mean(relaxTimes) << "\t" << stDev(relaxTimes) << endl;
		relaxIndex++;
	}

	cout << "\n\n <TSV> =\t" << mean(terminalServiceVolumes) << "\t +/t\t" << stDev(terminalServiceVolumes)
		<< "\n\nn\t <meas>\t +/-\t time\t +/-" << endl;
	for(int n(0); n < numNetworks; n++){
		cout << n << "\t" << mean(measures[n]) << "\t" << stDev(measures[n])
			<< "\t" << mean(times[n]) << "\t" << stDev(times[n])
			<< endl;
	}
}

void drawLegend(){
	int width(100), height(300);
	double lineRadius(2.0);
	string fn("outputs/legend.png");
	int red(111), green(55), blue(55), // for lines
		redTip(0), greenTip(111), blueTip(0),
		redNode(0), greenNode(0), blueNode(111),
		redHeart(111), greenHeart(0), blueHeart(0),
		standardRadius(width/10);
	vector<unsigned char> image;
	image.resize(width*height*4);
	for(int i(0); i < width*height*4; i++)
		image[i] = 255;
	drawLine(image, width, height, red, green, blue, width/4, 7*height/8, 3*width/4, 7*height/8, lineRadius);
	drawCircle(image, width, height, redNode, greenNode, blueNode, width/2, 5*height/8, standardRadius/2, lineRadius);
	drawSquare(image, width, height, redTip, greenTip, blueTip, width/2, 3*height/8, standardRadius, lineRadius);
	drawTriangle(image, width, height, redHeart, greenHeart, blueHeart, width/2, 1*height/8, 2*standardRadius, lineRadius);
	unsigned error = lodepng::encode(fn.c_str(), image, width, height);
	if(error)
		cout << "\nError snapshotPNG(): error " << error << " from lodepng" << endl;
}

// it would be better to have the configuration have no crossings because of optimality
//vector<BranchPoint2D> rerouteCrossings(const vector<BranchPoint2D> &b, int maxPasses = 1000){
//	int passCount(0), numRerouted(1);
//	while(numRerouted > 0 && passCount < maxPasses){
//
//	}
//}

//void randomlySelectedConfigurations(){
//
//}

void test3D(){
	//kisset (7446, 13195, 32375);
	string uniquer("test3D195_");
	bool doDraw(true);
	double volStart(90.0), volInc(10.0), volEnd(1.0*volStart),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0), //pow(3.0/4.0/acos(-1.0), 1.0/5.0)
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(9.0),
		maxPathLengthCoefficient(5.0);
	int numNetworks(100), numOpt(10), graftPasses(100), settleIts(10000), numBins(25),
		maxFailedInsertions(10000), width(1000), height(width), border(width/20);
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		<< "\n volStart = " << volStart << "\n volInc = " << volInc << "\n volEnd = " << volEnd
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses << "\n numBins = " << numBins
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<vector<double> > times, measures;
	vector<vector<vector<double> > >lambdaLs, gammas, tipRatios;
	vector<vector<unsigned int> > terminalServiceVolumes;
	vector<vector<int> > ignoreCounts, ignoreCountsGamma;
	for(double vol(volStart); vol < volEnd + volInc/2.0; vol += volInc){
		double side(pow(vol/4.0, 1.0/3.0)), heartx(side/2.0), hearty(nearestRadius/2.0), heartz(side/2.0);//side(sqrt(area))
		times.push_back(vector<double>());
		measures.push_back(vector<double>());
		lambdaLs.push_back(vector<vector<double> >());
		gammas.push_back(vector<vector<double> >());
		tipRatios.push_back(vector<vector<double> >());
		terminalServiceVolumes.push_back(vector<unsigned int>());
		ignoreCounts.push_back(vector<int>());
		ignoreCountsGamma.push_back(vector<int>());
		for(int n(0); n < numNetworks; n++){
			//kisset(24698, 17231, 19458); // first network does not settle
			//kisset(3139653558, 3603707825, 437701738) // does not settle, different threshes
			//kisset(4285140577, 1936677691, 1822336653) // left off
			kissprint();
			vector<Bifurcation3D> b(initializeEvenRandomSpacing(maxFailedInsertions, side, side, side, nearestRadius, heartx, hearty, heartz));
			//vector<Bifurcation2D> bCopy(b);
			//introduceBifurcationsFromSources(bCopy);
			//vector<BranchPoint2D> bGuess = consolidateDegenerates(bCopy, consolidateThresh);
			//double guessMeasure = tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			//cout << "\n\t\t\t\t\t\t Guess measure = " << guessMeasure << endl;
			//if(doDraw)
			//	drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_guess.png", bGuess, width, height, border);
			time_t tempStart(time(NULL));
			vector<BranchPoint3D> bGraft = optimizeByGlobalConsolidatedGrafting(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, graftPasses, settleIts);
			double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			times.back().push_back(difftime(time(NULL), tempStart));
			measures.back().push_back(graftMeasure);
			terminalServiceVolumes.back().push_back(b.size() - 1);
			cout << "\n\t\t Graft took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = " << graftMeasure << endl;
			if(doDraw)
				drawTree(uniquer + makeString((int)ceil(vol)) + "_n" + makeString(n) + "_graft.png", bGraft, width, height, border);
			ignoreCounts.back().push_back(0);
			ignoreCountsGamma.back().push_back(0);
			lambdaLs.back().push_back(validBranchingLengthRatios(bGraft, ignoreCounts.back().back(), lengthThresh));
			gammas.back().push_back(validParentChildLengthRatios(bGraft, ignoreCountsGamma.back().back(), lengthThresh));
			tipRatios.back().push_back(validTipCountRatios(bGraft, ignoreCounts.back().back(), lengthThresh)); // ignoreCounts should be the same as for lambdaL
		}
	}
	cout << "\n\n area\t <TSV>\t +/-\t fracIgn_rat\t +/-\t fracIgn_gam\t +/-\t measure\t +/-\t <time>\t +/-" << endl;
	for(unsigned int a(0); a < times.size(); a++){
		vector<double> ignoreFracs, ignoreFracsGamma;
		for(int n(0); n < numNetworks; n++){
			ignoreFracs.push_back(double(ignoreCounts[a][n])/double(terminalServiceVolumes[a][n] - 1));
			ignoreFracsGamma.push_back(double(ignoreCountsGamma[a][n])/double(2*terminalServiceVolumes[a][n] - 2));
		}
		cout << volStart + a*volInc
			<< "\t" << mean(terminalServiceVolumes[a]) << "\t" << stDev(terminalServiceVolumes[a])
			<< "\t" << mean(ignoreFracs) << "\t" << stDev(ignoreFracs)
			<< "\t" << mean(ignoreFracsGamma) << "\t" << stDev(ignoreFracsGamma)
			<< "\t" << mean(measures[a]) << "\t" << stDev(measures[a])
			<< "\t" << mean(times[a]) << "\t" << stDev(times[a])
			<< endl;
	}

	vector<vector<map<double, double> > > lambdaL, gamma, tipRatio; // single networks
	vector<map<double, double> > areaLambdaL, areaGamma, areaTipRatio; // within area, over networks
	vector<double> allLambdaLs, allGammas, allTipRatios; // over all areas
	for(unsigned int a(0); a < lambdaLs.size(); a++){
		lambdaL.push_back(vector<map<double, double> >());
		gamma.push_back(vector<map<double, double> >());
		tipRatio.push_back(vector<map<double, double> >());
		vector<double> areaLambdaLs, areaGammas, areaTipRatios; // within area, over networks
		for(int n(0); n < numNetworks; n++){
			areaLambdaLs.insert(areaLambdaLs.end(), lambdaLs[a][n].begin(), lambdaLs[a][n].end());
			areaGammas.insert(areaGammas.end(), gammas[a][n].begin(), gammas[a][n].end());
			areaTipRatios.insert(areaTipRatios.end(), tipRatios[a][n].begin(), tipRatios[a][n].end());
			lambdaL.back().push_back(probabilityDistribution(lambdaLs[a][n], numBins));
			gamma.back().push_back(probabilityDistribution(gammas[a][n], numBins));
			tipRatio.back().push_back(probabilityDistribution(tipRatios[a][n], numBins));
		}
		areaLambdaL.push_back(probabilityDistribution(areaLambdaLs, numBins));
		areaGamma.push_back(probabilityDistribution(areaGammas, numBins));
		areaTipRatio.push_back(probabilityDistribution(areaTipRatios, numBins));
		allLambdaLs.insert(allLambdaLs.end(), areaLambdaLs.begin(), areaLambdaLs.end());
		allGammas.insert(allGammas.end(), areaGammas.begin(), areaGammas.end());
		allTipRatios.insert(allTipRatios.end(), areaTipRatios.begin(), areaTipRatios.end());
	}
	map<double, double> allLambdaL(probabilityDistribution(allLambdaLs, numBins)),
		allGamma(probabilityDistribution(allGammas, numBins)),
		allTipRatio(probabilityDistribution(allTipRatios, numBins));

	printReport("lambdaL", numBins, lambdaL, areaLambdaL, allLambdaL, volStart, volInc, numNetworks);
	printReport("gamma", numBins, gamma, areaGamma, allGamma, volStart, volInc, numNetworks);
	printReport("tipRatio", numBins, tipRatio, areaTipRatio, allTipRatio, volStart, volInc, numNetworks);
}

void localSingleGraftTest(){
	//kisset (9267, 2654, 4259);
	string uniquer("outputs/localSingleGraftTest_c195_");
	bool doDraw(true);
	double area(10.0),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(9.0),
		maxPathLengthCoefficient(5.0),
		localRadStart(0.5*nearestRadius), localRadInc(0.5*nearestRadius), localRadEnd(2.0*nearestRadius);
	int numNetworks(5), numOpt(3), graftPasses(100), settleIts(1000),
		maxFailedInsertions(10000), width(500), height(width), border(width/20);
	cout << "\n uniquer = " << uniquer
		<< "\n area = " << area
		<< "\n settleThresh = " << settleThresh << "\n consolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;

	vector<double> times, measures;
	vector<vector<double> > timesSingle, measuresSingle;
	vector<int> terminalServiceVolumes;
	double side(sqrt(area/2.0)), heartx(side/2.0), hearty(nearestRadius/2.0);//side(sqrt(area))
	for(int n(0); n < numNetworks; n++){
		timesSingle.push_back(vector<double>());
		measuresSingle.push_back(vector<double>());
		//unsigned long ii, jj, kk;
		//kisseed(ii, jj, kk);
		//kisset(rand(), rand(), rand());
		vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, side, side, nearestRadius, heartx, hearty));
		terminalServiceVolumes.push_back(b.size());
		time_t tempStart(time(NULL));
		kisset(1, 222, 333333);
		vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGrafting(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, graftPasses, settleIts);
		double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
		times.push_back(difftime(time(NULL), tempStart));
		measures.push_back(graftMeasure);
		if(doDraw){
			drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_graft.png", bGraft, width, height, border);
		}
		for(double localRad(localRadStart); localRad < localRadEnd + localRadInc/2.0; localRad += localRadInc){
			tempStart = time(NULL);
			kisset(1, 222, 333333);
			vector<BranchPoint2D> bGraftSingle = optimizeByGlobalConsolidatedGraftingWithSingle(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, graftPasses, settleIts);
			double graftSingleMeasure = tripletMeasure(bGraftSingle, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			timesSingle.back().push_back(difftime(time(NULL), tempStart));
			measuresSingle.back().push_back(graftSingleMeasure);
			cout << "\n\t\t n = " << n << "\t localRad = " << localRad << " \tGraft took " << niceTime(times.back())
				<< ", Graft within locality took " << niceTime(timesSingle.back().back()) << " for " << b.size() - 1
				<< " TSVs\n\t\t Graft measure = " << graftMeasure
				<< "\n\t\t Graft within locality measure = " << graftSingleMeasure
				<< endl;
			if(doDraw){
				drawTree(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_r" + makeString<double>(localRad, 2) + "_graftSingle.png", bGraftSingle, width, height, border);
			}
		}
	}

	cout << "\n\n localRad\t <time>\t +/-" << endl;
	int localRadIndex(0);
	for(double localRad(localRadStart); localRad < localRadEnd + localRadInc/2.0; localRad += localRadInc){
		vector<double> relaxTimes;
		for(int n(0); n < numNetworks; n++)
			relaxTimes.push_back(timesSingle[n][localRadIndex]);
		cout << localRad << "\t" << mean(relaxTimes) << "\t" << stDev(relaxTimes) << endl;
		localRadIndex++;
	}

	cout << "\n\n <TSV> =\t" << mean(terminalServiceVolumes) << "\t +/t\t" << stDev(terminalServiceVolumes)
		<< "\n\n n\t no_loc_meas\t <meas>\t +/-\t no_loc_time\t time\t +/-" << endl;
	for(int n(0); n < numNetworks; n++){
		cout << n << "\t" << measures[n] << "\t" << mean(measuresSingle[n]) << "\t" << stDev(measuresSingle[n])
			<< "\t" << times[n] << "\t" << mean(timesSingle[n]) << "\t" << stDev(timesSingle[n])
			<< endl;
	}
}

void testRegraftSingle(){
	kisset (14297, 19881, 20173);
	double side(2.0), nearestRadius(1.0);
	int width(200), height(width), border(width/20);
	vector<Bifurcation2D> b(initializeEvenRandomSpacing(10000, side, side, nearestRadius, side/2.0, nearestRadius));
	introduceBifurcationsFromSources(b);
	vector<int> ps, cs;
	for(unsigned int i(1); i < b.size(); i++){
		cs.push_back(i);
		if(b[i].childIndex.size() > 0)
			ps.push_back(i);
	}
	drawBifurcationTree("outputs/original.png", b, width, height, border);
	cout << "\n original tree = " << makeStringTreeOneLine(b) << "\n full original tree:\n" << makeStringTree(b) << endl;
	for(unsigned int p(0); p < ps.size(); p++){
		for(unsigned int c(0); c < cs.size(); c++){
			cout << "\n\n attemping ps[p] = " << ps[p] << " and cs[c] = " << cs[c] << endl;
			vector<Bifurcation2D> bCopy(b);
			if(regraftSingle(bCopy, ps[p], cs[c])){
				settleBifurcations(bCopy);
				drawBifurcationTree("outputs/regraft_p" + makeString(ps[p]) + "_c" + makeString(cs[c]) + "_t.png", bCopy, width, height, border);
				cout << "\n ps[p] = " << ps[p] << "  cs[c] = " << cs[c] << " keep first tree = " << makeStringTreeOneLine(bCopy) << endl;
				cout << "FULL TREE: " << makeStringTree(bCopy) << endl;
				cout << "\t measure = " << tripletMeasure(bCopy, 1.0, 1.0, 1.0) << endl;
				bCopy = b;
			}else
				cout << "\n ps[p] = " << ps[p] << "  cs[c] = " << cs[c] << " keep first SKIP" ;
			if(regraftSingle(bCopy, ps[p], cs[c], false)){
				settleBifurcations(bCopy);
				drawBifurcationTree("outputs/regraft_p" +  makeString(ps[p]) + "_c" + makeString(cs[c]) + "_f.png", bCopy, width, height, border);
				cout << "\n ps[p] = " << ps[p] << "  cs[c] = " << cs[c] << " keep second tree = " << makeStringTreeOneLine(bCopy) << endl;
				cout << "FULL TREE: " << makeStringTree(bCopy) << endl;
				cout << "\t measure = " << tripletMeasure(bCopy, 1.0, 1.0, 1.0) << endl;
			}else
				cout << "\n ps[p] = " << ps[p] << "  cs[c] = " << cs[c] << " keep second SKIP" ;
		}
	}
}

void localSingleGraftDebug(){
	kisset (9267, 2654, 4259);
	string uniquer("outputs/localSingleGraftDebug_c195_");
	bool doDraw(true);
	double area(4.0),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(9.0),
		maxPathLengthCoefficient(5.0),
		localRadStart(2.0*nearestRadius), localRadInc(nearestRadius), localRadEnd(2.0*nearestRadius);
	int numNetworks(1), numOpt(1), graftPasses(100), settleIts(1000),
		maxFailedInsertions(10000), width(500), height(width), border(width/20);
	cout << "\n uniquer = " << uniquer
		<< "\n area = " << area
		<< "\n settleThresh = " << settleThresh << "\n consolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;

	double side(sqrt(area)), heartx(side/2.0), hearty(nearestRadius/2.0), localRad(1.0);
	cout << "\n Creating b..." << endl;
	vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, side, side, nearestRadius, heartx, hearty));
	cout << "\n Optimizing..." << endl;
	vector<BranchPoint2D> bGraftSingle = optimizeByGlobalConsolidatedGraftingNEW(b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, graftPasses, settleIts);
	cout << "\n DONE!" << endl;
}

void makeStringTest(){
	for(double d(0.001); d < 2.0; d += 0.008)
		cout << d << " gives " << makeString(d, 2) << endl;
}

void asymArea(){
	int numNetworks(1), maxConsecutiveFailedInsertions(10000), width(1000), height(width), border(10),
		numOpt(10), maxPasses(200), settleIts(2000);
	double minX(1.0), minY(minX), maxX(2.0*minX), maxY(2.0*minY), nearestRadius(0.5/sqrt(acos(-1.0))),
		totalLengthCoefficient(1.0), avePathLengthCoefficient(9.0), maxPathLengthCoefficient(0.0);
	double consolidateThresh(0.01), settleThresh(0.001), sepThresh(0.25*nearestRadius), heartx(nearestRadius), hearty((minY + maxY)/2.0);

	cout << "\n numNetworks = " << numNetworks << "\n maxConsecutiveFailedInsertions = " << maxConsecutiveFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border << "\n numOpt = " << numOpt
		<< "\n maxPasses = " << maxPasses << "\n settleIts = " << settleIts
		<< "\n minX = " << minX << "\n minY = " << minY << "\n maxX = " << maxX << "\n maxY = " << maxY
		<< "\n nearestRadius = " << nearestRadius << "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient << "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n consolidateThresh = " << consolidateThresh << "\n settleThresh = " << settleThresh << "\n sepThresh = " << sepThresh
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< endl;

	for(int n(0); n < numNetworks; n++){
		vector<Bifurcation2D> b(initializeEvenRandomSpacingSingleBend(maxConsecutiveFailedInsertions, minX, minY, maxX, maxY, nearestRadius, heartx, hearty));
		drawBentBifurcationTree("outputs/asymArea_" + makeString(n) + ".png", b, minX, minY, width, height, border);
		vector<Bifurcation2D> bCopy(b);
		introduceBifurcationsFromSources(bCopy);
		settleBifurcations(bCopy, false, settleThresh);
		drawBentBifurcationTree("outputs/asymArea_" + makeString(n) + "_built.png", bCopy, minX, minY, width, height, border);
		vector<BranchPoint2D> br(consolidateDegenerates(bCopy, consolidateThresh));
		drawBentTree("outputs/asymArea_" + makeString(n) + "_consol.png", br, minX, minY, width, height, border);
		vector<BranchPoint2D> bo(optimizeByGlobalConsolidatedGraftingWithSingle(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, sepThresh, maxPasses, settleIts));
		drawBentTree("outputs/asymArea_" + makeString(n) + "_opt.png", bo, minX, minY, width, height, border);
	}
}

void regularLattice(){
	int numOpt(10), width(500), height(500), border(5), maxPasses(100), settleIts(10000);
	double side(15.0), nearestRadius(1.0/sqrt(acos(-1.0)));
	double consolidateThresh(0.01), settleThresh(0.001), sepThresh(0.25*nearestRadius),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(9.0),
		maxPathLengthCoefficient(0.0);
	vector<Bifurcation2D> b;
	for(double x(0.0); x < side + 0.5; x += 1.0){
		for(double y(0.0); y < side + 0.5; y += 1.0){
			b.push_back(Bifurcation2D());
			b.back().position[0] = x;
			b.back().position[1] = y;
		}
	}
	drawBifurcationTree("outputs/symArea.png", b, width, height, border);
	vector<Bifurcation2D> bCopy(b);
	introduceBifurcationsFromSources(bCopy);
	settleBifurcations(bCopy, false, settleThresh);
	drawBifurcationTree("outputs/symArea_built.png", bCopy, width, height, border);
	vector<BranchPoint2D> br(consolidateDegenerates(bCopy, consolidateThresh));
	drawTree("outputs/symArea_consol.png", br, width, height, border);
	vector<BranchPoint2D> bo(optimizeByGlobalConsolidatedGraftingWithSingle(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, sepThresh, maxPasses, settleIts));
	drawTree("outputs/symArea_opt.png", bo, width, height, border);
}

void circularSym(){
	int numOpt(5), maxConsecutiveFailedInsertions(10000), width(500), height(500), border(width/25), maxPasses(100), settleIts(10000);
	double side(5.0), nearestRadius(1.0/sqrt(acos(-1.0)));
	double consolidateThresh(0.01), settleThresh(0.001), sepThresh(0.5*nearestRadius),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0);
	vector<Bifurcation2D> b(initializeEvenRandomSpacingCircle(maxConsecutiveFailedInsertions, side, nearestRadius));
	drawBifurcationTree("outputs/circArea.png", b, width, height, border);
	vector<Bifurcation2D> bCopy(b);
	introduceBifurcationsFromSources(bCopy);
	settleBifurcations(bCopy, false, settleThresh);
	drawBifurcationTree("outputs/circArea_built.png", bCopy, width, height, border);
	vector<BranchPoint2D> br(consolidateDegenerates(bCopy, consolidateThresh));
	drawTreeCirc("outputs/circArea_consol.png", br, width, height, border);
	vector<BranchPoint2D> bo(optimizeByGlobalConsolidatedGraftingWithSingle(numOpt, b, 1.0, 0.0, 0.0, settleThresh, consolidateThresh, sepThresh, maxPasses, settleIts));
	drawTreeCirc("outputs/circArea_opt100.png", bo, width, height, border);
	vector<BranchPoint2D> bo2(optimizeByGlobalConsolidatedGraftingWithSingle(numOpt, b, 1.0, 9.0, 0.0, settleThresh, consolidateThresh, sepThresh, maxPasses, settleIts));
	drawTreeCirc("outputs/circArea_opt190.png", bo, width, height, border);
}

void lambdaLGammaOptCirc(){
	kisset (7446, 13195, 32375);
	string uniquer("lambdaLGammaOptOneNetCirc_c20_");
	bool doDraw(true);
	// area here is diameter
	double areaStart(20.0), areaInc(areaStart), areaEnd(areaStart),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0),
		localRad(2.0*nearestRadius);
	int numNetworks(100), numOpt(10), graftPasses(200), settleIts(10000), numBins(50),
		maxFailedInsertions(10000), width(1000), height(1000), border(1000/25);
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		<< "\n note that area is really halved, and then circled"
		<< "\n areaStart = " << areaStart << "\n areaInc = " << areaInc << "\n areaEnd = " << areaEnd
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n localRad = " << localRad
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses << "\n numBins = " << numBins
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<vector<double> > times, measures;
	vector<vector<vector<double> > >lambdaLs, gammas, tipRatios;
	vector<vector<unsigned int> > terminalServiceVolumes;
	vector<vector<int> > ignoreCounts, ignoreCountsGamma;
	for(double area(areaStart); area < areaEnd + areaInc/2.0; area += areaInc){
		double side(sqrt(area/2.0));//side(sqrt(area))
		times.push_back(vector<double>());
		measures.push_back(vector<double>());
		lambdaLs.push_back(vector<vector<double> >());
		gammas.push_back(vector<vector<double> >());
		tipRatios.push_back(vector<vector<double> >());
		terminalServiceVolumes.push_back(vector<unsigned int>());
		ignoreCounts.push_back(vector<int>());
		ignoreCountsGamma.push_back(vector<int>());
		for(int n(0); n < numNetworks; n++){
			vector<Bifurcation2D> b(initializeEvenRandomSpacingCircle(maxFailedInsertions, side, nearestRadius));
			time_t tempStart(time(NULL));
			vector<Bifurcation2D> bCopy(b);
			introduceBifurcationsFromSources(bCopy);
			vector<BranchPoint2D> bBuilt(consolidateDegenerates(bCopy));
			if(doDraw)
				drawTreeCirc(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_built.png", bBuilt, width, height, border);
			vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingle(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, graftPasses, settleIts);
			double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			times.back().push_back(difftime(time(NULL), tempStart));
			measures.back().push_back(graftMeasure);
			terminalServiceVolumes.back().push_back(b.size() - 1);
			cout << "\n\t\t Graft took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = " << graftMeasure << endl;
			if(doDraw)
				drawTreeCirc(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_graft.png", bGraft, width, height, border);
			ignoreCounts.back().push_back(0);
			ignoreCountsGamma.back().push_back(0);
			lambdaLs.back().push_back(validBranchingLengthRatios(bGraft, ignoreCounts.back().back(), lengthThresh));
			gammas.back().push_back(validParentChildLengthRatios(bGraft, ignoreCountsGamma.back().back(), lengthThresh));
			tipRatios.back().push_back(validTipCountRatios(bGraft, ignoreCounts.back().back(), lengthThresh)); // ignoreCounts should be the same as for lambdaL
		}
	}
	cout << "\n\n area\t <TSV>\t +/-\t fracIgn_rat\t +/-\t fracIgn_gam\t +/-\t measure\t +/-\t <time>\t +/-" << endl;
	for(unsigned int a(0); a < times.size(); a++){
		vector<double> ignoreFracs, ignoreFracsGamma;
		for(int n(0); n < numNetworks; n++){
			ignoreFracs.push_back(double(ignoreCounts[a][n])/double(terminalServiceVolumes[a][n] - 1));
			ignoreFracsGamma.push_back(double(ignoreCountsGamma[a][n])/double(2*terminalServiceVolumes[a][n] - 2));
		}
		cout << areaStart + a*areaInc
			<< "\t" << mean(terminalServiceVolumes[a]) << "\t" << stDev(terminalServiceVolumes[a])
			<< "\t" << mean(ignoreFracs) << "\t" << stDev(ignoreFracs)
			<< "\t" << mean(ignoreFracsGamma) << "\t" << stDev(ignoreFracsGamma)
			<< "\t" << mean(measures[a]) << "\t" << stDev(measures[a])
			<< "\t" << mean(times[a]) << "\t" << stDev(times[a])
			<< endl;
	}

	vector<vector<map<double, double> > > lambdaL, gamma, tipRatio; // single networks
	vector<map<double, double> > areaLambdaL, areaGamma, areaTipRatio; // within area, over networks
	vector<double> allLambdaLs, allGammas, allTipRatios; // over all areas
	for(unsigned int a(0); a < lambdaLs.size(); a++){
		lambdaL.push_back(vector<map<double, double> >());
		gamma.push_back(vector<map<double, double> >());
		tipRatio.push_back(vector<map<double, double> >());
		vector<double> areaLambdaLs, areaGammas, areaTipRatios; // within area, over networks
		for(int n(0); n < numNetworks; n++){
			areaLambdaLs.insert(areaLambdaLs.end(), lambdaLs[a][n].begin(), lambdaLs[a][n].end());
			areaGammas.insert(areaGammas.end(), gammas[a][n].begin(), gammas[a][n].end());
			areaTipRatios.insert(areaTipRatios.end(), tipRatios[a][n].begin(), tipRatios[a][n].end());
			lambdaL.back().push_back(probabilityDistribution(lambdaLs[a][n], numBins));
			gamma.back().push_back(probabilityDistribution(gammas[a][n], numBins));
			tipRatio.back().push_back(probabilityDistribution(tipRatios[a][n], numBins));
		}
		areaLambdaL.push_back(probabilityDistribution(areaLambdaLs, numBins));
		areaGamma.push_back(probabilityDistribution(areaGammas, numBins));
		areaTipRatio.push_back(probabilityDistribution(areaTipRatios, numBins));
		allLambdaLs.insert(allLambdaLs.end(), areaLambdaLs.begin(), areaLambdaLs.end());
		allGammas.insert(allGammas.end(), areaGammas.begin(), areaGammas.end());
		allTipRatios.insert(allTipRatios.end(), areaTipRatios.begin(), areaTipRatios.end());
	}
	map<double, double> allLambdaL(probabilityDistribution(allLambdaLs, numBins)),
		allGamma(probabilityDistribution(allGammas, numBins)),
		allTipRatio(probabilityDistribution(allTipRatios, numBins));

	printReport("lambdaL", numBins, lambdaL, areaLambdaL, allLambdaL, areaStart, areaInc, numNetworks);
	printReport("gamma", numBins, gamma, areaGamma, allGamma, areaStart, areaInc, numNetworks);
	printReport("tipRatio", numBins, tipRatio, areaTipRatio, allTipRatio, areaStart, areaInc, numNetworks);
}

void lambdaLOptCircConst(){
	kisset (7446, 13195, 32375);
	string uniquer("lambdaLOptCircConst_c100_");
	bool doDraw(true);
	// area here is diameter
	double areaStart(100.0), areaInc(20.0), areaEnd(100.0),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(5.0),
		localRad(2.0*nearestRadius),
		capConst(nearestRadius);
	int numNetworks(10), numOpt(5), graftPasses(200), settleIts(10000), numBins(27),
		maxFailedInsertions(10000), width(1000), height(1000), border(1000/25);
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		<< "\n note that area is really halved"
		<< "\n areaStart = " << areaStart << "\n areaInc = " << areaInc << "\n areaEnd = " << areaEnd
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n capConst  = " << capConst
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n localRad = " << localRad
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses << "\n numBins = " << numBins
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<vector<double> > times, measures;
	vector<vector<vector<double> > >lambdaLs;//, gammas, tipRatios;
	vector<vector<unsigned int> > terminalServiceVolumes;
	vector<vector<int> > ignoreCounts;//, ignoreCountsGamma;
	for(double area(areaStart); area < areaEnd + areaInc/2.0; area += areaInc){
		double side(sqrt(area/2.0));//side(sqrt(area))
		times.push_back(vector<double>());
		measures.push_back(vector<double>());
		lambdaLs.push_back(vector<vector<double> >());
		//gammas.push_back(vector<vector<double> >());
		//tipRatios.push_back(vector<vector<double> >());
		terminalServiceVolumes.push_back(vector<unsigned int>());
		ignoreCounts.push_back(vector<int>());
		//ignoreCountsGamma.push_back(vector<int>());
		for(int n(0); n < numNetworks; n++){
			vector<Bifurcation2D> b(initializeEvenRandomSpacingCircle(maxFailedInsertions, side, nearestRadius));
			time_t tempStart(time(NULL));
			vector<Bifurcation2D> bCopy(b);
			introduceBifurcationsFromSources(bCopy);
			vector<BranchPoint2D> bBuilt(consolidateDegenerates(bCopy));
			if(doDraw)
				drawTreeCirc(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_built.png", bBuilt, width, height, border);
			vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingle(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, graftPasses, settleIts);
			double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			times.back().push_back(difftime(time(NULL), tempStart));
			measures.back().push_back(graftMeasure);
			terminalServiceVolumes.back().push_back(b.size() - 1);
			cout << "\n\t\t Graft took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = " << graftMeasure << endl;
			if(doDraw)
				drawTreeCirc(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_graft.png", bGraft, width, height, border);
			ignoreCounts.back().push_back(0);
			//ignoreCountsGamma.back().push_back(0);
			lambdaLs.back().push_back(validBranchingLengthRatiosWithConstant(bGraft, ignoreCounts.back().back(), capConst, lengthThresh));
			//gammas.back().push_back(validParentChildLengthRatios(bGraft, ignoreCountsGamma.back().back(), lengthThresh));
			//tipRatios.back().push_back(validTipCountRatios(bGraft, ignoreCounts.back().back(), lengthThresh)); // ignoreCounts should be the same as for lambdaL
		}
	}
	cout << "\n\n area\t <TSV>\t +/-\t fracIgn_rat\t +/-\t fracIgn_gam\t +/-\t measure\t +/-\t <time>\t +/-" << endl;
	for(unsigned int a(0); a < times.size(); a++){
		vector<double> ignoreFracs;//, ignoreFracsGamma;
		for(int n(0); n < numNetworks; n++){
			ignoreFracs.push_back(double(ignoreCounts[a][n])/double(terminalServiceVolumes[a][n] - 1));
			//ignoreFracsGamma.push_back(double(ignoreCountsGamma[a][n])/double(2*terminalServiceVolumes[a][n] - 2));
		}
		cout << areaStart + a*areaInc
			<< "\t" << mean(terminalServiceVolumes[a]) << "\t" << stDev(terminalServiceVolumes[a])
			<< "\t" << mean(ignoreFracs) << "\t" << stDev(ignoreFracs)
			//<< "\t" << mean(ignoreFracsGamma) << "\t" << stDev(ignoreFracsGamma)
			<< "\t" << mean(measures[a]) << "\t" << stDev(measures[a])
			<< "\t" << mean(times[a]) << "\t" << stDev(times[a])
			<< endl;
	}

	vector<vector<map<double, double> > > lambdaL;//, gamma, tipRatio; // single networks
	vector<map<double, double> > areaLambdaL;//, areaGamma, areaTipRatio; // within area, over networks
	vector<double> allLambdaLs;//, allGammas, allTipRatios; // over all areas
	for(unsigned int a(0); a < lambdaLs.size(); a++){
		lambdaL.push_back(vector<map<double, double> >());
		//gamma.push_back(vector<map<double, double> >());
		//tipRatio.push_back(vector<map<double, double> >());
		vector<double> areaLambdaLs;//, areaGammas, areaTipRatios; // within area, over networks
		for(int n(0); n < numNetworks; n++){
			areaLambdaLs.insert(areaLambdaLs.end(), lambdaLs[a][n].begin(), lambdaLs[a][n].end());
			//areaGammas.insert(areaGammas.end(), gammas[a][n].begin(), gammas[a][n].end());
			//areaTipRatios.insert(areaTipRatios.end(), tipRatios[a][n].begin(), tipRatios[a][n].end());
			lambdaL.back().push_back(probabilityDistribution(lambdaLs[a][n], numBins));
			//gamma.back().push_back(probabilityDistribution(gammas[a][n], numBins));
			//tipRatio.back().push_back(probabilityDistribution(tipRatios[a][n], numBins));
		}
		areaLambdaL.push_back(probabilityDistribution(areaLambdaLs, numBins));
		//areaGamma.push_back(probabilityDistribution(areaGammas, numBins));
		//areaTipRatio.push_back(probabilityDistribution(areaTipRatios, numBins));
		allLambdaLs.insert(allLambdaLs.end(), areaLambdaLs.begin(), areaLambdaLs.end());
		//allGammas.insert(allGammas.end(), areaGammas.begin(), areaGammas.end());
		//allTipRatios.insert(allTipRatios.end(), areaTipRatios.begin(), areaTipRatios.end());
	}
	map<double, double> allLambdaL(probabilityDistribution(allLambdaLs, numBins));//,
		//allGamma(probabilityDistribution(allGammas, numBins)),
		//allTipRatio(probabilityDistribution(allTipRatios, numBins));

	printReport("lambdaL", numBins, lambdaL, areaLambdaL, allLambdaL, areaStart, areaInc, numNetworks);
}

void lambdaLOptCircIgnoreTips(){
	//kisset (7446, 13195, 32375);
	string uniquer("lambdaLOptCircIgTips_c100_");
	bool doDraw(true);
	// area here is diameter
	double areaStart(100.0), areaInc(20.0), areaEnd(100.0),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0),
		localRad(2.0*nearestRadius),
		capConst(nearestRadius);
	int numNetworks(10), numOpt(5), graftPasses(200), settleIts(10000), numBins(27),
		maxFailedInsertions(10000), width(1000), height(1000), border(1000/25);
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		<< "\n note that area is really halved"
		<< "\n areaStart = " << areaStart << "\n areaInc = " << areaInc << "\n areaEnd = " << areaEnd
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n capConst  = " << capConst
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n localRad = " << localRad
		<< "\n numNetworks = " << numNetworks << "\n graftPasses = " << graftPasses << "\n numBins = " << numBins
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<vector<double> > times, measures;
	vector<vector<vector<double> > >lambdaLs;//, gammas, tipRatios;
	vector<vector<unsigned int> > terminalServiceVolumes;
	vector<vector<int> > ignoreCounts;//, ignoreCountsGamma;
	for(double area(areaStart); area < areaEnd + areaInc/2.0; area += areaInc){
		double side(sqrt(area/2.0));//side(sqrt(area))
		times.push_back(vector<double>());
		measures.push_back(vector<double>());
		lambdaLs.push_back(vector<vector<double> >());
		//gammas.push_back(vector<vector<double> >());
		//tipRatios.push_back(vector<vector<double> >());
		terminalServiceVolumes.push_back(vector<unsigned int>());
		ignoreCounts.push_back(vector<int>());
		//ignoreCountsGamma.push_back(vector<int>());
		for(int n(0); n < numNetworks; n++){
			vector<Bifurcation2D> b(initializeEvenRandomSpacingCircle(maxFailedInsertions, side, nearestRadius));
			time_t tempStart(time(NULL));
			vector<Bifurcation2D> bCopy(b);
			introduceBifurcationsFromSources(bCopy);
			vector<BranchPoint2D> bBuilt(consolidateDegenerates(bCopy));
			if(doDraw)
				drawTreeCirc(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_built.png", bBuilt, width, height, border);
			vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingle(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, graftPasses, settleIts);
			double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
			times.back().push_back(difftime(time(NULL), tempStart));
			measures.back().push_back(graftMeasure);
			terminalServiceVolumes.back().push_back(b.size() - 1);
			cout << "\n\t\t Graft took " << niceTime(times.back().back()) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = " << graftMeasure << endl;
			if(doDraw)
				drawTreeCirc(uniquer + makeString((int)ceil(area)) + "_n" + makeString(n) + "_graft.png", bGraft, width, height, border);
			ignoreCounts.back().push_back(0);
			//ignoreCountsGamma.back().push_back(0);
			lambdaLs.back().push_back(validBranchingLengthRatiosIgnoreTips(bGraft, ignoreCounts.back().back(), lengthThresh));
			//gammas.back().push_back(validParentChildLengthRatios(bGraft, ignoreCountsGamma.back().back(), lengthThresh));
			//tipRatios.back().push_back(validTipCountRatios(bGraft, ignoreCounts.back().back(), lengthThresh)); // ignoreCounts should be the same as for lambdaL
		}
	}
	cout << "\n\n area\t <TSV>\t +/-\t fracIgn_rat\t +/-\t fracIgn_gam\t +/-\t measure\t +/-\t <time>\t +/-" << endl;
	for(unsigned int a(0); a < times.size(); a++){
		vector<double> ignoreFracs;//, ignoreFracsGamma;
		for(int n(0); n < numNetworks; n++){
			ignoreFracs.push_back(double(ignoreCounts[a][n])/double(terminalServiceVolumes[a][n] - 1));
			//ignoreFracsGamma.push_back(double(ignoreCountsGamma[a][n])/double(2*terminalServiceVolumes[a][n] - 2));
		}
		cout << areaStart + a*areaInc
			<< "\t" << mean(terminalServiceVolumes[a]) << "\t" << stDev(terminalServiceVolumes[a])
			<< "\t" << mean(ignoreFracs) << "\t" << stDev(ignoreFracs)
			//<< "\t" << mean(ignoreFracsGamma) << "\t" << stDev(ignoreFracsGamma)
			<< "\t" << mean(measures[a]) << "\t" << stDev(measures[a])
			<< "\t" << mean(times[a]) << "\t" << stDev(times[a])
			<< endl;
	}

	vector<vector<map<double, double> > > lambdaL;//, gamma, tipRatio; // single networks
	vector<map<double, double> > areaLambdaL;//, areaGamma, areaTipRatio; // within area, over networks
	vector<double> allLambdaLs;//, allGammas, allTipRatios; // over all areas
	for(unsigned int a(0); a < lambdaLs.size(); a++){
		lambdaL.push_back(vector<map<double, double> >());
		//gamma.push_back(vector<map<double, double> >());
		//tipRatio.push_back(vector<map<double, double> >());
		vector<double> areaLambdaLs;//, areaGammas, areaTipRatios; // within area, over networks
		for(int n(0); n < numNetworks; n++){
			areaLambdaLs.insert(areaLambdaLs.end(), lambdaLs[a][n].begin(), lambdaLs[a][n].end());
			//areaGammas.insert(areaGammas.end(), gammas[a][n].begin(), gammas[a][n].end());
			//areaTipRatios.insert(areaTipRatios.end(), tipRatios[a][n].begin(), tipRatios[a][n].end());
			lambdaL.back().push_back(probabilityDistribution(lambdaLs[a][n], numBins));
			//gamma.back().push_back(probabilityDistribution(gammas[a][n], numBins));
			//tipRatio.back().push_back(probabilityDistribution(tipRatios[a][n], numBins));
		}
		areaLambdaL.push_back(probabilityDistribution(areaLambdaLs, numBins));
		//areaGamma.push_back(probabilityDistribution(areaGammas, numBins));
		//areaTipRatio.push_back(probabilityDistribution(areaTipRatios, numBins));
		allLambdaLs.insert(allLambdaLs.end(), areaLambdaLs.begin(), areaLambdaLs.end());
		//allGammas.insert(allGammas.end(), areaGammas.begin(), areaGammas.end());
		//allTipRatios.insert(allTipRatios.end(), areaTipRatios.begin(), areaTipRatios.end());
	}
	map<double, double> allLambdaL(probabilityDistribution(allLambdaLs, numBins));//,
		//allGamma(probabilityDistribution(allGammas, numBins)),
		//allTipRatio(probabilityDistribution(allTipRatios, numBins));

	printReport("lambdaL", numBins, lambdaL, areaLambdaL, allLambdaL, areaStart, areaInc, numNetworks);
}

void settleFig(){
	overrelaxationFactor = 0.0;
	kisset(1, 2, 3);
	double side(3.0), nearestRadius(2.0/acos(-1.0));
	int width(1000), height(width), border(width/25);
	vector<Bifurcation2D> b;//(initializeEvenRandomSpacing(10000, side, side, nearestRadius, 0.0, 0.0));
	/*drawBifurcationTree("outputs/dist.png", b, width, height, border);
	introduceBifurcationsFromSources(b);
	vector<Bifurcation2D> bOpt();
	bool changed(true);
	while(changed){
		changed = false;
		for(unsigned int i(0); i < b.size(); i++){
			if(b[i].childIndex.size() == 2){
				double x(b[i].position[0]), y(b[i].position[1]);
				b[i].position[0] = (b[b[i].childIndex[0]].position[0] + b[b[i].childIndex[1]].position[0])/2.0;
				b[i].position[1] = (b[b[i].childIndex[0]].position[1] + b[b[i].childIndex[1]].position[1])/2.0;
				changed = changed || x != b[i].position[0] || y!= b[i].position[1];
			}
		}
	}*/
	
	/*b.push_back(Bifurcation2D()); b.back().position[0] = 0.0; b.back().position[1] = 0.0;
	b.push_back(Bifurcation2D()); b.back().position[0] = 0.8; b.back().position[1] = 0.0;
	b.push_back(Bifurcation2D()); b.back().position[0] = 1.0; b.back().position[1] = 0.8;
	b.push_back(Bifurcation2D()); b.back().position[0] = 0.7; b.back().position[1] = 1.0;
	b.push_back(Bifurcation2D()); b.back().position[0] = 0.3; b.back().position[1] = 0.9;
	double lineFrac(0.9);
	double x(b[1].position[0] + lineFrac*(b[3].position[0] - b[1].position[0])), y(b[1].position[1] + lineFrac*(b[3].position[1] - b[1].position[1]));
	addParent(b, 1, 3, x, y);
	addParent(b, 5, 2);
	addParent(b, 4, 6);
	b[7].parentIndex = 0;
	b[0].childIndex.push_back(7);*/
	b.push_back(Bifurcation2D()); b.back().position[0] = 0.5; b.back().position[1] = 0.0;
	b.push_back(Bifurcation2D()); b.back().position[0] = 0; b.back().position[1] = 1.5;
	b.push_back(Bifurcation2D()); b.back().position[0] = 1.0; b.back().position[1] = 0.7;
	addParent(b, 1, 2);
	b[3].parentIndex = 0;
	b[0].childIndex.push_back(3);
	cout << "\n Fermat initial total length = " << tripletMeasure(b, 1.0, 0.0, 0.0);
	bool drawSettleFermat(false);
	settleBifurcations(b, drawSettleFermat, 0.001, 100);
	cout << "\n Fermat end total length = " << tripletMeasure(b, 1.0, 0.0, 0.0);
	drawTree("outputs/Fermat.png", consolidateDegenerates(b, 0.01), width, height, border);
	b[3].position[0] = 0.5;
	b[3].position[1] = 0.4;
	drawTree("outputs/wrongFermat.png", consolidateDegenerates(b, 0.01), width, height, border);
	cout << "\n wrongFermat total length = " << tripletMeasure(b, 1.0, 0.0, 0.0);
	
	b.clear();
	b.push_back(Bifurcation2D()); b.back().position[0] = 0.5; b.back().position[1] = 0.0;
	b.push_back(Bifurcation2D()); b.back().position[0] = 1.0; b.back().position[1] = 0.4;
	b.push_back(Bifurcation2D()); b.back().position[0] = 0.5; b.back().position[1] = 1.0;
	b.push_back(Bifurcation2D()); b.back().position[0] = 0.0; b.back().position[1] = 0.6;
	addParent(b, 2, 3);
	addParent(b, 1, 4);
	b[5].parentIndex = 0;
	b[0].childIndex.push_back(5);
	bool drawSettleGeomMedPrevious(false);
	settleBifurcations(b, drawSettleGeomMedPrevious, 0.001, 100);
	cout << "\n\n total length = " << totalLength(b)
		<< "\n average height = " << averageDistanceFromHeart(b)
		<< "\n max height = " << maxDistanceFromHeart(b)
		<< endl;
	drawTree("outputs/previous.png", consolidateDegenerates(b, 0.01), width, height, border);

	if(!regraftSingle(b, 4, 5, true))
		cout << "\n Could not regraft" << endl;
	drawTree("outputs/initial.png", consolidateDegenerates(b, 0.01), width, height, border);
	bool drawSettleGeomMed(true);
	settleBifurcations(b, drawSettleGeomMed, 0.001, 20);
	cout << "\n\n total length = " << totalLength(b)
		<< "\n average height = " << averageDistanceFromHeart(b)
		<< "\n max height = " << maxDistanceFromHeart(b)
		<< endl;
	drawTree("outputs/final.png", consolidateDegenerates(b, 0.01), width, height, border);

	if(!regraftSingle(b, 4, 1, true))
		cout << "\n Could not regraft" << endl;
	bool drawSettleGeomMedNext(false);
	settleBifurcations(b, drawSettleGeomMedNext, 0.001, 20);
	drawTree("outputs/third.png", consolidateDegenerates(b, 0.01), width, height, border);
	cout << "\n\n total length = " << totalLength(b)
		<< "\n average height = " << averageDistanceFromHeart(b)
		<< "\n max height = " << maxDistanceFromHeart(b)
		<< endl;

	if(!regraftSingle(b, 4, 3, true))
		cout << "\n Could not regraft" << endl;
	b[4].position[0] = 0.4; b[4].position[1] = 0.6;
	b[5].position[0] = 0.67; b[5].position[1] = 0.18;
	drawTree("outputs/impossibleGeomMed1.png", consolidateDegenerates(b, 0.01), width, height, border);

	b[4].position[0] = 0.5; b[4].position[1] = 0.4;
	b[5].position[0] = 0.65; b[5].position[1] = 0.25;
	drawTree("outputs/impossibleGeomMed2.png", consolidateDegenerates(b, 0.01), width, height, border);

	b[4].position[0] = 0.35; b[4].position[1] = 0.6;
	b[5].position[0] = 0.6; b[5].position[1] = 0.2;
	drawTree("outputs/impossibleGeomMed3.png", consolidateDegenerates(b, 0.01), width, height, border);

	b[4].position[0] = 0.4; b[4].position[1] = 0.55;
	b[5].position[0] = 0.55; b[5].position[1] = 0.4;
	drawTree("outputs/impossibleGeomMed4.png", consolidateDegenerates(b, 0.01), width, height, border);
	
}

void testSimilarity(){
	//kisset(24611, 19046, 7706);
	//kisset (24678, 4714, 15707); // takes a long time for 5
	kisset (5372, 4528, 29834); // only find 1 configuration for 5 TSVs because of single central hub
	//kisset (8855, 30979, 16034);
	int maxFailedInsertions(100000), targetServiceVolumes(6), numBest(945);
	bool doDraw(false);
	double nearestRadius(1.0),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0),
		settleThresh(0.001), consolidateThresh(0.01);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	string uniquer("testSimilarity_" + makeString(sizeX));
	int width((int)floor(sizeX)*10), height((int)floor(sizeY)*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300/2;
	border = 10/2;

	cout << "\n testSimilarity()\n uniquer = " << uniquer
		<< "\n sizeX = " << sizeX << "\n sizeY = " << sizeY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n numBest = " << numBest
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient" << maxPathLengthCoefficient
		<< "\n settleThresh = " << settleThresh
		<< "\n consolidateThresh = " << consolidateThresh
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< endl;
	vector<Bifurcation2D> b(initializeEvenRandomSpacingServiceVolumes(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty, targetServiceVolumes, maxFailedInsertions));
	drawBifurcationTree("outputs/" + uniquer + "_initial.png", b, width, height, border);
	vector<BranchPoint2D> *bestBs = new vector<BranchPoint2D>[numBest];
	double *bestBMeasures = new double[numBest];
	for(int i(0); i < numBest; i++)
		bestBMeasures[i] = -1.0;
	exhaustiveConsolidatedHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, true, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, true);
	
	int worstIndex(0);
	double *sims = new double[numBest];
	for(int i(0); i < numBest; i++){
		if(bestBMeasures[i] < 0.0)
			sims[i] = -1.0;
		else{
			sims[i] = subtreeSimilarity(bestBs[0], bestBs[i]);
			worstIndex = i;
		}
	}

	cout << "\n\n RANKED BY MEASURE\n bestBs\t measure\t rr\t normalized_measure\t brSim" << endl;
	for(int i(0); i < numBest; i++){
		if(bestBMeasures[i] >= 0.0){
			cout << makeStringTreeOneLine(bestBs[i]) << "\t" << bestBMeasures[i]
			<< "\t" << double(i)/double(worstIndex)
			<< "\t" << bestBMeasures[i]/bestBMeasures[0]
			<< "\t" << sims[i]
			<< endl;
			if(doDraw)
				drawTree("outputs/" + uniquer + "_ranked" + makeString(i) + ".png", bestBs[i], width, height, border);
		}
	}

	//re-order based on similarity
	bool changed(true);
	while(changed){
		changed = false;
		for(int i(0); i < numBest - 1; i++){
			if(sims[i] < sims[i + 1]){
				swap(sims[i], sims[i + 1]);
				swap(bestBs[i], bestBs[i + 1]);
				swap(bestBMeasures[i], bestBMeasures[i + 1]);
			}
		}
	}
	
	cout << "\n\n RANKED BY SIMILARITY\n bestBs\t measure\t rr\t normalized_measure\t brSim" << endl;
	for(int i(0); i < numBest; i++){
		if(bestBMeasures[i] >= 0.0){
			cout << makeStringTreeOneLine(bestBs[i]) << "\t" << bestBMeasures[i]
			<< "\t" << double(i)/double(worstIndex)
			<< "\t" << bestBMeasures[i]/bestBMeasures[0]
			<< "\t" << sims[i]
			<< endl;
			//drawTree("outputs/" + uniquer + "_ranked" + makeString(i) + ".png", bestBs[i], width, height, border);
		}
	}

	delete[] bestBs;
	delete[] bestBMeasures;
	delete[] sims;
}

void bestGuessRank(){
	kisset (5372, 4528, 29834); // only find 1 configuration for 5 TSVs because of single central hub
	//kisset (8855, 30979, 16034);
	int maxFailedInsertions(100000), targetServiceVolumes(4), numBest(15), numNetworks(1000);
	double nearestRadius(1.0),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(9.0),
		maxPathLengthCoefficient(0.0),
		settleThresh(0.001), consolidateThresh(0.01);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	string uniquer("bestGuessRank_" + makeString(sizeX));
	int width((int)floor(sizeX)*10), height((int)floor(sizeY)*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300/2;
	border = 10/2;

	cout << "\n bestGuessRank()\n uniquer = " << uniquer
		<< "\n sizeX = " << sizeX << "\n sizeY = " << sizeY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n numBest = " << numBest << "\n numNetworks = " << numNetworks
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient" << maxPathLengthCoefficient
		<< "\n settleThresh = " << settleThresh
		<< "\n consolidateThresh = " << consolidateThresh
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< endl;
	vector<int> guessIndex, maxIndex;
	vector<double> rescaledGuess;
	cout << "\n\n guessIndex\t maxIndex\t guessRescaled" << endl;
	for(int n(0); n < numNetworks; n++){
		vector<Bifurcation2D> b(initializeEvenRandomSpacingServiceVolumes(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty, targetServiceVolumes, maxFailedInsertions));
		vector<Bifurcation2D> bGuess(b);
		introduceBifurcationsFromSources(bGuess);
		vector<BranchPoint2D> *bestBs = new vector<BranchPoint2D>[numBest];
		double *bestBMeasures = new double[numBest];
		for(int i(0); i < numBest; i++)
			bestBMeasures[i] = -1.0;
		exhaustiveConsolidatedHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, true);
		int matchingIndex(-1), mxIndex(-1);
		for(int i(0); i < numBest && bestBMeasures[i] >= 0.0; i++){
			if(equivalentTrees(consolidateDegenerates(bGuess, consolidateThresh), bestBs[i], settleThresh))
				matchingIndex = i;
			mxIndex = i;
		}
		guessIndex.push_back(matchingIndex);
		maxIndex.push_back(mxIndex);
		rescaledGuess.push_back(double(matchingIndex)/double(mxIndex));
		cout << matchingIndex << "\t" << mxIndex << "\t" << double(matchingIndex)/double(mxIndex) << endl;
		delete[] bestBs;
		delete[] bestBMeasures;
	}
	cout << "\n\n <guessIndex>  = " << mean(guessIndex) << " +/- " << stDev(guessIndex)
		<< "\n <maxIndex> = " << mean(maxIndex) << " +/- " << stDev(maxIndex)
		<< "\n <rescaledGuess> = " << mean(rescaledGuess) << " +/- " << stDev(rescaledGuess)
		<< endl;

}

void bestSymGuessRank(){
	kisset (5372, 4528, 29834); // only find 1 configuration for 5 TSVs because of single central hub
	//kisset (8855, 30979, 16034);
	int maxFailedInsertions(100000), targetServiceVolumes(4), numBest(15), numNetworks(1000);
	double nearestRadius(1.0),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0),
		settleThresh(1.0e-6), consolidateThresh(0.01);
	double sizeX(nearestRadius*sqrt((2 + targetServiceVolumes)*acos(-1.0))/2.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	string uniquer("bestGuessRank_" + makeString(sizeX));
	int width((int)floor(sizeX)*10), height((int)floor(sizeY)*100), border(width/25);
	if(height < width)
		border = height/25;
	width = height = 300/2;
	border = 10/2;

	cout << "\n bestGuessRank()\n uniquer = " << uniquer
		<< "\n sizeX = " << sizeX << "\n sizeY = " << sizeY
		<< "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n numBest = " << numBest << "\n numNetworks = " << numNetworks
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient" << maxPathLengthCoefficient
		<< "\n settleThresh = " << settleThresh
		<< "\n consolidateThresh = " << consolidateThresh
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n nearestRadius = " << nearestRadius
		<< endl;
	vector<int> guessIndex, maxIndex;
	vector<double> rescaledGuess;
	cout << "\n\n guessIndex\t maxIndex\t guessRescaled\t guessTree" << endl;
	for(int n(0); n < numNetworks; n++){
		vector<Bifurcation2D> b(initializeEvenRandomSpacingServiceVolumes(maxFailedInsertions, sizeX, sizeY, nearestRadius, heartx, hearty, targetServiceVolumes, maxFailedInsertions));
		vector<Bifurcation2D> bGuess(b);
		introduceBifurcationsSymmetricallyFromHeart(bGuess, settleThresh);
		vector<BranchPoint2D> *bestBs = new vector<BranchPoint2D>[numBest];
		double *bestBMeasures = new double[numBest];
		for(int i(0); i < numBest; i++)
			bestBMeasures[i] = -1.0;
		exhaustiveConsolidatedHierarchyConfigurationSearch(b, numBest, bestBs, bestBMeasures, false, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, true);
		int matchingIndex(-1), mxIndex(-1);
		for(int i(0); i < numBest && bestBMeasures[i] >= 0.0; i++){
			if(matchingIndex < 0 && equivalentTrees(consolidateDegenerates(bGuess, consolidateThresh), bestBs[i], settleThresh))
				matchingIndex = i;
			//else
			///	cout << "\tbGuess " << n << " is not " << makeStringTreeOneLine(bestBs[i]) << endl;
			mxIndex = i;
		}
		guessIndex.push_back(matchingIndex);
		maxIndex.push_back(mxIndex);
		rescaledGuess.push_back(double(matchingIndex)/double(mxIndex));
		cout << matchingIndex << "\t" << mxIndex << "\t" << double(matchingIndex)/double(mxIndex)  << "\t" << makeStringTreeOneLine(bGuess) << endl;
		delete[] bestBs;
		delete[] bestBMeasures;
	}
	cout << "\n\n <guessIndex>  = " << mean(guessIndex) << " +/- " << stDev(guessIndex)
		<< "\n <maxIndex> = " << mean(maxIndex) << " +/- " << stDev(maxIndex)
		<< "\n <rescaledGuess> = " << mean(rescaledGuess) << " +/- " << stDev(rescaledGuess)
		<< endl;

}

struct PathOfBranchPoint2D{
	vector<BranchPoint2D> *b; // pointer to vector with bifurcations
	vector<int> path; // indices in b in the path
	double length; // total length of path
	int parentIndex; // index of parent in the context of the vector of PathOfBranchPoint2Ds
	vector<int> childIndex; // indices of children in the context of the vector of PathOfBranchPoint2Ds
};

void printPathOfBranchPoint2DVector(const vector<PathOfBranchPoint2D> &p){
	streamsize prec(cout.precision());
	cout.precision(3);
	for(unsigned int i(0); i < p.size(); i++){
		cout << "\n\t path " << i
			<< "\t\t parentIndex = " << p[i].parentIndex
			<< "\t childIndex = " << makeVectorString(p[i].childIndex)
			<< "\n\t\t length = " << p[i].length
			<< endl;
		int lastInParent(-1);
		if(p[i].parentIndex >= 0)
			lastInParent = p[p[i].parentIndex].path.back();
		cout << "\t\t last in parent = " << lastInParent << endl;
		for(unsigned int j(0); j < p[i].path.size(); j++){
			cout << "\t\t branchPoint2D at path step " << j
				<< ": b[" << p[i].path[j] << "]\t parentIndex = " << (*(p[i].b))[p[i].path[j]].parentIndex
				<< "\t childIndex = " << makeVectorString((*(p[i].b))[p[i].path[j]].childIndex)
				<< "\t position = (" <<  (*(p[i].b))[p[i].path[j]].position[0] << ", " << (*(p[i].b))[p[i].path[j]].position[1] << ")"
				<< endl;
		}
	}
	cout.precision(prec);
}

// default i = 0 assumes heart is 0
// i is the index of the first bifurcation point in b after the start (i's parent)
// b has been trimmed of degenerate sinks
void buildPathOfBranchPoint2D(vector<PathOfBranchPoint2D> &p, vector<BranchPoint2D> *bOrig, vector<BranchPoint2D> &b, int i = 0, int pIndex = -1){
	p.push_back(PathOfBranchPoint2D());
	p.back().b = bOrig;
	if(b[i].parentIndex >= 0){
		p.back().path.push_back(b[i].parentIndex);
		p.back().length = separation2D(b[i], b[b[i].parentIndex]);
	}else
		p.back().length = 0.0;
	p.back().path.push_back(i);
	p.back().parentIndex = pIndex;
	while(b[i].childIndex.size() == 1){
		i = b[i].childIndex[0];
		p.back().path.push_back(i);
		p.back().length += separation2D(b[i], b[b[i].parentIndex]);
	}
	if(b[i].childIndex.size() > 1){
		int newPIndex(p.size() - 1);
		for(unsigned int c(0); c < b[i].childIndex.size(); c++){
			p[newPIndex].childIndex.push_back(p.size());
			buildPathOfBranchPoint2D(p, bOrig, b, b[i].childIndex[c], newPIndex);
		}
	}
}

// assumes heart is b[0]
vector<PathOfBranchPoint2D> trimDegenerateSinks(vector<BranchPoint2D> &bOrig, double consolidateThresh){
	vector<BranchPoint2D> bCopy(bOrig);
	for(unsigned int i(0); i < bCopy.size(); i++){
		if(bCopy[i].parentIndex < 0)
			continue;
		if(bCopy[i].childIndex.size() < 1 && separation2D(bCopy[i], bCopy[bCopy[i].parentIndex]) < consolidateThresh)
			removeFromVector(i, bCopy[bCopy[i].parentIndex].childIndex);
	}
	vector<PathOfBranchPoint2D> p;
	buildPathOfBranchPoint2D(p, &bOrig, bCopy);
	return p;
}

void addBranchingLengthRatios(vector<double> &lambdaLs, const vector<PathOfBranchPoint2D> &p, int pIndex = 0){
	if(p[pIndex].childIndex.size() < 1)
		return;
	for(unsigned int c(0); c < p[pIndex].childIndex.size() - 1; c++){
		for(unsigned int cc(c + 1); cc < p[pIndex].childIndex.size(); cc++){
			double lambdaL(p[p[pIndex].childIndex[c]].length/p[p[pIndex].childIndex[cc]].length);
			if(lambdaL > 1.0)
				lambdaL = 1.0/lambdaL;
			lambdaLs.push_back(lambdaL);
			addBranchingLengthRatios(lambdaLs, p, p[pIndex].childIndex[c]);
		}
	}
	addBranchingLengthRatios(lambdaLs, p, p[pIndex].childIndex.back());
}

// assumes p[0] starts at the heart
vector<double> branchingLengthRatios(const vector<PathOfBranchPoint2D> &p){
	vector<double> lambdaLs;
	addBranchingLengthRatios(lambdaLs, p);
	return lambdaLs;
}

void addParentChildLengthRatios(vector<double> &gammas, const vector<PathOfBranchPoint2D> &p, int pIndex = 0){
	if(p[pIndex].childIndex.size() < 1)
		return;
	for(unsigned int c(0); c < p[pIndex].childIndex.size() - 1; c++){
		gammas.push_back(p[p[pIndex].childIndex[c]].length/p[pIndex].length);
		addParentChildLengthRatios(gammas, p, p[pIndex].childIndex[c]);
	}
}

vector<double> parentChildLengthRatios(const vector<PathOfBranchPoint2D> &p){
	vector<double> gammas;
	addParentChildLengthRatios(gammas, p);
	return gammas;
}

void testTrimDegenerateSinks(){
	kisset(4297, 8618, 4779);  // numOpt(10), graftPasses(200), settleIts(10000)
	int maxFailedInsertions(10000), numOpt(1), graftPasses(20), settleIts(1000), numBins(10);
	double nearestRadius(1.0), localRad(2.0*nearestRadius),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(9.0),
		maxPathLengthCoefficient(5.0),
		settleThresh(0.001), consolidateThresh(0.01);
	double sizeX(15.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	string uniquer("testTrimDegenerateSinks_" + makeString(sizeX));
	int width(1000), height(1000), border(20);

	cout << "\n testTrimDegenerateSinks()\n uniquer = " << uniquer
		<< "\n sizeX = " << sizeX //<< "\n sizeY = " << sizeY
		<< "\n heartx = " << heartx //<< "\n hearty = " << hearty
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient" << maxPathLengthCoefficient
		<< "\n settleThresh = " << settleThresh
		<< "\n consolidateThresh = " << consolidateThresh
		<< "\n nearestRadius = " << nearestRadius
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n numOpt = " << numOpt
		<< "\n graftPasses = " << graftPasses << "\n settleIts = " << settleIts
		<< "\n localRad = " << localRad << "\n numBins = " << numBins
		<< endl;
	vector<Bifurcation2D> b(initializeEvenRandomSpacingCircle(maxFailedInsertions, sizeX, nearestRadius));
	vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingle(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, graftPasses, settleIts);
	drawTreeCirc("outputs/" + uniquer + "_opt.png", bGraft, width, height, border);
	vector<PathOfBranchPoint2D> bPath(trimDegenerateSinks(bGraft, consolidateThresh));
	printPathOfBranchPoint2DVector(bPath);
	vector<double> lambdaLs(branchingLengthRatios(bPath));
	cout << "\n lambdaLs: " << makeVectorString(lambdaLs) << endl;
	map<double, double> lambdaLDist(probabilityDistribution(lambdaLs, numBins));
	printMap(lambdaLDist, "lambda_L");
}

void ensembleTrimDegenerateSinks(){
	kisset(4297, 8618, 4779);  // numOpt(10), graftPasses(200), settleIts(10000)
	int maxFailedInsertions(10000), numOpt(10), graftPasses(200), settleIts(10000), numBins(27), numNetworks(100);
	double nearestRadius(1.0), localRad(2.0*nearestRadius),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(9.0),
		maxPathLengthCoefficient(5.0),
		settleThresh(0.001), consolidateThresh(0.01);
	double sizeX(8.0), sizeY(sizeX);
	double heartx(sizeX/2.0), hearty(nearestRadius/2.0);
	string uniquer("ensembleTrimDegenerateSinks" + makeString(sizeX));
	int width(1000), height(1000), border(20);

	cout << "\n ensembleTrimDegenerateSinks()\n uniquer = " << uniquer
		<< "\n sizeX = " << sizeX //<< "\n sizeY = " << sizeY
		<< "\n heartx = " << heartx //<< "\n hearty = " << hearty
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient" << maxPathLengthCoefficient
		<< "\n settleThresh = " << settleThresh
		<< "\n consolidateThresh = " << consolidateThresh
		<< "\n nearestRadius = " << nearestRadius
		<< "\n maxFailedInsertions = " << maxFailedInsertions << "\n numOpt = " << numOpt
		<< "\n graftPasses = " << graftPasses << "\n settleIts = " << settleIts
		<< "\n localRad = " << localRad << "\n numBins = " << numBins
		<< "\n numNetworks = " << numNetworks
		<< endl;
	vector<double> allLambdaLs, allGammas;
	for(int n(0); n < numNetworks; n++){
		time_t tempStart(time(NULL));
		vector<Bifurcation2D> b(initializeEvenRandomSpacingCircle(maxFailedInsertions, sizeX, nearestRadius));
		cout << endl << b.size() - 1 << " TSVs" << endl;
		vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingle(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, graftPasses, settleIts);
		cout << " Optimization took " << niceTime(difftime(time(NULL), tempStart)) << endl;
		drawTreeCirc(uniquer + "_" + makeString(n) + "_opt.png", bGraft, width, height, border);
		vector<PathOfBranchPoint2D> bPath(trimDegenerateSinks(bGraft, consolidateThresh));
		//printPathOfBranchPoint2DVector(bPath);
		vector<double> lambdaLs(branchingLengthRatios(bPath));
		allLambdaLs.insert(allLambdaLs.end(), lambdaLs.begin(), lambdaLs.end());
		vector<double> gammas(parentChildLengthRatios(bPath));
		allGammas.insert(allGammas.end(), gammas.begin(), gammas.end());
		//cout << "\n lambdaLs: " << makeVectorString(lambdaLs) << endl;
	}
	map<double, double> lambdaLDist(probabilityDistribution(allLambdaLs, numBins));
	printMap(lambdaLDist, "lambda_L");
	map<double, double> gammaDist(probabilityDistribution(allGammas, numBins));
	printMap(gammaDist, "gamma");
}

/*
// b has no connections
vector<BranchPoint2D> optimizeByGlobalConsolidatedGraftingSingleAndRadius(int numOpt, vector<Bifurcation2D> b,
	double totalLengthCoefficient, double avePathLengthCoefficient, double maxPathLengthCoefficient,
	double settleThresh, double consolidateThresh, int graftPasses, int settleIts){
}

void optRadiusTuning(){
	int maxFailedInsertions(10000), numOpt(10), graftPasses(200), settleIts(10000);
	double sizeX(5.0), nearestRadius(1.0), settleThresh(0.01), consolidateThresh(0.1), 
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0);
	time_t tempStart(time(NULL));
	vector<Bifurcation2D> b(initializeEvenRandomSpacingCircle(maxFailedInsertions, sizeX, nearestRadius));
	vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGraftingSingleAndRadius(numOpt, b, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, graftPasses, settleIts);

}*/

double settleBifurcationsWithRadiusDownFrom(int i, vector<Bifurcation2D> &b, terseSeg **bt, vector<BranchPoint2D> &bc, double consolidateThresh, int *terseIndex, double rootRad, double betaRad, double capRad, double radThresh){
	double settleChange(0.0);
	for(unsigned int c(0); c < b[i].childIndex.size(); c++){
		int cIndex(b[i].childIndex[c]);
		if(b[cIndex].parentIndex >= 0 && b[cIndex].childIndex.size() > 0) // child not the heart, not a TSV
			settleChange += settleBifurcationWithRadius(cIndex, b, bt, bc, consolidateThresh, terseIndex, rootRad, betaRad, capRad, radThresh);
		//cout << "\n settleBifurcationsWithRadiusDownFrom() i = " << i << " c = " << c << " settleChange 0  = " << settleChange << endl;
		settleChange += settleBifurcationsWithRadiusDownFrom(cIndex, b, bt, bc, consolidateThresh, terseIndex, rootRad, betaRad, capRad, radThresh);
		//cout << "\n settleBifurcationsWithRadiusDownFrom() i = " << i << " c = " << c << " settleChange 1  = " << settleChange << endl;
	}
	return settleChange;
}

double settleBifurcationsWithRadiusUpFrom(int i, vector<Bifurcation2D> &b, terseSeg **bt, vector<BranchPoint2D> &bc, double consolidateThresh, int *terseIndex, double rootRad, double betaRad, double capRad, double radThresh){
	double settleChange(0.0);
	int pIndex(b[i].parentIndex);
	if(pIndex < 0)
		return 0.0;
	if(b[pIndex].parentIndex >= 0 && b[pIndex].childIndex.size() > 0){// parent not the heart, not a TSV
		settleChange = settleBifurcationWithRadius(pIndex, b, bt, bc, consolidateThresh, terseIndex, rootRad, betaRad, capRad, radThresh);
		cout << "\n settleBifurcationsWithRadiusUpFrom() i = " << i << " settleChange 0  = " << settleChange << endl;
		settleChange += settleBifurcationsWithRadiusUpFrom(pIndex, b, bt, bc, consolidateThresh, terseIndex, rootRad, betaRad, capRad, radThresh);
		cout << "\n settleBifurcationsWithRadiusUpFrom() i = " << i << " settleChange 1  = " << settleChange << endl;
		if(b[pIndex].childIndex.size() == 2){
			int cIndex(otherChildIndex(b, pIndex, i));
			settleChange += settleBifurcationWithRadius(cIndex, b, bt, bc, consolidateThresh, terseIndex, rootRad, betaRad, capRad, radThresh);
			cout << "\n settleBifurcationsWithRadiusUpFrom() i = " << i << " settleChange 2  = " << settleChange << endl;
			settleChange += settleBifurcationsWithRadiusDownFrom(cIndex, b, bt, bc, consolidateThresh, terseIndex, rootRad, betaRad, capRad, radThresh);
			cout << "\n settleBifurcationsWithRadiusUpFrom() i = " << i << " settleChange 3  = " << settleChange << endl;
		}
	}
	return settleChange;
}

double settleBifurcationsWithRadiusStartAt(vector<Bifurcation2D> &b, int startAt, double settleThresh, double consolidateThresh, double capRad, int settleIts){
	//cout << "\n settleBifurcationsWithRadiusStartAt(startAt = " << startAt << ")" << endl;
	double rootRad(0.01), betaRad(0.95), radThresh(capRad/2.0), settleChange(2.0*settleThresh);
	vector<BranchPoint2D> bc(consolidateDegenerates(b, consolidateThresh));
	terseSeg** bt = createTerseSegTree(b.size(), bc, rootRad, betaRad, capRad, radThresh);
	int *terseIndex = new int[b.size()];
	for(unsigned int i(0); i < b.size(); i++)
		terseIndex[i] = findCorrespondingTerseIndex(b, i, bc, consolidateThresh);
	//cout << "\n\t settleBifurcationsWithRadiusStartAt():\n\tb:\n" << makeStringTree(b) << endl;
	//cout << "\t bc:\n" << makeStringTree(bc) << endl;
	//cout << "\t terseIndex:\n" << makeArrayString(terseIndex, b.size()) << endl;
	//cout << "\t bt tree:\n" << makeStringTreeRadOnly(getRoot(bc, bt)) << endl;
	//cout << endl;
	for(int it(0); it < settleIts && settleChange > settleThresh; it++){
		settleChange = 0.0;
		if(b[startAt].parentIndex >= 0 && b[startAt].childIndex.size() > 0) // not the heart, not a TSV
			settleChange = settleBifurcationWithRadius(startAt, b, bt, bc, consolidateThresh, terseIndex, rootRad, betaRad, capRad, radThresh);
		double changeDown(settleBifurcationsWithRadiusDownFrom(startAt, b, bt, bc, consolidateThresh, terseIndex, rootRad, betaRad, capRad, radThresh));
		double changeUp(settleBifurcationsWithRadiusUpFrom(startAt, b, bt, bc, consolidateThresh, terseIndex, rootRad, betaRad, capRad, radThresh));
		settleChange += changeUp + changeDown;

		//cout << "\n b after Fermat during it = " << it << " gives:\n" << makeStringTree(b) << endl;

		// move gobs
		vector<vector<int> > gobs, boundaries;
		findDegenerateGobs(b, boundaries, gobs, consolidateThresh);
		for(unsigned int g(0); g < gobs.size(); g++){
			if(gobs[g].size() > 1){ // otherwise leave it to Fermat
				int closestBoundaryToHeart(boundaries[g][0]);
				unsigned int *lvlFromHeart(levelsFromHeart(b));
				for(unsigned int j(1); j < boundaries[g].size(); j++){
					if(lvlFromHeart[closestBoundaryToHeart] > lvlFromHeart[boundaries[g][j]])
						closestBoundaryToHeart = boundaries[g][j];
				}
				int closestGobToHeart(gobs[g][0]);
				for(unsigned int j(1); j < gobs[g].size(); j++){
					if(lvlFromHeart[closestGobToHeart] > lvlFromHeart[gobs[g][j]])
						closestGobToHeart = gobs[g][j];
				}
				delete[] lvlFromHeart;
				vector<double> w;
				for(unsigned int j(0); j < boundaries[g].size(); j++){
					if(boundaries[g][j] == closestBoundaryToHeart){
						if(terseIndex[closestGobToHeart] > 0)
							w.push_back(bt[terseIndex[closestGobToHeart]]->rad);
						else
							w.push_back(rootRad);
					}else{
						if(b[boundaries[g][j]].childIndex.size() > 0){
							if(terseIndex[boundaries[g][j]] < 1)
								w.push_back(rootRad);
							else
								w.push_back(bt[terseIndex[boundaries[g][j]]]->rad);
						}else
							w.push_back(capRad);
					}
				}
				// with gradual underrelaxation
				double prevX(b[gobs[g][0]].position[0]), prevY(b[gobs[g][0]].position[1]),
					underrelaxationFactor(0.5*exp(-it/double(b.size())) + 0.5);
				TorresWeiszfeld(b, boundaries[g], w, gobs[g], settleThresh);
				double dx(b[gobs[g][0]].position[0] - prevX), dy(b[gobs[g][0]].position[1] - prevY);
				double underX(prevX + underrelaxationFactor*dx), underY(prevY + underrelaxationFactor*dy);
				for(unsigned int j(0); j < gobs[g].size(); j++){
					b[gobs[g][j]].position[0] = underX;
					b[gobs[g][j]].position[1] = underY;
				}
				bc = consolidateDegenerates(b);
				recreateTerseSegTree(bt, bc, rootRad, betaRad, capRad, radThresh);
				for(unsigned int j(0); j < b.size(); j++)
					terseIndex[j] = findCorrespondingTerseIndex(b, j, bc, consolidateThresh);
				//cout << "\n updated bt to:" << endl;
				//for(unsigned int j(1); j < bc.size(); j++)
				//	cout << j << "\t" << bt[j]->rad << endl;
			}
		}

		//cout << "\n settleBifurcationsWithRadiusStartAt(): it = " << it << " \t settleChange = " << settleChange
		//	<< "\n\t\t\t\t changeDown = " << changeDown << "\n\t\t\t\t changeUp = " << changeUp
		//	<< endl;
		//drawTree("outputs/settleBifurcationsWithRadiusStartAt_" + makeString(it) + ".png", bc, 1000, 1000, 25);
	}
	return settleChange;
}

void testWeightedFermatSettle(){
	kisset (14804, 13905, 7418);
	string uniquer("outputs/testWeightedFermat");
	int maxFailedInsertions(10000), settleIts(10000), width(500), height(width), border(width/40);
	double sizeX(4.0), nearestRadius(1.0), settleThresh(0.01), consolidateThresh(0.1), capRad(0.005);
	time_t tempStart(time(NULL));
	vector<Bifurcation2D> b(initializeEvenRandomSpacingCircle(maxFailedInsertions, sizeX, nearestRadius));
	/*b.push_back(Bifurcation2D());
	b.back().position[0] = 0.0;
	b.back().position[1] = 0.0;
	b.push_back(Bifurcation2D());
	b.back().position[0] = 1.0;
	b.back().position[1] = 0.0;
	b.push_back(Bifurcation2D());
	b.back().position[0] = 0.0;
	b.back().position[1] = 1.0;*/
	drawBifurcationTree(uniquer + "_TSVs.png", b, width, height, border);
	introduceBifurcationsFromSources(b);
	drawTreeCirc(uniquer + "_NearestNeighbor.png", consolidateDegenerates(b, consolidateThresh), width, height, border);
	for(unsigned int i(0); i < 0 + 1*b.size(); i++){
		vector<Bifurcation2D> bCopy(b);
		double settleChange(settleBifurcationsWithRadiusStartAt(bCopy, i, settleThresh, consolidateThresh, capRad, settleIts));
		cout << "\n\n Settled with radius starting at " << i << " at time " << niceTime(difftime(time(NULL), tempStart))
			<< " with settleChange = " << settleChange << endl;
		if(settleChange > settleThresh)
			cout << "\t WARNING! settleChange = " << settleChange << " > " << settleThresh << " = " << settleThresh << endl;
		drawTreeCirc(uniquer + "_settled_startedAt" + makeString(i) + ".png", consolidateDegenerates(bCopy), width, height, border);
		cout << "\t\t Drew tree from starting at " << i << " at time " << niceTime(difftime(time(NULL), tempStart)) << endl;
		cout << "\t Resulting tree" << endl;
		for(unsigned int j(0); j < b.size(); j++)
			cout << j << "\t" << bCopy[j].position[0] << "\t" << bCopy[j].position[1] << endl;
	}
}

int rpp(int &x){return x++;}
int ppr(int &x){return ++x;}

void stupidTest(){
	int x(0);
	cout << "\n x = " << x << endl;
	int y(rpp(x));
	cout << "After rpp, x = " << x << " and y = " << y << endl;
	y = ppr(x);
	cout << "After ppr, x = " << x << " and y = " << y << endl;
	int *z = new int;
	*z = 5;
	cout << "\n *z = " << *z << endl;
	delete z;
	z = new int;
	*z = 4;
	cout << "\n *z = " << *z << endl;
}

void checkRegraftForSymmetry(vector<Bifurcation2D> &b, int *nTips, double forceSymThresh, double settleThresh, int settleIts, int i = 0){
	for(unsigned int c(0); c < b[i].childIndex.size(); c++)
		checkRegraftForSymmetry(b, nTips, forceSymThresh, settleThresh, settleIts, b[i].childIndex[c]);
	if(b[i].childIndex.size() == 2){
		double tipRatio(1.0);
		if(nTips[b[i].childIndex[0]] < nTips[b[i].childIndex[1]])
			tipRatio = double(nTips[b[i].childIndex[0]] + 1)/double(nTips[b[i].childIndex[1]]);
		else if(nTips[b[i].childIndex[0]] > nTips[b[i].childIndex[1]])
			tipRatio = double(nTips[b[i].childIndex[1]] + 1)/double(nTips[b[i].childIndex[0]]);
		if(tipRatio < forceSymThresh){
			int statBi(b[i].childIndex[0]), moBi(b[i].childIndex[1]);
			if(nTips[b[i].childIndex[0]] < nTips[b[i].childIndex[1]]){
				statBi = moBi;
				moBi = b[i].childIndex[0];
			}
			int changedIndex(0);
			if(b[statBi].childIndex.size() == 2){// regraft moBi to closer child of statBi
				//drawTree("outputs/moBi" + makeString(moBi) + ".png", consolidateDegenerates(b), 500, 500, 500/40);
				if(separation2D(b[moBi], b[b[statBi].childIndex[0]]) < separation2D(b[moBi], b[b[statBi].childIndex[1]])){ // child 0 is new parent of moBi
					//cout << "\nMoving " << moBi << " to be child of " << statBi << endl;
					if(!regraftSingle(b, statBi, moBi, true)){
						cout << "\n Error checkRegraftForSymmetry: cannot move b[moBi = " << moBi << "] to parent " << b[statBi].childIndex[0] << endl;
						cout << "\n\t (" << b[moBi].position[0] << ", " << b[moBi].position[1] << ") to parent ("
							<< b[b[statBi].childIndex[0]].position[0] << ", " << b[b[statBi].childIndex[0]].position[1] << ")" << endl;
						abort();
					}
				}else{ // child 1 is new parent of moBi
					changedIndex = 1;
					cout << "\nMoving " << moBi << " to be child of " << statBi << endl;
					if(!regraftSingle(b, statBi, moBi, false)){//b, moBi, b[statBi].childIndex[1], true
						cout << "\n Error checkRegraftForSymmetry: cannot move b[moBi = " << moBi << "] to parent " << b[statBi].childIndex[1] << endl;
						cout << "\n\t (" << b[moBi].position[0] << ", " << b[moBi].position[1] << ") to parent ("
							<< b[b[statBi].childIndex[1]].position[0] << ", " << b[b[statBi].childIndex[1]].position[1] << ")" << endl;
						abort();
					}
				}
				allNumTips(b, i, nTips);
				settleBifurcations(b, false, settleThresh, settleIts);
				//drawTree("outputs/statBi" + makeString(statBi) + "p" + makeString(moBi) + ".png", consolidateDegenerates(b), 500, 500, 500/40);
			}else{
				cout << "\n Warning checkRegraftForSymmetry: skipping b[statBi = " << statBi
					<< "] because it has childIndex.size() = " << b[statBi].childIndex.size() << endl;
			}
		}
	}
}

// uses adjustedMinimumTipRatio to determine if an update is needed
vector<Bifurcation2D> regraftForSymmetry(const vector<Bifurcation2D> &bOrig, double forceSymThresh, double settleThresh, int settleIts){
	vector<Bifurcation2D> b(bOrig);
	int *nTips(numTips(b)), maxRegraftIts(b.size());
	for(int it(0); it < maxRegraftIts && adjustedMinimumTipRatio(b, nTips) < forceSymThresh; it++){
		checkRegraftForSymmetry(b, nTips, forceSymThresh, settleThresh, settleIts);
	}
	if(adjustedMinimumTipRatio(b, nTips) < forceSymThresh)
		cout << "\n Warning regraftForSymmetry(): aMTR = " << adjustedMinimumTipRatio(b, nTips) << " after " << maxRegraftIts << " passes." << endl;
	return b;
}

void testRegraftForSymmetry(){
	kisset (7446, 13195, 32375);
	string uniquer("outputs/testRegraftForSymmetry");// "OneFour" for oneFour, "Square" for square
	bool doDraw(true);
	double forceSymThresh(0.1), side(3.0),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),heartx(side/2.0), hearty(nearestRadius/2.0),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0),
		localRad(2.0*nearestRadius),
		tipRatioThresh(0.01);
	int numOpt(1), graftPasses(100), settleIts(10000),
		maxFailedInsertions(10000), width(1000), height(1000), border(width/40);// width 250 for oneFour
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		<< "\n forceSymThresh = " << forceSymThresh
		<< "\n side = " << side << "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n localRad = " << localRad << "\n tipRatioThresh = " << tipRatioThresh
		<< "\n graftPasses = " << graftPasses
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;

	vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, 1.0*side, 1.0*side, nearestRadius, heartx, hearty));
	
	vector<Bifurcation2D> bNSC(b);
	introduceBifurcationsFromSources(bNSC);
	vector<BranchPoint2D> bNSCGuess = consolidateDegenerates(bNSC, consolidateThresh);
	cout << "\n bNSC mTR = " << minimumTipRatio(bNSC) << " \t aMTR = " << adjustedMinimumTipRatio(bNSC)
		<< "\n bNSCGuess mTR = " << minimumTipRatio(bNSCGuess)
		<< endl;
	double NSCGuessMeasure = tripletMeasure(bNSCGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	cout << "\n\t\t\t\t\t\t bNSCGuess measure = " << NSCGuessMeasure << endl;
	if(doDraw)
		drawTree(uniquer + makeString((int)ceil(side)) + "_NSCguess.png", bNSCGuess, width, height, border);

	vector<Bifurcation2D> bCopy(b);
	introduceBifurcationsSymmetricallyFromHeart(bCopy, settleThresh);
	vector<BranchPoint2D> bGuess = consolidateDegenerates(bCopy, consolidateThresh);
	cout << "\n bCopy mTR = " << minimumTipRatio(bCopy) << " \t aMTR = " << adjustedMinimumTipRatio(bCopy)
		<< "\n bGuess mTR = " << minimumTipRatio(bGuess)
		<< endl;
	double guessMeasure = tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	cout << "\n\t\t\t\t\t\t Guess measure = " << guessMeasure << endl;
	if(doDraw)
		drawTree(uniquer + makeString((int)ceil(side)) + "_guess.png", bGuess, width, height, border);



	//vector<Bifurcation2D> bSymBiNSC(regraftForSymmetry(bNSC, forceSymThresh, settleThresh, settleIts));
	//vector<BranchPoint2D> bSymNSC(consolidateDegenerates(bSymBiNSC));
	//cout << "\n bSymBiNSC has minimumTipRatio = " << minimumTipRatio(bSymBiNSC) 
	//	<< "\n bSymNSC has minimumTipRatio = " << minimumTipRatio(bSymNSC) << endl;
	//double symNSCMeasure = tripletMeasure(bSymNSC, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	//cout << "\n\t\t\t\t\t\t bSymNSC measure = " << symNSCMeasure << endl;
	//if(doDraw)
	//	drawTree(uniquer + makeString((int)ceil(side)) + "_NSCsym.png", bSymNSC, width, height, border);
	
	//vector<Bifurcation2D> bSymBi(regraftForSymmetry(bCopy, forceSymThresh, settleThresh, settleIts));
	//vector<BranchPoint2D> bSym(consolidateDegenerates(bSymBi));
	//cout << "\n bSymBi has minimumTipRatio = " << minimumTipRatio(bSymBi) 
	//	<< "\n bSym has minimumTipRatio = " << minimumTipRatio(bSym) << endl;
	//double symMeasure = tripletMeasure(bSym, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	//cout << "\n\t\t\t\t\t\t bSym measure = " << symMeasure << endl;
	//if(doDraw)
	//	drawTree(uniquer + makeString((int)ceil(side)) + "_sym.png", bSym, width, height, border);

	time_t tempStart(time(NULL));
	vector<Bifurcation2D> bBi;
	vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingleAndTipRatioThresh(numOpt, b, bBi, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, tipRatioThresh, graftPasses, settleIts);
	double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	cout << "\n\t\t Graft took " << niceTime(difftime(time(NULL), tempStart)) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = "
		<< graftMeasure << "\n\t\t minimumTipRatio(bBi) = " << minimumTipRatio(bBi)
		<< "\n\t\t minimumTipRatio(bGraft) = " << minimumTipRatio(bGraft) << endl;
	if(doDraw)
		drawTree(uniquer + makeString((int)ceil(side)) + "_graft.png", bGraft, width, height, border);
}

void regraftForSymmetryEnsemble(){
	kisset (7446, 13195, 32375);
	string uniquer("outputs/regraftForSymmetryEnsemble0p0_");// "OneFour" for oneFour, "Square" for square
	bool doDraw(true);
	double side(10.0), //forceSymThresh(0.2),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),heartx(side/2.0), hearty(nearestRadius/2.0),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(9.0),
		maxPathLengthCoefficient(0.0),
		localRad(2.0*nearestRadius),
		tipRatioThresh(0.0);
	int numOpt(10), graftPasses(100), settleIts(10000), numBins(20), numNetworks(20),
		maxFailedInsertions(10000), width(1000), height(1000), border(width/40);// width 250 for oneFour
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		//<< "\n forceSymThresh = " << forceSymThresh
		<< "\n side = " << side << "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n localRad = " << localRad << "\n tipRatioThresh = " << tipRatioThresh
		<< "\n graftPasses = " << graftPasses << "\n numBins = " << numBins
		<< "\n numNetworks = " << numNetworks
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<double> lambdaLs;
	for(int n(0); n < numNetworks; n++){
		vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, 1.0*side, 1.0*side, nearestRadius, heartx, hearty));
		vector<Bifurcation2D> bCopy(b);
		introduceBifurcationsSymmetricallyFromHeart(bCopy, settleThresh);
		vector<BranchPoint2D> bGuess = consolidateDegenerates(bCopy, consolidateThresh);
		//cout << "\n bCopy mTR = " << minimumTipRatio(bCopy) << " \t aMTR = " << adjustedMinimumTipRatio(bCopy)
		//	<< "\n bGuess mTR = " << minimumTipRatio(bGuess)
		//	<< endl;
		//double guessMeasure = tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
		//cout << "\n\t\t\t\t\t\t Guess measure = " << guessMeasure << endl;
		if(doDraw)
			drawTree(uniquer + makeString((int)ceil(side)) + "_" + makeString(n) + "_guess.png", bGuess, width, height, border);
		time_t tempStart(time(NULL));
		vector<Bifurcation2D> bBi;
		vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingleAndTipRatioThresh(numOpt, b, bBi, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, tipRatioThresh, graftPasses, settleIts);
		double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
		cout << "\n\t\t Graft took " << niceTime(difftime(time(NULL), tempStart)) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = "
			<< graftMeasure << "\n\t\t minimumTipRatio(bBi) = " << minimumTipRatio(bBi) << " \t aMTR(bBi) = " << adjustedMinimumTipRatio(bBi)
			<< "\n\t\t minimumTipRatio(bGraft) = " << minimumTipRatio(bGraft) << endl;
		if(doDraw)
			drawTree(uniquer + makeString((int)ceil(side)) + "_" + makeString(n) + "_graft.png", bGraft, width, height, border);
		int ignoreCount(0);
		vector<double> lambdaL(validBranchingLengthRatios(bGraft, ignoreCount, consolidateThresh));
		lambdaLs.insert(lambdaLs.begin(), lambdaL.begin(), lambdaL.end());
		cout << "\n" << n << " ignored " << ignoreCount << endl;
	}
	map<double, double> lambdaLDist(probabilityDistribution(lambdaLs, numBins));
	printMap(lambdaLDist, "lambda_L");
}

void assignHortonStrahlerLevel(const vector<BranchPoint2D> &b, unsigned int *HS, int i = 0){
	int highestLevel(1), highestLevelCount(0);
	for(unsigned int c(0); c < b[i].childIndex.size(); c++){
		assignHortonStrahlerLevel(b, HS, b[i].childIndex[c]);
		if(highestLevel == HS[b[i].childIndex[c]])
			highestLevelCount++;
		else if(highestLevel < (int)HS[b[i].childIndex[c]]){
			highestLevel = HS[b[i].childIndex[c]];
			highestLevelCount = 1;
		}
	}
	HS[i] = highestLevel;
	if(highestLevelCount > 1)
		HS[i]++;
}

unsigned int* HortonStrahlerLevels(const vector<BranchPoint2D> &b){
	unsigned int *HS = new unsigned int[b.size()];
	assignHortonStrahlerLevel(b, HS);
	return HS;
}

void levelsToHeart(unsigned int *LH, int level, int index, const vector<BranchPoint2D> &b){
	LH[index] = level;
	for(unsigned int c(0); c < b[index].childIndex.size(); c++)
		levelsToHeart(LH, level + 1, b[index].childIndex[c], b);
}


unsigned int* levelsToHeart(const vector<BranchPoint2D> &b){
	unsigned int *LH = new unsigned int[b.size()];
	levelsToHeart(LH, 0, 0, b);
	return LH;
}

void coutFullGamma(const vector<BranchPoint2D> &b){
	cout << "tree:\n" << recordString(b) << endl;
	// do more . . .
}

void testGammaAnalysis(){
	kisset (7446, 13195, 32375);
	string uniquer("outputs/testRegraftForSymmetry");// "OneFour" for oneFour, "Square" for square
	bool doDraw(true);
	double side(3.0),//forceSymThresh(0.1), 
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))), heartx(side/2.0), hearty(nearestRadius/2.0),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0);
		//localRad(2.0*nearestRadius),
		//tipRatioThresh(0.01);
	int numOpt(1), settleIts(10000), //graftPasses(100),
		maxFailedInsertions(10000), width(1000), height(1000), border(width/40);// width 250 for oneFour
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		//<< "\n forceSymThresh = " << forceSymThresh
		<< "\n side = " << side << "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		//<< "\n localRad = " << localRad << "\n tipRatioThresh = " << tipRatioThresh
		//<< "\n graftPasses = " << graftPasses
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;

	vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, 1.0*side, 1.0*side, nearestRadius, heartx, hearty));
	
	vector<Bifurcation2D> bNSC(b);
	introduceBifurcationsFromSources(bNSC);
	vector<BranchPoint2D> bNSCGuess = consolidateDegenerates(bNSC, consolidateThresh);
	cout << "\n bNSC mTR = " << minimumTipRatio(bNSC) << " \t aMTR = " << adjustedMinimumTipRatio(bNSC)
		<< "\n bNSCGuess mTR = " << minimumTipRatio(bNSCGuess)
		<< endl;
	double NSCGuessMeasure = tripletMeasure(bNSCGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
	cout << "\n\t\t\t\t\t\t bNSCGuess measure = " << NSCGuessMeasure << endl;
	if(doDraw)
		drawTree(uniquer + makeString((int)ceil(side)) + "_NSCguess.png", bNSCGuess, width, height, border);
	coutFullGamma(bNSCGuess);
}


void testRecordReadNetwork(){
	string uniquer("outputs/recordReadTest");
	double side(3.0), nearestRadius(1.0/sqrt(acos(-1.0))), heartx(side/2.0), hearty(nearestRadius/2.0);
	int maxFailedInsertions(10000), width(1000), height(1000), border(width/40);
	vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, 1.0*side, 1.0*side, nearestRadius, heartx, hearty));
	vector<Bifurcation2D> bNSC(b);
	introduceBifurcationsFromSources(bNSC);
	vector<BranchPoint2D> bcNSC(consolidateDegenerates(bNSC));
	drawTree(uniquer + makeString((int)ceil(side)) + "_bcNSC.png", bcNSC, width, height, border);
	string recTree(recordString(bNSC));
	cout << "\n\n recTree:\n" << recTree << endl;
	vector<Bifurcation2D> bTrans(readRecordedString(recTree));
	cout << "\n bTrans recorded =\n" << recordString(bTrans) << endl;
	cout << "\n oneLine before = " << makeStringTreeOneLine(bNSC) << endl;
	cout << "\n oneLine after = " << makeStringTreeOneLine(bTrans) << endl;
	vector<BranchPoint2D> bcTrans(consolidateDegenerates(bTrans));
	drawTree(uniquer + makeString((int)ceil(side)) + "_bcTrans.png", bcTrans, width, height, border);
}

void plainGammaDistribution(){
	kisset (7446, 13195, 32375);
	string uniquer("dists0p5_");// "OneFour" for oneFour, "Square" for square
	bool doDraw(true);
	double side(5.0), //forceSymThresh(0.2),
		settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))),heartx(side/2.0), hearty(nearestRadius/2.0),
		totalLengthCoefficient(1.0),
		avePathLengthCoefficient(0.0),
		maxPathLengthCoefficient(0.0),
		localRad(2.0*nearestRadius),
		tipRatioThresh(0.5);
	int numOpt(10), graftPasses(100), settleIts(10000), numBins(20), numNetworks(20),
		maxFailedInsertions(10000), width(1000), height(1000), border(width/40);// width 250 for oneFour
	cout << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		//<< "\n forceSymThresh = " << forceSymThresh
		<< "\n side = " << side << "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n localRad = " << localRad << "\n tipRatioThresh = " << tipRatioThresh
		<< "\n graftPasses = " << graftPasses << "\n numBins = " << numBins
		<< "\n numNetworks = " << numNetworks
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	vector<double> lambdaLs, gammas, allGammas;
	for(int n(0); n < numNetworks; n++){
		vector<Bifurcation2D> b(initializeEvenRandomSpacing(maxFailedInsertions, 1.0*side, 1.0*side, nearestRadius, heartx, hearty));
		vector<Bifurcation2D> bCopy(b);
		introduceBifurcationsSymmetricallyFromHeart(bCopy, settleThresh);
		vector<BranchPoint2D> bGuess = consolidateDegenerates(bCopy, consolidateThresh);
		//cout << "\n bCopy mTR = " << minimumTipRatio(bCopy) << " \t aMTR = " << adjustedMinimumTipRatio(bCopy)
		//	<< "\n bGuess mTR = " << minimumTipRatio(bGuess)
		//	<< endl;
		//double guessMeasure = tripletMeasure(bGuess, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
		//cout << "\n\t\t\t\t\t\t Guess measure = " << guessMeasure << endl;
		if(doDraw)
			drawTree(uniquer + makeString((int)ceil(side)) + "_" + makeString(n) + "_guess.png", bGuess, width, height, border);
		time_t tempStart(time(NULL));
		vector<Bifurcation2D> bBi;
		vector<BranchPoint2D> bGraft = optimizeByGlobalConsolidatedGraftingWithSingleAndTipRatioThresh(numOpt, b, bBi, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, tipRatioThresh, graftPasses, settleIts);
		double graftMeasure = tripletMeasure(bGraft, totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
		cout << "\n\t\t Graft took " << niceTime(difftime(time(NULL), tempStart)) << " for " << b.size() - 1 << " TSVs\t\t Graft measure = "
			<< graftMeasure << "\n\t\t minimumTipRatio(bBi) = " << minimumTipRatio(bBi) << " \t aMTR(bBi) = " << adjustedMinimumTipRatio(bBi)
			<< "\n\t\t minimumTipRatio(bGraft) = " << minimumTipRatio(bGraft) << endl;
		if(doDraw)
			drawTree(uniquer + makeString((int)ceil(side)) + "_" + makeString(n) + "_graft.png", bGraft, width, height, border);
		int ignoreCount(0);
		vector<double> lambdaL(validBranchingLengthRatios(bGraft, ignoreCount, consolidateThresh));
		vector<double> gamma(validParentChildLengthRatios(bGraft, ignoreCount, consolidateThresh)); // same ignore count
		vector<double> allGamma(allParentChildLengthRatios(bGraft));
		lambdaLs.insert(lambdaLs.begin(), lambdaL.begin(), lambdaL.end());
		gammas.insert(gammas.begin(), gamma.begin(), gamma.end());
		allGammas.insert(allGammas.begin(), allGamma.begin(), allGamma.end());
		cout << "\n" << n << " ignored " << ignoreCount << endl;
	}
	map<double, double> lambdaLDist(probabilityDistribution(lambdaLs, numBins)),
		gammaDist(probabilityDistribution(gammas, numBins)),
		allGammaDist(probabilityDistribution(allGammas, numBins));
	printMap(lambdaLDist, "lambda_L");
	printMap(gammaDist, "gamma");
	printMap(allGammaDist, "gamma(all)");
	ofstream outLambdaL(uniquer + "valid_lambda_L.dat");
	for(unsigned int i(0); i < lambdaLs.size(); i++)
		outLambdaL << lambdaLs[i] << endl;
	outLambdaL.close();
	ofstream outGamma(uniquer + "valid_gamma.dat");
	for(unsigned int i(0); i < gammas.size(); i++)
		outGamma << gammas[i] << endl;
	outGamma.close();
	ofstream outAllGamma(uniquer + "all_gamma.dat");
	for(unsigned int i(0); i < allGammas.size(); i++)
		outAllGamma << allGammas[i] << endl;
	outAllGamma.close();
}

void prebinnedDistributions(){
	vector<string> fns;
	string dir("outputs/gammas/");
	fns.push_back("allGamma_dicom.dat");
	fns.push_back("allGamma_mouse.dat");
	fns.push_back("allGamma_mouse825.dat");
	fns.push_back("c100/dists0p0_all_gamma.dat");
	fns.push_back("c100/dists0p0_valid_gamma.dat");
	fns.push_back("c100/dists0p1_all_gamma.dat");
	fns.push_back("c100/dists0p1_valid_gamma.dat");
	fns.push_back("c100/dists0p2_all_gamma.dat");
	fns.push_back("c100/dists0p2_valid_gamma.dat");
	fns.push_back("c100/dists0p3_all_gamma.dat");
	fns.push_back("c100/dists0p3_valid_gamma.dat");
	fns.push_back("c100/dists0p4_all_gamma.dat");
	fns.push_back("c100/dists0p4_valid_gamma.dat");
	fns.push_back("c100/dists0p5_all_gamma.dat");
	fns.push_back("c100/dists0p5_valid_gamma.dat");
	fns.push_back("c190/dists0p0_all_gamma.dat");
	fns.push_back("c190/dists0p0_valid_gamma.dat");
	fns.push_back("c190/dists0p1_all_gamma.dat");
	fns.push_back("c190/dists0p1_valid_gamma.dat");
	fns.push_back("c190/dists0p2_all_gamma.dat");
	fns.push_back("c190/dists0p2_valid_gamma.dat");
	fns.push_back("c190/dists0p3_all_gamma.dat");
	fns.push_back("c190/dists0p3_valid_gamma.dat");
	fns.push_back("c190/dists0p4_all_gamma.dat");
	fns.push_back("c190/dists0p4_valid_gamma.dat");
	fns.push_back("c190/dists0p5_all_gamma.dat");
	fns.push_back("c190/dists0p5_valid_gamma.dat");
	int numBins(100);
	double binEnd(5.0);
	double binWidth(binEnd/numBins);
	vector<unsigned int*> bins;
	vector<unsigned int> counts;
	for(unsigned int f(0); f < fns.size(); f++){
		ifstream inFile((dir + fns[f]).c_str());
		string ln;
		bins.push_back(new unsigned int[numBins + 1]);
		for(int i(0); i < numBins + 1; i++)
			bins.back()[i] = 0;
		counts.push_back(0);
		while(getline(inFile, ln)){
			double x(atof(ln.c_str()));
			counts.back()++;
			if(x < binWidth)
				bins.back()[0]++;
			else if(x > binEnd)
				bins.back()[numBins]++;
			else{
				bool found(false);
				for(int i(0); i < numBins && !found; i++){
					if(x < binWidth*(i + 1)){
						bins.back()[i]++;
						found = true;
					}
				}
				if(!found){
					cout << "\n Error prebinnedDistributions(): Did not find bin for x = " << x << " given binEnd = " << binEnd << "." << endl;
					counts.back()--;
				}
			}
		}
		inFile.close();
	}
	
	cout << "\nbinMiddle";
	for(unsigned int f(0); f < fns.size(); f++)
		cout << "\t" << fns[f];
	cout << endl;
	for(int i(0); i < numBins; i++){
		cout << (i + 0.5)*binWidth;
		for(unsigned int f(0); f < fns.size(); f++)
			cout << "\t" << (bins[f][i]/double(counts[f]))/binWidth;
		cout << endl;
	}
	cout << "\nfracOver";
	for(unsigned int f(0); f < fns.size(); f++)
		cout << "\t" << bins[f][numBins]/double(counts[f]);

	for(unsigned int f(0); f < fns.size(); f++)
		delete[] bins[f];
}

string paddedString(int x, int pad = 5){
	string s(makeString(x));
	s.insert(s.begin(), pad - s.size(), '0');
	return s;
}

void singleNetOptRec(double lengthCoef, double avePathCoef, double maxPathCoef, double tipRatioThresh, double side, int shape){
	string uniquer("saved_outputs//singleNetOptRec_c"+ makeString(lengthCoef, 2) + "-" + makeString(avePathCoef, 2) + "-" + makeString(maxPathCoef, 2)
		+ "_hierSym" + makeString(tipRatioThresh, 4) + "_side" + makeString(side));

	cout << uniquer << endl;
	if(shape == 0)
		uniquer += "circ";
	else if(shape == 1)
		uniquer += "sq";
	else if(shape == 2)
		uniquer += "elong";
	else{
		shape = 0;
		uniquer = "circ";
	}

	string enumeratorFn("enumerator.dat");
	ifstream inEnum(enumeratorFn.c_str());
	int count(0);
	if(inEnum.is_open())
		inEnum >> count;
	inEnum.close();
	uniquer += "_" + paddedString(count);
	ofstream outEnum(enumeratorFn.c_str());
	outEnum << count + 1;
	outEnum.close();

	bool doDraw(true);
	double settleThresh(0.001), consolidateThresh(0.01), lengthThresh(consolidateThresh),
		nearestRadius(1.0/sqrt(acos(-1.0))), heartx(side/2.0), hearty(nearestRadius/2.0),
		totalLengthCoefficient(lengthCoef),
		avePathLengthCoefficient(avePathCoef),
		maxPathLengthCoefficient(maxPathCoef),
		localRad(2.0*nearestRadius);
	if(shape == 2)
		 heartx = side/4.0;
	const int numOpt(10), graftPasses(100), settleIts(10000), numBins(20), numNetworks(20),
		maxFailedInsertions(10000);
	int width(1000), height(1000), border(width/40);// width 250 for oneFour
	if(shape == 2){
		width /= 2;
		height *= 2;
	}
	ofstream tempOut(uniquer + "_temp_out.txt");
	streambuf* old_buffer = cout.rdbuf(tempOut.rdbuf());

	kissprint();

	tempOut << "\n uniquer = " << uniquer << "\n overrelaxationFactor = " << overrelaxationFactor
		<< "\n side = " << side << "\n heartx = " << heartx << "\n hearty = " << hearty
		<< "\n settleThresh = " << settleThresh << "\nconsolidateThresh = " << consolidateThresh << "\n lengthThresh = " << lengthThresh
		<< "\n nearestRadius = " << nearestRadius << "\n numOpt = " << numOpt
		<< "\n totalLengthCoefficient = " << totalLengthCoefficient
		<< "\n avePathLengthCoefficient = " << avePathLengthCoefficient
		<< "\n maxPathLengthCoefficient = " << maxPathLengthCoefficient
		<< "\n localRad = " << localRad << "\n tipRatioThresh = " << tipRatioThresh
		<< "\n graftPasses = " << graftPasses << "\n numBins = " << numBins
		<< "\n numNetworks = " << numNetworks
		<< "\n settleIts = " << settleIts << "\n maxFailedInsertions = " << maxFailedInsertions
		<< "\n width = " << width << "\n height = " << height << "\n border = " << border
		<< endl;
	
	vector<Bifurcation2D> b;
	if(shape == 0)
		b = initializeEvenRandomSpacingCircle(maxFailedInsertions, side, nearestRadius); // <---
	else if(shape == 1)
		b = initializeEvenRandomSpacing(maxFailedInsertions, 1.0*side, 1.0*side, nearestRadius, heartx, hearty);
	else if(shape == 2)
		b = initializeEvenRandomSpacing(maxFailedInsertions, 0.5*side, 2.0*side, nearestRadius, heartx, hearty);
	vector<Bifurcation2D> bCopy(b);


	for (int i = 0; i < b.size(); i++)
	{
		cerr << b.at(i).position[0] << " " << b.at(i).position[1] << endl;
	}
	cerr << "Kappa" << endl;
	system("pause");
	/////////////////////////////////////////////////////////////////////////

	introduceBifurcationsSymmetricallyFromHeart(bCopy, settleThresh);

	vector<BranchPoint2D> bGuess= consolidateDegenerates(bCopy, consolidateThresh);
	drawTreeCirc(uniquer + "_guessSym.png", bGuess, width, height, border);


	if(doDraw){
		if(shape == 0)
			drawTreeCirc(uniquer + "_guessSym.png", bGuess, width, height, border);
		else
			drawTree(uniquer + "_guessSym.png", bGuess, width, height, border);
	}
	cout << "\n b has " << b.size() - 1 << " TSVs." << endl;
	//system("pause");
	vector<Bifurcation2D> bBi[numOpt];
	vector<BranchPoint2D> bGraft[numOpt];
	ofstream outAll(uniquer + "_allOpt.dat");

	for(int opt(0); opt < numOpt; opt++){
		time_t tempStart(time(NULL));
		bGraft[opt] = optimizeByGlobalConsolidatedGraftingWithSingleAndTipRatioThresh(1, b, bBi[opt], totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient, settleThresh, consolidateThresh, localRad, tipRatioThresh, graftPasses, settleIts);
		double graftMeasure = tripletMeasure(bGraft[opt], totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient);
		if(doDraw){
			if(shape == 0)
				drawTreeCirc(uniquer + "_optGraft" + paddedString(opt, 4) + ".png", bGraft[opt], width, height, border);
			else
				drawTree(uniquer + "_optGraft" + paddedString(opt, 4) + ".png", bGraft[opt], width, height, border);
		}
		outAll << recordString(bBi[opt]) << endl;
		cout << "\n\t\t Graft took " << niceTime(difftime(time(NULL), tempStart))
			<< " with opt = " << opt
			<< " for " << b.size() - 1 << " TSVs\t\t Graft measure = "
			<< graftMeasure << "\n\t\t minimumTipRatio(bBi) = " << minimumTipRatio(bBi[opt]) << " \t aMTR(bBi) = " << adjustedMinimumTipRatio(bBi[opt])
			<< "\n\t\t minimumTipRatio(bGraft) = " << minimumTipRatio(bGraft[opt]) << endl;
		
	}
	outAll.close();
	
	int optNet(0);
	for(int opt(1); opt < numOpt; opt++){
		if(tripletMeasure(bGraft[optNet], totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient) >
				tripletMeasure(bGraft[opt], totalLengthCoefficient, avePathLengthCoefficient, maxPathLengthCoefficient))
			optNet = opt;
	}
	cout << "\n optNet = " << optNet << endl;
	if(doDraw){
		if(shape == 0)
			drawTreeCirc(uniquer + "_optGraft.png", bGraft[optNet], width, height, border);
		else
			drawTree(uniquer + "_optGraft.png", bGraft[optNet], width, height, border);
	}
	ofstream outOpt(uniquer + "_opt.dat");
	outOpt << recordString(bBi[optNet]) << endl;
	
	outOpt.close();

	time_t endTime = time(NULL);
	cout << "\n\nRuntime: " << niceTime(difftime(endTime, startTime));
	cout << "\nProgram complete.\n";

	cout.rdbuf(old_buffer);
	tempOut.close();
}

int main0(int argc, char *argv[]){
	startTime = time(NULL);
	srand((int)startTime);
	kisset(rand(), rand(), rand());
	cout.precision(16);
	//system("del /Q .\\outputs\\*.png");

	//testCopyConstructor();
	//testSettle();
	//testInsertBifurcations();
	//testLengthAnalysis();
	//multipleLengths();
	//testHeirarchicalSwaps();
	//heirarchicalSwapAnalysis();
	//generalAnalysis();
	//testIntroduceBifurcationsSpeedup();
	//testSwapConditions();
	//swappingSequence();
	//gradualGrowthFromHeartWithSwaps();
	//asymTheory();
	//evenRandomSpacingTest();
	//annealingTest();
	//exhaustiveSearchTest();
	//exhaustiveTrendsMem();
	//exhaustiveProblemSettle();
	//rankedExchangeConnectivity();
	//rescaledDerivativesAndAnalysis();
	//rescaledDerivatives();
	//capillariesPerArea();
	//exhaustiveGreed();
	//greedySuccessRate();
	//randomlySelectedConfigurations();
	//exhaustiveOptimal();
	//tipSwapping();
	//testTorres();
	//testTerseSeg();
	//testRescaling();
	//testConsolidation();
	//testExhaustiveConsolidated();
	//greedyConsolidatedSuccessRate();
	//greedyConsolidatedCompetition();
	//exhaustiveConsolidatedMeasures(); // use for meas vs. rank
	//graftScaling();
	//lambdaLGammaOpt();
	//lambdaLGammaOptTips(); // with tipRatioThresh filter before settling
	//lambdaLGammaOptNoSearch();
	//optConsistency();
	//drawLegend();
	//testOverrelaxation();
	//test3D();
	//capillariesPerVolume();
	//localSingleGraftTest();
	//testRegraftSingle();
	//localSingleGraftDebug();
	//makeStringTest();
	//initializeEvenRandomSpacingSingleBendTest();
	//asymArea();
	//regularLattice();
	//circularSym();
	//lambdaLGammaOptCirc();
	//lambdaLOptCircConst();
	//lambdaLOptCircIgnoreTips();
	//settleFig();
	//lambdaL3D();
	//lambdaLSphere();
	//allSearchedConfig();
	//testSimilarity();
	//bestGuessRank();
	//bestSymGuessRank();
	//allSearchedConfigPathAndRandom(); // tracks paths and similarities; includes averaging
	//scalingDataForGammaNoSearch();
	//testTrimDegenerateSinks();
	//ensembleTrimDegenerateSinks();
	//optRadiusTuning();
	//testTerseSeg();
	//testWeightedFermatSettle();
	//stupidTest();
	//testRegraftForSymmetry();
	//regraftForSymmetryEnsemble();
	//testGammaAnalysis();
	//testRecordReadNetwork();
	//plainGammaDistribution();
	//prebinnedDistributions();
	//
	if (argc == 7)
	{
		cout << "Hello" << endl;
		singleNetOptRec(atof(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), atoi(argv[6])); //(1,9,(1),1,1,0)
	}
	else{
		cout << "\n Required double precision input for singleNetOptRec(lengthCoef, avePathCoef, maxPathCoef, tipRatioThresh, side, shape)" // side is 
			<< "\n\t with shape as 0=circle 1=square 2=elongated." << endl;
	}


	time_t endTime = time(NULL);
	cout << "\n\nRuntime: " << niceTime(difftime(endTime, startTime));
	cout << "\nProgram complete.\n";
	string ending;
	system("pause");
	return 0;
}
