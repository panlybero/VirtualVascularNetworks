#pragma once
#include<iostream>
struct Point
{
	double x;
	double y;
	double z=0;
	double rad = 0; // holds radius for vessel that ENDS there
	bool hasZ;

	Point()
	{

	}
	Point(Point* other)
	{
		x = other->x;
		y = other->y;
		rad = other->rad;
	}
	Point(double _x, double _y, double _z)
	{
		x = _x;
		y = _y;
		z = _z;
		hasZ = true;
	}
	Point(double _x, double _y)
	{
		x = _x;
		y = _y;
		hasZ = false;
		
	}
	double radius()
	{
		return rad;
	}

	void dump()
	{
		std::cout << "(" << x << "," << y << "," << rad << ")";
	}
};