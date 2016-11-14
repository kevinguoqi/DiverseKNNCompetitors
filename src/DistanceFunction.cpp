/*
 * DistanceFunction.cpp
 *
 *  Created on: Oct 4, 2016
 *      Author: zhengyuxin
 */

#include "DistanceFunction.h"
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <iostream>

using namespace std;

DistanceFunction::DistanceFunction() {
	// TODO Auto-generated constructor stub

}

DistanceFunction::~DistanceFunction() {
	// TODO Auto-generated destructor stub
}


double DistanceFunction::getSpatialDistance( const Point& p1, const Point& p2, const vector<int>& spatialIndex ) const
{
	double result = 0.0;
	for( int index: spatialIndex )
	{
		double diff = p1.attributes[index] - p2.attributes[index];
		result += diff * diff;
	}
	return pow(result, 0.5);
}


double DistanceFunction::getSpatialDistanceSquare( const Point& p1, const Point& p2, const vector<int>& spatialIndex ) const
{
	double result = 0.0;
	for( int index: spatialIndex )
	{
		double diff = p1.attributes[index] - p2.attributes[index];
		result += diff * diff;
	}
	return result;
}


double DistanceFunction::getDiverseDistance( const Point& p1, const Point& p2, const vector<int>& diverseIndex, const WeightedFunction& weightedFunction) const
{
	vector< double > diffVec;

	for( int index: diverseIndex )
	{
		diffVec.push_back( abs( p1.attributes[index] - p2.attributes[index] ) );
	}

	// sort the diff vector in descending order
	std::sort( diffVec.begin(), diffVec.end(), std::greater<double>() );

	// calculate the diverse distance
	double result = 0.0;
	for ( int i = 0; i < weightedFunction.dim; i++ )
	{
		result += diffVec[i] * weightedFunction.weights[i];
	}

	diffVec.clear();
	return result;
}


double DistanceFunction::getL2Norm( const vector<double>& vec ) const
{
	double result = 0.0;
	for( double v: vec )
	{
		result += v * v;
	}
	return pow(result, 0.5);
}

double DistanceFunction::dotProduct( const vector<double>& v1, const vector<double>& v2 ) const
{
	double result = 0.0;
	for ( unsigned int i = 0; i < v1.size(); i++ )
	{
		result += v1[i] * v2[i];
	}
	return result;
}

double DistanceFunction::cosine( const vector<double>& v1, const vector<double>& v2 ) const
{
	return dotProduct(v1, v2) / (getL2Norm(v1) * getL2Norm(v2));
}
