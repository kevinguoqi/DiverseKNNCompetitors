/*
 * DistanceFunction.h
 *
 *  Created on: Oct 4, 2016
 *      Author: zhengyuxin
 */

#ifndef DISTANCEFUNCTION_H_
#define DISTANCEFUNCTION_H_

#include "Point.h"
#include "WeightedFunction.h"
#include <vector>

using namespace std;

class DistanceFunction {
public:
	DistanceFunction();
	virtual ~DistanceFunction();


public:
	double getSpatialDistance( const Point& p1, const Point& p2, const vector<int>& spatialIndex ) const;
	double getSpatialDistanceSquare( const Point& p1, const Point& p2, const vector<int>& spatialIndex ) const;
	double getDiverseDistance( const Point& p1, const Point& p2, const vector<int>& diverseIndex, const WeightedFunction& weightedFunction) const;

	double getL2Norm( const vector<double>& vec ) const;
	double dotProduct( const vector<double>& v1, const vector<double>& v2 ) const;
	double cosine( const vector<double>& v1, const vector<double>& v2 ) const;

};

#endif /* DISTANCEFUNCTION_H_ */
