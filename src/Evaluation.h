/*
 * Evaluation.h
 *
 *  Created on: Oct 5, 2016
 *      Author: zhengyuxin
 */

#ifndef EVALUATION_H_
#define EVALUATION_H_

#include "Point.h"
#include <vector>
#include "DistanceFunction.h"


#define PI 3.14159265
#define EPSILON 0.0000000001

using namespace std;

class Evaluation {
public:
	Evaluation();
	virtual ~Evaluation();


public:
	double computeOriDiv( const Point& query, const vector<Point>& resultset, const DistanceFunction& func );
	double computeRel( const vector<double>& resultSetSpatialDist, const vector<double>& groundTruthSpatialDist );
	double computeDivRel( const double div, const double rel, const double alpha );
	double computeAngDiv( const Point& query, const vector<Point>& resultset, const DistanceFunction& func );

};

#endif /* EVALUATION_H_ */
