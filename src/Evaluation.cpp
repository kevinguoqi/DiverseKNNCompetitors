/*
 * Evaluation.cpp
 *
 *  Created on: Oct 5, 2016
 *      Author: zhengyuxin
 */

#include "Evaluation.h"
#include <cmath>
#include <limits>
#include <iostream>

using namespace std;




Evaluation::Evaluation() {
	// TODO Auto-generated constructor stub

}

Evaluation::~Evaluation() {
	// TODO Auto-generated destructor stub
}

double Evaluation::computeOriDiv( const Point& query, const vector<Point>& resultset, const DistanceFunction& func  )
{
	// initialize the direction for result
	int dim = query.attributes.size();
	vector<double> resultDirection( dim, 0.0 );

	for( Point p: resultset )
	{
		vector<double> direction( dim, 0.0 );
		for ( int i = 0; i < dim; i++ )
		{
			direction[i] = p.attributes[i] - query.attributes[i];
		}
		double norm = func.getL2Norm( direction );

		if ( norm < EPSILON )
		{
			continue;
		}

		for ( int i = 0; i < dim; i++ )
		{
			resultDirection[i] += direction[i] / norm;
		}
	}

	double score = 1 - ( func.getL2Norm( resultDirection ) / ( resultset.size() * 1.0 ) );
	resultDirection.clear();
	return score;
}



double Evaluation::computeRel( const vector<double>& resultSetSpatialDist, const vector<double>& groundTruthSpatialDist )
{
	int num = min( resultSetSpatialDist.size(), groundTruthSpatialDist.size() );
	double numerator = 0.0;
	double denominator = 0.0;

	for ( int i = 0; i < num; i++ )
	{
		numerator += groundTruthSpatialDist[i];
	}

	for ( int i = 0; i < num; i++ )
	{
		denominator += resultSetSpatialDist[i];
	}
	return numerator / denominator;
}


double Evaluation::computeDivRel( const double div, const double rel, const double alpha )
{
	return alpha * div + (1 - alpha) * rel;
}


double Evaluation::computeAngDiv( const Point& query, const vector<Point>& resultset, const DistanceFunction& func  )
{
	double maxCosine = -1;
	int dim = resultset[0].attributes.size();
	int size = resultset.size();
	for ( int i = 0; i < size - 1; i++ )
	{
		Point p1 = resultset[i];
		vector<double> v1( dim, 0.0 );
		for ( int d = 0; d < dim; d++ )
		{
			v1[d] = p1.attributes[d] - query.attributes[d];
		}

		for ( int j = i + 1; j < size; j++ )
		{
			Point p2 = resultset[j];
			vector<double> v2( dim, 0.0 );
			for ( int d = 0; d < dim; d++ )
			{
				v2[d] = p2.attributes[d] - query.attributes[d];
			}
			double cosV1V2 = func.cosine(v1, v2);

			if ( maxCosine < cosV1V2 )
			{
				maxCosine = cosV1V2;
			}
			v2.clear();
		}
		v1.clear();
	}

	if ( maxCosine < ( 1 + EPSILON ) && maxCosine > ( 1 - EPSILON )	)
	{
		return 0;
	}
	else
	{
		return acos(maxCosine) * 180.0 / PI;
	}
}
