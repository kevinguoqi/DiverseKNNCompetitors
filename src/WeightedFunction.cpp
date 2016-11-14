/*
 * WeightedFunction.cpp
 *
 *  Created on: Oct 4, 2016
 *      Author: zhengyuxin
 */

#include "WeightedFunction.h"
#include <math.h>

WeightedFunction::WeightedFunction( int dim, double alpha ) {

	this->dim = dim;
	this->alpha = alpha;

	for ( int i = 0; i < dim; i++ )
	{
		double weight = pow( alpha, i ) * ( 1 - alpha ) / ( 1 - pow( alpha, dim ) );
		weights.push_back( weight );
	}

}

WeightedFunction::~WeightedFunction() {
	// TODO Auto-generated destructor stub
}

