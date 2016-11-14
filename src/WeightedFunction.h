/*
 * WeightedFunction.h
 *
 *  Created on: Oct 4, 2016
 *      Author: zhengyuxin
 */

#ifndef WEIGHTEDFUNCTION_H_
#define WEIGHTEDFUNCTION_H_

#include <vector>
using namespace std;

class WeightedFunction {


public:
	int dim;
	double alpha;
	vector<double> weights;

public:
	WeightedFunction( int dim, double alpha );
	virtual ~WeightedFunction();


};

#endif /* WEIGHTEDFUNCTION_H_ */
