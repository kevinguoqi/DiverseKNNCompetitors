/*
 * Dataset.h
 *
 *  Created on: Oct 4, 2016
 *      Author: zhengyuxin
 */

#ifndef DATASET_H_
#define DATASET_H_

#include "Point.h"
#include <vector>
#include <string>
#include <climits>

using namespace std;




class Dataset {

public:

	int dim;
	vector<Point> dataset;
	vector<Point> queryset;

	vector<int> spatialIndex;
	vector<int> diverseIndex;


public:
	Dataset(int dim);
	virtual ~Dataset();


	void setDataSet( vector<Point>& dataset );
	void setQuerySet( vector<Point>& queryset );

	void setSpatialIndex( vector<int>& spatialIndex );
	void setDiverseIndex( vector<int>& diverseIndex );

	// read data set
	void readDataSet( string& inputFile, int num, char sep, bool skipFirstColumn );
	// read query set
	void readQuerySet( string& inputFile, int num, char sep, bool skipFirstColumn );


private:
	void readFile( string& inputFile, vector<Point>& set, int num, char sep );
	void readFileSkipFirstColumn( string& inputFile, vector<Point> &set, int num, char sep );

};

#endif /* DATASET_H_ */
