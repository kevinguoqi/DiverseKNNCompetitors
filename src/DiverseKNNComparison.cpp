//============================================================================
// Name        : DiverseKNNComparison.cpp
// Author      : zyx
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "Dataset.h"
#include "KNDN.h"
#include <string>
#include <math.h>


using namespace std;


void testCaldata()
{
	// set parameters
	int dim = 512;
	char seperator = ' ';

//	bool skipFirstColumn = true;
//	string datafile = "/mnt/hd4/home/zhengyuxin/workspace/diverseKnnDataSet/caldata/caldata_data";
//	string queryfile =	"/mnt/hd4/home/zhengyuxin/workspace/diverseKnnDataSet/caldata/caldata_query_test";

	bool skipFirstColumn = true;
//	string datafile = "/mnt/hd4/home/zhengyuxin/workspace/diverseKnnDataSet/motor_normalize.txt";
//	string queryfile = "/mnt/hd4/home/zhengyuxin/workspace/diverseKnnDataSet/motor_query_10.txt";

//	string datafile = "/mnt/hd4/home/zhengyuxin/workspace/diverseKnnDataSet/CALdataset_test.txt";
//	string queryfile = "/mnt/hd4/home/zhengyuxin/workspace/diverseKnnDataSet/CALquery.txt";

//	string datafile = "/mnt/hd4/home/guoqi/div_exp/dataset/motor_normalize.txt";
//	string queryfile = "/mnt/hd4/home/guoqi/div_exp/dataset/motor_query.txt";

//	string datafile = "/mnt/hd4/home/guoqi/div_exp/dataset/sun.ds_normalize";
//	string queryfile = "/mnt/hd4/home/guoqi/div_exp/dataset/sun_500.q_normalize";

//	string datafile = "/mnt/hd4/home/guoqi/div_exp/dataset/mnist.train_normalize";
//	string queryfile = "/mnt/hd4/home/guoqi/div_exp/dataset/mnist_500.test_normalize";

//	string datafile = "/mnt/hd4/home/guoqi/div_exp/dataset/sgp_noid_4f_uniq.txt";
//	string queryfile = "/mnt/hd4/home/guoqi/div_exp/dataset/sgp_query_500.txt";

	string datafile = "/mnt/hd4/home/guoqi/div_exp/dataset/coil100_gist.ds";
	string queryfile = "/mnt/hd4/home/guoqi/div_exp/dataset/coil100_gist.ds";

	double diverseThreshold = 0.1;
//	int topk[] = {5, 7, 9, 11, 13, 15};
//	int topk[] = {5, 10, 20, 40, 80, 160 };
	int topk[] = {10};
	double linearWeight = 0.5;
	double alpha = 0.1;
//	int KNN = topk * 5;

	Evaluation evaluation;
	DistanceFunction distanceFunction;
	WeightedFunction weightedFunction( dim, alpha );

	int dataNum = INT_MAX;
	int queryNum = INT_MAX;


	// load the dataset
	Dataset dataset(dim);
	dataset.readDataSet(datafile, dataNum, seperator, skipFirstColumn);
	dataset.readQuerySet(queryfile, queryNum, seperator, skipFirstColumn);

	// set index, just use the same index
	for (int i = 0; i < dim; i++) {
		dataset.spatialIndex.push_back(i);
		dataset.diverseIndex.push_back(i);
	}

	KNDN kndn;

	for(int i=0;i<sizeof(topk)/sizeof(int);i++){

		cout << "#################################################################################################" << endl;
		cout << "current topk is " << topk[i] << endl;
		cout << "#################################################################################################" << endl;

		int KNN = topk[i] * 5;
//		cout << "running greedy algorithm" << endl;
//		kndn.greedyAlgorithm( dataset, evaluation, distanceFunction, weightedFunction, diverseThreshold, topk[i], linearWeight );
//		cout << endl;

		cout << "running buffered algorithm" << endl;
		kndn.bufferedAlgorithm( dataset, evaluation, distanceFunction, weightedFunction, diverseThreshold, topk[i], linearWeight );
		cout << endl;

		cout << "running Gabriel Neighbor algorithm" << endl;
		kndn.GabrielNeighborAlgorithm( dataset, evaluation, distanceFunction, weightedFunction, diverseThreshold, topk[i], linearWeight );
		cout << endl;

		cout << "running MaxMinDist algorithm" << endl;
		kndn.maxMinDistAlgorithm( dataset, evaluation, distanceFunction, weightedFunction, diverseThreshold, topk[i], linearWeight, KNN );
		cout << endl;

	}

}


void testAcos()
{
	std::cout << "acos(-1) = " << acos(-1) << '\n'
	              << "acos(0.0) = " << acos(0.0) << " 2*acos(0.0) = " << 2*acos(0) << '\n'
	              << "acos(0.5) = " << acos(0.5) << " 3*acos(0.5) = " << 3*acos(0.5) << '\n'
	              << "acos(1) = " << acos(1.0) << '\n';

	std::cout << "acos(1.00000000001) = " << acos(1.00000000001) << '\n';
	std::cout << "acos(1.1) = " << acos(1.1) << '\n';
}


int main() {


//	testAcos();

	testCaldata();


	return 0;
}
