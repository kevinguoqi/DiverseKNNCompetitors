/*
 * KNDN.h
 *
 *  Created on: Oct 4, 2016
 *      Author: zhengyuxin
 */

#ifndef KNDN_H_
#define KNDN_H_


#include "Dataset.h"
#include "DistanceFunction.h"
#include "Evaluation.h"
#include <unordered_set>
#include <vector>
#include <unordered_map>

using namespace std;


class KNDN {

public:
	KNDN();
	virtual ~KNDN();


public:
	void greedyAlgorithm( Dataset& dataset,
			Evaluation& evaluation, DistanceFunction& distanceFunction, WeightedFunction& weightedFunction,
			double diverseThreshold, int topk, double alpha  );

	void bufferedAlgorithm( Dataset& dataset,
			Evaluation& evaluation, DistanceFunction& distanceFunction, WeightedFunction& weightedFunction,
			double diverseThreshold, int topk, double alpha );


	void GabrielNeighborAlgorithm( Dataset& dataset,
			Evaluation& evaluation, DistanceFunction& distanceFunction, WeightedFunction& weightedFunction,
			double diverseThreshold, int topk, double alpha );

	void maxMinDistAlgorithm( Dataset& dataset,
				Evaluation& evaluation, DistanceFunction& distanceFunction, WeightedFunction& weightedFunction,
				double diverseThreshold, int topk, double alpha, int candidateKNNSize );

private:
	unordered_set<int> optimizeLayer( unordered_map<int, Point>& dataMap, Point& query,
			unordered_set<int>& resultIDs, unordered_set<int>& thislayer, int topk,
			Evaluation& evaluation, DistanceFunction& distanceFunction  );

	double adjustDiverseThreshold( double diverseThreshold, Dataset& dataset, unordered_map<int, Point>& dataMap, vector<std::pair<int, double>>& knnDistVec,
			int topk, DistanceFunction& distanceFunction, WeightedFunction& weightedFunction  );


	int selectNextResultIndex( const Dataset& dataset, unordered_map<int, Point>& dataMap, const vector<std::pair<int, double>>& spatialDistVec,
			const unordered_set<int>& resultIDs, const int candidateKNNSize, const DistanceFunction& distanceFunction );

	void writeResult(const string& filename, const vector<int>& resultid);
};

#endif /* KNDN_H_ */
