/*
 * KNDN.cpp
 *
 *  Created on: Oct 4, 2016
 *      Author: zhengyuxin
 */

#include "KNDN.h"
#include <algorithm>
#include <iostream>
#include "IDDistance.h"
#include "Point.h"
#include "Timer.h"
#include <time.h>
#include <stdlib.h>
#include <fstream>

using namespace std;


#define double double

KNDN::KNDN() {
	// TODO Auto-generated constructor stub

}

KNDN::~KNDN() {
	// TODO Auto-generated destructor stub
}


bool comparePairs( const std::pair<int, double>& lhs, const std::pair<int, double>& rhs ){
	return lhs.second < rhs.second;
}


vector<std::pair<int, double>> computeDistanceAndSort( Dataset& dataset,  Point& query, DistanceFunction& distanceFunction )
{
	// calculate the spatial distance for every pair
	vector<std::pair<int, double>> spatialDistVec;
	int datapointid = 0;
	for ( Point p: dataset.dataset )
	{
		double spatialDistance = distanceFunction.getSpatialDistanceSquare( query, p, dataset.spatialIndex );
		std::pair<int, double> idDistance( datapointid, spatialDistance );
		spatialDistVec.push_back( idDistance );
		datapointid ++;
	}

	std::sort( spatialDistVec.begin(), spatialDistVec.end(), comparePairs );
	return spatialDistVec;
}

void KNDN::writeResult(const string& filename, const vector<int>& resultid){

	ofstream myfile( filename.c_str(), fstream::out | fstream::app );

	if(myfile.is_open()){
		myfile << endl;
		for(int pid: resultid){
			myfile << (pid+1) << "," ;
		}

		myfile.close();
	}


}


void KNDN::greedyAlgorithm( Dataset& dataset,
		Evaluation& evaluation, DistanceFunction& distanceFunction, WeightedFunction& weightedFunction,
		double diverseThreshold, int topk, double alpha )
{
	cout << "************************************ Greedy algorithm results ******************************************************" << endl;

	// build an index for the dataset
	unordered_map<int, Point> dataMap;
	int pointid = 0;
	for ( Point p: dataset.dataset )
	{
		dataMap[pointid] = p;
		pointid ++;
	}

	// get start time
	Timer timer;

	// store the final results
	vector<double> divVec;
	vector<double> relVec;
	vector<double> divRelVec;
	vector<double> angDivVec;

	// calculate for every query
	int query_count = 0;
	for ( Point query: dataset.queryset )
	{
		// maintain a result list for every query
		vector<Point> result;
		vector<int> resultid;

		for( int i = 0; i < query.attributes.size(); i++ )
		{
			dataset.dataset[query_count].attributes[i] = 1000;
		}

		// calculate the spatial distance for every pair
		vector<std::pair<int, double>> spatialDistVec = computeDistanceAndSort( dataset,  query, distanceFunction );


		// calculate the ground truth
		vector<double> groundTruthDist( topk, 0.0 );
		vector<double> resultSetDist;

		diverseThreshold = adjustDiverseThreshold( diverseThreshold, dataset, dataMap, spatialDistVec, topk, distanceFunction, weightedFunction );


		// calculate the diverse set
		// the 1NN must be in the result set
		int query_in_dataset = 0;
		result.push_back( dataMap[spatialDistVec[ (0+query_in_dataset) ].first] );
		resultSetDist.push_back( pow( spatialDistVec[ (0+query_in_dataset)].second, 0.5 ) );
		resultid.push_back( spatialDistVec[ (0+query_in_dataset) ].first );


		for ( unsigned int i = 1+query_in_dataset; i < spatialDistVec.size(); i++ )
		{
			pair<int, double> pair = spatialDistVec[i];
			int pid = pair.first;
			Point currentPoint = dataMap[pid];


			// a point is a diverse point if it's diversity is greater than the diverseThreshold to all tuples in the result set
			bool isResult = true;
			for ( Point resultPoint: result )
			{
				if ( distanceFunction.getDiverseDistance( resultPoint, currentPoint, dataset.diverseIndex, weightedFunction ) < diverseThreshold )
				{
					isResult = false;
					break;
				}
			}
			if ( isResult )
			{
				result.push_back( currentPoint );
				resultSetDist.push_back( pow( pair.second, 0.5 ));
				resultid.push_back( pid );
				// stop when we have enough results
				if ( result.size() >= topk )
				{
					break;
				}
			}
		}

		string resultprefix = "/mnt/hd4/home/guoqi/coil100/coil100_all/query" + to_string(query_count) + "_diverse.csv";
		writeResult(resultprefix, resultid);


		double div = evaluation.computeOriDiv( query, result, distanceFunction ) ;
		double rel = evaluation.computeRel( resultSetDist, groundTruthDist );
		double angDiv = evaluation.computeAngDiv( query, result, distanceFunction );

		divVec.push_back( div );
		relVec.push_back( rel );
		divRelVec.push_back( evaluation.computeDivRel( div, rel, alpha ) );
		angDivVec.push_back( angDiv );


		for( int i = 0; i < query.attributes.size(); i++ )
		{
			dataset.dataset[query_count].attributes[i] = query.attributes[i];
		}
		query_count++;

	}

	// calculate and print the average
	double avgDiv = 0.0;
	double avgRel = 0.0;
	double avgDivRel = 0.0;
	double avgAngDiv = 0.0;
	int querySize = dataset.queryset.size();

	for ( int i = 0; i < querySize; i++ )
	{
		avgDiv += divVec[i];
		avgRel += relVec[i];
		avgDivRel += divRelVec[i];
	    avgAngDiv += angDivVec[i];
	}


	double running_time = timer.elapsed();


	cout << "running time: " << running_time << endl;
	cout << "average div: " << avgDiv / (querySize * 1.0) << endl;
	cout << "average rel: " << avgRel / (querySize * 1.0) << endl;
	cout << "average divRel: " << avgDivRel / (querySize * 1.0) << endl;
	cout << "average angle div: " << avgAngDiv / (querySize * 1.0) << endl;

}

int getNumPoints( vector< vector<Point> >& followers )
{
	int count = 0;
	for ( unsigned int i = 0; i < followers.size(); i++ )
	{
		for ( unsigned int j = 0; j < followers[i].size(); j++ )
		{
			count++;
		}
	}
	return count;
}


void KNDN::bufferedAlgorithm( Dataset& dataset,
		Evaluation& evaluation, DistanceFunction& distanceFunction, WeightedFunction& weightedFunction,
		double diverseThreshold, int topk, double alpha )
{
	cout << "************************************ Buffered algorithm results ******************************************************" << endl;

	// build an index for the dataset
	unordered_map<int, Point> dataMap;
	int pointid = 0;
	for (Point p : dataset.dataset) {
		dataMap[pointid] = p;
		pointid++;
	}

	// start the timer
	Timer timer;

	// store the final results
	vector<double> divVec;
	vector<double> relVec;
	vector<double> divRelVec;
	vector<double> angDivVec;


	// compute the non-diverse radius
//	double nonDiverseRadius = distanceFunction.getNotDiverseRadius( weightedFunction, diverseThreshold );

	int query_count = 0;
	// calculate for every query
	for (Point query : dataset.queryset) {

		for( int i = 0; i < query.attributes.size(); i++ )
		{
				dataset.dataset[query_count].attributes[i] = 1000;
		}


//		cout << "query " << queryid++ << ": [" << query.attributes[0] << " ," << query.attributes[dataset.dim - 1 ] << "]" << endl;

		// maintain a result list for every query
		vector<Point> result;
		vector<int> resultid;
		// maintain a follower set for replacement
		vector< vector<Point> > followers;



		// calculate the spatial distance for every pair
		vector<std::pair<int, double>> spatialDistVec = computeDistanceAndSort( dataset, query, distanceFunction );
		// calculate the ground truth
		vector<double> groundTruthDist(topk, 0.0);
		vector<double> resultSetDist;

		for (int i = 0; i < topk; i++) {
			groundTruthDist[i] = pow( spatialDistVec[i].second, 0.5 );
		}

		diverseThreshold = adjustDiverseThreshold( diverseThreshold, dataset, dataMap, spatialDistVec, topk, distanceFunction, weightedFunction );

		// calculate the diverse set
		// the 1NN must be in the result set
		result.push_back( dataMap[spatialDistVec[0].first] );
		//resultid.push_back( query_count );
		followers.push_back( vector<Point>() );

		int knn = 1;
		while( result.size() < topk && knn < dataset.dataset.size() )
		{
			pair<int, double> pair = spatialDistVec[knn];
			int pid = pair.first;
			Point currentPoint = dataMap[pid];

//			cout << "k = " << knn << endl;
//			cout << "pid = " << pid << endl;
//			cout << "result size = " << result.size() << endl;
//			cout << "follower size = " << getNumPoints(followers) << endl << endl;


			// check which leaders is not diverse for the currentPoint
			vector<int> notDiverseLeaders;
			// a point is a diverse point if it's diversity is greater than the diverseThreshold to all tuples in the result set
			bool isResult = true;
			// get the original result size
			int resultSize = result.size();
			for ( unsigned int i = 0; i < resultSize; i++ )
			{
				if ( distanceFunction.getDiverseDistance( result[i], currentPoint, dataset.diverseIndex, weightedFunction ) < diverseThreshold )
				{
					isResult = false;
					notDiverseLeaders.push_back( i );
					if ( notDiverseLeaders.size() >= 2 )
					{
						break;
					}
				}
			}


			// if the currentPoint is diverse, put it in the result
			if ( isResult )
			{
				// the currentPoint becomes a new leader
				result.push_back( currentPoint );
				//resultid.push_back( pid );
				followers.push_back( vector<Point>() );

				// remove all the non-dedicated followers
				for ( int i = 0; i < resultSize; i++ )
				{
					for ( int j = ( followers[i].size() - 1 ); j >=0; j-- )
					{
						// if the follower is not diverse w.r.t. the new leader,
						// then it is not a dedicated follower, remove it
						if ( distanceFunction.getDiverseDistance( followers[i][j], currentPoint, dataset.diverseIndex, weightedFunction) < diverseThreshold )
						{
							followers[i].erase( followers[i].begin() + j );
						}
					}
				}
			}
			// else, check whether we can replace the leader from its followers
			// if the currentPoint is not diverse w.r.t. only one leader, add the currentPoint to the follower
			else if ( notDiverseLeaders.size() == 1 )
			{
				int replacedLeaderCandidateID = notDiverseLeaders[0];

				// since a new follower is added, check whether we can find exclusive followers
				bool foundExclusiveFollower = false;
				vector<Point> newLeaders;
				vector<Point> remainingFollowers;

				for ( unsigned int j = 0; j < followers[replacedLeaderCandidateID].size(); j++)
				{
					// if we can find another dedicated follower, put it into newLeader
					if ( distanceFunction.getDiverseDistance( followers[replacedLeaderCandidateID][j], currentPoint, dataset.diverseIndex, weightedFunction) > diverseThreshold )
					{
						foundExclusiveFollower = true;
						newLeaders.push_back( followers[replacedLeaderCandidateID][j] );
						newLeaders.push_back( currentPoint );

						followers[replacedLeaderCandidateID].erase( followers[replacedLeaderCandidateID].begin() + j );
						remainingFollowers = followers[replacedLeaderCandidateID];


						// delete the corresponding leader and follower
						result.erase( result.begin() + replacedLeaderCandidateID );
						//resultid.erase( resultid.begin() + replacedLeaderCandidateID );

						followers.erase( followers.begin() + replacedLeaderCandidateID );

						break;
					}
				}


				// if we find exclusiVe followers, do the replacement of leaders
				if ( foundExclusiveFollower )
				{
					// check whether the existing followers is diverse with the new leaders
					for ( unsigned int i = 0; i < followers.size(); i++)
					{
						for (int j = followers[i].size() - 1; j >= 0; j--)
						{
							for (Point newLeader : newLeaders)
							{
								// check whether the follower is not diverse with the new leader,
								// if it is, remove the follower
								if ( distanceFunction.getDiverseDistance( followers[i][j], newLeader, dataset.diverseIndex, weightedFunction) < diverseThreshold )
								{
									followers[i].erase( followers[i].begin() + j);
									break;
								}
							}
						}
					}
					// get the start index of the new leader
					int newLeaderStartIndex = result.size();

					// add the newHeaders to the result
					for (Point newLeader : newLeaders)
					{
						result.push_back(newLeader);
						followers.push_back(vector<Point>());
					}

					// check the remaining Followers
					for (int i = remainingFollowers.size() - 1; i >= 0;	i--)
					{
						// check the follower is dedicated or not
						Point currentFollower = remainingFollowers[i];

						vector<int> notDiverseLeaders4Follower;
						int index = newLeaderStartIndex;
						for (Point newLeader : newLeaders)
						{
							// check whether the follower is not diverse with the new leader,
							if ( distanceFunction.getDiverseDistance(currentFollower, newLeader, dataset.diverseIndex, weightedFunction) < diverseThreshold ) {
								notDiverseLeaders4Follower.push_back(index);
							}
							index++;
						}

						// if the remaining follower is a dedicated follower, add it to the corresponding leader
						if (notDiverseLeaders4Follower.size() == 1) {
							followers[notDiverseLeaders4Follower[0]].push_back(	currentFollower );
						}
						notDiverseLeaders4Follower.clear();
					}


				}
				// if we can not find another exclusive follower, add the currentPoint to the follower set
				else
				{
					if ( followers[replacedLeaderCandidateID].size() < topk )
					{
						followers[replacedLeaderCandidateID].push_back( currentPoint );
					}
				}

				// clear the newLeaders vector and the remainingFollowers vector
				newLeaders.clear();
				remainingFollowers.clear();
			}

			notDiverseLeaders.clear();
			knn++;
		}


		// calculate the relevance
		for ( int i = 0; i < result.size(); i++ )
		{
			resultSetDist.push_back( distanceFunction.getSpatialDistance( query, result[i], dataset.spatialIndex) );
		}

		double div = evaluation.computeOriDiv(query, result, distanceFunction);
		double rel = evaluation.computeRel(resultSetDist, groundTruthDist);
		double angDiv = evaluation.computeAngDiv(query, result, distanceFunction);

		divVec.push_back(div);
		relVec.push_back(rel);
		divRelVec.push_back(evaluation.computeDivRel(div, rel, alpha));
		angDivVec.push_back(angDiv);

//		cout << "size = " << result.size() << " , div = "<< div << " , rel = " << rel << " , angDiv = " << angDiv << endl;

		for (Point res: result){

			for(int i=0; i < pointid; i++){

				bool isCurrent = true;
				Point current = dataMap[i];

				for(int j=0; j < 512; j++){

					if(res.attributes[j] != current.attributes[j] ){
						isCurrent = false;
						break;
					}

				}

				if(isCurrent){

					//cout << i << ",";
					resultid.push_back(i);
				}

			}


		}

		string resultprefix = "/mnt/hd4/home/guoqi/coil100/coil100_all/query" + to_string(query_count) + "_diverse.csv";
		writeResult(resultprefix, resultid);



		// clear the results
		result.clear();
		followers.clear();
		spatialDistVec.clear();
		groundTruthDist.clear();
		resultSetDist.clear();

		for( int i = 0; i < query.attributes.size(); i++ )
		{
				dataset.dataset[query_count].attributes[i] = query.attributes[i];
		}
		query_count++;
	}

	// calculate and print the average
	double avgDiv = 0.0;
	double avgRel = 0.0;
	double avgDivRel = 0.0;
	double avgAngDiv = 0.0;
	int querySize = dataset.queryset.size();

	for (int i = 0; i < querySize; i++) {
		avgDiv += divVec[i];
		avgRel += relVec[i];
		avgDivRel += divRelVec[i];
		avgAngDiv += angDivVec[i];
	}

	double running_time = timer.elapsed();


	cout << "running time: " << running_time << endl;
	cout << "average div: " << avgDiv / (querySize * 1.0) << endl;
	cout << "average rel: " << avgRel / (querySize * 1.0) << endl;
	cout << "average divRel: " << avgDivRel / (querySize * 1.0) << endl;
	cout << "average angle div: " << avgAngDiv / (querySize * 1.0) << endl;
}



unordered_set<int> KNDN::optimizeLayer( unordered_map<int, Point>& dataMap, Point& query,
		unordered_set<int>& resultIDs, unordered_set<int>& thislayer, int topk,
		Evaluation& evaluation, DistanceFunction& distanceFunction  )
{
	unordered_set<int> selected;
	selected.clear();
	unordered_set<int>::iterator it;
	vector<Point> results;

	int requiredSize = topk - resultIDs.size();
	float* tmpdiv;
	float tmpdivi;
	// select new point from this layer until a total of k points are selected
	while ( selected.size() < requiredSize )
	{
		results.clear();

		for ( int id : resultIDs )
		{
			results.push_back( dataMap[id] );
		}
		for ( int id : selected )
		{
			results.push_back( dataMap[id] );
		}


		// calculate angular diversity of individual results
		tmpdiv = new float[thislayer.size()];
		int tmpind = 0;
		for ( it = thislayer.begin(); it != thislayer.end(); it++)
		{
			// add one and calculate the diverse score
			if ( it == thislayer.begin() )
			{
				results.push_back( dataMap[*it] );
			}
			else
			{
				// delete the previous one
				results.erase( results.end() - 1 );
				results.push_back( dataMap[*it] );
			}
			// calculate the diverse score
			tmpdiv[tmpind++] = evaluation.computeOriDiv( query, results, distanceFunction );
		}

		// greedy algorithm
		// select the one with highest diversity
		tmpdivi = tmpdiv[0];
		it = thislayer.begin();
		tmpind = *it;
		for (int j = 1; j < thislayer.size(); j++)
		{
			it++;
			if (tmpdivi < tmpdiv[j]) {
				tmpdivi = tmpdiv[j];
				tmpind = *it;
			}
		}

		// include it in result set
		selected.insert( tmpind );
		thislayer.erase( thislayer.find(tmpind) );

		// clear memory
		delete[] tmpdiv;
	}
	return selected;
}


void KNDN::GabrielNeighborAlgorithm( Dataset& dataset,
		Evaluation& evaluation, DistanceFunction& distanceFunction, WeightedFunction& weightedFunction,
		double diverseThreshold, int topk, double alpha )
{

	cout << "************************************ Gabriel Neighbor algorithm results ******************************************************" << endl;

	// build an index for the dataset
	unordered_map<int, Point> dataMap;
	int pointid = 0;
	for (Point p : dataset.dataset) {
		dataMap[pointid] = p;
		pointid++;
	}


	// start the timer
	Timer timer;

	// store the final results
	vector<double> divVec;
	vector<double> relVec;
	vector<double> divRelVec;
	vector<double> angDivVec;

	int numData = dataset.dataset.size();
//	cout << "numData:" << numData << endl;


	int queryid = 0;
	int query_count = 0;
	for ( Point query: dataset.queryset )
	{
		cout << "Query " << queryid++ << " :[" << query.attributes[0] << ", " << query.attributes[1] << "]"<< endl;

		for( int i = 0; i < query.attributes.size(); i++ )
		{
			dataset.dataset[query_count].attributes[i] = 1000;
		}

		// maintain a result list for every query
		unordered_set<int> resultIDs;

		// calculate the spatial distance for every pair
		vector<std::pair<int, double>> spatialDistVec = computeDistanceAndSort( dataset,  query, distanceFunction );

		// calculate the ground truth
		vector<double> groundTruthDist( topk, 0.0 );
		for ( int i = 0; i < topk; i++ )
		{
			groundTruthDist[i] = pow( spatialDistVec[i].second, 0.5 );
		}

		// insert the 1NN to this layer
		unordered_set<int> thislayer;
		thislayer.insert( spatialDistVec[0].first );
		resultIDs.insert( spatialDistVec[0].first );
		int previousResultSize = 0;

		while ( resultIDs.size() < topk && resultIDs.size() > previousResultSize )
		{
			// reset the previousResultSize
			previousResultSize = resultIDs.size();

			unordered_set<int> lastlayer = unordered_set<int>( thislayer );
			thislayer.clear();

			// for all points from last layer
			for ( unordered_set<int>::iterator it = lastlayer.begin(); it != lastlayer.end(); it++ )
			{
				// find all Gabriel neighbors
				for ( int j = 0; j < numData; ++j)
				{
					// ignore the point itself
					if ( *it == j )
					{
						continue;
					}

					double ijlength = distanceFunction.getSpatialDistanceSquare( dataMap[*it], dataMap[j], dataset.spatialIndex );
					bool isGGedge = true;
					for ( int kk = 0; kk < numData; ++kk)
					{
						if ( kk == *it || kk == j )
						{
							continue;
						}
						// check whether it is a gabriel graph edge
						if ( ( distanceFunction.getSpatialDistanceSquare( dataMap[*it], dataMap[kk], dataset.spatialIndex ) +
								distanceFunction.getSpatialDistanceSquare( dataMap[kk], dataMap[j], dataset.spatialIndex ) ) < ijlength + EPSILON )
						{
							isGGedge = false;
							break;
						}
					}
					if ( isGGedge && resultIDs.count(j) == 0 && thislayer.count(j) == 0)
					{
						if ( distanceFunction.getSpatialDistanceSquare( query, dataMap[j], dataset.spatialIndex ) > EPSILON )
						{
							thislayer.insert(j);
						}
					}
				}
			}

			// add layer to the result set
			if ( ( resultIDs.size() + thislayer.size() ) <= topk )
			{
				// add all
				for ( int id : thislayer )
				{
					resultIDs.insert( id );
				}
			}
			else
			{
				// optimize result set
				unordered_set<int> selected = optimizeLayer( dataMap, query, resultIDs, thislayer, topk, evaluation, distanceFunction  );
				for ( int id : selected )
				{
					resultIDs.insert( id );
				}
			}
		}

		vector<Point> resultPoints;
		vector<int> resultid;
		for ( int id : resultIDs )
		{
			//cout << id << ",";
			resultPoints.push_back( dataMap[id] );
			resultid.push_back( id );
		}


		// calculate the relevance
		vector<double> resultSetDist;
		for ( int i = 0; i < resultPoints.size(); i++ )
		{
			resultSetDist.push_back( distanceFunction.getSpatialDistance( query, resultPoints[i], dataset.spatialIndex) );
		}

		double div = evaluation.computeOriDiv(query, resultPoints, distanceFunction);
		double rel = evaluation.computeRel(resultSetDist, groundTruthDist);
		double angDiv = evaluation.computeAngDiv(query, resultPoints, distanceFunction);

		divVec.push_back(div);
		relVec.push_back(rel);
		divRelVec.push_back(evaluation.computeDivRel(div, rel, alpha));
		angDivVec.push_back(angDiv);

		cout << "div = "<< div << " , rel = " << rel << " , angDiv = " << angDiv << endl;

		string resultprefix = "/mnt/hd4/home/guoqi/coil100/coil100_all/query" + to_string(queryid-1) + "_diverse.csv";
		writeResult(resultprefix, resultid);

		for( int i = 0; i < query.attributes.size(); i++ )
		{
			dataset.dataset[query_count].attributes[i] = query.attributes[i];
		}
		query_count++;

	}

	// calculate and print the average
	double avgDiv = 0.0;
	double avgRel = 0.0;
	double avgDivRel = 0.0;
	double avgAngDiv = 0.0;
	int querySize = dataset.queryset.size();

	for ( int i = 0; i < querySize; i++ )
	{
		avgDiv += divVec[i];
		avgRel += relVec[i];
		avgDivRel += divRelVec[i];
	    avgAngDiv += angDivVec[i];
	}


	double running_time = timer.elapsed();

	cout << "running time: " << running_time << endl;
	cout << "average div: " << avgDiv / (querySize * 1.0) << endl;
	cout << "average rel: " << avgRel / (querySize * 1.0) << endl;
	cout << "average divRel: " << avgDivRel / (querySize * 1.0) << endl;
	cout << "average angle div: " << avgAngDiv / (querySize * 1.0) << endl;


}

double KNDN::adjustDiverseThreshold( double diverseThreshold, Dataset& dataset, unordered_map<int, Point>& dataMap, vector<std::pair<int, double>>& knnDistVec,
		int topk, DistanceFunction& distanceFunction, WeightedFunction& weightedFunction  )
{
	Point onenn = dataMap[knnDistVec[0].first];
	double knnMaxDivDist = 0.0;
	double knnMeanDivDist = 0.0;

	for(int i = 1; i < topk; i++)
	{
		Point inn = dataMap[knnDistVec[i].first];
		double innDivDist = distanceFunction.getDiverseDistance( onenn, inn, dataset.diverseIndex, weightedFunction );
		if ( knnMaxDivDist < innDivDist )
		{
			knnMaxDivDist = innDivDist;
		}
		knnMeanDivDist += innDivDist;
	}

	if ( diverseThreshold > knnMaxDivDist)
	{
		diverseThreshold = knnMeanDivDist / (topk * 1.0);
	}
	return diverseThreshold;
}


int KNDN::selectNextResultIndex( const Dataset& dataset, unordered_map<int, Point>& dataMap, const vector<std::pair<int, double>>& spatialDistVec,
		const unordered_set<int>& resultIDs, const int candidateKNNSize, const DistanceFunction& distanceFunction )
{
	int pIndex = -1;
	double minDist = 0;

	for ( int i = 0; i < candidateKNNSize; i++ )
	{
		int currentPID = spatialDistVec[i].first;

		// if the current point is not in the result set;
		if ( resultIDs.find( currentPID ) == resultIDs.end() )
		{
			Point currentPoint = dataMap[currentPID];
			double currentMin =  std::numeric_limits<double>::max();

			// compare the currentPoint with each point in the result set
			for ( int resultId: resultIDs )
			{
				// compute the distance and update the currentMin distance
				double dist = distanceFunction.getSpatialDistanceSquare( currentPoint, dataMap[resultId], dataset.spatialIndex );
				if ( currentMin > dist )
				{
					currentMin = dist;
				}
			}

			// update the global min, we choose a point with maximized global min
			if ( minDist < currentMin )
			{
				minDist = currentMin;
				pIndex = i;
			}
		}
	}
	return pIndex;
}


void KNDN::maxMinDistAlgorithm( Dataset& dataset,
				Evaluation& evaluation, DistanceFunction& distanceFunction, WeightedFunction& weightedFunction,
				double diverseThreshold, int topk, double alpha, int candidateKNNSize )
{
	cout << "************************************ MaxMinDist algorithm results ******************************************************" << endl;

	// initialize random seed
	srand (time(NULL));

		// build an index for the dataset
		unordered_map<int, Point> dataMap;
		int pointid = 0;
		for ( Point p: dataset.dataset )
		{
			dataMap[pointid] = p;
			pointid ++;
		}

		// get start time
		Timer timer;

		// store the final results
		vector<double> divVec;
		vector<double> relVec;
		vector<double> divRelVec;
		vector<double> angDivVec;

		int query_count = 0;
		// calculate for every query
		for ( Point query: dataset.queryset )
		{
//			cout << "Query :[" << query.attributes[0] << ", " << query.attributes[1] << "]"<< endl;

			for( int i = 0; i < query.attributes.size(); i++ )
			{
				dataset.dataset[query_count].attributes[i] = 1000;
			}

			// maintain a result list for every query
			vector<Point> result;
			unordered_set<int> resultIDs;


			// calculate the spatial distance for every pair
			vector<std::pair<int, double>> spatialDistVec = computeDistanceAndSort( dataset,  query, distanceFunction );


			// calculate the ground truth
			vector<double> groundTruthDist( topk, 0.0 );
			vector<double> resultSetDist;
			for ( int i = 0; i < topk; i++ )
			{
				groundTruthDist[i] = pow( spatialDistVec[i].second, 0.5 );
			}

			diverseThreshold = adjustDiverseThreshold( diverseThreshold, dataset, dataMap, spatialDistVec, topk, distanceFunction, weightedFunction );



			// randomly select a point
			 int pIdx = rand() % candidateKNNSize;
			 resultIDs.insert( spatialDistVec[pIdx].first );
			 result.push_back( dataMap[spatialDistVec[pIdx].first] );
			 resultSetDist.push_back( pow( spatialDistVec[pIdx].second, 0.5 ) );


			 for ( int i = 1; i < topk; i++ )
			 {
				 int nextPIndex = selectNextResultIndex( dataset, dataMap, spatialDistVec, resultIDs, candidateKNNSize, distanceFunction );
				 resultIDs.insert( spatialDistVec[nextPIndex].first );
				 result.push_back( dataMap[spatialDistVec[nextPIndex].first] );
				 resultSetDist.push_back( pow( spatialDistVec[nextPIndex].second, 0.5 ) );
			 }

			 vector<int> resultid;
			 for(int id: resultIDs){

				 resultid.push_back(id);
			 }

			 string resultprefix = "/mnt/hd4/home/guoqi/coil100/coil100_all/query" + to_string(query_count) + "_diverse.csv";
			 writeResult(resultprefix, resultid);


			double div = evaluation.computeOriDiv( query, result, distanceFunction ) ;
			double rel = evaluation.computeRel( resultSetDist, groundTruthDist );
			double angDiv = evaluation.computeAngDiv( query, result, distanceFunction );

			divVec.push_back( div );
			relVec.push_back( rel );
			divRelVec.push_back( evaluation.computeDivRel( div, rel, alpha ) );
			angDivVec.push_back( angDiv );

			for( int i = 0; i < query.attributes.size(); i++ )
			{
				dataset.dataset[query_count].attributes[i] = query.attributes[i];
			}
			query_count++;

//			cout << "div = "<< div << " , rel = " << rel << " , angDiv = " << angDiv << endl;
		}

		// calculate and print the average
		double avgDiv = 0.0;
		double avgRel = 0.0;
		double avgDivRel = 0.0;
		double avgAngDiv = 0.0;
		int querySize = dataset.queryset.size();

		for ( int i = 0; i < querySize; i++ )
		{
			avgDiv += divVec[i];
			avgRel += relVec[i];
			avgDivRel += divRelVec[i];
		    avgAngDiv += angDivVec[i];
		}


		double running_time = timer.elapsed();


		cout << "running time: " << running_time << endl;
		cout << "average div: " << avgDiv / (querySize * 1.0) << endl;
		cout << "average rel: " << avgRel / (querySize * 1.0) << endl;
		cout << "average divRel: " << avgDivRel / (querySize * 1.0) << endl;
		cout << "average angle div: " << avgAngDiv / (querySize * 1.0) << endl;
}


