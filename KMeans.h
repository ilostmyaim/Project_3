#pragma once


#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <list>
#include <map>
#include "hash.h"
#include "LSH.h"
#include "Cluster.h"
#include "CUBE.h"

typedef struct { 
	int K; //number of clusters
	int number_of_hash_functions; //default: L=4
	int number_of_hash_tables; //default: L=5
	int max_iterations;
	std::string input_file;
	std::string output_file;
	std::string conf_file;
	std::string met; //metric

	int init_choice, assign_choice, update_choice,probes,MC;
	bool complete=false;
}init_params_t;



class KMeans {
private: 
	int _K; // number of cluster(subsets) for K-Means
	int _max_iter;
	long int _totalItems; //number of items read from dataset
	int _num_hash_functions;
	int _num_hashtables;
	int _max_iterations;
	std::vector<item_t> _items; //points read from input file
	std::vector<Cluster> _clusters;
	std::string _inputFile;
	std::string _outputFile;
	std::string _confFile;
	std::string _met; //metric
	double _avgSilhouette;
	vector_t _bi;
	vector_t _ai;
	vector_t _si; //average si for each cluster
	std::chrono::duration<double> _duration;// execution time of kmeans

	LSH *_LSHObject;
	CUBE *_CUBEObject;

	int _initChoice;
	int _assignChoice;
	int _updateChoice;
	int _MC;
	int _probes;
	bool _complete;

public:
	KMeans(init_params_t init_params);
	void randomInitialization();
	void printClusters();
	bool executeKMeans();
	bool lloydsAssignment();
	bool LSHAssignment();
	bool CUBEAssignment();
	double initialRangeLSH(Metric metric);
	int getNearestCluster(item_t item,Metric metric);
	bool updateMeans();
	void computeMean(const std::multimap<int, const item_t*> &multimap, int clusterID, Cluster *cluster);
	bool improvedLloydsUpdate();
	double computeSilhouette();
};

void initParametersKMeans(init_params_t *init_params, int argc, char **argv);