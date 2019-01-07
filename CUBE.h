#pragma once

#include<string>
#include<vector>
#include "LSH.h"
enum Metric;

class Hash;
typedef struct item item_t;
typedef std::vector<double> vector_t;

class CUBE {

private:
	int _k; //the dimension in which the points are projected
	int _MC; //max number of points to be checked
	int _probes; //max number of vertices to be checked
	std::string _inputFile;
	std::string _outputFile;
	std::string _queryFile;
	Hash *_hashTable;
	int _hashTableSize;

public:
	CUBE(int k=4, int MC=10, int probes=2, std::string inputFile="input_file", std::string outputFile="output_file", std::string queryFile="query_file",Metric = euclidean);
	void setQueryFileName(std::string newName);
	int get_k();
	int get_M();
	int get_probes();
	std::string get_inputFile();
	std::string get_outputFile();
	void displayCUBE();
	void insertAllItems(std::vector<item_t> &items,Metric metric);
	std::vector< std::vector<item_t> >rangeSearch(vector_t q,int cluster_id ,double R, double C, Metric metric);
	double nearestNeighbor(vector_t, Metric metric);

};


void initParametersCube(int* k, int* m, int* probes, std::string &input_file, std::string & output_file, std::string & query_file,std::string &met,int argc, char** argv);

