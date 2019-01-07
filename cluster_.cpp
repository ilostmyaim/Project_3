#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <chrono>
#include "LSH.h"
#include "CUBE.h"
#include "hash.h"
#include "Cluster.h"
#include "KMeans.h"


using namespace std;

int main(int argc, char **argv) {
	vector_t vec;
	string line;

	init_params_t init_params;
	Metric metric;
	string input;
	//initialize parameters
	initParametersKMeans(&init_params, argc, argv);

	ofstream outputFile(init_params.output_file);
	streambuf *coutbuf = cout.rdbuf();
	cout.rdbuf(outputFile.rdbuf());

	KMeans KMeansObject(init_params);
	if(init_params.init_choice == 1) { 
		KMeansObject.randomInitialization();
	}
	else{
		//kmeans++ initialization
	}
	
	KMeansObject.executeKMeans();
	KMeansObject.computeSilhouette();
	KMeansObject.printClusters();

	cout.rdbuf(coutbuf);

}