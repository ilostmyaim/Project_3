#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include "CUBE.h"

using namespace std;

int main(int argc, char **argv) {
	vector_t vec;

	int d,m,probes;
	string input_file, output_file, query_file, met;
	Metric metric;
	string input;
	/*initialize parameters*/
	initParametersCube(&d, &m, &probes, input_file, output_file, query_file, met, argc, argv);
	if(met.compare("cosine") == 0) { 
		metric = cosine;
	}
	else{ 
		metric = euclidean;
	}

	/*create a Hybercube object*/
	//metric = cosine;
	CUBE cubeObject(d, m, probes, input_file, output_file, query_file, metric);

	ofstream outputFile(output_file);
	streambuf *coutbuf = cout.rdbuf();
	cout.rdbuf(outputFile.rdbuf());

	
	cubeObject.executeCUBE(metric);

	while(true) {

		cout.rdbuf(coutbuf);

		cout << "Do you want to use a different query file?(y/n): ";
		getline(cin, input);
		if(input.compare("y") == 0){
			cout << "Enter name/path of query file: ";
			getline(cin, input);

			ofstream outputFile(output_file);
			streambuf *coutbuf = cout.rdbuf();
			cout.rdbuf(outputFile.rdbuf());

			cubeObject.setQueryFileName(input);
			cubeObject.executeCUBE(metric);
			continue;
		}
		else{
			break;
		}
	}
	
	//cubeObject.displayCUBE();
	
	cout.rdbuf(coutbuf);


}	
