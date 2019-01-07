#pragma once

#include<string>
#include<vector>


class Hash;
typedef struct item item_t;
typedef std::vector<double> vector_t;


enum Metric { euclidean , cosine } ;


#define M 60000
const int w = 4;
class LSH {

private:

	int _k; //number of locality-sensitive functions hi
	int _L; //number of hash tables
	std::string _inputFile;
	std::string _outputFile;
	std::string _queryFile;
	Hash **_arrayOfHashTables; //L-Hash objects
	int _hashTableSize;



public:

	LSH(int k=4, int L=5, std::string inputFile = "input_file", std::string outputFile="output_file", std::string queryFile="query_file", Metric metric = euclidean); //constructor
	~LSH();
	int get_k();
	int get_L();
	std::string get_inputFile();
	std::string get_outputFile();
	void setQueryFileName(std::string);
	void executeLSH(Metric);
	void insertAllItems(std::vector<item_t> &items,Metric metric);
	void displayLSH();
	std::vector< std::vector<item_t> > rangeSearch(vector_t q, int cluster_id, double R, double C, Metric metric);
	double nearestNeighbor(vector_t, Metric metric);
	int sizeofLSH();

};

/*parse command line parameters*/
void initParameters(int* k, int* L, std::string & input_file, std::string & output_file, std::string & query_file,std::string &met,int argc, char** argv);





