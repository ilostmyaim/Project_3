#include<iostream>
#include<vector>
#include<string>
#include<cstring>
#include<string.h>
#include<fstream>
#include<sstream>
#include<math.h>
#include<algorithm>
#include<chrono>
#include<map>


using namespace std;

int main(int argc, char **argv) {

	string input_file = "vader_lexicon.csv";
	string line;
	multimap<string, const double> vaderLexiconMap;
	ifstream inFile(input_file);
	string word;
	double score;
	if(inFile.is_open()) {
		while(getline(inFile, line)) {
			//line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
			if(line.empty())
				continue;
			istringstream ss(line);
			ss >> word >> score;
			//cout << score << endl;
			auto delimiterPos = line.find(" ");
			auto newPair = std::pair<string, const double>(word, score);
			vaderLexiconMap.insert(newPair);
		}
	}

	cout<< "10q = " << vaderLexiconMap.find("grrr")->second << "\n";

	return 0;
}