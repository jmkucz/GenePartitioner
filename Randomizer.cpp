#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <cstring>
#include <string.h>
#include <algorithm>
#include <vector>
using namespace std;

int main() {
    string name;
    string gene;
    vector<string> names = { "Extra_Back.txt", "Extra_Front.txt", "n3B30.txt",
    "n3F30.txt", "n4F30.txt", "n4B30.txt", "nnF30.txt", "nnB30.txt" };
    
    for (int i = 0; i < 8; i++) {
        string inFile = "Rand_9_ERBB2_3U_" + names[i];
        string outFile = "Rand_10_ERBB2_3U_" + names[i];
        ifstream ifs;
        ifs.open(inFile);
        if (!ifs) {
            cout << "didn't work" << endl;
        }

        ofstream ofs;
        ofs.open(outFile);
        if (!ifs) {
            cout << "didn't work" << endl;
        }

        while (getline(ifs, name)) {
            getline(ifs, gene);
            random_shuffle(gene.begin(), gene.end());
            ofs << name << endl;
            ofs << gene << endl;
        }
        ifs.close();
        ofs.close();
    }
    return 0;
}

