#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <cstring>
#include <sstream>
#include <map>
#include <string.h>
#include <vector>
#include <algorithm>
#include <iterator>
using namespace std;

string newFile(const char *file, int length);
void copyString(const char *tempfile, char *file, int length);
int size(const char *tempfile);
int numOfSpaces(const char *tempfile);
void addFrontNucleotidesToFile(string target, int range, const char *file, int fileLength, string nameOfFile);
void addBackNucleotidesToFile(string target, int range, const char *file, int fileLength, string nameOfFile);
void addBackFirstThirtyToFile(string target, int range, const char *file, int fileLength, string nameOfFile);
void addFrontFirstThirtyToFile(string target, int range, const char *file, int fileLength, string nameOfFile);
void addFrontSecondThirtyToFile(string target, int range, const char *file, int fileLength, string nameOfFile);
void addFrontThirdThirtyToFile(string target, int range, const char *file, int fileLength, string nameOfFile);
void addBackSecondThirtyToFile(string target, int range, const char *file, int fileLength, string nameOfFile);
void addBackThirdThirtyToFile(string target, int range, const char *file, int fileLength, string nameOfFile);
int firstLocationOnGene(const char *file, int fileLength, string target, int start);
int lastLocationOnGene(const char *file, int fileLength, string target, int start);
string frontThirty(int positionOfTarget, int range, const char *file, int fileLength, string target);
string backThirty(int positionOfTarget, int range, const char *file, int fileLength, string target);
string frontNucleotides(int positionOfStart, int positionOfPreviousTarget, const char *file, int fileLength, int rangePlusTarget, int range);
string backNucleotides(int positionOfStart, char *positionOfNextTarget, const char *file, int fileLength);
void thirtyOfString(string group, int place, ofstream &ofs, string nameOfFile, int start, int fileLength);
void thirtyOfStringFront(string group, int place, ofstream &ofs, string nameOfFile, int start, int fileLength);

//Looks through each nucleotide to find the squence of the target gene. 
int firstLocationOnGene(const char *file, int fileLength, string target, int start) {
    int i = start;
    for (const char *ptr = file + start; ptr < file + fileLength; ptr++) {
        if (target[0] == *ptr && target[1] == *(ptr + 1) && ptr < file + fileLength - size(target)) {
            int z = 0;
            for (int j = 0; j < size(target); j++) {
                
                if (target[j] != *(ptr + j)) {
                    j = size(target);
                }
                else {
                    z += 1;
                    if (z == size(target)) {
                        return i;
                    }
                }
            }
        }
        i++;
    }
    return -1;
}

int lastLocationOnGene(const char *file, int fileLength, string target, int start) {
    if (firstLocationOnGene(file, fileLength, target, start) != -1) {
        return firstLocationOnGene(file, fileLength, target, start) + size(target) - 1;
    }
    return -1;
}

string frontThirty(int positionOfTarget, int range, const char *file, int fileLength, string target) {
    string firstThirty = "";
    if (positionOfTarget <= -1) {
        return "";
    }
    if (positionOfTarget > range) {
        for (int i = range; i >= 1; i--) {
            firstThirty += *(file + positionOfTarget - i);
        }
    }
    else {
        for (int j = positionOfTarget; j >= 1; j--) {
            firstThirty += *(file + positionOfTarget - j);
        }
    }
    return firstThirty;
}

string backThirty(int positionOfEndTarget, int range, const char *file, int fileLength, string target) {
    string lastThirty = "";
    if (positionOfEndTarget == -1) {
        return "";
    }
    if (fileLength - positionOfEndTarget < range) {
        for (int i = 1; i < fileLength - positionOfEndTarget; i++) {
            lastThirty = *(file + positionOfEndTarget + i);
        }
    }
    else {
        for (int i = 1; i < range + 1; i++) {
            lastThirty += *(file + positionOfEndTarget + i);
        }
    }
    return lastThirty;
}


void addFrontFirstThirtyToFile(string target, int range, const char *file, int fileLength, string nameOfFile) {
    int i = 0;
    ofstream ofs;
    ofs.open("NM_005228_3U_nnF30.txt", ofs.out | ofs.app);
    if (!ofs.is_open()) {
        cout << "not writing" << endl;
    }

    if (frontThirty(firstLocationOnGene(file, fileLength, target, i), range, file, fileLength, target) != "") {
        
        for (const char *ptr = file; ptr < file + fileLength; ptr++) {
            if (frontThirty(firstLocationOnGene(file, fileLength, target, i), range, file, fileLength, target) != "") {
                int length = range - 1;
                int end = firstLocationOnGene(file, fileLength, target, i);
                int start = end - length;
                if (end < range) {
                    start = 1;
                }
                ofs << ">" << nameOfFile << ":nnF30=" << start << "-" << end << '\n';
                //cout << frontThirty(firstLocationOnGene(file, fileLength, target, i), range, file, fileLength, target) << '\n';
                ofs << frontThirty(firstLocationOnGene(file, fileLength, target, i), range, file, fileLength, target) << '\n';
                i = firstLocationOnGene(file, fileLength, target, i);
            }
            i++;
        }
    }
    ofs.close();

}

void addFrontSecondThirtyToFile(string target, int range, const char *file, int fileLength, string nameOfFile) {
    int i = 0;
    ofstream ofs;
    ofs.open("NM_005228_3U_n3F30.txt", ofs.out | ofs.app);
    if (!ofs.is_open()) {
        cout << "not writing" << endl;
    }

    if (frontThirty(firstLocationOnGene(file, fileLength, target, i) - range, range, file, fileLength, target) != "") {
        
        for (const char *ptr = file; ptr < file + fileLength; ptr++) {
            if (frontThirty(firstLocationOnGene(file, fileLength, target, i) - range, range, file, fileLength, target) != "") {
                if (i != -1) {
                    int length = range - 1;
                    int end = firstLocationOnGene(file, fileLength, target, i) - range;
                    int start = end - length;
                    if (end < range) {
                        start = 1;
                    }
                    ofs << ">" << nameOfFile << ":n3F30=" << start << "-" << end << '\n';
                    //cout << frontThirty(firstLocationOnGene(file, fileLength, target, i) - range, range, file, fileLength, target) << '\n';
                    ofs << frontThirty(firstLocationOnGene(file, fileLength, target, i) - range, range, file, fileLength, target) << '\n';
                    i = firstLocationOnGene(file, fileLength, target, i);
                }
                if (i == -1) {
                    return;
                }
            }
            i++;
        }
    }
    ofs.close();

}

void addFrontThirdThirtyToFile(string target, int range, const char *file, int fileLength, string nameOfFile) {
    int i = 0;
    ofstream ofs;
    ofs.open("NM_005228_3U_n4F30.txt", ofs.out | ofs.app);
    if (!ofs.is_open()) {
        cout << "not writing" << endl;
    }

    if (frontThirty(firstLocationOnGene(file, fileLength, target, i) - range * 2, range, file, fileLength, target) != "") {
        
        for (const char *ptr = file; ptr < file + fileLength; ptr++) {
            if (frontThirty(firstLocationOnGene(file, fileLength, target, i) - range * 2, range, file, fileLength, target) != "") {
                if (i != -1) {
                    int length = range - 1;
                    int end = firstLocationOnGene(file, fileLength, target, i) - range * 2;
                    int start = end - length;
                    if (end < range) {
                        start = 1;
                    }
                    ofs << ">" << nameOfFile << ":n4F30=" << start << "-" << end << '\n';
                    //cout << frontThirty(firstLocationOnGene(file, fileLength, target, i) - range * 2, range, file, fileLength, target) << '\n';
                    ofs << frontThirty(firstLocationOnGene(file, fileLength, target, i) - range * 2, range, file, fileLength, target) << '\n';
                    i = firstLocationOnGene(file, fileLength, target, i);
                }
                if (i == -1) {
                    return;
                }
            }
            i++;
        }
    }
    ofs.close();

}

void addBackFirstThirtyToFile(string target, int range, const char *file, int fileLength, string nameOfFile) {
    int i = 0;
    ofstream ofs;
    ofs.open("NM_005228_3U_nnB30.txt", ofs.out | ofs.app);
    if (!ofs.is_open()) {
        cout << "not writing" << endl;
    }

    if (backThirty(lastLocationOnGene(file, fileLength, target, i), range, file, fileLength, target) != "") {
        
        for (const char *ptr = file; ptr < file + fileLength; ptr++) {
            if (backThirty(lastLocationOnGene(file, fileLength, target, i), range, file, fileLength, target) != "") {
                int length = range - 1;
                int start = lastLocationOnGene(file, fileLength, target, i) + 2;
                if (fileLength - start < range) {
                    length = fileLength - start + 1;
                }
                int end = start + length;
                ofs << ">" << nameOfFile << ":nnB30=" << start << "-" << end << '\n';
                //cout << backThirty(lastLocationOnGene(file, fileLength, target, i), range, file, fileLength, target) << '\n';
                ofs << backThirty(lastLocationOnGene(file, fileLength, target, i), range, file, fileLength, target) << '\n';
                i = lastLocationOnGene(file, fileLength, target, i);
            }
            i++;
        }
    }
    ofs.close();

}

void addBackSecondThirtyToFile(string target, int range, const char *file, int fileLength, string nameOfFile) {
    int i = 0;
    ofstream ofs;
    ofs.open("NM_005228_3U_n3B30.txt", ofs.out | ofs.app);
    if (!ofs.is_open()) {
        cout << "not writing" << endl;
    }

    if (backThirty(lastLocationOnGene(file, fileLength, target, i) + range, range, file, fileLength, target) != "" && lastLocationOnGene(file, fileLength, target, i) != -1) {
        
        for (const char *ptr = file; ptr < file + fileLength; ptr++) {
            if (backThirty(lastLocationOnGene(file, fileLength, target, i) + range, range, file, fileLength, target) != "") {
                if (lastLocationOnGene(file, fileLength, target, i) != -1) {
                    int length = range - 1;
                    int start = lastLocationOnGene(file, fileLength, target, i) + range + 2;
                    if (fileLength - start < range) {
                        length = fileLength - start + 1;
                    }
                    int end = start + length;
                    ofs << ">" << nameOfFile << ":n3B30=" << start << "-" << end << '\n';
                    //cout << backThirty(lastLocationOnGene(file, fileLength, target, i) + range, range, file, fileLength, target) << '\n';
                    ofs << backThirty(lastLocationOnGene(file, fileLength, target, i) + range, range, file, fileLength, target) << '\n';
                    i = lastLocationOnGene(file, fileLength, target, i);
                }
                if (i == -1) {
                    return;
                }
            }
            i++;
        }
    }
    ofs.close();

}

void addBackThirdThirtyToFile(string target, int range, const char *file, int fileLength, string nameOfFile) {
    int i = 0;
    int positionOfPreviousTarget = 0;
    ofstream ofs;
    ofs.open("NM_005228_3U_n4B30.txt", ofs.out | ofs.app);
    if (!ofs.is_open()) {
        cout << "not writing" << endl;
    }
    if (backThirty(lastLocationOnGene(file, fileLength, target, i) + range * 2, range, file, fileLength, target) != "" && lastLocationOnGene(file, fileLength, target, i) != -1) {
        
        for (const char *ptr = file; ptr < file + fileLength; ptr++) {
            if (backThirty(lastLocationOnGene(file, fileLength, target, i) + range * 2, range, file, fileLength, target) != "") {
                if (lastLocationOnGene(file, fileLength, target, i) != -1) {
                    int length = range - 1;
                    int start = lastLocationOnGene(file, fileLength, target, i) + range * 2 + 2;
                    if (fileLength - start < range) {
                        length = fileLength - start + 1;
                    }
                    int end = start + length;
                    ofs << ">" << nameOfFile << ":n4B30=" << start << "-" << end << '\n';
                    //cout << backThirty(lastLocationOnGene(file, fileLength, target, i) + range * 2, range, file, fileLength, target) << '\n';
                    ofs << backThirty(lastLocationOnGene(file, fileLength, target, i) + range * 2, range, file, fileLength, target) << '\n';
                    i = lastLocationOnGene(file, fileLength, target, i);
                }
                if (i == -1) {
                    return;
                }
            }
            i++;
        }
    }
    ofs.close();

}

string frontNucleotides(int positionOfStart, int positionOfPreviousTarget, const char *file, int fileLength, int rangePlusTarget, int range) {
    string frontNucleotides = "";
    if (positionOfStart == -1 || positionOfPreviousTarget > positionOfStart || (positionOfStart + range*3) <= rangePlusTarget) {
        return "";
    }

    for (const char *ptr = file + positionOfPreviousTarget; ptr < file + positionOfStart - 1; ptr++) {
        frontNucleotides += *ptr;
    }
    return frontNucleotides;
}


void addFrontNucleotidesToFile(string target, int range, const char *file, int fileLength, string nameOfFile) {
    int positionOfPreviousTarget = 0;
    int frontPositionOfPreviousTargetArr[10000] = {};
    int endPositionOfTwoTargetsAgoArr[10000] = {};
    endPositionOfTwoTargetsAgoArr[0] = 0;
    int i = 0;
    ofstream ofs;
    int rangePlusTarget = range * 3 + size(target);
    ofs.open("NM_005228_3U_Extra_Front.txt", ofs.out | ofs.app);
    if (!ofs.is_open()) {
        cout << "not writing" << endl;
    }
    string one = frontNucleotides(firstLocationOnGene(file, fileLength, target, i) - range * 3 + 1, positionOfPreviousTarget, file, fileLength, rangePlusTarget, range);
    string two = frontNucleotides(firstLocationOnGene(file, fileLength, target, i) + 1, positionOfPreviousTarget, file, fileLength, rangePlusTarget, range);
    if (one != "" || (two != file && two != "")) {
        if (frontNucleotides(firstLocationOnGene(file, fileLength, target, i) - range * 3 + 1, positionOfPreviousTarget, file, fileLength, rangePlusTarget, range) == "" && frontNucleotides(firstLocationOnGene(file, fileLength, target, i) + 1, positionOfPreviousTarget, file, fileLength, rangePlusTarget, range) != "") {
            positionOfPreviousTarget = firstLocationOnGene(file, fileLength, target, i) + size(target) + range * 3;
        }
        
        
        for (const char *ptr = file; ptr < file + fileLength; ptr++) {
            if (frontNucleotides(firstLocationOnGene(file, fileLength, target, i) - range * 3 + 1, positionOfPreviousTarget, file, fileLength, rangePlusTarget, range) != "") {
                
                string group = frontNucleotides(firstLocationOnGene(file, fileLength, target, i) - range * 3 + 1, positionOfPreviousTarget, file, fileLength, rangePlusTarget, range);
                //cout << frontNucleotides(firstLocationOnGene(file, fileLength, target, i) - range * 3 + 1, positionOfPreviousTarget, file, fileLength, rangePlusTarget, range) << endl;
                int x = positionOfPreviousTarget;
                int f = 0;
                reverse(begin(group), end(group));
                for (int j = 0; j < ceil(size(group) / 30) + 1; j++) {
                    thirtyOfStringFront(group, j, ofs, nameOfFile, positionOfPreviousTarget + size(group) + f, fileLength);
                    f -= 30;
                }
                i = firstLocationOnGene(file, fileLength, target, i);
                if (i > 0) {
                    positionOfPreviousTarget = i + range * 3 + size(target);

                }
            }
            i++;
        }
    }
    ofs.close();
}

void thirtyOfStringFront(string group, int place, ofstream &ofs, string nameOfFile, int start, int fileLength) {
    place = place * 30;
    string thirty = "";
    int range = 30;
    if (size(group) - place < 30) {
        range = size(group) - place;
    }

    int ending = start;
    start = ending - range + 1;
    if (ending == fileLength - 1) {
        ending++;
    }
    for (int i = place; i < place + range; i++) {
        thirty += group[i];
    }
    reverse(begin(thirty), end(thirty));
    ofs << ">" << nameOfFile << ":Front_Extra=" << start << "-" << ending << '\n';
    //cout << thirty << endl;
    ofs << thirty << endl;
}

void thirtyOfString(string group, int place, ofstream &ofs, string nameOfFile, int start, int fileLength) {
    place = place * 30;
    string thirty = "";
    int range = 30;
    if (size(group) - place < 30) {
        range = size(group) - place;
    }
    
    int end = start + range - 1;
    if (end == fileLength + 1) {
        end--;
    }
    for (int i = place; i < place + range; i++) {
        thirty += group[i];
    }
    ofs << ">" << nameOfFile << ":Back_Extra=" << start << "-" << end << '\n';
    //cout << thirty << endl;
    ofs << thirty << endl;
}

string backNucleotides(int positionOfStart, const char *positionOfNextTarget, const char *file, int fileLength) {
    string backNucleotides = "";
    if (positionOfStart == -1) {
        return "";
    }
    for (const char *ptr = file + positionOfStart; ptr <= positionOfNextTarget; ptr++) {
        backNucleotides += *ptr;
    }
    return backNucleotides;
}

void addBackNucleotidesToFile(string target, int range, const char *file, int fileLength, string nameOfFile) {
    const char *positionOfNextTarget = file + fileLength - 1;
    int i = 0;
    int frontPositionOfTwoTargetsFromNowArr[10000] = {};
    int endPositionOfNextTargetArr[10000] = {};
    frontPositionOfTwoTargetsFromNowArr[0] = fileLength - 1;
    ofstream ofs;
    ofs.open("NM_005228_3U_Extra_Back.txt", ofs.out | ofs.app);
    if (!ofs.is_open()) {
        cout << "not writing" << endl;
    }
    if (backNucleotides(lastLocationOnGene(file, fileLength, target, i) + range * 3 + 1, positionOfNextTarget, file, fileLength) != "" && lastLocationOnGene(file, fileLength, target, i) != -1) {
       
        for (const char *ptr = file; ptr < file + fileLength; ptr++) {
            if (backNucleotides(lastLocationOnGene(file, fileLength, target, i), positionOfNextTarget, file, fileLength) != "") {
                if (firstLocationOnGene(file, fileLength, target, lastLocationOnGene(file, fileLength, target, i)) != -1) {
                    positionOfNextTarget = file + firstLocationOnGene(file, fileLength, target, lastLocationOnGene(file, fileLength, target, i + range)) - range * 3;
                }
                string group = backNucleotides(lastLocationOnGene(file, fileLength, target, i) + range * 3 + 1, positionOfNextTarget, file, fileLength);
                //cout << backNucleotides(lastLocationOnGene(file, fileLength, target, i) + range * 3 + 1, positionOfNextTarget, file, fileLength) << endl;

                int x = lastLocationOnGene(file, fileLength, target, i) + range * 3 + 1;

                int f = 1;
                for (int j = 0; j < ceil(size(group) / 30) + 1; j++) {
                    thirtyOfString(group, j, ofs, nameOfFile, x + f, fileLength);
                    f += 30;
                }
                i = lastLocationOnGene(file, fileLength, target, i);
                if (i > 0) {
                    if (lastLocationOnGene(file, fileLength, target, lastLocationOnGene(file, fileLength, target, i + 1)) == -1) {
                        positionOfNextTarget = file + fileLength;
                    }
                    else {
                        positionOfNextTarget = file + i + range * 3 + size(target);
                    }

                }
            }
            i++;
        }
    }
    ofs.close();
}

void copyString(const char *tempfile, char *file, int length) {
    char *filePtr = file;
    for (const char *ptr = tempfile; ptr < tempfile + length; ptr++) {
        *filePtr = *ptr;
        filePtr++;
    }
}

int size(const char *tempfile) {
    int num = 0;
    const char *ptr = tempfile;
    while (*ptr != '\0') {
        num++;
        ptr++;
    }
    return num;
}

int numOfSpaces(const char *tempfile) {
    int num = 0;
    const char *ptr = tempfile;
    while (*ptr != '\0') {
        if (*ptr == '\n') {
            num++;
        }
        ptr++;
    }
    return num;
}

string newFile(const char *file, int length) {
    string newFile = "";
    for (const char *ptr = file; ptr < file + length; ptr++) {
        if (*ptr != '\n') {
            newFile += *ptr;
        }
    }
    return newFile;
}

int main() {
    string input = "";
    int spot = 0;
    map<string, int> names;
    while (getline(cin, input)) {
        istringstream stream(input);
        vector<string> vec;
        copy(istream_iterator<string>(stream), istream_iterator<string>(), back_inserter(vec));
        vec[spot].erase(0, 4);
        
        string nameOfFile = "NM_005228:" + vec[spot];
        auto it = names.find(nameOfFile);
        
        spot++;
        string target = vec[spot];
        spot = 0;
        int range = 30;
        string tempfile = "";
        ifstream ifs;
        ifs.open("EGFR_NM_005228_3U.txt");
        if (!ifs) {
            cout << "didn't work" << endl;
        }
        string sFile((istreambuf_iterator<char>(ifs)),
            istreambuf_iterator<char>());
        const char *oldFile = sFile.c_str();
        int sizeFile = size(oldFile);
        int numSpaces = numOfSpaces(oldFile);
        int sizeNew = sizeFile - numSpaces;
        tempfile = newFile(oldFile, sizeFile);
        const char *file = tempfile.c_str();
        int fileLength = size(file);
        if (it == names.end()) {
            if (firstLocationOnGene(file, fileLength, target, 0) != -1) {
                names.insert({ nameOfFile, 1 });
            }
        }
        else {
            nameOfFile += "-";
            nameOfFile += (char)(it->second + 48);
            it->second++;
        }
        //cout << "front first" << endl;
        addFrontFirstThirtyToFile(target, range, file, fileLength, nameOfFile);
        //cout << endl;
        //cout << "front second" << endl;
        addFrontSecondThirtyToFile(target, range, file, fileLength, nameOfFile);
        //cout << endl;
        //cout << "front third" << endl;
        addFrontThirdThirtyToFile(target, range, file, fileLength, nameOfFile);
        //cout << "Back first" << endl;
        addBackFirstThirtyToFile(target, range, file, fileLength, nameOfFile);
        //cout << endl;
        //cout << "back second" << endl;
        addBackSecondThirtyToFile(target, range, file, fileLength, nameOfFile);
        //cout << endl;
        //cout << "Back third" << endl;
        addBackThirdThirtyToFile(target, range, file, fileLength, nameOfFile);
        //cout << endl;
        //cout << "front extra" << endl;
        addFrontNucleotidesToFile(target, range, file, fileLength, nameOfFile);
        //cout << endl;
        //cout << "back extra" << endl;
        addBackNucleotidesToFile(target, range, file, fileLength, nameOfFile);
        ifs.close();
    }    
    return 0;
}