#include <fstream>
#include <algorithm>
#include <string>
#include <unordered_map>
#include "translator.h"

#ifndef default_translate_codons
#define READ_DEFAULT() readDefault()
#define default_translate_codons READ_DEFAULT()
#endif


string translator::defaultCodon = "";
unordered_map<string, string> translator::codonTable;
vector<char> translator::bases;


bool translator::init(string &codonfile) {
    bool initialized = false;
    initAllowedBases();
    if (!codonfile.empty()) {
        ifstream file(codonfile);
        if (file.good()) {
            string line;
            while (getline(file, line)) {
                translator::addTranslationUnit(line);
            }
            file.close();
            initialized = true;
        }
    } else {
        size_t pos = 0;
        string delimiter = ";";
        string token;
        while ((pos = translator::defaultCodon.find(delimiter)) != string::npos) {
            token = translator::defaultCodon.substr(0, pos);
            translator::addTranslationUnit(token);
            translator::defaultCodon.erase(0, pos + delimiter.length());
        }
        initialized = true;
    }

    return initialized;
}

void translator::initAllowedBases() {
    bases.push_back('A');
    bases.push_back('C');
    bases.push_back('G');
    bases.push_back('T');
    bases.push_back('U');
}

void translator::addTranslationUnit(string &unit) {
    if (unit[0] != '>') {
        size_t pos = 0;
        string token;
        string delimiter = "=";
        string triplet;
        string amino;

        while ((pos = unit.find(delimiter)) != string::npos) {
            triplet = unit.substr(0, pos);
            amino = unit.substr(pos+1, unit.length());
            unit.erase(0, pos + delimiter.length());
        }

        if (translator::checkUnit(triplet)) {
            codonTable[triplet] = amino;
            cout << "Base " << triplet << " Amino " << codonTable[triplet] << endl;
        }else{
            cerr << "Unknown base in triplet ->" << unit;
        }


    }
}

string translator::getTranslatedAminoAcid(string &unit) {
    // To imitate the biological process. Can be omitted in favour of the runtime.
    std::replace( unit.begin(), unit.end(), 'T', 'U');
    return translator::codonTable.find(unit)->second;
}

void translator::readDefault() {
    initAllowedBases();
    translator::defaultCodon = "D";
}

bool translator::checkUnit(string &basicString) {
    bool ok = true;

    // unit has to be 3 chars
    ok = basicString.length() == 3;
    if(ok){
        // only valid DNA bases are allowed
        ok = translator::isBase(basicString[0]);
        ok = ok && translator::isBase(basicString[1]);
        ok = ok && translator::isBase(basicString[2]);
    }

    return ok;
}

/**
 * This function checks if the character at the given position is allowed.
 * @param pos position in str
 * @param str the current part of the sequence
 * @return true if allowed, false otherwise
 */
bool translator::isBase(char &base) {
    bool allowed = false;

    for(int i = 0; i < translator::bases.size() && !allowed; i++){
        allowed =  translator::bases.at(i) == base;
    }

    return allowed;
}

string translator::translate(string &line) {
    string translated;
    //cout << "Zeile ist " << line.length() << " lang";


    for(auto pos = 0; pos+3 < line.length(); pos = pos +3){
        string unit = line.substr(pos, 3);
        translated+= translator::getTranslatedAminoAcid(unit);
    }

    return translated;
}
