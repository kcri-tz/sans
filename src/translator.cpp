#include "translator.h"


string translator::defaultCodon = "";
unordered_map<string, string> translator::codonTable;
unordered_set<char> translator::bases;


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
        } else {
            cerr << "Cannot find file " << codonfile << ", using default translation" << endl;
            readDefault();
            initialized = true;
        }

    } else {
        readDefault();
        initialized = true;
    }

    return initialized;
}

void translator::initAllowedBases() {
    bases.insert('A');
    bases.insert('C');
    bases.insert('G');
    bases.insert('T');
    bases.insert('U');
}

void translator::addTranslationUnit(string &unit) {
    if (unit[0] != '>') {
        size_t pos = 0;
        string searchString;
        string token;
        string delimiter = "=";
        string triplet;
        string amino;

        searchString.assign(unit);
        if ((pos = searchString.find(delimiter)) != string::npos) {
            triplet = searchString.substr(0, pos);
            amino = searchString.substr(pos+1, searchString.length());
            searchString.erase(0, pos + delimiter.length());
        }

        if (translator::checkUnit(triplet)) {
            codonTable[triplet] = amino;
           // cout << "Base " << triplet << " Amino " << codonTable[triplet] << endl;
        }else{
            cerr << "Unknown base in triplet ->" << unit;
        }


    }
}

string translator::getTranslatedAminoAcid(string &unit) {
    string translated;
    // To imitate the biological process. Can be omitted in favour of the runtime.
    std::replace( unit.begin(), unit.end(), 'T', 'U');
    auto iterator = translator::codonTable.find(unit);

    if (iterator == translator::codonTable.end()) {
        cerr << "Could not find translation for " << unit << ", skipping sequence" << endl;
        translated = "";
    } else {
        translated = iterator->second;

        if (translated.find('#') != string::npos) {
            cerr << "Found stop codon for sequence " << unit << ", skipping sequence" << endl;
            translated = "";
        }
    }

    return  translated;
}

void translator::readDefault() {
    //we could read the default data from a remote-stream later
    translator::defaultCodon = "UUU=F;UUC=F;UUA=L;UUG=L;UCU=S;UCC=S;UCA=S;UCG=S;UAU=Y;UAC=Y;UAA=#;UAG=#;UGU=C;UGC=C;UGA=#;UGG=W;CUU=L;CUC=L;CUA=L;CUG=L;CCU=P;CCC=P;CCA=P;CCG=P;CAU=H;CAC=H;CAA=Q;CAG=Q;CGU=R;CGC=R;CGA=R;CGG=R;AUU=I;AUC=I;AUA=I;AUG=M;ACU=T;ACC=T;ACA=T;ACG=T;AAU=N;AAC=N;AAA=K;AAG=K;AGU=S;AGC=S;AGA=R;AGG=R;GUU=V;GUC=V;GUA=V;GUG=V;GCU=A;GCC=A;GCA=A;GCG=A;GAU=D;GAC=D;GAA=D;GAG=D;GGU=G;GGC=G;GGA=G;GGG=G";
    size_t pos = 0;
    string delimiter = ";";
    string token;
    string searchString;
    searchString.assign(translator::defaultCodon);
    while ((pos = searchString.find(delimiter)) != string::npos) {
        token = searchString.substr(0, pos);
        translator::addTranslationUnit(token);
        searchString.erase(0, pos + delimiter.length());
    }
    translator::addTranslationUnit(searchString);
}

bool translator::checkUnit(string &basicString) {
    bool ok = true;

    // unit has to be 3 chars
    util::trim(basicString);
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
    unordered_set<char>::const_iterator got = translator::bases.find (base);
    return got != translator::bases.end();
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
