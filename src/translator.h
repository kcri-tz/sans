#ifndef SANS_TRANSLATOR_H
#define SANS_TRANSLATOR_H


#include <string>
#include <iostream>
#include <vector>

using namespace std;


class translator {

private:
    static unordered_map<string, string> codonTable;
    static string defaultCodon;
    static void readDefault();
    static vector<char> bases;

    static bool checkUnit(string &basicString);
    static bool isBase(char &base);
public:
    static bool init(string& codonfile);
    static void addTranslationUnit(string& unit);
    static string getTranslatedAminoAcid(string& unit);
    static string translate(string& line);

protected:
    static void initAllowedBases();
};


#endif //SANS_TRANSLATOR_H
