#ifndef BIFROST_TOPTRIE_H
#define BIFROST_TOPTRIE_H

#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>

using namespace std;

class TopTrie {

protected:
    struct node {
        friend class TopTrie;

    private:
        unsigned int weight;
        unsigned int inverse;
        double value;

        unordered_map<string, node *> children;

    public:
        node();
        ~node();
    };

private:
    unsigned int max;
    unsigned int num;

    node *root;

public:
    TopTrie(unsigned int size);
    ~TopTrie();

    map<double, vector<string>, greater<double>> list;

    void addNodes(vector<string> &objects, unsigned int weight, unsigned int inverse, const SANS_opt& opt);

    void filterStrict();
    void filterWeakly();

    static bool compatible(vector<string>&, vector<vector<string>>&);
};

#endif //BIFROST_TOPTRIE_H
