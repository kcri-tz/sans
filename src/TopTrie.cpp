#include "TopTrie.h"

TopTrie::node::node() : weight(0), inverse(0), value(0) {
    children = unordered_map<string, node *>();
}

TopTrie::node::~node() {
    //delete children;
}

TopTrie::TopTrie(unsigned int size) : max(size), num(0) {
    list = map<double, vector<string>, greater<double>>();
    root = new node();
}

TopTrie::~TopTrie() {
    //delete list;
    delete root;
}

void TopTrie::addNodes(vector<string> &objects, unsigned int weight, unsigned int inverse, const SANS_opt& opt) {

    node *subtrie = root;
    for (string &obj : objects) {
        if (subtrie->children.count(obj) == 0) {
            subtrie->children[obj] = new node();
        }
        subtrie = subtrie->children[obj];
    }

    subtrie->weight += weight;
    subtrie->inverse += inverse;
    
    auto new_value = 0.0;
    if (opt.allow_asym) {
        new_value = sqrt(subtrie->weight+1) * sqrt(subtrie->inverse+1);
    } else {
        new_value = sqrt(subtrie->weight) * sqrt(subtrie->inverse);
    }

    if (new_value > 0) {
        if (list.erase(subtrie->value) == 1) {
            subtrie->value = new_value;
            while (!list.emplace(subtrie->value, objects).second) {
                subtrie->value -= 0.000001;
            }
            return;
        }

        if (max > 0 && num == max) {
            double minimum = list.rbegin()->first;
            if (minimum < new_value && list.erase(minimum) == 1) {
                subtrie->value = new_value;
                while (!list.emplace(subtrie->value, objects).second) {
                    subtrie->value -= 0.000001;
                }
            }
            return;
        }

        ++num;
        subtrie->value = new_value;
        while (!list.emplace(subtrie->value, objects).second) {
            subtrie->value -= 0.000001;
        }
    }
}
