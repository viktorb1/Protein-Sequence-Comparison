#ifndef _BLOSUM_H_
#define _BLOSUM_H_

#include <unordered_map>

class Blosum {
    const static int num_aa = 20;

    // stores the blosum62 matrix
    double blosum[num_aa][num_aa];

    // used to retrieve column/row for amino acid
    unordered_map<char, int> aarowcol;


public:
    Blosum();
    double getValue(char aa1, char aa2);
};

#endif