#ifndef _SEQUENCES_H_
#define _SEQUENCES_H_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <string>
#include "Blosum.h"


class SeqComparison {
    string seq1, seq2;
    string seq1Name, seq2Name;
    double K3Hat, distance;
    Blosum blosum62;

    // kmax stores the min(seq1Length, seq2Length);
    int kmax;
    double beta;

    // read AA sequence files
    void readFile(string& seq, string& name, int filenum);
    void readToString(ifstream& inf, string& seq, string& name);

    // partial sum of equation 5
    double calculateK3_1(vector<double>& product, string& s1, string& s2);

    // partial sum of equation 5
    double calculateK3_k(vector<double>& product, int k, string& s1, string& s2);

    //equation 6
    double calculateK3(string& s1, string& s2, int length);


public:
    SeqComparison();
    double calculateK3Hat();
    double calculateDistance();
    void printResults();
};

#endif