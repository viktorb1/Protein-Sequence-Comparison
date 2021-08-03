using namespace std;

#include "Blosum.h"
#include <fstream>
#include <iostream>
#include <unordered_map>

Blosum::Blosum()
{
    char aa;
    string temp;
    ifstream inf("BLOSUM62-2.txt");

    // skip first row in blosum62 file
    getline(inf, temp);
    
    for(int i = 0; i < num_aa; i++)
    {
        // read in amino acid order from first column
        inf >> aa;
        aarowcol[aa] = i;

        // add blosum62 data to 2D array
        for(int j = 0; j < num_aa; j++)
            inf >> blosum[i][j];
    }

    inf.close();
}


// returns Blosum62 data given two amino acids
double Blosum::getValue(char aa1, char aa2)
{
    return blosum[ aarowcol[aa1] ][ aarowcol[aa2] ];

}