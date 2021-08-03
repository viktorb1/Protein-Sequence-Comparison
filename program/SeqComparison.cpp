using namespace std;

#include "SeqComparison.h"
#include <math.h>
#include <iomanip>


SeqComparison::SeqComparison()
{
    readFile(seq1, seq1Name, 1);
    readFile(seq2, seq2Name, 2);
    kmax = min(seq1.length(), seq2.length());

    cout << "3) Enter the value of Beta: ";
    cin >> beta;
    calculateDistance();
}


void SeqComparison::readFile(string& seq, string& name, int filenum)
{
    ifstream inf;
    string f;

    cout << filenum << ") Enter the filename of sequence " << filenum << ": ";
    cin >> f;
    inf.open(f.c_str());

    // if the file could not be read properly
    while (!inf)
    {
        cout << "   File not found. Please try again: ";
        cin.clear();
        cin >> f;
        inf.open(f.c_str());
    }

    readToString(inf, seq, name);
}


void SeqComparison::readToString(ifstream& inf, string& seq, string& name)
{
    string temp;
    char c;

    // if there is a description line in the file
    if (inf.peek() == '>')
    {
        inf >> name;
        // remove '>' from name of amino acid discription
        name.erase(0, 1);

        // skip the rest of the first line
        getline(inf, temp);
    }

    while (inf.get(c))
        if(c != '\n')
            seq += c;
}


// calculate K3 for each single AA
double SeqComparison::calculateK3_1(vector<double>& product, string& s1, string& s2)
{
    double sum = 0;

    for (size_t i = 0; i < s1.length(); i++)
        for (size_t j = 0; j < s2.length(); j++)
        {
            double value = pow(blosum62.getValue(s1[i], s2[j]), beta);
            product.push_back(value);
            sum += value;
        }

    return sum;
}


// calculates K3 for each substring where the length of the substring is 
// 1 < |AA| <= kmax
//
// this function uses the previous multiplication calculations to
// recursively build up the solution
double SeqComparison::calculateK3_k(vector<double>& product, int k, string& s1, string& s2)
{
    double sum = 0;

    for (size_t i = 0; i < s1.length() - k + 1; i++)
        for (size_t j = 0; j < s2.length() - k + 1; j++)
        {
            double value = pow(blosum62.getValue(s1[i+k-1], s2[j+k-1]), beta);
            product[i*s2.length()+j] *= value;
            sum += product[i*s2.length()+j];
        }

    return sum;
}


// returns the un-normalized string kernel K3 (equation 5)
double SeqComparison::calculateK3(string& s1, string& s2, int length)
{
    vector<double> product;
    double totalsum = calculateK3_1(product, s1, s2);

    for (int k = 2; k <= length; k++)
        totalsum += calculateK3_k(product, k, s1, s2);

    return totalsum;
}


// returns the normalized string kernel (equation 6)
double SeqComparison::calculateK3Hat()
{
    double K3_ST = calculateK3(seq1, seq2, kmax),
           K3_SS = calculateK3(seq1, seq1, seq1.length()),
           K3_TT = calculateK3(seq2, seq2, seq2.length());

    return K3Hat = K3_ST / pow(K3_SS*K3_TT, 0.5);
}


// returns the distance between two AA sequences (equation 7)
double SeqComparison::calculateDistance()
{
    calculateK3Hat();
    return distance = pow(2*(1-K3Hat), 0.5);
}


void SeqComparison::printResults()
{
    cout << "   ---------------------------------------------------------\n";
    cout << "   |  Identifier                    |  Number of Residues  |\n";
    cout << "   |--------------------------------|----------------------|\n";
    cout << "   |  " << setw(28) << left <<  seq1Name << "  |         " 
         << setw(13) << seq1.length() << "|\n";
    cout << "   |  " << setw(28) << left <<  seq2Name << "  |         " 
         << setw(13) << seq2.length() << "|\n";
    cout << "   |-------------------------------------------------------|\n";

    cout << "   |                    DISTANCE: " << setfill(' ') << left << setw(16) << distance << "         |\n";
    cout << "   ---------------------------------------------------------\n\n";
}