using namespace std;

#include <iostream>
#include "Blosum.h"
#include "Menu.h"
#include "SeqComparison.h"


int main()
{
    printMenu();
    bool loop = true;

    while (loop)
    {
        char YorN;
        SeqComparison seq; // constructor calculates distance
        seq.printResults();

        cout << "Would you like to try another sequence? (Y = new sequence, N = exit) ";
        cin >> YorN;

        if (YorN != 'Y' && YorN != 'y')
            break;

        cout << '\n';
    }

    return 0;
}