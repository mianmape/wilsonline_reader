#include <iostream>
#include <string>
#include "ipglasma.hpp"

using namespace std;

int main(int argc, char* argv[])
{
    // Arguments: filaname
    string fname = argv[1];
    IPGlasma ipglasma(fname, 0, BINARY);

    // Get Wilson line at x=1 GeV^-1, y=3 GeV^-1
    WilsonLine w = ipglasma.GetWilsonLine(1,3);
    cout << w.Element(1,2) << endl; // Print element [1,2], note that it is a complex number


    return 0;
}   