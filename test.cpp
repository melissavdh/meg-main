//g++ profile.cpp -o profile -I/usr/include/cfitsio -lcfitsio
// histogram - sum up column values in rectangle and plot those?

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include "fitness.cpp"
#include "minimise.cpp"

    using namespace std;
    
    int main(int argc, char *argv[]) {

        vector<double> p(1);
        p[0] = 600;
        //p[1] = 1;
        Powell<double(vector<double>&)> powell(function);
        p = powell.minimise(p);
        
        cout << "p[0] = " << p[0] << /*", p[1] = " << p[1] <<*/ endl;
        
        return 0;
        
    }
