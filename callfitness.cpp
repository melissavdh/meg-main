//g++ profile.cpp -o profile -I/usr/include/cfitsio -lcfitsio
// histogram - sum up column values in rectangle and plot those?

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <CCfits/CCfits>

    using namespace std;
    
    typedef vector<vector<double> > pixelArray;
    extern pixelArray actual_values;
    extern pixelArray coords;
    
    double fitness(vector<double>& estimates) {
        
        double b = estimates[1];
                    
        // calculate size of region and central co-ordinates in PIXEL CO-ORDS        
        int xrange = (coords[1][0] - coords[0][0])+1;
        int yrange = (coords[1][1] - coords[0][1])+1;
        
        pixelArray sersic_values = newProfile(xrange,yrange,estimates);
        
        // calculate fitness statistic
        double statistic = fitnessStatistic(actual_values,sersic_values,xrange,yrange,b);
        cout << "fitness = " << statistic << endl;
        
        return statistic;
        
    }
