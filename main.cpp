//g++ profile.cpp -o profile -I/usr/include/cfitsio -lcfitsio
// histogram - sum up column values in rectangle and plot those?

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <CCfits/CCfits>
#include "estimates.cpp"
#include "readin.cpp"
#include "callfitness.cpp"
#include "minimise.cpp"
//#include "histogram.cpp"

    using namespace std;
        
    typedef vector<vector<double> > pixelArray;
    pixelArray actual_values;
    pixelArray coords;

    int main(int argc, char *argv[]) {
                
        // create gamma function lookup table
        createLookUp();
        
        //histogram(argv[1]);
        //cout << "histo executed" << endl;
               
        double b = 1700;
        double n = 0.5;
        
        // read in coordinates and put in 2d array coords to be passed to other
        // functions called from here.
        cout << "argv[2] = " << argv[2] << endl;
        coords = readin(argv[2]);
        if (coords[0][0] == coords[1][1]) return 1;
            
        // calculate size of region and central co-ords in PIXEL CO-ORDS        
        int xrange = (coords[1][0] - coords[0][0])+1;
        int yrange = (coords[1][1] - coords[0][1])+1;
        int xx = floor(xrange/2 + coords[0][0]);
        int yy = floor(yrange/2 + coords[0][1]);
        
        // calculate central co-ords in terms of ARRAY INDICES
        int array_x = xx - coords[0][0];
        int array_y = yy - coords[0][1];
        
        // get pixel values from region of FITS file
        actual_values = read_values(argv[1],coords,xrange,yrange);
        
        // estimate total light using sky estimate and actual pixel values
        double totalLight = total_light(b,actual_values,xrange,yrange);

        // estimate rs by finding half light and corresponding distance
        double Rs = estimate_rs(array_x,array_y,actual_values,totalLight);
                
        // estimate L0 using total light and rs
        double start_light = estimate_l0(totalLight,Rs,n);
        
        // calculate fitness
        vector<double> estimates(6);
        estimates[0] = n;
        estimates[1] = b;
        estimates[2] = Rs;
        estimates[3] = start_light;
        estimates[4] = (double)array_x;
        estimates[5] = (double)array_y;        
        Powell<double(vector<double>&)> powell(fitness);
        estimates = powell.minimise(estimates);
        
        cout << "n = " << estimates[0] << ", b = " << estimates[1] << ", Rs = " << estimates[2] << ", L0 = " << estimates[3] << endl;
        cout << "array_x = " << estimates[4] << ", array_y = " << estimates[5] << endl;
        cout << "Number of iterations =  " << powell.iter << endl;
            
        pixelArray sersic_values = newProfile(xrange,yrange,estimates);
        
        pixelArray new_region = createRegion(actual_values,sersic_values,xrange,yrange);
               
        newImage(argv[1], new_region, xrange, yrange, coords);
               
        return 0;
           
    }
	
