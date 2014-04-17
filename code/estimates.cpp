#include "estimates.h"
#include "globals.h"

    using namespace std;

    typedef vector<vector<double> > pixelArray;
    
    // histogram of pixel values to estimate b
    double estimate_b(pixelArray all_values) {
        
        double back = 0;
        int bins = 200;
        int mode = 0;
        
        double min = 10000;
        double max = 0;
        
        int yrange = all_values.size();
        int xrange = all_values[0].size();
        
        for (int i = 0; i < yrange; i++) {
            for (int j = 0; j < xrange; j++) {
                if ( (all_values[i][j] != 0) && !(isnan(all_values[i][j])) && !(isinf(all_values[i][j])) ) {
                    if ( all_values[i][j] < min ) min = all_values[i][j];
                    else if ( all_values[i][j] > max ) max = all_values[i][j];
                }
            }
        }
        
        double width = (max - min)/bins;
        width = ceil(width);
        
        vector<int> bin_values(bins);
        
        for (int i = 0; i < yrange; i++) {
            for (int j = 0; j < xrange; j++) {
                int index = floor((all_values[i][j]-min)/width);
                bin_values[index]++;
            }
        }
        for (int i = 0; i < bins; i++) {
            if (bin_values[i] > mode) {
                mode = bin_values[i];
                back = ((min+(i*width))+(min+((i+1)*width)))/2;
            }
        }

        return back;
        
    }

    // estimate scale radius: distance at which the central brightness has dropped by a factor e
    double estimate_rs(int x,int y,pixelArray actual_values,double l0,double q,double b) {
    
        double rs = 0;
        double drop = l0/exp(1);
        double current = l0;
        int newx = x;
        
        while ( current > drop ) {
            current = actual_values[y][newx] - b;
            newx++;
        }
        
        return rs = sqrt( pow((newx-2-x)/q,2) );
        
    }
   
    // estimate r_e from r_s
    /*double estimate_re(double r_s, double n) {
        
        double k = (2*n) - (1/3) + (0.009876/n);
        double r_e = r_s*(pow(k,n));
        
        return r_e;
        
    }*/


