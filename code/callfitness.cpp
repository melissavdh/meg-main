#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include "estimates.h"
#include "profile.h"
#include "globals.h"
#include "callfitness.h"

    using namespace std;
    
    typedef vector<vector<double> > pixelArray;
    int count = 1;
    
    // calculate fitness statistic
    double fitnessStatistic(pixelArray actual_values,pixelArray sersic_values,vector<double>& estimates,int xrange,int yrange) {
        
        double statistic = 0;
        
        double B = estimates[1];
        
        for (int y = 0; y < yrange; y++) {
            for (int x = 0; x < xrange; x++) {
                double num = pow((actual_values[y][x] - (sersic_values[y][x]+B)),2);
                double den = actual_values[y][x];
                statistic += num/den;
            }
        }
       
        statistic /= ((xrange*yrange) - estimates.size());
        
        return statistic;
        
    }
    
    double fitness(vector<double>& estimates) {
        
        // translate incorrect estimates to correct (eg rotate 180 degrees to stay within boundary)
        // to prevent NaNs etc
        if (estimates[6] < 0) estimates[6] = M_PI - estimates[6];
        if (estimates[6] >= M_PI) estimates[6] -= M_PI;
        if (estimates[8] < 0) estimates[8] = 0.001;
        for (int i = 0; i < 10; i++) if (estimates[i] < 0) estimates[i] = 0.0000000000000001;
        
        // calculate size of region and central co-ordinates in PIXEL CO-ORDS        
        int xrange = (coords[1][0] - coords[0][0])+1;
        int yrange = (coords[1][1] - coords[0][1])+1;
        
        // calculate new Sersic model based on current set of estimates
        pixelArray sersic_values = newProfile(xrange,yrange,estimates);
        
        int array_x =  (int)(estimates[3]+0.5);
        int array_y = (int)(estimates[4]+0.5);
        
        //cout << "actual values" << endl;
        //for (int i = -50; i <= 50; i++) cout << i << " " << actual_values[array_y+i][array_x]-estimates[1] << endl;
        
        // calculate fitness statistic for this new Sersic model
        double statistic = fitnessStatistic(actual_values,sersic_values,estimates,xrange,yrange);
        
        //cout << "sersic values" << endl;
        //for (int i = -50; i <= 50; i++) cout << i << " " << sersic_values[array_y+i][array_x] << endl;
        
        // penalty constraints
        //if ( (estimates[5] > 1) || (estimates[5] <= 0) ) statistic += 100000000;
        //if (estimates[7] > 100000) statistic += 100000000;
        //if (estimates[9] > 0.2) statistic += 100000000;
        
        if (statistic < 100000) cout << count << " " << statistic << endl;
        count++;
        
        return statistic;
        
    }
    
    // calculate fitness statistic
    double psf_fitnessStatistic(pixelArray psf_values,pixelArray psf,vector<double>& psf_estimates,int psf_xrange,int psf_yrange) {
        
        double psf_statistic = 0;
        
        double b = psf_estimates[3];
        
        for (int y = 0; y < psf_yrange; y++) {
            for (int x = 0; x < psf_xrange; x++) {
                double num = pow((psf_values[y][x] - (psf[y][x]+b)),2);
                double den = psf_values[y][x];
                psf_statistic += num/den;
            }
        }
        
        psf_statistic /= ((psf_xrange*psf_yrange) - psf_estimates.size());
        
        return psf_statistic;
        
    }
    
    double PSFfitness(vector<double>& psf_estimates) {
                                 
        // calculate size of region and central co-ordinates in PIXEL CO-ORDS        
        int psf_xrange = (psf_coords[1][0] - psf_coords[0][0])+1;
        int psf_yrange = (psf_coords[1][1] - psf_coords[0][1])+1;
        
        // calculate new psf profile based on current set of estimates
        pixelArray psf = psf_profile(psf_xrange,psf_yrange,psf_estimates);
        
        // calculate fitness statistic for this new Sersic model
        double psf_statistic = psf_fitnessStatistic(psf_values,psf,psf_estimates,psf_xrange,psf_yrange);
        
        return psf_statistic;
        
    }
