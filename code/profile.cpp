#include <stdlib.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include "FFT.h"
#include "profile.h"

    using namespace std;
    
    typedef vector<vector<double> > pixelArray;
    
    // apply delta function to pixel based on distance from centre
    double delta(double height,double distance) {
        
        double norm = 1.0/(height * sqrt(M_PI));
        double power = (pow(distance,2)/pow(height,2));
        double new_value = norm * exp(-power);
        
        return new_value;
        
    }

    // calculate sersic value of pixel based on distance from centre
    double contsersic(double rs, double n, double distance, double l0) {
        
        double contvalue = 0.0;

        contvalue = l0 * exp( (pow( (distance/rs),(1.0/n) ))*(-1.0) );
        
        return contvalue;
        
    }
    
    // calculate sersic profile in 360 degrees
    pixelArray newProfile(int xrange, int yrange, vector<double> estimates) {
        
        pixelArray sersic_values;
        sersic_values.resize(yrange);
        for (int i = 0; i < yrange; i++) {
            sersic_values[i].resize(xrange);
        }
        int array_x =  (int)(estimates[3]+0.5);
        int array_y = (int)(estimates[4]+0.5);
        double q = estimates[5];
        double PA = estimates[6];

        double distance = 0;
        for (int y = 0; y < yrange; y++) {
            for (int x = 0; x < xrange; x++) {
                double new_x = (x-array_x)*cos(PA) + (y-array_y)*sin(PA) + array_x;
                double new_y = (y-array_y)*cos(PA) - (x-array_x)*sin(PA) + array_y;
                distance = sqrt(pow((new_x-array_x),2) + pow(((new_y-array_y)/q),2));
                sersic_values[y][x] = contsersic(estimates[2],estimates[0],distance,estimates[7]);
                //sersic_values[y][x] += delta(estimates[9],distance);
            }
        }
        //cout << "before delta" << endl;
        //for (int i = -50; i <= 50; i++) cout << i << " " << sersic_values[array_y+i][array_x] << endl;
        
        //sersic_values[array_y][array_x] += delta(estimates[9],0);
        
        for (int y = 0; y < yrange; y++) {
            for (int x = 0; x < xrange; x++) {
                double new_x = (x-array_x)*cos(PA) + (y-array_y)*sin(PA) + array_x;
                double new_y = (y-array_y)*cos(PA) - (x-array_x)*sin(PA) + array_y;
                distance = sqrt(pow((new_x-array_x),2) + pow(((new_y-array_y)/q),2));
                sersic_values[y][x] += delta(estimates[9],distance);
            }
        }

        //cout << "after delta" << endl;
        //for (int i = -50; i <= 50; i++) cout << i << " " << sersic_values[array_y+i][array_x] << endl;
        
        // now convolve with PSF
        sersic_values = convolve(sersic_values,xrange,yrange,estimates[8]);
        //cout << "after convolution" << endl;
        //for (int i = -50; i <= 50; i++) cout << i << " " << sersic_values[array_y+i][array_x] << endl;
        
        return sersic_values;
        
    }
    
    // calculate new image region by subtracting sersic values from actual values
    pixelArray createRegion(pixelArray actual_values, pixelArray sersic_values, int xrange, int yrange) {
        
        pixelArray new_region;
        new_region.resize(yrange);
        for (int i = 0; i < yrange; i++) {
            new_region[i].resize(xrange);
        }
        
        for (int y = 0; y < yrange; y++) {
            for (int x = 0; x < xrange; x++) {
                new_region[y][x] = actual_values[y][x] - sersic_values[y][x];
            }
        }
        
        return new_region;
        
    }
    
    // calculate profile value of psf pixel based on distance from centre
    double contpsf(double sigma, double distance, double l0) {
        
        double contvalue = 0;

        contvalue = l0 * exp( ((pow(distance,2))/( 2*pow(sigma,2) )) * (-1) );
        
        return contvalue;
        
    }
    
    // calculate profile of psf star
    pixelArray psf_profile(int psf_xrange, int psf_yrange, vector<double> psf_estimates) {
        
        pixelArray psf_profile;
        psf_profile.resize(psf_yrange);
        for (int i = 0; i < psf_yrange; i++) {
            psf_profile[i].resize(psf_xrange);
        }
        int psf_x = (int)(psf_estimates[0]+0.5);
        int psf_y = (int)(psf_estimates[1]+0.5);

        double distance = 0;
        for (int y = 0; y < psf_yrange; y++) {
            for (int x = 0; x < psf_xrange; x++) {
                distance = sqrt(pow((x-psf_x),2) + pow((y-psf_y),2));
                psf_profile[y][x] = contpsf(psf_estimates[4],distance,psf_estimates[2]);
            }
        }
        
        return psf_profile;
        
    }
