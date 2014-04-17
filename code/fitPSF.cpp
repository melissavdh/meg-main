#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include "estimates.h"
#include "megio.h"
#include "callfitness.h"
#include "psfminimise.h"
#include "profile.h"
#include "fitPSF.h"


    using namespace std;
        
    typedef vector<vector<double> > pixelArray;
    pixelArray psf_values;
    pixelArray psf_coords;

    double fit_psf(char *argv[]) {
                
        // read in coordinates and put in 2d array coords to be passed to other
        // functions called from here.
        psf_coords = readin(argv[2]);
        if (psf_coords[0][0] == psf_coords[1][1]) {
            cerr << "ERROR: Problem with co-ordinates file, check that the co-ordinates describe a sensible region." << endl;
            exit(1);
        }
            
        // calculate size of region and central co-ords in PIXEL CO-ORDS        
        int psf_xrange = (psf_coords[1][0] - psf_coords[0][0])+1;
        int psf_yrange = (psf_coords[1][1] - psf_coords[0][1])+1;
        int psf_xx = floor(psf_xrange/2 + psf_coords[0][0]);
        int psf_yy = floor(psf_yrange/2 + psf_coords[0][1]);
        
        // calculate central co-ords in terms of ARRAY INDICES
        int psf_x = psf_xx - psf_coords[0][0];
        int psf_y = psf_yy - psf_coords[0][1];

        // get pixel values from region of FITS file
        psf_values = read_values(argv[0],psf_coords,psf_xrange,psf_yrange);

        // estimate b
        //pixelArray all_values = read_all_values(argv[0]);
        double psf_b = estimate_b(psf_values);
        
        // get l0 by finding pixel value at centre
        double psf_l0 = psf_values[psf_y][psf_x] - psf_b;
        
        double sum = 0;
        for (int i = 0; i < psf_yrange; i++) {
            for (int j = 0; j < psf_xrange; j++) {
                sum += (psf_values[i][j] - psf_b);
            }
        }
        //for (int i = -15; i <= 15; i++) cout << i << " " << (psf_values[16][16+i]-psf_b)/sum << endl;
        
        // estimate sigma
        double sigma_est = 1;
        
        vector<double> psf_estimates(5);
        psf_estimates[0] = psf_x;
        psf_estimates[1] = psf_y;
        psf_estimates[2] = psf_l0;
        psf_estimates[3] = psf_b;
        psf_estimates[4] = sigma_est;
        
        //cout << "x = " << psf_estimates[0] << ", y = " << psf_estimates[1] << ", l0 = " << psf_estimates[2] << endl;
        //cout << "b = " << psf_estimates[3] << ", sigma = " << psf_estimates[4] << endl;
        
        psf_estimates = psf_callPowell(psf_estimates);
        
        //cout << "x = " << psf_estimates[0] << ", y = " << psf_estimates[1] << ", l0 = " << psf_estimates[2] << endl;
        //cout << "b = " << psf_estimates[3] << ", sigma = " << psf_estimates[4] << endl;
        
        //pixelArray profile_values = psf_profile(psf_xrange,psf_yrange,psf_estimates);
        // for sersic values only
        //sersicImage(argv[1], profile_values, psf_xrange, psf_yrange);
        
        return psf_estimates[4];

           
    }
	
