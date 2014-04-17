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
#include "callminimise.h"
#include "profile.h"
#include "globals.h"
#include "fitGAL.h"


    using namespace std;
        
    typedef vector<vector<double> > pixelArray;
    pixelArray actual_values;
    pixelArray coords;
    pixelArray psfkernel;

    int fit_galaxy(char *argv[], double sigma) {
        
        // set Sersic index to average value (between 1 and 4)
        double n = 2.5;
        
        // read in coordinates and put in 2d array coords to be passed to other
        // functions called from here.
        coords = readin(argv[1]);
        if (coords[0][0] == coords[1][1]) {
            cerr << "ERROR: Problem with co-ordinates file, check that the co-ordinates describe a sensible region." << endl;
            exit(1);
        }
            
        // calculate dimensions of region       
        int xrange = (coords[1][0] - coords[0][0])+1;
        int yrange = (coords[1][1] - coords[0][1])+1;
        
        // calculate central co-ords of region
        double array_x = (xrange-1)/2;
        double array_y = (yrange-1)/2;
        
        // calculate central co-ords of array
        int x =  (int)(array_x + 0.5);
        int y = (int)(array_y + 0.5);
        
        // get pixel values from region of FITS file
        actual_values = read_values(argv[0],coords,xrange,yrange);
        
        // estimate b
        double b = estimate_b(actual_values);
        
        // l0 is pixel value at central co-ords
        double l0 = (actual_values[y][x] - b);
        
        // set axis ratio b/a (assume 1 at start)
        double q = 1.0;

        // estimate rs by finding radius at which light has dropped by a factor e
        double r_s = estimate_rs(x,y,actual_values,l0,q,b);
                
        // set position angle theta (assume 0 at start) (in radians)
        double theta = 0.0;
        
        // set height of delta function (assume 0.5 at start)
        double delta = 0.1;
        
        // if PSF is input as image, read in image
        if ( psfimage == 0 ) psfkernel = read_all_values(argv[2]);
        
        // calculate fitness
        vector<double> estimates(10);
        estimates[0] = n;
        estimates[1] = b;
        estimates[2] = r_s;
        estimates[3] = array_x;
        estimates[4] = array_y;   
        estimates[5] = q;
        estimates[6] = theta;
        estimates[7] = l0;
        estimates[8] = sigma;
        estimates[9] = delta;
        
        /******************** FINAL PARAMETER VALUES ******************** delta function and convolution
        Fitness of model          =  10.1171
        Sersic index              =  1.38226
        Background level          =  1880.29
        Scale radius              =  4.34478
        Central co-ordinates      =  (73.0701,67.1897)
        Axis ratio                =  0.997556
        Position angle (degrees)  =  2.78799
        Central brightness        =  28062.3
        PSF sigma                 =  1.5
        Height of delta function  =  0.999137
        ****************************************************************/
        /******************** FINAL PARAMETER VALUES ******************** delta function no convolution
        Fitness of model          =  3.70023
        Sersic index              =  1.67435
        Background level          =  1903.94
        Scale radius              =  2.49611
        Central co-ordinates      =  (73.2988,67.0977)
        Axis ratio                =  0.851595
        Position angle (degrees)  =  58.5764
        Central brightness        =  49294.7
        PSF sigma                 =  2.31157
        Height of delta function  =  0.999311
        ****************************************************************/
        
        /******************** FINAL PARAMETER VALUES ******************** fpC no convolution
        Fitness of model          =  0.121091
        Sersic index              =  1.49405
        Background level          =  1160.88
        Scale radius              =  5.86691
        Central co-ordinates      =  (94.6468,96.7897)
        Axis ratio                =  0.933541
        Position angle (degrees)  =  31.719
        Central brightness        =  2577.17
        PSF sigma                 =  3.46692
        Height of delta function  =  1
        ****************************************************************/
        /******************** FINAL PARAMETER VALUES ******************** fpC convolution
        Fitness of model          =  0.0707993
        Sersic index              =  2.194
        Background level          =  1154.59
        Scale radius              =  1.53685
        Central co-ordinates      =  (94.5,96.5805)
        Axis ratio                =  0.94916
        Position angle (degrees)  =  32.5894
        Central brightness        =  6337.82
        PSF sigma                 =  3.07585
        Height of delta function  =  0.996107
        ****************************************************************/
        /******************** FINAL PARAMETER VALUES ******************** amoeba, bobyqa
        Fitness of model          =  8.17722
        Sersic index              =  1.34671
        Background level          =  1932.47
        Scale radius              =  4.53047
        Central co-ordinates      =  (73.4928,67.5)
        Axis ratio                =  0.893054
        Position angle (degrees)  =  59.9024
        Central brightness        =  28940.2
        PSF sigma                 =  1.66809
        Height of delta function  =  1.66149
        ****************************************************************/
        
        /*estimates[0] = 1.59732;
        estimates[1] =  1912;
        estimates[2] =  2.81635;
        estimates[3] = 72.731;
        estimates[4] = 66.532;
        estimates[5] =  0.846082;
        estimates[6] =  0.9455879728;
        estimates[7] =  45690.1;
        estimates[8] =  1.70901;
        estimates[9] =  3.90346e-06;*/
        
        /*estimates[0] = 1.84675;
        estimates[1] = 1904.94;
        estimates[2] = 1.71001;
        estimates[3] = 72.659;
        estimates[4] = 66.911;   
        estimates[5] = 0.845812;
        estimates[6] = 0.33098666*M_PI;
        estimates[7] = 66056.2;
        estimates[8] = 1.61598;
        estimates[9] = 9.60272e-06;*/
        
        if ( string(argv[3]) == "amoeba" ) estimates = callAmoeba(estimates);
        else if ( string(argv[3]) == "powell" ) estimates = callPowell(estimates);
        else if ( string(argv[3]) == "bobyqa" ) estimates = callBobyqa(estimates);
        else if ( string(argv[3]) == "lm" ) estimates = callLM(estimates);
        else if ( string(argv[3]) == "none" ) cout << "No minimisation chosen.  Modelling initial estimates." << endl;
        else {
            cerr << "ERROR: Minimisation method \"" << argv[3] << "\" must be one of amoeba / powell / bobyqa / lm / none." << endl;
            exit(1);
        }
        
        // write ROI back into image and save as new file
        double final_statistic = fitness(estimates);
        pixelArray sersic_values = newProfile(xrange,yrange,estimates);
        pixelArray new_region = createRegion(actual_values,sersic_values,xrange,yrange);
        newImage(argv[0], new_region, xrange, yrange, coords);
        
        // write Sersic model into new file
        if ( string(argv[4]) != "n" ) sersicImage(argv[0], sersic_values, xrange, yrange);
        
        if ( string(argv[3]) == "none" ) {
            
        cout << endl << "******************* INITIAL PARAMETER VALUES *******************" << endl;
        cout << "Fitness of model           =  " << final_statistic << endl;
        cout << "Sersic index               =  " << estimates[0] << endl;
        cout << "Background level           =  " << estimates[1] << endl;
        cout << "Scale radius               =  " << estimates[2] << endl;
        cout << "Central co-ordinates       =  (" << (coords[0][0] + estimates[3]) << "," << (coords[0][1] + estimates[4]) << ")" << endl;
        cout << "Axis ratio                 =  " << estimates[5] << endl;
        cout << "Position angle (degrees)   =  " << estimates[6]*(180/M_PI) << endl;
        cout << "Central brightness         =  " << estimates[7] << endl;
        cout << "PSF sigma                  =  " << estimates[8] << endl;
        cout << "Delta height               =  " << estimates[9] << endl;
        cout << "****************************************************************" << endl << endl;
        
        } else {
        
        cout << endl << "******************** FINAL PARAMETER VALUES ********************" << endl;
        cout << "Fitness of model           =  " << final_statistic << endl;
        cout << "Sersic index               =  " << estimates[0] << endl;
        cout << "Background level           =  " << estimates[1] << endl;
        cout << "Scale radius               =  " << estimates[2] << endl;
        cout << "Central co-ordinates       =  (" << (coords[0][0] + estimates[3]) << "," << (coords[0][1] + estimates[4]) << ")" << endl;
        cout << "Axis ratio                 =  " << estimates[5] << endl;
        cout << "Position angle (degrees)   =  " << estimates[6]*(180/M_PI) << endl;
        cout << "Central brightness         =  " << estimates[7] << endl;
        cout << "PSF sigma                  =  " << estimates[8] << endl;
        cout << "Delta height               =  " << estimates[9] << endl;
        cout << "****************************************************************" << endl << endl;
        
        }
        
       
        return 0;
           
    }
	
