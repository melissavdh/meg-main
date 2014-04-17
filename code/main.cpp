#include <iostream>
#include <string>
#include <stdio.h>
#include "fitGAL.h"
#include "fitPSF.h"
#include "globals.h"


    using namespace std;
    
    int psfimage;
    
    int main(int argc, char *argv[]) {
        
        char *args[5];
        double sigma = 0;
        
        if ( string(argv[1]) == "--help" ) {
            cout << endl << "Welcome to meg.  Usage: "<< endl << endl;
            cout << "-i image        : FITS file containing the galaxy to be modelled" << endl;
            cout << "-c coords       : text file with co-ords defining ROI in image" << endl;
            cout << "-p image        : FITS file containing PSF" << endl;
            cout << "    OR " << endl;
            cout << "   coords       : text file with co-ords of PSF region in image" << endl;
            cout << "-m minimisation : choice of amoeba / powell / bobyqa / lm / none" << endl;
            cout << "-o output model : \"n\" if Sersic model should not be saved into new image" << endl;
            cout << "                  (defaut = yes)" << endl << endl;
            cout << "Please note that any new images output will overwrite images with the same name." << endl;
            cout << "Residual image saved as image_residual.fit and model saved as image_model.fit" << endl;
            cout << "where image is the original name of the input image." << endl << endl;
            return 0;
        }
        
        if ( string(argv[1]) == "--version" ) {
            cout << "This is meg, version 1.0 (May 2011)"<< endl;
            return 0;
        }
        
        if ( (argc < 9) || (argc > 11) ) {
            cerr << "ERROR: invalid number of arguments.  See --help" << endl;
            return 1;
        }
        
        // put args into correct place in args array
        for (int i = 1; i < argc; i+=2) {
            if ( string(argv[i]) == "-i" ) args[0] = argv[i+1];
            else if ( string(argv[i]) == "-c" ) args[1] = argv[i+1];
            else if ( string(argv[i]) == "-p" ) args[2] = argv[i+1];
            else if ( string(argv[i]) == "-m" ) args[3] = argv[i+1];
            else if ( string(argv[i]) == "-o" ) args[4] = argv[i+1];
            else {
                cout << "ERROR: invalid input, check your flags.  See --help." << endl;
                return 1;
            }
        }
        
        // set psfimage based on type of PSF input, find sigma if PSF not input as image
        size_t fit = string(args[2]).find(".fit");
        if (fit == string::npos) {
            psfimage = 1;
            sigma = fit_psf(args);
            cout << "sigma = " << sigma << endl;
        } else psfimage = 0;
        
        // fit the galaxy
        if ( fit_galaxy(args,sigma) != 0 ) {
            cerr << "ERROR: problem fitting galaxy." << endl;
            return 1;
        }
    
        return 0;
    }
	
