#include <stdlib.h>
#include <limits>
#include <iostream>
#include "callfitness.h"
#include "powell.cpp"
#include "amoeba.cpp"
#include "quasi.cpp"
#include "psfminimise.h"


    using namespace std;
    
    double psf_fitnessStat;
    int psf_its = 0;
    int psf_itmax = 5;
    
    vector<double> psf_callAmoeba(vector<double> psf_estimates) {
        
        Amoeba am(numeric_limits<double>::epsilon());
        vector<double> del(5);
        del[0] = 0.1;
        del[1] = 0.1;
        del[2] = 0.1;
        del[3] = 0.1;
        del[4] = 0.1;
        while ( psf_its < psf_itmax ) {
            psf_estimates = am.minimise(psf_estimates,del,PSFfitness);
            psf_fitnessStat = fitness(psf_estimates);
            psf_its++;
        }
        psf_its = 0;
        return psf_estimates;
        
    }
    
    vector<double> psf_callPowell(vector<double> psf_estimates) {
        
        Powell<double(vector<double>&)> powell(PSFfitness,0.1);
        while ( psf_its < psf_itmax ) {
            psf_estimates = powell.minimise(psf_estimates);
            psf_fitnessStat = PSFfitness(psf_estimates);
            psf_its++;
        }
        psf_its = 0;
        return psf_estimates;
            
    }     
