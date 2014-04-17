#ifndef __CALLFITNESS_H_INCLUDED__
#define __CALLFITNESS_H_INCLUDED__

#include <vector>
    
    typedef std::vector<std::vector<double> > pixelArray;
    
    // calculate fitness statistic
    double fitnessStatistic(pixelArray actual_values,pixelArray sersic_values,std::vector<double>& estimates,int xrange,int yrange);
    
    double fitness(std::vector<double>& estimates);
    
    // calculate psf fitness statistic
    double psf_fitnessStatistic(pixelArray psf_values,pixelArray psf,std::vector<double>& psf_estimates,int psf_xrange,int psf_yrange);
    
    double PSFfitness(std::vector<double>& psf_estimates);
    
#endif
