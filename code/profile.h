#ifndef __PROFILE_H_INCLUDED__
#define __PROFILE_H_INCLUDED__
    
    typedef std::vector<std::vector<double> > pixelArray;
    
    // calculate new value after applying delta function
    double delta(double height,double distance);

    // calculate sersic profile based on distance from centre
    double contsersic(double rs, double n, double distance, double l0);
    
    // calculate sersic profile in 360 degrees
    pixelArray newProfile(int xrange, int yrange, std::vector<double> estimates);
    
    // calculate new image region by subtracting sersic values from actual values
    pixelArray createRegion(pixelArray actual_values, pixelArray sersic_values, int xrange, int yrange);
    
    // calculate profile value of psf pixel based on distance from centre
    double contpsf(double sigma, double distance, double l0);
    
    // calculate profile of psf star
    pixelArray psf_profile(int psf_xrange, int psf_yrange, std::vector<double> psf_estimates);
    
#endif
