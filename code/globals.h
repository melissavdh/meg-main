#ifndef __GLOBALS_H_INCLUDED__
#define __GLOBALS_H_INCLUDED__

#include <vector>

    typedef std::vector<std::vector<double> > pixelArray;
    
    extern pixelArray actual_values;
    extern pixelArray coords;
    extern pixelArray psf_values;
    extern pixelArray psf_coords;
    extern int psfimage;
    extern pixelArray psfkernel;

#endif
