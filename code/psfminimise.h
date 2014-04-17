#ifndef __PSFMINIMISE_H_INCLUDED__
#define __PSFMINIMISE_H_INCLUDED__

#include <vector>
    
    std::vector<double> psf_callAmoeba(std::vector<double> estimates);
    
    std::vector<double> psf_callPowell(std::vector<double> estimates);
    
#endif
