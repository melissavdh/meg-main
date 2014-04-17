#ifndef __ESTIMATES_H_INCLUDED__
#define __ESTIMATES_H_INCLUDED__

#include <vector>
#include <cmath>

    typedef std::vector<std::vector<double> > pixelArray;
    
    // set up histogram to estimate b
    double estimate_b(pixelArray all_values);
    
    // estimate half light radius by stepping 1 pixel at a time
    double estimate_rs(int x,int y,pixelArray actual_values,double l0,double q,double b); 
       
    // estimate r_e from r_s
    //double estimate_re(double r_s, double n);

#endif
