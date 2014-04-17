#ifndef __FINDIFF_H_INCLUDED__
#define __FINDIFF_H_INCLUDED__

#include <vector>

    struct Ftor {
        
        double operator ()(std::vector<double> & x) const;

    };
    
    
    template <class T>
    struct Funcd {
        
        double EPS;
        T &func;
        double f;
        Funcd(T &funcc) : EPS(1.0e-8), func(funcc) {}
        
        double operator() (std::vector<double> &x);
        
        void df(std::vector<double> &x, std::vector<double> &df);
        
    };
    
#endif
