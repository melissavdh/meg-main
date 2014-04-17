#ifndef __POWELL_H__
#define __POWELL_H__

#include <vector>
#include "matrix.h"
    
    struct Bracketmethod {
        
        double ax,bx,cx,fa,fb,fc;
        template <class T>    
        void bracket(const double a, const double b, T &func); 

    };
    
    
    struct Brent : Bracketmethod {
        
        double xmin, fmin;
        const double tol;
        Brent(const double toll);
        template <class T>
        double minimise(T &func);
        
    };
    
    
    template <class T>
    struct F1dim {
        
        const std::vector<double> &p;
        const std::vector<double> &xi;
        int n;
        T &func;
        std::vector<double> xt;
        F1dim(std::vector<double> &pp, std::vector<double> &xii, T &funcc);
        double operator() (const double x);
        
    };
    
    
    template <class T>
    struct Linemethod {
        
        std::vector<double> p;
        std::vector<double> xi;
        T &func;
        int n;
        Linemethod(T &funcc);        
        double linmin();
        
    };
    
    
    template <class T>
    struct Powell : Linemethod<T> {
        
        int iter;
        double fret;
        using Linemethod<T>::func;
        using Linemethod<T>::linmin;
        using Linemethod<T>::p;
        using Linemethod<T>::xi;
        const double ftol;
        Powell(T &func, const double ftoll);
        std::vector<double> minimise(std::vector<double> &pp);
        std::vector<double> minimise(std::vector<double> &pp, MatDoub_IO &ximat);
        
    };
    
#endif
