#ifndef __LINEMETHOD_H_INCLUDED__
#define __LINEMETHOD_H_INCLUDED__

#include <vector>
       
    struct Bracketmethod {
        
        double ax,bx,cx,fa,fb,fc;
        template <class T>    
        void bracket(const double a, const double b, T &func);
        
    };


    struct Brent : Bracketmethod {
        
        double xmin, fmin;
        const double tol;
        Brent(const double toll = /*0.000000000000001yes*//*0.0000000000000001abort*/0.00000003) : tol(toll) {}
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
        F1dim(std::vector<double> &pp, std::vector<double> &xii, T &funcc) : p(pp), xi(xii), n(pp.size()), func(funcc), xt(n) {}
        double operator() (const double x);
        
    };
    
    template <class T>
    struct Linemethod {
        
        std::vector<double> p;
        std::vector<double> xi;
        T &func;
        int n;
        Linemethod(T &funcc) : func(funcc) {}
        double linmin();
        
    };
    
    
#endif    
