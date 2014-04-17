#ifndef __QUASI_H_INCLUDED__
#define __QUASI_H_INCLUDED__

#include <vector>

    //int qnfunc;
    
    struct Ftor {
        
        double operator ()(std::vector<double> & x) const;
        
    };
    
    
    template <class T>
    struct Funcd {
        
        double EPS;
        T &func;
        double f;
        Funcd(T &funcc);
        
        double operator() (std::vector<double> &x);
        void df(std::vector<double> &x, std::vector<double> &df);
        
    };

    
    template <class T>
    void lnsrch(std::vector<double> &xold, const double fold, std::vector<double> &g,
        std::vector<double> &p, std::vector<double> &x, double &f, const double stpmax,
        bool &check, T &func);
    
    
    template <class T>
    void dfpmin(std::vector<double> &p, const double gtol, int &iter, double &fret, T &funcd);
    
#endif
