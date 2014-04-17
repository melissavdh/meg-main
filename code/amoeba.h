#ifndef __AMOEBA_H_INCLUDED__
#define __AMOEBA_H_INCLUDED__

#include <vector>
#include "matrix.h"


    struct Amoeba {
        
        const double ftol;
        int nfunc;
        int mpts;
        int ndim;
        double fmin;
        std::vector<double> y;
        MatDoub p;
        
        Amoeba(const double ftoll);
        
		template <class T>
		std::vector<double> minimise(std::vector<double> &point, const double del, T &func);
		
		template <class T>
		std::vector<double> minimise(std::vector<double> &point, std::vector<double> &dels, T &func);
		
		template <class T>
		std::vector<double> minimise(MatDoub &pp, T &func);
		
		inline void get_psum(MatDoub &p, std::vector<double> &psum);
		
		template <class T>
		double amotry(MatDoub &p, std::vector<double> &y, std::vector<double> &psum, const int ihi, const double fac, T &func);
    
    };
        
#endif
        
