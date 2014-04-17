#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <CCfits/CCfits>
#include "amoeba.h"
#include "min_funcs.h"

    using namespace std;
    
	Amoeba::Amoeba(const double ftoll) : ftol(ftoll) {}
	
	
	template <class T>
	vector<double> Amoeba::minimise(vector<double> &point, const double del, T &func) {
		
		vector<double> dels(point.size(),del);
		return minimise(point,dels,func);
		
	}
	
	
	template <class T>
	vector<double> Amoeba::minimise(vector<double> &point, vector<double> &dels, T &func) {
		
		int ndim = point.size();
		MatDoub pp(ndim + 1, ndim);
		for (int i = 0; i < ndim + 1; i++) {
			for (int j = 0; j < ndim; j++) pp[i][j] = point[j];
			if ( i != 0 ) pp[i][i-1] += dels[i-1];
		}
		
		return minimise(pp,func);
		
	}
	
	
	template <class T>
	vector<double> Amoeba::minimise(MatDoub &pp, T &func) {
		
		const int NMAX = 5000;
		const double TINY = 1.0e-10;
		int ihi,ilo,inhi;
		
		mpts = pp.nrows();
		ndim = pp.ncols();
		
		vector<double> psum(ndim),pmin(ndim),x(ndim);
		p = pp;
		y.resize(mpts);
		for (int i = 0; i < mpts; i++) {
			for (int j = 0; j < ndim; j++) x[j] = p[i][j];
			y[i] = func(x);
		}
		
		nfunc = 0;
		get_psum(p,psum);
		for (;;) {
			ilo = 0;
			ihi = y[0]>y[1] ? (inhi = 1,0) : (inhi = 0,1);
			for (int i = 0; i < mpts; i++) {
				if ( y[i] <= y[ilo] ) ilo = i;
				if ( y[i] > y[ihi] ) {
					inhi = ihi;
					ihi = i;
				} else if ( (y[i] > y[inhi]) && (i != ihi) ) inhi = i;
			}
			
			double rtol = 2.0*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo])+TINY);
			
			if ( rtol < ftol ) {
				swap(y[0],y[ilo]);
				for (int i = 0; i < ndim; i++) {
					swap(p[0][i],p[ilo][i]);
					pmin[i] = p[0][i];
				}
				fmin = y[0];
				return pmin;
			}
			
			if ( nfunc >= NMAX ) throw("NMAX exceeded");
			nfunc += 2;
			double ytry = amotry(p,y,psum,ihi,-1.0,func);
			if ( ytry <= y[ilo] ) ytry = amotry(p,y,psum,ihi,2.0,func);
			else if ( ytry >= y[inhi] ) {
				double ysave = y[ihi];
				ytry = amotry(p,y,psum,ihi,0.5,func);
				if ( ytry >= ysave ) {
					for (int i = 0; i < mpts; i++) {
						if ( i != ilo ) {
							for (int j = 0; j < ndim; j++) p[i][j] = psum[j] = 0.5*(p[i][j] + p[ilo][j]);
							y[i] = func(psum);
						}
					}
					nfunc += ndim;
					get_psum(p,psum);
				}
			} else --nfunc;
		}
	}
	
	
	inline void Amoeba::get_psum(MatDoub &p, vector<double> &psum) {
		
		for (int j = 0; j < ndim; j++) {
			double sum = 0.0;
			for (int i = 0; i < mpts; i++) sum += p[i][j];
			psum[j] = sum;
		}
	}
	
	
	template <class T>
	double Amoeba::amotry(MatDoub &p, vector<double> &y, vector<double> &psum, const int ihi, const double fac, T &func) {
		
		vector<double> ptry(ndim);
		double fac1 = (1.0-fac)/ndim;
		double fac2 = fac1 - fac;
		for (int j = 0; j < ndim; j++) ptry[j] = psum[j]*fac1 - p[ihi][j]*fac2;
		
		double ytry = func(ptry);
		if ( ytry < y[ihi] ) {
			y[ihi] = ytry;
			for (int j = 0; j < ndim; j++) {
				psum[j] += ptry[j] - p[ihi][j];
				p[ihi][j] = ptry[j];
			}
		}
		
		return ytry;
		
	} 
        
        
