#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <limits>
#include "minfuncs.h"
#include "matrix.h"
#include "linemethod.h"

    using namespace std;

    template <class T>
    struct Powell : Linemethod<T> {
        
        int iter;
        double fret;
        using Linemethod<T>::func;
        using Linemethod<T>::linmin;
        using Linemethod<T>::p;
        using Linemethod<T>::xi;
        const double ftol;
        Powell(T &func, const double ftoll = 0.00000003) : Linemethod<T>(func), ftol(ftoll) {}
        
        vector<double> minimise(vector<double> &pp) {
            
            //cout << "powell" << endl;
            
            int n = pp.size();
            MatDoub ximat(n,n,0.0);
            for (int i = 0; i < n; i++) ximat[i][i] = 1.0;
            return minimise(pp,ximat);
            
        }

        vector<double> minimise(vector<double> &pp, MatDoub_IO &ximat) {
            
            const int ITMAX = 200;
            const double TINY = 0.0000000000000000000000001;
            
            double fptt;
            int n = pp.size(); // 4
            p = pp;
            vector<double> pt(n);
            vector<double> ptt(n);
            xi.resize(n);
            fret = func(p); // function to be minimised
            
            for (int i = 0; i < n; i++) pt[i] = p[i]; // save the initial point
            
            for (iter = 0;;++iter) {
                //cout << "iter = " << iter << endl;
                double fp = fret;
                int ibig = 0;
                double del = 0.0;
                    for (int k = 0; k < n; k++) {
                       for (int l = 0; l < n; l++) xi[l] = ximat[l][k]; // copy the direction
                    fptt = fret;
                    fret = linmin();
                    if ( (fptt - fret) > del ) {
                       del = fptt - fret;
                       ibig = k + 1;
                    }
                }
                if ( 2.0*(fp-fret) <= (ftol*(abs(fp) + abs(fret) )/* + TINY*/) ) {
                    //cout << "returning estimates" << endl;
                    return p;
                }
                if ( iter == ITMAX ) throw("powell exceeding maximum iterations.");
                for (int j = 0; j < n; j++) {
                    ptt[j] = 2.0*p[j] - pt[j];
                    xi[j] = p[j] - pt[j];
                    pt[j] = p[j];
                }
                fptt = func(ptt);
                if ( fptt < fp ) {
                    double t = 2.0 * (fp - 2.0*fret + fptt) * pow((fp - fret - del),2) - del*pow((fp - fptt),2);
                    if ( t < 0.0 ) {
                        fret = linmin();
                        for (int m = 0; m < n; m++) {
                            ximat[m][ibig - 1] = ximat[m][n-1];
                            ximat[m][n-1] = xi[m];
                        }
                    }
                }
            }
        }
    };
