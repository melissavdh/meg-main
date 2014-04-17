//g++ profile.cpp -o profile -I/usr/include/cfitsio -lcfitsio
// histogram - sum up column values in rectangle and plot those?

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <CCfits/CCfits>
#include "F1dim.cpp"
#include "brent.cpp"

    using namespace std;
    
    template <class T>
    struct Linemethod {
        
        vector<double> p;
        vector<double> xi;
        T &func;
        int n;
        Linemethod(T &funcc) : func(funcc) {}
        
        double linmin() {
            
            //cout << "linmin" << endl;
            
            double ax,xx,xmin;
            n = p.size();
            F1dim<T> f1dim(p,xi,func);
            ax = 0.0;
            xx = 1.0;
            Brent brent;
            brent.bracket(ax,xx,f1dim);
            xmin = brent.minimise(f1dim);
            for (int i = 0; i < n; i++) {
                xi[i] *= xmin;
                p[i] += xi[i];
            }
            return brent.fmin;
        }
    };
            
            
            
            
