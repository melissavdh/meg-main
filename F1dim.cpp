#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <CCfits/CCfits>

    using namespace std;
    
    template <class T>
    struct F1dim {
        
        const vector<double> &p;
        const vector<double> &xi;
        int n;
        T &func;
        vector<double> xt;
        
        F1dim(vector<double> &pp, vector<double> &xii, T &funcc) : p(pp), xi(xii), n(pp.size()), func(funcc), xt(n) {}
        
        double operator() (const double x) {
            
            //cout << "F1dim" << endl;
            
            for (int i = 0; i < n; i++) {
                xt[i] = p[i] + x*xi[i];
            }
            return func(xt);
            
        }
        
    };
