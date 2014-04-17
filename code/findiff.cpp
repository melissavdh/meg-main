#include <vector>
#include <cmath>
#include "callfitness.h"
#include "findiff.h"

    using namespace std;
       
    double Ftor::operator ()(vector<double> & x) const {

        return fitness(x);
        
    }
    
    template <class T>
    double Funcd<T>::operator() (vector<double> &x) {
        //qnfunc++;
        return f = func(x);
    }
    
    template <class T>
    void Funcd<T>::df(vector<double> &x, vector<double> &df) {
        
        int n = x.size();
        vector<double> xh = x;
        double fold = f;
        for (int j = 0; j < n; j++) {
            double temp = x[j];
            double h = EPS*abs(temp);
            if ( h == 0.0 ) h = EPS;
            xh[j] = temp + h;
            h = xh[j] - temp;
            double fh = operator()(xh);
            xh[j] = temp;
            df[j] = (fh - fold)/h;
        }
    }

