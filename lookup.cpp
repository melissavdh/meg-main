#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

    using namespace std;
    
    double q[7] = {
        75122.6331530,
        80916.6278952,
        36308.2951477,
        8687.24529705,
        1168.92649479,
        83.8676043424,
        2.50662827511
    };
    
    double value[1000001];
    
    double GAMMA(double z) {
        
        double num = 0;
        double den = 1;
        
        for (int i = 0; i < 7; i++) {
            num += (q[i])*(pow(z,i));
            den *= (z + i);
        }
        
        double frac = num/den;
        double b = pow((z + 5.5),(z + 0.5));
        double c = exp((-1)*(z + 5.5));
        
        double gamma = frac * b * c;
        
        return gamma;
        
    }
    
    
    void fillTable() {
        
        for (int i = 0; i < 1000001; i++) {
            double n = i*0.0001;
            double insert = GAMMA(n);
            value[i] = insert;
        }
        
    }
    
    
    double gamma(double x) {
        
        double result = x/0.0001;
        int val = (int)(result+0.5);
        return value[val];
        
    }

//http://www.rskey.org/gamma.htm
    

