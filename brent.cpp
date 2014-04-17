#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <limits>
#include "bracket.cpp"

    using namespace std;
    
    struct Brent : Bracketmethod {
        
        double xmin, fmin;
        const double tol;
        Brent(const double toll = /*0.000000000000001yes*//*0.0000000000000001abort*/0.00000003) : tol(toll) {}
        
        template <class T>
        double minimise(T &func) {
            
            //cout << "brent" << endl;
            
            const int ITMAX = 100;
            const double CGOLD = 0.3819660;
            const double ZEPS = numeric_limits<double>::epsilon()*0.001;
            
            double a,b,d=0.0,etemp,fu,fv,fw,fx;
            double p,q,r,tol1,tol2,u,v,w,x,xm;
            double e=0.0;
            
            a = (ax < cx ? ax : cx);
            b = (ax > cx ? ax : cx);
            x = w = v = bx;
            fw = fv = fx = func(x);
            for (int iter = 0; iter < ITMAX; iter++) {
                xm = 0.5 * (a+b);
                tol2 = 2.0 * (tol1=tol*abs(x)+ZEPS);
                if ( abs(x-xm) <= (tol2 - 0.5*(b-a)) ) {
                    fmin = fx;
                    return xmin = x;
                }
                if ( abs(e) > tol1 ) {
                    
                    r = (x-w)*(fx-fv);
                    q = (x-v)*(fx-fw);
                    p = (x-v)*q - (x-w)*r;
                    q = 2.0 * (q-r);
                    if ( q > 0.0 ) p = -p;
                    q = abs(q);
                    etemp = e;
                    e = d;
                    if ( abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x) ) {
                        d = CGOLD*(e = (x >= xm ? a-x : b-x) );
                    } else {
                        d = p/q;
                        u = x + d;
                        if ( u-a < tol2 || b-u < tol2 ) d = SIGN(tol1,xm-x);
                    }
                } else {
                    d = CGOLD*(e = (x >= xm ? a-x : b-x));
                }
                u = ( abs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
                fu = func(u);
                if ( fu <= fx ) {
                    if ( u >= x) a = x;
                    else b = x;
                    shft3(v,w,x,u);
                    shft3(fv,fw,fx,fu);
                } else {
                    if ( u < x ) a = u;
                    else b = u;
                    if ( fu <= fw || w == x ) {
                        v = w;
                        w = u;
                        fv = fw;
                        fw = fu;
                    } else if ( fu <= fv || v == x || v == w ) {
                        v = u;
                        fv = fu;
                    }
                }
            }
            throw("Too many iterations in brent");
        }
    };
                    
                    
            
