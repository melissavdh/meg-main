#include <cmath>
#include <CCfits/CCfits>
#include <limits>
#include "powell.h"
#include "min_funcs.h"

    using namespace std;

    template <class T>    
    void Bracketmethod::bracket(const double a, const double b, T &func) {
        
        //cout << "bracket" << endl;
        
        const double GOLD = 1.618034;
        const double GLIMIT = 100.0;
        const double TINY = 0.00000000000000000001;
        
        ax = a;
        bx = b;
        double fu;
        fa = func(ax);
        fb = func(bx);
        
        if ( fb > fa ) {
            SWAP(ax,bx);
            SWAP(fb,fa);
        }
        
        cx = bx+GOLD*(bx-ax);
        fc = func(cx);
        
        while ( fb > fc ) {
            double r = (bx-ax)*(fb-fc);
            double q = (bx-cx)*(fb-fa);
            double u = bx-((bx-cx)*q-(bx-ax)*r)/(2.0*SIGN(MAX(abs(q-r),TINY),q-r));
            double ulim = bx+GLIMIT*(cx-bx);
            
            if ( (bx - u)*(u - cx) > 0.0 ) {
                fu = func(u);
                if ( fu < fc ) {
                    ax = bx;
                    bx = u;
                    fa = fb;
                    fb = fu;
                    return;
                } else if ( fu > fb ) {
                    cx = u;
                    fc = fu;
                    return;
                }
                u = cx + GOLD * (cx - bx);
                fu = func(u);
            } else if ( (cx - u)*(u - ulim) > 0.0 ) {
                fu = func(u);
                if ( fu < fc ) {
                    shft3(bx,cx,u,u+GOLD*(u-cx));
                    shft3(fb,fc,fu,func(u));
                }
            } else if ( (u - ulim)*(ulim - cx) >= 0.0 ) {
                u = ulim;
                fu = func(u);
            } else {
                u = cx + GOLD*(cx-bx);
                fu = func(u);
            }
            shft3(ax,bx,cx,u);
            shft3(fa,fb,fc,fu);
        }
    }


    Brent::Brent(const double toll = /*0.000000000000001yes*//*0.0000000000000001abort*/0.00000003) : tol(toll) {}
        
    template <class T>
    double Brent::minimise(T &func) {
        
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
    
    
    template <class T>
    F1dim<T>::F1dim(vector<double> &pp, vector<double> &xii, T &funcc) : p(pp), xi(xii), n(pp.size()), func(funcc), xt(n) {}
    
    template <class T>
    double F1dim<T>::operator() (const double x) {
        
        //cout << "F1dim" << endl;
        
        for (int i = 0; i < n; i++) {
            xt[i] = p[i] + x*xi[i];
        }
        return func(xt);
        
    }
    
    
    template <class T>
    Linemethod<T>::Linemethod(T &funcc) : func(funcc) {}
    
    template <class T>
    double Linemethod<T>::linmin() {
        
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
    
    
    template <class T>
    Powell<T>::Powell(T &func, const double ftoll = 0.00000003) : Linemethod<T>(func), ftol(ftoll) {}
    
    template <class T>
    vector<double> Powell<T>::minimise(vector<double> &pp) {
            
        //cout << "powell" << endl;
        
        int n = pp.size();
        MatDoub ximat(n,n,0.0);
        for (int i = 0; i < n; i++) ximat[i][i] = 1.0;
        return minimise(pp,ximat);
        
    }

    template <class T>
    vector<double> Powell<T>::minimise(vector<double> &pp, MatDoub_IO &ximat) {
            
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
        
