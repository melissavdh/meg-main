#include <vector>
#include <cmath>
#include <algorithm>
#include "minfuncs.h"
#include "linemethod.h"

    using namespace std;
    
    struct Bracketmethod {
        
        double ax,bx,cx,fa,fb,fc;
   
        template <class T>    
        void bracket(const double a, const double b, T &func) {
            
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
    };


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

    
    template <class T>
    struct Frprmn : Linemethod<T> {
        
        int iter;
        double fret;
        using Linemethod<T>::func;
        using Linemethod<T>::linmin;
        using Linemethod<T>::p;
        using Linemethod<T>::xi;
        const double ftol;
        
        Frprmn(T &funcd, const double ftoll = 3.0e-8) : Linemethod<T>(funcd), ftol(ftoll) {}
        
        vector<double> minimise(vector<double> &pp) {
            
            const int ITMAX = 200;
            const double EPS = 1.0e-18;
            const double GTOL = 1.0e-8;
            double gg, dgg;
            int n = pp.size();
            p = pp;
            vector<double> g(n), h(n);
            
            xi.resize(n);
            double fp = func(p);
            func.df(p,xi);
            
            for (int j = 0; j < n; j++) {
                g[j] = -xi[j];
                xi[j] = h[j] = g[j];
            }
            
            for (int its = 0; its < ITMAX; its++) {
                iter = its;
                fret = linmin();
                if ( 2.0*abs(fret-fp) <= ftol*(abs(fret) + abs(fp) + EPS) ) return p;
                fp = fret;
                func.df(p,xi);
                double test = 0.0;
                double den = max(fp,1.0);
                
                for (int j = 0; j < n; j++) {
                    double temp = abs(xi[j]) * max(abs(p[j]),1.0)/den;
                    if ( temp > test ) test = temp;
                }
                
                if ( test < GTOL ) return p;
                dgg = gg = 0.0;
                
                for (int j = 0; j < n; j++) {
                    gg += g[j]*g[j];
                    dgg += xi[j]*xi[j];
                }
                
                if ( gg == 0.0 ) return p;
                double gam = dgg/gg;
                for (int j = 0; j < n; j++) {
                    g[j] = -xi[j];
                    xi[j] = h[j] = g[j] + gam*h[j];
                }
            }
            throw("Too many iterations in frprmn");
        }
    };
            
            