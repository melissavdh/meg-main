#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include "min_funcs.cpp"

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
            
