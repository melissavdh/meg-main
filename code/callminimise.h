#ifndef __CALLMINIMISE_H_INCLUDED__
#define __CALLMINIMISE_H_INCLUDED__

#include <vector>
    
    std::vector<double> callAmoeba(std::vector<double> estimates);
    
    std::vector<double> callPowell(std::vector<double> estimates);
    
    void callQuasi(std::vector<double> estimates);
    
    std::vector<double> callBobyqa(std::vector<double> estimates);
    
    std::vector<double> callLM(std::vector<double> estimates);
    
#endif
