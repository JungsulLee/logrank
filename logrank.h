#ifndef __DEFINITION_OF_LOGRANK__
#define __DEFINITION_OF_LOGRANK__

#include <vector>

// [group1] and [group2] : time_censoring. 0 (false) if censored or lost follow-up. 1 (true) if event occured. 
bool logrank(const std::vector<std::pair<double,bool> >& group1, const std::vector<std::pair<double,bool> >& group2, double* const z_score, double* const pvalue, double* const chisquare); 


#endif
