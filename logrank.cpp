#include <util.h>

#include <math.h>

typedef struct {
	int event; // 1 if occurred, 0 if censored. 
	double event_time; 
	int group; // group index, 0 or 1.
} events; 

typedef struct {
	int event; // 1 if occurred, 0 otherwise including censoring. 
	double event_time; 
	double value; // by which subject will be sorted in ascending order. 
	int group; // group index, 0 or 1. -1 if invalid. 
} subject; 

bool insert(const std::vector<std::pair<double,bool> >& group, const int& group_index, unsigned int* const valid_n, std::vector<events>* const dest)
{

	unsigned int n = 0; 
	std::vector<std::pair<double,bool> >::const_iterator pos;
	for(pos = group.begin(); pos != group.end(); pos++){
		if(pos->first < 1e-10) continue; 
		n++; 
		events e; 
		e.event = (int)(pos->second); 
		e.event_time = pos->first; 
		e.group = group_index;
		dest->push_back(e); 
	}
	*valid_n = n; 

	return true;
}

bool compare_with_event_time(const events& op1, const events& op2)
{
	if(fabs(op1.event_time - op2.event_time) < 1e-20){
		return op1.event < op2.event; 
	}
	return op1.event_time < op2.event_time;
}


// [group1] and [group2] : time_censoring. 0 if censored. 
bool logrank(const std::vector<std::pair<double,bool> >& group1, const std::vector<std::pair<double,bool> >& group2, double* const z_score, double* const pvalue, double* const chisquare)
{
	if(group1.empty() == true || group2.empty() == true) return false; 

	std::vector<events> pool; 
	unsigned int remained1 = 0; 
	unsigned int remained2 = 0; 
	insert(group1, 0, &remained1, &pool); 
	insert(group2, 1, &remained2, &pool); 
	
	if(remained1 == 0 || remained2 == 0) return false;	

	std::sort(pool.begin(), pool.end(), compare_with_event_time); 
	events invalid_event = pool[pool.size()-1]; 
	invalid_event.event_time = invalid_event.event_time + 1.0; 
	invalid_event.event = 0; 
	pool.push_back(invalid_event); 

	double O = 0.0; double E = 0.0; double V = 0.0; 
	unsigned int event1_at_j = 0; 
	unsigned int event2_at_j = 0; 


	double previous_event_time = pool[0].event_time;
	std::vector<events>::const_iterator pos;
	int tied1 = 0; int tied2 = 0; 
	for(pos = pool.begin(); pos != pool.end()-1; pos++){
		if(pos->event == 0){
			if(fabs(pos->event_time - ((pos+1)->event_time)) > 1e-20){
				if(pos->group == 0) remained1--; 
				else if(pos->group == 1) remained2--; 
				remained1 -= tied1; 
				remained2 -= tied2; 
				tied1 = 0; 
				tied2 = 0; 
			}
			else{
				if(pos->group == 0) tied1++; 
				else if(pos->group == 1) tied2++; 
			}
			continue; 
		}
		// Now that event is 1. 

		if(pos->group == 0){
			event1_at_j++; 
		}
		else if(pos->group == 1){
			event2_at_j++; 
		}
		
		if(fabs(pos->event_time - ((pos+1)->event_time)) < 1e-20){
			if(pos->group == 0){
				tied1++; 
			}
			else if(pos->group == 1){
				tied2++; 
			}
		}
		else{
			int r0j = remained1 - event1_at_j; 
			int r1j = remained2 - event2_at_j; 

			int total_event_at_j = event1_at_j + event2_at_j; 
			int total_remain = remained1 + remained2; 

			double o = 0.0; double e = 0.0; double v = 0.0; 
			if(total_remain > 1){
				o = event2_at_j; 
				O += o; 
				e = remained2*total_event_at_j/(double)total_remain; 
				E += e; 
				v = remained1*remained2*total_event_at_j*(total_remain - total_event_at_j)/(double)(total_remain * total_remain * (total_remain - 1.0) ); 
				V += v; 
			}

			//std::cout << pos->group << '\t' << pos->event_time << '\t' << pos->event << '\t' << event1_at_j << '\t' << r0j << '\t' << remained1 << '\t' << event2_at_j << '\t' << r1j << '\t' << remained2 << '\t' << total_event_at_j << '\t' << total_remain << '\t' << o << '\t' << e << '\t' << v << std::endl; 

			remained1 -= tied1; 
			remained2 -= tied2; 

			event1_at_j = 0; event2_at_j = 0; 
			if(pos->group == 0) remained1--;
			else if(pos->group == 1) remained2--; 

			tied1 = 0; tied2 = 0; 
		}
	}


	//std::cout << O << '\t' << E << '\t' << V << std::endl; 
	double Z = (O - E)/sqrt(V); 
	
	if(z_score != NULL) *z_score = (-1.0)*Z;
	if(chisquare != NULL) *chisquare = Z*Z; 
	if(pvalue != NULL) *pvalue = 1.0 - erf(fabs(Z)/sqrt(2.0)); 

	//std::cout << *chisquare << '\t' << *pvalue << std::endl; 
	
	return true;
}
