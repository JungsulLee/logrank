#include "logrank.h"


#include <set>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <math.h>
#include <dirent.h>
#include <stdlib.h>

std::string strip(const std::string& str, const std::string& wc)
{
	int length = str.size();
	int i = 0; 
	int first_pos = 0;
	int last_pos = length - 1;
	std::set<char> white_char; 
	for(i = 0; i<wc.size(); i++){
		white_char.insert(wc[i]); 
	}

	for(i = 0; i<length; i++){
		if(white_char.find(str[i]) != white_char.end()) continue; 
		if(white_char.find(str[i]) == white_char.end()){
			first_pos = i;
			break;
		}
	}
	for(i = length - 1; i>=0; i--){
		if(white_char.find(str[i]) != white_char.end()) continue; 
		if(white_char.find(str[i]) == white_char.end()){
			last_pos = i; 
			break; 
		}
	}

	if(first_pos > last_pos){
		first_pos = 0; 
		last_pos = length - 1;
	}

	return str.substr(first_pos, last_pos - first_pos + 1); 
}

std::vector<std::string> split(const std::string& src, const char& delimit)
{
	std::vector<std::string> result; 
	std::string line = strip(src, "\r\n"); 
	std::string str(1, delimit); 
	line += delimit;
	int length = line.size(); 
	int i = 0; 
	int pre_pos = 0; 
	for(i = 0; i<length; i++){
		if(line[i] == delimit || line[i] == '\r' || line[i] == '\n'){
			std::string str = line.substr(pre_pos, i - pre_pos); 
			pre_pos = i + 1; 
			result.push_back(str); 
		}
	}

	return result; 
}

bool get_file_list(const std::string& init_dir, const std::string& filter, std::vector<std::string>* dest)
{
	struct dirent **namelist; 
	int n = 0; 
	n = scandir(init_dir.c_str(),&namelist, NULL, alphasort); 
	if(n <= 0) return false; 

	while(n--){
		std::string file_name = namelist[n]->d_name; 
		if(file_name.size() <= filter.size()){
			free(namelist[n]); 
			continue; 
		}
		if(filter[0] == '*' && filter[filter.size()-1] != '*'){
			std::string new_filter = filter.substr(1,filter.size()-1); 
			if(file_name.substr(file_name.size()-new_filter.size()) == new_filter){
				dest->push_back(file_name); 
			}
		}
		else if(filter[0] == '*' && filter[filter.size()-1] == '*'){
			std::string new_filter = filter.substr(1, filter.size()-2); 
			if(file_name.find(new_filter) != std::string::npos){
				dest->push_back(file_name);
			}
		}
		free(namelist[n]); 
	}
	free(namelist);


	return true;
}


int main(int argc, const char* argv[])
{

	std::vector<std::string> files; 
	get_file_list("./case/", "*.txt", &files); 
	std::sort(files.begin(), files.end()); 

	std::vector<std::string>::const_iterator fpos;
	for(fpos = files.begin(); fpos != files.end(); fpos++){
		std::string file_name = "./case/" + *fpos; 
		std::ifstream ifile(file_name.c_str()); 
		std::string line; 
		double solution_chsq = -1.0; 
		double solution_p = -1.0; 
		std::vector<std::pair<double,bool> > low_group;
		std::vector<std::pair<double,bool> > high_group;
		bool passed_stat_section = false; 
		while(std::getline(ifile, line)){
			if(line == "" || line[0] == '#') continue; 
			if(passed_stat_section == false){
				if(line.find("Comparison of Survival Curves") != std::string::npos) passed_stat_section = true; 
			}
			std::vector<std::string> v = split(line, '\t'); 
			if(passed_stat_section == false){// && v.size() == 3){
				double tvalue = atof(v[0].c_str()); 
				if(fabs(tvalue) < 1e-50) continue; 
				//if(v[1] == ""){
				if(v.size() == 3){
					if(v[1] == ""){
						high_group.push_back(std::make_pair<double,bool>(tvalue, (bool)(v[2] == "1"))); 
					}
					else{
						low_group.push_back(std::make_pair<double,bool>(tvalue, (bool)(v[1] == "1"))); 
					}
				}
				//else if(v[2] == ""){
				else if(v.size() == 2){
					low_group.push_back(std::make_pair<double,bool>(tvalue, (bool)(v[1] == "1"))); 
				}
			}

			if(passed_stat_section == false) continue; 
			/*
				if(v[0] == "INFO"){
					if(v[1] == "MANTEL-COX"){
						if(v[2] == "CHISQUARE") solution_chsq = atof(v[3].c_str()); 
						else if(v[2] == "PVALUE") solution_p = atof(v[3].c_str()); 
					}
				}
			*/
			if(line.find("Gehan-Breslow-Wilcoxon") != std::string::npos) break; 
			if(line.find("Chi square") != std::string::npos) solution_chsq = atof(v[1].c_str()); 
			if(line.find("P value summary") == std::string::npos && line.find("P value") != std::string::npos) solution_p = atof(v[1].c_str()); 
		}

		double p = 1.0; 
		double chsq = 0.0; 
		//std::cout << low_group.size() << '\t' << high_group.size() << std::endl;

		double z = 0.0; 
		logrank(low_group, high_group, &z, &p, &chsq); //&chsq);

		double chsq_diff = fabs(solution_chsq - chsq)/solution_chsq; 
		std::string tag = "SMALL"; 
		if(chsq_diff > 0.01) tag = "BIG"; 

		std::cout << file_name << '\t' << tag ;
		std::cout << "\tChiSq.diff\t" << chsq_diff << '\t' << solution_chsq << '\t' << chsq << '\t' << fabs(solution_chsq - chsq) ;
		std::cout << "\tP.diff\t" << fabs(solution_p - p) << '\t' << solution_p << '\t' << p  << std::endl; 

		//std::cout << chsq << '\t' << p << std::endl; 
	}



	return 0;
}
