#pragma once
#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <sstream>
#include <string>
#include <cstdint>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <set> 
#include <fstream>
#include <thread>
#include <unordered_map>
#include <unordered_set>


// using namespace std;
using uint = unsigned int;
using uchar = unsigned char;
using V_string = std::vector<std::string>; 
using VV_string = std::vector<V_string>; 
using UMap_str_VV_string = std::unordered_map<std::string, VV_string>;
using VV_float = std::vector<std::vector<float>>;
using VVV_float = std::vector<VV_float>;
using V_bgen = std::vector<std::vector<std::vector<std::vector<float>>>>;
using Set_string = std::set<std::string>;
typedef std::numeric_limits<double> dbl;


struct GEMOptions 
{
    std::string pheno_file;
    std::string cov_file;
    char delim_pheno;
    char delim_cov;
    std::string bgen_file;
    std::string sample_file;
    bool do_filters = false;
    bool use_sample_file =false;
    std::string includeVariantFile;
    int stream_snps = 1;
    std::string sampleid_header_name;
    std::string random_slope_header_name;
    std::string missing_key = "NA";
    std::vector<std::string> covariates;
    std::vector<std::string> exposures;
    std::vector<std::string> interactions;
    int threads;
    std::string out_file;
};

struct CovariateReadResult {
    std::vector<std::string> sampleID_list;
    UMap_str_VV_string covMap;
    int samSize;
    bool cov_is_duplicated;
    int numSelCol;
};

char resolve_delim(const std::string& s);
// Read options from argv
GEMOptions get_options(int argc, const char* argv[]);

CovariateReadResult read_covariate_data(
    const std::string& cov_file,
    std::vector<std::string> covSelHeadersName,
    const std::vector<std::string>& expCovSelHeadersName,
    const std::vector<std::string>& intCovSelHeadersName,
    const std::string& sampleIDHeaderName,
    const std::string& randomSlopeHeaderName,
    char delim_cov = ',',
    const std::string missing_key = "NA") ;

void process_phenotype_file(
    std::string const& filename, 
    std::string sample_id_hdr,
    V_string sampleID_list,
    std::set<int> &valid_indices,
    V_string &column_names,
    VV_string &phenotype_data,
    char delim =  ',',
    std::string missing_key= "NA"
    ); 

// void update_covMap(UMap_str_VV_string &covMap, std::set<int> pheno_valid_indices);

void clean_covMap_by_invalid_indices(
    std::vector<std::string> const& sampleID_list,
    std::set<int> const& pheno_valid_indices,
    UMap_str_VV_string& covMap
);

void write_bgen_result_to_file(const V_bgen& results, const std::string& filename,
                                char delimiter = ',');
