#pragma once
#include<iostream>
#include<sstream>
#include<string>
#include<cstdint>
#include <limits>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<set> 
#include<fstream>
#include<thread>
#include<unordered_map>
#include<unordered_set>


using namespace std;
using uint = unsigned int;
using uchar = unsigned char;
using V_string = std::vector<std::string>; 
using VV_string = std::vector<V_string>; 
using UMap_str_VV_string = std::unordered_map<std::string, VV_string>;
using V_bgen = std::vector<std::vector<std::vector<std::vector<double>>>>;
typedef std::numeric_limits<double> dbl;

namespace extTypes 
{
    using StringSet = std::set<std::string>;
    using VVString = std::vector<std::vector<uint>>;
    using SampleUnmap = std::unordered_map<std::string, std::vector<std::vector<std::string>>>;

}


struct GEMOptions 
{
    std::string phenoFile;
    std::string covFile;
    std::string bgenFile;
    std::string sampleFile;
    std::string sampleIDHeaderName;
    std::string randomSlopeHeaderName;
    std::string missingKey;
    std::string outFile;
    int threads;
    char delim_pheno;
    char delim_cov;
    std::vector<std::string> covariates;
    std::vector<std::string> exposures;
    std::vector<std::string> interactions;
};

struct CovariateReadResult {
    std::vector<std::string> sampleID_list;
    UMap_str_VV_string covMap;
    int samSize;
    bool cov_is_duplicated;
    int numSelCol;
};

char resolve_delim(const std::string& s);
GEMOptions getOptions(int argc, char* argv[]);

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
