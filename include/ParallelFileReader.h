#pragma once
#include "ReadFiles.h"
/**
 * @brief Struct to hold the result of reading covariate data.
 * 
 */
struct GEMOptions 
{
    std::string pheno_file;
    std::string cov_file;
    char delim_pheno = ',';
    char delim_cov = ',';
    std::string geno_file;
    std::string sample_file = "";
    bool do_filters = false;
    bool use_sample_file = false;
    std::string includeVariantFile = "";
    int stream_snps = 1;
    std::string sampleid_header_name;
    std::string random_slope_header_name = "";
    std::ext::V_string covariates;
    std::ext::V_string exposures;
    std::ext::V_string interactions;
    std::string missing_key = "NA";
    std::string kin_path;
    char delim_k = ',';
    double kin_diag = 1;
    int threads;
    int num_chunks = 0;
    std::string outfile = "Bgen_dosage.out";
    // Determine type of genotype file
    GEMOptions(); 
};

struct CovariateReadResult
 {
    std::ext::V_string sampleID_list;
    std::ext::UMap_str_VV_string covMap;
    int samSize;
    bool cov_is_duplicated;
    int numSelCol;
    /**
     * @brief Indices of lines to keep.
     */
    std::ext::V_int valid_indices;
};

void write_bgen_result_to_file(const std::ext::V_bgen& results, const std::string& filename,
                                char delimiter = ',');

void output_file_generator(int threads, std::string outfile);
