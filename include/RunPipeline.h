#pragma once
#include "ReadBGEN.h"
#include "FitNullModel.h" 
#include "GEMConfig.h" 

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

class GEMRunner 
{
private:
    // Internal helpers
    CovariateReadResult read_covariate_data() ;
    // void process_phenotype_file();     
    void clean_covMap_by_invalid_indices();
    
    public:
    // Options passed by user
    GEMOptions opt;
    Bgen bgen;
    std::ext::V_string bgen_sample_id;
    // Covariates and phenotype data
    CovariateReadResult shared_cov_result;
    bool kin_flag = false;
    bool is_dup_id = false;
    std::string genofile_type; //To be filled by geno_file_type
    std::ext::V_double c2_values;
    // Constructor
    explicit GEMRunner(const GEMOptions& user_opt);
    void find_genofile_type();
    void run_fit_nullmodel(); 
    void check_kinship_usage();
};
