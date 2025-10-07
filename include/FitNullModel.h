#pragma once
#include "GMMAT.h"
#include "Cov.h"
#include "Declars.h"
#include "FileIO.h" 
#include "GEMConfig.h" 
#include "MatrixUtils.h"
#include "ReadBGEN.h"
#include "TimeUtils.h"

void center(int center, int scale, int samSize, int numSelCol, std::ext::V_double covdata, std::ext::V_double* covdata_ret);
void fitNullModel2(int samSize, int numSelCol, int phenoType, double epsilon, 
                    int robust, std::ext::V_string covSelHeadersName, std::ext::V_double phenodata, 
                    std::ext::V_double covdata, std::ext::V_double* XinvXTX_ret, std::ext::V_double* miu_ret, 
                    std::ext::V_double* resid_ret, double* sigma2_ret, std::ext::V_double& beta_ret,
                    std::ext::V_double& Xbeta_ret);

void printCovVarMat(int numCovs, std::ext::V_string covNames, double* covVarMat, double* beta, int phenoType, int samSize);

class NullModel 
{
    private:
        GEMOptions opt;
        
    public:
        explicit NullModel(GEMOptions const& user_opt);
        void fit_nullmodel(bool kin_flag,
            std::ext::V_string& bgen_sample_id,
            bool dup_id,
            std::ext::V_string& shared_colnames,
            std::ext::VV_string& shared_phenotype_data,
            std::vector<std::set<int>>& pheno_valid_indices,
            std::ext::V_double& c2_out);
        void process_gmmat(const std::string kin_add, 
                            const std::string cov_add, const char kin_delim, 
                            const double kin_diag, const char cov_delim, 
                            const std::string &sampleid_header_name, const std::ext::V_string &cov_headers, 
                            std::ext::V_string &bgen_sample_id, const std::string missing_key,
                            std::ext::V_string const& pheno_column_names,
                            std::ext::VV_string const& phenotype_data,
                            std::vector<std::set<int>>& pheno_valid_indices,
                            int num_threads,
                            std::ext::FitNull_f const& fitNullModel2,
                            std::ext::V_string const& covariates,
                            std::string const& random_slope_header_name,
                            std::string const& output,
                            std::ext::V_double& c2_out);
        void print_res(std::string output, std::ext::V_string const& column_names,
                        std::ext::V_double const& c2, 
                        std::ext::V_string const&bgen_sample_id,
                        std::ext::VV_string const& id_include, 
                        std::ext::VV_double const& output_matrix);


};

