#pragma once
#include "GMMAT.h"
#include "Declars.h"
#include "ParallelFileReader.h" 
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
        explicit NullModel(const GEMOptions& user_opt);
        void fit_nullmodel(bool kin_flag,
            std::ext::V_string& bgen_sample_id,
            CovariateReadResult& shared_cov_result,
            std::set<int>& shared_pheno_valid_indices,
            std::ext::V_string& shared_colnames,
            std::ext::VV_string& shared_phenotype_data,
            std::ext::V_double& c2_out);
        void process_gmmat(const std::ext::V_string& column_names,
            const std::ext::VV_string& phenotype_data,
            int num_threads,
            SparseInverse& sp,
            const std::ext::FitNull_f& fitNullModel2,
            const std::ext::V_string& covariates,
            const std::string& random_slope_header_name,
            const std::string& output,
            std::ext::V_double& c2_out);
        void print_res(std::string output, std::ext::V_string const& column_names,
            std::ext::V_double const& c2, std::ext::V_string const& id_include,
            std::ext::VV_double const& output_matrix);


};

