/*  GEM : Gene-Environment interaction analysis for Millions of samples
 *  Copyright (C) 2018-2025  Liang Hong, Han Chen, Duy Pham, Cong Pan, Samaneh Salehi Nasab
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

 /*
   software description and change log

   version 1: GWAS of snp one by one
   version 2: streaming multi-snps and GWAS analysis of multi-snps simultaneously
   version bgen: multi-snps and bgen
   version Log:  logistic regression for binary phenodata
   7/27/18:  input column header name in pheno file not index number
   8/01/18:  sample ID matching
   8/30/18:  initial sample Size is not necessary equal
   8/30/18:  hash table of genoUnMatchID for significant unmatching numbers;
   12/10/18: aritrary stream_snps implemented, omp parallel with > 20 covariates (HC)
   2/7/19:   speed up reading genotype data without looking up in genoUnMatchID

   To-Do List:
   1. OOP
 */

#include "FitNullModel.h" 

NullModel::NullModel(const GEMOptions& user_opt): opt(user_opt){}
 
void NullModel::process_gmmat(const std::string kin_add, 
                            const std::string cov_add, const char kin_delim, 
                            const double kin_diag, const char cov_delim, 
                            const std::string &sampleid_header_name, const std::ext::V_string &cov_headers, 
                            std::ext::V_string &bgen_sample_id, const std::string missing_key, 
                            const std::ext::V_string& pheno_column_names,
                            const std::ext::VV_string& phenotype_data,
                            std::vector<std::set<int>>& pheno_valid_indices,
                            int num_threads,
                            const std::ext::FitNull_f& fitNullModel2,
                            const std::ext::V_string& covariates,
                            const std::string& random_slope_header_name,
                            const std::string& output,
                            std::ext::V_double& c2_out)
{
    int col = 0;
    int num_columns = pheno_column_names.size();
    std::ext::VV_double output_matrix;
    size_t res_size;
    bool sampleid_filled = false;    
    std::ext::V_string id_include;
    const int pheno_columns = num_columns - 2;
    const int block_size = (pheno_columns + num_threads - 1) / num_threads;  // ceil division

    std::vector<std::thread> threads;
    //glmmkin_residuals is struct return type by GMMAT
    std::ext::VV_string id_include_vec(pheno_columns);
    std::map<std::string, glmmkin_residuals> residual_map;
    std::vector<std::map<std::string, glmmkin_residuals>> thread_local_maps(num_threads);

    Cov cov;
    cov.m_sam_id_hdr = sampleid_header_name;
    cov.m_v_hdrs = cov_headers;
    cov.m_data_frame.m_geno_ids = bgen_sample_id;
    cov.m_data_frame.m_missing_key = missing_key;
    cov.set_path(cov_add);


    auto path = cov.get_path();
    cov.read_file(path, cov_delim, cov.m_v_hdrs);
    //Match genofile sample IDs
    fmt::println("Number of observation in covariate file before matching IDs with genotype IDS is: {}", cov.m_data_frame.n_rows());
    //Remove lines with missing data from cov data based on missing value in cov and missing sampleID in genotype file
    cov.m_data_frame.match_genoids(cov.m_sam_id_hdr, cov.m_v_hdrs);
    fmt::println("Number of observation in covariate file after matching IDs with genotype IDs is: {}", cov.m_data_frame.n_rows());
    fmt::println("****************************************************************************");
    //Map cov sample ids to int to be used as matrix indices

    for (int t = 0; t < num_threads; ++t) 
    {
        int start_col = t * block_size;
        int end_col = std::min(start_col + block_size, pheno_columns);

        if (start_col >= end_col) break;  // no more work

        threads.emplace_back(
            [cov_copy = cov, &kin_add, &kin_delim, &kin_diag,
            &cov_delim, &bgen_sample_id, 
            &missing_key, start_col, t, end_col, &fitNullModel2, 
            &phenotype_data, &pheno_valid_indices, &covariates, &random_slope_header_name, 
            &pheno_column_names, &thread_local_maps] () mutable
            {
                auto& local_map = thread_local_maps[t];
                
                for (int this_col = start_col; this_col < end_col; ++this_col) 
                {
                    GMMAT gmmat;
                    
                    glmmkin_residuals residuals = gmmat.glmmkin_init(cov_copy,
                        kin_add, kin_delim, kin_diag,
                        cov_delim, bgen_sample_id, 
                        missing_key, fitNullModel2, phenotype_data[this_col], 
                        pheno_valid_indices[this_col], 
                        covariates, random_slope_header_name, "", "REML", "AI",
                        500, 1e-5, 1e-5, 1e+5, 10
                    );

                    local_map[pheno_column_names[this_col + 2]] = std::move(residuals);
                }
            }
        );
    }

    for (auto& th : threads) 
    {
        th.join();
    }

    for (const auto& local_map : thread_local_maps)
    {
        for (const auto& [key, value] : local_map)
        {
            residual_map[key] = value;
        }
    }

    // if (id_include.size() <= 0)
    // {
        // id_include = residual_map[pheno_column_names[2]].id_include;
    // }
    
    std::ext::V_double c2;
    for (int col = 0; col < pheno_columns; ++col) 
    {
        auto& residual = residual_map[pheno_column_names[col + 2]].scaled_residuals_c1;
        id_include_vec.emplace_back(residual_map[pheno_column_names[col + 2]].id_include);
        output_matrix.emplace_back(std::move(residual));
        c2.push_back(residual_map[pheno_column_names[col + 2]].c2);
    }

    // Return C2 values
    c2_out = c2;

    print_res(output, pheno_column_names, c2, bgen_sample_id, id_include_vec, output_matrix);
}

 
void NullModel::fit_nullmodel(bool kin_flag,
        std::ext::V_string& bgen_sample_id,
        bool dup_id,
        std::ext::V_string& shared_pheno_colnames,
        std::ext::VV_string& shared_phenotype_data,
        std::vector<std::set<int>>& pheno_valid_indices,
        std::ext::V_double& c2_out)
{
    if (kin_flag || dup_id)
    {
        auto start_time_gmmat = std::chrono::high_resolution_clock::now();
        vector <string> cov_headers(opt.covariates);
        cov_headers.insert(cov_headers.begin(), opt.sampleid_header_name);
        if(std::find(cov_headers.begin(), cov_headers.end(), opt.random_slope_header_name) == cov_headers.end())
        {
            cov_headers.insert(cov_headers.end(), opt.random_slope_header_name);
        }           
     
        process_gmmat(opt.kin_add, opt.cov_add, opt.kin_delim,
                        opt.kin_diag, opt.cov_delim, opt.sampleid_header_name, 
                        cov_headers, bgen_sample_id, opt.missing_key,
                        shared_pheno_colnames, shared_phenotype_data, 
                        pheno_valid_indices, opt.threads, fitNullModel2, 
                        opt.covariates, opt.random_slope_header_name, 
                        opt.outfile, c2_out);
        cout << "\nEnd of association test\n";
        cout << "****************************************************************************\n";
        cout << "calculating the duration of association test...\n";
        auto end_time_gmmat = std::chrono::high_resolution_clock::now();
        printExecutionTime(start_time_gmmat, end_time_gmmat);
        cout << std::flush;
    }  
    else
    {
        std::cerr << "Please make sure you have repetaed measure data or define a kinship\n";
        std::exit(EXIT_FAILURE);
        
    } 
}
 

/**
 * @brief perfom centering
 * 
 * @param center 
 * @param scale 
 * @param samSize 
 * @param numSelCol 
 * @param covdata 
 * @param covdata_ret 
 */
void center(int center, int scale, int samSize, int numSelCol, std::ext::V_double covdata, std::ext::V_double* covdata_ret) 
{
    std::ext::V_double tmp1(samSize, 1);
    double* tmpMean = new double[numSelCol + 1];
    std::ext::V_double tmpSD(numSelCol + 1);
    if (center) 
    {
        matmatprod(&tmp1[0], &covdata[0], tmpMean, 1, samSize, numSelCol + 1);
        if (!scale) {
            cout << "Centering without rescaling..." << endl;
            for (int i = 1; i < numSelCol + 1; i++) {
                tmpMean[i] /= double(samSize * 1.0);
                tmpSD[i] = 1.0;
            }
        }
        else {
            cout << "Centering and rescaling..." << endl;
            for (int i = 1; i < numSelCol + 1; i++) {
                tmpMean[i] /= double(samSize * 1.0);
            }
            for (int i = 0; i < samSize; i++) {
                for (int j = 1; j < numSelCol + 1; j++) {
                    tmpSD[j] += pow(covdata[i * (numSelCol + 1) + j] - tmpMean[j], 2.0);
                }
            }
            for (int i = 1; i < numSelCol + 1; i++) {
                tmpSD[i] = sqrt(tmpSD[i] / double(samSize * 1.0 - 1.0));
            }
        }

        for (int i = 0; i < samSize; i++) {
            for (int j = 1; j < numSelCol + 1; j++) {
                covdata[i * (numSelCol + 1) + j] = (covdata[i * (numSelCol + 1) + j] - tmpMean[j]) / tmpSD[j];
            }
        }

    }
    else 
    {
        if (scale) {
            cout << "Scaling ALL exposures and covariates..." << endl;
            matmatprod(&tmp1[0], &covdata[0], tmpMean, 1, samSize, numSelCol + 1);
            for (int i = 1; i < numSelCol + 1; i++) {
                tmpMean[i] /= double(samSize * 1.0);
            }
            for (int i = 0; i < samSize; i++) {
                for (int j = 1; j < numSelCol + 1; j++) {
                    tmpSD[j] += pow(covdata[i * (numSelCol + 1) + j] - tmpMean[j], 2.0);
                }
            }
            for (int i = 1; i < numSelCol + 1; i++) {
                tmpSD[i] = sqrt(tmpSD[i] / double(samSize * 1.0 - 1.0));
            }
            for (int i = 0; i < samSize; i++) {
                for (int j = 1; j < numSelCol + 1; j++) {
                    covdata[i * (numSelCol + 1) + j] /= tmpSD[j];
                }
            }
        }
    }
    delete[] tmpMean;
    *covdata_ret = covdata;
} 
     

void printCovVarMat(int numCovs, std::ext::V_string covNames, double* covVarMat, double* beta, int phenoType, int samSize) 
{
    covNames.insert(covNames.begin(), "Intercept");
    boost::math::chi_squared chisq_dist_M(1);

    cout << "\nCoefficients: \n";
    cout << boost::format("%-26s %-17s %-22s %-19s %-15s\n") % "" % "Estimate" % "Std. Error" % "Z-value" % "P-value";
    for (int i = 0; i < numCovs; i++) 
    {
        double stdError = sqrt(covVarMat[i * numCovs + i]);
        double zvalue = beta[i] / stdError;
        double pr = (isnan(zvalue)) ? NAN : boost::math::cdf(complement(chisq_dist_M, (beta[i] * beta[i]) / covVarMat[i * numCovs + i]));
        cout << boost::format("%+15s %19.6e %19.6e %19.6e %19.6e\n") % covNames[i] % beta[i] % stdError % zvalue % pr;
    }

    cout << "\nVariance-Covariance Matrix: \n";
    cout << boost::format("%+35s") % covNames[0];
    for (int i = 1; i < numCovs; i++) {
        cout << boost::format("%+20s") % covNames[i];
    }
    cout << "\n";
    for (int i = 0; i < numCovs; i++) {
        cout << boost::format("%+15s") % covNames[i];
        for (int j = 0; j < numCovs; j++) {
            cout << boost::format("%20.6e") % covVarMat[j * numCovs + i];
        }
        cout << "\n";
    }
    cout << "\n";
}

 
void fitNullModel2(int samSize, int numSelCol, int phenoType, double epsilon, 
                int robust, std::ext::V_string covariates, std::ext::V_double phenodata, 
                std::ext::V_double covdata, std::ext::V_double* XinvXTX_ret, std::ext::V_double* miu_ret, 
                std::ext::V_double* resid_ret, double* sigma2_ret, std::ext::V_double& beta_ret,
                std::ext::V_double& Xbeta_ret)
{
    double* phenoY = &phenodata[0];
    double* covX = &covdata[0];
    vector <double> residvec(samSize);
    // for logistic regression
    vector <double> miu(samSize);
    int Check = 1; // convergence condition of beta^(i+1) - beta^(i)
    int iter = 1;


    cout << "Precalculations and fitting null model..." << endl;
    //auto start_time = std::chrono::high_resolution_clock::now();
    // transpose(X) * X
    double* XTransX = new double[(numSelCol + 1) * (numSelCol + 1)];
    matTmatprod(covX, covX, XTransX, samSize, numSelCol + 1, numSelCol + 1);
    // invert (XTransX)
    matInv(XTransX, numSelCol + 1);
    // transpose(X) * Y
    double* XTransY = new double[(numSelCol + 1)];
    matTvecprod(covX, phenoY, XTransY, samSize, numSelCol + 1);
    // beta = invert(XTransX) * XTransY
    double* beta = new double[(numSelCol + 1)];
    beta_ret.resize(numSelCol + 1);
    matvecprod(XTransX, XTransY, beta, numSelCol + 1, numSelCol + 1);

    // logistic regression
    while ((phenoType == 1) && (Check != (numSelCol + 1))) 
    {
        iter++;
        // X * beta
        double* XbetaFL = new double[samSize];
        matvecprod(covX, beta, XbetaFL, samSize, numSelCol + 1);
        double* Yip1 = new double[samSize];
        // W * X and W * Y
        double* WX = new double[samSize * (numSelCol + 1)];
        double* WYip1 = new double[samSize];
        for (int i = 0; i < samSize; i++) {
            miu[i] = exp(XbetaFL[i]) / (1.0 + exp(XbetaFL[i]));
            Yip1[i] = XbetaFL[i] + (phenoY[i] - miu[i]) / (miu[i] * (1 - miu[i]));
            WYip1[i] = miu[i] * (1 - miu[i]) * Yip1[i];
            for (int j = 0; j < numSelCol + 1; j++) {
                WX[i * (numSelCol + 1) + j] = miu[i] * (1 - miu[i]) * covX[i * (numSelCol + 1) + j];
            }
        }
        // transpose(X) * WX
        matTmatprod(covX, WX, XTransX, samSize, numSelCol + 1, numSelCol + 1);
        // invert (XTransX)
        matInv(XTransX, numSelCol + 1);
        // transpose(X) * WYip1
        matTvecprod(covX, WYip1, XTransY, samSize, numSelCol + 1);
        // beta = invert(XTransX) * XTransY
        double* betaT = new double[(numSelCol + 1)];
        matvecprod(XTransX, XTransY, betaT, numSelCol + 1, numSelCol + 1);
        Check = 0;
        for (int i = 0; i < numSelCol + 1; i++) {
            if (std::abs(betaT[i] - beta[i]) <= epsilon) Check++;
            beta[i] = betaT[i];
        }

        delete[] Yip1;
        delete[] WYip1;
        delete[] WX;
        delete[] XbetaFL;
        delete[] betaT;
    }

    // X * beta
    double* Xbeta = new double[samSize];
    Xbeta_ret.resize(samSize);
    matvecprod(covX, beta, Xbeta, samSize, numSelCol + 1);

    for(int i{0}; i < samSize; ++i)
    {
        Xbeta_ret[i] = Xbeta[i];
    }

    // X*[invert (XTransX)]
    std::ext::V_double XinvXTXvec(samSize * (numSelCol + 1));
    double* XinvXTX = &XinvXTXvec[0];
    if (phenoType == 1) {
        double* WX = new double[samSize * (numSelCol + 1)];
        for (int i = 0; i < samSize; i++) {
            miu[i] = exp(Xbeta[i]) / (1.0 + exp(Xbeta[i]));
            Xbeta[i] = miu[i];
            for (int j = 0; j < numSelCol + 1; j++) {
                WX[i * (numSelCol + 1) + j] = miu[i] * (1.0 - miu[i]) * covX[i * (numSelCol + 1) + j];
            }
        }
        // transpose(X) * WX
        matTmatprod(covX, WX, XTransX, samSize, numSelCol + 1, numSelCol + 1);
        // invert (XTransX)
        matInv(XTransX, numSelCol + 1);
        matmatprod(WX, XTransX, XinvXTX, samSize, numSelCol + 1, numSelCol + 1);
        delete[] WX;

        cout << "Logistic regression reaches convergence after " << iter << " steps...\n";
    }
    else {
        matmatprod(covX, XTransX, XinvXTX, samSize, numSelCol + 1, numSelCol + 1);
    }


    // residual = Y - X * beta
    double sigma2 = 0;
    for (int i = 0; i < samSize; i++) {
        residvec[i] = phenoY[i] - Xbeta[i];
        sigma2 += residvec[i] * residvec[i];
    }
    double* resid = &residvec[0];

    // sqr(sigma) = transpose(resid)*resid/[samSize-(numSelCol+1)]
    sigma2 = sigma2 / (samSize - (numSelCol + 1));
    if (phenoType == 1) sigma2 = 1.0;


    std::ext::V_double XR2vec;
    if (!robust) 
    {
        for (int i = 0; i < (numSelCol + 1) * (numSelCol + 1); i++) {
            XTransX[i] = XTransX[i] * sigma2;
        }
        printCovVarMat(numSelCol + 1, covariates, XTransX, beta, phenoType, samSize);
    }
    else {
        std::ext::V_double XR2vec = covdata;
        for (int i = 0; i < samSize; i++) {
            for (int j = 0; j < numSelCol + 1; j++) {
                XR2vec[i * (numSelCol + 1) + j] = XR2vec[i * (numSelCol + 1) + j] * resid[i] * resid[i];
            }
        }

        double* XR2 = &XR2vec[0];
        double* XR2tX = new double[(numSelCol + 1) * (numSelCol + 1)];
        matTmatprod(XR2, covX, XR2tX, samSize, numSelCol + 1, numSelCol + 1);
        double* XTransXtXR2tX = new double[(numSelCol + 1) * (numSelCol + 1)];
        matmatTprod(XR2tX, XTransX, XTransXtXR2tX, numSelCol + 1, numSelCol + 1, numSelCol + 1);
        double* XTransXR2 = new double[(numSelCol + 1) * (numSelCol + 1)];
        matmatTprod(XTransX, XTransXtXR2tX, XTransXR2, numSelCol + 1, numSelCol + 1, numSelCol + 1);
        printCovVarMat(numSelCol + 1, covariates, XTransXR2, beta, phenoType, samSize);
        delete[] XR2tX;
        delete[] XTransXtXR2tX;
        delete[] XTransXR2;
    }

    //filling beta_ret (alpha) and Xbeta(eta),  needed for calculations in GMMAT 
    for(size_t i{0}; i < numSelCol + 1; ++i)
    {
        beta_ret[i] = beta[i];
    } 
    
    delete[] XTransX;
    XTransX = nullptr; 
    delete[] XTransY;
    XTransY = nullptr; 
    delete[] beta;
    beta = nullptr; 
    delete[] Xbeta;
    Xbeta = nullptr; 
    if (phenoType == 1)
    {
        *miu_ret = miu;
    }
    else
    {
        *miu_ret = Xbeta_ret;
    }
    
    *sigma2_ret = sigma2;
    *resid_ret = residvec;
    *XinvXTX_ret = XinvXTXvec;
}



void NullModel::print_res(
    std::string output,
    std::ext::V_string const& pheno_column_names,
    std::ext::V_double const& c2,
    std::ext::V_string const& bgen_sample_id, 
    std::ext::VV_string const& id_include,
    std::ext::VV_double const& output_matrix)
{
    std::ofstream out(output);

    // Header line: phenotype names
    for (size_t i = 0; i < pheno_column_names.size(); ++i) 
    {
        out << pheno_column_names[i];
        if (i != pheno_column_names.size() - 1) out << '\t';
    }
    out << '\n';

    // Correlation row
    out << '#' << '\t' << '#';
    for (auto cor : c2) 
    {
        out << '\t' << cor;
    }
    out << '\n';

    // --- build fast lookup for id_include --- to accomodate each ID seperately
    std::vector<std::unordered_map<std::string, size_t>> id_lookup(id_include.size());
    for (size_t ph = 0; ph < id_include.size(); ++ph) 
    {
        for (size_t idx = 0; idx < id_include[ph].size(); ++idx) 
        {
            id_lookup[ph][id_include[ph][idx]] = idx;
        }
    }

    // Main loop over all samples
    for (size_t row = 0; row < bgen_sample_id.size(); ++row) 
    {
        out << bgen_sample_id[row] << '\t' << bgen_sample_id[row];

        // check each phenotype
        for (size_t ph = 0; ph < pheno_column_names.size(); ++ph) 
        {
            auto it = id_lookup[ph].find(bgen_sample_id[row]);
            if (it != id_lookup[ph].end()) 
            {
                // found → use matching value
                out << '\t' << output_matrix[ph][it->second];
            } 
            else 
            {
                // not found → missing, write zero
                out << '\t' << 0;
            }
        }
        out << '\n';
    }
}

