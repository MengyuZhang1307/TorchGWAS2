#include "ReadBGEN.h"

// Run pipeline function
V_bgen run_bgen(const GEMOptions& opt) 
{
    // Step 1: Read covariates
    CovariateReadResult cov_result = read_covariate_data(
        opt.cov_file,
        opt.covariates,
        opt.exposures,
        opt.interactions,
        opt.sampleid_header_name,
        opt.random_slope_header_name,
        opt.delim_cov,
        opt.missing_key
    );

    // Step 2: Read phenotype
    std::set<int> pheno_valid_indices;
    V_string colnames;
    VV_string phenotype_data;

    process_phenotype_file(opt.pheno_file, opt.sampleid_header_name,cov_result.sampleID_list,
                                pheno_valid_indices,
                                colnames,
                                phenotype_data,
                                opt.delim_pheno,
                                opt.missing_key
                        );

    // Step 3: Clean covMap
    clean_covMap_by_invalid_indices( cov_result.sampleID_list, pheno_valid_indices, 
                                        cov_result.covMap);

 
    Bgen bgen;

    // Step 5: Parallel dosage calculation
    bgen.processBgenHeaderBlock(opt.bgen_file);
    bgen.processBgenSampleBlock(bgen, opt.sample_file.c_str(), opt.use_sample_file, cov_result.covMap,
                                    opt.missing_key, cov_result.numSelCol, cov_result.samSize);
    bgen.getPositionOfBgenVariant(bgen, opt.threads, opt.includeVariantFile, opt.do_filters);
    int stream_snps =1;
    
    V_bgen results(opt.threads);

    auto worker = [&](int tid) 
    {
        auto thread_results = calcDosage(opt.bgen_file, opt.stream_snps, tid, bgen); // returns nested result
        results[tid] = std::move(thread_results);
    };

    std::vector<std::thread> pool;
    for (int i = 0; i < opt.threads; ++i) 
    {
        pool.emplace_back(worker, i);
    }

    for (auto& t : pool) t.join();

    return results;
}
