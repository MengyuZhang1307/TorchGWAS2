#include "RunPipeline.h"

GEMRunner::GEMRunner(const GEMOptions& user_opt) : opt(user_opt) 
{
    find_genofile_type();
    check_kinship_usage();
    // Step 1:  Read covariate file
    shared_cov_result = read_covariate_data();

    // Step 2 run BGEN metods
    bgen.process_bgen_header_block(opt.geno_add);
    bgen.process_bgen_sample_block(opt.sample_add.c_str(), opt.use_sample_file, 
                                    shared_cov_result.covMap, opt.missing_key, 
                                    shared_cov_result.numSelCol, 
                                    shared_cov_result.samSize);   
    // bgen.get_position_bgen_variant(opt.num_chunks, opt.includeVariantFile,
    //                                          opt.do_filters);
    bgen_sample_id = bgen.sampleID;
    bgen.filterVariants = opt.do_filters;
    is_dup_id = shared_cov_result.cov_is_duplicated;
    // free heavy members
    shared_cov_result.sampleID_list.clear();
    shared_cov_result.covMap.clear();
    shared_cov_result.valid_indices.clear();
    shared_cov_result.samSize = 0;
    shared_cov_result.numSelCol = 0;
}

/**
 * @brief find genotype file type
 * 
 * @return std::string 
 */

void GEMRunner::find_genofile_type()
{
    std::string ext = fs::path(opt.geno_add).extension().string();
        if (ext == ".bgen") 
        {
            genofile_type = "BGEN";
        } 
        else if (ext == ".bed") 
        {
            genofile_type = "BED";
        } 
        else if (ext == ".pgen") 
        {
            genofile_type = "PGEN";
        } 
        else 
        {
            genofile_type = "UNKNOWN";
            std::cerr << "The extextion of genotype should be bgen, pgen or bed \n";
            std::exit(EXIT_FAILURE);
        }
}


/**
 * @brief reade covariate file
 * 
 * @param cov_add
 * @param opt.covariates 
 * @param opt.exposures 
 * @param opt.interactions 
 * @param opt.sampleid_header_name 
 * @param opt.random_slope_header_name 
 * @param opt.cov_delim
 * @param opt.missing_key 
 * @return CovariateReadResult: Structure containing parsed data.
 */
CovariateReadResult GEMRunner::read_covariate_data() 
{
    CovariateReadResult result;
    bool& cov_is_duplicated = result.cov_is_duplicated;
    cov_is_duplicated = false;

    int numExpSelCol = opt.exposures.size();
    int numIntSelCol = opt.interactions.size();

    for (int i = numIntSelCol - 1; i >= 0; --i)
        opt.covariates.insert(opt.covariates.begin(), opt.interactions[i]);
    for (int i = numExpSelCol - 1; i >= 0; --i)
        opt.covariates.insert(opt.covariates.begin(), opt.exposures[i]);

    result.numSelCol = opt.covariates.size() - numExpSelCol - numIntSelCol;
    std::ext::V_int colSelVec(opt.covariates.size());

    std::ifstream fincov(opt.cov_add);
    if (fincov.fail()) 
    {
        throw std::runtime_error("ERROR: Cannot open covariate file");
    }

    std::unordered_map<std::string, int> colNames;
    std::string line, headerName;
    std::getline(fincov, line);
    std::istringstream issHead(line);
    int header_i = 0;
    while (std::getline(issHead, headerName, opt.cov_delim)) 
    {
        headerName.erase(std::remove(headerName.begin(), headerName.end(), '\r'), headerName.end());
        headerName.erase(std::remove(headerName.begin(), headerName.end(), '"'), headerName.end());
        if (colNames.count(headerName)) 
        {
            throw std::runtime_error("ERROR: Duplicate header in covariate file: " + headerName);
        }
        colNames[headerName] = header_i++;
    }

    if (!colNames.count(opt.sampleid_header_name))
        throw std::runtime_error("ERROR: Sample ID column not found");
    int SamIDCol = colNames[opt.sampleid_header_name];

    if (!opt.random_slope_header_name.empty()) 
    {
        if (!colNames.count(opt.random_slope_header_name))
            throw std::runtime_error("ERROR: Random slope column not found");
    }

    for (const auto& h : opt.exposures)
        if (!colNames.count(h)) throw std::runtime_error("ERROR: Exposure column not found: " + h);
    for (const auto& h : opt.interactions)
        if (!colNames.count(h)) throw std::runtime_error("ERROR: Interaction column not found: " + h);
    for (size_t i = 0; i < opt.covariates.size(); ++i) 
    {
        if (!colNames.count(opt.covariates[i]))
        {
            throw std::runtime_error("ERROR: Covariate column not found: " + opt.covariates[i]);
        }
        colSelVec[i] = colNames[opt.covariates[i]];
    }

    int nrows = 0;
    while (std::getline(fincov, line)) ++nrows;

    result.samSize = nrows;
    std::cout << "****************************************************************************\n";
    std::cout << "Number of samples in covariate file: " << nrows << "\n";
    std::cout << "****************************************************************************\n";

    fincov.clear();
    fincov.seekg(0);
    std::getline(fincov, line); // skip header
    
    for (int r = 0; r < result.samSize; ++r) 
    {
        std::getline(fincov, line);
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        std::istringstream iss(line);
        std::string value;
        std::ext::V_string values;
        while (std::getline(iss, value, opt.cov_delim)) values.push_back(value);

        if (values.size() == colNames.size() - 1) values.push_back("");
        if (values.size() != colNames.size())
        {
            throw std::runtime_error("ERROR: Column mismatch at row " + std::to_string(r));
        }
        if (!cov_is_duplicated && result.covMap.count(values[SamIDCol]))
        {
            cov_is_duplicated = true;
        }

        bool has_missing;
        std::ext::V_string entry;
        for (int c : colSelVec) 
        {
            std::string val = values[c];

            if (val == "" || val == opt.missing_key)
            {
                has_missing = true;
            }
            val.erase(std::remove(val.begin(), val.end(), '"'), val.end());
            entry.push_back(val);
        }
        // Return lines number to keep
        if (!has_missing) 
        {
            result.valid_indices.push_back(r);
        }

        result.sampleID_list.push_back(values[SamIDCol]);
        result.covMap[values[SamIDCol]].push_back(entry);
    }

    fincov.close();

    // Check for potential categorical variables by counting unique values
    std::cout << "Checking for potential categorical variables...\n";    
    const int CATEGORICAL_THRESHOLD = 10; // Consider as categorical if <= 10 unique values
    
    // For each covariate column, count unique values
    for (size_t col_idx = 0; col_idx < opt.covariates.size(); ++col_idx) 
    {
        std::unordered_set<std::string> unique_values;
        int non_missing_count = 0;
        bool exceeds_threshold = false;
        
        // Collect all non-missing values for this column with early exit optimization
        for (const auto& sample_entry : result.covMap) 
        {
            for (const auto& measurement : sample_entry.second) 
            {
                if (col_idx < measurement.size()) 
                {
                    const std::string& val = measurement[col_idx];
                    if (val != opt.missing_key && val != "") 
                    {
                        unique_values.insert(val);
                        non_missing_count++;
                        
                        // Early exit if we already exceed the categorical threshold
                        if (unique_values.size() > CATEGORICAL_THRESHOLD) 
                        {
                            exceeds_threshold = true;
                            break;
                        }
                    }
                }
            }
            if (exceeds_threshold) break; // Break outer loop too
        }
        
        int unique_count = unique_values.size();
        
        // Only process and warn if within categorical threshold
        if (unique_count <= CATEGORICAL_THRESHOLD && non_missing_count > 0) 
        {
            std::cout << "WARNING: Covariate '" << opt.covariates[col_idx] 
                      << "' has only " << unique_count << " unique values out of " 
                      << non_missing_count << " non-missing observations.\n";
            std::cout << "         Unique values: ";
            
            // Convert to vector and sort for consistent output
            std::vector<std::string> sorted_values(unique_values.begin(), unique_values.end());
            std::sort(sorted_values.begin(), sorted_values.end(), 
                     [](const std::string& a, const std::string& b) {
                         // Try to sort numerically if possible, otherwise lexicographically
                         try {
                             double num_a = std::stod(a);
                             double num_b = std::stod(b);
                             return num_a < num_b;
                         } catch (...) {
                             return a < b;
                         }
                     });
            
            for (size_t i = 0; i < sorted_values.size(); ++i) 
            {
                std::cout << sorted_values[i];
                if (i < sorted_values.size() - 1) std::cout << ", ";
            }
            std::cout << "\n";
            
            if (unique_count == 2) 
            {
                std::cout << "         This appears to be a binary categorical variable.\n";
                std::cout << "         RECOMMENDATION: Convert to one-hot encoding using (k-1) dummy variables.\n";
                std::cout << "         For 2 categories, create 1 dummy variable (0/1 encoding).\n";
            }
            else 
            {
                std::cout << "         This may be a categorical variable with multiple levels.\n";
                std::cout << "         RECOMMENDATION: Convert to one-hot encoding if nominal categorical.\n";
                std::cout << "         For " << unique_count << " categories, create " << (unique_count - 1) << " dummy variables.\n";
            }           
            std::cout << "\n";
        }
    }
    std::cout << "****************************************************************************\n";

    return result;
}

/**
 * @brief Check if user provided the kinship file
 * 
 */

void GEMRunner::check_kinship_usage() 
{
    kin_flag = !opt.kin_add.empty();  // sets true if user provided kin_path
}

/**
 * @brief Function to call null model from python
 * 
 */
void GEMRunner::run_fit_nullmodel() 
{
    NullModel model(opt);
    model.fit_nullmodel(kin_flag,
                        bgen_sample_id,
                        is_dup_id,                       
                        c2_values // Pass C2 values by reference
    );
}
