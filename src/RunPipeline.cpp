#include "RunPipeline.h"

GEMRunner::GEMRunner(const GEMOptions& user_opt) : opt(user_opt) 
{
    find_genofile_type();
    check_kinship_usage();
    // Step 1:  Read covariate file
    shared_cov_result = read_covariate_data();

    // Step 2: Read phenotype file
    process_phenotype_file();

    // Step 3: Clean covariate map based on valid phenotype samples
    clean_covMap_by_invalid_indices();
    // Step 4 run BGEN metods
    bgen.process_bgen_header_block(opt.geno_file);
    bgen.process_bgen_sample_block(opt.sample_file.c_str(), opt.use_sample_file, 
                                    shared_cov_result.covMap, opt.missing_key, 
                                    shared_cov_result.numSelCol, 
                                    shared_cov_result.samSize);   
    // bgen.get_position_bgen_variant(opt.num_chunks, opt.includeVariantFile,
    //                                          opt.do_filters);
    bgen_sample_id = bgen.sampleID;
    bgen.filterVariants = opt.do_filters;
    dup_id = shared_cov_result.cov_is_duplicated;
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
    std::string ext = fs::path(opt.geno_file).extension().string();
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
 * @param cov_file 
 * @param opt.covariates 
 * @param opt.exposures 
 * @param opt.interactions 
 * @param opt.sampleid_header_name 
 * @param opt.random_slope_header_name 
 * @param opt.delim_cov 
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

    std::ifstream fincov(opt.cov_file);
    if (fincov.fail()) 
    {
        throw std::runtime_error("ERROR: Cannot open covariate file");
    }

    std::unordered_map<std::string, int> colNames;
    std::string line, headerName;
    std::getline(fincov, line);
    std::istringstream issHead(line);
    int header_i = 0;
    while (std::getline(issHead, headerName, opt.delim_cov)) 
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
        while (std::getline(iss, value, opt.delim_cov)) values.push_back(value);

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
    std::cout << "\n****************************************************************************\n";
    std::cout << "Checking for potential categorical variables...\n";
    std::cout << "****************************************************************************\n";
    
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
 * @brief Read phenotype file
 * 
 * @param opt.pheno_file 
 * @param opt.sampleid_header_name 
 * @param shared_cov_result.sampleID_list 
 * @param valid_indices 
 * @param shared_colnames 
 * @param shared_phenotype_data 
 * @param opt.delim_pheno 
 * @param opt.missing_key 
 */

void GEMRunner::process_phenotype_file() 
{
    std::unordered_set<std::string> seen;
    std::ifstream file(opt.pheno_file);

    if (!file.is_open()) 
    {
        std::cerr << "Error opening file: " << opt.pheno_file << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Read header (column names)
    std::string header_line;
    std::getline(file, header_line);
    std::stringstream ss(header_line);
    std::string col_name;
    int row_indx = 0;
    int col_indx = 0;
    int hdr_id_indx;
    
    while (std::getline(ss, col_name, opt.delim_pheno)) 
    {
        if (seen.insert(col_name).second)
        {
            shared_colnames.push_back(col_name);
            if(col_name == opt.sampleid_header_name)
            {
                hdr_id_indx = col_indx;
            }
        }
        else
        {
            std::cerr << "ERROR: there are repeated columns'name in the phenotype file please check your file.\n";
            exit(EXIT_FAILURE); 
        }
        ++col_indx;
    }
    
    int num_columns = shared_colnames.size();
    std::cout << "Total columns: " << num_columns << "\n"; 
    std::cout << "****************************************************************************\n";
    shared_phenotype_data.resize(num_columns - 1);
    // The first two cols are FID and IID
    if(num_columns < 3)
    {
        std::cerr << "ERROR: number of columns in phenotype file at least should be 3.\n";
        exit(EXIT_FAILURE);
    }

    std::string line;
    while(getline(file, line))
    {
        std::stringstream ss(line);
        std::string value;
        std::ext::V_string values;
        while(getline(ss, value, opt.delim_pheno))
        {
            values.push_back(value);
        }
                    
        if (!line.empty() && line.back() == opt.delim_pheno) 
        {
            values.push_back(opt.missing_key);
        }

        if (values.size() != num_columns) 
        {
            std::cerr << "ERROR: expect: " << num_columns << " columns at row: " << row_indx + 1 << " , while there is: " << values.size() << " columns."<< '\n';
            std::cerr << "If delimiter is space check for extra spaces in line \n";
            exit(EXIT_FAILURE);
        }

        if (row_indx >= shared_cov_result.sampleID_list.size())
        {
            std::cerr << "ERROR: Sample IDs in pheno file are more than covariate file " << '\n';
            exit(EXIT_FAILURE);
        }

        if (values[hdr_id_indx] != shared_cov_result.sampleID_list[row_indx]) 
        {
            std::cerr << "ERROR: Sample ID mismatch at line " << row_indx + 1
                    << ". Expected: " << shared_cov_result.sampleID_list[row_indx]
                    << ", Found: " << values[hdr_id_indx] << '\n';
            exit(EXIT_FAILURE);
        }
        
        bool invalid_indices = false;

        for(int i = 2; i < num_columns; i++)
        {
            if(values[i] == opt.missing_key || values[i].empty())
            {
                shared_phenotype_data[i - 2].push_back(opt.missing_key);
                invalid_indices = true;
            }
            else
            {
                shared_phenotype_data[i - 2].push_back(values[i]);
            }
        }

        if(!invalid_indices)
        {
            shared_pheno_valid_indices.insert(row_indx);
        }
        row_indx++;
    }
    
    if (row_indx < shared_cov_result.sampleID_list.size())
        {
            std::cerr << "ERROR: Sample IDs in covariate file are more than pheno file " << '\n';
            exit(EXIT_FAILURE);
        }
}


/**
 * @brief Remove lines with missing data coming from phenotype file
 * 
 * @param shared_cov_result.sampleID_list
 
 * @param shared_pheno_valid_indices 
 * @param shared_cov_result.covMap 
 */
void GEMRunner::clean_covMap_by_invalid_indices() 
{
    // Keep each sample ID and it is line number
    std::unordered_map<std::string, int> seen_count;

    // First pass: track which positions to delete per ID
    std::unordered_map<std::string, std::ext::V_int> delete_positions;

    for (size_t i = 0; i < shared_cov_result.sampleID_list.size(); ++i) 
    {
        // Find the occurance of line inside covData map
        const std::string& id = shared_cov_result.sampleID_list[i];
        int pos = seen_count[id]++;
        if (!shared_pheno_valid_indices.count(i)) 
        {
            delete_positions[id].push_back(pos);
        }
    }

    // Second pass: remove entries in reverse to preserve indexing
    for (auto& [id, positions] : delete_positions) 
    {
        auto it = shared_cov_result.covMap.find(id);
        if (it == shared_cov_result.covMap.end()) continue;

        auto& vecs = it->second;

        // Sort in descending order to erase from back to front
        std::sort(positions.rbegin(), positions.rend());
        // Remove the missing line from cov data 
        for (int pos : positions) 
        {
            if (pos >= 0 && pos < static_cast<int>(vecs.size())) 
            {
                vecs.erase(vecs.begin() + pos);
            }
        }

        if (vecs.empty()) 
        {
            shared_cov_result.covMap.erase(id);
        }
    }
}


/**
 * @brief Check if user provided the kinship file
 * 
 */

void GEMRunner::check_kinship_usage() 
{
    kin_flag = !opt.kin_path.empty();  // sets true if user provided kin_path
}


/**
 * @brief calculate dosages for bgen allesss per each sample
 * 
 * @param start_chunck 
 * @param chunks_to_read 
 * @return py::array_t<float> 
 */

/**
 * @brief Function to call null model from python
 * 
 */
void GEMRunner::run_fit_nullmodel() 
{
    NullModel model(opt);
    model.fit_nullmodel(kin_flag,
                        bgen_sample_id,
                        dup_id,
                        shared_pheno_valid_indices,
                        shared_colnames,
                        shared_phenotype_data,
                        c2_values // Pass C2 values by reference
    );
}
