#include "ReadFiles.h"

char resolve_delim(const std::string& s)
{
    if (s == "\\t") return '\t';
    if (s == "\\n") return '\n';
    if (s == "\\s") return ' ';  // optional
    return s[0];  // just return first char
}

std::vector<std::string> split(const std::string& s, char delim)  
{
        std::vector<std::string> out;
        std::istringstream ss(s);
        std::string token;
        while (std::getline(ss, token, delim)) {
            out.push_back(token);
        }
        return out;
}

GEMOptions get_options(int argc, const char* argv[]) 
{
    if (argc < 14) 
    {
        std::cerr << "Usage: <program> <pheno.txt> <cov.txt> <kin.txt> <bgen.bgen> <sample.sample> "
                     "<sample_id_header> <random_slope_header> <missing_key> <output.txt> <threads> "
                     "<delim_pheno> <delim_cov> <delim_kin> <cov1,cov2,...> <exp1,exp2,...> <int1,int2,...>\n";
        exit(EXIT_FAILURE);
    }

    GEMOptions opt;
    opt.pheno_file  = argv[1];
    opt.cov_file  = argv[2];
    opt.delim_pheno  = resolve_delim(argv[3]);
    opt.delim_cov = resolve_delim(argv[4]);
    opt.bgen_file = argv[5];
    opt.sample_file = argv[6];
    opt.do_filters = std::string(argv[7]) == "true";
    opt.use_sample_file = std::string(argv[8]) == "true";
    opt.includeVariantFile = argv[9];
    opt.stream_snps = std::stoi(argv[10]);
    opt.sampleid_header_name = argv[11];
    opt.random_slope_header_name= argv[12];
    opt.covariates = split(argv[13], ',');
    opt.exposures = split(argv[14], ',');
    opt.interactions = split(argv[15], ',');
    opt.missing_key = argv[16];
    opt.threads = std::stoi(argv[17]);
    opt.out_file = argv[18];
    return opt;
}


CovariateReadResult read_covariate_data(
    const std::string& cov_file,
    std::vector<std::string> cov_sel_headers_name,
    const std::vector<std::string>& exp_cov_sel_headers_name,
    const std::vector<std::string>& int_cov_sel_headers_name,
    const std::string& sampleid_header_name,
    const std::string& random_slope_header_name,
    char delim_cov,
    const std::string missing_key) 
{
    CovariateReadResult result;
    bool& cov_is_duplicated = result.cov_is_duplicated;
    cov_is_duplicated = false;

    int numExpSelCol = exp_cov_sel_headers_name.size();
    int numIntSelCol = int_cov_sel_headers_name.size();

    for (int i = numIntSelCol - 1; i >= 0; --i)
        cov_sel_headers_name.insert(cov_sel_headers_name.begin(), int_cov_sel_headers_name[i]);
    for (int i = numExpSelCol - 1; i >= 0; --i)
        cov_sel_headers_name.insert(cov_sel_headers_name.begin(), exp_cov_sel_headers_name[i]);

    result.numSelCol = cov_sel_headers_name.size() - numExpSelCol - numIntSelCol;
    std::vector<int> colSelVec(cov_sel_headers_name.size());

    std::ifstream fincov(cov_file);
    if (fincov.fail()) 
    {
        throw std::runtime_error("ERROR: Cannot open covariate file");
    }

    std::unordered_map<std::string, int> colNames;
    std::string line, headerName;
    std::getline(fincov, line);
    std::istringstream issHead(line);
    int header_i = 0;
    while (std::getline(issHead, headerName, delim_cov)) 
    {
        headerName.erase(std::remove(headerName.begin(), headerName.end(), '\r'), headerName.end());
        headerName.erase(std::remove(headerName.begin(), headerName.end(), '"'), headerName.end());
        if (colNames.count(headerName)) 
        {
            throw std::runtime_error("ERROR: Duplicate header in covariate file: " + headerName);
        }
        colNames[headerName] = header_i++;
    }

    if (!colNames.count(sampleid_header_name))
        throw std::runtime_error("ERROR: Sample ID column not found");
    int SamIDCol = colNames[sampleid_header_name];

    if (!random_slope_header_name.empty()) 
    {
        if (!colNames.count(random_slope_header_name))
            throw std::runtime_error("ERROR: Random slope column not found");
    }

    for (const auto& h : exp_cov_sel_headers_name)
        if (!colNames.count(h)) throw std::runtime_error("ERROR: Exposure column not found: " + h);
    for (const auto& h : int_cov_sel_headers_name)
        if (!colNames.count(h)) throw std::runtime_error("ERROR: Interaction column not found: " + h);
    for (size_t i = 0; i < cov_sel_headers_name.size(); ++i) 
    {
        if (!colNames.count(cov_sel_headers_name[i]))
        {
            throw std::runtime_error("ERROR: Covariate column not found: " + cov_sel_headers_name[i]);
        }
        colSelVec[i] = colNames[cov_sel_headers_name[i]];
    }

    int nrows = 0;
    while (std::getline(fincov, line)) ++nrows;
    result.samSize = nrows;
    fincov.clear();
    fincov.seekg(0);
    std::getline(fincov, line); // skip header

    for (int r = 0; r < result.samSize; ++r) {
        std::getline(fincov, line);
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        std::istringstream iss(line);
        std::string value;
        std::vector<std::string> values;
        while (std::getline(iss, value, delim_cov)) values.push_back(value);
        if (values.size() == colNames.size() - 1) values.push_back("");
        if (values.size() != colNames.size())
            throw std::runtime_error("ERROR: Column mismatch at row " + std::to_string(r));

        if (!cov_is_duplicated && result.covMap.count(values[SamIDCol]))
            cov_is_duplicated = true;

        std::vector<std::string> entry;
        for (int c : colSelVec) 
        {
            std::string val = values[c];
            val.erase(std::remove(val.begin(), val.end(), '"'), val.end());
            entry.push_back(val);
        }

        result.sampleID_list.push_back(values[SamIDCol]);
        result.covMap[values[SamIDCol]].push_back(entry);
    }

    fincov.close();


    return result;
}


void process_phenotype_file(
    std::string const& filename, 
    std::string sample_id_hdr,
    V_string sampleID_list,
    std::set<int> &valid_indices,
    V_string &column_names,
    VV_string &phenotype_data,
    char delim,
    std::string missing_key
    ) 
{
    std::unordered_set<std::string> seen;
    std::ifstream file(filename);

    if (!file.is_open()) 
    {
        std::cerr << "Error opening file: " << filename << std::endl;
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
    
    while (std::getline(ss, col_name, delim)) 
    {
        if (seen.insert(col_name).second)
        {
            column_names.push_back(col_name);
            if(col_name == sample_id_hdr)
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
    
    int num_columns = column_names.size();
    std::cout << "Total columns: " << num_columns << "\n"; 
    std::cout << "****************************************************************************\n";
    phenotype_data.resize(num_columns - 1);

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
        V_string values;
        while(getline(ss, value, delim))
        {
            values.push_back(value);
        }
                    
        if (!line.empty() && line.back() == delim) 
        {
            values.push_back(missing_key);
        }

        if (values.size() != num_columns) 
        {
            std::cerr << "ERROR: expect: " << num_columns << " columns at row: " << row_indx + 1 << " , while there is: " << values.size() << " columns."<< '\n';
            exit(EXIT_FAILURE);
        }

        if (row_indx >= sampleID_list.size())
        {
            std::cerr << "ERROR: Sample IDs in pheno file are more than covariate file " << '\n';
            exit(EXIT_FAILURE);
        }

        if (values[hdr_id_indx] != sampleID_list[row_indx]) 
        {
            std::cerr << "ERROR: Sample ID mismatch at line " << row_indx + 1
                    << ". Expected: " << sampleID_list[row_indx]
                    << ", Found: " << values[hdr_id_indx] << '\n';
            exit(EXIT_FAILURE);
        }
        
        bool invalid_indices = false;

        for(int i = 2; i < num_columns; i++)
        {
            if(values[i] == missing_key || values[i].empty())
            {
                phenotype_data[i - 2].push_back(missing_key);
                invalid_indices = true;
            }
            else
            {
                phenotype_data[i - 2].push_back(values[i]);
            }
        }

        if(!invalid_indices)
        {
            valid_indices.insert(row_indx);
        }
        row_indx++;
    }
    
    if (row_indx < sampleID_list.size())
        {
            std::cerr << "ERROR: Sample IDs in covariate file are more than pheno file " << '\n';
            exit(EXIT_FAILURE);
        }

}

// void update_covMap(UMap_str_VV_string &covMap, std::set<int> pheno_valid_indices)
// {

//     std::vector<std::string> keys_to_remove;
//     size_t i = 0;


//     for (const auto& [key, val] : covMap)
//     {
//         if (!pheno_valid_indices.count(i)) 
//         {
//             keys_to_remove.push_back(key);
//             std::cout << "\n inside readfile \n" << i << "  "  << key << '\n';
//         }
//         ++i;
//     }

//     for (const auto& key : keys_to_remove)
//     {
//         covMap.erase(key);
//     }

// }


void clean_covMap_by_invalid_indices(
    std::vector<std::string> const& sampleID_list,
    std::set<int> const& pheno_valid_indices,
    UMap_str_VV_string& covMap
) 
{
    std::unordered_map<std::string, int> seen_count;

    // First pass: track which positions to delete per ID
    std::unordered_map<std::string, std::vector<int>> delete_positions;

    for (size_t i = 0; i < sampleID_list.size(); ++i) 
    {
        const std::string& id = sampleID_list[i];
        int pos = seen_count[id]++;
        if (!pheno_valid_indices.count(i)) 
        {
            delete_positions[id].push_back(pos);
        }
    }

    // Second pass: remove entries in reverse to preserve indexing
    for (auto& [id, positions] : delete_positions) 
    {
        auto it = covMap.find(id);
        if (it == covMap.end()) continue;

        auto& vecs = it->second;

        // Sort in descending order to erase from back to front
        std::sort(positions.rbegin(), positions.rend());

        for (int pos : positions) 
        {
            if (pos >= 0 && pos < static_cast<int>(vecs.size())) 
            {
                vecs.erase(vecs.begin() + pos);
            }
        }

        if (vecs.empty()) 
        {
            covMap.erase(id);
        }
    }
}

// Optional: pass delimiter (default is comma)
void write_bgen_result_to_file(const V_bgen& results, const std::string& filename, char delimiter) 
{
    std::ofstream out(filename);
    if (!out.is_open()) 
    {
        std::cerr << "Error: could not open file " << filename << " for writing.\n";
        return;
    }

    for (const auto& thread_vec : results) 
    {
        for (const auto& snp_vec : thread_vec) 
        {
            for (const auto& sample_vec : snp_vec) 
            {
                for (size_t i = 0; i < sample_vec.size(); ++i) 
                {
                    out << sample_vec[i];
                    if (i != sample_vec.size() - 1)
                        out << delimiter;
                }
                out << '\n';  // new line for each innermost vector
            }
        }
    }

    out.close();
    std::cout << "\u2705 BGEN results written to: " << filename << '\n';
}
