// #include "../include/ReadFiles.h"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include "ReadFiles.h"
#include <unordered_set>

bool DataFrame::isDataFrameCreated() const 
{
    return (m_nrows != 0 && m_ncols != 0);
}


int DataFrame::n_rows() const
{ 
    if (!isDataFrameCreated()) 
    {
        std::cerr << "Error: DataFrame has not been created.\n";
        exit(EXIT_FAILURE); 
    }
    return m_nrows; 
}


int DataFrame::n_cols() const
{
    if (!isDataFrameCreated()) 
    {
        std::cerr << "Error: DataFrame has not been created.\n";
        exit(EXIT_FAILURE); 
    }
    return m_ncols; 
}


std::ext::V_string  DataFrame::read_lines(std::string_view path)
{
    std::ifstream ifs(path.data());

    if(ifs.fail())
    {
        fmt::print("Error in reading file {}.\nPlease check your inputs.\n", path);
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::ext::V_string v_strs;

    while(std::getline(ifs, line))
    {
        line.erase(std::remove(line.begin(), line.end(), '\r'),line.end());
        v_strs.emplace_back(line);
    }

    return v_strs;
}

void DataFrame::fill_data(std::ext::V_string const& lines, char delim)
{
    std::ext::VV_string vv_strs;
    for(auto const& line : lines)
    {
        std::istringstream iss(line);
        std::string cell;
        std::ext::V_string v_str_tmp;
        while(std::getline(iss, cell, delim))
        {
            cell.erase(std::remove(cell.begin(), cell.end(), '\"'), cell.end());
            // v_str_tmp.emplace_back(cell);
            if (cell.empty())
            {   
                v_str_tmp.emplace_back(m_missing_key);
            }
                else
            {
                v_str_tmp.emplace_back(cell);
            }
        }
        vv_strs.emplace_back(v_str_tmp);
    }
    m_nrows = vv_strs.size();
    m_ncols = vv_strs[0].size();

    m_headers.resize(m_ncols);

    for(int i{0}; i < m_ncols; ++i)
    {
        m_headers[i] = vv_strs[0][i];
        std::ext::V_string values;

        for(int j{1}; j < m_nrows; ++j)
        {
            if(vv_strs[j].size() != m_ncols)
            {
                values.emplace_back(m_missing_key);
                continue;
            }
            
            values.emplace_back(vv_strs[j][i]);
        }
        m_data[m_headers[i]] = values;
    }
    m_nrows = m_data[m_headers[0]].size();
}

void DataFrame::fill_data(std::ext::V_string const& lines, char delim,
                          std::ext::V_string const& cov_col_names)
{
    std::ext::VV_string vv_strs;
    for(auto const& line : lines)
    {
        std::istringstream iss(line);
        std::string cell;
        std::ext::V_string v_str_tmp;
        while(std::getline(iss, cell, delim))
        {
            cell.erase(std::remove(cell.begin(), cell.end(), '\"'), cell.end());
            if (cell.empty())
                v_str_tmp.emplace_back(m_missing_key);
            else
                v_str_tmp.emplace_back(cell);
        }
        vv_strs.emplace_back(v_str_tmp);
    }

    // all headers
    std::ext::V_string all_headers = vv_strs[0];

    // find indices of headers to keep
    std::ext::V_int cov_col_indices;
    for (int i = 0; i < all_headers.size(); ++i) 
    {
        if (std::find(cov_col_names.begin(), cov_col_names.end(), all_headers[i]) != cov_col_names.end()) 
        {
            cov_col_indices.push_back(i);
        }
    }

    m_ncols = cov_col_indices.size();
    m_headers.resize(m_ncols);

    // fill data only for kept headers
    for (int idx = 0; idx < cov_col_indices.size(); ++idx) 
    {
        int col = cov_col_indices[idx];
        m_headers[idx] = all_headers[col];
        std::ext::V_string values;

        for (int j = 1; j < vv_strs.size(); ++j) 
        {
            if (vv_strs[j].size() <= col) 
            {
                values.emplace_back(m_missing_key);
            } 
            else 
            {
                values.emplace_back(vv_strs[j][col]);
            }
        }
        m_data[m_headers[idx]] = values;
    }

    m_nrows = m_data[m_headers[0]].size();
}



void DataFrame::head(int n)
{
        for(auto const& curr_hdr : m_headers)
        {
            std::cout << curr_hdr << " ";
        }

        std::cout << "\n";

        for(int i{0}; i <= n; ++i)
        {
            for(auto const& curr_hdr : m_headers)
            {
                std::cout << m_data[curr_hdr][i] << " ";
            }

            std::cout << "\n";
        }

}

void DataFrame::read_file(std::string_view path, char delim)
{
    auto v_strs = read_lines(path);
    fill_data(v_strs, delim);
}

void DataFrame::read_file(std::string_view path, char delim, 
                std::ext::V_string const& cov_col_names)
{
    auto v_strs = read_lines(path);
    fill_data(v_strs, delim, cov_col_names);
}

std::ext::V_string DataFrame::get_header(std::string const& hdr) const
{
    auto it = m_data.find(hdr);

    if (it != m_data.end()) 
    {
        return it->second;
    } 
    else 
    {
        fmt::print("Error: Header {} not found. It must be one of the variables in the submitted file\n", hdr);
        exit(EXIT_FAILURE); // Return an empty vector
    }
}

void DataFrame:: remove_missing(std::ext::V_string const& v_hdrs, std::ext::V_int const& valid_indices)
{
    for (const auto& hdr : v_hdrs)
    {
        auto it = m_data.find(hdr);
        if (it != m_data.end())
        {
            std::ext::V_string filtered_data;
            filtered_data.reserve(valid_indices.size()); // Reserve space for efficiency

            for (int idx : valid_indices)
            {
                if (idx < it->second.size())
                {
                    filtered_data.push_back(it->second[idx]);
                }
            }

            // Swap the filtered data back into the original data structure
            it->second.swap(filtered_data);
        }
    }
    m_nrows = m_data[m_headers[0]].size();
}
//Remove rows with missing data and match phenoIDs with genofile IDs
void DataFrame::match_genoids(std::string hdr_id, std::ext::V_string const& v_hdrs)
{
    // Convert m_geno_IDs(bgenIDs) to an unordered set for fast lookup
    std::unordered_set<std::string> geno_id_set(m_geno_ids.begin(), m_geno_ids.end());
    const std::ext::V_string& sam_ids = m_data[hdr_id];
    std::vector<int> keep_indices;

    for (int i = 0; i < sam_ids.size(); ++i)
    {
        if (geno_id_set.count(sam_ids[i]) > 0)
        {
            keep_indices.push_back(i);
        } 
    }

    for (int idx : keep_indices)
    {
        bool is_valid_row = true;

        for (const auto& hdr : v_hdrs)
        {
            auto it = m_data.find(hdr);
            if (it != m_data.end())
            {
                if (idx < it->second.size() && it->second[idx] == m_missing_key) //check missing value as well
                {
                    // std::cerr << "Warning: missing value at row: " << idx + 1 << " for header: " << hdr << "\n";
                    is_valid_row = false; // If there's any missing value, mark the row as invalid
                    break; // No need to check further headers for this row
                }
            }
        }

        if (is_valid_row)
        {
            m_valid_indices.push_back(idx);
        }
    }
    remove_missing(v_hdrs, m_valid_indices);
}


DataFrame DataFrame::copy_by_hdrs(std::ext::V_string const& v_hdrs)
{
    DataFrame new_DataFrame;
    new_DataFrame.m_ncols = v_hdrs.size();
    new_DataFrame.m_nrows = m_data[v_hdrs[0]].size();
    new_DataFrame.m_headers = v_hdrs;
    new_DataFrame.m_missing_key = m_missing_key;
    new_DataFrame.m_geno_ids = m_geno_ids;

    if (v_hdrs.empty()) 
    {
        throw std::invalid_argument("Header vector is empty");
    }
    
    for(unsigned int i{0}; i < v_hdrs.size(); ++i)
    {
        new_DataFrame.m_data[v_hdrs[i]] = m_data[v_hdrs[i]];
    }
    return new_DataFrame;
}


size_t DataFrame::size_wo_duplicates(std::string const& hdr)
{
    std::unordered_set<std::string> tmp_st;

    for(auto elm : m_data[hdr])  
    {
        tmp_st.insert(elm);
        
    } 

    return tmp_st.size();
}


bool DataFrame::any_duplicated(std::string const& hdr)
{
    std::unordered_set<std::string> tmp_st;
    for(auto elm : m_data[hdr])  
    {
        auto [it, success] = tmp_st.insert(elm);
        if(!success)
        {
            return true;
        }
    } 
    return false;
}


std::set<std::string> DataFrame::list_duplicates(std::string const& hdr)
{
    std::unordered_set<std::string> unique_id;
    std::ext::V_string list_dup;

    for(auto elm : m_data[hdr])
    {
        auto [it, success] = unique_id.insert(elm);
        if(!success)
        {
            list_dup.push_back(elm);
        }
    }
    std::set<std::string> unique_dup(list_dup.begin(), list_dup.end());
    return unique_dup;
}


std::ext::V_string DataFrame:: uniq_ids(std::string const& hdr) 
{
    std::unordered_set<std::string> seen; 
    std::ext::V_string list_uniq;

    for(auto elm : m_data[hdr])
    {
        auto [it, success] = seen.insert(elm);
        if(success)
        {
            list_uniq.push_back(elm);
        }
    }
    
    return list_uniq;
}
