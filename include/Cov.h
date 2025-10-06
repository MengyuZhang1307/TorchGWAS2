#pragma once

#include "ReadFiles.h"

class Cov
{
    public:
        DataFrame m_data_frame;
        std::string m_sam_id_hdr;
        std::ext::V_string m_v_hdrs;
        std::set<int> m_pheno_valid_indices;

        // std::string m_phenoMissingKey;
        /**
         * @brief a function to get the size of cov data
         * 
         * @return unsigned int 
         */
        unsigned int size();
        /**
         * @brief Set the path
         * 
         * @param path 
         */
        void set_path(std::string_view path);
        /**
         * @brief Get the path
         * 
         * @return std::string 
         */
        std::string get_path() const;
        /**
         * @brief A function to read covariate file with a given delimeter (default ",") and path
         * 
         * @param delim 
         */
        void read_file(std::string_view, char delim = ','); 
        void read_file(std::string_view path, char delim, 
            std::ext::V_string const& cov_col_names);
        std::pair<std::string, std::string> check_binary(std::ext::V_double const& ph_column);
    private:
        std::string m_path;
};