#include "Cov.h"
#include <iostream>
#include <iterator>
#include <unordered_set>


void Cov::read_file(std::string_view path, char delim, 
    std::ext::V_string const& cov_col_names)
{
    m_data_frame.read_file(path, delim, cov_col_names);
    //std::ext::V_string tmp_hdrs = {m_sam_id};
    // m_data_frame = m_data_frame.copy_by_hdrs(m_v_hdrs);
}

std::pair<std::string, std::string> Cov::check_binary(std::ext::V_double const& ph_column)
{
    std::unordered_set<double> s_tmp;
    for(auto ph : ph_column)
    {
        s_tmp.insert(ph);
        if(s_tmp.size() > 2)
        {
            return {"gaussian", "identity"};
        }
    }
    return {"binomial", "logit"};
}

unsigned int Cov::size()
{
    return m_data_frame.m_data[m_sam_id].size();
}


void Cov::set_path(std::string_view path)
{
    m_path = path.data();
}


std::string Cov::get_path() const
{
    return m_path;
}
