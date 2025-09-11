#pragma once
#include <armadillo>
#include <cholmod.h>
#include <cs.h> 
#include "Cov.h"
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "Kinship.h"

namespace std
{
    namespace ext
    {
        using Index_map = unordered_map<string, V_int>;
        using Index_map_kin = unordered_map<string, int>;
        using Triplet_d = Eigen::Triplet<double>;
        using Triplet_i = Eigen::Triplet<int>;
        using VecTuples4spmat = vector<Triplet_d>;
        using VecTriple_i = vector<Triplet_i>;


        // a template helper function and functions for printing tuples
        template<typename Tuple, size_t... Indices>
        void tuplePrint(const Tuple& t, std::index_sequence<Indices...>)
        {
            ((cout << std::get<Indices>(t) << " "), ...);
        }

        template<typename... Args>
        void tuplePrint(const std::tuple<Args...>& t)
        {
            tuplePrint(t, std::index_sequence_for<Args...>());
            cout << "\n";
        }
        
    }
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (std::pair<T1,T2> const& pair) const 
    {
        auto h1 = std::hash<T1>{}(pair.first);
        auto h2 = std::hash<T2>{}(pair.second);
        return h1 ^ h2; 
    }
};

using SpaMat = Eigen::SparseMatrix<double, Eigen::ColMajor>;
// using SpaMat = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
using Mat = Eigen::MatrixXd;
using DensMat = Eigen::MatrixXd;
using DensVec = Eigen::VectorXd;
using DensVecInt = Eigen::VectorXi;
using DenseMatInt = Eigen::MatrixXi;


class SparseInverse
{
    public:
        Cov cov;
        Kinship kin;
        char cov_delim;
        char kin_delim; 

        SparseInverse() =  default;
        SparseInverse(const std::string kin_add, const std::string cov_add, const char kin_delim,
                      const double kin_diag, const char cov_delim, const std::string &m_sam_id,
                      const std::ext::V_string &m_v_hdrs, std::ext::V_string &bgen_sample_id,
                      const std::string missing_key, std::set<int> pheno_valid_indices);
        void set_idx_mp(std::ext::V_string v_strs);
        //to create kinship from unique IDs
        void set_idx_mp_uniqkin(std::ext::V_string v_strs);
        std::ext::Index_map get_idx_mp();  
        std::ext::Index_map_kin get_idx_mp_uniqkin();       
        void set_spmat();
        void set_spmat(SpaMat sm);
        SpaMat& get_spmat();
        SpaMat& get_uniqkin();
        SpaMat inv_spamat();
        //define static as we might don't create an object
        static SpaMat inv_spamat(SpaMat const& sm);
        static DensMat inv(DensMat const &dm);

    private:
       std::ext::Index_map m_idx_mp;
       std::ext::Index_map_kin m_idx_mp_uniq;
       SpaMat m_spmat;
       //if there are duplicated IDs create a kin from unique IDS
       SpaMat m_uniqkinmat;//to generate N in N kinship
       [[nodiscard]] std::ext::VecTuples4spmat create_tuple4spmat();
       [[nodiscard]] void fill_vt4spmat(std::ext::VecTuples4spmat &vt4spmat);
       [[nodiscard]] std::ext::VecTuples4spmat create_tuple4_uniqkin();
       [[nodiscard]] void fill_vt4uniqkin(std::ext::VecTuples4spmat &vt4spmat);
       [[nodiscard]] bool is_missing(int id1, int id2);
};