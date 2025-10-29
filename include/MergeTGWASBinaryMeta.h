// #pragma once
// #include "GEMConfig.h"
// #include <nlohmann/json.hpp>

// #pragma pack(push, 1) //prevent add padding by compiler
// struct VariantMeta 
// {
//     char snpid[32];
//     char rsid[32];
//     int32_t chr;[2]
//     int32_t pos;
//     char non_effect_allele[8];
//     char effect_allele[8];
//     int32_t n_samples;
//     float af;
//     float gv;
// };
// #pragma pack(pop)

// /**
//  * @brief class to do post processing. Add meta data from genotype file.
//  * 
//  */
// class MergeMetaData
// {
//     private:
//         GEMOptions opt;
//         Bgem bgen;
//     public:
//         explicit MergeMetaData(GEMOptions const& user_opt);
//         void read_bgen_data();
//         void combine_metadata();
// };

// // helper functions signatures
// size_t count_snps(std::string const& filename);
// nlohmann::json meta_json read_json(std::string const& filename,
//                                     size_t snp_txt_count);