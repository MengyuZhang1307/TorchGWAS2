#include <gtest/gtest.h>
#include "../ReadBGEN.h" 
#include "../ReadFiles.h"

// Optional: mock data helpers
// std::string write_temp_file(const std::string& content, const std::string& filename) {
//     std::ofstream out(filename);
//     out << content;
//     out.close();
//     return filename;
// }

// Test getOptions
TEST(GEMOptionsTest, ParsesCorrectly) {
    char* argv[] = {
        "program",
        "pheno.txt", "cov.txt", "bgen.bgen", "sample.sample",
        "eid", "sex, coll", "-9", "out.txt", "4",
        "\t", "\t",
        "age,sex", "bmi", "bmi*sex"
    };
    int argc = sizeof(argv) / sizeof(argv[0]);

    GEMOptions opt = getOptions(argc, argv);
    EXPECT_EQ(opt.phenoFile, "pheno.txt");
    EXPECT_EQ(opt.threads, 4);
    EXPECT_EQ(opt.delim_cov, '\t');
    EXPECT_EQ(opt.covariates.size(), 2);
    EXPECT_EQ(opt.interactions[0], "bmi*sex");
}

TEST(CovariateTest, ParsesCovariateCorrectly) {
    std::string cov_file = "example.cov-2id-sorted";  // your real file path

    std::vector<std::string> covariates = {"cov1", "cov2"};
    std::vector<std::string> exposures = {"cov3"};
    std::vector<std::string> interactions = {};
    std::string sampleID = "sampleid";
    std::string slope = "";
    char delim = ',';  // comma-separated
    std::string missing = "NA";

    CovariateReadResult res = read_covariate_data(
        cov_file, covariates, exposures, interactions,
        sampleID, slope, delim, missing
    );

    EXPECT_EQ(res.sampleID_list.size(), 2500);
    EXPECT_EQ(res.sampleID_list[0], "sample_001");
    EXPECT_EQ(res.covMap["sample_001"][0][0], "-1.049");       // cov1
    // EXPECT_EQ(res.covMap["sample_356"][0][2], "0.402");   // cov3
    // EXPECT_EQ(res.covMap["sample_436"][0][1], "NA");      // cov2
}


TEST(PhenotypeProcessingTest, HandlesCleanPhenotypeFile) {
    // Assume your file is saved at this path
    std::string pheno_path = "example.pheno2-2id-sorted";
    std::string sample_id_hdr = "ID";
    char delim = ',';
    std::string missing = "NA";

    // Simulate a covMap for all sample IDs
    UMap_str_VV_string covMap = {
    {"sample_001", {{"1", "2", "3"}, {"1", "2", "3"}}},
    {"sample_002", {{"5", "6", "7"}, {"5", "6", "7"}}}
};

    V_string sampleID_list = {"sample_001","sample_001","sample_001","sample_001","sample_001", "sample_002", "sample_002", "sample_002", "sample_002", "sample_002"};

    std::set<int> valid_indices;
    V_string column_names;
    VV_string phenotype_data;

    process_phenotype_file(pheno_path, sample_id_hdr, sampleID_list, valid_indices, column_names, phenotype_data, delim, missing);

    // Check valid samples (no missing)
    EXPECT_EQ(valid_indices.size(), 9);
    // EXPECT_TRUE(valid_indices.count(0));
    EXPECT_TRUE(valid_indices.count(1));
    EXPECT_TRUE(valid_indices.count(2));

    for(auto ind : valid_indices)
    {
        std::cout << ind << '\n';
    }

    // Check column names
    EXPECT_EQ(column_names[0], "ID");
    EXPECT_EQ(column_names[1], "sampleid");
    EXPECT_EQ(column_names[2], "pheno1");

    // // Check actual phenotype data
    EXPECT_EQ(phenotype_data[0][0], "NA"); // pheno2 for sample_500
    EXPECT_EQ(phenotype_data[0][1], "1.457"); // pheno2 for sample_498
    EXPECT_EQ(covMap.size(), 2);
    std::cout << covMap.size();
    // // Now test covMap cleanup
    // update_covMap(covMap, valid_indices);
    clean_covMap_by_invalid_indices(sampleID_list, valid_indices, covMap);

    std::cout << covMap.size();
    EXPECT_EQ(covMap.size(), 2);
    for(auto& [key, val] : covMap)
    {
        std::cout << '\n' << "key in testfile " << '\n' << key << '\n';
    }
    std::cout << '\n' << covMap["sample_001"].size();
    std::cout << '\n' << covMap["sample_002"].size() << '\n';
    EXPECT_TRUE(covMap.count("sample_001"));
    EXPECT_TRUE(covMap.count("sample_002"));
}
