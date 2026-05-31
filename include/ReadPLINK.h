#pragma once
#include "../thirdparty/plink-2.0/pgenlib_read.h"
#include "ReadFiles.h"
#include "BoundedQueue.h"

using uchar = unsigned char;

/**
 * @brief Class for reading PLINK format genotype files (BED/BIM/FAM and PGEN/PVAR/PSAM)
 *
 * This class provides functionality to read both PLINK 1.x (BED/BIM/FAM) and
 * PLINK 2.0 (PGEN/PVAR/PSAM) format files. It uses the pgenlib library for
 * efficient reading and handles sample/variant filtering and ID matching.
 */
class Plink
{
    public:
        // File format type: "BED" for PLINK 1.x, "PGEN" for PLINK 2.0
        std::string format_type;

        // File paths
        std::string pgen_path;      // .pgen or .bed file path
        std::string pvar_path;      // .pvar or .bim file path
        std::string psam_path;      // .psam or .fam file path

        // pgenlib structures for reading
        plink2::PgenFileInfo pgfi;
        plink2::PgenReader pgr;

        // Sample information (from .fam/.psam)
        uint32_t raw_sample_ct;     // Total samples in file
        int new_samSize;             // Samples after filtering/matching
        std::ext::V_string sampleID;      // Matched sample IDs (in output order)
        std::ext::V_string sampleID_all;  // All sample IDs before matching
        std::ext::V_int plink_to_out;     // Mapping: PLINK sample idx -> output idx (-1 if excluded)
        std::ext::V_double new_covdata;   // Updated covariate data (for collinearity check)

        // Variant information (from .bim/.pvar)
        uint32_t raw_variant_ct;    // Total variants in file
        std::vector<std::string> variant_ids;      // Variant IDs (rsID or snpID)
        std::vector<std::string> chromosome;       // Chromosome names
        std::vector<uint32_t> base_pair_pos;       // Base pair positions
        std::vector<std::string> allele_ref;       // Reference alleles
        std::vector<std::string> allele_alt;       // Alternate alleles

        // For variant filtering
        bool filterVariants;
        std::vector<long int> include_idx;              // Indices of variants to include
        std::vector<uint32_t> includeVariantIndex;      // Variant indices after filtering
        std::vector<std::vector<uint32_t>> keepVariants; // Variants to keep per thread block

        // For covariate collinearity checking
        int numIntSelCol_new;
        int numExpSelCol_new;
        int numSelCol_new;
        std::ext::V_int excludeCol;

        // For multi-threading
        uint32_t threads;
        std::vector<uint32_t> variant_begin;  // Start variant index for each thread block
        std::vector<uint32_t> variant_end;    // End variant index for each thread block

        // Phenotype type (for future use)
        int phenoType;

        // pgenlib initialization state
        bool pgfi_inited = false;
        unsigned char* pgfi_alloc = nullptr;
        uint32_t header_ctrl = 0;
        uint32_t max_vrec_width = 0;
        uintptr_t pgr_alloc_cacheline_ct = 0;

        // Variant filtering (set from GEMOptions before calling calc_dosage_plink)
        std::string includeVariantFile;

        /**
         * @brief Initialize PLINK file reader and read header information
         *
         * Detects file format (BED vs PGEN) and reads the corresponding header.
         * For BED files, also reads BIM and FAM files.
         * For PGEN files, also reads PVAR and PSAM files.
         *
         * @param pgen_or_bed_file Path to .pgen or .bed file
         */
        void process_plink_header_block(std::string const& pgen_or_bed_file);

        /**
         * @brief Read sample information and perform ID matching
         *
         * Reads sample IDs from .fam/.psam file and matches them with covariate file.
         * Handles missing covariate values and creates sample mapping.
         *
         * @param fam_or_psam_file Path to .fam or .psam file (can be empty if auto-detected)
         * @param use_fam_psam Whether to use provided fam/psam file path
         * @param covmap Covariate data map (sample_id -> covariate values)
         * @param pheno_missing_key String representing missing values
         * @param numSelCol Number of selected covariate columns
         * @param sam_size Number of samples in covariate file
         * @param id_path Path to ID file for matching (optional)
         * @param match_ids Whether to perform ID matching
         */
        void process_plink_sample_block(const char fam_or_psam_file[300],
                                       bool use_fam_psam,
                                       std::ext::UMap_str_VV_string covmap,
                                       std::string pheno_missing_key,
                                       int numSelCol,
                                       int sam_size,
                                       std::string id_path = "",
                                       bool match_ids = false);

        /**
         * @brief Read variant information and prepare for multi-threaded processing
         *
         * Reads variant IDs, positions, and alleles from .bim/.pvar file.
         * Optionally filters variants based on an include list.
         * Divides variants into blocks for multi-threaded reading.
         *
         * @param threads Number of threads for parallel processing
         * @param includeVariantFile Path to file containing variants to include (optional)
         * @param do_filters Whether to apply variant filtering
         */
        void get_variant_positions(int threads,
                                  std::string includeVariantFile,
                                  bool do_filters);

        /**
         * @brief Destructor - cleanup pgenlib resources
         */
        ~Plink();

    private:
        /**
         * @brief Read .fam file (PLINK 1.x format)
         *
         * @param fam_file Path to .fam file
         * @return Vector of sample IDs (FID_IID format)
         */
        std::ext::V_string read_fam_file(std::string const& fam_file);

        /**
         * @brief Read .psam file (PLINK 2.0 format)
         *
         * @param psam_file Path to .psam file
         * @return Vector of sample IDs
         */
        std::ext::V_string read_psam_file(std::string const& psam_file);

        /**
         * @brief Read .bim file (PLINK 1.x format)
         *
         * @param bim_file Path to .bim file
         */
        void read_bim_file(std::string const& bim_file);

        /**
         * @brief Read .pvar file (PLINK 2.0 format)
         *
         * @param pvar_file Path to .pvar file
         */
        void read_pvar_file(std::string const& pvar_file);

        /**
         * @brief Detect file paths for BIM and FAM files based on BED file path
         *
         * @param bed_file Path to .bed file
         */
        void detect_bed_companion_files(std::string const& bed_file);

        /**
         * @brief Detect file paths for PVAR and PSAM files based on PGEN file path
         *
         * @param pgen_file Path to .pgen file
         */
        void detect_pgen_companion_files(std::string const& pgen_file);
};

/**
 * @brief Stream genotype dosages from PLINK file in chunks
 *
 * Multi-threaded function to read genotype data from PLINK files and convert to dosages.
 * Dosages are calculated as: dosage = 0*P(AA) + 1*P(Aa) + 2*P(aa)
 * where AA is homozygous reference, Aa is heterozygous, aa is homozygous alternate.
 *
 * @param plinkFile Path to .pgen or .bed file
 * @param plink Reference to Plink object with initialized file reader
 * @param queue Bounded queue for streaming chunks to consumer
 * @param threads Number of threads for parallel reading
 * @param snps_per_chunk Number of SNPs to read per chunk
 */
void calc_dosage_plink(std::string const& plinkFile,
                      Plink& plink,
                      BoundedChunkQueue& queue,
                      int threads,
                      int snps_per_chunk);
