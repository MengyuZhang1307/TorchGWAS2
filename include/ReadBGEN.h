#pragma once
#include "../thirdparty/zstd-1.5.5/lib/zstd.h"
#include "../thirdparty/libdeflate-1.18/libdeflate.h"
#include "zlib.h"
#include "ReadFiles.h"

class Bgen 
{
    public:

        // For file
        FILE* fin;

        // For BGEN offset
        uint32_t offset;

        // For BGEN header block
        uint32_t Mbgen;
        uint32_t Nbgen;
        uint32_t CompressedSNPBlocks;
        uint32_t Layout;

        // For BGEN header-flag block;
        uint32_t SampleIdentifiers;


        // For ID matching
        int new_samSize;
        std::vector<std::string>   sampleID;
        //AllsampleIDs before matching
        std::vector<std::string>   sampleID_all;
        std::vector<double>   new_covdata;
        std::vector<double>   new_phenodata;
        std::vector<long int> include_idx;
        std::vector <long int> variant_pos;
        std::vector<unsigned int> includeVariantIndex;
        // For check of co-linear relations between covX;
        int numIntSelCol_new;
        int numExpSelCol_new;
        int numSelCol_new;
        std::vector<int> excludeCol;
        // For multithreading BGEN file
        int phenoType;
        uint32_t threads;
        bool filterVariants;
        std::vector<uint32_t> Mbgen_begin;
        std::vector<uint32_t> Mbgen_end;
        std::vector<long long unsigned int> bgenVariantPos;
        std::vector<std::vector<uint32_t> > keepVariants;

        void processBgenHeaderBlock(std::string bgenfile);
        void processBgenSampleBlock(Bgen bgen, const char sample_file[300], bool use_sample, UMap_str_VV_string phenomap, std::string pheno_missing_key, int numSelCol, int sam_size);
        void getPositionOfBgenVariant(Bgen bgen, int threads, std::string includeVariantFile, bool do_filters);
};

// void gemBGEN(int thread_num, double sigma2, double* resid, double* XinvXTX, vector<double> miu, BinE binE, Bgen bgen, CommandLine cmd);
void Bgen13GetTwoVals(const unsigned char* prob_start, uint32_t bit_precision, uintptr_t offset, uintptr_t* first_val_ptr, uintptr_t* second_val_ptr);
VVV_float  calcDosage(std::string bgenFile, int stream_snps, int thread_num, Bgen bgen);



