#include "ReadBGEN.h"

int main(int argc, char *argv[])
{
    // BGEN bgen;
    std::cout << "hello\n";
    std::tuple<std::string, int> options;
    GEMOptions opt = getOptions(argc, argv);
    CovariateReadResult result = read_covariate_data(
    opt.covFile,
    opt.covariates,
    opt.exposures,
    opt.interactions,
    opt.sampleIDHeaderName,
    opt.randomSlopeHeaderName,
    opt.delim_cov,
    opt.missingKey
    );
    // Bgen bgen;
    // bgen.processBgenHeaderBlock(bgenFile);
    // bgen.processBgenSampleBlock(bgen, samplefile, useSampleFile, phenomap, phenoMissingKey, numSelCol, samSize);
    // bgen.getPositionOfBgenVariant(bgen);
    // int stream_snps =1;
    // std::vector<std::vector<std::vector<float>>> dosageWhole;

    // if(threadsNum > 1)
    // {
    //     std::vector<std::thread> threads;
    //     for(int t  = 0; t < threadsNum; ++t)
    //     {
    //         threads.emplace_back([&](){
    //         std::vector<std::vector<float>> dosages =  calcDosage(bgenFile, stream_snps,  threadsNum,  bgen);
    //         dosageWhole[t] = std::move(dosages);
    //         });            
    //     }

    // }
    // else
    // {
    //     std::vector<std::vector<float>> dosages =  calcDosage(bgenFile, stream_snps,  threadsNum,  bgen);
    //     dosageWhole[0] = std::move(dosages);

    // }
    return 0;
}


