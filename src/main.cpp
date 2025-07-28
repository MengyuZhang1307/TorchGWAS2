#include "ReadBGEN.h"
#include "RunPipeline.h"
// To run use sample
// ./build/GEM example/example.pheno2-2id example/example.cov-2id , , example/example.bgen example/example.sample false true "" 1 sampleid  "" cov1,cov3 "" time NA 2 out1.txtC
int main(int argc, const char *argv[])
{
    // BGEN bgen;
    std::cout << "Reading BGEN file in multithread..................\n";

    GEMOptions opt = get_options(argc, argv);

    auto results = run_bgen(opt);
    write_bgen_result_to_file(results, opt.out_file);
    return 0;
}


