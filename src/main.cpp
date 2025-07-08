#include "ReadBGEN.h"
#include "RunPipeline.h"

int main(int argc, const char *argv[])
{
    // BGEN bgen;
    std::cout << "Reading BGEN file in multithread..................\n";

    GEMOptions opt = getOptions(argc, argv);

    auto results = run_bgen(opt);
    write_bgen_result_to_file(results, opt.out_file);
    return 0;
}


