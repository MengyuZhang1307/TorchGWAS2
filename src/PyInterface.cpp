#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "RunPipeline.h"
//Module to be called in python
PYBIND11_MODULE(Mygen, m) 
{
    m.doc() = R"doc(
        GEM interface for BGEN dosage calculation.

        Command sample:
        ./build/GEM pheno_file 
                    cov_file pheno_del cov_del bgen_file sample_bgen_file
                    do_filter use_sample_file includeVariantFile stream_snps sampleid
                    randomslope covariates interactions missing_key threads out_file
        Example:           
        ./build/GEM example/example.pheno2-2id example/example.cov-2id , , example/example.bgen 
        example/example.sample false true "" 1 sampleid  "" cov1,cov3 "" time NA 2 out1.txt
        Note:
        Please use "" if you do not want to pass a value for a specific argument.
        )doc";
    pybind11::class_<GEMOptions>(m, "GEMOptions")
        .def(pybind11::init<>())
        .def_readwrite("pheno_file", &GEMOptions::pheno_file)
        .def_readwrite("cov_file", &GEMOptions::cov_file)
        .def_readwrite("delim_pheno", &GEMOptions::delim_pheno)
        .def_readwrite("delim_cov", &GEMOptions::delim_cov)
        .def_readwrite("bgen_file", &GEMOptions::bgen_file)
        .def_readwrite("sample_file", &GEMOptions::sample_file)
        .def_readwrite("do_filters", &GEMOptions::do_filters)
        .def_readwrite("use_sample_file",  &GEMOptions::use_sample_file)
        .def_readwrite("includeVariantFile", &GEMOptions::includeVariantFile)
        .def_readwrite("stream_snps",  &GEMOptions::stream_snps)
        .def_readwrite("sampleid_header_name", &GEMOptions::sampleid_header_name)
        .def_readwrite("random_slope_header_name", &GEMOptions::random_slope_header_name)
        .def_readwrite("covariates", &GEMOptions::covariates)
        .def_readwrite("exposures", &GEMOptions::exposures)
        .def_readwrite("interactions", &GEMOptions::interactions)
        .def_readwrite("missing_key", &GEMOptions::missing_key)
        .def_readwrite("threads", &GEMOptions::threads)
        .def_readwrite("out_file", &GEMOptions::out_file);

    // 🔗 Bind the function that uses GEMOptions
    m.def("run_bgen", &run_bgen,
        pybind11::arg("opts"),
        "Run dosage processing pipeline using GEMOptions object");
}
