#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "RunPipeline.h"
PYBIND11_MODULE(mygen, m) 
{
    m.doc() = "GEM interface for BGEN dosage calculation";

    pybind11::class_<GEMOptions>(m, "GEMOptions")
        .def(pybind11::init<>())
        .def_readwrite("pheno_file", &GEMOptions::pheno_file)
        .def_readwrite("cov_file", &GEMOptions::cov_file)
        .def_readwrite("bgen_file", &GEMOptions::bgen_file)
        .def_readwrite("sample_file", &GEMOptions::sample_file)
        .def_readwrite(" do_filters", &GEMOptions::do_filters)
        .def_readwrite("use_sample_file",  &GEMOptions::use_sample_file)
        .def_readwrite("stream_snps",  &GEMOptions::stream_snps)
        .def_readwrite("sampleid_header_name", &GEMOptions::sampleid_header_name)
        .def_readwrite("random_slope_header_name", &GEMOptions::random_slope_header_name)
        .def_readwrite("missing_key", &GEMOptions::missing_key)
        .def_readwrite("out_file", &GEMOptions::out_file)
        .def_readwrite("threads", &GEMOptions::threads)
        .def_readwrite("delim_pheno", &GEMOptions::delim_pheno)
        .def_readwrite("delim_cov", &GEMOptions::delim_cov)
        .def_readwrite("covariates", &GEMOptions::covariates)
        .def_readwrite("exposures", &GEMOptions::exposures)
        .def_readwrite("interactions", &GEMOptions::interactions);

    // 🔗 Bind the function that uses GEMOptions
    m.def("run_from_options", &run_bgen,
        pybind11::arg("opts"),
        "Run dosage processing pipeline using GEMOptions object");
}
