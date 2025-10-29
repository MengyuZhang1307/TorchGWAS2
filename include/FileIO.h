#pragma once
#include "ReadFiles.h"
void write_bgen_result_to_file(const std::ext::V_bgen& results, const std::string& filename,
                                char delimiter = ',');

void output_file_generator(int threads, std::string outfile);
