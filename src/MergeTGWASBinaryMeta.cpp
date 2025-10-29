// #include "MergeTGWASBinaryMeta.h"

// //helper functions
// /**
//  * @brief count number of snps in file
//  * 
//  * @param filename 
//  * @return size_t 
//  */
// size_t count_snpss(std::string const& filename) 
// {
//     std::ifstream file(filename);
//     if (!file) 
//     {
//         std::cerr << "Failed to open file: " << filename << std::endl;
//         return 0;
//     }

//     std::string line;
//     size_t count = 0;

//     // Skip header line
//     if (std::getline(file, line)) 
//     {
//         // now count the rest
//         while (std::getline(file, line))
//         {
//             ++count;
//         }
//     }
//     return count;
// }


// /**
//  * @brief Read JSON metadata
//  * 
//  * @param filename 
//  * @param snp_txt_count 
//  * @return nlohmann::json 
//  */
// nlohmann::json meta_json read_json(std::string const& filename,
//                                     size_t snp_txt_count)
// {
//     std::ifstream json_file(filename);
//     if (!json_file) 
//     {
//         std::cerr << "Cannot open JSON file: " << json_path << "\n";
//         std::cout << "****************************************************************************\n";
//         exit(EXIT_FAILURE);
//     }

//     nlohmann::json meta_json;
//     json_file >> meta_json;
//     json_file.close();

//     // Extract rows from JSON
//     if (!meta_json.contains("rows")) 
//     {
//         std::cerr << "JSON missing 'rows' field.\n";
//         std::cout << "****************************************************************************\n";
//         exit(EXIT_FAILURE);
//     }

//     size_t snp_json_count = meta_json["rows"];

//     // Compare json file with GEM meta data
//     if (snp_txt_count != snp_json_count) 
//     {
//         std::cerr << "Mismatch: GEM snps count=" << snp_txt_count
//                   << " vs JSON=" << snp_json_count << "\n";
//         std::cout << "****************************************************************************\n";
//         exit(EXIT_FAILURE);
//     }
//     return meta_json;
// }

// // constructor
// MergeMetaData::MergeMetaData(GEMOptions const& user_opt) : opt(user_opt){}
// /**
//  * @brief read meta data from bgen file
//  * 
//  */
// void MergeMetaData::read_bgen_data()
// {
//     opt.find_genofile_type();
//     opt.check_kinship_usage();
//     // Step 1:  Read covariate file
//     opt.shared_cov_result = opt.read_covariate_data();
//     bgen.process_bgen_header_block(opt.geno_add);
//     bgen.process_bgen_sample_block(opt.sample_add.c_str(), opt.use_sample_file, 
//                                     opt.shared_cov_result.covMap, opt.missing_key, 
//                                     opt.shared_cov_result.numSelCol, 
//                                     opt.shared_cov_result.samSize);   
//     bgen.get_position_bgen_variant(opt.num_chunks, opt.includeVariantFile,
//                                              opt.do_filters);
//     bgen_sample_id = bgen.sampleID;
//     bgen.filterVariants = opt.do_filters;
//     std::vector<std::thread> threads;
//     if (opt.threads > 1)
//     {
//         std::cout << "Running multithreading...\n";
//         for (int t = 0; t < opt.threads; t++)
//         {
//             threads.emplace_back(&Read_bgen_file, opt.geno_add, bgen, t, opt.stream_snps, "bgen_meta.txt");
//         }
//         std::cout << "Continuing GEI test... joining threads...\n";
//         for (auto& thread : threads) 
//         {
//             thread.join(); 
//         }
//         combine_metadata();
//     }
//     else if (opt.threads == 1)
//     {
//         std::cout << "Running with single thread...\n";
//         Read_bgen_file(opt.geno_add, bgen, t, opt.stream_snps, "bgen_meta.txt");
//     }
// }

// /**
//  * @brief combine multiple output data from bgen file
//  * 
//  */
// void MergeMetaData::combine_metadata()
// {
//     cout << "Combining results... \n";
//     std::ofstream results(opt.outFile, std::ios::binary | std::ios_base::app);
//     results << "SNPID" << "\tRSID\t" : "\t" << "CHR" << "\t" << "POS" << "\t" << "Non_Effect_Allele" << "\t" << "Effect_Allele" << "\t" << "N_Samples" << "\t" << "AF" << "\t" << "GV" << "\t"; 

//     for (int t = 0; t < opt.thread; t++)
//     {
//         input = opt.output + "_bin_" + std::to_string(thread_num) + ".tmp";
//         std::ifstream(input);
//         if (input.peek() != std::ifstream::traits_type::eof())
//         {
//             results << input.rdbuf();
//         }
//         input.close();
//         input.remove();
//     }

//     results.close();
// }

// void MergeMetaData::combine_metadata_TGWAS()
// {
//     std::ifstream input("intermediate" + opt.output); //TGWAS output binary file
//     std::istringstream iss;
//     std::string line;

//     size_t snp_txt_count = count_snps("intermediate" + opt.output + ".meta.json"); //GEM output metadata
//     std::cout << "SNP count (excluding header): " << snp_txt_count << "\n";
//     nlohmann::json meta_json = read_json("intermediate" + opt.output + ".meta.json", snp_txt_count); //Json output TGWAS
//     size_t rows = meta_json["rows"];
//     size_t cols = meta_json["cols"];
//     std::string dtype = meta_json["dtype"];
//     std::ext::V_string meta_headers = meta_json["headers"].get<std::ext::V_string>();

// }


// void merge_chunk(long long start_row,
//                  long long end_row,
//                  int total_cols)
// {
//     // --- open binary and text files ---
//     std::ifstream data("results_" + opt.output, , std::ios::binary); //TGWAS output binary file
//     std::ifstream meta_txt(meta_txt_path);
//     std::ifstream (data_path, std::ios::binary);
//     std::ofstream out(out_path, std::ios::binary | std::ios::in | std::ios::out);

//     size_t meta_size = sizeof(VariantMeta);
//     size_t row_bytes = total_cols * sizeof(float);

//     // --- skip lines in meta.txt until start_row ---
//     std::string dummy;
//     for (long long i = 0; i < start_row && std::getline(meta_txt, dummy); ++i) {}

//     // --- seek data and output to correct byte offsets ---
//     data.seekg(start_row * row_bytes);
//     out.seekp(start_row * (meta_size + row_bytes));

//     std::vector<float> buffer(total_cols);
//     VariantMeta m;
//     std::string line;

//     for (long long i = start_row; i < end_row && std::getline(meta_txt, line); ++i) {
//         std::istringstream ss(line);
//         std::string snpid, allele1, allele2;
//         ss >> snpid >> m.chr >> m.pos >> allele1 >> allele2 >> m.af;

//         // copy strings into fixed-size fields
//         std::snprintf(m.snpid, sizeof(m.snpid), "%s", snpid.c_str());
//         std::snprintf(m.allele1, sizeof(m.allele1), "%s", allele1.c_str());
//         std::snprintf(m.allele2, sizeof(m.allele2), "%s", allele2.c_str());

//         data.read(reinterpret_cast<char*>(buffer.data()), row_bytes);

//         out.write(reinterpret_cast<char*>(&m), meta_size);
//         out.write(reinterpret_cast<char*>(buffer.data()), row_bytes);
//     }
// }
