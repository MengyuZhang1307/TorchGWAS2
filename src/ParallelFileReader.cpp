#include "ParallelFileReader.h"

GEMOptions::GEMOptions() 
{
    unsigned int hw_threads = std::thread::hardware_concurrency();
    if (hw_threads == 0) hw_threads = 1; // fallback
    threads = hw_threads;
    num_chunks = hw_threads;
}

/**
 * @brief Write bgen dosages to file
 * 
 * @param results 
 * @param filename 
 * @param delimiter 
 */
void write_bgen_result_to_file(const std::ext::V_bgen& results, const std::string& filename, char delimiter) 
{
    std::ofstream out(filename);
    if (!out.is_open()) 
    {
        std::cerr << "Error: could not open file " << filename << " for writing.\n";
        return;
    }

    for (const auto& thread_vec : results) 
    {
        for (const auto& snp_vec : thread_vec) 
        {
            for (const auto& sample_vec : snp_vec) 
            {
                for (size_t i = 0; i < sample_vec.size(); ++i) 
                {
                    out << sample_vec[i];
                    if (i != sample_vec.size() - 1)
                        out << delimiter;
                }
                out << '\n';  // new line for each innermost vector
            }
        }
    }

    out.close();
    std::cout << "\u2705 BGEN results written to: " << filename << '\n';
    std::cout << "****************************************************************************\n";
}

/**
 * @brief  Write all thredead's output files in one file
 * 
 * @param threads 
 * @param outfile 
 */
void output_file_generator(int threads, std::string outfile)
{
    if (fs::exists(outfile)) 
    {
        fs::remove(outfile);  // Delete file
    }

    std::ofstream results(outfile, std::ios::binary | std::ios_base::app);
    if (!results.is_open()) 
    {
        std::cerr << "Failed to open " << outfile << "\n";  // Append file content
        exit(EXIT_FAILURE);
    } 

    for (int thread = 0; thread < threads; ++thread)
    {
        std::string thread_output = outfile + std::to_string(thread) + ".tmp";
        std::ifstream thread_output_file(thread_output);
        results << thread_output_file.rdbuf();
        thread_output_file.close();
        std::remove(thread_output.c_str());
    }
    
    results.close();

    std::cout << "\u2705 BGEN results written to: " << outfile << '\n';
    std::cout << "****************************************************************************\n";
}
