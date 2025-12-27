#include "ReadBGEN.h"
#include <thread>
#include <atomic>
#include <limits>
#include <memory>
#include <cstring>
#include <mutex>
#include <condition_variable>

/**************************************
Helper function to remove extra bytes
*************************************/
auto trim_null = [](const char* s, size_t maxlen) -> std::string {
    return std::string(s, strnlen(s, maxlen));
};
/**************************************
This function is revised based on the Parse function in BOLT-LMM v2.3 source code
*************************************/

/**
 * @brief process Bgen Header Block
 * 
 * @param bgenfile 
 */
void Bgen::process_bgen_header_block(std::string bgenfile) 
{
    char genofile[300];
    strcpy(genofile, bgenfile.c_str());
    std::string genopath(genofile);
    if (genopath.substr(genopath.length() - 5, 5) != ".bgen") 
    {
        std::cout << "\nERROR: " << genopath << " does not have a .bgen extension. \n\n";
        exit(1);
    }

    fin = fopen(genofile, "rb");
    if (fin == 0) 
    {
        std::cerr << "\nERROR: BGEN file could not be opened.\n\n";
        exit(1);
    }

    std::cout << "General information of BGEN file: \n";
    if (!fread(&offset, 4, 1, fin)) 
    {
        std::cerr << "\nERROR: Cannot read BGEN header (offset).\n\n";
        exit(1);
    }

    uint L_H; 
    if (!fread(&L_H, 4, 1, fin)) 
    {
        std::cerr << "\nERROR: Cannot read BGEN header (LH).\n\n";
        exit(1);
    }

    if (fread(&Mbgen, 4, 1, fin)) 
    { 
        if (Mbgen <= 0) {
            std::cerr << "\nERROR: The number of variants in the BGEN file is 0.\n\n";
            exit(1);
        }
        std::cout << "Number of variants: " << Mbgen << '\n';
    } 
    else 
    { 
        std::cerr << "\nERROR: Cannot read BGEN header (M).\n\n"; 
        exit(1);
    }

    if (fread(&Nbgen, 4, 1, fin)) 
    {
        if (Nbgen <= 0) 
        {
            std::cerr << "\nERROR: The number of samples in the BGEN file is 0.\n\n";
            exit(1);
        }
        std::cout << "Number of samples: " << Nbgen << '\n';
    }
    else 
    {
        std::cerr << "\nERROR: Cannot read BGEN header (N). \n\n";
        exit(1);
    }


    char magic[5]; 
    if (!fread(magic, 1, 4, fin)) 
    { 
        std::cerr << "\nERROR: Cannot read BGEN header (magic bytes). \n\n"; 
        exit(1);
    }
    magic[4] = '\0';
    if (!(magic[0] == 'b' && magic[1] == 'g' && magic[2] == 'e' && magic[3] == 'n')) 
    {
        std::cerr << "\nERROR: BGEN file's four magic number bytes does not match 'b' 'g' 'e' 'n'.\n\n";
        exit(1);
    }

    fseek(fin, L_H - 20, SEEK_CUR);         

    uint flags; 
    if (!fread(&flags, 4, 1, fin)) 
    {
        std::cerr << "\nERROR: Cannot read BGEN header (flags). \n\n";
        exit(1); 
    }


    // The header block - flag definitions
    CompressedSNPBlocks = flags & 3;
    
    switch (CompressedSNPBlocks)
    {
    case 0:
        std::cout << "Genotype Block Compression Type: Uncompressed\n";
        break;
    case 1:
        std::cout << "Genotype Block Compression Type: Zlib\n";
        break;
    case 2:
        std::cout << "Genotype Block Compression Type: Zstd\n";
        break;
    default:
        std::cout << "\nERROR: BGEN compression flag must be 0 (uncompressed), 1 (zlib compression), or 2 (zstd compression). Value in file: " << CompressedSNPBlocks << ".\n\n";
        exit(1);
    }

    Layout = (flags >> 2) & 0xf; 
    std::cout << "Layout: " << Layout << '\n';
    if (Layout != 1U && Layout != 2U) 
    {
        std::cerr << "\nERROR: BGEN layout flag must be 1 or 2.\n\n";
        exit(1);
    }

    SampleIdentifiers = flags >> 31; 
    std::cout << "Sample Identifiers Present: ";  
    SampleIdentifiers == 0 ? std::cout << "False \n" : std::cout << "True \n";
    if (SampleIdentifiers != 0 && SampleIdentifiers != 1) 
    {
        std::cerr << "\nERROR: BGEN sample identifier flag must be 0 or 1.\n\n";
        exit(1);
    }
    
}


/**
 * @brief Calculate snp dosage per sample 
 * 
 * @param bgnefile
 * @param Bgen
 * @param queue
 * @param snp_per_chunk
 */
// Orchestrator: determine thread count from queue capacity and stream batched chunks without extra copies
void calc_dosage(const std::string& bgenFile, Bgen &bgen, BoundedChunkQueue& queue, int snps_per_chunk)
{
    const int samSize = bgen.new_samSize;
    if (snps_per_chunk <= 0) snps_per_chunk = 1000;

    int threads = static_cast<int>(queue.capacity());
    if (threads <= 0) threads = 1;

    // Partition work across threads using existing helper
    bgen.get_position_bgen_variant(threads, std::string(), bgen.filterVariants);

    std::atomic<bool> stop{false};
    // Enforce ordered pushing: only thread `next_thread_to_push` may push.
    std::atomic<int> next_thread_to_push{0};
    std::mutex order_mtx;
    std::condition_variable order_cv;
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(threads));

    for (int t = 0; t < threads; ++t) 
    {
        workers.emplace_back([&, t]() 
        {
            double MAF = 0.001;
            double maxMAF = 1 - MAF;
            auto start_time = std::chrono::high_resolution_clock::now();

            const uint Nbgen  = bgen.Nbgen;
            const uint Layout = bgen.Layout;
            const uint CompressedSNPBlocks = bgen.CompressedSNPBlocks;
            const bool filterVariants = bgen.filterVariants;
            const std::vector<long int>& include_idx = bgen.include_idx;

            constexpr uint maxLA = 65536;
            std::vector<char> snpID(maxLA + 1), rsID(maxLA + 1), chrStr(maxLA + 1), allele1(maxLA + 1), allele0(maxLA + 1);
            // Work buffers
            std::vector<uchar> zBuf, shortBuf, zBuf1;
            std::vector<uint16_t> shortBuf1;
            const uLongf destLen1 = 6 * Nbgen; // for Layout 1
            if (Layout == 1) 
            {
                if (CompressedSNPBlocks == 0) zBuf1.resize(destLen1); else shortBuf1.resize(destLen1);
            }
            
            // Decompressors and file
            struct DecompDel 
            { 
                void operator()(libdeflate_decompressor* p) const noexcept 
                { 
                    if (p) libdeflate_free_decompressor(p); 
                } 
            };
            std::unique_ptr<libdeflate_decompressor, DecompDel> decompressor(libdeflate_alloc_decompressor());
            std::unique_ptr<FILE, decltype(&fclose)> fin(fopen(bgenFile.c_str(), "rb"), &fclose);
            if (!fin) { stop = true; return; }
            
            // Seek to thread's start
            fseek(fin.get(), static_cast<long>(bgen.bgenVariantPos[t]), SEEK_SET);
            uint snploop = bgen.Mbgen_begin[t];
            const uint end = bgen.Mbgen_end[t];
            int keepIndex = 0;
            int ret;
            
            // Allocate initial chunk buffer
            std::shared_ptr<float> chunk_buf(new (std::nothrow) float[static_cast<size_t>(snps_per_chunk) * static_cast<size_t>(samSize)], std::default_delete<float[]>());
            if (!chunk_buf) { stop = true; return; }
            int row_in_chunk = 0;
            const float nanv = std::numeric_limits<float>::quiet_NaN();
            
            Chunk current_chunk;
            current_chunk.data = chunk_buf;
            current_chunk.cols = static_cast<std::size_t>(samSize);
            current_chunk.snpid.reserve(snps_per_chunk);
            current_chunk.rsid.reserve(snps_per_chunk);
            current_chunk.chr.reserve(snps_per_chunk);
            current_chunk.pos.reserve(snps_per_chunk);
            current_chunk.allele1.reserve(snps_per_chunk);
            current_chunk.allele0.reserve(snps_per_chunk);
            current_chunk.n_samples.reserve(snps_per_chunk);

            while (!stop && snploop <= end) 
            {
                double gsqmean = 0;
                int stream_i = 0;
                if (Layout == 1) 
                {
                    uint Nrow; ret = fread(&Nrow, 4, 1, fin.get());
                    if (Nrow != Nbgen) { stop = true; break; }
                }
                
                ushort LS; ret = fread(&LS, 2, 1, fin.get()); ret = fread(snpID.data(), 1, LS, fin.get()); snpID[LS] = '\0';
                ushort LR; ret = fread(&LR, 2, 1, fin.get()); ret = fread(rsID.data(), 1, LR, fin.get()); rsID[LR] = '\0';
                ushort LC; ret = fread(&LC, 2, 1, fin.get()); ret = fread(chrStr.data(), 1, LC, fin.get()); chrStr[LC] = '\0';
                uint32_t physpos; 
                std::vector <double> AF(snps_per_chunk);
                std::string physpos_tmp;
                ret = fread(&physpos, 4, 1, fin.get()); 
                physpos_tmp = std::to_string(physpos);
                
                if (Layout == 2) {
                    uint16_t LKnum; ret = fread(&LKnum, 2, 1, fin.get());
                    if (LKnum != 2) { stop = true; break; }
                }
                
                uint32_t LA; ret = fread(&LA, 4, 1, fin.get()); ret = fread(allele1.data(), 1, LA, fin.get()); allele1[LA] = '\0';
                uint32_t LB; ret = fread(&LB, 4, 1, fin.get()); ret = fread(allele0.data(), 1, LB, fin.get()); allele0[LB] = '\0';
                
                // Write directly into the current row
                float* row_ptr = chunk_buf.get() + static_cast<size_t>(row_in_chunk) * static_cast<size_t>(samSize);
                std::size_t out_cols = 0;
                uint nmiss = 0;

                if (Layout == 1) 
                {
                    uint16_t* probs_start;
                    if (CompressedSNPBlocks == 1) 
                    {
                        uint zLen; ret = fread(&zLen, 4, 1, fin.get()); zBuf1.resize(zLen);
                        ret = fread(&zBuf1[0], 1, zLen, fin.get());
                        if (libdeflate_zlib_decompress(decompressor.get(), &zBuf1[0], zLen, &shortBuf1[0], destLen1, NULL) != LIBDEFLATE_SUCCESS) { stop = true; break; }
                        probs_start = &shortBuf1[0];
                    } 
                    else 
                    {
                        ret = fread(&zBuf1[0], 1, destLen1, fin.get()); probs_start = reinterpret_cast<uint16_t*>(&zBuf1[0]);
                    }
                    const double scale = 1.0 / 32768; int idx_k = 0;
                    for (uint i = 0; i < Nbgen && idx_k < samSize; i++) 
                    {
                        if (include_idx[idx_k] == static_cast<long int>(i)) 
                        {
                            double p11 = probs_start[3 * i] * scale, p10 = probs_start[3 * i + 1] * scale, p00 = probs_start[3 * i + 2] * scale;
                            
                            if (!(p11 == 0 && p10 == 0 && p00 == 0)) 
                            {
                                double pTot = p11 + p10 + p00;
                                auto dosage = static_cast<float>((2 * p00 + p10) / pTot);
                                row_ptr[out_cols++] = dosage;
                                AF[stream_i] += dosage;
                                gsqmean += dosage * dosage;
                            } 
                            else 
                            {
                                row_ptr[out_cols++] = nanv;
                                nmiss++;
                            }
                            idx_k++;
                        }
                    }
                } 
                else 
                { 
                    // Layout 2
                    uint zLen; ret = fread(&zLen, 4, 1, fin.get());
                    // Filtering semantics: use +1 offset like the original calc_dosage
                    if (filterVariants && bgen.keepVariants.size() > static_cast<size_t>(t) && keepIndex < static_cast<int>(bgen.keepVariants[t].size()) && bgen.keepVariants[t][keepIndex] + 1 != snploop) 
                    {
                        // skip block payload
                        if (CompressedSNPBlocks > 0) fseek(fin.get(), 4 + zLen - 4, SEEK_CUR); else fseek(fin.get(), zLen, SEEK_CUR);
                        snploop++;
                        continue;
                    }

                    uint DLen; uchar* bufAt;
                    if (CompressedSNPBlocks == 1) 
                    {
                        zBuf.resize(zLen - 4);
                        ret = fread(&DLen, 4, 1, fin.get()); 
                        ret = fread(&zBuf[0], 1, zLen - 4, fin.get());
                        shortBuf.resize(DLen); 
                        uLongf destLen = DLen;
                        if (libdeflate_zlib_decompress(decompressor.get(), &zBuf[0], zLen - 4, &shortBuf[0], destLen, NULL) != LIBDEFLATE_SUCCESS) { stop = true; break; }
                        bufAt = &shortBuf[0];
                    } 
                    else if (CompressedSNPBlocks == 2) 
                    {
                        zBuf.resize(zLen - 4);
                        ret = fread(&DLen, 4, 1, fin.get()); 
                        ret = fread(&zBuf[0], 1, zLen - 4, fin.get());
                        shortBuf.resize(DLen); 
                        uLongf destLen = DLen;
                        size_t dret = ZSTD_decompress(&shortBuf[0], destLen, &zBuf[0], zLen - 4);
                        if (ZSTD_isError(dret)) { stop = true; break; }
                        bufAt = &shortBuf[0];
                    } 
                    else 
                    {
                        zBuf.resize(zLen); ret = fread(&zBuf[0], 1, zLen, fin.get()); bufAt = &zBuf[0];
                    }

                    uint32_t N; std::memcpy(&N, bufAt, sizeof(int32_t)); if (N != Nbgen) { stop = true; break; }
                    uint16_t K; std::memcpy(&K, &(bufAt[4]), sizeof(int16_t)); if (K != 2U) { stop = true; break; }
                    const uint32_t min_ploidy = bufAt[6]; if (min_ploidy != 2U) { stop = true; break; }
                    const uint32_t max_ploidy = bufAt[7]; if (max_ploidy != 2U) { stop = true; break; }

                    const unsigned char* missing_and_ploidy_info = &(bufAt[8]);
                    const unsigned char* probs_start = &(bufAt[10 + N]);
                    const uint32_t is_phased = probs_start[-2]; if (is_phased != 1 && is_phased != 0) { stop = true; break; }
                    const uint32_t bit_precision = probs_start[-1]; if (bit_precision != 8 && bit_precision != 16 && bit_precision != 24 && bit_precision != 32) { stop = true; break; }
                    const uintptr_t numer_mask = (1U << bit_precision) - 1; const uintptr_t probs_offset = bit_precision / 8;

                    int idx_k = 0;
                    if (!is_phased) 
                    {
                        for (uint32_t i = 0; i < N && idx_k < samSize; i++) 
                        {
                            const uint32_t mp = missing_and_ploidy_info[i]; uintptr_t numer_aa, numer_ab;
                            if (mp == 2) 
                            { 
                                bgen13_get_two_vals(probs_start, bit_precision, probs_offset, &numer_aa, &numer_ab); probs_start += (probs_offset * 2); 
                            }
                            else if (mp == 130) 
                            { 
                                if (include_idx[idx_k] == static_cast<long int>(i)) 
                                { 
                                    row_ptr[out_cols++] = nanv; idx_k++; 
                                } 
                                probs_start += (probs_offset * 2); 
                                continue; 
                            }
                            else 
                            { 
                                stop = true; 
                                break; 
                            }
                            if (include_idx[idx_k] == static_cast<long int>(i)) 
                            { 
                                double p11 = numer_aa / double(numer_mask); 
                                double p10 = numer_ab / double(numer_mask); 
                                double dosage = static_cast<float>(2 * (1 - p11 - p10) + p10); 
                                AF[stream_i] += dosage;
                                gsqmean += dosage * dosage;
                                row_ptr[out_cols++] = dosage;
                                idx_k++; 
                            }
                        }
                    } 
                    else 
                    {
                        for (uint32_t i = 0; i < N && idx_k < samSize; i++) 
                        {
                            const uint32_t mp = missing_and_ploidy_info[i]; uintptr_t numer_aa, numer_ab;
                            if (mp == 2) { bgen13_get_two_vals(probs_start, bit_precision, probs_offset, &numer_aa, &numer_ab); probs_start += (probs_offset * 2); }
                            else if (mp == 130) 
                            { 
                                if (include_idx[idx_k] == static_cast<long int>(i)) 
                                { 
                                    row_ptr[out_cols++] = nanv; idx_k++; 
                                } 
                                probs_start += (probs_offset * 2); 
                                continue; 
                            }
                            else { stop = true; break; }
                            if (include_idx[idx_k] == static_cast<long int>(i)) 
                            { 
                                double p11 = numer_aa / double(numer_mask); 
                                double p10 = numer_ab / double(numer_mask); 
                                
                                double dosage = static_cast<float>(2 - (p11 + p10)); 
                                row_ptr[out_cols++] = dosage;
                                AF[stream_i] += dosage;
                                gsqmean += dosage * dosage;
                                idx_k++; 
                            }
                        }
                    }
                    if (filterVariants) keepIndex++;
                }

                double gmean  = AF[stream_i] / double(samSize - nmiss);
                gsqmean /= static_cast<double>(samSize - nmiss);
                double cur_AF = gmean / 2.0 ;
                double gvar = (gsqmean - gmean * gmean) * static_cast<double>(samSize - nmiss) / static_cast<double>(samSize - nmiss - 1);

                if ((cur_AF < MAF || cur_AF > maxMAF) ) 
                { 
                    AF[stream_i] = 0.0;
                    snploop++;
                    stream_i++;
                    continue;
                }
                else 
                {
                    AF[stream_i] = cur_AF;
                }

                // Pad remaining columns with NaN to fixed width
                for (std::size_t c = out_cols; c < static_cast<std::size_t>(samSize); ++c) row_ptr[c] = nanv;
                row_in_chunk++;
                
                current_chunk.snpid.push_back(
                    LS > 0 ? trim_null(snpID.data(), 65536) : "NA"
                );
                current_chunk.rsid.push_back(
                    LR > 0 ? trim_null(rsID.data(), 65536) : "NA"
                );
                current_chunk.chr.push_back(
                    LC > 0 ? trim_null(chrStr.data(), 65536) : "-1"
                );
                current_chunk.pos.push_back(physpos_tmp);
                current_chunk.allele1.push_back(
                    trim_null(allele1.data(), 65536)
                );
                current_chunk.allele0.push_back(
                    trim_null(allele0.data(), 65536)
                );
                current_chunk.n_samples.push_back(std::to_string(samSize - nmiss));
                current_chunk.af.push_back(cur_AF);
                current_chunk.gv.push_back(gvar);


                if (row_in_chunk == snps_per_chunk) 
                {
                    // Wait until it's this thread's turn to push
                    {
                        std::unique_lock<std::mutex> lk(order_mtx);
                        order_cv.wait(lk, [&]{ return stop || next_thread_to_push.load() == t; });
                    }
                
                    current_chunk.data = chunk_buf; 
                    current_chunk.rows = row_in_chunk; 
                    current_chunk.cols = static_cast<std::size_t>(samSize);
                    if (!queue.push(std::move(current_chunk))) { stop = true; break; }
                    chunk_buf.reset(new (std::nothrow) float[static_cast<size_t>(snps_per_chunk) * static_cast<size_t>(samSize)], std::default_delete<float[]>());
                    if (!chunk_buf) { stop = true; break; }
                    row_in_chunk = 0;

                    current_chunk.snpid.clear();
                    current_chunk.rsid.clear();
                    current_chunk.chr.clear();
                    current_chunk.pos.clear();
                    current_chunk.allele0.clear();
                    current_chunk.allele1.clear();
                    current_chunk.n_samples.clear();
                    current_chunk.af.clear();
                    current_chunk.gv.clear();
                }
                
                snploop++;
                stream_i++;
            }

            // Flush remaining rows
            if (!stop && row_in_chunk > 0) 
            {
                {
                    std::unique_lock<std::mutex> lk(order_mtx);
                    order_cv.wait(lk, [&]{ return stop || next_thread_to_push.load() == t; });
                }
                
                current_chunk.data = chunk_buf; 
                current_chunk.rows = row_in_chunk; 
                current_chunk.cols = static_cast<std::size_t>(samSize);
                (void)queue.push(std::move(current_chunk));
            }

            // Handoff to next thread when this thread finished pushing
            {
                std::unique_lock<std::mutex> lk(order_mtx);
                // If this thread produced no chunks, it still must wait for its turn
                order_cv.wait(lk, [&]{ return stop || next_thread_to_push.load() == t; });
                next_thread_to_push.fetch_add(1);
            }
            order_cv.notify_all();

            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            std::cout << "Thread " << t << " finished in ";
            std::cout << "Elapsed time: " << duration.count() << " ms\n";
        });
    }

    for (auto& th : workers) th.join();
    queue.close();
}



void Bgen13GetTwoVals(const unsigned char* prob_start, uint32_t bit_precision, uintptr_t offset, uintptr_t* first_val_ptr, uintptr_t* second_val_ptr) {

    switch (bit_precision) {
    case 8:
        *first_val_ptr = prob_start[0];
        prob_start += offset;
        *second_val_ptr = prob_start[0];
        break;
    case 16:
        *first_val_ptr = prob_start[0] | (prob_start[1] << 8);
        prob_start += offset;
        *second_val_ptr = prob_start[0] | (prob_start[1] << 8);
        break;
    case 24:
        *first_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16);
        prob_start += offset;
        *second_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16);
        break;
    case 32:
        *first_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16) | (prob_start[3] << 24);
        prob_start += offset;
        *second_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16) | (prob_start[3] << 24);
        break;
    }

}


/**
 * @brief This functions reads the sample block of BGEN v1.1, v1.2, 
 * and v1.3. Also finds which samples to remove if they have missing values in the covariate file.
 * 
 * @param samplefile 
 * @param useSample 
 * @param covmap 
 * @param MissingKey 
 * @param numSelCol 
 * @param samSize 
 */
/**********************************************************************************
This function is revised based on the Parse function in BOLT-LMM v2.3 source code
***********************************************************************************/

// This functions reads the sample block of BGEN v1.1, v1.2, and v1.3. Also finds which samples to remove if they have missing values in the covariate file.
void Bgen::process_bgen_sample_block(const char samplefile[300], bool useSample, std::ext::UMap_str_VV_string covmap, std::string MissingKey, int numSelCol, int samSize) 
{
    int k = 0;
    std::unordered_set<int> genoUnMatchID;
    std::ext::V_string tempID;

    std::ext::V_double new_covdata_orig(samSize * (numSelCol+1));
    if ((SampleIdentifiers == 0) || useSample) 
    {
        if (SampleIdentifiers == 0 && !useSample) {
            std::cerr << "\nERROR: BGEN file does not contain sample identifiers. A .sample file is required. \n"
                << "       See https://www.well.ox.ac.uk/~gav/qctool/documentation/sample_file_formats.html for .sample file format. \n\n";
            exit(1);
        }

        std::ifstream fIDMat;
        fIDMat.open(samplefile);
        if (!fIDMat.is_open()) {
            std::cerr << "\nERROR: Sample file could not be opened.\n\n";
            exit(1);
        }

        std::string IDline;
        std::getline(fIDMat, IDline);
        std::getline(fIDMat, IDline);
        uint nSamples = 0;
        while (getline(fIDMat, IDline)) {
            nSamples++;
        }
        if (nSamples != Nbgen) {
            std::cout << "\nERROR: Number of sample identifiers in .sample file (" << nSamples << ") does not match the number of samples specified in BGEN file (" << Nbgen << ").\n\n";
            exit(1);
        }
        else {
            fIDMat.clear();
            fIDMat.seekg(0, fIDMat.beg);
            getline(fIDMat, IDline);
            getline(fIDMat, IDline);
        }

        for (uint m = 0; m < Nbgen; m++) 
        {
            // IDMatching
            getline(fIDMat, IDline);
            std::istringstream iss(IDline);
            std::string strtmp;
            iss >> strtmp;
            //AllsampleIDs before matching
            sampleID_all.push_back(strtmp);
            int itmp = k;
            
            if (covmap.find(strtmp) != covmap.end()) 
            {
                auto& tmp_valvecs = covmap[strtmp]; 
                bool ID_added = false;

                for (const auto& tmp_valvec : tmp_valvecs) {
                    // Check for missing covariate values in the current vector
                    if (find(tmp_valvec.begin(), tmp_valvec.end(), MissingKey) == tmp_valvec.end() &&
                        find(tmp_valvec.begin(), tmp_valvec.end(), "") == tmp_valvec.end()) 
                    {     
                        new_covdata_orig[k * (numSelCol + 1)] = 1.0;
                        for (int c = 0; c < numSelCol; c++) 
                        {
                            sscanf(tmp_valvec[c].c_str(), "%lf", &new_covdata_orig[k * (numSelCol + 1) + c + 1]);
                        }
                        // sampleID.push_back(strtmp);
                        // k++;
                        if(!ID_added)
                        {   
                            sampleID.push_back(strtmp);
                            k++;
                            ID_added = true;
                        }
                    }
                }
            } 
            // save the index with unmatched ID into genoUnMatchID.
            if (itmp == k) 
            {
                genoUnMatchID.insert(m);
            }
        } 

        fIDMat.close();
    }

    if ((SampleIdentifiers == 1) && !useSample) {

        uint maxLA = 65536;
        char* samID = new char[maxLA + 1];

        uint LS1;  
        if (!fread(&LS1, 4, 1, fin)) {
            std::cerr << "\nERROR: Cannot read BGEN sample block (LS).\n\n";
            exit(1);
        }

        uint Nrow; 
        if (!fread(&Nrow, 4, 1, fin)) {
            std::cerr << "\nERROR: Cannot read BGEN sample block (N).\n\n";
            exit(1);
        }
        if (Nrow != Nbgen) {
            std::cerr << "\nERROR: Number of sample identifiers (" << Nrow << ") does not match number of samples specified in BGEN file (" << Nbgen << ").\n\n";
            exit(1);
        }


        for (uint m = 0; m < Nbgen; m++) {
            ushort LSID; 
            if (!fread(&LSID, 2, 1, fin)) { 
                std::cerr << "\nERROR: Cannot read BGEN sample block (LSID).\n\n"; 
                exit(1); 
            }
            if (!fread(samID, 1, LSID, fin)) {
                std::cerr << "\nERROR: Cannot read BGEN sample block (sample id).\n\n";
                exit(1);
            }
            samID[LSID] = '\0';

            std::string strtmp(samID);
            sampleID_all.push_back(strtmp);
            int itmp = k;
            
            if (m < 5) {
                tempID.push_back(strtmp);
            }

            if (covmap.find(strtmp) != covmap.end()) 
            {
                auto& tmp_valvecs = covmap[strtmp]; 
                bool ID_added = false;

                for (const auto& tmp_valvec : tmp_valvecs) 
                {
                    // Check for missing covariate values in the current vector
                    if (find(tmp_valvec.begin(), tmp_valvec.end(), MissingKey) == tmp_valvec.end() &&
                        find(tmp_valvec.begin(), tmp_valvec.end(), "") == tmp_valvec.end()) 
                    {     
                        new_covdata_orig[k * (numSelCol + 1)] = 1.0;
                        for (int c = 0; c < numSelCol; c++) 
                        {
                            sscanf(tmp_valvec[c].c_str(), "%lf", &new_covdata_orig[k * (numSelCol + 1) + c + 1]);
                        }
                        if(!ID_added)
                        {   
                            sampleID.push_back(strtmp);
                            k++;
                            ID_added = true;
                        }
                    }
                }
            }

            if (itmp == k) 
            {
                genoUnMatchID.insert(m);
            }
        }
        delete[] samID;
    } // end SampleIdentifiers == 1


    // After IDMatching, resizing covdata and covdata, and updating samSize;

    new_covdata_orig.resize(k * (numSelCol + 1));
    samSize = k;

    if (samSize == 0) 
    {
        std::cerr << "\nERROR: Sample size changed from " << samSize + genoUnMatchID.size() << " to " << samSize << ".\n\n";
        if (SampleIdentifiers == 1 && !useSample) {
            int print_i = 5;
            if (Nbgen < 5)
            { 
                print_i = Nbgen; 
            }
            std::cout << "ID matching was done using the BGEN sample identifier block. \nHere are the first " << print_i << " sample identifiers in BGEN file: \n";
            for (int i = 0; i < print_i; i++) 
            {
                std::cout << " " << tempID[i] << "\n";
            }
        }

        if (SampleIdentifiers == 0 || useSample) 
        {
            std::cout << "Check if sample IDs are consistent between the covariate file and sample file, or check if (--sampleid-name) is specified correctly. \n\n";
        }
        exit(1);
    }


    int ii = 0;
    include_idx.resize(samSize);
    for (uint i = 0; i < Nbgen; i++) {
        if (genoUnMatchID.find(i) == genoUnMatchID.end()) 
        {
             include_idx[ii] = i;
              ii++;
        }
    }

    std::cout << "****************************************************************************\n";
    if (genoUnMatchID.empty()) 
    {
        std::cout << "After processes of sample IDMatching and checking missing values, the sample size does not change.\n\n";
    }
    else 
    {
        std::cout << "After processes of sample IDMatching and checking missing values, the sample size changes from "
            << samSize + genoUnMatchID.size() << " to " << samSize << ".\n\n";
    }
    std::cout << "Sample IDMatching and checking missing values processes have been completed.\n";
    std::cout << "****************************************************************************\n";


    new_samSize = samSize;
    if (new_samSize<(numSelCol+1) || new_samSize == (numSelCol+1))
    {
        std::cout << "\nERROR: The sample size should be greater than the number of predictors!" <<std::endl;
        exit(1);
    }

    // //the first column of matcovX is Y
    // MatrixXd matcovX (samSize,(numSelCol+1));
    // for (int i=0; i<samSize; i++){    
    //     for (int j=0; j<(numSelCol+1); j++) {
    //       matcovX(i,j) =new_covdata_orig [i * (numSelCol+1) +j];
    //     }
    // }
    // Eigen::HouseholderQR<MatrixXd> qr;
    // qr.compute(matcovX);
    // Eigen::MatrixXd R = qr.matrixQR();
    // int colR=R.cols();
    // VectorXd diagR (colR);
    // for (int i=0; i<colR; i++){
    //     diagR(i)=abs(R(i,i));
    // }

    // double sqrtEps =sqrt(std::numeric_limits<double>::epsilon());
    // double maxdiag = *std::max_element( diagR.begin(), diagR.end() ) ;
    // double colinear_cut = abs(maxdiag * sqrtEps);
    // for (int i=0; i<colR; i++){
    //     if (abs(diagR(i)) < colinear_cut){
    //         excludeCol.push_back(i);    
    //     }
    // }
    // matcovX.resize(0,0);
    // R.resize(0,0);

    // int NumExcludeCol = excludeCol.size();
    // if (excludeCol.size()>0){        
    //     std::vector <int> remove_colinear;
    //     for (int i=0; i<excludeCol.size(); i++){
    //         for (int j=0; j<samSize; j++) {
    //             remove_colinear.push_back(j * (numSelCol+1) + excludeCol[i]);
    //         }
    //     }

    //     numSelCol=numSelCol- excludeCol.size();
    //     new_covdata.resize(samSize * (numSelCol+1));
    //     std::vector<double> temp;
    //     for (int i=0; i<new_covdata_orig.size(); i++)
    //     {
    //         if (std::find(remove_colinear.begin(), remove_colinear.end(), i) == remove_colinear.end())
    //         {
    //             temp.push_back(new_covdata_orig[i]);
                
    //         }
    //     }
    //     new_covdata = temp;
    // } 
    // else 
    // {
            // new_covdata.resize(samSize * (numSelCol+1));
            // new_covdata = new_covdata_orig;
    // }

}



/***********************************************************************************
This function contains code that is revised based on BOLT-LMM v2.3 source code
************************************************************************************/

// This function reads just the variant block for BGEN files version v1.1, v1.2, and v1.3 and is used to grab the byte where the variant begins.
//    Necesary when there's no bgen index file.
/**
 * @brief This function reads just the variant block for BGEN files version v1.1, v1.2, and v1.3 and 
 * is used to grab the byte where the variant begins.
 * 
 * @param threads 
 * @param includeVariantFile 
 * @param doFilters 
 */
void Bgen::get_position_bgen_variant(int threads, std::string includeVariantFile, bool doFilters) 
{
    int count = 0;
    uint CompressedSNPBlocks = this->CompressedSNPBlocks;
    uint Layout = this->Layout;
    uint offset = this->offset;
    uint Mbgen = this->Mbgen;
    uint nSNPS = this->Mbgen;
    uint Nbgen = this->Nbgen;
    uint maxLA = 65536;
    char* snpID   = new char[maxLA + 1];
    char* rsID    = new char[maxLA + 1];
    char* chrStr  = new char[maxLA + 1];
    char* allele1 = new char[maxLA + 1];
    char* allele0 = new char[maxLA + 1];
    std::string IDline;

    std::ext::Set_string includeVariant;
    std::vector<std::vector<uint>> includeVariantIndex;
    bool checkSNPID = false;
    bool checkRSID = false;
    bool checkInclude = false;
    int ret;

    if (doFilters) 
    {
        filterVariants = true;
        if (!includeVariantFile.empty()) 
        {
            checkInclude = true;

            std::ifstream fInclude;
            fInclude.open(includeVariantFile);
            if (!fInclude.is_open()) {
                std::cerr << "\nERROR: The file (" << includeVariantFile << ") could not be opened.\n\n";
                exit(1);
            }

            std::string vars;
            getline(fInclude, IDline);
            std::transform(IDline.begin(), IDline.end(), IDline.begin(), ::tolower);
            IDline.erase(std::remove(IDline.begin(), IDline.end(), '\r'), IDline.end());
            if (IDline == "snpid") {
                std::cout << "An include snp file was detected... \nIncluding SNPs for analysis based on their snpid... \n";
                checkSNPID = true;
            }
            else if (IDline == "rsid") {
                std::cout << "An include snp file was detected... \nIncluding SNPs for analysis based on their rsid... \n";
                checkRSID = true;
            }
            else {
                std::cerr << "\nERROR: Header name of " << includeVariantFile << " must be snpid or rsid.\n\n";
                exit(1);
            }

            while (fInclude >> vars) 
            {
                if (includeVariant.find(vars) != includeVariant.end()) 
                {
                    std::cout << "\nERROR: " << vars << " is a duplicate variant in " << includeVariantFile << ".\n\n";
                    exit(1);
                }
                includeVariant.insert(vars);
                count++;
            }
            nSNPS = count;
            std::cout << "Detected " << nSNPS << " variants to be used for analysis... \nAll other variants will be excluded.\n\n\n";
            std::cout << "Detected" << std::thread::hardware_concurrency() << " available thread(s)...\n"; 
            if (nSNPS < threads) {
                std::cout << "Number of variants (" << nSNPS << ") is less than the number of specified threads (" << threads << ")...\n";
                threads = nSNPS;
                std::cout << "Dividing to " << threads << " chunks for finding SNPs positions... \n\n";
            }
            else {
                std::cout << "Dividing to " << threads << " chunks for finding SNPs positions... \n\n";
            }
        }

        std::cout << "Dividing BGEN file into " << threads << " block(s)...\n";
        std::cout << "Identifying start position of each block...\n";
        std::vector<uint> endIndex(threads);
        int nBlocks = ceil(nSNPS / threads);
        uint index = 0;
        uint k = 0;
        uint sucessCount = 0;
        Mbgen_begin.resize(threads);
        Mbgen_end.resize(threads);
        bgenVariantPos.resize(threads);
        keepVariants.resize(threads);

        for (uint t = 0; t < threads; t++) 
        {
            endIndex[t] = ((t + 1) == threads) ? nSNPS - 1 : floor(((nSNPS / threads) * (t + 1)) - 1);
        }

        FILE* fin = this->fin;
        fseek(fin, offset + 4, SEEK_SET);

        for (uint snploop = 0; snploop < Mbgen; snploop++) 
        {
            long long unsigned int prev = ftell(fin);

            uint Nrow;
            if (Layout == 1) {
                ret = fread(&Nrow, 4, 1, fin); 
                if (Nrow != Nbgen) 
                {
                    std::cerr << "\nERROR: Number of samples (" << Nrow << ") with genotype probabilities does not match number of samples specified in BGEN file (" << Nbgen << ").\n\n";
                    exit(1);
                }
            }

            ushort LS; 
            ret = fread(&LS, 2, 1, fin);

            ret = fread(snpID, 1, LS, fin);
            snpID[LS] = '\0';
            if (checkSNPID) {
                if ((checkInclude) && (includeVariant.find(snpID) != includeVariant.end())) 
                {
                    sucessCount++;
                    keepVariants[k].push_back(snploop);
                    if (index == (nBlocks * k)) 
                    {
                        Mbgen_begin[k] = snploop;
                        long long int curr = ftell(fin);
                        bgenVariantPos[k] = curr - (curr - (prev));
                    }
                    if (index == endIndex[k]) 
                    {
                        Mbgen_end[k] = snploop;
                        k++;
                        if (k == threads) {break;}
                    }
                    index++;
                }
            }

            ushort LR; 
            ret = fread(&LR, 2, 1, fin);

            ret = fread(rsID, 1, LR, fin); 
            rsID[LR] = '\0';
            if (checkRSID) {
                if (checkInclude && (includeVariant.find(rsID) != includeVariant.end())) 
                {
                    sucessCount++;
                    keepVariants[k].push_back(snploop);
                    if (index == (nBlocks * k)) 
                    {
                        Mbgen_begin[k] = snploop;
                        long long unsigned int curr = ftell(fin);
                        bgenVariantPos[k] = curr - (curr - (prev));
                    }
                    if (index == endIndex[k]) 
                    {
                        Mbgen_end[k] = snploop;
                        k++;
                        if (k == threads) {break;}
                    }
                    index++;
                }
            }

            ushort LC; 
            ret = fread(&LC, 2, 1, fin);

            ret = fread(chrStr, 1, LC, fin); 
            chrStr[LC] = '\0';

            uint32_t physpos; 
            ret = fread(&physpos, 4, 1, fin);

            uint16_t LKnum;
            if (Layout == 2) {
                ret = fread(&LKnum, 2, 1, fin);
                if (LKnum != 2) {
                    std::cerr << "\nERROR: " << std::string(snpID) << " is a non-bi-allelic variant with " << LKnum << " alleles. Please filter these variants for now.\n\n";
                    exit(1);
                }
            }

            uint32_t LA; 
            ret = fread(&LA, 4, 1, fin);
            ret = fread(allele1, 1, LA, fin); 
            allele1[LA] = '\0';

            uint32_t LB; 
            ret = fread(&LB, 4, 1, fin);
            ret = fread(allele0, 1, LB, fin); 
            allele0[LB] = '\0';


            if (Layout == 2) {
                if (CompressedSNPBlocks > 0) 
                {
                    uint zLen; 
                    ret = fread(&zLen, 4, 1, fin);
                    fseek(fin, 4 + zLen - 4, SEEK_CUR);

                }
                else 
                {
                    uint zLen; 
                    ret = fread(&zLen, 4, 1, fin);
                    fseek(fin, zLen, SEEK_CUR);
                }
            }
            else 
            {
                if (CompressedSNPBlocks == 1) 
                {
                    uint zLen; 
                    ret = fread(&zLen, 4, 1, fin);
                    fseek(fin, zLen, SEEK_CUR);

                }
                else 
                {
                    fseek(fin, 6 * Nbgen, SEEK_CUR);
                }
            }
        }

        if (sucessCount != nSNPS) 
        {
            std::cerr << "\nERROR: There are one or more SNPs in BGEN file with " << IDline << " not in " << includeVariantFile << ".\n\n";
            exit(1);
        }
    }
    else 
    {
        filterVariants = false;
        std::cout << "Detected " << std::thread::hardware_concurrency() << " available thread(s)...\n";
        if (Mbgen < threads) 
        {
            threads = Mbgen;
            std::cout << "Number of variants (" << Mbgen << ") is less than the number of specified threads (" << threads << ")...\n";
            std::cout << "dividing to " << threads << " chunks for finding SNPs positions... \n\n";
        }
        else 
        {
            std::cout << "dividing to  " << threads << " chunks for finding SNPs positions... \n\n";
        }

        std::cout << "Dividing BGEN file into " << threads << " block(s)..." << std::endl;
        Mbgen_begin.resize(threads);
        Mbgen_end.resize(threads);
        bgenVariantPos.resize(threads);
        std::cout << std::flush;
        keepVariants.resize(threads);
        std::cout << std::flush;
        for (uint t = 0; t < threads-1; t++) 
        {
            Mbgen_begin[t] = floor((Mbgen / threads) * t);
            Mbgen_end[t] = floor(((Mbgen / threads) * (t + 1)) - 1);
        }

        Mbgen_begin[threads-1] = floor((Mbgen / threads) * (threads - 1));
        Mbgen_end[threads-1] = Mbgen - 1;

        uint t = 0;
        FILE* fin = this->fin;
        fseek(fin, offset + 4, SEEK_SET);

        for (uint snploop = 0; snploop < Mbgen; snploop++) 
        {
            if (snploop == Mbgen_begin[t]) 
            {
                bgenVariantPos[t] = ftell(fin);
                t++;
                if (t == (Mbgen_begin.size())) 
                {
					break;
				}
            }

            uint Nrow;
            if (Layout == 1) 
            {
                ret = fread(&Nrow, 4, 1, fin);
                if (Nrow != Nbgen) 
                {
                    std::cerr << "\nERROR: Number of samples (" << Nrow << ") with genotype probabilities does not match number of samples specified in BGEN file (" << Nbgen << ").\n\n";
                    exit(1);
                }
            }

            ushort LS; 
            ret = fread(&LS, 2, 1, fin);
            ret = fread(snpID, 1, LS, fin); 
            snpID[LS] = '\0';
            ushort LR; 
            ret = fread(&LR, 2, 1, fin);
            ret = fread(rsID, 1, LR, fin); 
            rsID[LR] = '\0';

            ushort LC; 
            ret = fread(&LC, 2, 1, fin);
            ret = fread(chrStr, 1, LC, fin); 
            chrStr[LC] = '\0';

            uint32_t physpos; 
            ret = fread(&physpos, 4, 1, fin);

            uint16_t LKnum;
            if (Layout == 2) 
            {
                ret = fread(&LKnum, 2, 1, fin);
                if (LKnum != 2) 
                {
                    std::cerr << "\nERROR: " << std::string(snpID) << " is a non-bi-allelic variant with " << LKnum << " alleles. Please filter these variants for now.\n\n";
                    exit(1);
                }
            }

            uint32_t LA; 
            ret = fread(&LA, 4, 1, fin);
            ret = fread(allele1, 1, LA, fin); 
            allele1[LA] = '\0';

            uint32_t LB; 
            ret = fread(&LB, 4, 1, fin);
            ret = fread(allele0, 1, LB, fin); 
            allele0[LB] = '\0';
            // Seeks past the uncompressed genotype.
            if (Layout == 2) 
            {
                if (CompressedSNPBlocks > 0) 
                {
                    uint zLen;  
                    ret = fread(&zLen, 4, 1, fin);
                    ret = fseek(fin, 4 + zLen - 4, SEEK_CUR);

                }
                else 
                {
                    uint zLen; 
                    ret = fread(&zLen, 4, 1, fin);
                    ret = fseek(fin, zLen, SEEK_CUR);
                }
            }
            else 
            {
                if (CompressedSNPBlocks == 1) 
                {
                    uint zLen;  
                    ret = fread(&zLen, 4, 1, fin);
                    ret = fseek(fin, zLen, SEEK_CUR);

                }
                else 
                {
                    ret = fseek(fin, 6 * Nbgen, SEEK_CUR);
                }
            }
        }        
    }
    
    std::cout << std::flush;
    (void)ret;
    delete[] snpID;
    delete[] rsID;
    delete[] chrStr;
    delete[] allele1;
    delete[] allele0;
}

/***********************************************************************************
bgen13_get_two_vals function return probs
************************************************************************************/

/**
 * @brief Reads two encoded genotype probabilities (integers).Interprets 
 * them according to the bit_precision. Stores the results in 
 * first_val_ptr and second_val_ptr.
 * 
 * @param prob_start 
 * @param bit_precision 
 * @param offset 
 * @param first_val_ptr 
 * @param second_val_ptr 
 */
void bgen13_get_two_vals(const unsigned char* prob_start, uint32_t bit_precision, uintptr_t offset, uintptr_t* first_val_ptr, uintptr_t* second_val_ptr) {

    switch (bit_precision) {
    case 8:
        *first_val_ptr = prob_start[0];
        prob_start += offset;
        *second_val_ptr = prob_start[0];
        break;
    case 16:
        *first_val_ptr = prob_start[0] | (prob_start[1] << 8);
        prob_start += offset;
        *second_val_ptr = prob_start[0] | (prob_start[1] << 8);
        break;
    case 24:
        *first_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16);
        prob_start += offset;
        *second_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16);
        break;
    case 32:
        *first_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16) | (prob_start[3] << 24);
        prob_start += offset;
        *second_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16) | (prob_start[3] << 24);
        break;
    }

}