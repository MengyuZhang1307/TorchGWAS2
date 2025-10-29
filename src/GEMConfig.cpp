#include "GEMConfig.h"

GEMOptions::GEMOptions() 
{
    unsigned int hw_threads = std::thread::hardware_concurrency();
    if (hw_threads == 0) hw_threads = 1; // fallback
    threads = hw_threads;
    num_chunks = hw_threads;
}