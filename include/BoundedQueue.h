#pragma once

#include <cstddef>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <deque>

struct Chunk {
    std::shared_ptr<float> data; // contiguous buffer of size rows*cols
    std::size_t rows = 0;
    std::size_t cols = 0;
};

class BoundedChunkQueue {
public:
    explicit BoundedChunkQueue(std::size_t capacity);
    BoundedChunkQueue(const BoundedChunkQueue&) = delete;
    BoundedChunkQueue& operator=(const BoundedChunkQueue&) = delete;

    bool push(Chunk&& item);
    bool pop(Chunk& out);
    void close();
    bool closed() const;
    std::size_t size() const;
    std::size_t capacity() const { return cap_; }

private:
    const std::size_t cap_;
    mutable std::mutex m_;
    std::condition_variable not_empty_;
    std::condition_variable not_full_;
    bool closed_ = false;
    std::deque<Chunk> q_;
};
