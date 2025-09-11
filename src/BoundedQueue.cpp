#include "BoundedQueue.h"

BoundedChunkQueue::BoundedChunkQueue(std::size_t capacity)
    : cap_(capacity) {}

bool BoundedChunkQueue::push(Chunk&& item) {
    std::unique_lock<std::mutex> lk(m_);
    not_full_.wait(lk, [&]{ return q_.size() < cap_ || closed_; });
    if (closed_) return false;
    q_.emplace_back(std::move(item));
    not_empty_.notify_one();
    return true;
}

bool BoundedChunkQueue::pop(Chunk& out) {
    std::unique_lock<std::mutex> lk(m_);
    not_empty_.wait(lk, [&]{ return !q_.empty() || closed_; });
    if (q_.empty()) return false; // closed and empty
    out = std::move(q_.front());
    q_.pop_front();
    not_full_.notify_one();
    return true;
}

void BoundedChunkQueue::close() {
    std::lock_guard<std::mutex> lk(m_);
    closed_ = true;
    not_empty_.notify_all();
    not_full_.notify_all();
}

bool BoundedChunkQueue::closed() const {
    std::lock_guard<std::mutex> lk(m_);
    return closed_;
}

std::size_t BoundedChunkQueue::size() const {
    std::lock_guard<std::mutex> lk(m_);
    return q_.size();
}
