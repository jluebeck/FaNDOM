#include <queue>
#include <mutex>

using namespace std;

#ifndef FANDOM_THREADSAFE_QUEUE_H
#define FANDOM_THREADSAFE_QUEUE_H

template<typename Data>

class threadsafe_queue {
private:
    queue<Data> q;
    mutable mutex mt;

public:
    Data pop() {
        lock_guard<mutex> mlock(mt);
        if (q.empty()) {
            return NULL;
        }
        Data &item = q.front();
        q.pop();
        return item;
    }

    void push(const Data& val) {
        lock_guard<mutex> lock(mt);
        q.push(val);
    }

    size_t size() {
        lock_guard<mutex> mlock(mt);
        return q.size();
    }

//    bool empty() {
//        lock_guard<mutex> mlock(mt);
//        return q.empty();
//    }
};

#endif //FANDOM_THREADSAFE_QUEUE_H
