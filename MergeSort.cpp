#include <iostream>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <cstring>
#include <chrono>
#include <cilk/opadd_reducer.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#define MERGESIZE (2048)

const int64_t INIT_CHUNK = 100000;

/* Validate returns a count of the number of times the i-th number is less than (i-1)th number.
   This version is parallel using cilk_for and cilk reducer. */
int64_t ValidateParallel(const std::vector<int64_t> & input) {

    int64_t n = input.size();
    if (n <= 1) return 0;
    
    cilk::opadd_reducer<int64_t> count(0);

    cilk_for (uint64_t i = 1; i < n; i++) {
        if (input[i-1] > input[i]) {
            count += 1;
        }
    }
    return count;
}

/* SerialMerge merges two sorted arrays vec[Astart..Aend] and vec[Bstart..Bend] into vector tmp starting at index tmpIdx. */

void SerialMerge(const std::vector<int64_t> & vec,  int64_t Astart,  int64_t Aend, int64_t Bstart,
                 int64_t Bend, std::vector<int64_t> & tmp, int64_t tmpIdx){
    if (Astart <= Aend && Bstart <= Bend) {
        while(true) {
            if (vec[Astart] < vec[Bstart]) {
                tmp[tmpIdx++] = vec[Astart++];
                if (Astart > Aend) break;
            } else {
                tmp[tmpIdx++] = vec[Bstart++];
                if (Bstart > Bend) break;
            }
        }
    }
    if (Astart > Aend) {
        if (Bstart <= Bend)
            memcpy(&tmp[tmpIdx], &vec[Bstart], sizeof(int64_t) * (Bend - Bstart + 1));
    } else {
        if (Astart <= Aend)
            memcpy(&tmp[tmpIdx], &vec[Astart], sizeof(int64_t) * (Aend - Astart + 1));
    }
}

// ParallelMerge: merges two sorted arrays vec[Astart..Aend] and vec[Bstart..Bend] into tmp[tmpIdx...],
//  makes Halfway split the mid element in its correct position and spawns 2 More 'll merges 

void ParallelMerge(const std::vector<int64_t> & vec, int64_t Astart, int64_t Aend, int64_t Bstart,
                 int64_t Bend, std::vector<int64_t> & tmp, int64_t tmpIdx) {

    // Base Cases
    if (Astart > Aend) {
        if (Bstart <= Bend)
            memcpy(&tmp[tmpIdx], &vec[Bstart], sizeof(int64_t) * (Bend - Bstart + 1));
        return;
    }
    if (Bstart > Bend) {
        if (Astart <= Aend)
            memcpy(&tmp[tmpIdx], &vec[Astart], sizeof(int64_t) * (Aend - Astart + 1));
        return;
    }

    int64_t n1 = Aend - Astart + 1;
    int64_t n2 = Bend - Bstart + 1;

    // small case , jst merge serially
    if ((n1 + n2) <= MERGESIZE) {
        SerialMerge(vec, Astart, Aend, Bstart, Bend, tmp, tmpIdx);
        return;
    }

    // pick the mid ele from the larger of the two arrays to balance recursion.
    if (n1 >= n2) {

        int64_t midA = Astart + (n1 / 2);
        int64_t val = vec[midA];    // find mid value

        // find first pos in B where >= mid val
        int64_t Bk = std::lower_bound( (const int64_t*)(&vec[Bstart]), (const int64_t*)(&vec[Bend]) + 1, val )
                    - (const int64_t*)(&vec[0]);
        // std::lower_bound returns pointer so converting it to index by - &vec[0].
        
        int64_t leftCount = (midA - Astart) + (Bk - Bstart);
        int64_t pos = tmpIdx + leftCount;
        tmp[pos] = val;

       

        // here chosing vec[mid] divides the problem into 2 discreete problems , 
        // 1-> where elements are less than vec[mid]
        // 2-> where elements are grater than vec[mid]
        // we can assign them to different workers 

        cilk_spawn ParallelMerge(vec, Astart, midA - 1, Bstart, Bk - 1, tmp, tmpIdx);
        ParallelMerge(vec, midA + 1, Aend, Bk, Bend, tmp, pos + 1);
        cilk_sync;

    } else { // th mirror case 

        // n2 > n1: pick median from B (for balenced resursion)
        int64_t midB = Bstart + (n2 / 2);
        int64_t val = vec[midB];

        int64_t Ak = std::lower_bound( (const int64_t*)(&vec[Astart]), (const int64_t*)(&vec[Aend]) + 1, val )
                    - (const int64_t*)(&vec[0]);
        ;
        int64_t leftCount = (Ak - Astart) + (midB - Bstart);
        int64_t pos = tmpIdx + leftCount;
        tmp[pos] = val;

        cilk_spawn ParallelMerge(vec, Astart, Ak - 1, Bstart, midB - 1, tmp, tmpIdx);
        ParallelMerge(vec, Ak, Aend, midB + 1, Bend, tmp, pos + 1);
        cilk_sync;
    }
}
void MergeSortSerial(std::vector<int64_t> & vec, int64_t start, std::vector<int64_t> & tmp, 
                    int64_t tmpStart, int64_t sz, int64_t depth, int64_t limit){
    
    // Note: also need to use depth and limit to control the granularity in the parallel case.
    if (sz < MERGESIZE) {
        sort(vec.begin() + start,vec.begin() + start + sz);
        return;
    }
    depth++;
    
    auto half = sz >> 1;
    auto A = start;
    auto tmpA = tmpStart;
    auto B = start + half;
    auto tmpB = tmpStart + half;
    // Sort left half.
    MergeSortSerial(vec, A, tmp, tmpA, half, depth, limit);
    // Sort right half.
    MergeSortSerial(vec, B, tmp, tmpB, start + sz - B, depth, limit);
    // More for paralle case?
    
    // Merge sorted parts into tmp.
    SerialMerge(vec, A, A + half-1, B, start + sz -1, tmp, tmpStart);
    // Copy result back to vec.
    memcpy(&vec[A], &tmp[tmpStart],sizeof(int64_t) * (sz));
}

// MergeSort sorts vec[start..start+sz-1] using tmp[] as a temporary array. if depth >= limit then switch to serial.
void MergeSortParallel(std::vector<int64_t> & vec, int64_t start, std::vector<int64_t> & tmp, 
                        int64_t tmpStart, int64_t sz, int64_t depth, int64_t limit) {

    if (sz <= 1) return;

    // If reached cutoff for parallelism, or the piece is small, do serial.
    if (depth >= limit || sz <= MERGESIZE) {
        MergeSortSerial(vec, start, tmp, tmpStart, sz, depth, limit);
        return;
    }

    int64_t half = sz >> 1;

    int64_t A = start;
    int64_t tmpA = tmpStart;
    int64_t B = start + half;
    int64_t tmpB = tmpStart + half;

    // Recurse in parallel on the two halves
    cilk_spawn MergeSortParallel(vec, A, tmp, tmpA, half, depth + 1, limit);
    MergeSortParallel(vec, B, tmp, tmpB, sz - half, depth + 1, limit);
    cilk_sync;

    // Merge sorted parts into tmp (in parallel)
    ParallelMerge(vec, A, A + half - 1, B, start + sz - 1, tmp, tmpStart);

    // Copy merged result back into vec.
    memcpy(&vec[A], &tmp[tmpStart], sizeof(int64_t) * sz);
}

// Initialize vec[start..end) with random numbers in parallel chunks.
void InitParallel(std::vector<int64_t> & vec) {
    const int64_t n = vec.size();
    const int64_t nchunks = (n + INIT_CHUNK - 1) / INIT_CHUNK;

    cilk_for (int64_t c = 0; c < nchunks; ++c) {
        int64_t start = c * INIT_CHUNK;
        int64_t end   = std::min(n, start + INIT_CHUNK);

        struct drand48_data buffer;
        srand48_r(static_cast<unsigned>(time(NULL) ^ (start << 6)), &buffer);

        for (int64_t i = start; i < end; ++i) {
            long int v;
            lrand48_r(&buffer, &v);
            vec[i] = v;
        }
    }
}


int main(int argc, char* argv[]) {

    if (argc != 3) {
        std:: cout << "Usage " << argv[0] << " <vector size> <cutoff>\n";
        std:: cout << " cutoff is 0-based recursion depth at which parallelism stops (0 => serial)\n";
        return -1;
    }

    int64_t sz = atol(argv[1]);
    if (sz <= 0) {
        std::cerr << "enter valid vector size (> 0) !" << std::endl;
        return -1;
    }

    int64_t cutoff = atol(argv[2]);

    std::vector<int64_t> input(sz);
    std::vector<int64_t> tmp(sz);

    // Parallel initialization in chunks to allow drand48_r per-chunk.
    InitParallel(input);


    auto start = std::chrono::system_clock::now();
    MergeSortParallel(input, 0, tmp, 0, sz, /*depth=*/0 , /*limit=*/cutoff);
    auto end = std::chrono::system_clock::now();

    std::cout<<"\nruntime = "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()<<" ms\n";

    auto mistakes = ValidateParallel(input);
    std::cout << "Mistakes=" << mistakes << "\n";

    assert ( (mistakes == 0)  && " Validate() failed");

    return 0;
}