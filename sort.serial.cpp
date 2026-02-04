/*
 Merge sort.
 Author: Milind Chabbi
 Date: Jan/8/2022
 */

#include<iostream>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <cstring>
#include <chrono>
#ifdef CILK
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#else
#define cilk_spawn
#define cilk_sync
#define cilk_for for
#endif

// Tune this value suitably (2 is definately not good for performance)
#define MERGESIZE (2)
using namespace std;

/* Validate returns a count of the number of times the i-th number is less than (i-1)th number. */
int64_t Validate(const vector<int64_t> & input) {
    int64_t count = 0;
    for (int64_t i = 1; i < input.size(); i++) {
        if (input[i] < input[i-1]) {
            count++;
        }
    }
    return count;
}

/* SerialMerge merges two sorted arrays vec[Astart..Aend] and vec[Bstart..Bend] into vector tmp starting at index tmpIdx. */
void SerialMerge(vector<int64_t> & vec,  int64_t Astart,  int64_t Aend, int64_t Bstart, int64_t Bend, vector<int64_t> & tmp, int64_t tmpIdx){
    if (Astart <= Aend && Bstart <= Bend) {
        for (;;) {
            if (vec[Astart] < vec[Bstart]) {
                tmp[tmpIdx] = vec[Astart];
                tmpIdx++;
                Astart++;
                if (Astart > Aend)
                    break;
            } else {
                tmp[tmpIdx] = vec[Bstart];
                Bstart++;
                tmpIdx++;
                if (Bstart > Bend)
                    break;
            }
        }
    }
    if (Astart > Aend) {
        memcpy(&tmp[tmpIdx], &vec[Bstart], sizeof(int64_t) * (Bend - Bstart + 1));
    } else {
        memcpy(&tmp[tmpIdx], &vec[Astart], sizeof(int64_t) * (Aend - Astart + 1));
    }
}


// NOTE: reimplement the below function into a parallel MergeSort.
// MergeSort sorts vec[start..start+sz-1] using tmp[] as a temporary array.
void MergeSort(vector<int64_t> & vec, int64_t start, vector<int64_t> & tmp, int64_t tmpStart, int64_t sz, int64_t depth, int64_t limit){
    
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
    MergeSort(vec, A, tmp, tmpA, half, depth, limit);
    // Sort right half.
    MergeSort(vec, B, tmp, tmpB, start + sz - B, depth, limit);
    // More for paralle case?
    
    // Merge sorted parts into tmp.
    SerialMerge(vec, A, A + half-1, B, start + sz -1, tmp, tmpStart);
    // Copy result back to vec.
    memcpy(&vec[A], &tmp[tmpStart],sizeof(int64_t) * (sz));
}

// Initialize vec[start..end) with random numbers.
void Init(vector<int64_t> & vec, int64_t start, int64_t end){
    struct drand48_data buffer;
    srand48_r(time(NULL), &buffer);
    for(int64_t i = start; i < end; i++){
        lrand48_r(&buffer, &vec[i]);
    }
}

int main(int argc, char**argv) {
    if (argc != 3) {
        cout << "Usage " << argv[0] << " <vector size> <cutoff>\n";
        return -1;
    }
    int64_t sz = atol(argv[1]);
    cout<<"\nvector size=" << sz << "\n";
    if (sz == 0) {
        cout << "Usage " << argv[0] << " <vector size> <cutoff>\n";
        return -1;
    }
    int64_t cutoff = atol(argv[2]);
    vector<int64_t> input(sz);
    vector<int64_t> tmp(sz);
    
    const int64_t SZ = 100000;
    for(int64_t i = 0; i < input.size(); i+=SZ){
        Init(input, i, min(static_cast<int64_t>(input.size()), static_cast<int64_t>(i+SZ)));
    }
    auto now = chrono::system_clock::now();
    MergeSort(input, 0,  tmp, 0,  input.size(),  /*depth=*/0 , /*limit=*/cutoff);
    cout<<"\nmillisec="<<std::chrono::duration_cast<std::chrono::milliseconds>(chrono::system_clock::now() - now).count()<<"\n";
    auto mistakes = Validate(input);
    cout << "Mistakes=" << mistakes << "\n";
    assert ( (mistakes == 0)  && " Validate() failed");
    return 0;
}

