#include <iostream>
#include <vector>
#include <set>
#include <atomic>
#include <cstdint>
#include <algorithm>
#include <cstdlib>
#include <chrono>

#include <cilk/cilk.h>
#include <cilk/cilk_api.h> 
#include <cilk/opadd_reducer.h>


// reducer to store all solutions,
struct solution_store {
    std::vector<std::vector<int>> sols;
    solution_store() = default;
};

// identity to initialize a reducer in-place 
void identity_all_sol(void *view){
    new (view) solution_store();
}

// merge: append right.sols into left.sols 
void merge_all_sol(void *left, void *right){

    solution_store* L = reinterpret_cast<solution_store*>(left);
    solution_store* R = reinterpret_cast<solution_store*>(right);
    
    // Append R->sols into L->sols
    L->sols.insert(L->sols.end(), R->sols.begin(), R->sols.end());
}

// reducer to store unique solutions
struct unique_store{
    std::set<std::vector<int>> uniq;
    unique_store() = default;
};

void identity_uniq_sol(void* view){
    new (view) unique_store();
}

void merge_uniq_sol(void* left,void* right){
    unique_store* L = reinterpret_cast<unique_store*>(left);
    unique_store* R = reinterpret_cast<unique_store*>(right);

    // deterministically merge left and right sols
    L->uniq.insert(R->uniq.begin(), R->uniq.end()); 
}


// reducer variable to store all solutions
solution_store cilk_reducer(identity_all_sol, merge_all_sol) all_solutions;

// reducer variable to store uniq solutions
unique_store cilk_reducer(identity_uniq_sol,merge_uniq_sol) uniq_solutions;

// Atomic counter to support 'stopafter' 
std::atomic<long long> found_count(0);

/* 
    try to increment found_count only if stopafter permits.

   -> returns true if the caller may record the solution.
*/
bool try_to_record(long long stopafter){
    if (stopafter == -1) {
        found_count.fetch_add(1, std::memory_order_relaxed);
        return true;
    }

    //lock-free retry logic (since locking and unlocking is overhead)
    long long cur = found_count.load(std::memory_order_relaxed);
    while (cur < stopafter) {
        // implementing Compare and swap Atomic operation to update found_count. 
        if (found_count.compare_exchange_weak(cur, cur + 1, std::memory_order_relaxed)){
            return true;
        }
    }
    return false;
}

/* helpers for 8 symmetries */

// rotate 90 degrees clockwise 
// each queen is at (i,v[i]).
// rotation rule : roatate N x N by 90deg => (row,col) --> (col,N-1-row)

std::vector<int> rotate90(const std::vector<int>& v){
    int N = (int)v.size();
    std::vector<int> r(N);

    // for each queen placed
    for (int row = 0; row < N; ++row){
        int col = v[row];
       
        r[col] = N - 1 - row;
    }
    return r;
}

// reflects the board horizontally
// reflection rule for N x N : (row,col) --> (row,N-1-col)

std::vector<int> reflect(const std::vector<int>& v){
    int N = (int)v.size();
    std::vector<int> r(N);

    // for each queen
    for(int i = 0; i < N; ++i) 
        r[i] = N - 1 - v[i];

    return r;
}

// given one solution, generate all symmetric equivalents and pick one representative.

std::vector<int> canonical(const std::vector<int>& sol){

    std::vector<std::vector<int>> cand;
    std::vector<int> cur = sol;

    for (int k = 0; k < 4; k++) {
        cand.push_back(cur);
        cand.push_back(reflect(cur));
        cur = rotate90(cur);
    }

    // smallest as representative
    auto it = std::min_element(cand.begin(), cand.end());
    
    return *it;
}

// helper for n-queens using bit masks 
int count_trailing_zeros(uint64_t bit) {
    int count = 0;
    while ((bit & 1ULL) == 0) {
        bit >>= 1;
        count++;
    }
    return count;
} 

/* 
   n-queens recursion using bitmasks.
   board is passed by value so each child gets its own copy.
*/
void solve(int N, int row,
           uint64_t cols, uint64_t d1, uint64_t d2, /*bitmasks for columns,main diagonals,anti diagonals*/
           std::vector<int> board, int cutoff, long long stopafter){

    // stop if we've already reached stopafter
    if (stopafter != -1 && found_count.load(std::memory_order_relaxed) >= stopafter)
        return;

    if (row == N){
        // found a solution; only record if try_to_record allows

        if (try_to_record(stopafter)){
            // this is safe each worker has his private reducer view
            all_solutions.sols.push_back(board);

            // also add uniq sol to unique_store reducer
            std::vector<int> canon = canonical(board);
            uniq_solutions.uniq.insert(canon);
        }
        return;
    }

    // compute available columns as bits 0..N-1

    // if col 0 is available then avail = 0001 if N is 4

    uint64_t all_ones = (1ULL << N) - 1ULL;

    uint64_t avail = (~(cols | d1 | d2)) & all_ones;

    while (avail) {

        uint64_t bit = avail & -avail; // rightmost set bit
        avail -= bit; // remove the rightmost set bit for next ittr

        int col = count_trailing_zeros(bit); // count the number of zeros on right of rightmost set bit
                                        // ---> col index

        board[row] = col;

        // spawn if allowed by cutoff
        if (cutoff == -1 || row < cutoff) {
            cilk_spawn solve(N, row + 1, cols | bit, (d1 | bit) << 1, (d2 | bit) >> 1, board, cutoff, stopafter);
        } else {
            // serial
            solve(N, row + 1, cols | bit, (d1 | bit) << 1, (d2 | bit) >> 1, board, cutoff, stopafter);
        }
    }
    cilk_sync;
}

int main(int argc, char* argv[]) {

    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <N> <cutoff> <stopafter>\n";
        return 1;
    }

    int N = atoi(argv[1]);
    int cutoff = atoi(argv[2]);
    long long stopafter = atoll(argv[3]);

    if (N <= 0) {
        std::cerr << "N must be positive\n";
        return 1;
    }

    // initialize
    found_count.store(0);
    std::vector<int> board(N, -1);

    // start 
    auto start = std::chrono::system_clock::now();
    solve(N, 0, 0ULL, 0ULL, 0ULL, board, cutoff, stopafter);
    auto end = std::chrono::system_clock::now();

    // At this point, reduction has happened and all_solutions.sols contains all collected solutions
    std::vector<std::vector<int>> allsols = all_solutions.sols;

    // Print All solutions
    std::cout << "All solutions:" << std::endl;
    for (size_t i = 0; i < allsols.size(); ++i) {
        std::cout << "solution " << i << ":";
        for (int x : allsols[i]) 
            std::cout << " " << x;
        std::cout << std::endl;
    }

    int idx = 0;
    for(const auto& sol : uniq_solutions.uniq){
        std::cout << "unique solution " << idx++ << ":";
        for(int x : sol)
            std::cout << " " << x;
        std::cout << std::endl;
    }

    std::cout << "Number of solutions: " << allsols.size() << "\n";
    std::cout << "Number of unique solutions: " << uniq_solutions.uniq.size() << "\n";
    std::cout<<"\nruntime = "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()<<" ms\n";

    
    return 0;
}
