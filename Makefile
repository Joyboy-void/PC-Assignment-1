CXX=clang++

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CXXFLAGS=-O0 -g 
else
	CXXFLAGS=-O3
endif

CILK_FLAGS=-fopencilk -DCILK
DR_FLAGS=-fsanitize=cilk
BM_FLAGS=-fcilktool=cilkscale-benchmark
SCALE_FLAGS=-fcilktool=cilkscale

OUT=bin

SRC_1=src/MergeSort.cpp
EXECUTABLE_1 =  sortSerial sortCilk sortSan sortBenchmark sortScale

SRC_2=src/NQueens.cpp
EXECUTABLE_2 =  queenSerial queenCilk queenSan queenScale queenBenchmark

all: $(EXECUTABLE_1) $(EXECUTABLE_2)
	@echo "==========================================="
	@echo "Targets             : Variables Used"
	@echo "-------------------------------------------"
	@echo "runSortSerial       : SIZE"
	@echo "runSortCilk         : SIZE CUTOFF CILK_NWORKERS"
	@echo "runSortSan          : SIZE CUTOFF"
	@echo "runSortScale        : SIZE CUTOFF"
	@echo "runQueenSerial      : N REPORTAFTER"
	@echo "runQueenCilk        : N CUTOFF REPORTAFTER CILK_NWORKERS"
	@echo "runQueenSan         : N CUTOFF REPORTAFTER"
	@echo "==========================================="

sortSerial: $(SRC_1)
	@echo "strating build..."
	@echo 
	$(CXX) $(CXXFLAGS) $< -o $(OUT)/$@

# as cutoff = 0 implies serial exicution 
runSortSerial: sortSerial
	./$< $(SIZE) 0 

sortCilk: $(SRC_1)
	$(CXX) $(CXXFLAGS) $(CILK_FLAGS) $< -o $(OUT)/$@ 

runSortCilk: sortCilk
	CILK_NWORKERS=$(CILK_NWORKERS) ./$< $(SIZE) $(CUTOFF)

sortSan: $(SRC_1)
	$(CXX) $(CXXFLAGS) $(CILK_FLAGS) $(DR_FLAGS) $^ -o $(OUT)/$@

runSortSan: sortSan
	CILK_NWORKERS=1 ./$< $(SIZE) $(CUTOFF)

sortBenchmark: $(SRC_1)
	$(CXX) $(CXXFLAGS) $(BM_FLAGS) $(CILK_FLAGS) $^ -o $(OUT)/$@

sortScale: $(SRC_1)
	$(CXX) $(CXXFLAGS) $(SCALE_FLAGS) $(CILK_FLAGS) $^ -o $(OUT)/$@

runSortScale: sortScale sortBenchmark
	python3 /home/apps/opencilk/cilktools/Cilkscale_vis/cilkscale.py -c ./sortScale -b ./sortBenchmark -cpus 1,2,3,4,8,16,32  --output-csv report.csv --output-plot report.pdf -a $(SIZE) $(CUTOFF)

queenSerial: $(SRC_2)
	$(CXX) $(CXXFLAGS) $(CILK_FLAGS) $< -o $(OUT)/$@

#setting cutoff = 0 runs the NQueens code serially.
runQueenSerial:	queenSerial
	./$< $(N) 0 $(REPORTAFTER)

queenCilk: $(SRC_2)
	$(CXX) $(CXXFLAGS) $(CILK_FLAGS) $< -o $(OUT)/$@

runQueenCilk: queenCilk
	CILK_NWORKERS=$(CILK_NWORKERS) ./$< $(N) $(CUTOFF) $(REPORTAFTER)

queenSan: $(SRC_2)
	$(CXX) $(CXXFLAGS) $(CILK_FLAGS) $(DR_FLAGS) $^ -o $(OUT)/$@

runQueenSan: queenSan
	CILK_NWORKERS=1 ./$< $(N) $(CUTOFF) $(REPORTAFTER)

queenBenchmark: $(SRC_2)
	$(CXX) $(CXXFLAGS) $(BM_FLAGS) $(CILK_FLAGS) $^ -o $(OUT)/$@

queenScale: $(SRC_2)
	$(CXX) $(CXXFLAGS) $(SCALE_FLAGS) $(CILK_FLAGS) $^ -o $(OUT)/$@

#TODO
runQueenScale: queenScale queenBenchmark
	
clean:
	rm -f bin/* && touch bin/.gitkeep