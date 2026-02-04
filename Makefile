CXX=/home/apps/opencilk/build/bin/clang++

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CXXFLAGS=-O0 -g -std=c++11 
else
	CXXFLAGS=-O3 -g -std=c++11
endif

CILK_FLAGS=-fopencilk -DCILK
DR_FLAGS=-fsanitize=cilk
BM_FLAGS=-fcilktool=cilkscale-benchmark
SCALE_FLAGS=-fcilktool=cilkscale

EXECUTABLE =  sortSerial sortCilk sortSan sortBenchmark sortScale
#SRC=sort.serial.cpp
SRC=sort.parallel2.cpp 
all: $(EXECUTABLE)

sortSerial: $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^

runSortSerial: sortSerial
	./$< $(SIZE) $(CUTOFF)

sortCilk: $(SRC)
	$(CXX) $(CXXFLAGS) $(CILK_FLAGS) -o $@ $^

runSortCilk: sortCilk
	CILK_NWORKERS=$(CILK_NWORKERS) ./$< $(SIZE) $(CUTOFF)

sortSan: $(SRC)
	$(CXX) $(CXXFLAGS) $(CILK_FLAGS) $(DR_FLAGS) -o $@ $^

runSortSan: sortSan
	CILK_NWORKERS=1 ./$< $(SIZE) $(CUTOFF)

sortBenchmark: $(SRC)
	$(CXX) $(CXXFLAGS) $(BM_FLAGS) $(CILK_FLAGS) -o $@ $^

sortScale: $(SRC)
	$(CXX) $(CXXFLAGS) $(SCALE_FLAGS) $(CILK_FLAGS) -o $@ $^

runSortScale: sortScale sortBenchmark
	python3 /home/apps/opencilk/cilktools/Cilkscale_vis/cilkscale.py -c ./sortScale -b ./sortBenchmark -cpus 1,2,3,4,8,16,32  --output-csv report.csv --output-plot report.pdf -a $(SIZE) $(CUTOFF)


clean:
	rm -f *.o *~ core $(EXECUTABLE)
