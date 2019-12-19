CC=g++
CPPFLAGS=-I
LDFLAGS=-g
BOOST_ROOT=/import/sigmod00/1/scratch/z5119798/multikcore/boost_1_68_0/
DEPS = bigraph.h dynamic.h experiment.h gephi.h incremental_test.h mischievous.h paper.h preprocessing.h random.h static.h utility.h Multithread.h dyn_rebuild.h
OBJ = bigraph.o dynamic.o experiment.o gephi.o incremental_test.o main.o mischievous.o paper.o preprocessing.o random.o static.o utility.o Multithread.o dyn_rebuild.o

%.o: %.cpp $(DEPS)
	$(CC) -std=c++1y $(LDFLAGS) -c -O3 -o $@ $< $(CPPFLAGS)$(BOOST_ROOT)

abcore: $(OBJ)
	$(CC) -std=c++1y -pthread $(LDFLAGS) -O3 -o $@ $^ $(CPPFLAGS)$(BOOST_ROOT)
