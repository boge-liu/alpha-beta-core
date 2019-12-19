CC=g++
CPPFLAGS=-I
LDFLAGS=-g
BOOST_ROOT=PATH_TO_BOOST_ROOT
DEPS = bigraph.h dynamic.h experiment.h gephi.h incremental_test.h mischievous.h paper.h preprocessing.h random.h static.h utility.h Multithread.h dyn_rebuild.h
OBJ = bigraph.o dynamic.o experiment.o gephi.o incremental_test.o main.o mischievous.o paper.o preprocessing.o random.o static.o utility.o Multithread.o dyn_rebuild.o

%.o: %.cpp $(DEPS)
	$(CC) -std=c++1y $(LDFLAGS) -c -O3 -o $@ $< $(CPPFLAGS)$(BOOST_ROOT)

abcore: $(OBJ)
	$(CC) -std=c++1y -pthread $(LDFLAGS) -O3 -o $@ $^ $(CPPFLAGS)$(BOOST_ROOT)
