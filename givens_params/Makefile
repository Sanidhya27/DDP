all: mimo_ofdm_precoding givens_decomposition_tests mimo_ofdm_precoding_tdl
CFLAGS ?= -Wall
CXXFLAGS ?= -Wall
CXXFLAGS=-std=c++11 
givens_decomposition_tests: givens_decomposition_tests.o givens_decomposition.o
	$(CXX) -pthread -o $@ $^ -L . -litpp -lgtest -lgtest_main

mimo-tdl-channel.o: mimo-tdl-channel.hpp mimo-tdl-channel.cpp

mimo_ofdm_precoding: mimo_ofdm_precoding.o givens_decomposition.o
	$(CXX) -o $@ $^ -litpp

mimo_ofdm_precoding_tdl: mimo_ofdm_precoding_tdl.o givens_decomposition.o mimo-tdl-channel.o
	$(CXX) -o $@ $^ -litpp

mimo_ofdm_precoding_tdl.o: mimo_ofdm_precoding_tdl.cpp givens_decomposition.hpp mimo-tdl-channel.hpp
	$(CXX) -std=c++11 -c -o $@ $<

mimo_ofdm_precoding.o: mimo_ofdm_precoding.cpp givens_decomposition.hpp
	$(CXX) -std=c++11 -c -o $@ $<

givens_decomposition_tests.o: givens_decomposition_tests.cpp givens_decomposition.hpp
	$(CXX) -c -o $@ $<

givens_decomposition.o: givens_decomposition.cpp givens_decomposition.hpp
	$(CXX) -c -o $@ $<

.PHONY: clean test

test: givens_decomposition_tests
	./givens_decomposition_tests

clean:
	$(RM) *.o mimo_ofdm_precoding givens_decomposition_tests
