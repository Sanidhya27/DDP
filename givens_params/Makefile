all: mimo_ofdm_precoding givens_decomposition_tests

givens_decomposition_tests: givens_decomposition_tests.o givens_decomposition.o
	cp /usr/lib/x86_64-linux-gnu/libgtest.a .
	cp /usr/lib/x86_64-linux-gnu/libgtest_main.a .
	$(CXX) -pthread -o $@ $^ -L . -litpp -lgtest -lgtest_main

mimo_ofdm_precoding: mimo_ofdm_precoding.o givens_decomposition.o
	$(CXX) -o $@ $^ -litpp

mimo_ofdm_precoding.o: mimo_ofdm_precoding.cpp givens_decomposition.hpp
	$(CXX) -c -o $@ $<

givens_decomposition_tests.o: givens_decomposition_tests.cpp givens_decomposition.hpp
	$(CXX) -c -o $@ $<

givens_decomposition.o: givens_decomposition.cpp givens_decomposition.hpp
	$(CXX) -c -o $@ $<

.PHONY: clean test

test: givens_decomposition_tests
	./givens_decomposition_tests

clean:
	$(RM) *.o mimo_ofdm_precoding givens_decomposition_tests
