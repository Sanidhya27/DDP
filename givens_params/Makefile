all: mimo_ofdm_precoding

mimo_ofdm_precoding: mimo_ofdm_precoding.o givens_decomposition.o
	$(CXX) -o $@ $^ -litpp

mimo_ofdm_precoding.o: mimo_ofdm_precoding.cpp givens_decomposition.hpp
	$(CXX) -c -o $@ $<

givens_decomposition.o: givens_decomposition.cpp givens_decomposition.hpp
	$(CXX) -c -o $@ $<

.PHONY: clean
clean:
	$(RM) *.o mimo_ofdm_precoding
