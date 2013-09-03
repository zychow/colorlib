CXX = g++
CXXFLAGS = -ggdb3 -Wall -lgsl -lgslcblas

testprogram: main.o colorCalculation.o
	$(CXX) main.o colorCalculation.o -o testprogram $(CXXFLAGS)
main.o:main.cpp
	$(CXX) -c main.cpp -o main.o $(CXXFLAGS)
colorCalculation.o:colorCalculation.cpp
	$(CXX) -c colorCalculation.cpp -o colorCalculation.o $(CXXFLAGS)

clean:
	rm -f *.o
