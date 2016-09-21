CXX=g++
CXXFLAGS=-O3 -Wall -g -mfpmath=sse -m32 -msse2

#merlinscatter_interface.a: merlinscatter_interface.o testclass.o
#	
merlinscatter_interface.o: #testclass.o
	$(CXX) $(CXXFLAGS) -c merlinscatter_interface.cpp

clean:
	@rm *.o *~
