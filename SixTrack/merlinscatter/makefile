CXX=g++
CXXFLAGS=-O3 -Wall -Wextra -g -mfpmath=sse -m32 -msse2 -std=c++11

merlinscatter.a: ElasticScatter.o DiffractiveScatter.o HelperFunctions.o merlinscatter_interface.o RandomNG.o RNG.o ACG.o Landau.o MLCG.o Normal.o Uniform.o LandauEnergyLoss.o
		ar rcs merlinscatter.a RandomNG.o ElasticScatter.o DiffractiveScatter.o HelperFunctions.o merlinscatter_interface.o RNG.o ACG.o Landau.o MLCG.o Normal.o Uniform.o LandauEnergyLoss.o

merlinscatter_interface.o:
	$(CXX) $(CXXFLAGS) -c merlinscatter_interface.cpp

ElasticScatter.o:
	$(CXX) $(CXXFLAGS) -c ElasticScatter.cpp

DiffractiveScatter.o:
	$(CXX) $(CXXFLAGS) -c DiffractiveScatter.cpp

HelperFunctions.o:
	$(CXX) $(CXXFLAGS) -c HelperFunctions.cpp

LandauEnergyLoss.o:
	$(CXX) $(CXXFLAGS) -c LandauEnergyLoss.cpp

#RNG files
RandomNG.o:
	$(CXX) $(CXXFLAGS) -c RandomNG/RandomNG.cpp
RNG.o:
	$(CXX) $(CXXFLAGS) -c RandomNG/RNG.cpp
ACG.o:
	$(CXX) $(CXXFLAGS) -c RandomNG/ACG.cpp
Landau.o:
	$(CXX) $(CXXFLAGS) -c RandomNG/Landau.cpp
MLCG.o:
	$(CXX) $(CXXFLAGS) -c RandomNG/MLCG.cpp
Normal.o:
	$(CXX) $(CXXFLAGS) -c RandomNG/Normal.cpp
Uniform.o:
	$(CXX) $(CXXFLAGS) -c RandomNG/Uniform.cpp
clean:
	@rm -f *.o *~ *.a
