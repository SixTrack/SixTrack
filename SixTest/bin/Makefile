SHELL = /bin/sh

FC = gfortran
#FC = g77
FCF = -c -O -m32 -frecord-marker=4 
FCL = -static -m32 
#FCF = -c -O -m32
#FCL = -static -m32 

all: verify10 checkf10 compf10 checkf110 checkf1014

verify10: verify10.o
	 ${FC} -o verify10 ${FCL} verify10.o 

verify10.o: verify10.f
	${FC} -c ${FCF} verify10.f

checkf10: checkf10.o
	 ${FC} -o checkf10 ${FCL} checkf10.o 

checkf10.o: checkf10.f
	${FC} -c ${FCF} checkf10.f

checkf110: checkf110.o
	 ${FC} -o checkf110 ${FCL} checkf110.o 

checkf110.o: checkf110.f
	${FC} -c ${FCF} checkf110.f

checkf1014: checkf1014.o
	 ${FC} -o checkf1014 ${FCL} checkf1014.o 

checkf1014.o: checkf1014.f
	${FC} -c ${FCF} checkf1014.f

compf10: compf10.o
	 ${FC} -o compf10 ${FCL} compf10.o 

compf10.o: compf10.f
	${FC} -c ${FCF} compf10.f

clean:
	rm -f *.o verify10 checkf10 compf10 checkf110 checkf1014
