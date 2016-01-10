SHELL = /bin/sh

FC = gfortran
FCF =  -frecord-marker=4 -O4 -m32 -g -fno-second-underscore -funroll-loops -mfpmath=sse -msse2
FCL = -static -m32 -L/usr/lib/gcc/i686-linux-gnu/5.2.1
FCLD = -m32 -L/usr/lib/gcc/i686-linux-gnu/5.2.1
FCLMAP = -Wl,--print-map -static -m32 -L/usr/lib/gcc/i686-linux-gnu/5.2.1
CRLIBM = 
CERNLIB = 
gfortran_lib = 
BOINC = 
NAGLIBC = 
BEAMGAS = 
SIXTRACK = 

%.o: %.f
	$(FC) -c $(FCF) $<

%.o: %.f90
	$(FC) -c $(FCF) $<

all: 

dafor: dafor.o
	$(FC) $(FCLD) -o dafor dafor.o

astuce: astuce.o
	$(FC) $(FCLD) -o astuce astuce.o $(NAGLIBC)

myboinc.o: myboinc.f
	$(FC) -c $(FCF) myboinc.f

: $(MYHDF5LIB) track.o sixve.o sixvefox.o dabnews.o lielib.o $(BEAMGAS) $(BOINC) $(CRLIBM)
	$(FC) $(FCLMAP) -o  track.o sixve.o sixvefox.o dabnews.o lielib.o $(BEAMGAS) \
	$(BOINC) $(CERNLIB) $(CRLIBM) $(gfortran_lib) > sixtrack.map 2>&1

: sixda.o sixsc.o sixscfox.o dabnew.o lielib.o $(CRLIBM)
	$(FC) $(FCLMAP) -o  sixda.o sixsc.o sixscfox.o dabnew.o lielib.o \
	$(NAGLIB) $(NAGLIBC) $(CRLIBM) $(gfortran_lib) > sixtrack_da.map 2>&1

clean: 
	rm -f *.ast *.s *.c *.o  
