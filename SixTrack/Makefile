SHELL = /bin/sh

%.o: %.f
	$(FC) -c $(FCF) $<

dafor: dafor.o
	$(FC) $(FCL) -o dafor dafor.o

astuce: astuce.o
	$(FC) $(FCL) -o astuce astuce.o $(NAGLIBC)

SixTrack: track.o sixve.o sixvefox.o dabnews.o lielib.o $(CRLIBM)
	$(FC) $(FCL) -o SixTrack track.o sixve.o sixvefox.o dabnews.o lielib.o \
	$(CERNLIB) $(CRLIBM)

SixTrack_da: sixda.o sixsc.o sixscfox.o dabnew.o lielib.o $(CRLIBM)
	$(FC) $(FCL) -o SixTrack_da sixda.o sixsc.o sixscfox.o dabnew.o lielib.o \
	$(NAGLIB) $(CRLIBM)

# make clean will execute this: 
clean: 
	rm -f track.f sixve.f sixvefox.f sixda.f sixsc.f sixscfox.f dabnew.f dabnews.f lielib.f
	rm -f SixTrack SixTrack_da dafor astuce
	rm -f crlibm/*.o crlibm/crlibm.a
	rm -f *.o *~
