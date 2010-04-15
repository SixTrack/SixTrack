SHELL = /bin/sh

%.o: %.f
	$(FC) -c $(FCF) $<

dafor: dafor.o
	$(FC) -m32 -o dafor dafor.o

astuce: astuce.o
	$(FC) -m32 -o astuce astuce.o $(NAGLIBC)

%%SIXTRACK%%: track.o sixve.o sixvefox.o dabnews.o lielib.o $(BOINC) $(CRLIBM)
	$(FC) $(FCLMAP) -o %%SIXTRACK%% track.o sixve.o sixvefox.o dabnews.o lielib.o \
	$(BOINC) $(CERNLIB) $(CRLIBM) $(gfortran_lib) > sixtrack.map 2>&1

%%SIXTRACK_DA%%: sixda.o sixsc.o sixscfox.o dabnew.o lielib.o $(CRLIBM)
	$(FC) $(FCLMAP) -o %%SIXTRACK_DA%% sixda.o sixsc.o sixscfox.o dabnew.o lielib.o \
	$(NAGLIB) $(CRLIBM) $(gfortran_lib) > sixtrack_da.map 2>&1

clean: 
	rm -f *.ast *.s *.f *.c *.o %%SIXTRACK%% %%SIXTRACK_DA%%
