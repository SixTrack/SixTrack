SHELL = /bin/sh

CC = gcc
OBJS = addition_scs.o multiplication_scs.o zero_scs.o \
     division_scs.o double2scs.o scs2double.o zero_scs.o \
     trigo_fast.o tan.o sine.o cosine.o \
     exp.o logsix.o exp_fast.o log_fast.o \
     atan.o atan_fast.o division_scs.o \
     ericc.o csh_fast.o log10.o rem_pio2.o \
     dtoa_c.o dtoaf.o round_near.o disable_xp.o enable_xp.o
all: crlibm.a
crlibm.a: $(OBJS)
	ar -rv crlibm.a $(OBJS)
AM1_CFLAGS = -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith  -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline -pedantic -fno-strict-aliasing
DEFS = -DLINUX_INLINE -DHAVE_CONFIG_H -I.
CFLAGS = -g -O2 -m32
COMPILE1 = $(CC) $(DEFS) $(AM1_CFLAGS) $(CFLAGS)
csh_fast.o: csh_fast.c crlibm.h crlibm_private.h scs.h \
        scs_config.h scs_private.h crlibm_config.h \
        csh_fast.h
	$(COMPILE1) -c csh_fast.c
rem_pio2.o: rem_pio2.c crlibm.h crlibm_private.h scs.h \
        scs_config.h scs_private.h crlibm_config.h \
        coefpi2.h 
	$(COMPILE1) -c rem_pio2.c
trigo_fast.o: trigo_fast.c crlibm.h crlibm_private.h scs.h \
        scs_config.h scs_private.h crlibm_config.h \
        trigo_fast.h 
	$(COMPILE1) -c trigo_fast.c
tan.o: tan.c crlibm.h crlibm_private.h scs.h \
        scs_config.h scs_private.h crlibm_config.h \
        tan.h 
	$(COMPILE1) -c tan.c
sine.o: sine.c crlibm.h crlibm_private.h scs.h \
        scs_config.h scs_private.h crlibm_config.h \
        sine.h 
	$(COMPILE1) -c sine.c
cosine.o: cosine.c crlibm.h crlibm_private.h scs.h \
        scs_config.h scs_private.h crlibm_config.h \
        cosine.h 
	$(COMPILE1) -c cosine.c
exp.o: exp.c crlibm.h crlibm_private.h scs.h \
	scs_config.h scs_private.h crlibm_config.h \
	exp.h
	$(COMPILE1) -c exp.c
exp_fast.o: exp_fast.c crlibm.h crlibm_private.h scs.h \
	scs_config.h scs_private.h crlibm_config.h \
	exp_fast.h
	$(COMPILE1) -c exp_fast.c
logsix.o: logsix.c logsix.h crlibm.h crlibm_private.h scs.h \
	scs_config.h scs_private.h crlibm_config.h \
        logsix.h
	$(COMPILE1) -c logsix.c
log_fast.o: log_fast.c crlibm.h crlibm_private.h scs.h \
	scs_config.h scs_private.h crlibm_config.h \
        log_fast.h
	$(COMPILE1) -c log_fast.c
log10.o: log10.c log10.h crlibm.h crlibm_private.h scs.h \
	scs_config.h scs_private.h crlibm_config.h \
        log10.h
	$(COMPILE1) -c log10.c
atan.o: atan.c atan.h crlibm.h crlibm_private.h scs.h \
        scs_config.h scs_private.h crlibm_config.h \
        atan_fast.h
	$(COMPILE1) -c atan.c
atan_fast.o: atan_fast.c crlibm.h crlibm_private.h scs.h \
        scs_config.h scs_private.h crlibm_config.h \
        atan_fast.h
	$(COMPILE1) -c atan_fast.c
ericc.o: ericc.c crlibm.h crlibm_private.h 
	$(COMPILE1) -c ericc.c
disable_xp.o: disable_xp.c
	$(COMPILE1) -c disable_xp.c
enable_xp.o: enable_xp.c
	$(COMPILE1) -c enable_xp.c
round_near.o: round_near.c
	$(CC) -m32 -std=c99 -W -Wall -pedantic -c round_near.c
dtoa_c.o: dtoa_c.c
	$(CC) -m32 -std=c99 -W -Wall -pedantic -c dtoa_c.c
dtoaf.o: dtoaf.c
	$(CC) -m32 -std=c99 -W -Wall -pedantic -c dtoaf.c

AM2_CFLAGS = -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline 
COMPILE2 = $(CC) $(DEFS) $(AM2_CFLAGS) $(CFLAGS)

addition_scs.o: addition_scs.c scs.h scs_config.h scs_private.h
	$(COMPILE2) -c addition_scs.c
double2scs.o: double2scs.c scs.h scs_config.h scs_private.h
	$(COMPILE2) -c double2scs.c
multiplication_scs.o: multiplication_scs.c scs.h scs_config.h \
	scs_private.h
	$(COMPILE2) -c multiplication_scs.c
division_scs.o: division_scs.c scs.h scs_config.h \
	scs_private.h
	$(COMPILE2) -c division_scs.c
scs2double.o: scs2double.c scs.h scs_config.h scs_private.h
	$(COMPILE2) -c scs2double.c
zero_scs.o: zero_scs.c scs.h scs_config.h scs_private.h
	$(COMPILE2) -c zero_scs.c

clean:
	rm -f $(OBJS) crlibm.a
