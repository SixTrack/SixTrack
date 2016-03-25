# O---------------------------------------------------------------------o
# |
# | sixtrack makefile - management of options
# |
# |  TODO:
# |  hdf5
# |  g77, g95, lf95, nagfor
# |  logic for 64bit
# |
# |
# o---------------------------------------------------------------------o
# |

SHELL:=/bin/bash

# Library Path configuration
# Mac
#MAC_BOINC_LIB_PATH     := ../../../*.a
#MAC_LIBS_PATH          := /opt/local/lib/gcc48/i386/libstdc++.a /opt/local/libexec/llvm-3.3/lib/clang/3.3/lib/darwin/libclang_rt.eprintf.a /opt/local/lib/gcc48/gcc/x86_64-apple-darwin12/4.8.3/i386/libgcc_eh.a

#Windows
Windows_BOINC_LIB_PATH := C:\andrei_6track\win32\libxp2\*.lib
Windows_MS_SDK_PATH    := "C:\Program Files\Microsoft SDKs\Windows\v7.0A\Lib\User32.lib" "C:\Program Files\Microsoft SDKs\Windows\v7.0A\Lib\psapi.lib" "C:\Program Files\Microsoft SDKs\Windows\v7.0A\Lib\AdvAPI32.Lib"

# default options, can be overriden by the command line
#OPTIONS := fast tilt tracking crlibm ifort
OPTIONS := fast tilt tracking crlibm gfortran O3
#OPTIONS := crlibm bnlelens gfortran boinc api

NOPT  := $(strip $(subst O,,$(lastword $(filter O3 O2 O1 O0, $(OPTIONS)))))



ALL_OPTIONS := $(sort api beamgas big bignblz bnlelens boinc bpm cernlib collimat cpss cr crlibm da debug fast fio g77 g95 gfortran hdf5 ifort lf95 nagfor naglib small sse2 sse3 sse41 sse42 tilt tracking vvector windows O0 O1 O2 O3 O4)

# list of valid user options
USER_OPTIONS := $(sort api beamgas bignblz bnlelens boinc bpm cernlib collimat cpss cr crlibm da debug fast fio g77 g95 gfortran hdf5 ifort lf95 nagfor naglib sse2 sse3 sse41 sse42 tilt tracking windows O0 O1 O2 O3 O4)

ALL_OPTIONS := $(sort api beamgas big bignblz bnlelens boinc bpm cernlib collimat cpss cr crlibm da debug fast fio g77 g95 gfortran hdf5 ifort lf95 nagfor naglib small sse2 sse3 tilt tracking vvector windows)

# list of valid compiler options
COMPILER_OPTIONS := g77 g95 gfortran ifort lf95 nagfor


# options dependencies (+)enable, -disable, !exclude

beamgas    := +collimat
boinc      := +cpss
bpm        := !da !collimat
collimat   := -crlibm !da
cpss       := +cr !da !collimat
da         := -naglib !collimat
hdf5       := +collimat
windows    := +boinc +crlibm

# dictionary between options and patterns

# decks options to be checked in the .ast file
DECKS_OPT := collimat cr  bnlelens debug
# if disabled <lhs>  delete deck <rhs> 
collimat := nwrtcoll
bnlelens := nwrtbnl
cr       := checkpt
debug    := dumps

# if enabled <lhs>  delete deck <rhs> 
DECKS_OPT_INV := cernlib naglib
cernlib  := plotdumy
naglib   := nagdumy

# clean-up, expand and classify OPTIONS

expand      = $(foreach o,$1,$o $(if $($o),( $(call expand,$(patsubst +%,%,$($o))) ),))
PROVIDED   := $(sort $(patsubst +%,%,$(OPTIONS)))
#EXPANDED   := $(strip $(call expand,$(PROVIDED)))
EXPANDED   := $(PROVIDED)
ENABLED    := $(sort $(filter-out -% !% ( ),$(EXPANDED)))
ENABLED    := $(filter-out O0 O1 O2 O3,$(ENABLED))
DISABLED   := $(sort $(patsubst -%,%,$(filter -%,$(EXPANDED))))
EXCLUDED   := $(sort $(patsubst !%,%,$(filter !%,$(EXPANDED))))
UNSELECTED := $(filter-out $(ENABLED),$(USER_OPTIONS))
CONFLICT_D := $(filter $(DISABLED),$(ENABLED))
CONFLICT_E := $(filter $(EXCLUDED),$(ENABLED))


# Macro of enabling more options
is_equ = $(if $(patsubst $1,,$2),,t)
is_ena = $(filter $1,$(ENABLED))
if_any = $(if $(filter $1,$2),$3,$4)
if_all = $(if $(call is_equ,$(sort $1),$(filter $1,$(ENABLED))),$2,$3)

# sixtrack DESTDIR
space  =
space +=

# Build version of SixTrack
VERSION:= $(shell grep "data version" sixtrack.s | sed -e "s/.*'\([0-9]\).\([0-9]*\).\([0-9]*\).*/\1\2\3/")

DESTDIR := SixTrack_$(VERSION)_$(subst $(space),_,$(ENABLED))_O$(NOPT)
$(info Sixtrack build dir: $(DESTDIR))

# options to enable from other combinations
$(call enable,if_all,collimat da,,tracking)

# diagnostics
DEBUG:=yes
ifneq ($(CONFLICT_D)$(CONFLICT_E)$(DEBUG),no)
$(info Options provided: $(PROVIDED))
$(info Options expanded: $(EXPANDED))
$(info Options enabled:  $(ENABLED))
$(info Options disabled: $(DISABLED))
$(info Options excluded: $(EXCLUDED))
$(info Conflict with disabled: $(CONFLICT_D))
$(info Conflict with excluded: $(CONFLICT_E))
$(if $(CONFLICT_D)$(CONFLICT_E),$(error ERROR: please select a consistent set of options),)
endif

# compilers configuration
# currently supported ifort and gfortran
# detecting the selected compiler. If none selected, gfortran is default
COMPILER := $(call if_any,$(COMPILER_OPTIONS),$(ENABLED),$(filter $(COMPILER_OPTIONS),$(ENABLED)),gfortran)

# ifort configuration
FC_ifort   := ifort
FCF_ifort := -mia32 -fp-model source -fp-model strict -fp-model no-except
sseopt=$(filter sse2 sse3 sse41 see42,$(ENABLED))
ifneq ($(sseopt),)
  sse2:=sse2
  sse3:=sse3
  sse41:=sse4.1
  sse42:=sse4.2
  FCF_ifort:= $(FCF_ifort) -m$($(sseopt))
endif
FCF_ifort := $(FCF_ifort) -O$(NOPT)
FCLD_ifort := -m32
FCL_ifort  := -static -m32

# gfortran configuration
FC_gfortran   := gfortran -m32 
FCF_gfortran  := -frecord-marker=4 -g -fno-second-underscore -funroll-loops -mfpmath=sse -msse2
ifneq ($(NOPT),)
  FCF_gfortran:=$(FCF_gfortran) -O$(NOPT)
endif
FCLD_gfortran := 
FCL_gfortran  := -static

# nagfor configuration
NAGLIBC :=-L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/32 -lg2c
FC_nagfor  := nagfor $(NAGLIBC)
FCF_nagfor := -O$(NOPT) -abi=32 -dusty -maxcontin=60 -ieee=full -g90
FCL_nagfor := -Bstatic -abi=32
FCLD_nagfor:= -abi=32

# pgf90 configuration
FC_pgf90   :=pgf90
FCF_pgf90  :=-m32 -Kieee
ifneq ($(NOPT),)
  FCF_pgf90:= $(FCF_pgf90)-fast
endif
FCL_pgf90  :=-Bstatic -m32
FCLD_pgf90 :=-m32

# g95 configuration
#g95libc :=-L/usr/lib/ -lc -L/usr/lib/gcc/i386-redhat-linux/4.1.1 -lgcc_eh
#g95libc :=-L/usr/lib/ -lc -L/usr/lib/gcc/i586-redhat-linux/4.4.0 -lgcc_eh
#g95libc :=-L/usr/lib/ -lc -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/32 -lgcc_eh
FC_g95  :=g95 $(g95libc)
FCF_g95 :=-Wa,--32 -g -fno-second-underscore -funroll-loops
FCL_g95 :=-static -m32 -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/32 -lg2c
FCLD_g95:=-m32 -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/32 -lg2c

# lf95 configuration
FC_lf95  :=lf95
FCF_lf95 :=--o0 --tp -Wa,--32 #o>0 produce wrong results
FCL_lf95 :=-static -m32
FCLD_lf95:=-m32



# selecting the enabled compiler
FC   := $(FC_$(COMPILER))    #compiler
FCF  := $(FCF_$(COMPILER))   #complile options
FCLD := $(FCLD_$(COMPILER))  #dynamic link options
FCL  := $(FCL_$(COMPILER))   #static link options

# if ifort is not sourced displaying a messsage and exiting
ifeq ($(filter ifort,$(ENABLED)),ifort)
  ifeq ($(INTEL_LICENSE_FILE),)
    $(info *** Error: Intel compiler not set up. If on lxplus call:)
    $(info   source /afs/cern.ch/sw/IntelSoftware/linux/setup.sh)
    $(info   source /afs/cern.ch/sw/IntelSoftware/linux/x86_64/xe2011/composerxe-2011.2.137/bin/ifortvars.sh ia32)
    $(info   source /afs/cern.ch/sw/IntelSoftware/linux/x86_64/xe2011/composerxe-2011.2.137/bin/iccvars.sh ia32)
    $(error "Intel compiler not set up")
  endif
endif
# end compiler configuration

# add specific rules external library options
# these are default settings for Linux. Will be later overriden for Windows and Mac

# crlibm
CRLIBM_O := $(call if_any,crlibm,$(ENABLED),$(DESTDIR)/crlibm/crlibm.a,)

# cernlib
CERNLIB  := $(call if_any,cernlib,$(ENABLED),-L./lib32 -lmathlib -lgraflib -lgrafX11 -lpacklib_noshift -lkernlib_noshift -lg2c -lX11 -lxcb -lXau -lXdmcp -ldl -lpthread -lgcc_eh,)

# boinc
ifeq ($(filter boinc,$(ENABLED)),boinc)
  BOINCROOT= /afs/cern.ch/project/sixtrack/boinc_linux32
  BOINCAPI := $(BOINCROOT)
  BOINCLIB := $(BOINCROOT)
  BOINCZIP := $(BOINCROOT)/zip
  BOINC    :=  $(BOINCAPI)/libboinc_api.a $(BOINCLIB)/libboinc.a $(BOINCZIP)/libboinc_zip.a -ldl -lnsl -lstdc++ -lpthread -lm -lc  -lgcc_eh
  BOINC_O_FILES  :=  boinc_api_fortran.o boinc_zipitall.o
endif

# prepare AST rule files for ASTUCE and DAFOR
MASKS_DAFOR := sixvefox sixscfox sixda

#different masks depending on DA
MASKS_DA   := lielib dabnews sixscfox sixsc sixda
MASKS_OT   := track sixve sixvefox dabnews lielib beamgas
MASKS      := $(call if_any,da,$(ENABLED),$(MASKS_DA),$(MASKS_OT))

# preparing fortran files to be preprocessed with dafor. 
# variable is empty if to be preprocessed by dafor, else contains a comment
# added the d to the name to not conflict with options management
sixvefox_dafor  := 
sixscfox_dafor  := 
sixda_dafor     := 
lielib_dafor    := \#
dabnews_dafor   := \#
sixsc_dafor     := \#
track_dafor     := \#
sixve_dafor     := \#
beamgas_dafor   := \#


#detecting operating system
OPSYS := $(shell uname -s)

# Linux
# various options used in commands that vary on different OSs
/            := /
;            := ;
output       := -o
objects      := .o
obj_out      := -o
c            := -c
p            := -p
"            := "
DAFOR        := dafor
ASTUCE       := astuce
SIXTRACK     := sixtrack
WIN_CRLIBM_O :=

#linking options on Linux
SX_LINUX     := $(FCL) $(BOINC_O_FILES:%=$(DESTDIR)$/%) $(CERNLIB)

# Darwin
ifeq ($(OPSYS),Darwin)
  $(info Building on Mac)
  # Mac specific option
  FC +=     -mmacosx-version-min=10.6
  # removing -static keywoard from linking step
  SX_LINUX := $(filter-out -static,$(SX_LINUX))
  # path to BOINC library
  ifeq ($(filter boinc,$(ENABLED)),boinc)
    BOINC     := $(MAC_BOINC_LIB_PATH)
  endif
  #Mac specific libraries
  SX_DARWIN := $(MAC_LIBS_PATH)
endif

# Windows
# OPSYS variable not working on Windows. Running uname -s from cmd return Windows
# running the same command from make returns a different value
# better OS detection for Windows needed
ifeq ($(OS),Windows_NT)
  $(info Building on windows)
  / := $(strip \ )
  ; := &
  c := /c
  " := #"
  p := 
  SX_LINUX        := 
  FCF             := 
  BOINC_O_FILES   := 
  output          := /exe:
  objects         := .obj
  obj_out         := /Fo
  # SixTrack wiki shows compilation for ifort only
  FC              := ifort
  # path to BOINC Library
  BOINC           := $(Windows_BOINC_LIB_PATH)
  O_FILES         := $(patsubst %.o,%$(objects),$(O_FILES))
  WIN_ARCH        := /arch:$(call if_any,sse2 sse3,$(ENABLED),$(filter sse2 sse3,$(ENABLED)),IA32)
  WIN_IFORT       := /fp:source /fp:strict /fp:except-
  WIN_OPTIONS     := /O2 /MT /names:lowercase /assume:underscore
  ICL_WIN_OPTIONS := /I. /D_CRT_SECURE_NO_WARNINGS /Qstd=c99
  WIN_CRLIBM      := addition_scs.c atan.c atan_fast.c cosine.c csh_fast.c disable_xp.c division_scs.c double2scs.c dtoaf.c dtoa_c.c enable_xp.c ericc.c exp.c exp_fast.c log10.c log.c log_fast.c multiplication_scs.c rem_pio2.c round_near.c scs2double.c sine.c tan.c trigo_fast.c zero_scs.c
  # windows crlibm objects
  WIN_CRLIBM_O    := $(patsubst %.c,%$(objects),$(WIN_CRLIBM))
  # setting Linux crlibm objects to blank to remove any dependency in the rules
  CRLIBM_O        :=
  MS_SDK          := $(Windows_MS_SDK_PATH)
  # Windows linking options
  SX_WINDOWS      := /VERBOSE:LIB $(WIN_OPTIONS)  $(MS_SDK) /link
  ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
     $(info AMD64 achitecture)
  endif

  ifeq ($(PROCESSOR_ARCHITECTURE),x86)
     $(info x86 Architecture)
  endif
endif

# preparing proper file names with extension
MASKS:= $(addprefix $(DESTDIR)$/,$(MASKS))
O_FILES   := $(MASKS:=.o)
O_FILES   := $(O_FILES) plato_seq.o
AST_FILES := $(MASKS:=.ast)
FN_FILES   := $(MASKS:=n.f)
F_FILES   := $(MASKS:=.f)
# .s files to be copied
S_FILES := dabnew.s lielib.s sixtrack.s

# separate variable for source required on Windows
SRC_AST := .$/ast_mask$/$/

# separate variable for source required on Windows
SRC_CRLIBM := .$/crlibm$/$/

# SED command to remove inactive options and decks. Handling 4 situations when option to be removed is first, middle, last, or only
ASTUCE_INACTIVE := $(filter-out $(ENABLED),$(ALL_OPTIONS))

# adjustment from the make_six file. Doesn't fit our current logic.
# ifort doesn't compile with this options
ASTUCE_INACTIVE := $(filter-out small,$(ASTUCE_INACTIVE))
DECKS_INACTIVE  := $(filter-out $(ENABLED),$(DECKS_OPT))
DECKS_INACTIVE  := $(DECKS_INACTIVE) $(filter $(ENABLED),$(DECKS_OPT_INV))



# preparing expression to convert .ast files (remove unselected options and unused decks)
SED_EXP := $(foreach elem,$(ASTUCE_INACTIVE),s/\(^df .*\), *$(elem) *,\(.*\)/\1,\2/;s/\(^df .*\), *$(elem) *$$/\1/;s/\(^df  *$(elem) *,\)\(.*\)/df \2/;s/\(^df  *$(elem) *$$\)/df/; )
SED_DCK := $(foreach elem,$(DECKS_INACTIVE),/^e  *$($(elem))/d;)
SED_EXP := $(SED_EXP)$(SED_DCK)

#CRLIBM variables
# include compiler options and CRLIBM makefile
CRLIBM_SRC := .$/crlibm
CRLIBM_DST := $(DESTDIR)$/crlibm

CC = gcc
OBJS := addition_scs.o multiplication_scs.o zero_scs.o \
     division_scs.o double2scs.o scs2double.o zero_scs.o \
     trigo_fast.o tan.o sine.o cosine.o \
     exp.o logsix.o exp_fast.o log_fast.o \
     atan.o atan_fast.o division_scs.o \
     ericc.o csh_fast.o log10.o rem_pio2.o \
     dtoa_c.o dtoaf.o round_near.o disable_xp.o enable_xp.o

OBJ1:= csh_fast.o trigo_fast.o tan.o sine.o cosine.o exp.o exp_fast.o logsix.o log_fast.o log10.o atan.o atan_fast.o
OBJ2:= addition_scs.o double2scs.o multiplication_scs.o division_scs.o scs2double.o zero_scs.o
OBJS := $(OBJ1) $(OBJ2) ericc.o disable_xp.o enable_xp.o round_near.o dtoaf.o dtoa_c.o rem_pio2.o

OBJ1:= $(addprefix $(CRLIBM_DST)$/,$(OBJ1))
OBJ2:= $(addprefix $(CRLIBM_DST)$/,$(OBJ2))
OBJS:= $(addprefix $(CRLIBM_DST)$/,$(OBJS))

HEAD1:= $(addprefix $(CRLIBM_SRC)$/,crlibm.h crlibm_private.h scs.h scs_config.h scs_private.h crlibm_config.h)
HEAD2:= $(addprefix $(CRLIBM_SRC)$/,scs.h scs_config.h scs_private.h)

AM1_CFLAGS = -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith  -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline -pedantic -fno-strict-aliasing
AM2_CFLAGS = -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline
DEFS = -DLINUX_INLINE -DHAVE_CONFIG_H -I.
CFLAGS = -g -O2 -m32

COMPILE1 = $(CC) $(DEFS) $(AM1_CFLAGS) $(CFLAGS)
COMPILE2 = $(CC) $(DEFS) $(AM2_CFLAGS) $(CFLAGS)

#################
#
#     Rules
#
#################


.PHONY: clean sixtrack help crlibm

sixtrack: $(DESTDIR)$/sixtrack

crlibm: $(CRLIBM_DST)$/crlibm.a

$(DESTDIR)$/sixtrack: $(O_FILES) $(BOINC_O_FILES) $(CRLIBM_O) $(WIN_CRLIBM_O)
	$(info Building: $(DESTDIR)$/sixtrack)
	@$(FC)  $(output)$@ $(O_FILES) $(BOINC) $(CRLIBM_O) $(WIN_CRLIBM_O) $(SX_LINUX) $(SX_WINDOWS) $(SX_DARWIN)


# sed cand handle at most 192 words in an expression (experimentally determined) (expr for one option is 12 words long)
# splitting original expression into 2 sections and storing them to files
$(DESTDIR): $(S_FILES)
	$(info Preparing $@)
	@mkdir $(DESTDIR)
	@cp $(S_FILES) $(DESTDIR)
	@echo $"$(wordlist 1,192,$(SED_EXP))$" > $(DESTDIR)$/sed_expr_1
	@echo $"$(wordlist 193,384,$(SED_EXP))$" > $(DESTDIR)$/sed_expr_2

$(CRLIBM_DST): | $(DESTDIR)
	$(info Preparing $@)
	@mkdir  $(CRLIBM_DST)

$(CRLIBM_DST)$/crlibm.a: |$(CRLIBM_DST)

$(CRLIBM_DST)$/crlibm.a: $(CRLIBM_DST) $(OBJS)
	$(info Building $(@F))
	@ar -rv $(CRLIBM_DST)/crlibm.a $(OBJS) &>$@.log
	@if [[ `wc -l <$@.log` -gt 0 ]] $; then wc -l  $@.log $; fi

$(OBJ1): $(CRLIBM_DST)$/%.o: $(CRLIBM_SRC)$/%.c $(CRLIBM_SRC)/%.h $(HEAD1)
	$(info Compiling $(@F))
	@$(COMPILE1) -o $@ -c $< &>$@.log
	@if [[ `wc -l <$@.log` -gt 0 ]] $; then wc -l  $@.log $; fi

$(OBJ2): $(CRLIBM_DST)$/%.o: $(CRLIBM_SRC)$/%.c $(HEAD2)
	$(info Compiling $(@F))
	@$(COMPILE2) -o $@ -c $< &>$@.log
	@if [[ `wc -l <$@.log` -gt 0 ]] $; then wc -l  $@.log $; fi

$(CRLIBM_DST)$/rem_pio2.o: $(CRLIBM_SRC)$/rem_pio2.c $(CRLIBM_SRC)/coefpi2.h $(HEAD1)
	$(info Compiling $(@F))
	@$(COMPILE1) -o $@ -c $< &>$@.log
	@if [[ `wc -l <$@.log` -gt 0 ]] $; then wc -l  $@.log $; fi

$(CRLIBM_DST)/ericc.o: $(CRLIBM_SRC)$/ericc.c $(CRLIBM_SRC)$/crlibm.h $(CRLIBM_SRC)$/crlibm_private.h
	$(info Compiling $(@F))
	@$(COMPILE1) -o $@ -c $< &>$@.log
	@if [[ `wc -l <$@.log` -gt 0 ]] $; then wc -l  $@.log $; fi
$(CRLIBM_DST)/disable_xp.o: $(CRLIBM_SRC)$/disable_xp.c
	$(info Compiling $(@F))
	@$(COMPILE1) -o $@ -c $< &>$@.log
	@if [[ `wc -l <$@.log` -gt 0 ]] $; then wc -l  $@.log $; fi
$(CRLIBM_DST)/enable_xp.o: $(CRLIBM_SRC)$/enable_xp.c
	$(info Compiling $(@F))
	@$(COMPILE1) -o $@ -c $< &>$@.log
	@if [[ `wc -l <$@.log` -gt 0 ]] $; then wc -l  $@.log $; fi
$(CRLIBM_DST)/round_near.o: $(CRLIBM_SRC)$/round_near.c
	$(info Compiling $(@F))
	@$(CC) -m32 -std=c99 -W -Wall -pedantic -o $@ -c $< &>$@.log
	@if [[ `wc -l <$@.log` -gt 0 ]] $; then wc -l  $@.log $; fi
$(CRLIBM_DST)/dtoa_c.o: $(CRLIBM_SRC)$/dtoa_c.c
	$(info Compiling $(@F))
	@$(CC) -m32 -std=c99 -W -Wall -pedantic -o $@ -c $< &>$@.log
	@if [[ `wc -l <$@.log` -gt 0 ]] $; then wc -l  $@.log $; fi
$(CRLIBM_DST)/dtoaf.o: $(CRLIBM_SRC)$/dtoaf.c
	$(info Compiling $(@F))
	@$(CC) -m32 -std=c99 -W -Wall -pedantic -o $@ -c $< &>$@.log
	@if [[ `wc -l <$@.log` -gt 0 ]] $; then wc -l  $@.log $; fi

$(AST_FILES): $(DESTDIR)$/%.ast: $(SRC_AST)%.ast |$(DESTDIR)
	$(info Copy $(@F))
	@sed -f $(DESTDIR)$/sed_expr_1 $< > $(DESTDIR)$/temp_sed $; sed -f $(DESTDIR)$/sed_expr_2 $(DESTDIR)$/temp_sed > $@


$(FN_FILES): %n.f : %.ast $(DESTDIR)$/astuce $(S_FILES)
	$(info Applying astuce on $(@F))
	@cd $(DESTDIR) $; .$/astuce < $(<F) >$(*F).astuce_log

$(F_FILES): %.f : %n.f $(DESTDIR)$/dafor
	$(info Applying dafor on $(@F))
	@cp $< $(DESTDIR)/fort.1
	@cd $(DESTDIR) $; .$/dafor >$(*F).dafor_log
	@cp $(DESTDIR)/fort.2 $@

$(O_FILES): %$(objects): %.f
	$(info Building $(@F))
	@$(FC) $(FCF) $(WIN_ARCH) $(WIN_IFORT) $(WIN_OPTIONS) $(obj_out)$@ $(c) $< &>$@.log
	@if [[ `wc -l <$@.log` -gt 0 ]] $; then wc -l  $@.log $; fi

$(BOINC_O_FILES): %.o : %.cpp
	$(info Building $@)
	@gcc -o $@ -m32 -c -I$(BOINCROOT) -I$(BOINCAPI) -I$(BOINCLIB) -I$(BOINCZIP) $< &>$@.log
	@if [[ `wc -l <$@.log` -gt 0 ]] $; then wc -l  $@.log $; fi

$(DESTDIR)$/astuce : astuce.f
	$(info Building astuce)
	@$(FC) $(FCF) $(output)$(DESTDIR)$/astuce $<

$(DESTDIR)$/dafor : dafor.f
	$(info Building dafor)
	@$(FC) $(FCF) $(output)$(DESTDIR)$/dafor $<

$(WIN_CRLIBM_O): %.obj: $(SRC_CRLIBM)%.c
	$(info Building $(@F))
	@icl $(WIN_ARCH) $(WIN_IFORT) $(WIN_OPTIONS) $(ICL_WIN_OPTIONS) /Fo$(DESTDIR)$/crlibm$/$@ /c $<

help:
	@echo ***Usage example:
	@echo Default OPTIONS: $(OPTIONS)
	@echo Default COMPILER: gfortran
	@echo Usage: make or
	@echo        make OPTIONS=[list of options]
	@echo known options: $(ALL_OPTIONS)

clean:
	$(info Removing $(DESTDIR))
	@if [[ -d $(DESTDIR) ]] $;then rm -r $(DESTDIR) $; else echo "No $(DESTDIR) found"$; fi

# end of makefile
