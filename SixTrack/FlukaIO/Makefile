#Set this to @ to keep the makefile quiet
SILENCE = @

#---- Outputs ----#
COMPONENT_NAME = FlukaIO
TARGET_LIB = \
	lib/lib$(COMPONENT_NAME).a
TARGET_SHARED_LIB = \
	lib/lib$(COMPONENT_NAME).so

TEST_TARGET = \
	$(COMPONENT_NAME)_tests

ifeq ($(BUILD64), Y)
	TARGET_SUFFIX = 64
	CPPUTEST_CPPFLAGS += -m64
	G77FLAGS += -m64
endif
ifeq ($(BUILD32), Y)
	TARGET_SUFFIX = 32
	CPPUTEST_CPPFLAGS += -m32
	CPPUTEST_LDFLAGS += -m32 -O2
	G77FLAGS += -m32
endif

COMPONENT_NAME = FlukaIO
TARGET_LIB = \
	lib/lib$(COMPONENT_NAME)$(TARGET_SUFFIX).a
TARGET_SHARED_LIB = \
	lib/lib$(COMPONENT_NAME)$(TARGET_SUFFIX).so

TEST_TARGET = \
	$(COMPONENT_NAME)$(TARGET_SUFFIX)_tests

#--- Inputs ----#
PROJECT_HOME_DIR = .
#CPPUTEST_HOME = /home/dsinuela/src/cpputest
CPP_PLATFORM = Gcc

# in case of Ubuntu 16.04, comment out -Werror
WARNINGFLAGS = -Wall -Wswitch-default -Wswitch-enum -Werror 
#CPPUTEST_CFLAGS are set to override malloc and free to get memory leak detection in C programs
#CPPUTEST_CFLAGS = -Dmalloc=cpputest_malloc -Dfree=cpputest_free
#GCOVFLAGS = -fprofile-arcs -ftest-coverage

CPPUTEST_CPPFLAGS += $(WARNINGFLAGS)
CPPUTEST_CPPFLAGS += -g -O2
CPPUTEST_CPPFLAGS += -fpic
#CPPUTEST_CPPFLAGS += -DFLUKAIO_NODELAY

G77 = gfortran
G77FLAGS += $(WARNINGFLAGS)
G77FLAGS += -g -O2

#SRC_DIRS is a list of source directories that make up the target library
#If test files are in these directories, their IMPORT_TEST_GROUPs need
#to be included in main to force them to be linked in.  By convention
#put them into an AllTests.h file in each directory
SRC_DIRS = \
	src

#TEST_SRC_DIRS is a list of directories including 
# - A test main (AllTests.cpp by conventin)
# - OBJ files in these directories are included in the TEST_TARGET
# - Consequently - AllTests.h containing the IMPORT_TEST_GROUPS is not needed
# - 
TEST_SRC_DIRS = \
	tests \
	tests/fakes

#includes for all compiles
INCLUDES =\
	-I.\
	-Iinclude
#	-I$(CPPUTEST_HOME)/include\

#Flags to pass to ld
USER_LIBS =
LD_LIBRARIES += -lstdc++

include ComponentMakefile

all: $(TARGET_LIB) $(TARGET_SHARED_LIB)

libs: $(TARGET_LIB) $(TARGET_SHARED_LIB)

$(TARGET_SHARED_LIB): $(OBJ)
	$(SILENCE)echo Building archive $@
	$(SILENCE)mkdir -p lib
	$(SILENCE)$(CC) $(LDFLAGS) -shared -o $@ $^

STUFF_TO_CLEAN += $(TARGET_SHARED_LIB)

# Samples

.PHONY: clean_samples

clean_samples:
	rm -rf samples/{ServerTest,ClientTest,fserver,fclient}

samples: ServerTest ClientTest fserver fclient

ServerTest: samples/ServerTest.c $(TARGET_LIB)
	$(CC) $^ $(CPPFLAGS) -I./include -o samples/$@

ClientTest: samples/ClientTest.c $(TARGET_LIB)
	$(CC) $^ $(CPPFLAGS) -I./include -o samples/$@

fserver: samples/fserver.f $(TARGET_LIB)
	$(G77) $^ $(G77FLAGS) -o samples/$@

fclient: samples/fclient.f $(TARGET_LIB)
	$(G77) $^ $(G77FLAGS) -o samples/$@

STUFF_TO_CLEAN += samples/ServerTest
STUFF_TO_CLEAN += samples/ClientTest
STUFF_TO_CLEAN += samples/fserver
STUFF_TO_CLEAN += samples/fclient

# Documentation

.PHONY: doc

doc:
	@echo Making documentation
	$(SILENCE)doxygen doc/Doxyfile 1> doxygen.log
	#@$(RM) doxygen.log
	cd doc; make

STUFF_TO_CLEAN += doc/html doxygen.log

show:
	@echo "CPPFLAGS=$(CPPFLAGS)"
	@echo "LDFLAGS=$(LDFLAGS)"

.PHONY: clean_all

clean_all:
	$(SILENCE)echo "Making clean_all"
	$(SILENCE)make clean BUILD32=Y > /dev/null
	$(SILENCE)make clean BUILD64=Y > /dev/null

