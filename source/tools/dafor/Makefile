# |
# o---------------------------------------------------------------------o
# |
# | astuce makefile
# |
# o---------------------------------------------------------------------o
# |

###################
# Project settings

PROJECT := dafor

#################
# Build settings
#

# architecture bit: detect/32/64 (default is detect)
ARCH    := detect

# debugging mode: yes/no (default is no)
DEBUG   := no

# profiling mode: yes/no (default is no)
PROFILE := no

# link with static libs: yes/no (default is no)
STATIC  := no

#############################
# Compilers/Linkers settings

# C compiler (default is ifort)
FC := ifort

# Linker (default is fortran compiler)
LD  = $(FC)

####################
# Includes settings

FILE_F90 := Makefile_f90

####################
# Makefile includes

makedir := ../../make
include $(makedir)/make.inc

# end of makefile
