#.PHONY : clean

FC = gfortran
FCFLAGS = -O3 -W -Wall -fcheck=all

SRCS = cube_neighbor.f90 cube_partition_sfc.f90 cube_partition_band.f90
OBJS = $(SRCS:.f90=.o)
MODS = $(SRCS:.f90=.mod)
LIBNAME = libshared


.SUFFIXES:
.SUFFIXES:.f90 .o

build : $(OBJS)
	$(FC) -shared $^ -o $(LIBNAME).so

clean :
	rm -f $(OBJS) $(MODS) $(LIBNAME).so

.f90.o :
	$(FC) $(FCFLAGS) -c -fPIC $<
