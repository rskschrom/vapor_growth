# compile options
FC = gfortran
FCFLAGS = -Ofast
#FLIBS = -lnetcdff -lnetcdf
FLIBS = 
LPATH = -L/usr/local/lib
IPATH = -I/usr/local/include
OBJS = thermo.o micro.o
MAIN = test.f90

# main program compilation
test : $(MAIN) $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $^ $(IPATH) $(FLIBS)

# object files
%.o : %.f90
	$(FC) $(FCFLAGS) -c $< $(IPATH) $(FLIBS)

# clean
clean :
	rm *.o *.mod
