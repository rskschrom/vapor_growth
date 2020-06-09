# compile options
FC = gfortran
FCFLAGS = -Ofast
#FLIBS = -lnetcdff -lnetcdf
FLIBS = 
LPATH = -L/usr/local/lib
IPATH = -I/usr/local/include
OBJS = thermo.o micro.o solver.o util.o
MAIN = box_model.f90

# main program compilation
box_model : $(MAIN) $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $^ $(IPATH) $(FLIBS)

# object files
%.o : %.f90
	$(FC) $(FCFLAGS) -c $< $(IPATH) $(FLIBS)

# clean
clean :
	rm *.o *.mod
