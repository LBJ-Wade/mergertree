
EXEC   = mergertree

OBJS   = allvars.o  cosmo.o  funcs_durham.o  init.o  io.o  main.o  misc.o  timestep.o  test.o


INCL   = allvars.h  cosmo.h  proto.h  Makefile


#OPT   += -DUPDATETYPETWO  #  This updates the positions of type 2 galaxies when the galaxies are written to file


CC       =    gcc          # sets the C-compiler (default)
OPTIMIZE =   -O3 -Wall    # optimization and warning flags (default)

GSL_LIBS = -L/opt/local/lib

LIBS   =   -g -lm  $(GSL_LIBS) -lgsl -lgslcblas 

CFLAGS =   -g $(OPTIONS) $(OPT) $(OPTIMIZE) $(GSL_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS)

tidy:
	rm -f $(OBJS) .$(EXEC)
