CC	= gfortran#ifort #gfortran
CFLAGS	=  -Ofast #-g -fbounds-check -Wall -ffpe-trap=invalid

LDFLAGS	=
EXEC	= main

OBJS = main.o init.o solver.o wavefunctions.o fields.o

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)

init.o: init.f90

wavefunctions.o: wavefunctions.f90 init.o

fields.o: fields.f90 init.o wavefunctions.o

solver.o: solver.f90 init.o fields.o wavefunctions.o

main.o: main.f90 init.o wavefunctions.o fields.o solver.o

.SUFFIXES : .f90 .o

%.o : %.mod

.f90.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC) *.mod
