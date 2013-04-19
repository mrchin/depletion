# compiler
FC := ifort 
#
# compile flags
# FCFLAGS = -g -c -kind=byte
FCFLAGS = -g -c
# link flags
FLFLAGS =
#
# source files and objects
SRCS = $(patsubst %.F90, %.o, $(wildcard *.F90))
#
# program name
PROGRAM = sparse_tester
#
all: $(PROGRAM)
#
$(PROGRAM): $(SRCS)
	$(FC) $(FLFLAGS) -o $@ $^
#
depletion.o: depletion_header.o sparse_header.o
sparse.o: sparse_header.o   
#
%.o: %.F90 
	$(FC) $(FCFLAGS) -o $@ $<
#
%.mod: %.h
	$(FC) $(FCLAGS) -o $@ $<

clean:
	rm -f *.o *.mod
