# compiler
FC := $(F90)
#
# compile flags
# FCFLAGS = -g -c -kind=byte
ifeq ($(F90),nagfor)
FCFLAGS = -g -c -kind=byte
else
FCFLAGS = -g -c
endif
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
main.o: sparse_header_debug.o
#
%.o: %.F90 
	$(FC) $(FCFLAGS) -o $@ $<
#
%.mod: %.h
	$(FC) $(FCLAGS) -o $@ $<

clean:
	rm -f *.o *.mod
