# Edmanuel Torres
# eetorers@gmail.com
# https://github.com/eetorres/mdargon

CXX = g++

PROG = argon_gl

SCRS = $(PROG).cxx
OBJS = $(SCRS:.cxx=.o)

MDIR = $(shell sh -c pwd)

OPTIM = -O2 -Wall

CFLAGS = $(OPTIM) -c  -I.
LFLAGS = $(OPTIM) -lm -lGL -lGLU -lglut

.SUFFIXES:    .cxx .o .h

.cxx.o:
	$(CXX) $(CFLAGS) $< -o $@

all:    clean $(OBJS)
	echo "==== Compiling  ====";
	$(CXX) -o $(PROG) $(OBJS) $(LFLAGS)

clean:
	echo "==== Cleaning ====";
	rm -rf *.o *.xyz *.dat $(OBJS) $(PROG)

