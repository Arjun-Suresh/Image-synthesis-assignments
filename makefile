# Please see the header information in the supplied .ccp file for more
# information on this template project. For more information on C programming
# in Linux, please see the excellent introduction to makefile structure and
# the gcc compiler here:
#
# http://www.cs.txstate.edu/labs/tutorials/tut_docs/Linux_Prog_Environment.pdf

CC		= g++
LDFLAGS 	= -lglut -lGL -lm
PROJECT		= pr04
FILES		= pr04.cpp
PROJECT2	= subSurface
FILES2		= subSurface.cpp

subSurface: subSurface.o pr04.o
	${CC} ${CFLAGS} -o subSurface.o pr04.o

subSurface.o:subSurface.cpp
	${CC} ${CFLAGS} -c subSurface.cpp ${LDFLAGS}

pr04.o:pr04.cpp
	${CC} ${CFLAGS} -c pr04.cpp ${LDFLAGS}
	

