OPTIM = -O0 -fno-loop-optimize -fno-aggressive-loop-optimizations
DEBUG = -g
WARN_FLAGS = -Wall -W -Wwrite-strings -Winline -Wstrict-prototypes -Wnested-externs -Wpointer-arith -Wcast-align -Wcast-qual -Wshadow -Werror=vla
DEFINES = -DTESTING -DCA_N_U=20 -DCA_N_V=6 -DTESTING_LSQ_PACKAGE

CC = gcc
CC_FLAGS = -fstack-usage -fwrapv -fPIC ${OPTIM} ${DEBUG} ${WARN_FLAGS} # never optimise main
INCLUDE_FLAGS = -I../package/common
LINK_FLAGS = -lm

MAKE = make
RM = rm -f

.DEFAULT_GOAL = all

all : main.c package
	$(CC) ${CC_FLAGS} -o main.o main.c ../package/bin/solveActiveSet.so $(DEFINES) $(INCLUDE_FLAGS) $(LINK_FLAGS)

package:
	@cd ../package && make

clean :
	$(RM) *.o
	${RM} *.su
	@cd ../package && make clean

