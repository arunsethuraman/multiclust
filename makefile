# type
# 'make' to compile mixture model as run_mix
# 'make admix' to compile mixture model as run_admix

CC=gcc
CFLAGS=-std=c99 -Wall -W -O3 -ffast-math -fexpensive-optimizations -funroll-all-loops -finline-functions -pedantic
#CFLAGS=-std=c99 -Wall -W -pedantic -pg -O3
#CFLAGS=-std=c99 -Wall -W -pedantic -g
LDFLAGS=-lm #-fleading-underscore -I/home/tuf29140/lapack/lapack-3.4.2/lapacke/include/ -L/home/tuf29140/lapack/lapack-3.4.2 #-L/usr/lib64/atlas/ -I/usr/include -llapack -L/home/tuf29140/lapack/CBLAS/lib/#-I/usr/local/include -L/usr/local/lib -lgsl -l gslcblas

ifdef OLDWAY
	CFLAGS := $(CFLAGS) -DOLDWAY
endif

ifdef LAPACK
	CFLAGS := $(CFLAGS) -DLAPACK
	LDFLAGS := $(LDFLAGS) -llapack
endif


# Local variables
SRC = $(wildcard *.c)
HDR = $(wildcard *.h)
OBJ = $(SRC:.c=.o)
DEP = $(SRC:.c=.d) 
EXE = multiclust

multiclust : $(OBJ)
	$(CC) $(CFLAGS) -o multiclust $(OBJ) $(LDFLAGS)

test : $(OBJ)
	$(CC) $(CFLAGS) -o multiclust.test $(OBJ) $(LDFLAGS)

include $(DEP)

%.d : %.c
	-@$(SHELL) -ec '$(CC) -MM $(CFLAGS) $(IFLAGS) $< \
		| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@'

.PHONY : clean cleanall

clean:
	@rm -f *.o *.d 2> /dev/null

cleanall:
	@rm -f *.o *.d $(EXE) 2> /dev/null
