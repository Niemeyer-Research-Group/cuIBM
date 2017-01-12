.PHONY: lib
.SUFFIXES: .cpp

CC = /usr/local/cuda-7.5/bin/nvcc ${OSXOPTS}
OPTS = -arch=compute_20 -O3 -I../include

SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)

LIBNAME=libyaml-cpp.a

lib: $(OBJ)
	$(AR) ruv  ../$(LIBNAME) $?
	ranlib ../$(LIBNAME)

clean:
	@rm -f *.o

.cpp.o:
	$(CC) $(OPTS) -c $< -o $@


