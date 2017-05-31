CC = gcc
CFLAGS = -Wall -g 
LIB = -lgmp
OBJ = RDRgmp.o

all: wNAF.exe

wNAF.exe: wNAF.c $(OBJ)
	$(CC) $(CFLAGS) wNAF.c $(OBJ) -o wNAF.exe $(LIB)

RDRgmp.o: RDRgmp.c

clean:
	rm *.o
	rm *~

