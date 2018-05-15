SHELL = /bin/sh
.SUFFIXES: .o .f .F .c .f90 .F90


CC = g++
CFLAGS  = -ggdb

OBJS    =  main.o
test.exe: 	$(OBJS)	
	$(CC) -o main.exe $(CFLAGS) $(OBJS) 
.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.exe *.o

