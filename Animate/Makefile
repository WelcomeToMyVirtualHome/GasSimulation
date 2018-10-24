CC=g++
CFLAGS=-Wall
LIBS= -lGL -lGLU -lglut
all:    animateAr
animateAr:    main.o  
	$(CC) -o animateAr main.o $(CFLAGS) $(LIBS)
main.o: main.cpp
	$(CC) -o main.o -c main.cpp $(CFLAGS)
clean:
	rm -f *.o
