CC=g++ -std=c++11 -fopenmp -O3

all: rz_pic

rz_pic: main.o PoissonSolver.o Matrix.o Pusher.o Particles.o Grid.o Interpolation.o
	$(CC) -o rz_pic main.o PoissonSolver.o Matrix.o Pusher.o Particles.o Grid.o Interpolation.o
	rm -rf *.o

main.o: main.cpp
	$(CC) -c main.cpp

PoissonSolver.o: Field/PoissonSolver.cpp
	$(CC) -c Field/PoissonSolver.cpp

Matrix.o: Tools/Matrix.cpp
	$(CC) -c Tools/Matrix.cpp

Pusher.o: Particles/Pusher.cpp
	$(CC) -c Particles/Pusher.cpp

Particles.o: Particles/Particles.cpp
	$(CC) -c Particles/Particles.cpp

Grid.o: Grid/Grid.cpp
	$(CC) -c Grid/Grid.cpp

Interpolation.o: Particles/Interpolation.cpp
	$(CC) -c Particles/Interpolation.cpp

clean:
	rm -rf *.o rz_pic
