CC = g++
OBJS = main.o amoeba.o callfitness.o callminimise.o estimates.o FFT.o fitGAL.o fitPSF.o megio.o psfminimise.o powell.o profile.o quasi.o

default: all

all: $(OBJS)
	$(CC) $(OBJS) -I/usr/include/cfitsio -lcfitsio -lfftw3 -L/usr/lib -llevmar -lblas -llapack -lm fftw++.cc /usr/lib/libf2c.a -o meg
	
main: main.cpp fitGAL.h fitPSF.h globals.h
	$(CC) main.cpp

fitGAL: fitGAL.cpp fitGAL.h estimates.h megio.h callfitness.h callminimise.h profile.h globals.h
	$(CC) fitGAL.cpp
	
fitPSF: fitPSF.cpp fitPSF.h megio.h callfitness.h psfminimise.h profile.h estimates.h
	$(CC) fitPSF.cpp
	
callfitness: callfitness.h callfitness.cpp estimates.h profile.h globals.h FFT.h
	$(CC) callfitness.cpp
	
estimates: estimates.h estimates.cpp globals.h
	$(CC) estimates.cpp
    
megio: megio.h megio.cpp globals.h fitsio.h
	$(CC) megio.cpp -I/usr/local/include
	
FFT: FFT.h FFT.cpp megio.h fftw++.cc globals.h
	$(CC) FFT.cpp -I/usr/include -lfftw3 -lm
    
callminimise: callminimise.h callminimise.cpp amoeba.cpp powell.cpp quasi.cpp callfitness.h dlib/all/source.cpp lmbc_core.c levmar.h
	$(CC) callminimise.cpp -I/home/melissa/Documents/Y4/Project/dlib-17 -L/usr/lib -llevmar -lm /usr/lib/libf2c.a -lblas -llapack
	
psfminimise: psfminimise.h psfminimise.cpp amoeba.cpp powell.cpp quasi.cpp callfitness.h
	$(CC) psfminimise.cpp
    
profile: profile.h profile.cpp
	$(CC) profile.cpp
    
amoeba: amoeba.cpp matrix.h minfuncs.h
	$(CC) amoeba.cpp
    
powell: powell.cpp matrix.h minfuncs.h
	$(CC) powell.cpp
    
quasi: quasi.cpp matrix.h minfuncs.h callfitness.h
	$(CC) quasi.cpp

clean:
	\rm *.o meg

