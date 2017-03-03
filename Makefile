## Set which compiler to use by defining CCCOM:
##GNU GCC compiler
#CCCOM=g++ -m64 -std=c++11 
##Clang compiler (good to use on Mac OS)
#CCCOM=clang++ -std=c++1y
##Intel C++ compiler (good to use with Intel MKL if available)
CCCOM=mpic++ -std=c++11 -g
#########


## Flags to give the compiler for "release mode"



#LIBFLAGS = -larmadillo
LIBSPECTRA = -I/media/xuejian/WORK/spectra/spectra-0.2.0/spectra-0.2.0/include/ -I/media/xuejian/WORK/spectra/eigen-eigen-07105f7124f9/







obj=main.o Parameter.o JC_Parameter.o OP.o Hamiltanian.o Evolution.o#Sub.o QWave.o Super.o DMRGP.o Corr.o
tevolve:$(obj)
	$(CCCOM) -o tevolve $(obj)  $(LIBSPECTRA)
main.o:main.cpp Mat.h CalQ.h
	$(CCCOM) -c main.cpp $(LIBSPECTRA)
Parameter.o:Parameter.cpp Parameter.h
	$(CCCOM) -c Parameter.cpp $(LIBSPECTRA)
JC_Parameter.o:JC_Parameter.cpp JC_Parameter.h Parameter.h
	$(CCCOM) -c JC_Parameter.cpp $(LIBSPECTRA)
OP.o:OP.cpp OP.h
	$(CCCOM) -c OP.cpp $(LIBSPECTRA)
Hamiltanian.o:Hamiltanian.cpp Hamiltanian.h
	$(CCCOM) -c Hamiltanian.cpp $(LIBSPECTRA)
Evolution.o:Evolution.cpp Evolution.h
	$(CCCOM) -c Evolution.cpp $(LIBSPECTRA)
.PHONY:clean
clean:
	rm -f tevolve $(obj)















