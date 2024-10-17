# hpc 
#CC = gcc -DH5_USE_16_API
#CFLAGS =  -I$(CPATH) -std=gnu99 -O3 -march=native -mtune=native -flto -fopenmp -funroll-loops -g
#LDFLAGS = -L$(LD_LIBRARY_PATH) -lm -lpthread -lz -ldl -lgsl -lgslcblas -lhdf5 -lhdf5_hl -lmpi

# pc
CC = gcc -DH5_USE_16_API
CFLAGS =  -I/usr/include/gsl/ -I/usr/include/hdf5/serial/ -std=gnu99 -O3 -march=native -mtune=native -flto -fopenmp -funroll-loops -g
LDFLAGS = -L/usr/lib/x86_64-linux-gnu -lm -lpthread -lz -ldl \
	  -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_hl -lmpi -lgsl -lgslcblas

## choose model
MODEL=iharm3d
#MODEL=analytic

#ipolarray.o
OBJIPO = \
main.o geodesics.o image.o ipolarray.o radiation.o tetrads.o geometry.o \
model_tetrads.o model_radiation.o model_geodesics.o \
model_$(MODEL).o

HEADERS = decs.h params.h defs.h constants.h

#rule to make them all
all: $(OBJIPO) ipole

# rule to make object files
main.o: main.c $(HEADERS)
	$(CC) $(CFLAGS) -c main.c

geodesics.o: geodesics.c 
	$(CC) $(CFLAGS) -c geodesics.c

geometry.o: geometry.c 
	$(CC) $(CFLAGS) -c geometry.c

image.o: image.c 
	$(CC) $(CFLAGS) -c image.c

ipolarray.o: ipolarray.c
	$(CC) $(CFLAGS) -c ipolarray.c

model_geodesics.o: model_geodesics.c 
	$(CC) $(CFLAGS) -c model_geodesics.c

model_radiation.o: model_radiation.c 
	$(CC) $(CFLAGS) -c model_radiation.c

model_tetrads.o: model_tetrads.c
	$(CC) $(CFLAGS) -c model_tetrads.c

radiation.o: radiation.c 
	$(CC) $(CFLAGS) -c radiation.c

tetrads.o: tetrads.c 
	$(CC) $(CFLAGS) -c tetrads.c

model_$(MODEL).o: model/model_$(MODEL).c
	$(CC) $(CFLAGS) -c model/model_$(MODEL).c


#rule to make executable from object files, here no header files necessery?
ipole: $(OBJIPO)
	$(CC) $(CFLAGS) -o ipole $(OBJIPO) $(LDFLAGS)


# rules to clean up
clean:
	rm *.o model/*.o
cleanup:
	rm output_ipole/ipole*.ppm output_ipole/ipole.dat



