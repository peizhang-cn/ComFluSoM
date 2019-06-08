CC = g++

CFLAGS = -O3 -Wall -std=c++11

LFLAGS = -lhdf5_serial -lhdf5_cpp -fopenmp

INCLUDES = -I /usr/include/hdf5/serial/ -I $(ComFluSoM)/Library/DEM -I $(ComFluSoM)/Library/LBM -I $(ComFluSoM)/Library/DELBM -I /usr/include/eigen3/

TARGET = t_delbm001

all: $(TARGET)

$(TARGET) : $(TARGET).cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(TARGET).cpp $(LFLAGS)

clean:
	$(RM) $(TARGET) *.h5 *.xmf *.res
