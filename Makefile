CC = g++
OPT = -O3
#OPT = -g
WARN = -Wall
CFLAGS = $(OPT) $(WARN) $(INC) $(LIB)

# Source file
SIM_SRC = Cache_simulator.cpp      

# Object file (named according to the source file)
SIM_OBJ = Cache_simulator.o

#################################

# Default rule to create sim_cache
all: sim_cache
	@echo "my work is done here..."

# Rule for making sim_cache
sim_cache: $(SIM_OBJ)
	$(CC) -o sim_cache $(CFLAGS) $(SIM_OBJ) -lm
	@echo "-----------DONE WITH SIM_CACHE-----------"

# Generic rule for compiling .cpp files into .o files
.cpp.o:
	$(CC) $(CFLAGS) -c $<

# Rule to remove all .o files plus the sim_cache binary
clean:
	rm -f *.o sim_cache

