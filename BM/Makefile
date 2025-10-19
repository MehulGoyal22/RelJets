# Fortran compiler
FC = gfortran
FFLAGS = -O2 -Wall

# Target executable
TARGET = bm_solver

# Source files
SRC_MODULE = bm76_module.F90
SRC_MAIN   = test_bm76.F90

# Object files
OBJ_MODULE = bm76_module.o
OBJ_MAIN   = test_bm76.o

# Build rule
all: $(TARGET)

$(TARGET): $(OBJ_MODULE) $(OBJ_MAIN)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJ_MODULE) $(OBJ_MAIN)

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f *.o *.mod $(TARGET)
