# Define the compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -I/opt/gurobi/gurobi1100/linux64/include/
LDFLAGS = -L/opt/gurobi/gurobi1100/linux64/lib -lgurobi_c++ -lgurobi110

# Target name
#TARGET = deg_xy
TARGET = deg
# Source file
#SRC = deg_aradi_div_prop_product_xy.cpp
SRC = deg_aradi_div_prop.cpp

# Compile and link
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

# Clean up build files
clean:
	rm -f $(TARGET)

