# Name of the executable
EXE = TestGenWF

# Source files
SRC = TestGenWF.C

# Compiler
CXX = g++

# Compiler & linker flags from ROOT
ROOTCFLAGS   = $(shell root-config --cflags)
ROOTLIBS     = $(shell root-config --libs)

# Optimization / warnings
CXXFLAGS = -O2 -Wall -fPIC $(ROOTCFLAGS) -std=c++17
LDFLAGS  = $(ROOTLIBS)

all: $(EXE)

$(EXE): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(EXE) *.o
