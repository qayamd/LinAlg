# Compiler settings
CXX = g++
CXXFLAGS = -Wall -std=c++11

# Target executable name
TARGET = demo

# Source files
SOURCES = demo.cpp linalg.cpp

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Default target
all: $(TARGET)

# Rule to link the program
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to compile the source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

# Clean the build
clean:
	rm -f $(OBJECTS) $(TARGET)

# Dependencies
demo.o: demo.cpp linalg.h
linalg.o: linalg.cpp linalg.h
