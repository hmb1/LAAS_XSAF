CXX = g++ -std=c++11 -fopenmp -march=native
CXXFLAGS = -Wall -c -I/usr/local/include
LDFLAGS = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -lSystem
SOURCES = main_gnu.cpp XSAF_C_gnu.cpp
EXECUTABLE = main_gnu.exe
