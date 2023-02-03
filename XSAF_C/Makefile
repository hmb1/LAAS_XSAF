ifneq ($(MAKECMDGOALS),clean)
 include $(MAKECMDGOALS).make
endif

OBJECTS = $(SOURCES:.cpp=.o)

$(MAKECMDGOALS): clean $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -f *.o $(EXECUTABLE)
