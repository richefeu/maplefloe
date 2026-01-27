
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  CXX = g++-15
  CXXFLAGS = -O3 -Wall -Wextra -pedantic -std=c++17
  GLLINK = `pkg-config --libs gl glu glut`
	GLFLAGS = `pkg-config --cflags gl glu glut`
	# on apple, use brew to install freeglut and mesa-glu

else
  CXX = g++
  CXXFLAGS = -O3 -Wall -std=c++17
  GLLINK = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXext -lX11
endif

# Ensure the toofus directory exists
ifeq ($(wildcard ./toofus),)
$(info Cloning ToOfUs repository)
$(shell git clone https://github.com/richefeu/toofus.git > /dev/null 2>&1)
endif

# The list of source files for the library SPICE
SOURCES = MappleFloe.cpp FloeElement.cpp Interaction.cpp MFloeAssembler.cpp

# Each cpp file listed below corresponds to an object file
OBJECTS = $(SOURCES:%.cpp=%.o)

.PHONY: all clean clone_toofus

all: pac run see

clean:
	@echo "\033[0;32m-> Remove object files\033[0m"
	rm -f *.o
	@echo "\033[0;32m-> Remove compiled applications\033[0m"
	rm -f pac run see
	@echo "\033[0;32m-> Remove libMappleFloe.a\033[0m"
	rm -f libMappleFloe.a

clean+: clean
	@echo "\033[0;32m-> Remove local folder toofus\033[0m"
	rm -rf ./toofus

clone_toofus:
	@if [ ! -d "./toofus" ]; then \
		@echo "\033[0;32m-> CLONING ToOfUs (Tools often used)\033[0m"; \
		git clone https://github.com/richefeu/toofus.git; \
	fi

%.o: %.cpp
	@echo "\033[0;32m-> COMPILING OBJECT" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o $@

libMappleFloe.a: $(OBJECTS)
	@echo "\033[0;32m-> BUILDING LIBRARY" $@ "\033[0m"
	ar rcs $@ $^

run: run.cpp libMappleFloe.a
	@echo "\033[0;32m-> BUILDING RUN APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -I toofus/toofus-gate -c $< -o run.o
	$(CXX) -o $@ run.o libMappleFloe.a

see: see.cpp libMappleFloe.a
	@echo "\033[0;32m-> BUILDING SEE APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o see.o $(GLFLAGS)
	$(CXX) -o $@ see.o libMappleFloe.a $(GLLINK)

pac: pac.cpp libMappleFloe.a
	@echo "\033[0;32m-> BUILDING PAC APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -I toofus/toofus-gate -c $< -o pac.o $(GLFLAGS)
	$(CXX) -o $@ pac.o libMappleFloe.a $(GLLINK)
