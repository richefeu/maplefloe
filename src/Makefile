# ============================================================================
#
# ============================================================================

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)

  CXX = g++-15
  CXXFLAGS = -O3 -Wall -Wextra -pedantic -std=c++17 -D DEBUG
	
  GLLINK = `pkg-config --libs gl glu glut`
	GLFLAGS = `pkg-config --cflags gl glu glut`
	# on apple, use brew to install freeglut and mesa-glu
	
	GLFWLINK = `pkg-config --libs glfw3`
	GLFWFLAGS = `pkg-config --cflags glfw3`
	
else

  CXX = g++
  CXXFLAGS = -O3 -Wall -std=c++17
  GLLINK = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXext -lX11
	
	# pour le moment on a pas encore tester de compiler sur un PC linux
	
endif

# Ensure the toofus directory exists
ifeq ($(wildcard ./toofus),)
$(info Cloning ToOfUs repository)
$(shell git clone https://github.com/richefeu/toofus.git > /dev/null 2>&1)
endif

# The list of source files for the library
SOURCES = MapleFloe.cpp FloeElement.cpp Interaction.cpp Driving.cpp MFloeAssembler.cpp

# Each cpp file listed below corresponds to an object file
OBJECTS = $(SOURCES:%.cpp=%.o)

.PHONY: all clean clean+ clone_toofus

all: pac run see

clean:
	@echo "\033[0;32m-> Remove object files\033[0m"
	rm -f *.o
	@echo "\033[0;32m-> Remove compiled applications\033[0m"
	rm -f pac run seeold see
	@echo "\033[0;32m-> Remove libMapleFloe.a\033[0m"
	rm -f libMapleFloe.a

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

libMapleFloe.a: $(OBJECTS)
	@echo "\033[0;32m-> BUILDING LIBRARY" $@ "\033[0m"
	ar rcs $@ $^

run: run.cpp libMapleFloe.a
	@echo "\033[0;32m-> BUILDING RUN APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o run.o
	$(CXX) -o $@ run.o libMapleFloe.a

seeold: seeold.cpp libMapleFloe.a
	@echo "\033[0;32m-> BUILDING SEEOLD APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o seeold.o $(GLFLAGS)
	$(CXX) -o $@ seeold.o libMapleFloe.a $(GLLINK)
	
see: see.cpp libMapleFloe.a
	@echo "\033[0;32m-> BUILDING SEE APPLICATION (GLFW version)" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o see.o $(GLFWFLAGS)
	$(CXX) -o $@ see.o libMapleFloe.a $(GLFWLINK) -framework OpenGL

pac: pac.cpp MFloeAssembler.cpp
	@echo "\033[0;32m-> BUILDING PAC APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c MFloeAssembler.cpp -o MFloeAssembler.o
	$(CXX) $(CXXFLAGS) -c $< -o pac.o
	$(CXX) -o $@ pac.o MFloeAssembler.o
