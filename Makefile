SOURCES   = Principal.cpp luv.cpp
PROGNAME  = app
LIBRARIES = -lm -lGL -lGLU -lglut -lX11 -lcv -lhighgui -lcvaux
INCLUDES = /usr/include/opencv/

PROGNAME : Principal.cpp luv.cpp 
	g++ -o app Principal.cpp luv.cpp -I/usr/include/opencv/ $(LIBRARIES)

all:	$(PROGNAME) 

clean:
	rm -f *.o core* *~ $(PROGNAME)
