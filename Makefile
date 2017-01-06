SOURCES   = main.cpp luv.cpp
PROGNAME  = app
LIBRARIES = -lm -lGL -lGLU -lglut -lX11 -lcv -lhighgui -lcvaux
INCLUDES = /usr/include/opencv/

PROGNAME : main.cpp luv.cpp 
	g++ -g $(DFLAGS) -o app $^ -I/usr/include/opencv/ $(LIBRARIES)

clean:
	rm -f *.o core* *~ $(PROGNAME)
