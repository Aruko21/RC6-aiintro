# Program sources
SRCS=src/* main.cpp

# Program headers
HDRS=include/*

# Include directory
INC=-Iinclude

# Gnuplot
GNPLT_GIF=gnuplot_gif.gp
GNPLT_WINDOW=gnuplot_window.gp

# Compiler
CC=g++

CXX_FLAGS=-std=c++11

EXEC=lab3

.PHONY: clean

$(EXEC): $(SRCS) $(HDRS)
	$(CC) $(INC) -o $(EXEC) $(SRCS) $(CXX_FLAGS)

plot_gif:
	gnuplot gnuplot_gif.gp

plot_window:
	gnuplot gnuplot_window.gp

clean:
	rm $(EXEC) *.gif *.dat
