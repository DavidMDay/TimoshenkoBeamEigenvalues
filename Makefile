CC=c++
STANDARD=-std=c++11
COMPILE= $(STANDARD) -c -g -Wall
LINK= $(STANDARD)
all: rect pipe plotabkj

objects = timoshenko_wave_numbers.o wave_frequency_equation.o wave_number.o \
	wave_number_gradients.o frequency_equation.o frequency_equation_cc.o \
	frequency_equation_cf.o secular.o pseudoarclength.o

rect : rectangular_timoshenko.o $(objects)
	$(CC) $(LINK) rectangular_timoshenko.o $(objects) -o rect
pipe : pipe_timoshenko.o $(objects)
	$(CC) $(LINK) pipe_timoshenko.o $(objects) -o pipe
plotabjk: arclengthfig.o $(objects)
	$(CC) $(LINK) arclengthfig.o $(objects) -o figure
arclengthfig.o : arclengthfig.C
	$(CC) $(COMPILE) arclengthfig.C -I.
pipe_timoshenko.o : pipe_timoshenko.C
	$(CC) $(COMPILE) pipe_timoshenko.C -I.
rectangular_timoshenko.o : rectangular_timoshenko.C
	$(CC) $(COMPILE) rectangular_timoshenko.C -I.
frequency_equation_cc.o : frequency_equation_cc.C
	$(CC) $(COMPILE) frequency_equation_cc.C -I.
frequency_equation_cf.o : frequency_equation_cf.C
	$(CC) $(COMPILE) frequency_equation_cf.C -I.
frequency_equation.o : frequency_equation.C frequency_equation_cc.C frequency_equation_cf.C
	$(CC) $(COMPILE) frequency_equation.C -I.
pseudoarclength.o : pseudoarclength.C
	$(CC) $(COMPILE) pseudoarclength.C -I.
secular.o : secular.C
	$(CC) $(COMPILE) secular.C -I.
wave_number.o : wave_number.C 
	$(CC) $(COMPILE) wave_number.C -I. 
wave_number_gradients.o : wave_number_gradients.C
	$(CC) $(COMPILE) wave_number_gradients.C -I.
wave_frequency_equation.o : wave_frequency_equation.C
	$(CC) $(COMPILE) wave_frequency_equation.C -I.
timoshenko_wave_numbers.o: timoshenko_wave_numbers.C
	$(CC) $(COMPILE) timoshenko_wave_numbers.C -I.
clean:
	rm arclengthfig.o frequency_equation.o frequency_equation_cc.o \
	frequency_equation_cf.o pseudoarclength.o secular.o \
	timoshenko_wave_numbers.o wave_frequency_equation.o wave_number.o \
	wave_number_gradients.o
