CC=c++
STANDARD=-std=c++11
COMPILE= $(STANDARD) -c -g -Wall
LINK= $(STANDARD)
all: rect pipe plotabkj

objects = timoshenko_wave_numbers.o wave_frequency_equation.o wave_number.o \
	wave_number_gradients.o frequency_equation.o clamped_clamped.o \
	clamped_free.o secular.o arclengthcontinuation.o

rect : rectangular_cross_section.o $(objects)
	$(CC) $(LINK) rectangular_cross_section.o $(objects) -o rect
pipe : circular_cross_section.o $(objects)
	$(CC) $(LINK) circular_cross_section.o $(objects) -o pipe
plotabkj: arclengthfig.o $(objects)
	$(CC) $(LINK) arclengthfig.o $(objects) -o plotabkj

rectangular_cross_section.o : rectangular_cross_section.cxx
	$(CC) $(COMPILE) rectangular_cross_section.cxx -I.
circular_cross_section.o : circular_cross_section.cxx
	$(CC) $(COMPILE) circular_cross_section.cxx -I.
arclengthfig.o : arclengthfig.cxx
	$(CC) $(COMPILE) arclengthfig.C -I.

clamped_clamped.o : clamped_clamped.C
	$(CC) $(COMPILE) clamped_clamped.C -I.
clamped_free.o : clamped_free.C
	$(CC) $(COMPILE) clamped_free.C -I.
frequency_equation.o : frequency_equation.C clamped_clamped.C clamped_free.C
	$(CC) $(COMPILE) frequency_equation.C -I.
arclengthcontinuation.o : arclengthcontinuation.C
	$(CC) $(COMPILE) arclengthcontinuation.C -I.
secular.o : secular.C
	$(CC) $(COMPILE) secular.C -I.
timoshenko_wave_numbers.o: timoshenko_wave_numbers.C
	$(CC) $(COMPILE) timoshenko_wave_numbers.C -I.
wave_number.o : wave_number.C 
	$(CC) $(COMPILE) wave_number.C -I. 
wave_number_gradients.o : wave_number_gradients.C
	$(CC) $(COMPILE) wave_number_gradients.C -I.
wave_frequency_equation.o : wave_frequency_equation.C
	$(CC) $(COMPILE) wave_frequency_equation.C -I.
clean:
	rm arclengthfig.o rectangular_cross_section.o circular_cross_section.o \
	frequency_equation.o clamped_clamped.o clamped_free.o \
	arclengthcontinuation.o secular.o \
	timoshenko_wave_numbers.o wave_number.o \
	wave_frequency_equation.o wave_number_gradients.o

clean_exec:
	rm plotabkj pipe rect
