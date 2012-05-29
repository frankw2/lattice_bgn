CC = g++
CFLAGS = -g -Wall -O2
LDFLAGS = -lntl -lgmp -lm -lgsl

SOURCES = test.c bgn.c
TARGET = test
OBJECTS =  $(SOURCES:.c=.o)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJECTS)
	
.PHONY: clean
clean:
	rm -fr $(TARGET) $(OBJECTS)
