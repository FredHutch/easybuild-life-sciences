PROG=../../bin/vicuna_analysis

CC=g++
CFLAGS=-c -Wall -O3  
LDFLAGS= #
SRCS=$(wildcard *.cpp) 
OBJS=$(SRCS:.cpp=.o)

all: $(SRCS) $(PROG)
		
$(PROG): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $@


$(OBJS): %.o : %.cpp
	$(CC) $(CFLAGS) $<   

clean:
	rm -rf $(PROG) $(OBJS)
