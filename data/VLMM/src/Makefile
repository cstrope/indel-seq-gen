
# c-code makefile

CC = gcc
CFLAGS = -O2
# add -Wall to see most warning messages
#     -g    to add debugging information
#     -O2   to add optimization

SRCS = *.c
OBJS = pattern_match.o priority_queue.o pst.o file_handle.o numrec.o pfa.o learn_algs.o

.PHONY: all clean wipe depend

all: train predict emit 2pfa

train: $(OBJS) train.o
	$(CC) $(CFLAGS) $(OBJS) train.o -o train -lm

predict: $(OBJS) predict.o
	$(CC) $(CFLAGS) $(OBJS) predict.o -o predict -lm

emit: $(OBJS) emit.o
	$(CC) $(CFLAGS) $(OBJS) emit.o -o emit -lm

2pfa: $(OBJS) 2pfa.o
	$(CC) $(CFLAGS) $(OBJS) 2pfa.o -o 2pfa -lm


clean:
	-\rm -f *~
	-\rm -f *.o
	-\rm -f core

wipe: clean
	-\rm -f train predict emit 2pfa


