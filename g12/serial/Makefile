CC = gcc
LD = gcc

SRC = src
LIB = lib
OUT = build

DMSG = 0

CFLAGS = -std=c17 -I. -pedantic-errors -Werror -Wall -Wextra -DMSG_LEVEL=$(DMSG) -O3

.PHONY: prepare clean program remake runall validate
remake: clean prepare program

clean:
	rm -rf $(OUT)

prepare:
	mkdir -p $(OUT)

program: $(OUT)/matrix.o $(OUT)/tsp.o $(OUT)/queue.o
	$(LD) -o tsp $(OUT)/matrix.o $(OUT)/tsp.o $(OUT)/queue.o -fopenmp

# Files
build/matrix.o: $(SRC)/matrix.c
	$(CC) $(CFLAGS) -o $(OUT)/matrix.o -c $(SRC)/matrix.c

build/tsp.o: $(SRC)/tsp.c
	$(CC) $(CFLAGS) -o $(OUT)/tsp.o -c $(SRC)/tsp.c -fopenmp

build/queue.o: $(LIB)/nqueue/queue.c
	$(CC) $(CFLAGS) -o $(OUT)/queue.o -c $(LIB)/nqueue/queue.c

runall:
	for filename in ../tests/*.in; do \
                fileout=`echo $$filename | cut -d'-' -f1`; \
                inputing=`echo $$filename | cut -d'-' -f2`; \
                inputing=`echo $$inputing | cut -d'.' -f1`; \
                echo $${filename} $${inputing}; \
                echo -n '' > $${fileout}.serial; \
                for number in 1 2 3 4 5; do \
                         echo -n "$${number} "; \
                         /usr/bin/time -f "%M KB" ./tsp $${filename} $${inputing} 2>> $${fileout}.serial >> $${fileout}.serial ;\
		done; \
                echo " "; \
        done

validate:
	for filename in tests/*.out; do \
		filediff=`echo $$filename | cut -d'.' -f1`; \
		echoing=`echo $$filediff | cut -d'/' -f2`; \
		echo $${echoing}; \
		diff $${filename} $${filediff}.test; \
	done
