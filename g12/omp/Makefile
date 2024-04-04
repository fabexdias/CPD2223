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

program: $(OUT)/matrix.o $(OUT)/tsp-omp.o $(OUT)/queue.o
	$(LD) -o tsp-omp $(OUT)/matrix.o $(OUT)/tsp-omp.o $(OUT)/queue.o -fopenmp

# Files
build/matrix.o: $(SRC)/matrix.c
	$(CC) $(CFLAGS) -o $(OUT)/matrix.o -c $(SRC)/matrix.c

build/tsp-omp.o: $(SRC)/tsp-omp.c
	$(CC) $(CFLAGS) -o $(OUT)/tsp-omp.o -c $(SRC)/tsp-omp.c -fopenmp

build/queue.o: $(LIB)/nqueue/queue.c
	$(CC) $(CFLAGS) -o $(OUT)/queue.o -c $(LIB)/nqueue/queue.c

runall:
	make; \
	for filename in ../tests/gen19-23.in ../tests/gen20-5000.in ../tests/gen22-25000.in ../tests/gen24-50000.in ../tests/gen26-50000.in ../tests/gen30-5000.in ; do \
		fileout=`echo $$filename | cut -d'-' -f1`; \
		inputing=`echo $$filename | cut -d'-' -f2`; \
		inputing=`echo $$inputing | cut -d'.' -f1`; \
		echo $${filename} $${inputing}; \
		for threads in 1 2 4 8; do \
			export OMP_NUM_THREADS=$${threads}; \
			echo "====$${OMP_NUM_THREADS}===="; \
			echo -n '' > $${fileout}.$${OMP_NUM_THREADS}; \
			for number in 1 2 3 4 5 6 7 8 9 10; do \
				echo -n "$${number} "; \
				/usr/bin/time -f "%M KB" ./tsp $${filename} $${inputing} 2>> $${fileout}.$${OMP_NUM_THREADS} >> $${fileout}.$${OMP_NUM_THREADS}; \
			done; \
			echo " "; \
		done \
	done

validate:
	for filename in tests/*.out; do \
		filediff=`echo $$filename | cut -d'.' -f1`; \
		echoing=`echo $$filediff | cut -d'/' -f2`; \
		echo $${echoing}; \
		diff $${filename} $${filediff}.test; \
	done
