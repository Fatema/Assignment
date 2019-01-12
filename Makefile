CFLAGS = -O3 --std=c++11
CC = icpc

reset: result-*.vtp result.pvd
	rm result-*.vtp result.pvd

%.o: %.c
	$(CC) $(CFLAGS) -o $@ $<
