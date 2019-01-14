CFLAGS = -O3 --std=c++11 -debug full
CC = icpc
OMP = -qopenmp
VECREPORT = -qopt-report-phase=vec -qopt-report=5

reset: result-*.vtp result.pvd
	rm result-*.vtp result.pvd

%.o: %.c
	$(CC) $(CFLAGS) $(OMP) $(VECREPORT) -o $@ $<

step6-%: solution-step6.o
	./solution-step6.o `python2 params.py $*`
