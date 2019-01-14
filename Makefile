CFLAGS = -O3 --std=c++11
CC = icpc
OMP = -qopenmp
VECREPORT = -qopt-report-phase=vec -qopt-report=5

reset: result-*.vtp result.pvd
	rm result-*.vtp result.pvd

%.o: %.c
	$(CC) $(CFLAGS) $(OMP) $(VECREPORT) -o $@ $<

step6-%: solution-step6.o
	./solution-step6.o `python3 params.py $*`
