.jac:
	nvc++ step_jac.cpp -acc -Minfo=accel -o stepjac

.sor:
	nvc++ step_sor.cpp -acc -Minfo=accel -o stepsor

all: .jac .sor

jac: .jac

sor: .sor

clean:
	rm stepjac stepsor