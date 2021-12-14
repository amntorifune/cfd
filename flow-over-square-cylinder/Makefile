.jac:
	nvc++ fosc_jac.cpp -acc -Minfo=accel -o foscjac

.sor:
	nvc++ fosc_sor.cpp -acc -Minfo=accel -o foscsor

all: .jac .sor

jac: .jac

sor: .sor

clean:
	rm foscjac foscsor