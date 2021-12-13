.jac128:
	nvc++ ldc128jac.cpp -acc -Minfo=accel -o ldc128jac

.sor128:
	nvc++ ldc128sor.cpp -acc -Minfo=accel -o ldc128sor

.jac256:
	nvc++ ldc256jac.cpp -acc -Minfo=accel -o ldc256jac

.sor256:
	nvc++ ldc256sor.cpp -acc -Minfo=accel -o ldc256sor

all: .jac128 .sor128 .jac256 .sor256

jac: .jac128 .jac256

sor: .sor128 .sor256

128: .jac128 .sor128

256: .jac256 .sor256

clean:
	rm ldc128jac ldc128sor ldc256jac ldc256sor