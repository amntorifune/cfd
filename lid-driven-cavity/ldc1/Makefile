.128:
	nvc++ ns2d.cpp -acc -Minfo=accel -o ns2d
.256:
	nvc++ ns2d256.cpp -acc -Minfo=accel -o ns2d256
all: .128 .256
128: .128
256: .256
clean:
	rm ns2d
	rm ns2d256