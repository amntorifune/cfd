.128:
	nvc++ ns2dccg.cpp -acc -Minfo=accel -o ns2dccg
.256:
	nvc++ ns2dccg256.cpp -acc -Minfo=accel -o ns2dccg256
all: .128 .256
128: .128
256: .256
clean:
	rm ns2dccg
	rm ns2dccg256