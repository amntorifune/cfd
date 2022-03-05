#!/usr/bin/bash

date > build.txt

echo nvc++ -c geom.cpp    -acc -Minfo=accel -o ./obj/geom.o    >> build.txt
     nvc++ -c geom.cpp    -acc -Minfo=accel -o ./obj/geom.o   &>> build.txt
echo ___________________________________________________       >> build.txt
echo nvc++ -c init.cpp    -acc -Minfo=accel -o ./obj/init.o    >> build.txt
     nvc++ -c init.cpp    -acc -Minfo=accel -o ./obj/init.o   &>> build.txt
echo ___________________________________________________       >> build.txt
echo nvc++ -c bcu.cpp     -acc -Minfo=accel -o ./obj/bcu.o     >> build.txt
     nvc++ -c bcu.cpp     -acc -Minfo=accel -o ./obj/bcu.o    &>> build.txt
echo ___________________________________________________       >> build.txt
echo nvc++ -c bcp.cpp     -acc -Minfo=accel -o ./obj/bcp.o     >> build.txt
     nvc++ -c bcp.cpp     -acc -Minfo=accel -o ./obj/bcp.o    &>> build.txt
echo ___________________________________________________       >> build.txt
echo nvc++ -c solid.cpp   -acc -Minfo=accel -o ./obj/solid.o   >> build.txt
     nvc++ -c solid.cpp   -acc -Minfo=accel -o ./obj/solid.o  &>> build.txt
echo ___________________________________________________       >> build.txt
echo nvc++ -c fs1.cpp     -acc -Minfo=accel -o ./obj/fs1.o     >> build.txt
     nvc++ -c fs1.cpp     -acc -Minfo=accel -o ./obj/fs1.o    &>> build.txt
echo ___________________________________________________       >> build.txt
echo nvc++ -c contra.cpp  -acc -Minfo=accel -o ./obj/contra.o  >> build.txt
     nvc++ -c contra.cpp  -acc -Minfo=accel -o ./obj/contra.o &>> build.txt
echo ___________________________________________________       >> build.txt
echo nvc++ -c diver.cpp   -acc -Minfo=accel -o ./obj/diver.o   >> build.txt
     nvc++ -c diver.cpp   -acc -Minfo=accel -o ./obj/diver.o  &>> build.txt
echo ___________________________________________________       >> build.txt
echo nvc++ -c sor.cpp     -acc -Minfo=accel -o ./obj/sor.o     >> build.txt
     nvc++ -c sor.cpp     -acc -Minfo=accel -o ./obj/sor.o    &>> build.txt
echo ___________________________________________________       >> build.txt
echo nvc++ -c jacob.cpp   -acc -Minfo=accel -o ./obj/jacob.o   >> build.txt
     nvc++ -c jacob.cpp   -acc -Minfo=accel -o ./obj/jacob.o  &>> build.txt
echo ___________________________________________________       >> build.txt
echo nvc++ -c avp0.cpp    -acc -Minfo=accel -o ./obj/avp0.o    >> build.txt
     nvc++ -c avp0.cpp    -acc -Minfo=accel -o ./obj/avp0.o   &>> build.txt
echo ___________________________________________________       >> build.txt
echo nvc++ -c fs2.cpp     -acc -Minfo=accel -o ./obj/fs2.o     >> build.txt
     nvc++ -c fs2.cpp     -acc -Minfo=accel -o ./obj/fs2.o    &>> build.txt
echo ___________________________________________________       >> build.txt
echo nvc++ -c fio.cpp     -acc -Minfo=accel -o ./obj/fio.o     >> build.txt
     nvc++ -c fio.cpp     -acc -Minfo=accel -o ./obj/fio.o    &>> build.txt
echo ___________________________________________________       >> build.txt
echo nvc++ -c main.cpp    -acc -Minfo=accel -o ./obj/main.o    >> build.txt
     nvc++ -c main.cpp    -acc -Minfo=accel -o ./obj/main.o   &>> build.txt
echo ___________________________________________________       >> build.txt
echo nvc++ ./obj/geom.o ./obj/init.o ./obj/bcu.o ./obj/bcp.o ./obj/solid.o ./obj/fs1.o ./obj/contra.o ./obj/diver.o ./obj/sor.o ./obj/jacob.o ./obj/avp0.o ./obj/fs2.o ./obj/fio.o ./obj/main.o -acc -Minfo=accel -o 3dturb  >> build.txt
nvc++ ./obj/geom.o ./obj/init.o ./obj/bcu.o ./obj/bcp.o ./obj/solid.o ./obj/fs1.o ./obj/contra.o ./obj/diver.o ./obj/sor.o ./obj/jacob.o ./obj/avp0.o ./obj/fs2.o ./obj/fio.o ./obj/main.o -acc -Minfo=accel -o 3dturb &>> build.txt