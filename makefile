 
all: main rank select access memory 

main: main.o hybridId.o leafId.o hybridBV.o staticBV.o leafBV.o basics.o
	gcc -O9 -o main main.o hybridId.o leafId.o hybridBV.o staticBV.o leafBV.o basics.o

main.o: main.c hybridId.h leafId.h hybridBV.h staticBV.h leafBV.h
	gcc -O9 -c main.c

rank: rank.o hybridBV.o staticBV.o leafBV.o basics.o
	gcc -O9 -o rank rank.o hybridBV.o staticBV.o leafBV.o basics.o

rank.o: rank.c hybridBV.h staticBV.h leafBV.h
	gcc -O9 -c rank.c

select: select.o hybridBV.o staticBV.o leafBV.o basics.o
	gcc -O9 -o select select.o hybridBV.o staticBV.o leafBV.o basics.o

select.o: select.c hybridId.h leafId.h hybridBV.h staticBV.h leafBV.h
	gcc -O9 -c select.c

access: access.o hybridBV.o staticBV.o leafBV.o basics.o
	gcc -O9 -o access access.o hybridBV.o staticBV.o leafBV.o basics.o

access.o: access.c hybridBV.h staticBV.h leafBV.h
	gcc -O9 -c access.c

memory: memory.o hybridBV.o staticBV.o leafBV.o basics.o
	gcc -O9 -o memory memory.o hybridBV.o staticBV.o leafBV.o basics.o

memory.o: memory.c hybridBV.h staticBV.h leafBV.h
	gcc -O9 -c memory.c

hybridId.o: hybridId.c hybridId.h leafId.h basics.h
	gcc -O9 -c hybridId.c

leafId.o: leafId.c leafId.h basics.h
	gcc -O9 -c leafId.c

hybridBV.o: hybridBV.c hybridBV.h staticBV.h leafBV.h basics.h
	gcc -O9 -c hybridBV.c

staticBV.o: staticBV.c staticBV.h basics.h
	gcc -O9 -c staticBV.c

leafBV.o: leafBV.c leafBV.h basics.h
	gcc -O9 -c leafBV.c

basics.o: basics.c basics.h
	gcc -O9 -c basics.c

