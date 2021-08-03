main: main.cpp Menu.o SeqComparison.o Blosum.o
	g++ -std=c++11 main.cpp Blosum.o Menu.o SeqComparison.o -o prog.out -Wall -Wextra -O2

Blosum.o: Blosum.cpp Blosum.h
	g++ -std=c++11 -c Blosum.cpp -Wall -Wextra -O2


SeqComparison.o: SeqComparison.cpp SeqComparison.h
	g++ -std=c++11 -c SeqComparison.cpp -Wall -Wextra -O2


Menu.o: Menu.cpp Menu.h
	g++ -std=c++11 -c Menu.cpp -Wall -Wextra -O2


clean: 
	rm prog.out Menu.o SeqComparison.o Blosum.o