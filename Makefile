all: logrank_test logrank.o

logrank_test: logrank_test.o logrank.o
	g++ -o logrank_test logrank_test.o logrank.o
logrank_test.o: logrank_test.cpp
	g++ -c logrank_test.cpp

logrank.o: logrank.cpp
	g++ -c logrank.cpp
