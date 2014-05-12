CC = gcc
CFLAGS = -w -Wall -O3 -Wno-deprecated -static
OBJ  = ./main.o ./ctbparser.o ./crfparser.o ./crfparser_thread.o ./templet_feature.o ./base_feature.o ./user.o ./crf.o ./crf_thread.o ./lbfgs.o ./fun.o ./normalstr.o ./trie.o
LINKOBJ  = ./main.o ./ctbparser.o ./crfparser.o ./crfparser_thread.o ./templet_feature.o ./base_feature.o ./user.o ./crf.o ./crf_thread.o ./lbfgs.o ./fun.o ./normalstr.o ./trie.o
BIN  = ./ctbparser


LIBRARIES = -lm -lnsl -dl  -lstdc++ -lpthread

.PHONY: all all-before all-after clean clean-custom

all: ./ctbparser


clean: clean-custom
		rm -f $(OBJ) $(BIN)

$(BIN): $(LINKOBJ)
		$(CC) $(CFLAGS) $(LINKOBJ) -o $(BIN) $(LIBRARIES)

./main.o: main.cpp ctbparser.h
		$(CC) $(CFLAGS) -c ./main.cpp -o ./main.o

./ctbparser.o: ctbparser.cpp ctbparser.h const.h fun.h
		$(CC) $(CFLAGS) -c ./ctbparser.cpp -o ./ctbparser.o

./crfparser.o: crfparser.cpp crfparser.h const.h lbfgs.h
		$(CC) $(CFLAGS) -c ./crfparser.cpp -o ./crfparser.o

./crfparser_thread.o: crfparser_thread.cpp crfparser_thread.h const.h fun.h
		$(CC) $(CFLAGS) -c ./crfparser_thread.cpp -o ./crfparser_thread.o

./templet_feature.o: templet_feature.cpp templet_feature.h const.h fun.h dat.h
		$(CC) $(CFLAGS) -c ./templet_feature.cpp -o ./templet_feature.o

./base_feature.o: base_feature.cpp base_feature.h
		$(CC) $(CFLAGS) -c ./base_feature.cpp -o ./base_feature.o

./user.o: user.cpp crfparser.h templet_feature.h
		$(CC) $(CFLAGS) -c ./user.cpp -o ./user.o

./crf_thread.o: crf_thread.cpp crf_thread.h fun.h
		$(CC) $(CFLAGS) -c ./crf_thread.cpp -o crf_thread.o

./crf.o: crf.cpp crf.h fun.h lbfgs.h const.h dat.h
		$(CC) $(CFLAGS) -c crf.cpp -o crf.o

./lbfgs.o: lbfgs.cpp lbfgs.h
		$(CC) $(CFLAGS) -c lbfgs.cpp -o lbfgs.o

./fun.o: fun.cpp fun.h const.h
		$(CC) $(CFLAGS) -c fun.cpp -o fun.o

./normalstr.o: normalstr.cpp normalstr.h const.h
		$(CC) $(CFLAGS) -c ./normalstr.cpp -o ./normalstr.o

./trie.o: trie.cpp trie.h
		$(CC) $(CFLAGS) -c ./trie.cpp -o ./trie.o

