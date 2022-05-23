#compiler
CXX=g++ -std=c++0x

#compiler flags
COPTI= -O3
CDEBG= -g
CWARN= -Wall -Wextra -Wno-sign-compare

#linker flags
LOPTI= -O3
LDEBG= -g
LSTDD= -lm -lboost_iostreams -lboost_program_options -lz -lbz2 -lpthread -llzma

#executable file
EFILE= bin/toy

#header files
HFILE= $(shell find src -name *.h)

#source files
CFILE= $(shell find src -name *.cpp)

#source path
VPATH= $(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#include path
ISTDD= -Isrc -I./ -I/home/theo/commands/htslib-1.9

IBAMP= -Ilib

#object files
OFILE= $(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)

#default
all: dynamic

#dynamic release
dynamic: CFLAG=$(COPTI) $(CWARN)
dynamic: LFLAG=$(LOPTI) $(LSTDD)
dynamic: IFLAG=$(ISTDD)
dynamic: $(EFILE)


$(EFILE): $(OFILE)
	$(CXX) $^ $(HOME)/commands/htslib-1.9/libhts.a -o $@ $(LFLAG)

obj/%.o: %.cpp $(HFILE)
	$(CXX) -o $@ -c $< $(CFLAG) $(IFLAG)

clean: 
	rm -f obj/*.o $(EFILE)

oxford:
	cp $(EFILE) ~/bin/.

install:
	cp $(EFILE) /usr/local/bin/.

