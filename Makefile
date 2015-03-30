# Standard Linux, performance
CC = g++
CFLAGS = -g -lm -Wall -O3 -Wextra -pedantic -std=c++11 -msse -msse2 -msse3 -mfpmath=sse -march=native -mtune=native -funroll-loops -flto #-m64

# Mac OS
# CC = g++-4.9
# CFLAGS = -lm -Wall -O3 -std=c++11 


# Windows
# CC = g++
# CFLAGS = x86_64-w64-mingw32-g++ -lm -Wall -O3 -msse -msse2 -msse3 -mfpmath=sse -funroll-loops -flto -m64 -std=gnu++11 -static -static-libgcc -static-libstdc++

# Debugging
# CC = g++
# CFLAGS = g++ -g -lm -Wall -Wextra -pedantic -std=c++11 -msse -msse2 -msse3


HEADER = $(wildcard src/*.h)
CPP = $(wildcard src/*.cpp)
OBJ = $(addprefix obj/,$(notdir $(CPP:.cpp=.o)))

default: bin/neurofield

bin/neurofield: $(OBJ)
	@mkdir -p bin
	@echo "$(CC) $(CFLAGS) $(OBJ) -o $@"
	@$(CC) $(CFLAGS) $(OBJ) -o $@ || (echo "mycommand failed $$?"; exit 1)
	@echo "====="
	@cat license.txt
	@echo "====="
	@echo "USE OF NEUROFIELD CONSTITUTES ACCEPTANCE OF THE LICENSE CONDITIONS ABOVE"

-include $(OBJ:.o=.d)

obj/%.o: src/%.cpp
	@mkdir -p obj
	@$(CC) $(CFLAGS) -c $< -o $@
	@$(CC) -MM $(CFLAGS) $< > obj/$*.d
	@echo "CC $<"
	@mv -f obj/$*.d obj/$*.d.tmp
	@sed -e 's|.*:|$@:|' < obj/$*.d.tmp > obj/$*.d
	@sed -e 's/.*://' -e 's/\\$$//' < obj/$*.d.tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> obj/$*.d
	@rm -f obj/$*.d.tmp

Documentation/user.pdf: Documentation/user.tex
	cd Documentation && pdflatex user && pdflatex user

Documentation/developer.pdf: Documentation/developer.tex
	cd Documentation && pdflatex developer && pdflatex developer

Paper/neurofield.pdf: Paper/neurofield.tex
	cd Paper && pdflatex neurofield #&& pdflatex neurofield

.PHONY: clean doc paper

doc: Documentation/user.pdf Documentation/developer.pdf

paper: Paper/neurofield.pdf

clean:
	@-rm -rf bin obj Documentation/{user,developer}.{aux,log,out,toc} Documentation/x.log
