CC=g++-9 -std=c++17 -fopenmp
DIR_INC=
CFLAGS=-Wall -Wno-unknown-pragmas -O3 $(DIR_INC)
DIR_SRC= ../src
CPP_FILES=$(wildcard ${DIR_SRC}/*.cpp) $(wildcard ${DIR_SRC}/*/*.cpp) $(wildcard ${DIR_SRC}/*/*/*.cpp)
SRC=$(CPP_FILES)
OBJ=$(SRC:.cpp=.o)

TARGET=../input/ILFT

defaut: $(TARGET)
	-rm $(OBJ)

$(TARGET): $(OBJ)
	$(CC) -o $@ $(OBJ)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY : clean
clean:
	rm $(TARGET) $(OBJ)