# Makefile by Cris Vergara
# Adapted from the one I made for school years ago

CPP=nvcc
OBJ_DIR=obj
INC_DIR=include
BIN_DIR=bin
SRC_DIR=src
CXXFLAGS=-I$(INC_DIR) -arch sm_52
CPPFLAGS=

SRC=$(shell ls $(SRC_DIR))
CUOBJ=$(SRC:.cu=.o)
OBJ=$(CUOBJ:.cpp=.o)

%.o: $(SRC_DIR)/%.cu
	if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi; $(CPP) $(CXXFLAGS) $(CPPFLAGS) $< -o $(OBJ_DIR)/$@ -dc

main: $(OBJ)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi; $(CPP) -lgmp $(CXXFLAGS) $(CPPFLAGS) $(addprefix $(OBJ_DIR)/,$(OBJ)) -o $(BIN_DIR)/$@

docs:
	doxygen;

cleandocs:
	if [ -d doc ]; then rm -rfv doc; fi;

clean:
	if [ -d $(BIN_DIR) ]; then rm -rfv $(BIN_DIR); fi; if [ -d $(OBJ_DIR) ]; then rm -rfv $(OBJ_DIR); fi;

all: main
