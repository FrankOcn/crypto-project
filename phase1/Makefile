# Makefile by Cris Vergara
# Adapted from the one I made for school years ago

CC=gcc
OBJ_DIR=obj
INC_DIR=include
BIN_DIR=bin
SRC_DIR=src
CCFLAGS=-g -lgmp -lpthread -I$(INC_DIR)

SRC=$(shell ls $(SRC_DIR))
OBJ=$(SRC:.c=.o)

%.o: $(SRC_DIR)/%.c
	if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi; $(CC) $(CCFLAGS) -c $< -o $(OBJ_DIR)/$@

main: $(OBJ)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi; $(CC) $(CCFLAGS) $(addprefix $(OBJ_DIR)/,$(OBJ)) main.c -o $(BIN_DIR)/$@

docs:
	doxygen;

cleandocs:
	if [ -d doc ]; then rm -rfv doc; fi;

clean:
	if [ -d $(BIN_DIR) ]; then rm -rfv $(BIN_DIR); fi; if [ -d $(OBJ_DIR) ]; then rm -rfv $(OBJ_DIR); fi;

all: main
