CC = mpic++
FC = mpif90

CFLAGS   := -g -O3
LDFLAGS  :=
INCLUDES := -I/usr/local/include -I./../../../test_thermal/submodule/monolis/include
LIBS     := -L/usr/local/lib -lm -lstdc++ -L../../../test_thermal/submodule/monolis/lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -lscalapack -llapack -lblas

INCLUDES_BB = -I../../../test_thermal/include
LIBS_BB     = -L../../../test_thermal/lib -lBBFE_std -lBBFE_sys -lBBFE_elemmat -lBBFE_manusol -lBB

INCLUDES_ROM = -I./../../include/rom_BB -I./../../include/rom_std -I./../../include/rom_sys -I./diff_core

LIBS_ROM     = -L./../../lib -lrom_BB -lrom_std -lrom_sys

SRC_DIR := ./diff_core
SRC     := $(wildcard $(addsuffix /*.c,$(SRC_DIR)))
OBJ1    := $(patsubst %.c,%.o,$(SRC)) ./diff_main/main.o

TARGET1 := hlpod_diff

# Default make target
all: build clean_data

# Build targets
build: $(TARGET1)

$(TARGET1): $(OBJ1)
	$(FC) $(CFLAGS) $(INCLUDES) $(INCLUDES_BB) $(INCLUDES_ROM) $(LDFLAGS) -o $@ $^ -Wl,--start-group $(LIBS_BB) $(LIBS_ROM) $(LIBS) -Wl,--end-group

# Compilation rule for .c to .o
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) $(INCLUDES_BB) $(INCLUDES_ROM) -c $< -o $@

.PHONY clean:;
	$(RM) -f $(OBJ1) $(TARGET1)

clean_data:
	$(RM) -f $(OBJ1)
