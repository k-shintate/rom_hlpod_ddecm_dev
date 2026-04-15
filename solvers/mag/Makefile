CC = mpic++
FC = mpif90

CFLAGS   := -g -O3
LDFLAGS  :=
INCLUDES := -I/usr/local/include -I./../../../test_thermal/submodule/monolis/include
LIBS     := -L/usr/local/lib -lm -lstdc++ -L../../../test_thermal/submodule/monolis/lib -lmonolis_solver -lgedatsu -lmonolis_utils -lmetis -lscalapack -llapack -lblas

INCLUDES_BB = -I../../../test_thermal/include
LIBS_BB     = -L../../../test_thermal/lib -lBBFE_std -lBBFE_sys -lBBFE_elemmat -lBBFE_manusol -lBB

INCLUDES_ROM = -I./../../include/rom_BB -I./../../include/rom_std -I./../../include/rom_sys -I./../../include/rom_ecm -I./mag_core -I./mag_matvec

LIBS_ROM     = -L./../../lib -lrom_BB -lrom_std -lrom_sys -lrom_ecm

SRC_DIR := ./mag_core ./mag_matvec
SRC     := $(wildcard $(addsuffix /*.c,$(SRC_DIR)))
OBJ1    := $(patsubst %.c,%.o,$(SRC)) ./mag_main/mag.o
OBJ2    := $(patsubst %.c,%.o,$(SRC)) ./mag_main/mag_ns.o
OBJ3    := $(patsubst %.c,%.o,$(SRC)) ./mag_main/mag_threephase.o
OBJ4    := $(patsubst %.c,%.o,$(SRC)) ./mag_main/mag_threephase_online.o

TARGET1 := mag
TARGET2 := mag_ns
TARGET3 := mag_threephase
TARGET4 := mag_threephase_online
# Default make target
all: build clean_data

# Build targets
build: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

$(TARGET1): $(OBJ1)
	$(FC) $(CFLAGS) $(INCLUDES) $(INCLUDES_BB) $(INCLUDES_ROM) $(LDFLAGS) -o $@ $^ -Wl,--start-group $(LIBS_BB) $(LIBS_ROM) $(LIBS) -Wl,--end-group

$(TARGET2): $(OBJ2)
	$(FC) $(CFLAGS) $(INCLUDES) $(INCLUDES_BB) $(INCLUDES_ROM) $(LDFLAGS) -o $@ $^ -Wl,--start-group $(LIBS_BB) $(LIBS_ROM) $(LIBS) -Wl,--end-group

$(TARGET3): $(OBJ3)
	$(FC) $(CFLAGS) $(INCLUDES) $(INCLUDES_BB) $(INCLUDES_ROM) $(LDFLAGS) -o $@ $^ -Wl,--start-group $(LIBS_BB) $(LIBS_ROM) $(LIBS) -Wl,--end-group

$(TARGET4): $(OBJ4)
	$(FC) $(CFLAGS) $(INCLUDES) $(INCLUDES_BB) $(INCLUDES_ROM) $(LDFLAGS) -o $@ $^ -Wl,--start-group $(LIBS_BB) $(LIBS_ROM) $(LIBS) -Wl,--end-group

# Compilation rule for .c to .o
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) $(INCLUDES_BB) $(INCLUDES_ROM) -c $< -o $@

.PHONY clean:;
	$(RM) -f $(OBJ1) $(TARGET1) $(OBJ2) $(TARGET2) $(OBJ3) $(TARGET3) $(OBJ4) $(TARGET4)

clean_data:
	$(RM) -f $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4)
