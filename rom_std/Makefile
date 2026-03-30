
CC = mpic++
CFLAGS = -Wall -Wextra -pedantic -O3
LDFLAGS =
INCLUDES = -I/usr/local/include -I../../test_thermal/include -I../../test_thermal/submodule/monolis/include -I../rom_BB
LIBS = -L/usr/local/lib -lm -L../../test_thermal/lib/BB -lBB -lBBFE_std

TARGET_LIB = librom_std.a
OBJS = hlpod_pre.o hlpod_comm.o hlpod_matvec.o hlpod_core.o set_modes.o monowrap.o write_std.o read_std.o

$(TARGET_LIB): $(OBJS)
	ar rcs $@ $?
	rm -f $(OBJS)

clean:
	-rm -f $(OBJS) $(TARGET_LIB)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(INCLUDES)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(INCLUDES)

.PHONY: clean;
	$(RM) -f $(OBJS) $(TARGET_LIB)

clean_data:
	$(RM) -f $(OBJS)
