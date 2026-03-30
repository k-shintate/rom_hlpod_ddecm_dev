
CC = mpic++
CFLAGS = -Wall -Wextra -pedantic -O3
LDFLAGS =
INCLUDES = -I/usr/local/include -I../../test_thermal/include -I../../test_thermal/submodule/monolis/include -I../rom_BB -I../rom_std -I../rom_sys
LIBS = -L/usr/local/lib -lm -L../../test_thermal/lib/BB -lBB -lBBFE_std -lBBFE_sys -lBBFE_elemmat

TARGET_LIB = librom_ecm.a
OBJS = HR.o DDHR.o DDHR_para.o DDHR_para_lb.o ecm_write.o inc_svd.o

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
