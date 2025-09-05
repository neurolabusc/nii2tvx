CC = gcc
SRC = nii2tvx.c
OUT = nii2tvx
LIBS = -lz -lm

all:
	$(CC) -O3 $(SRC) -o $(OUT) $(LIBS)

sanitize:
	$(CC) -O1 -g -fsanitize=address -fno-omit-frame-pointer $(SRC) -o $(OUT) $(LIBS)

clean:
	rm -f $(OUT)
