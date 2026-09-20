CC = gcc
SRC = nii2tvx.c
LIBS = -lz -lm

.PHONY: all sanitize clean

all:
	$(CC) -O3 $(SRC) -o nii2tvx $(LIBS)

sanitize:
	$(CC) -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer $(SRC) -o nii2tvx_asan $(LIBS)

clean:
	rm -f nii2tvx nii2tvx_asan
