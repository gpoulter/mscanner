.PHONY: all clean

all: cscore

cscore: cscore.c

cscore2.dll: cscore2.c
	gcc -o cscore2.dll -mno-cygwin -shared cscore2.c

clean:
	rm -f cscore.exe cscore
