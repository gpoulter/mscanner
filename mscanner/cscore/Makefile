.PHONY: all clean

all: cscore cscore2.so cscore2.dll

clean:
	rm -f cscore cscore.exe cscore2.so cscore2.dll

cscore: cscore.c

cscore2.so: cscore2.c
	gcc -o cscore2.so -shared cscore2.c

cscore2.dll: cscore2.c
	gcc -o cscore2.dll -mno-cygwin -shared cscore2.c
