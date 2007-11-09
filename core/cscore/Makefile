.PHONY: all clean

all: cscore cscore.so cscore.dll

clean:
	rm -f cscore cscore.exe cscore2.so cscore2.dll

cscore: cscore.c
	gcc -o cscore -DCSCORE cscore.c

cscore.so: cscore.c
	gcc -o cscore.so -shared cscore.c

cscore.dll: cscore.c
	gcc -o cscore.dll -mno-cygwin -shared cscore.c
