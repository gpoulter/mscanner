.PHONY: all clean

F = _FeatureCounter
C = _ScoreCalculator

GCC = gcc -lm
BINARIES = $(C) $(F)
WINBINS = $(BINARIES:%=%.exe)

all: $(BINARIES)

clean:
	rm -f $(BINARIES) $(WINBINS)
	
## Unix pipes approach

$(F): $(F).c
	$(GCC) -o $@ $<

$(C): $(C).c
	$(GCC) -DCSCORE -o $@ $<
	
## Shared library approach (but ctypes broken on Solaris)

#$(C).so: $(C).c
#	$(GCC) -shared -o $@ $<

#$(C).dll: $(C).c
#	$(GCC) -mno-cygwin -shared -o $@ $<

