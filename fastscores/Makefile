.PHONY: all clean

F = _FeatureCounter
C = _ScoreCalculator

GCC = gcc 
BINARIES = $(C) $(F)
WINBINS = $(BINARIES:%=%.exe)

all: $(BINARIES)

clean:
	rm -f $(BINARIES) $(WINBINS)
	
## Unix pipes approach

$(F): $(F).c
	$(GCC) -o $@ $< -lm

$(C): $(C).c
	$(GCC) -DCSCORE -o $@ $< -lm
	
## Shared library approach (but ctypes broken on Solaris)

#$(C).so: $(C).c
#	$(GCC) -shared -o $@ $< -lm

#$(C).dll: $(C).c
#	$(GCC) -mno-cygwin -shared -o $@ $< -lm

