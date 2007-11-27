.PHONY: all clean

F = _FeatureCounter
C = _ScoreCalculator

# Set up above using 16 and 32 bit features
F16 = $(F)16
F32 = $(F)32
C16 = $(C)16
C32 = $(C)32

GCC = gcc -lm
GCCL = $(GCC) -DLONGFEATS
BINARIES = $(C16) $(C32) $(F16) $(F32)
WINBINS = $(BINARIES:%=%.exe)

all: $(BINARIES)

clean:
	rm -f $(BINARIES) $(WINBINS)

$(F16): $(F).c
	$(GCC) -o $@ $<

$(F32): $(F).c	
	$(GCCL) -o $@ $<

$(C16): $(C).c
	$(GCC) -DCSCORE -o $@ $<

$(C32): $(C).c
	$(GCCL) -DCSCORE -o $@ $<

#$(C).so: $(C).c
#	$(GCC) -shared -o $@ $<

#$(C).dll: $(C).c
#	$(GCC) -mno-cygwin -shared -o $@ $<

