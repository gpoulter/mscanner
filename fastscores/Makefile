.PHONY: all clean

F = _FeatureCounter
C = _ScoreCalculator

BINARIES = $(C) $(C).so $(C).dll $(F)

all: $(BINARIES)

clean:
	rm -f $(BINARIES) $(C).exe $(F).exe


$(F): $(F).c
	gcc -o $@ $<	

$(C): $(C).c
	gcc -DCSCORE -o $@ $<

$(C).so: $(C).c
	gcc -shared -o $@ $<

$(C).dll: $(C).c
	gcc -mno-cygwin -shared -o $@ $<

