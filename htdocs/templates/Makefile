.PHONY: all clean

all: front.py form.py status.py output.py

%.py: %.tmpl
	${CHEETAH} compile $<

front.py: page.py

form.py: page.py

status.py: page.py

output.py: page.py

clean:
	rm -f *.bak *.pyc
