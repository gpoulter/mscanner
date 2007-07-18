.PHONY: all clean

all: contact.py form.py front.py output.py status.py

%.py: %.tmpl
	${CHEETAH} compile $<

contact.py: page.py
form.py: page.py
front.py: page.py
output.py: page.py
status.py: page.py

clean:
	rm -f *.bak *.pyc
