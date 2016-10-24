PROJECT = sample

.PHONY: all dist test version install clean

all:
	cd $PROJECT; make; cd ..

dist:
	python setup.py sdist

test:
	python -m unittest discover -s tests

install: version
	python setup.py install --home=~

clean:
	rm -f $PROJECT/*pyc  $PROJECT/*/*pyc tests/*pyc
