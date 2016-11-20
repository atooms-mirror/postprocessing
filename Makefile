PROJECT = postprocessing

.PHONY: all dist test version install clean

dist:
	python setup.py sdist

test:
	python -m unittest discover -s tests

install: version
	python setup.py install --user

develop:
	python setup.py develop --user

clean:
	rm -f $PROJECT/*pyc  $PROJECT/*/*pyc tests/*pyc
