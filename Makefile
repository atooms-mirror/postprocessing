PROJECT = postprocessing

.PHONY: all dist test install clean

dist:
	python setup.py sdist

test:
	python -m unittest discover -s tests

wrap:
	cd postprocessing; f2py -c -m neighbors_wrap neighbors.f90

install: wrap
	python setup.py install --user

develop:
	python setup.py develop --user

clean:
	rm -f ${PROJECT}/*pyc  ${PROJECT}/*/*pyc tests/*pyc
