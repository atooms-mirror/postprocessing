PROJECT = postprocessing

.PHONY: all dist test install clean

dist:
	python setup.py sdist

test:
	python -m unittest discover -s tests

install:
	python setup.py config_fc --opt '-O3 -funroll-loops' install --user

# This does not work with distutils
# develop:
# 	python setup.py develop --user

clean:
	rm -f ${PROJECT}/*pyc  ${PROJECT}/*/*pyc tests/*pyc data/*.pp.*
