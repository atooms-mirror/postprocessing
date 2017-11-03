PROJECT = postprocessing
COMMIT = $$(git describe --abbrev=6 --always 2>/dev/null || echo 0)
COMMIT_DIRTY = $$(git describe --abbrev=6 --dirty --always 2>/dev/null || echo 0)
DATE=$$(git show -s --format=%ci $(COMMIT) | cut -d ' ' -f 1)

.PHONY: all dist test install version clean

dist: version
	python setup.py sdist

test:	version
	python -m unittest discover -s tests

install: version
	python setup.py config_fc --opt '-O3 -funroll-loops' install --user

# This does not work with distutils
# develop:
# 	python setup.py develop --user

version:
	@echo __commit__ = \'$(COMMIT_DIRTY)\' > ${PROJECT}/_commit.py
	@echo __date__ = \'$(DATE)\' >> ${PROJECT}/_commit.py

clean:
	rm -f ${PROJECT}/*pyc  ${PROJECT}/*/*pyc tests/*pyc data/*.pp.*
