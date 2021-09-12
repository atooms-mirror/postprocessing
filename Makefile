PROJECT=atooms/postprocessing
PACKAGE=atooms.postprocessing

.PHONY: all version test coverage autopep8 flake8 install debug clean

all: install

# Version
COMMIT = $$(git describe --abbrev=6 --always 2>/dev/null || echo 0)
COMMIT_DIRTY = $$(git describe --abbrev=6 --dirty --always 2>/dev/null || echo 0)
DATE=$$(git show -s --format=%ci $(COMMIT) | cut -d ' ' -f 1)

version:
	@echo __commit__ = \'$(COMMIT_DIRTY)\' > ${PROJECT}/_commit.py
	@echo __date__ = \'$(DATE)\' >> ${PROJECT}/_commit.py

# Tests with fortran extension
test:
	mv $(PROJECT) $(PROJECT).tmp
	python -m unittest discover -s tests; mv $(PROJECT).tmp $(PROJECT)

coverage: install
	mv $(PROJECT) $(PROJECT).tmp
	coverage run --source $(PACKAGE) -m unittest discover -s tests; coverage report; mv $(PROJECT).tmp $(PROJECT)

# PEP8
autopep8:
	autopep8 -r -i $(PROJECT)
	autopep8 -r -i tests

flake8:
	flake8 $(PROJECT)

# Install with fortran extension
install: version
	python setup.py config_fc --quiet --opt '-O3 -funroll-loops' install

debug: version
	python setup.py config_fc --quiet --opt '-O3 -funroll-loops -fbounds-check' install

clean:
	find $(PROJECT) tests -name '*.pyc' -name '*.so' -exec rm '{}' +
	find $(PROJECT) tests -name -name '__pycache__' -exec rm -r '{}' +
	rm -rf build/ dist/

