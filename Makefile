PY := .venv/bin/python
export PYTHONPATH := src

.PHONY: help setup test verify clean

help:
	@echo "make setup    create the virtual environment and install the requirements"
	@echo "make test     run the tests"
	@echo "make verify   run the benchmark: does the ranking recover answers we already know?"

setup:
	python3 -m venv .venv
	$(PY) -m pip -q install --upgrade pip
	$(PY) -m pip -q install -r requirements.txt

test:
	$(PY) -m pytest tests -q

verify:
	$(PY) -m abx.benchmark

clean:
	rm -rf .pytest_cache src/abx/__pycache__ tests/__pycache__
