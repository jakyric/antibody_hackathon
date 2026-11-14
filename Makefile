PY := .venv/bin/python
export PYTHONPATH := src

.PHONY: help setup test verify clean

help:
	@echo "make setup    create the virtual environment and install the requirements"
	@echo "make test     run the tests"
	@echo "make verify   run the benchmark: does the ranking recover answers we already know?"

# On a stock Debian/Ubuntu python, `python3 -m venv` creates the environment
# but leaves out pip: the distribution ships ensurepip in a separate
# python3-venv package. Every later step then fails with "No module named
# pip". Rather than tell the reader to install a system package, bootstrap pip
# into the environment when it is missing. Verified from a clean clone on a
# machine without python3.12-venv.
setup:
	python3 -m venv .venv 2>/dev/null || python3 -m venv --without-pip .venv
	@$(PY) -m pip --version >/dev/null 2>&1 || { \
	  echo "this python has no ensurepip; bootstrapping pip into .venv"; \
	  curl -sSf https://bootstrap.pypa.io/get-pip.py -o .venv/get-pip.py && \
	  $(PY) .venv/get-pip.py -q; }
	$(PY) -m pip -q install --upgrade pip
	$(PY) -m pip -q install -r requirements.txt

test:
	$(PY) -m pytest tests -q

# A complete run, end to end: propose variants, apply the entry rules against
# every provided sequence, drop what the model rates unlikely, rank on the one
# signal that transferred to measured affinity, and pick a diverse set. Writes
# to results/demo/ so it never overwrites the submitted designs.
# Needs the two Chinery tables in data/ -- see "Getting the data" in the README.
demo:
	$(PY) -m abx.pipeline --candidates 5000 --submit 10 --seed 7 \
	  --out results/demo/designs_demo.csv
	@echo
	@echo "Ranked candidates with their rationale: results/demo/designs_demo.csv"
	@echo "Then check them against the entry rules:"
	@echo "  $(PY) scripts/requirement_check.py --submission <submission.csv>"

verify:
	$(PY) -m abx.benchmark

clean:
	rm -rf .pytest_cache src/abx/__pycache__ tests/__pycache__
