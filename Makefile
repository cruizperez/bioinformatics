SHELL=/bin/bash

setup:
	poetry install \
	&& poetry run pip install -e .

lockfile:
	poetry lock --no-update

sortimport:
	isort --profile black src

style:
	black --line-length 99 src

test/unit:
	coverage run -m pytest -vv -r s src/tests/unit; coverage report -m

test/lint:
	flake8 src

test/sort:
	isort src --check --diff

test/type:
	mypy src/bioinformatics

test/style:
	black --line-length 120 --check src

ci: \
	test/unit \
	test/style \
	test/lint \
	test/sort \
	test/type
