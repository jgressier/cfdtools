SRC := .
VENV_DIR := $(SRC)/venv
REQUIREMENTS_FILE := $(SRC)/requirements.txt
PKG := cfdtools

#env:
#	mkvirtualenv --python=$(which python3.7) $(PKG)

.PHONY: help venv

help: ## print this help
	@echo "\n$$(poetry version): use target ; available Makefile following targets"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' Makefile | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-16s\033[0m %s\n", $$1, $$2}'

venv: $(REQUIREMENTS_FILE)
	@echo "Creating virtual environment in $(VENV_DIR)..."
	$(PYTHON) -m venv $(VENV_DIR)
	@echo "Installing dependencies from $(REQUIREMENTS_FILE)..."
	$(VENV_DIR)/bin/pip install --upgrade pip
	$(VENV_DIR)/bin/pip install -r $(REQUIREMENTS_FILE)
	$(VENV_DIR)/bin/pip install -e .
	@echo "Virtual environment is ready."

.PHONY: install
install: ## install minimum required packages and flowdyn to local
	pip install -r $(SRC)/requirements.txt
	pip install $(SRC)

install_dev: install ## install package for development and testing
	#<missing># pip install -r $(SRC)/requirements-dev.txt
	#<missing># pip install -r $(SRC)/docs/requirements.txt

check_pyproject: ## check all requirements are defined in pyproject.toml
	cat requirements.txt | grep -E '^[^# ]' | cut -d= -f1 | xargs -n 1 poetry add
	cat requirements-dev.txt | grep -E '^[^# ]' | cut -d= -f1 | xargs -n 1 poetry add -D
	cat docs/requirements.txt | grep -E '^[^# ]' | cut -d= -f1 | xargs -n 1 poetry add -D

test: install_dev ## run tests with pytest
	pytest

cov_run:
	poetry run pytest --cov-report=xml

cov_publish: .codecov_token
	CODECOV_TOKEN=$$(cat .codecov_token)  bash <(curl -s https://codecov.io/bash)

clean: ## clean all unnecessary files
	find . -name "__pycache__" -exec rm -rf {} +
	find . -name ".mypy_cache" -exec rm -rf {} +
	find . -name ".pytest_cache" -exec rm -rf {} +
	find . -name ".coverage" -exec rm -f {} +
	find . -name ".ipynb_checkpoints" -exec rm -f {} +

clean_notebooks: ## remove ouputs in Jupyter notebooks files
	find lessons -name \*.ipynb -exec python3 scripts/remove_output.py {} +

clean_all: clean_notebooks clean ## run clean and clean_notebooks

