PYTHON ?= python3
VENV := venv
VENV_PYTHON := $(VENV)/bin/python

run: $(VENV_PYTHON)
	$(VENV_PYTHON) main.py

i: $(VENV_PYTHON)
	$(VENV_PYTHON) -m pip install -r requirements.txt

$(VENV_PYTHON):
	$(PYTHON) -m venv $(VENV)

# Activation must happen in the user's current shell. Run:
#   source venv/bin/activate
venv: $(VENV_PYTHON)
	@echo "Run: source $(VENV)/bin/activate"

.PHONY: run i venv
