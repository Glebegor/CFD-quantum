PYTHON ?= python3
VENV := venv
VENV_PYTHON := $(VENV)/bin/python
LOCAL ?= false

ifeq ($(LOCAL),true)
RUN_ARGS := -local
else
RUN_ARGS :=
endif

i: $(VENV_PYTHON)
	$(VENV_PYTHON) -m pip install -r requirements.txt

$(VENV_PYTHON):
	$(PYTHON) -m venv $(VENV)

pf:
	$(VENV_PYTHON) -m pip freeze > requirements.txt

run: $(VENV_PYTHON)
	$(VENV_PYTHON) src/main.py -local

run-local: $(VENV_PYTHON)
	$(VENV_PYTHON) src/main.py -local

run-qpu: $(VENV_PYTHON)
	$(VENV_PYTHON) src/main.py

clean-series:
	rm -rf output/time_series

run-series-local: clean-series $(VENV_PYTHON)
	MPLBACKEND=Agg $(VENV_PYTHON) src/main.py -local --time-series

run-series-qpu: clean-series $(VENV_PYTHON)
	MPLBACKEND=Agg $(VENV_PYTHON) src/main.py --time-series

animation: $(VENV_PYTHON)
	$(VENV_PYTHON) src/animation.py output/time_series --time-step 0.1

run-cfd: $(VENV_PYTHON)
	MPLBACKEND=Agg $(VENV_PYTHON) src/lbm_solver.py

run-hybrid-aer: $(VENV_PYTHON)
	MPLBACKEND=Agg $(VENV_PYTHON) src/hybrid_qlbm.py

run-hybrid-qpu: $(VENV_PYTHON)
	MPLBACKEND=Agg $(VENV_PYTHON) src/hybrid_qlbm.py --qpu

show-circuits: $(VENV_PYTHON)
	$(VENV_PYTHON) src/gates.py

test-qlbm: $(VENV_PYTHON)
	$(VENV_PYTHON) -m unittest discover -s tests -v

.PHONY: i pf run run-local run-qpu clean-series run-series-local run-series-qpu animation run-cfd run-hybrid-aer run-hybrid-qpu show-circuits test-qlbm
