.PHONY: test regression

PYTHON ?= python3
# Assemblies for the end-to-end check; they are not in the repository.
GENOMES ?= T_oshimai.fasta T_thermophilus_HB27.fasta

test:
	@$(PYTHON) -m unittest discover -s tests

regression:
	@tests/regression.sh $(GENOMES)
