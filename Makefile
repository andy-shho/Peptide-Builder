.PHONY: environment visualize

ENVIRONMENT=peptide_builder
VISUALIZE_SCRIPT = visualize.py

environment:
	conda remove --name $(ENVIRONMENT) --all --yes || true
	conda create -n $(ENVIRONMENT) "python=3.11" --yes
	conda install -c conda-forge numpy networkx gemmi matplotlib --name $(ENVIRONMENT) --yes

visualize:
	python ${VISUALIZE_SCRIPT} GIGAVLKVLTTGLPALISWIKRKRQQ melittin