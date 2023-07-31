################## WORK IN PROGRESS - probably not going to use it anymore ######################


# PDB_DIR = "../../datasets/pdb_structures_2"
# PDB_LST = $(PDB_DIR)/pdb_lst.txt
# #DIFF = $(diff <(ls $(PDB_DIR)) $(PDB_LST))
# CMAP_DIR = "../../datasets/cmaps"
# #SHELL := /bin/bash

# if the .checkpoint is here, the script doens't run if it doesn't find it
# I get nothing to be done even if it does things
# basically I could use a checkpoint as target




PDB_DIR := "pdb_structures_3"
AF_STR_CHKP := $(PDB_DIR)/.checkpoint

all: $(AF_STR_CHKP)

$(AF_STR_CHKP): $(PDB_DIR)
	@echo ""
	@echo "Rule 2"
	# Outputs all prerequisites
	@echo $^
	# Outputs all prerequisites newer than the target
	@echo $?
	@touch $@
	@echo ""

$(PDB_DIR): 
	@echo ""
	@echo "Rule 1"
	@mkdir -p $@ 
	@echo ""

# .PHONY: all



# PDB_DIR = "../../datasets/pdb_structures_2"
# SCRIPT_DIR = .
# AF_STR_CHKP = $(PDB_DIR)/.checkpoint


# $(AF_STR_CHKP): $(PDB_DIR)
# 	echo "Rule2"
# #	bash $(SCRIPT_DIR)/get_structures.sh $(PDB_DIR)
# 	touch $@

# $(PDB_DIR): 
# 	echo "Rule1"
# 	mkdir -p $@ 


# all: $(AF_STR_CHKP)

# .PHONY: all





# .PHONY: check_structures get_structures get_cmaps

# $(AF_STR_DIR): 
#    mkdir -p @



# get_structures: #$(PDB_DIR)/.checkpoint                                    
# 	@if [ ! -f $(PDB_DIR)/.checkpoint ]; then \
# 		bash get_structures.sh $(PDB_DIR); \
# 	fi

# check_structures:
# 	echo "pd"\;
# 	@if diff <(ls $(PDB_DIR)/*.pdb) $(PDB_LST) >/dev/null; then \ 
# 		echo "The files are the same!"; \
# 	else \
# 		echo "The files are different!"; \
# 		rm $(CMAP_DIR)/.checkpoint; \
# 		ls $(PDB_DIR)/*.pdb > $(PDB_DIR)/pdb_lst.txt
# 	fi

# get_cmaps: checkpoint_cmaps
# 	cd 
# 	runni $(SCRIPT_DIR)/cmaps script
# 	crea checkpoint_cmaps



# # Get structures from AF
# # Generate cmaps (depend on AF structures)
# # Generate df of DNA sequences (depend on AF structures and uniprot_gene_mapping)


# # Run
# run: datasets/cmaps/ datasets/af_uniprot_to_gene_id.json
# 	cd scripts
# 	python3 main.py -i ../tests/input/HARTWIG_WGS_PANCREAS_2020.in.maf -o ../tests/output/ \
# 	-p ../tests/input/HARTWIG_WGS_PANCREAS_2020.mutrate.json -H 0 -t PANCREAS -C HARTWIG_WGS_PANCREAS_2020