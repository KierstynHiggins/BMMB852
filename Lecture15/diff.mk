# BMMB 852: Simulate and run RNA-Seq analysis

usage:
	@echo "# make activate (activates bioinfo)"
	@echo "# make tools (downloads SRC tools)"
	@echo "# make stat  (activates stats)"
	@echo "# make simulate (simulate design and counts)
	@echo "# make edger    (make edger file)"
	@echo "# make pca   (generate PCA plot)"
	@echo "# make heat    (make heat map)"
	@echo "# make all   (run all targets)"
	
activate:
	# Activate bioinfo environment.
	conda activate bioinfo

tools:
	# Ensure toolbox is downloaded.
	bio code

stat:
	# Activate stats environment.
	conda activate stats
simulate:
	# Simulate a design and counts file.
	Rscript src/r/simulate_counts.r

edger:
	# Run edger.
	Rscript src/r/edger.r
	Rscript  src/r/evaluate_results.r  -a counts.csv -b edger.csv

pca:
	# Generate pca plot.
	src/r/plot_pca.r -c edger.csv

heat:
	# Generate heat map.
	src/r/plot_heatmap.r -c edger.csv

all: activate tools stat simulate edger pca heat
