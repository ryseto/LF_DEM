test:
		@echo "\n Building...\n"
		(cd .. ; make -B dev)

		@echo "\n Running test simulations...\n"
		../LF_DEM -f -r 0.2r -n confs/D3N200VF0.45Bidi1.4_0.5Cubic_1.bin repulsive.txt > repulsive_rate.log 2>&1
		../LF_DEM -f -s 5r 	 confs/D3N200VF0.45Bidi1.4_0.5Cubic_1.dat repulsive_stress.txt > repulsive_stress.log 2>&1
		@echo "\n Diffing...\n"
		@for f in data_*.dat par_*dat int_*dat st_*dat; do \
			echo "diff $$f reference/$$f"; \
			diff  -I '^# LF_DEM.*' $$f reference/$$f; \
		done
