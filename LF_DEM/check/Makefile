test:
		@echo "\n Building...\n"
		(cd .. ; make -B dev)

		@echo "\n Running test simulations...\n"
		../LF_DEM -p 10 -k knkt_mu1.dat -n confs/D3N500VF0.45Bidi1.4_0.5Cubic_1_colloid_determine_knkt_0.04_p10.bin colloid_determine_knkt_0.04.txt > knkt_det.log 2>&1
		../LF_DEM -r 0.2 -n confs/D3N500VF0.45Bidi1.4_0.5Cubic_1_repulsive_r0.2.bin -k knkt_mu1.dat repulsive.txt > repulsive_rate.log 2>&1
		../LF_DEM -s 5 -k knkt_mu1.dat -n confs/D3N500VF0.45Bidi1.4_0.5Cubic_1_repulsive_s5.bin repulsive.txt > repulsive_stress.log 2>&1
		../LF_DEM -p 100 -k knkt_mu1_peclet_ov0.02.txt -n confs/D3N500VF0.45Bidi1.4_0.5Cubic_1_colloid_p100.bin colloid.txt > brownian.log 2>&1
		../LF_DEM -p 0.1 -k knkt_mu1_peclet_ov0.02.txt -n confs/D3N500VF0.45Bidi1.4_0.5Cubic_1_colloid_p0.1.bin colloid.txt > brownian_lowPe.log 2>&1
		@echo "\n Diffing...\n"
		@for f in rheo*; do \
			echo "diff $$f reference/$$f"; \
			diff  -I '^# LF_DEM.*' $$f reference/$$f; \
		done

