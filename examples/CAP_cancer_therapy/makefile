build:
	cd ../.. && make fresh -j10

# run:
# 	../../build/./ABM4bio input.csv 1234

# tumor:
# 	../../build/./ABM4bio input_copy.csv 1234

# invivo:
# 	../../build/./ABM4bio input_inVivo.csv 1234

invitro:
	rm -rf results_inVitro
	OMP_NUM_THREADS=80 ../../build/./ABM4bio input_inVitro.csv 1234

# analyze_invivo:
# 	python3 analyze_invivo_tumor.py

clean:
	rm -rf results results_inVitro /home/aiwsif/Desktop/ABM4bio/build
