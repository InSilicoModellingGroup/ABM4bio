build:
	cd ../.. && make fresh -j10

invitro:
	rm -rf results_inVitro
	OMP_NUM_THREADS=80 ../../build/./ABM4bio input_inVitro.csv 1234

# analyze_invivo:
# 	python3 analyze_invivo_tumor.py

clean:
	rm -rf results results_inVitro /home/aiwsif/Desktop/ABM4bio/build
