build:
	cd ../.. && make fresh

run:
	../../build/./ABM4bio input.csv 1234

invivo:
	../../build/./ABM4bio input_inVivo.csv 1234

analyze_invivo:
	python3 analyze_invivo_tumor.py

clean:
	rm -rf results results_inVivo
