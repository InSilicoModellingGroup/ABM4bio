
ABM4bio = $(shell pwd)
#BIODYNAMO_VERSION = 25dc979062d7b4367537c821e80b95e5727f668a
#BIODYNAMO_FOLDER  = biodynamo-v1.05.120
BIODYNAMO_VERSION = a9d3c90e97164660d0ce567a357eb1cfe38035aa
BIODYNAMO_FOLDER  = biodynamo-v1.05.143
CC  = /usr/bin/gcc
CXX = /usr/bin/g++

all:
	@echo "\n *** Please provide a specific rule to execute *** \n"
whereami:
	@echo " ...over here:" $(ABM4bio)
setup_biodynamo:
	rm -Rf $(ABM4bio)/libs/biodynamo; \
	git clone https://github.com/BioDynaMo/biodynamo.git $(ABM4bio)/libs/biodynamo && \
	cd $(ABM4bio)/libs/biodynamo && \
	git checkout $(BIODYNAMO_VERSION)
install_biodynamo:
	cd $(ABM4bio)/libs/biodynamo && \
	./prerequisites.sh all && \
	rm -rf build; mkdir build; cd build; \
	cmake -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX) \
        -Dparaview=off -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$(ABM4bio)/libs .. && \
	make -j 4 && make install
biodynamo:
	make setup_biodynamo && make install_biodynamo
source_biodynamo:
	/bin/bash -c 'source $(ABM4bio)/libs/$(BIODYNAMO_FOLDER)/bin/thisbdm.sh'
clean:
	biodynamo clean && \
	rm -rf $(ABM4bio)/build
fresh:
	make clean; \
	mkdir $(ABM4bio)/build; cd $(ABM4bio)/build; \
	cmake -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX) .. && \
	make -j 4
%:
	echo "\n\n *****************"; \
	echo     " *** Running test: "$@; \
	echo     " ***************** "; \
	cd $(ABM4bio)/examples/$@; rm -rf results; make run
tests:
	make clean; \
	make --silent fresh                2>errors; \
	make --silent blood_vessel_CTC     2>errors; \
	make --silent cancer_angiogenesis  2>errors; \
	make --silent cancer_radiation     2>errors; \
	make --silent HeLa_cells           2>errors; \
	make --silent neuron_astrocyte     2>errors; \
	make --silent neuron_bipolar       2>errors; \
	make --silent obstacle_demo_1      2>errors; \
	make --silent obstacle_demo_2      2>errors; \
	make --silent tumour_spheroid      2>errors; \
	make --silent vasculogenesis       2>errors; \
	make --silent wound_assay          2>errors; \
	echo "\n\n ***************************"; \
	echo     " *** All tests completed ***"; \
	echo     " ***************************\n"
