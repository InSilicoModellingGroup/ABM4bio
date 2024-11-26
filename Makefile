
ABM4bio = $(shell pwd)
BIODYNAMO_VERSION = 25dc979062d7b4367537c821e80b95e5727f668a
BIODYNAMO_FOLDER  = biodynamo-v1.05.120
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
	make clean; make fresh; \
	make blood_vessel_CTC; \
	make cancer_angiogenesis; \
	make cancer_radiation; \
	make HeLa_cells; \
	make neuron_astrocyte; \
	make neuron_bipolar; \
	make obstacle_demo_1; \
	make obstacle_demo_2; \
	make tumour_spheroid; \
	make vasculogenesis; \
	make wound_assay; \
	echo "\n\n ***************************"; \
	echo     " *** All tests completed ***"; \
	echo     " ***************************\n"
