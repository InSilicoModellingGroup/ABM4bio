#!/bin/bash
mkdir -p /tmp/calib_bin
ln -sf /usr/bin/python3.10 /tmp/calib_bin/python3
export PATH="/tmp/calib_bin:$PATH"
cd /home/aiwsif/Desktop/ABM4bio/examples/CAP_cancer_therapy
make clean_calibration && make calibrate_control
