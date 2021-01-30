#!/usr/bin/env bash
python3 setup.py install --record files.txt
cd test
python3 run_all_unit_tests.py
cd ..
