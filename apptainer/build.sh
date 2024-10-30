#!/bin/bash

apptainer build contam.sif contam.def \
  && scp contam.sif mdz@ziemann-lab.net:~/public/contam/

