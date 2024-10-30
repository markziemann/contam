#!/bin/bash

apptainer build --disable-cache contam.sif contam.def \
  && scp contam.sif mdz@ziemann-lab.net:~/public/contam/

