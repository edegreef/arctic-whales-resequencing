#!/bin/bash

# bowhead
sample=88_Pang

# narwhal
#sample=ARRB_99_1026

/home/degreefe/programs/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o $sample.N25.t15.r5.psmc $sample.psmcfa
