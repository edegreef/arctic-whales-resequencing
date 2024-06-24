#!/bin/bash

sample=88_Pang_S51.rightwhale

# Bowhead
# 88_Pang
# BM_RB_2005_001

# Narwhal
# ARCR_07_1065 - BI
# ARGF_07_1127 - CHA
# ARRB_99_1026 - NHB

/home/degreefe/programs/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o $sample.N25.t15.r5.psmc $sample.psmcfa
