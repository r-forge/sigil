#!/bin/sh

cwb-scan-corpus -o bnc_adj_n.tsv BNC text_id+0 'hw+0=/([a-z]+-)*[a-z]+/' '?class+0=/ADJ/' 'hw+1=/([a-z]+-)*[a-z]+/' '?pos+1=/NN.*/'
