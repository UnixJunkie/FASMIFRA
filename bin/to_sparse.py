#!/usr/bin/python3
#
# encode '^bitstring:string val:float$' lines to a sparse
# .AP fingerprint format

import sys

input_fn = sys.argv[1]

with open(input_fn, 'r') as input:
    for i, line in enumerate(input.readlines()):
        stripped = line.strip()
        bitstring_target = stripped.split()
        bitstring = bitstring_target[0]
        target_val = bitstring_target[1]
        bits = list(bitstring)
        num_bits = len(bits)
        # assert(len(bits) == 1024) # current use case
        print('mol_%d,%s,[' % (i, target_val), file=sys.stdout, end='')
        # use 1-based indices; because of liblinear down the line
        start = True
        for i in range(0, num_bits):
            b = bits[i]
            assert(b == '1' or b == '0')
            if b == '1':
                if start:
                    print('%d:1' % (i+1), file=sys.stdout, end='')
                    start = False
                else:
                    print(';%d:1' % (i+1), file=sys.stdout, end='')
        print(']', file=sys.stdout) # EOL
