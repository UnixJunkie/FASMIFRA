#!/usr/bin/python3
#
# encode '^bitstring:string val:float$' lines to a sparse
# .AP fingerprint format

import sys

input_fn = sys.argv[1]

with open(input_fn, 'w') as input:
    for i, line in enumerate(input.readlines()):
        stripped = line.strip()
        bitstring, target_val = stripped.split(' ')
        bits = list(bitstring)
        num_bits = len(bits)
        assert(len(bits) == 1024)
        print('mol_%d,%s,[' % (i, target_val, file=sys.stdout, end=''))
        # use 1-based indices; because of liblinear down the line
        for i in range(0, num_bits):
            b = bits[i]
            assert(b == '1' or b == '0')
            if b == 1:
                if i > 0:
                    print(';%d:1' % (i+1), file=sys.stdout, end='')
                else: # i == 0
                    print('%d:1' % (i+1), file=sys.stdout, end='')
        print(']' % (i, target_val, file=sys.stdout)) # EOL
