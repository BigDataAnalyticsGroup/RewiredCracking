#!env python3

#=======================================================================================================================
#
#   Generate the tables for the AVX permutations.
#
#   This script generates C arrays for permuting elements in AVX registers.  The index into the array is the mask,
#   obtained by an AVX comparison operation, and converted to an integer via `movemask`.
#   Since SSE/AVX only support `cmpeq` and `cmpgt` for comparison, we evaluate `pivot > key`.  Hence, a `1` bit in the
#   mask means the element belongs to the lower partition.
#   The mask is processed from LSB to MSB.
#
#=======================================================================================================================


from array import array


def make_packable(permutation):
    return permutation[0:2] + permutation[4:6] + permutation[2:4] + permutation[6:8]


def get_lo_permutation_from_mask(lanes, mask):
    permutation = array('i', [0] * lanes)
    lo = bin(mask).count('1')
    hi = 0
    for i in range(0, lanes):
        if mask & 0x1:
            permutation[hi] = i
            hi += 1
        else:
            permutation[lo] = i
            lo += 1
        mask >>= 1
    return permutation.tolist()


def get_hi_permutation_from_mask(lanes, mask):
    permutation = array('i', [0] * lanes)
    lo = 0
    hi = lanes - bin(mask).count('1')
    for i in range(0, lanes):
        if mask & 0x1:
            permutation[hi] = i
            hi += 1
        else:
            permutation[lo] = i
            lo += 1
        mask >>= 1
    return permutation.tolist()


def gen_LUT(lanes, PACKABLE=False):
    size = int(2**lanes)
    element_size = int(256 / lanes)
    print('const uint8_t permutation{lanes}{packable}_table[{size}][8] = {{'\
            .format(size=size, lanes=lanes,packable='_packable' if PACKABLE else ''))
    if lanes == 8:
        for i in range(0, size):
            permutation = get_lo_permutation_from_mask(lanes, i)
            if PACKABLE:
                permutation = make_packable(permutation)
            print('    /* 0b{i:0{lanes}b} */ {{ {perm} }},'\
                    .format(i=i, lanes=lanes, perm=', '.join(str(x) for x in permutation)))
    elif lanes == 4:
        for i in range(0, size):
            permutation = get_lo_permutation_from_mask(lanes, i)
            permutation = list(map(lambda x : 2*x, permutation))
            permutation = [ permutation[0], permutation[0] + 1, \
                            permutation[1], permutation[1] + 1, \
                            permutation[2], permutation[2] + 1, \
                            permutation[3], permutation[3] + 1 ]
            if PACKABLE:
                permutation = make_packable(permutation)
            print('    /* 0b{i:0{lanes}b} */ {{ {perm} }},'\
                    .format(i=i, lanes=lanes, perm=', '.join(str(x) for x in permutation)))
    print('};')


if __name__ == '__main__':
    print('#pragma once\n')
    print('#include <cstdint>')
    print('#include <immintrin.h>\n')
    gen_LUT(4); # 4 * 64bit = 256bit
    print()
    gen_LUT(8); # 8 * 32bit = 256bit
    print()
    gen_LUT(4, True); # 4 * 64bit = 256bit
    print()
    gen_LUT(8, True); # 8 * 32bit = 256bit
