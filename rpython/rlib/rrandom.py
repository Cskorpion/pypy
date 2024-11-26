# this code is a version of the mersenne twister random number generator which
# is supposed to be used from RPython without the Python interpreter wrapping
# machinery etc.

# this is stolen from CPython's _randommodule.c

from rpython.rlib.rarithmetic import r_uint, intmask, r_uint32, r_uint64, r_int32

N = 624
M = 397
MATRIX_A = r_uint(0x9908b0df) # constant vector a
UPPER_MASK  = r_uint(0x80000000) # most significant w-r bits
LOWER_MASK = r_uint(0x7fffffff) # least significant r bits
MASK_32 = r_uint(0xffffffff)
TEMPERING_MASK_A = r_uint(0x9d2c5680)
TEMPERING_MASK_B = r_uint(0xefc60000)
MAGIC_CONSTANT_A = r_uint(1812433253)
MAGIC_CONSTANT_B = r_uint(19650218)
MAGIC_CONSTANT_C = r_uint(1664525)
MAGIC_CONSTANT_D = r_uint(1566083941)


class Random(object):
    def __init__(self, seed=r_uint(0)):
        self.state = [r_uint(0)] * N
        self.index = 0
        self.init_genrand(seed)

    def init_genrand(self, s):
        mt = self.state
        mt[0]= s & MASK_32
        for mti in range(1, N):
            mt[mti] = (MAGIC_CONSTANT_A *
                           (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + r_uint(mti))
            # See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
            # In the previous versions, MSBs of the seed affect
            # only MSBs of the array mt[].
            # for >32 bit machines 
            mt[mti] &= MASK_32
        self.index = N

    def init_by_array(self, init_key):
        key_length = len(init_key)
        mt = self.state
        self.init_genrand(MAGIC_CONSTANT_B)
        i = 1
        j = 0
        if N > key_length:
            max_k = N
        else:
            max_k = key_length
        for k in range(max_k, 0, -1):
            mt[i] = ((mt[i] ^
                         ((mt[i - 1] ^ (mt[i - 1] >> 30)) * MAGIC_CONSTANT_C))
                     + init_key[j] + r_uint(j)) # non linear
            mt[i] &= MASK_32 # for WORDSIZE > 32 machines
            i += 1
            j += 1
            if i >= N:
                mt[0] = mt[N - 1]
                i = 1
            if j >= key_length:
                j = 0
        for k in range(N - 1, 0, -1):
            mt[i] = ((mt[i] ^
                        ((mt[i - 1] ^ (mt[i - 1] >> 30)) * MAGIC_CONSTANT_D))
                     - i) # non linear
            mt[i] &= MASK_32 # for WORDSIZE > 32 machines
            i += 1
            if (i>=N):
                mt[0] = mt[N - 1]
                i = 1
        mt[0] = UPPER_MASK

    def _conditionally_apply(self, val, y):
        if y & r_uint(1):
            return val ^ MATRIX_A
        return val

    def genrand32(self):
        mt = self.state
        if self.index >= N:
            for kk in range(N - M):
                y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK)
                mt[kk] = self._conditionally_apply(mt[kk+M] ^ (y >> 1), y)
            for kk in range(N - M, N - 1):
                y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK)
                mt[kk] = self._conditionally_apply(mt[kk + (M - N)] ^ (y >> 1),
                                                   y)
            y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK)
            mt[N - 1] = self._conditionally_apply(mt[M - 1] ^ (y >> 1), y)
            self.index = 0
        y = mt[self.index]
        self.index += 1
        y ^= y >> 11
        y ^= (y << 7) & TEMPERING_MASK_A
        y ^= (y << 15) & TEMPERING_MASK_B
        y ^= (y >> 18)
        return y

    def random(self):
        a = intmask(self.genrand32() >> 5)
        b = intmask(self.genrand32() >> 6)
        return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0)

    def jumpahead(self, n):
        mt = self.state
        for i in range(N - 1, 1, -1):
            j = n % i
            mt[i], mt[j] = mt[j], mt[i]
        nonzero = False
        for i in range(1, N):
            mt[i] += r_uint(i + 1)
            mt[i] &= r_uint(0xffffffff)
            nonzero |= bool(mt[i])
        # Ensure the state is nonzero: in the unlikely event that mt[1] through
        # mt[N-1] are all zero, set the MSB of mt[0] (see issue #14591). In the
        # normal case, we fall back to the pre-issue 14591 behaviour for mt[0].
        if nonzero:
            mt[0] += r_uint(1)
            mt[0] &= r_uint(0xffffffff)
        else:
            mt[0] = r_uint(0x80000000)
        self.index = N

class Pcg32nogcRandom(object):
    """ https://www.pcg-random.org/download.html 
        Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
        * Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org """
    
    _alloc_flavor_ = 'raw'

    def __init__(self):
        self.state = r_uint64(0x853c49e6748fea9b)
        self.inc = r_uint64(0xda3e39cb94b95bdb)

    def seed_from_time(self):
        import time
        self.state = r_uint64(time.time() * 1000)    
    
    def randbelow(self, bound):
        assert bound > 0, "bound must be positive"

        bound = r_uint64(bound)

        assert bound < 2**32, "bound must be a 32 bit u_int"

        # To avoid bias, we need to make the range of the RNG a multiple of
        # bound, which we do by dropping output less than a threshold.
        # A naive scheme to calculate the threshold would be to do
        #
        #     uint32_t threshold = 0x100000000ull % bound;
        #
        # but 64-bit div/mod is slower than 32-bit div/mod (especially on
        # 32-bit platforms).  In essence, we do
        #
        #     uint32_t threshold = (0x100000000ull-bound) % bound;
        #
        # because this version will calculate the same modulus, but the LHS
        # value is less than 2^32.


        #threshold = r_uint32((2**32 - intmask(bound)) % intmask(bound))
        threshold = (r_uint64(2**32) - bound) % bound
    
        # Uniformity guarantees that this loop will terminate.  In practice, it
        # should usually terminate quickly; on average (assuming all bounds are
        # equally likely), 82.25% of the time, we can expect it to require just
        # one iteration.  In the worst case, someone passes a bound of 2^31 + 1
        # (i.e., 2147483649), which invalidates almost 50% of the range.  In 
        # practice, bounds are typically small and only a tiny amount of the range
        # is eliminated.

        while True:
            r = r_uint64(self.genrand32())
            if r >= threshold:
                return r_uint32(r % bound)

    def genrand32(self):
        oldstate = self.state
        # Advance internal state
        self.state = oldstate * r_uint64(6364136223846793005) + (self.inc | 1)
        # Calculate output function (XSH RR), uses old state for max ILP
        xorshifted = ((oldstate >> 18) ^ oldstate) >> 27
        rot = oldstate >> 59
        return r_uint32((xorshifted >> rot) | (xorshifted << ((-rot) & 31)))
