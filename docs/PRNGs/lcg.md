# Linear Congruential Generator (LCG)

!!! note "See Also"
    * [Wikipedia Article](https://en.wikipedia.org/wiki/Linear_congruential_generator)

A linear congruential generator (LCG, sometimes linear congruential random number generator, LCRNG) is a pseudo-random number generator defined by the recurrence relation:

$$
X_{n+1} = (aX_n + c) \bmod m
$$

for a modulus $m$, multiplier $a$, increment $c$, and internal state $X_n$.

## State Transition Function

The state transition function can be implemented as:

```cpp
// 32-bit LCG
uint32_t next_state(uint32_t state) {
    return (MULTIPLIER * state + INCREMENT) % MODULUS;
}
// 64-bit LCG
uint64_t next_state(uint64_t state) {
    return (MULTIPLIER * state + INCREMENT) % MODULUS;
}
```

or, in the case where the modulus $m = 2^{32}, 2^{64}$ the modulo operation becomes redundant due to wrapping addition:

```cpp
// 32-bit LCG
uint32_t next_state(uint32_t state) {
    return MULTIPLIER * state + INCREMENT;
}
// 64-bit LCG
uint64_t next_state(uint64_t state) {
    return MULTIPLIER * state + INCREMENT;
}
```

## Output Functions

Generally, the output function first truncates the state to just the higher-order bits (most commonly the top half of the state). This is because the higher bits have a longer period (see: [Distance Function](#distance-function)). Some exceptions are found in DPPt/HGSS/BDSP's Group Seed systems.

```cpp
// 32-bit LCG
uint16_t next_uint(uint32_t &state) {
    state = next_state(state);
    return state >> 16;
}
// 64-bit LCG
uint32_t next_uint(uint64_t &state) {
    state = next_state(state);
    return state >> 32;
}
```

### Modulo

The modulo output function simply reduces ``next_uint`` $\bmod n$ where $n$ is the range of the output. This is used throughout RSE, FRLG, Colosseum, XD, HGSS, and most instances of LCGs in the series.

```cpp
// 32-bit LCG
uint16_t next_rand(uint32_t &state, uint16_t n) {
    return next_uint(state) % n;
}
uint16_t next_rand(uint32_t &state, uint16_t min, uint16_t max) {
    return min + next_rand(state, max - min + 1);
}
```

This distribution is slightly uneven as it results in outputs $\in [0, \text{UINT_MAX} \bmod n]$ being weighted slightly higher than the rest of the potential outputs. For example, a 32-bit LCG with $n = 100$ (assuming uniform randomness), results in outputs $\in [0, 35]$ having a weight of 656 and outputs $\in [36, 99]$ having weight 655.

### Reciprocal Division

The reciprocal division output function takes the floor division of ``next_uint`` by a value $p = \lfloor \text{UINT_MAX} \div n \rfloor + 1$ where $n$ is the range of the output. This has the effect of partitioning the range of ``next_uint`` into $n - 1$ contiguous buckets of size $p$ and a final bucket of size $\text{UINT_MAX} \bmod p$. This is primarily seen in DPPt's main RNG.

```cpp
// 32-bit LCG
uint16_t next_rand(uint32_t &state, uint16_t n) {
    return next_uint(state) / ((0xFFFF / n) + 1);
}
uint16_t next_rand(uint32_t &state, uint16_t min, uint16_t max) {
    return min + next_rand(state, max - min + 1);
}
```

This distribution is slightly uneven as the last bucket (and thus the weight of output $n-1$) is smaller than the rest. For example, a 32-bit LCG with $n = 100$ (assuming uniform randomness), results in outputs $\in [0, 98]$ having a weight of 656 and output $99$ having weight 592. This issue could have been easily rectified by using $p = \lfloor \text{UINT_MAX} \div n \rfloor$ and rejection sampling out of range values (see: [Paper Mario](https://github.com/pmret/papermario/blob/8a0ea06aa996ecb1fd0545c272039d223b1a7ee9/src/43F0.c#L503-L524)).

### Multiply-Shift

The multiply-shift output function multiplies ``next_uint`` by $n$ (the range of the output) and right shifts by the size in bits of ``next_uint``. As right shifting is equivalent to floor division, this approximates multiplying by $n \div (\text{UINT_MAX}+1)$ or dividing by $(\text{UINT_MAX}+1) \div n$. This has the effect of partitioning the range of ``next_uint`` into $n$ contiguous buckets of roughly equal size $p \approx (\text{UINT_MAX}+1) \div n$. This is used in BW/B2W2's 64-bit LCG.

```cpp
// 64-bit LCG
uint32_t next_rand(uint64_t &state, uint32_t n) {
    return ((uint64_t)next_uint(state) * (uint64_t)n) >> 32;
}
uint32_t next_rand(uint64_t &state, uint32_t min, uint32_t max) {
    return min + next_rand(state, max - min + 1);
}
```

This distribution is slightly uneven as it results in $\text{UINT_MAX + 1} \bmod n$ of the outputs being weighted slightly higher than the rest (these are evenly distributed throughout the range). For example, a 32-bit LCG with $n = 100$ (assuming uniform randomness), results in 36 outputs having a weight of 656 and the other 64 having weight 655. This is a similar distribution to [modulo](#modulo), but is less lopsided and does not require expensive division.

### Random Float

To generate a uniform float $\in [0,1)$, ``next_uint`` is simply divided by ``UINT_MAX + 1``.

```cpp
// 32-bit LCG
float next_float(uint32_t &state) {
    return next_uint(state) / 65536.0f;
}
float next_float(uint32_t &state, float min, float max) {
    return min + next_float(state) * (max - min);
}
```

## Stepping Backwards

Stepping backwards requires solving the recurrence relation for $X_n$ as follows:

$$
\begin{aligned}
X_{n+1} &= (aX_n + c) \bmod m \\
X_{n+1} &\equiv aX_n + c \pmod m \\
X_{n+1} - c &\equiv aX_n \pmod m \\
X_n &\equiv a^{-1}(X_{n+1} - c) \pmod m
\end{aligned}
$$

This requires that the modular multiplicative inverse $a^{-1}$ exists (it does, as $a$ is chosen to be coprime to $m$). Notably, this relation can itself be simplified into an LCG:

$$
\begin{aligned}
X_n &\equiv a^{-1}(X_{n+1} - c) \pmod m \\
&\equiv a^{-1}X_{n+1} - a^{-1}c \pmod m \\
c^{-1} &= - a^{-1}c \bmod m \\
X_n &\equiv a^{-1}X_{n+1} + c^{-1} \pmod m \\
X_n &= (a^{-1}X_{n+1} + c^{-1}) \bmod m
\end{aligned}
$$

$a^{-1}$ and therefore $c^{-1}$ can be found via the extended Euclidean algorithm as follows:

```python
def extended_euclid(a, b): 
    if a == 0:
        return 0, 1

    x1, y1 = extended_euclid(b % a, a)
    x = y1 - (b // a) * x1
    y = x1
    
    return x, y

reverse_mult, _ = extended_euclid(mult, mod)
reverse_add = (reverse_mult * -add) % mod
```

## Jumping

Efficiently jumping an LCG requires finding some easy to calculate $f(X_0, n) = X_n$ for arbitrary $n$. The naive implementation is to simply apply the state transition function $n$ times which is $O(n)$ time complexity. A much more efficient constant-time approach can be found as follows:

First, notice that the recurrence relation is easily composable and reduced to another LCG:

$$
\begin{align}
X_{n+1} &= (aX_n + c) \bmod m \\
X_{n+2} &= (aX_{n+1} + c) \bmod m \\
&= (a(aX_n + c) + c) \bmod m \\
&= (a^2X_n + ac + c) \bmod m \\
c_2 &= (ac + c) \bmod m \\
&= (a^2X_n + c_2) \bmod m \\
&\implies \\
X_{n} &= (a^nX_0 + c_n) \bmod m
\end{align}
$$

so $f(X_0, n)$ is trivial to calculate if $a^n$ and $c_n$ are known.

Then, note that $f$ is composable such that $f(f(X_0, n_1), n_2) = f(X_0, n_1 + n_2)$. This means that if $n$ is decomposable into a sum $\sum_{i=1}^{N} n_i$ where each $a^i$ and $c_i$ are known, $f(X_0, n)$ is computable via $N < n$ LCG state transitions with $O(N)$ time complexity.

Given that the period of the LCG is $m$, any $n \geq m$ can simply be reduced $\bmod m$ and only $n \in [0, m)$ must be considered.

An easy to calculate and effective decomposition of $n$ is given by the binary representation of $n = \sum_{i=1}2^{i-1}b_i$. In the case where $m$ is a power of 2, the number of terms and thus the $N$ that the time complexity is dependent on is the exponent of $m$ and is constant. This means a 32-bit LCG can be jumped by an arbitrary $n$ in only 32 state transitions (and a 64-bit LCG in 64) as follows:

```cpp
// precomputed mults and adds
constexpr uint32_t MULTIPLIERS[32] = ...; // M[i] = a^{i}
constexpr uint32_t INCREMENTS[32] = ...; // I[i] = c_{i}

// 32-bit LCG
void jump(uint32_t &state, uint32_t n) {
    // instead of checking i < 32, this could check n != 0
    // to exit early when the most significant set bit is reached
    for (int i = 0; i < 32; i++) {
        // extract the i-th bit
        if (n & 1) {
            state = MULTIPLIERS[i] * state + INCREMENTS[i];
        }
        n >>= 1;
    }
}
```

Another method worth noting is jumping by representing the state transition function as multiplication by a particular matrix. The LCG recurrence relation is equivalent to the following:

$$
\begin{align}
\begin{bmatrix}
X_{n+1} \\ 1
\end{bmatrix} &= \begin{bmatrix}
a & c \\
0 & 1 \\
\end{bmatrix}\begin{bmatrix}
X_n \\ 1
\end{bmatrix} \bmod m \\
&\implies \\
\begin{bmatrix}
X_{n} \\ 1
\end{bmatrix} &= \begin{bmatrix}
a & c \\
0 & 1 \\
\end{bmatrix}^n\begin{bmatrix}
X_0 \\ 1
\end{bmatrix} \bmod m
\end{align}
$$

which allows efficient jumping if the matrix can be efficiently calculated via a method like exponentiation by squares.

## Distance Function

The problem of finding the distance between two LCG states can be thought of as the inverse of the problem of [jumping](#jumping). It requires some easy to calculate $d(X_n, X_m) = m - n$ when only $X_n$ and $X_m$ are known. Given the matrix representation of the state transition function above, it is clear that this is a particular case of the discrete logarithm problem. In the case of an LCG with a prime power modulus, $d$ can be calculated in constant-time as follows:

First, starting with the trivial fact that $(N \bmod p^u) \bmod p^v = N \bmod p^v \quad \forall v < u$, notice that reducing the state of the LCG by a smaller prime power produces a smaller LCG:

$$
\begin{align}
X_{n+1} &= (aX_n + c) \bmod p^u \\
X_{n+1} \bmod p^v &= ((aX_n + c) \bmod p^u) \bmod p^v \\
&= (aX_n + c) \bmod p^v \\
&= ((a\bmod p^v)X_n + (c\bmod p^v)) \bmod p^v \\
(X_v)_{n+1} &= (a_v(X_v)_n + c_v) \bmod p^v
\end{align}
$$

Given that $a$ and $c$ satisfy the properties for a full period LCG, $a_v$ and $c_v$ do as well (for the smaller modulus). This means that the sequence given by reducing the state of the LCG by $p^v$ has a period of $p^v$ which implies the distance function applied to this sequence $d_v = d((X_v)_n, (X_v)_m) = d \bmod p^v$.

Similar logic can be applied to other values of $p$, but for the case of $p = 2$, $d$ can be computed bit by bit as follows:

Starting with the least significant bit, the sequence $X_1$ has period $2^1 = 2$ and the distance between $(X_1)_n$ and $(X_1)_m$ ($d_1$) must be either 0 or 1. This is easy to compute as a distance of 0 implies equality and a distance of 1 implies inequality.

Next, the sequence $X_2$ has period $2^2 = 4$ and the distance between $(X_2)_n$ and $(X_2)_m$ ($d_2$) must be $\in [0, 3]$. From the prior point, $d_2 \bmod 2 = d_1 \implies d_2 = 2q + d_1$ where $q$ is either 0 or 1. Apply the state transition function to $(X_2)_n$ $d_1$ times to obtain a state $(X_2)_{n+d_1}$ that must be a distance $2q$ from $(X_2)_m$. $q$ is then trivial to compute as a distance of 0 implies equality and a distance of 2 implies inequality.

With $d_2$ known, $d_3$ and all subsequent partial distances can be computed via the same procedure. This process can be concluded when the state $(X_i)_{n+d_{i-1}} = X_m$ as this implies $d = d_{i-1}$ and the desired distance is acquired. In the case of a 32-bit LCG, the distance between two states $X_n$ and $X_m$ can be computed as follows:

```cpp

// 32-bit LCG
uint32_t distance(uint32_t state_n, uint32_t state_m) {
    uint32_t d_i = 0;
    uint32_t i = 1;
    // after each iteration, state_n becomes X_{n+d_{i-1}}
    while (state_n != state_m) {
        uint32_t state_i_n = i < 32 ? state_n % (1 << i) : state_n;
        uint32_t state_i_m = i < 32 ? state_m % (1 << i) : state_m;
        // q == 1
        if (state_i_n != state_i_m) {
            d_i += 1 << (i - 1);
            jump(state_n, 1 << (i - 1));
        }
        i++;
    }
    return d_i;
}

// simplified
uint32_t distance(uint32_t state_n, uint32_t state_m) {
    uint32_t d_i = 0;
    uint32_t mask = 1;
    while (state_n != state_m) {
        // because each jump won't modify the lower-order bits, we only need to check the next bit
        if ((state_n ^ state_m) & mask) {
            jump(state_n, mask);
            d_i += mask;
        }
        mask <<= 1;
    }
    return d_i;
}

// inlining jump logic
uint32_t distance(uint32_t state_n, uint32_t state_m) {
    uint32_t d_i = 0;
    uint32_t mask = 1;
    uint32_t i = 0;
    while (state_n != state_m) {
        if ((state_n ^ state_m) & mask) {
            // we only ever jump a power of 2, which we have precomputed mult/inc for
            state_n = MULTIPLIERS[i] * state_n + INCREMENTS[i];
            d_i += mask;
        }
        mask <<= 1;
        i++;
    }
    return d_i;
}
```

## Variants

| Name | Multiplier | Increment | Size | Usage |
| ---- | ---------- | --------- | ---- | ----- |
| PokeRNG | 0x41C64E6D | 0x6073 | 32-bit | RSE/FRLG/DPPt/HGSS main RNG |
| PokeRNGR | 0xEEB9EB65 | 0xA3561A1 | 32-bit | ^ stepped backwards |
| ARNG | 0x6C078965 | 0x1 | 32-bit | DPPt/HGSS group seed & pid rerolling |
| ARNGR | 0x9638806D | 0x69C77F93 | 32-bit | ^ stepped backwards |
| XDRNG | 0x343FD | 0x269EC3 | 32-bit | Colosseum/XD/Channel/Ageto main RNG |
| XDRNGR | 0xB9B33155 | 0xA170F641 | 32-bit | ^ stepped backwards |
| MRNG | 0x41C64E6D | 0x3039 | 32-bit | DPPt lotto, Rumble main RNG |
| MRNGR | 0xEEB9EB65 | 0xFC77A683 | 32-bit | ^ stepped backwards |
| BWRNG | 0x5d588b656c078965 | 0x269ec3 | 64-bit | BW/B2W2 main RNG |
| BWRNGR | 0xdedcedae9638806d | 0x9b1ae6e9a384e6f9 | 64-bit | ^ stepped backwards |

## References

* [Wikipedia](https://en.wikipedia.org/wiki/Linear_congruential_generator)
