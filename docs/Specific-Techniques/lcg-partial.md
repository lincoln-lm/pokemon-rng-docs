# LCG State Recovery From Two Consecutive Partial Outputs

Given two consecutive [PokeRNG](../../PRNGs/lcg/#variants) outputs $O_n$ and $O_{n+1}$, the goal is to find all internal states $X_{n-1}$ that produce them. This is particularly useful when trying to find all RNG states that produce a specific PID or set of IVs generated via Method 1/2/4 in Generations 3 and 4, but can be applied whenever consecutive RNG outputs are known.

Trivially, $X_n$ must be of the form $O_n \cdot 2^{16} + L_n$ and $X_{n+1}$ of the form $O_{n+1} \cdot 2^{16} + L_{n+1}$ where $O_n, O_{n+1}, L_n, L_{n+1} < 2^{16}$.

$L_{n+1}$ can be directly computed from $L_{n}$ due to the property that [the lower order bits of an LCG constitute their own smaller LCG](../../PRNGs/lcg/#distance-function).

Restated, given $O_n$, $O_{n+1}$, find all $L_n$ such that:

$$
O_{n+1} \cdot 2^{16} + L_{n+1} \equiv a(O_n \cdot 2^{16} + L_n) + c \pmod m
$$

for some $L_{n+1}$.

For simplicity, call the ordered sequence of solutions $S_i$, where each $S_i$ is a valid value for $L_n$ that satisfies the congruence.

Manipulate the congruence to obtain:

$$
\begin{align}
O_{n+1} \cdot 2^{16} + L_{n+1} &\equiv a(O_n \cdot 2^{16} + S_i) + c \pmod m \\
O_{n+1} \cdot 2^{16} + L_{n+1} &\equiv aO_n \cdot 2^{16} + aS_i + c \pmod m \\
(O_{n+1} - aO_n) \cdot 2^{16} + L_{n+1} &\equiv aS_i + c \pmod m \\
\end{align}
$$

Take any two solutions $S_0$ and $S_1$ such that $S_0 < S_1$ that satisfy the recurrence for $L$ and $L'$ respectively:

$$
\begin{align}
(O_{n+1} - aO_n) \cdot 2^{16} + L &\equiv aS_0 + c \pmod m \\
(O_{n+1} - aO_n) \cdot 2^{16} + L' &\equiv aS_1 + c \pmod m \\
\implies (L' - L) &\equiv a(S_1 - S_0) \pmod m
\end{align}
$$

Given that $L, L' < 2^{16}$, their difference lies in $(-2^{16}, 2^{16})$ and the difference between any two solutions is constrained as follows:

$$
\begin{align}
a(S_1 - S_0) \bmod m &< 2^{16} \quad \text{or} \\
a(S_1 - S_0) \bmod m &> -2^{16} \bmod m
\end{align}
$$

All valid differences can be found as follows:

```python
m = 2**32
a = 0x41C64E6D
for difference in range(1, 2**16):
    result = (a * difference) % m
    if result < 2**16 or result > m - 2**16:
        print(hex(difference))
```

which produces ``0x67D3`` and ``0xCFA6`` (``0x67D3 * 2``).

This result implies that, given the smallest potential solution $S_0$, $S_1 = S_0 + \text{0x67D3}$, $S_2 = S_0 + 2\cdot\text{0x67D3}$. Importantly, these are not all guaranteed to be actual solutions, so they must be manually tested.

To find a potential solution to work from, start with this form of the congruence:

$$
(O_{n+1} - aO_n) \cdot 2^{16} + L_{n+1} \equiv aS_i + c \pmod{2^{32}}
$$

Written as an equation:

$$
\begin{align}
(O_{n+1} - aO_n) \cdot 2^{16} + L_{n+1} + 2^{32}k &= aS_i + c \\
(O_{n+1} - aO_n + 2^{16}k) \cdot 2^{16} + L_{n+1} &= aS_i + c
\end{align}
$$

Multiplying both sides by $q = \text{0x67D3}$ allows the earlier property (``0x67D3*0x41C64E6D = 0x1aad00007ED7``) to be applied for the substitution $qa = 2^{32}p + r$

$$
\begin{align}
q(O_{n+1} - aO_{n} + 2^{16}k) \cdot 2^{16} + qL_{n+1} &= qaS_i + qc \\
q(O_{n+1} - aO_{n} + 2^{16}k) \cdot 2^{16} + qL_{n+1} &= (p\cdot2^{32} + r)S_i + qc \\
q(O_{n+1} - aO_{n} + 2^{16}k) + \frac{qL_{n+1}}{2^{16}} &= p\cdot2^{16}\cdot S_i + \frac{rS_i+qc}{2^{16}} \\
q(O_{n+1} - aO_{n} + 2^{16}k) + \frac{qL_{n+1} - rS_i - qc}{2^{16}} &= p\cdot2^{16}\cdot S_i
\end{align}
$$

Calling the fractional term $F$ and isolating $S_i$:

$$
\begin{align}
qk + \lfloor\frac{q(O_{n+1} - aO_{n}) + F}{2^{16}}\rfloor = p\cdot S_i \\
S_i \equiv \lfloor\frac{q(O_{n+1} - aO_{n}) + F}{2^{16}}\rfloor \cdot p^{-1} \pmod q
\end{align}
$$

Since this congruence is already $\pmod q$, the smallest potential solution $S_0$ is just the principal solution given by:

$$
S_0 = \lfloor\frac{q(O_{n+1} - aO_{n}) + F}{2^{16}}\rfloor \cdot p^{-1} \bmod q
$$

$F$ will vary depending on the value of $S_i$ and $L_{n+1}$ (though since these directly determine each other it can be thought of as only depending on one), but any value of $F$ less than $2^{16}$ above the true value will produce the correct $S_i$ due to the floored division.

The range of values for $F$ can be tested as follows:

```python
q = 0x67D3
r = 0x7ED7
a = 0x41C64E6D
c = 0x6073

min_F = float("inf")
max_F = float("-inf")
for S_i in range(0x10000):
    L = (a * S_i + c) & 0xFFFF
    F_numerator = q * L - r * S_i - q * c

    assert (F_numerator & 0xFFFF) == 0
    F = F_numerator >> 16
    if F < min_F:
        min_F = F
    if F > max_F:
        max_F = F

print(hex(min_F), hex(max_F))
```

which reveals $F \in [\text{-0xA584}, \text{0x4034}]$. Using a value of ``0x4034`` for $F$ will be between ``0x4034 - 0x4034 = 0x0`` and ``0x4034 - (-0xA584) = 0xE5B8`` above the true value which is always in the margin of error that allows for the computation of $S_0$. ``0x4034`` is the minimum such value, and anything up to and including ``0x5A7B``(``0xFFFF - 0xA584``) will work.

Computing $S_0$ is easily implemented as follows:

```python
q = 0x67D3
r = 0x7ED7
p = 0x1AAD
p_inv = pow(p, -1, q) # 0xD3E
a = 0x41C64E6D
c = 0x6073
F = 0x4034 # [0x4034, 0x5A7B]

def principal_candidate_solution(output_0, output_1):
    # floored division implemented as right shift
    return (((q * (output_1 - a * output_0) + F) >> 16) * p_inv) % q
```

and the full computation of all valid states that produce two consecutive outputs is implemented as:

```python
m = 2**32
a_inv = pow(a, -1, m) # 0xEEB9EB65

def test_solution(state, output_1):
    # output_0 need not be tested as it should just be the upper half of ``state``
    return (((state * a + c) >> 16) & 0xFFFF) == output_1

def solve(output_0, output_1):
    # at most 3
    solutions = []
    S_0 = principal_candidate_solution(output_0, output_1)
    for i in range(3):
        S = S_0 + q * i
        X_n = S | (output_0 << 16)
        if test_solution(X_n, output_1):
            # step backwards to get the state X_{n-1} that produces the outputs sequentially
            solutions.append(((X_n - c) * a_inv) % m)
    return solutions

```

## Attributions

- This technique can be attributed to user Parzival/StarfBerry who first shared it in a ([now deleted](https://github.com/StarfBerry/poke-scripts/blob/master/RNG/LCG_Reversal.py)) GitHub repository
