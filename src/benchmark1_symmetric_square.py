"""Symmetry-aware reversible truncated square.

Instead of generating both ordered partial products x_i*x_j and x_j*x_i,
combine them into a single contribution at the next binary position.
The dropped-pair cutoff is preserved exactly, including omitted carries.
"""
def square_symmetric(builder, x, dest, drop, inverse=False):
    n = len(x)
    if not (0 <= drop < 2 * n) or len(dest) != 2 * n - drop:
        raise ValueError("Incorrect symmetric-square dimensions")
    contributions = []
    for i in range(n):
        if 2 * i >= drop:
            contributions.append((2 * i - drop, i, i))
        for j in range(i + 1, n):
            if i + j >= drop:
                contributions.append((i + j + 1 - drop, i, j))

    for pos, i, j in (reversed(contributions) if inverse else contributions):
        if i == j:
            flag = x[i]
            allocated = False
        else:
            flag = builder.alloc_qubits(1)[0]
            builder.init_and(x[i], x[j], flag)
            allocated = True
        if inverse:
            builder.isub(1, target=dest[pos:], control=flag)
        else:
            builder.iadd(1, target=dest[pos:], control=flag)
        if allocated:
            builder.del_and(x[i], x[j], flag)
            builder.free(flag)
