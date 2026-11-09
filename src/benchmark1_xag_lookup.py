"""Generic reversible XOR/AND-graph (ANF) table lookup in native kickmix.

For <=9 address bits, compute the Boolean algebraic-normal-form (ANF)
of a packed classical table. Traverse monomials in depth-first order,
sharing prefix ANDs and measurement-uncomputing them in reverse order.
Every phase correction is handled by kickmix.init_and / del_and.

This is a research primitive, NOT the optimized 172-AND lookup reported
by other competitors. A full low-workspace Horner phase oracle is still
needed to make 30x viable.

The lookup XORs table[address] into target; calling it a second time
restores target without a separate high-workspace del_lookup.
"""
from __future__ import annotations


def anf_coefficients(table, address_bits: int) -> tuple[int, ...]:
    """Packed output-bit coefficients over GF(2); no numerical approximation."""
    if not 1 <= address_bits <= 9:
        raise ValueError("Supports 1..9 address bits")
    if len(table) != 1 << address_bits:
        raise ValueError("Table length must be exactly 2**address_bits")
    if any(v < 0 or not isinstance(v, int) for v in table):
        raise ValueError("Table entries must be nonnegative Python integers")
    result = list(table)
    for j in range(address_bits):
        for mask in range(1 << address_bits):
            if mask & (1 << j):
                result[mask] ^= result[mask ^ (1 << j)]
    return tuple(result)


def apply_anf_lookup(builder, address, target, table):
    """Reversible XOR lookup using depth-first shared AND monomials.

    Measurements occur only in kickmix.del_and, which also emits its
    required measurement-conditioned Clifford correction.
    """
    n = len(address)
    if not 1 <= n <= 9:
        raise ValueError("ANF lookup supports 1..9 address qubits")
    if len(table) != 1 << n:
        raise ValueError("Table size mismatch")
    if any(int(v).bit_length() > len(target) for v in table):
        raise ValueError("Target too narrow for table values")
    coefficients = anf_coefficients([int(v) for v in table], n)

    def xor_mask(control, word):
        k = 0
        while word:
            if word & 1:
                if control is None:
                    builder.x(target[k])
                else:
                    builder.cx(control, target[k])
            k += 1
            word >>= 1

    xor_mask(None, coefficients[0])

    def visit(start, prefix, mask):
        for j in range(start, n):
            child_mask = mask | (1 << j)
            if prefix is None:
                term = address[j]
                scratch = False
            else:
                term = builder.alloc_qubits(1)[0]
                builder.init_and(prefix, address[j], term)
                scratch = True
            xor_mask(term, coefficients[child_mask])
            visit(j + 1, term, child_mask)
            if scratch:
                builder.del_and(prefix, address[j], term)
                builder.free(term)

    visit(0, None, 0)


def comparison_circuits(table, address_bits, output_bits):
    """Return native generic-QROM and ANF-XAG *compute/erase* circuits."""
    import kickmix as km
    def circuit(method):
        b = km.CircuitBuilder()
        address = b.create_quantum_register(address_bits, "address")
        target = b.create_quantum_register(output_bits, "output")
        if method == "anf":
            apply_anf_lookup(b, address, target, table)
            apply_anf_lookup(b, address, target, table)
        elif method == "kickmix":
            b.init_lookup(table=table, address=address, target=target)
            b.del_lookup(table=table, address=address, target=target)
        else:
            raise ValueError("bad method")
        return b.finish_circuit()
    return circuit("kickmix"), circuit("anf")
