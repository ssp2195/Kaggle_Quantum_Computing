"""Low-workspace, measurement-corrected binary exponent normalization.

The older barrel normalizer retains ten prefix-AND flags until the end
of the oracle. Here only the four-qubit exponent register survives the
barrel shift:

  compute prefix chain on ORIGINAL x
  increment binary exponent for each true prefix
  measure-and-correct erase prefix AND temporaries while x remains original
  rotate x by exponent

After the phase:
  inverse rotate x using exponent
  recompute prefix chain on recovered ORIGINAL x
  subtract the same exponent contributions
  measure-and-correct erase the chain and free exponent

All nonlinear flag measurements have available live AND controls.
This is NOT free computation: each repeated prefix requires an extra
forward Toffoli, while del_and correction itself requires only Clifford
operations. The trade should lower peak qubits but add some gates.
"""
from __future__ import annotations
from .priority1_normalize import _zero_pair, _erase_zero_pair, _lanes


def _compute_prefixes(b,x):
    prefixes=[]
    for pair_index in range(10,0,-1):
        z=_zero_pair(b,x,pair_index)
        if prefixes:
            prefix=b.alloc_qubits(1)[0]
            b.init_and(prefixes[-1],z,prefix)
            _erase_zero_pair(b,x,pair_index,z)
        else:
            prefix=z
        prefixes.append(prefix)
    return prefixes


def _erase_prefixes(b,x,prefixes):
    if len(prefixes)!=10:
        raise ValueError("Expected ten original prefix flags")
    for idx in range(len(prefixes)-1,-1,-1):
        prefix=prefixes[idx]
        pair_index=10-idx
        if idx:
            z=_zero_pair(b,x,pair_index)
            b.del_and(prefixes[idx-1],z,prefix)
            b.free(prefix)
            _erase_zero_pair(b,x,pair_index,z)
        else:
            _erase_zero_pair(b,x,pair_index,prefix)


def normalize_compact(builder,x):
    """Interface-compatible return (None, binary_exponent_register)."""
    if len(x)!=22:
        raise ValueError("Expected 22 input qubits")
    b=builder
    shift=b.alloc_qubits(4)
    b.reset(shift)
    prefixes=_compute_prefixes(b,x)
    for prefix in prefixes:
        b.iadd(1,target=shift,control=prefix)
    _erase_prefixes(b,x,prefixes)
    even,odd=_lanes(x)
    b.left_rotate(even,shift=shift)
    b.left_rotate(odd,shift=shift)
    return None,shift


def denormalize_compact(builder,x,_unused,shift):
    """Unshift before recomputing prefixes to ensure causal correction."""
    b=builder
    even,odd=_lanes(x)
    b.right_rotate(even,shift=shift)
    b.right_rotate(odd,shift=shift)
    prefixes=_compute_prefixes(b,x)
    for prefix in reversed(prefixes):
        b.isub(1,target=shift,control=prefix)
    _erase_prefixes(b,x,prefixes)
    b.free(shift)
