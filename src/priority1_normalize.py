"""Alternative power-of-four normalization using a pairwise barrel rotation.

Computes the number of leading zero bit pairs *before* permuting x.
Rotating even and odd lanes separately by that amount is equivalent
to moving 2*s input bits left in a 22-bit circular register. Each
temporary predicate is cleared with its matching MBUC correction.
"""
from __future__ import annotations


def _zero_pair(builder, x, pair_index):
    high=2*pair_index+1
    low=2*pair_index
    builder.x(x[high]);builder.x(x[low])
    tmp=builder.alloc_qubits(1)[0]
    builder.init_and(x[high],x[low],tmp)
    builder.x(x[high]);builder.x(x[low])
    return tmp


def _erase_zero_pair(builder,x,pair_index,tmp):
    high=2*pair_index+1
    low=2*pair_index
    builder.x(x[high]);builder.x(x[low])
    builder.del_and(x[high],x[low],tmp)
    builder.x(x[high]);builder.x(x[low])
    builder.free(tmp)


def _lanes(x):
    import kickmix as km
    return (km.array([x[i] for i in range(0,22,2)]),
            km.array([x[i] for i in range(1,22,2)]))


def normalize_barrel(builder,x):
    if len(x)!=22:
        raise ValueError("Normalization expects the Benchmark I 22-bit register")
    scale=builder.alloc_qubits(4)
    builder.reset(scale)
    prefixes=[]
    for pair_index in range(10,0,-1):
        z=_zero_pair(builder,x,pair_index)
        if prefixes:
            prefix=builder.alloc_qubits(1)[0]
            builder.init_and(prefixes[-1],z,prefix)
            _erase_zero_pair(builder,x,pair_index,z)
        else:
            prefix=z
        prefixes.append(prefix)
        builder.iadd(1,target=scale,control=prefix)
    even,odd=_lanes(x)
    builder.left_rotate(even,shift=scale)
    builder.left_rotate(odd,shift=scale)
    return prefixes,scale


def denormalize_barrel(builder,x,prefixes,scale):
    even,odd=_lanes(x)
    builder.right_rotate(even,shift=scale)
    builder.right_rotate(odd,shift=scale)
    for idx in range(len(prefixes)-1,-1,-1):
        prefix=prefixes[idx]
        pair_index=10-idx
        builder.isub(1,target=scale,control=prefix)
        if idx:
            z=_zero_pair(builder,x,pair_index)
            builder.del_and(prefixes[idx-1],z,prefix)
            builder.free(prefix)
            _erase_zero_pair(builder,x,pair_index,z)
        else:
            _erase_zero_pair(builder,x,pair_index,prefix)
    builder.free(scale)
