"""Two-stage coefficient-QROM grouping for reversible direct-phase oracles.

Group A+B then C, or A then B+C, instead of keeping A+B+C live or
reloading each coefficient separately. All variants implement exactly
the same quantized phase and dropped partial-product arithmetic.
"""
from __future__ import annotations
from .priority1_model import Approximation, coefficient_matrix
from .priority1_exact import full_cut_screen
from .priority1_truncated import _mul_cut
from .benchmark1_symmetric_square import square_symmetric
from .priority1_normalize import normalize_barrel, denormalize_barrel
from .benchmark1_compact_normalize import normalize_compact, denormalize_compact
from .hybrid_kickmix import _normalize, _denormalize, _phase
from .hybrid_model import NORMAL_BITS


def build_grouped(config: Approximation, drops=(10,12,11), *,
                  grouping="ab_c", normalization="serial", square_method="rows"):
    import kickmix as km
    if config.degree!=2 or grouping not in ("ab_c","a_bc"):
        raise ValueError("Grouping must be ab_c or a_bc for a quadratic")
    if normalization not in ("serial","barrel","compact"):
        raise ValueError("Normalization must be serial, barrel or compact")
    if square_method not in ("rows", "symmetric"):
        raise ValueError("square_method must be rows or symmetric")
    dl,ds,dq=(int(v) for v in drops)
    cert=full_cut_screen(config,(dl,ds,dq))
    if not cert["passes"] or cert["worst_error_over_tolerance"]>=.85:
        raise ValueError("Phase approximation lacks a full-domain safety margin")
    mat=coefficient_matrix(config)
    width=[max(1,int(mat[:,i].max()).bit_length()) for i in range(3)]
    lo=config.low_bits
    b=km.CircuitBuilder()
    x=b.create_quantum_register(NORMAL_BITS,name="x")
    norm,undo=(
        (_normalize,_denormalize) if normalization=="serial"
        else (normalize_compact,denormalize_compact) if normalization=="compact"
        else (normalize_barrel,denormalize_barrel))
    flags,shift=norm(b,x)

    def load(indices):
        table=[]
        for row in mat:
            val,pos=0,0
            for i in indices:
                val|=int(row[i])<<pos
                pos+=width[i]
            table.append(val)
        storage=b.alloc_qubits(sum(width[i] for i in indices))
        b.init_lookup(table=table,address=x[lo:],target=storage)
        outputs=[]
        pos=0
        for i in indices:
            outputs.append(storage[pos:pos+width[i]])
            pos+=width[i]
        return table,storage,outputs

    def release(table,storage):
        b.del_lookup(table=table,address=x[lo:],target=storage)
        b.free(storage)

    def phase_constant(coeff):
        _phase(b,coeff,shift,config.precision,+1)

    def phase_linear(coeff):
        prod=b.alloc_qubits(len(coeff)+lo-dl)
        b.reset(prod)
        _mul_cut(b,coeff,x[:lo],prod,dl)
        _phase(b,prod,shift,config.precision+lo-dl,-1)
        _mul_cut(b,coeff,x[:lo],prod,dl,inverse=True)
        b.free(prod)

    def phase_quadratic(coeff):
        sq=b.alloc_qubits(2*lo-ds)
        b.reset(sq)
        if square_method == "symmetric":
            square_symmetric(b,x[:lo],sq,ds)
        else:
            _mul_cut(b,x[:lo],x[:lo],sq,ds)
        prod=b.alloc_qubits(len(coeff)+len(sq)-dq)
        b.reset(prod)
        _mul_cut(b,coeff,sq,prod,dq)
        _phase(b,prod,shift,config.precision+2*lo-ds-dq,+1)
        _mul_cut(b,coeff,sq,prod,dq,inverse=True)
        b.free(prod)
        if square_method == "symmetric":
            square_symmetric(b,x[:lo],sq,ds,inverse=True)
        else:
            _mul_cut(b,x[:lo],x[:lo],sq,ds,inverse=True)
        b.free(sq)

    if grouping=="ab_c":
        with b.mark("coefficient_A_B"):
            table,store,(a,beta)=load((0,1))
            phase_constant(a)
            phase_linear(beta)
            release(table,store)
        with b.mark("coefficient_C"):
            table,store,(gamma,)=load((2,))
            phase_quadratic(gamma)
            release(table,store)
    else:
        with b.mark("coefficient_A"):
            table,store,(a,)=load((0,))
            phase_constant(a)
            release(table,store)
        with b.mark("coefficient_B_C"):
            table,store,(beta,gamma)=load((1,2))
            phase_linear(beta)
            phase_quadratic(gamma)
            release(table,store)
    undo(b,x,flags,shift)
    return b.finish_circuit()
