"""Continuous cold-sheet equation residuals from the existing actual cache.

Every binary arithmetic interval operation is outward rounded. The force
remainder uses the proved Helmholtz/homotopy inequalities; no new force solve,
Fourier transform, particle deposition or trajectory advancement is performed.
"""
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from mpmath import iv

from p9_cic_replacement_validation import original_trajectory, WORK
from cosmic_cg.pm_certificate import coefficient_intervals
from p9_continuous_force_bound import nonlinear_coefficient


COUNTS = Counter()


def exact(x):
    v=np.asarray(x,dtype=float)
    return v,v


def down(x):
    a=np.asarray(x,dtype=float)
    COUNTS['directed_roundings']+=a.size
    return np.nextafter(a,-math.inf)


def up(x):
    a=np.asarray(x,dtype=float)
    COUNTS['directed_roundings']+=a.size
    return np.nextafter(a,math.inf)


def add(a,b):
    n=np.broadcast_shapes(a[0].shape,b[0].shape)
    COUNTS['real_additions']+=2*math.prod(n)
    return down(a[0]+b[0]),up(a[1]+b[1])


def sub(a,b):
    n=np.broadcast_shapes(a[0].shape,b[0].shape)
    COUNTS['real_additions']+=2*math.prod(n)
    return down(a[0]-b[1]),up(a[1]-b[0])


def mul(a,b):
    if b[0].ndim==0 and b[1].ndim==0 and b[0]==b[1]:
        scalar=float(b[0]);COUNTS['comparisons']+=1
        COUNTS['real_multiplications']+=2*a[0].size
        if scalar>=0:
            return down(a[0]*scalar),up(a[1]*scalar)
        return down(a[1]*scalar),up(a[0]*scalar)
    if a[0].ndim==0 and a[1].ndim==0 and a[0]==a[1]:
        return mul(b,a)
    values=[x*y for x in a for y in b]
    n=values[0].size
    COUNTS['real_multiplications']+=4*n
    COUNTS['comparisons']+=6*n
    return down(np.minimum.reduce(values)),up(np.maximum.reduce(values))


def div_positive(a,b):
    assert np.all(b[0]>0)
    if b[0].ndim==0 and b[1].ndim==0 and b[0]==b[1]:
        COUNTS['comparisons']+=1
        COUNTS['real_divisions']+=2*a[0].size
        return down(a[0]/b[0]),up(a[1]/b[1])
    values=[x/y for x in a for y in b]
    n=values[0].size
    COUNTS['real_divisions']+=4*n
    COUNTS['comparisons']+=6*n
    return down(np.minimum.reduce(values)),up(np.maximum.reduce(values))


def square_root(a):
    lo=np.maximum(a[0],0.)
    hi=np.maximum(a[1],0.)
    COUNTS['comparisons']+=2*lo.size
    COUNTS['square_roots']+=2*lo.size
    return np.maximum(0.,down(np.sqrt(lo))),up(np.sqrt(hi))


def sum_all(a):
    low,high=[v.reshape(-1) for v in a]
    while low.size>1:
        pairs=low.size//2
        COUNTS['real_additions']+=2*pairs
        next_low=down(low[:2*pairs:2]+low[1:2*pairs:2])
        next_high=up(high[:2*pairs:2]+high[1:2*pairs:2])
        if low.size%2:
            next_low=np.concatenate((next_low,low[-1:]))
            next_high=np.concatenate((next_high,high[-1:]))
        low,high=next_low,next_high
    return np.asarray(low[0]),np.asarray(high[0])


def roll(a,shift,axis):
    return tuple(np.roll(v,shift,axis) for v in a)


def select(a,index):
    return tuple(v[index] for v in a)


def mass(a,axes):
    out=a
    for axis in axes:
        out=div_positive(add(add(roll(out,1,axis),mul(exact(4.),out)),
                             roll(out,-1,axis)),exact(6.))
    return out


def inner(a,b):
    return div_positive(sum_all(mul(a,b)),exact(float(a[0].size)))


def field_norm(a):
    energy=exact(0.)
    for b in range(3):
        f=select(a,b)
        energy=add(energy,inner(f,mass(f,range(3))))
    return square_root(energy)


def particle_rms(a):
    return square_root(div_positive(sum_all(mul(a,a)),exact(float(a[0].shape[0]))))


def grad_curl_norms(a,h):
    """Exact Q1 weak-derivative quadratic forms, enclosed by intervals.

    Curl energy = gradient energy - divergence energy on the periodic box.
    Mixed divergence terms use centred differences in the two differentiated
    axes and the Q1 mass operator in the remaining axis.
    """
    gradient=exact(0.);diagonal=[];central=[]
    for b in range(3):
        f=select(a,b)
        for d in range(3):
            difference=div_positive(sub(roll(f,-1,d),f),exact(h))
            value=inner(difference,mass(difference,[e for e in range(3) if e!=d]))
            gradient=add(gradient,value)
            if b==d:diagonal.append(value)
        central.append(div_positive(sub(roll(f,-1,b),roll(f,1,b)),exact(2*h)))
    divergence=exact(0.)
    for value in diagonal:divergence=add(divergence,value)
    for d in range(3):
        for e in range(d+1,3):
            cross=inner(central[d],mass(central[e],[f for f in range(3) if f not in (d,e)]))
            divergence=add(divergence,mul(exact(2.),cross))
    return square_root(gradient),square_root(sub(gradient,divergence))


def refine(a,centred=False):
    out=a
    for axis in range(1,4):
        shape=list(out[0].shape);shape[axis]*=2
        low=np.empty(shape);high=np.empty(shape)
        even=[slice(None)]*4;odd=even.copy()
        even[axis]=slice(0,None,2);odd[axis]=slice(1,None,2)
        average=div_positive(add(out,roll(out,1 if centred else -1,axis)),exact(2.))
        for array,plain,mid in zip((low,high),out,average):
            array[tuple(odd if centred else even)]=plain
            array[tuple(even if centred else odd)]=mid
        out=low,high
    return out


def magnitude_upper(a):
    COUNTS['absolute_values']+=2*a[0].size
    COUNTS['comparisons']+=a[0].size
    return np.maximum(abs(a[0]),abs(a[1]))


def max_vector_norm(a):
    mag=magnitude_upper(a)
    sq=mul(exact(mag),exact(mag))
    total=exact(np.zeros_like(mag[0]))
    for b in range(3):total=add(total,select(sq,b))
    COUNTS['comparisons']+=total[1].size-1
    maximum=float(np.max(total[1]))
    return float(square_root(exact(maximum))[1])


def mean_norm(a):
    means=[div_positive(sum_all(select(a,b)),exact(float(a[0][b].size))) for b in range(3)]
    vec=tuple(np.asarray([m[k] for m in means]) for k in range(2))
    return float(square_root(sum_all(mul(vec,vec)))[1])


def stored_interval(value):
    return float(value.a),float(value.b)


def growing_momentum_factor(a0,params):
    """Alternating binomial integral, with explicit next-term enclosure."""
    iv.dps=40
    omega=iv.mpf(params.omega_b)+iv.mpf(params.omega_sigma)
    a=iv.mpf(a0);z=iv.mpf(params.omega_lambda)*a**3/omega
    assert float(z.a)>0 and float(z.b)<.02
    coefficient=iv.mpf(1);total=iv.mpf(0)
    for n in range(20):
        term=coefficient/(iv.mpf(2.5)+3*n)
        total=total+(-term if n%2 else term)
        coefficient=coefficient*(iv.mpf(1.5)+n)*z/(n+1)
    # n=20 is positive and bounds the alternating remainder.
    next_term=coefficient/(iv.mpf(2.5)+60)
    integral=iv.mpf([total.a,(total+next_term).b])
    f=-iv.mpf(1.5)/(1+z)+1/(iv.sqrt((1+z)**3)*integral)
    factor=iv.sqrt(omega*a*(1+z))*f
    lower=float(np.nextafter(float(factor.a),-math.inf))
    upper=float(np.nextafter(float(factor.b),math.inf))
    return exact(lower)[0],exact(upper)[1]


def state_shape(a,n):
    return tuple(v.reshape(n,n,n,3).transpose(3,0,1,2) for v in a)


def particle_shape(a):
    return tuple(v.transpose(1,2,3,0).reshape(-1,3) for v in a)


def main():
    COUNTS.clear()
    cfg,params,times,x,p,g,D,s,target,_=original_trajectory()
    n=cfg.mesh_size;h=cfg.box_size_mpch/n
    old=json.loads((WORK/'p9_time_spatial_bounds_20261004.json').read_text(encoding='utf-8'))
    assert old['initial_state_sha256']==target['initial_state_sha256']
    assert old['operator_sha256']==target['finite_operator_data_sha256']
    L=old['uniform_continuous_reconstruction_Lip_upper']
    iv.dps=40
    A=float(np.nextafter(float((1/iv.sqrt((1-iv.mpf(L))**3)).b),math.inf))
    C=nonlinear_coefficient(L)
    inverse_k=float(np.nextafter(float((iv.mpf(cfg.box_size_mpch)/(2*iv.pi)).b),math.inf))
    source=mul(exact(D),exact(s))
    centres=source
    for d in range(1,4):centres=div_positive(add(centres,roll(centres,-1,d)),exact(2.))
    interpolation_error=float(field_norm(sub(refine(source),refine(centres,True)))[1])
    source_norm=float(field_norm(source)[1])
    source_gradient,source_curl=grad_curl_norms(source,h)
    source_amplitude=max_vector_norm(source)
    f0=growing_momentum_factor(float(times[0]),params)
    source_p_error=float(add(mul(exact(float(f0[1])),exact(interpolation_error)),
        particle_rms(sub(mul(f0,particle_shape(centres)),exact(p[0]))))[1])
    increments=[];sampled=[];norms=[]
    for j in range(len(times)):
        increment=sub(exact(x[j]),exact(x[0]))
        increments.append(increment)
        sampled.append(add(increment,particle_shape(centres)))
        field=state_shape(increment,n)
        grad,curl=grad_curl_norms(field,h)
        norm=float(add(exact(source_norm),exact(float(field_norm(field)[1])))[1])
        gn=float(add(exact(float(source_gradient[1])),exact(float(grad[1])))[1])
        cn=float(add(exact(float(source_curl[1])),exact(float(curl[1])))[1])
        amplitude=float(add(exact(source_amplitude),exact(max_vector_norm(field)))[1])
        mean=mean_norm(add(source,field))
        transverse=float(add(exact(mean),mul(exact(inverse_k),exact(cn)))[1])
        remainder=float(mul(exact(A),add(mul(exact(C),exact(norm)),
                                        mul(exact(amplitude),exact(gn))))[1])
        norms.append({'a':float(times[j]),'psi_L2_upper':norm,'gradient_L2_upper':gn,
                      'curl_L2_upper':cn,'psi_Linfinity_Euclidean_upper':amplitude,
                      'mean_upper':mean,'transverse_upper':transverse,
                      'nonlinear_and_composition_upper':remainder})
    coefficients=[coefficient_intervals(float(a),params) for a in times]
    rows=[]
    for j in range(len(times)-1):
        dt=sub(exact(times[j+1]),exact(times[j]))
        dx=div_positive(sub(exact(x[j+1]),exact(x[j])),dt)
        dp=div_positive(sub(exact(p[j+1]),exact(p[j])),dt)
        alpha=(np.asarray(coefficients[j+1][0][0]),np.asarray(coefficients[j][0][1]))
        beta=(np.asarray(coefficients[j+1][1][0]),np.asarray(coefficients[j][1][1]))
        rx=max(float(particle_rms(sub(dx,mul(alpha,exact(p[k]))))[1]) for k in (j,j+1))
        rx=float(add(exact(rx),mul(exact(float(alpha[1])),exact(source_p_error)))[1])
        local=max(float(particle_rms(sub(dp,mul(beta,sampled[k])))[1]) for k in (j,j+1))
        local=float(add(exact(local),mul(exact(float(beta[1])),exact(interpolation_error)))[1])
        # Each base norm is convex along the affine time reconstruction.
        # Products of distinct norms must use their separate endpoint maxima.
        # Endpoint maxima of the already multiplied expressions are not enough.
        upper=lambda key:max(norms[k][key] for k in (j,j+1))
        transverse=float(add(exact(upper('mean_upper')),
                            mul(exact(inverse_k),exact(upper('curl_L2_upper'))))[1])
        nonlinear=float(mul(exact(A),add(mul(exact(C),exact(upper('psi_L2_upper'))),
            mul(exact(upper('psi_Linfinity_Euclidean_upper')),
                exact(upper('gradient_L2_upper')))))[1])
        extra=float(add(exact(transverse),exact(nonlinear))[1])
        rp=float(add(exact(local),mul(exact(float(beta[1])),exact(extra)))[1])
        rows.append({'a_interval':[float(times[j]),float(times[j+1])],
            'r_X_L2_Lagrangian_upper':rx,'r_P_L2_Lagrangian_upper':rp,
            'momentum_local_polynomial_upper':local,
            'force_projection_composition_upper':extra,
            'alpha_interval':[float(v) for v in alpha],'beta_interval':[float(v) for v in beta]})
    integral_x=exact(0.);integral_p=exact(0.)
    for row in rows:
        lo,hi=row['a_interval'];dt=sub(exact(hi),exact(lo))
        integral_x=add(integral_x,mul(dt,exact(row['r_X_L2_Lagrangian_upper'])))
        integral_p=add(integral_p,mul(dt,exact(row['r_P_L2_Lagrangian_upper'])))
    result={
        'scope':'Lagrangian L2 defects of original-initial-data continuous reconstruction in the unsoftened periodic Newton equations',
        'reconstruction_X':'q + original nodal psi0(q) + I_c[x(a)-x(a0)](q)',
        'reconstruction_P':'f0 psi0(q) + I_c[p(a)-p(a0)](q)',
        'force_normalization':'F_N=-grad Phi, Laplacian Phi=rho-1, zero mean force',
        'norm':'Euclidean vectors, box-volume-normalized L2 in continuous Lagrangian q',
        'time_interval':[float(times[0]),float(times[-1])],
        'initial_state_sha256':target['initial_state_sha256'],
        'operator_sha256':target['finite_operator_data_sha256'],
        'cache_sha256':hashlib.sha256((WORK/'p9_actual_trajectory_20261004.npz').read_bytes()).hexdigest(),
        'spatial_certificate_sha256':hashlib.sha256((WORK/'p9_time_spatial_bounds_20261004.json').read_bytes()).hexdigest(),
        'continuous_initial_momentum_factor_interval':[float(v) for v in f0],
        'initial_source_interpolation_L2_upper':interpolation_error,
        'initial_momentum_interpolation_L2_upper':source_p_error,
        'reconstruction_Lipschitz_Euclidean_upper':L,
        'Jacobian_pullback_factor_upper':A,
        'max_r_X_L2_upper':max(row['r_X_L2_Lagrangian_upper'] for row in rows),
        'max_r_P_L2_upper':max(row['r_P_L2_Lagrangian_upper'] for row in rows),
        'integrated_r_X_L2_upper':float(integral_x[1]),
        'integrated_r_P_L2_upper':float(integral_p[1]),
        'states':norms,'intervals':rows,
        'work':{
            'trajectory_advancements':0,'new_FFTs':0,'new_depositions':0,
            'previous_spatial_certificate_recomputed':False,
            'binary_primitive_operation_counts':dict(COUNTS),
            'coefficient_endpoint_evaluations':len(times),
            'additional_coefficient_arithmetic_upper':256*len(times),
            'additional_scalar_arithmetic_upper':4096,
            'interval_preprocessing':{
                'precision_decimal_digits':40,'growth_series_terms':20,
                'additional_interval_arithmetic_upper':512,'square_roots_upper':8,
                'pi_interval_reads':1,
                'transcendental_force_evaluations':0},
            'peak_spatial_refinement':'one initial 32^3 Q1 field pair; trajectory states remain at 16^3',
            'index_copies_and_hash_bytes_counted_as_memory_work_not_real_arithmetic':True},
        'Newtonian_evolved_trajectory_certified':False,
        'proof_boundary':'Residuals are continuous-function L2 bounds; a dynamical stability and density continuation argument is separate.',
        'all_bounds_finite':True,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    assert all(math.isfinite(result[key]) and result[key]>=0 for key in
               ('max_r_X_L2_upper','max_r_P_L2_upper','integrated_r_X_L2_upper','integrated_r_P_L2_upper'))
    path=WORK/'p9_continuous_residual_certificate_20261004.json'
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in (
        'max_r_X_L2_upper','max_r_P_L2_upper','integrated_r_X_L2_upper','integrated_r_P_L2_upper','work')},indent=2))


if __name__=='__main__':main()
