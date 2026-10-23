"""All-frequency Newton force bound on the actual continuous reconstruction.

The reconstruction is anchored to the original cold sheet, with displacement
increments taken from the real 4096-particle trajectory. This certifies a
spatial observable of that reconstructed state, not the Newtonian evolution.
"""
import itertools
import json
import math
from pathlib import Path

import numpy as np
from mpmath import iv

from p9_cic_replacement_validation import (original_trajectory,
    continuous_reconstruction_metrics, up, positive_sum, WORK,
    _spectral_setup, _multipliers, certified_force, bounded_fft, ArithmeticWork, U)
from cosmic_cg.pm_certificate import _add, _mul, _div, _down
from cosmic_cg.pm_spatial_certificate import (
    ComplexBall, ball_add, ball_mul, interval_ball,
    _vector_ball_bounds, _sum_square_bounds,
)


def positive_square_root_upper(value, work=None):
    if value < 0:
        raise ValueError('A nonnegative square-root argument is required')
    if work is not None:work.add(square_roots=1,directed_roundings=1)
    return up(math.sqrt(value))


def positive_square_root_lower(value, work=None):
    if value < 0:
        raise ValueError('A nonnegative square-root argument is required')
    if work is not None:work.add(square_roots=1)
    return max(0.,_down(math.sqrt(value),work))


def lower_difference(first, second, work=None):
    if work is not None:work.add(real_additions=1)
    return max(0.,_down(first-second,work))


def tail_norm_upper(norm_upper, retained_energy_lower, work=None):
    energy_upper=_mul(norm_upper,norm_upper,work)
    if work is not None:work.add(real_additions=1,directed_roundings=1)
    tail_energy=max(0.,up(energy_upper-retained_energy_lower))
    return positive_square_root_upper(tail_energy,work)


def sum_of_squared_norms_bounds(intervals,work=None):
    rows=[{'norm':value} for value in intervals]
    return _sum_square_bounds(rows,'norm',work)


def spline_norm_upper(fields,radius=0.,work=None):
    """L^-3/2 L2 norm of the periodic trilinear interpolant.

    The exact mass matrix is the tensor product of (1,4,1)/6.
    Its positive averaging and the final dot product are enclosed explicitly.
    """
    fields=np.asarray(fields,float)
    total=0.
    for field in fields:
        averaged=field.copy()
        for axis in range(3):
            averaged=(np.roll(averaged,1,axis)+4*averaged+np.roll(averaged,-1,axis))/6
        term=float(np.mean(field*averaged))
        magnitude=up(float(np.max(np.abs(field))))
        roundoff=_mul(64*U*field.size,_mul(magnitude,magnitude,work),work)
        roundoff=_add(roundoff,64*field.size*np.nextafter(0.,1.),work)
        total=_add(total,_add(term,roundoff,work),work)
        if work is not None:
            work.add(real_multiplications=4*field.size,
                real_additions=7*field.size-1,real_divisions=3*field.size,
                binary_scalings=1,absolute_values=field.size,comparisons=field.size-1)
    return _add(positive_square_root_upper(max(0.,total),work),
                _mul(positive_square_root_upper(3.,work),radius,work),work)


def interval_constant(value):
    result=interval_ball(value)
    return result.centre.real,result.radius


def spectral_factors(n,M):
    iv.dps=40
    factors={}
    for m in range(-M,M+1):
        theta=iv.pi*m/n
        sinc=iv.mpf(1) if m==0 else iv.sin(theta)/theta
        w,wr=interval_constant(sinc*sinc)
        re,rr=interval_constant(iv.cos(theta))
        im,ir=interval_constant(-iv.sin(theta))
        factors[m]=(w,wr,complex(re,im),up(rr+ir))
    return factors


def multiply_ball(z,r,f,fr,work=None):
    result=ball_mul(ComplexBall(complex(z),r),ComplexBall(complex(f),fr),work)
    return result.centre,result.radius


def vector_norm_interval(vector,radii,work=None):
    return _vector_ball_bounds([ComplexBall(complex(z),r) for z,r in zip(vector,radii)],work)


def nonlinear_coefficient(L):
    iv.dps=40
    ell=iv.mpf(L)
    result=(iv.mpf('1.5')*ell+ell**2+ell**3/4)/iv.sqrt((1-ell)**3)
    result+=ell*iv.sqrt((1+ell)**3)/2
    return up(float(result.b))


def validate_homotopy_bound():
    """Independent exact 1D Newton solution pulled back by x=q+eps sin(q)."""
    from numpy.polynomial.legendre import leggauss
    node,weight=leggauss(192)
    q=(node+1)*math.pi
    checks=[]
    for amp in (.003,.03,.15,.4):
        psi=amp*np.sin(q);jac=1+amp*np.cos(q)
        exact_difference=psi-amp*np.sin(q+psi)
        norm=math.sqrt(float(np.sum(weight*exact_difference**2*jac)/2))
        upper=nonlinear_coefficient(amp)*amp/math.sqrt(2)
        assert norm<upper
        checks.append({'amplitude':amp,'independent_force_difference_L2':norm,'analytic_bound':upper})
    return checks


def main():
    cfg,params,times,x,p,g,D,s,target,_=original_trajectory()
    reconstruction=continuous_reconstruction_metrics(times,x,p,D,s,cfg,target)
    n=cfg.mesh_size;G=n**3
    source=np.asarray(D,float)*s
    increment=(x[-1]-x[0]).reshape(n,n,n,3).transpose(3,0,1,2)
    setup=_spectral_setup(cfg);multipliers=_multipliers(setup,n)
    work=ArithmeticWork()
    _,grid_radius,grids=certified_force(x[-1]%cfg.box_size_mpch,cfg,multipliers,work)
    grids=np.asarray(grids)
    # Separate component splines share the same periodic norm; their
    # triangle bound is valid for the two different mesh origins.
    psi_norm=up(spline_norm_upper(source)+spline_norm_upper(increment))
    force_norm=spline_norm_upper(grids,grid_radius)
    L=reconstruction['uniform_reconstruction_Lip_Euclidean_upper']
    nonlinear=up(nonlinear_coefficient(L)*psi_norm)
    transforms=[];radii=[]
    for fields,radius in ((source,0.),(increment,0.),(grids,grid_radius)):
        modes=[];rs=[]
        for field in fields:
            mode,rad=bounded_fft(field,radius,work=work)
            modes.append(mode/G);rs.append(up(rad/G))
        transforms.append(np.asarray(modes));radii.append(rs)
    rows=[]
    # Full integer frequency cubes include both signs of Nyquist. The
    # trilinear interpolants have genuine aliases at every integer mode.
    for M in (2,4,8):
        factors=spectral_factors(n,M)
        error_lower=[];error_upper=[];psi_lower=[];force_lower=[];linear_force_lower=[]
        for m in itertools.product(range(-M,M+1),repeat=3):
            if m==(0,0,0):
                continue
            window,window_radius=1.,0.
            phase,phase_radius=1.+0j,0.
            for component in m:
                w,wr,z,zr=factors[component]
                window,window_radius=multiply_ball(window,window_radius,w,wr)
                phase,phase_radius=multiply_ball(phase,phase_radius,z,zr)
            index=tuple(component%n for component in m)
            psi=[];pr=[];force=[];fr=[]
            for b in range(3):
                delta,dr=multiply_ball(transforms[1][(b,*index)],radii[1][b],phase,phase_radius)
                total=transforms[0][(b,*index)]+delta
                tr=up(radii[0][b]+dr+16*U*max(1.,abs(total)))
                value,rad=multiply_ball(total,tr,window,window_radius)
                psi.append(value);pr.append(rad)
                value,rad=multiply_ball(transforms[2][(b,*index)],radii[2][b],window,window_radius)
                force.append(value);fr.append(rad)
            m2=sum(component*component for component in m)
            projected=[];projected_r=[]
            for b in range(3):
                value=sum((m[b]*m[d]/m2)*psi[d] for d in range(3))
                rad=up(sum(abs(m[b]*m[d]/m2)*pr[d] for d in range(3))+64*U*max(1.,abs(value)))
                projected.append(value);projected_r.append(rad)
            lo,hi=vector_norm_interval([projected[b]-force[b] for b in range(3)],
                                      [up(projected_r[b]+fr[b]+16*U) for b in range(3)])
            error_lower.append(lo*lo);error_upper.append(up(hi*hi))
            lo,_=vector_norm_interval(psi,pr);psi_lower.append(float(np.nextafter(lo*lo,-math.inf)))
            lo,_=vector_norm_interval(force,fr);force_lower.append(float(np.nextafter(lo*lo,-math.inf)))
            lo,_=vector_norm_interval(projected,projected_r);linear_force_lower.append(float(np.nextafter(lo*lo,-math.inf)))
        errlow=max(0.,float(np.nextafter(math.fsum(error_lower)*(1-16*U),-math.inf)))
        errhigh=positive_sum(error_upper)
        psi_energy_low=max(0.,float(np.nextafter(math.fsum(psi_lower)*(1-16*U),-math.inf)))
        force_energy_low=max(0.,float(np.nextafter(math.fsum(force_lower)*(1-16*U),-math.inf)))
        psi_tail=tail_norm_upper(psi_norm,psi_energy_low)
        force_tail=tail_norm_upper(force_norm,force_energy_low)
        tail_sum=_add(psi_tail,force_tail)
        whole_linear=positive_square_root_upper(_add(errhigh,_mul(tail_sum,tail_sum)))
        full_upper=_add(whole_linear,nonlinear)
        full_lower=lower_difference(positive_square_root_lower(errlow),nonlinear)
        linear_energy_lower=max(0.,_down(math.fsum(linear_force_lower)*(1-16*U)))
        newton_lower=lower_difference(positive_square_root_lower(linear_energy_lower),nonlinear)
        newton_upper=up(psi_norm+nonlinear)
        rows.append({'Fourier_cube':M,'modes':(2*M+1)**3-1,
            'finite_linear_discrepancy_L2_interval':[math.sqrt(errlow),math.sqrt(errhigh)],
            'unretained_displacement_L2_upper':psi_tail,'unretained_CIC_force_L2_upper':force_tail,
            'all_frequency_Newton_CIC_force_error_interval':[full_lower,full_upper],
            'Newton_force_norm_interval':[newton_lower,newton_upper],
            'relative_error_upper':None if newton_lower==0 else up(full_upper/newton_lower)})
    result={'target':'Original continuous cold sheet plus actual saved particle-displacement increments, final a=0.4',
        'initial_state_sha256':target['initial_state_sha256'],'upstream_run_id':params.source_run_id,
        'continuous_reconstruction_initial_error':0,
        'displacement_Lip_Euclidean_upper':L,'displacement_L2_upper':psi_norm,
        'density_lower':float(np.nextafter(float((1/(1+iv.mpf(L))**3).a),-math.inf)),
        'density_upper':up(float((1/(1-iv.mpf(L))**3).b)),
        'nonlinear_homotopy_remainder_L2_upper':nonlinear,
        'all_frequency_spatial_observable_enclosed':True,
        'Newtonian_evolved_trajectory_certified':False,'rows':rows,
        'independent_1d_Newton_force_checks':validate_homotopy_bound(),
        'work':{'additional_force_evaluations':1,'additional_3d_complex_ffts':13,
            'explicit_fft_and_force_arithmetic':work.counts,'finite_observable_mode_visits':sum(r['modes'] for r in rows),
            'scalar_observable_arithmetic_upper':8192*sum(r['modes'] for r in rows),
            'spline_norm_real_operations_upper':128*G,
            'interval_trigonometric_preprocessing_max_calls':6*(2*8+1),
            'scope':'one final reconstructed spatial state; spectral modes are reused, no new trajectory'} }
    output=WORK/'p9_continuous_force_bound_20261004.json'
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'output':str(output),'L':L,'psi_norm':psi_norm,'nonlinear_remainder':nonlinear,'rows':rows},indent=2))


if __name__=='__main__':
    main()
