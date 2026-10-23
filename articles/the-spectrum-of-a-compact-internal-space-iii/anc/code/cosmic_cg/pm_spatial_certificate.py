"""Initial Newtonian force observables of the generated continuous cold sheet.

The initial map and PM dynamics are unchanged.  This calculation encloses
Fourier projections of the continuous Newtonian force and of the actual
CIC force field.  It does not infer an all-mode or evolved continuum bound.
"""
from dataclasses import dataclass
from itertools import product
import math
import numpy as np
from mpmath import iv

from cosmic_cg.pm_certificate import (
    ArithmeticWork, U, ETA, _add, _mul, _div, _up, _down, _record,
    bounded_fft,
)


@dataclass(frozen=True)
class ComplexBall:
    centre: complex
    radius: float


def _norm1(z,work):
    _record(work,absolute_values=2)
    return _add(abs(z.real),abs(z.imag),work)


def ball_add(a,b,work):
    z=a.centre+b.centre
    _record(work,real_additions=2)
    magnitude=_add(_norm1(a.centre,work),_norm1(b.centre,work),work)
    radius=_add(_add(a.radius,b.radius,work),_mul(16*U,magnitude,work),work)
    return ComplexBall(z,_add(radius,16*ETA,work))


def ball_mul(a,b,work):
    z=a.centre*b.centre
    _record(work,real_multiplications=4,real_additions=2)
    aa,bb=_norm1(a.centre,work),_norm1(b.centre,work)
    radius=_add(_mul(aa,b.radius,work),_mul(bb,a.radius,work),work)
    radius=_add(radius,_mul(a.radius,b.radius,work),work)
    radius=_add(radius,_mul(32*U,_mul(aa,bb,work),work),work)
    return ComplexBall(z,_add(radius,32*ETA,work))


def interval_ball(value):
    lo,hi=float(value.a),float(value.b)
    centre=(lo+hi)/2
    # Conversion of interval endpoints and their midpoint to binary64.
    radius=_up(abs(hi-lo)+8*U*max(abs(lo),abs(hi))+16*ETA)
    return ComplexBall(complex(centre),radius)


def phase_ball(angle):
    c,s=interval_ball(iv.cos(angle)),interval_ball(iv.sin(angle))
    return ComplexBall(complex(c.centre.real,s.centre.real),_add(c.radius,s.radius))


def _sqrt_up(value,work):
    _record(work,square_roots=1)
    return _up(math.sqrt(value),work)


def _vector_norm_upper(array,work):
    """Uniform Euclidean norm for a vector-valued real grid."""
    squared=np.sum(array*array,axis=0)
    _record(work,real_multiplications=array.size,
            real_additions=2*squared.size,
            comparisons=squared.size-1)
    value=_mul(float(np.max(squared)),1+32*U,work)
    return _sqrt_up(value,work)


def _abs_bounds(z,work):
    # Scaling prevents squared small coefficients from underflowing.
    big=max(abs(z.real),abs(z.imag))
    small=min(abs(z.real),abs(z.imag))
    if big==0:
        return 0.,0.
    ratio=small/big
    _record(work,real_multiplications=2,real_additions=1,real_divisions=1,
            square_roots=1,absolute_values=4,comparisons=2)
    value=big*math.sqrt(1+ratio*ratio)
    allowance=_add(_mul(16*U,_norm1(z,work),work),16*ETA,work)
    return max(0.,_down(value-allowance,work)),_add(value,allowance,work)


def _vector_ball_bounds(balls,work):
    lower,upper=0.,0.
    for ball in balls:
        lo,hi=_abs_bounds(ball.centre,work)
        lo=max(0.,_down(lo-ball.radius,work))
        hi=_add(hi,ball.radius,work)
        lower=max(0.,_down(lower+_down(lo*lo,work),work))
        upper=_add(upper,_mul(hi,hi,work),work)
        _record(work,real_multiplications=1,real_additions=2)
    _record(work,square_roots=2)
    return max(0.,_down(math.sqrt(lower),work)),_up(math.sqrt(upper),work)


def _sum_square_bounds(rows,key,work):
    lo,hi=0.,0.
    for row in rows:
        lower,upper=row[key]
        lo=max(0.,_down(lo+_down(lower*lower,work),work))
        hi=_add(hi,_mul(upper,upper,work),work)
        _record(work,real_additions=1,real_multiplications=1)
    _record(work,square_roots=2)
    return max(0.,_down(math.sqrt(lo),work)),_up(math.sqrt(hi),work)


def initial_force_certificate(initial_sheet,config,multipliers,density_spectrum,*,bands=(1,2)):
    """Enclose the actual initial force discrepancy in finite Fourier bands.

    Three linear and forty-eight quadratic transforms are shared by every
    requested mode.  The initial density transform is reused from the
    particle-trajectory certificate.
    """
    n=config.mesh_size
    if not bands or any(type(b) is not int or not 0<b<n//2 for b in bands):
        raise ValueError("The observable bands must be strictly below Nyquist")
    iv.dps=40
    work=ArithmeticWork()
    volume=n**3
    scale=float(initial_sheet["displacement_scale"])
    source=np.asarray(initial_sheet["displacement_grid"],dtype=float)
    if source.shape!=(3,n,n,n) or not np.all(np.isfinite(source)) or not math.isfinite(scale):
        raise ValueError("A finite three-component displacement grid is required")
    if not np.array_equal(source,source.astype(np.float32).astype(float)) or scale!=float(np.float32(scale)):
        raise ValueError("Use the binary32 source grid and normalization of the actual initializer")
    # Source values and scale are binary32; their products fit binary64.
    psi=source*scale
    _record(work,real_multiplications=psi.size)
    amplitude=_vector_norm_upper(psi,work)
    squares=np.sum(psi*psi,axis=0)
    _record(work,real_multiplications=psi.size,real_additions=3*volume-1,binary_scalings=1)
    second_moment=_mul(float(np.sum(squares)/volume),1+16*volume*U,work)
    edge=[]
    for d in range(3):
        difference=psi-np.roll(psi,1,axis=d+1)
        _record(work,real_additions=psi.size)
        edge.append(_mul(_vector_norm_upper(difference,work),1+16*U,work))

    linear=[]
    for component in psi:
        centre,radius=bounded_fft(component,work=work)
        linear.append((centre/volume,_div(radius,volume,work)))
        _record(work,binary_scalings=2*volume)
    root3=iv.sqrt(3)
    gauss_iv=(iv.mpf('0.5')-1/(2*root3),iv.mpf('0.5')+1/(2*root3))
    gauss=[interval_ball(node).centre.real for node in gauss_iv]
    quadratics=[]
    for node_ids in product((0,1),repeat=3):
        t=[gauss[j] for j in node_ids]
        interpolated=np.zeros_like(psi)
        for shift in product((0,1),repeat=3):
            coefficient=math.prod(t[d] if shift[d] else 1-t[d] for d in range(3))
            interpolated+=coefficient*np.roll(psi,tuple(-v for v in shift),axis=(1,2,3))
            _record(work,real_multiplications=3*volume+2,real_additions=3*volume+3)
        # Eight positive trilinear weights and enclosed Gauss nodes.
        interpolation_radius=_mul(512*U,max(amplitude,ETA),work)
        product_radius=_add(_mul(2*amplitude,interpolation_radius,work),
                            _mul(interpolation_radius,interpolation_radius,work),work)
        product_radius=_add(product_radius,_mul(16*U,_mul(amplitude,amplitude,work),work),work)
        coefficients={}
        for i in range(3):
            for j in range(i,3):
                centre,radius=bounded_fft(interpolated[i]*interpolated[j],product_radius,work=work)
                coefficients[i,j]=(centre/volume,_div(radius,volume,work))
                _record(work,real_multiplications=volume,binary_scalings=2*volume)
        quadratics.append((node_ids,coefficients))

    rows=[]
    maximum=max(bands)
    for mode in product(range(-maximum,maximum+1),repeat=3):
        if mode==(0,0,0):
            continue
        index=tuple(m%n for m in mode)
        theta=[2*iv.pi*m/n for m in mode]
        k_iv=[2*iv.pi*m/config.box_size_mpch for m in mode]
        k=[interval_ball(value) for value in k_iv]
        k2=sum(value*value for value in k_iv)
        k_norm=interval_ball(iv.sqrt(k2))
        kup=_add(abs(k_norm.centre),k_norm.radius,work)
        window=iv.mpf(1)
        for m in mode:
            if m:
                angle=iv.pi*m/n
                window*=((iv.sin(angle)/angle)**2)
        window_ball=interval_ball(window)
        linear_density=ComplexBall(0j,0.)
        for d in range(3):
            fft,radius=linear[d]
            term=ball_mul(k[d],ComplexBall(complex(fft[index]),radius),work)
            linear_density=ball_add(linear_density,term,work)
        linear_density=ball_mul(ComplexBall(-1j,0.),ball_mul(window_ball,linear_density,work),work)
        quadratic=ComplexBall(0j,0.)
        for ids,coefficients in quadratics:
            angle=-sum(theta[d]*gauss_iv[ids[d]] for d in range(3))
            phase=phase_ball(angle)
            local=ComplexBall(0j,0.)
            for i in range(3):
                for j in range(i,3):
                    fft,radius=coefficients[i,j]
                    factor=ball_mul(k[i],k[j],work)
                    if i!=j:
                        factor=ball_mul(ComplexBall(2.,0.),factor,work)
                    local=ball_add(local,ball_mul(factor,ComplexBall(complex(fft[index]),radius),work),work)
            quadratic=ball_add(quadratic,ball_mul(ComplexBall(0.125,0.),ball_mul(phase,local,work),work),work)
        density=ball_add(linear_density,ball_mul(ComplexBall(-0.5,0.),quadratic,work),work)

        # Third-order displacement remainder, averaged before bounding.
        remainder=_div(_mul(_mul(_mul(kup,kup,work),kup,work),
                            _mul(amplitude,second_moment,work),work),6,work)
        gauss_error=0.
        zmax=_mul(kup,amplitude,work)
        for d in range(3):
            angle=interval_ball(theta[d])
            th=_add(abs(angle.centre),angle.radius,work)
            th2=_mul(th,th,work)
            th3=_mul(th2,th,work)
            th4=_mul(th2,th2,work)
            dz=_mul(kup,edge[d],work)
            term1=_mul(th4,_mul(zmax,zmax,work),work)
            term2=_mul(8,_mul(th3,_mul(zmax,dz,work),work),work)
            term3=_mul(12,_mul(th2,_mul(dz,dz,work),work),work)
            gauss_error=_add(gauss_error,_div(_add(_add(term1,term2,work),term3,work),4320,work),work)
        density=ComplexBall(density.centre,_add(density.radius,
                    _add(remainder,_mul(0.5,gauss_error,work),work),work))

        modes=density_spectrum["modes"]
        mr=density_spectrum["radius"]
        pm_density=ComplexBall(complex(modes[index]/volume),_div(mr,volume,work))
        newton,pm,discrepancy=[],[],[]
        for d in range(3):
            coefficient=interval_ball(k_iv[d]/k2)
            exact_force=ball_mul(ComplexBall(1j,0.),ball_mul(coefficient,density,work),work)
            stored=complex(multipliers[d][index])
            stored_radius=_mul(16*U,_norm1(stored,work),work)
            discrete_force=ball_mul(window_ball,ball_mul(ComplexBall(stored,stored_radius),pm_density,work),work)
            newton.append(exact_force)
            pm.append(discrete_force)
            discrepancy.append(ball_add(discrete_force,ball_mul(ComplexBall(-1.,0.),exact_force,work),work))
        rows.append({"mode":list(mode),"force_error_norm":_vector_ball_bounds(discrepancy,work),
                     "Newton_force_norm":_vector_ball_bounds(newton,work),
                     "Newton_density_coefficient":{"real":density.centre.real,
                         "imaginary":density.centre.imag,"radius":density.radius},
                     "density_displacement_remainder_upper":remainder,
                     "density_Gauss_remainder_upper":_mul(0.5,gauss_error,work)})

    summaries=[]
    for band in sorted(set(bands)):
        selected=[row for row in rows if max(abs(m) for m in row["mode"])<=band]
        difference=_sum_square_bounds(selected,"force_error_norm",work)
        newton=_sum_square_bounds(selected,"Newton_force_norm",work)
        relative=None
        if newton[0]>0:
            lower=max(0.,_down(difference[0]/newton[1],work))
            upper=_div(difference[1],newton[0],work)
            relative=[lower,upper]
        summaries.append({"mode_box_radius":band,"nonzero_modes":len(selected),
                          "normalized_L2_force_error":difference,
                          "normalized_L2_Newton_force":newton,
                          "relative_L2_force_error":relative})
    return {"target":"initial Newtonian and CIC force Fourier projections of the original continuous cold sheet",
            "initial_map":"X(q)=q+s I_h D(q), q uniform in the periodic box",
            "force_comparison":"continuous Newtonian force versus the trilinearly interpolated PM force from the actual particles",
            "norm":"volume-normalized L2 of the stated nonzero Fourier box",
            "displacement_amplitude_upper":amplitude,
            "displacement_second_moment_upper":second_moment,
            "cell_edge_displacement_upper":edge,
            "bands":summaries,"mode_enclosures":rows,
            "work":{"additional_3d_complex_FFTs":51,"reused_density_transforms":1,
                    "new_particle_trajectories":0,"explicit_arithmetic_counts":work.counts,
                    "remaining_scalar_arithmetic_upper":4096*len(rows),
                    "interval_preprocessing":{"decimal_digits":40,
                        "sin_cos_enclosures":16*len(rows)+sum(sum(m!=0 for m in row["mode"]) for row in rows),
                        "square_root_enclosures":len(rows)+1},
                    "retained_quadratic_spectral_arrays":48},
            "all_mode_spatial_error_certified":False,
            "evolved_Newtonian_trajectory_certified":False,
            "all_checks_passed":all(math.isfinite(b["normalized_L2_force_error"][1]) for b in summaries)}
