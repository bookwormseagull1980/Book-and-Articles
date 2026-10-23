"""A posteriori enclosure of the finite CIC particle--mesh equations.

The target is the periodic CIC operator with the published spectral
multipliers.  Binary64 elementary operations and square roots obey the
IEEE round-to-nearest model.  The enclosure includes an absolute underflow
allowance; nonfinite results are rejected.  Fourier roundoff is bounded
for the explicit radix-two butterflies below.  The initial stored
particle state is the initial datum.  No continuum extrapolation enters.
"""
from __future__ import annotations

import math
import hashlib
from dataclasses import dataclass, field
from functools import lru_cache
import numpy as np
from mpmath import iv

from cosmic_cg.cosmological_dynamics import expansion, drift_kick_integrals
from cosmic_cg.pm_nbody import (
    _initial_state, _spectral_setup, _particle_acceleration,
)

U = np.finfo(float).eps / 2
ETA = float(np.nextafter(0.0,1.0))


@dataclass
class ArithmeticWork:
    """Counts for explicit real arithmetic, independent of elapsed time."""
    counts: dict = field(default_factory=dict)

    def add(self, **counts):
        for key, value in counts.items():
            self.counts[key] = self.counts.get(key,0)+int(value)


def _record(work, **counts):
    if work is not None:
        work.add(**counts)


def _up(x, work=None):
    _record(work,directed_roundings=1)
    return float(np.nextafter(float(x), math.inf))


def _down(x, work=None):
    _record(work,directed_roundings=1)
    return float(np.nextafter(float(x), -math.inf))


def _add(a,b,work=None):
    _record(work,real_additions=1)
    return _up(a+b,work)


def _mul(a,b,work=None):
    _record(work,real_multiplications=1)
    return _up(a*b,work)


def _div(a,b,work=None):
    _record(work,real_divisions=1)
    return _up(a/b,work)


def _magnitude(z,work=None):
    """Complex modulus majorant using only absolute values and addition."""
    _record(work,real_additions=z.size,absolute_values=2*z.size,comparisons=max(0,z.size-1))
    return _up(float(np.max(np.abs(z.real)+np.abs(z.imag))),work)


def coefficient_intervals(a,params,work=None):
    """Outward intervals for alpha and beta from positive operations."""
    # alpha=(Omega_m a^3+Omega_Lambda a^6)^(-1/2),
    # beta=(3/2) Omega_m a alpha.  All data are stored binary numbers.
    def pair_mul(x,y):
        _record(work,real_multiplications=2)
        return _down(x[0]*y[0],work),_up(x[1]*y[1],work)
    _record(work,real_additions=2)
    omega=(_down(params.omega_b+params.omega_sigma,work),
           _up(params.omega_b+params.omega_sigma,work))
    aa=(a,a)
    a3=pair_mul(pair_mul(aa,aa),aa)
    a6=pair_mul(a3,a3)
    first=pair_mul(omega,a3)
    second=pair_mul((params.omega_lambda,params.omega_lambda),a6)
    _record(work,real_additions=2,square_roots=2,real_divisions=2)
    low=_down(first[0]+second[0],work)
    high=_up(first[1]+second[1],work)
    alpha=(_down(1/_up(math.sqrt(high),work),work),
           _up(1/_down(math.sqrt(low),work),work))
    beta=pair_mul(pair_mul(pair_mul((1.5,1.5),omega),aa),alpha)
    return alpha,beta


def exponential_upper(z,work=None):
    """Sixteen Taylor terms with an outward geometric remainder."""
    if not 0 <= z < 1:
        raise ValueError("Subdivide until the step stability product is below one")
    term,total=1.,1.
    for q in range(1,17):
        term=_div(_mul(term,z,work),q,work)
        total=_add(total,term,work)
    numerator=_div(_mul(term,z,work),17,work)
    ratio=_div(z,18,work)
    _record(work,real_additions=1)
    denominator=_down(1-ratio,work)
    return _add(total,_div(numerator,denominator,work),work)


@lru_cache(maxsize=None)
def _twiddles(n, inverse):
    # Interval trigonometry is performed once for each mesh, not per step.
    iv.dps = 40
    out, radius = [], 0.0
    for j in range(n // 2):
        angle = (2 if inverse else -2) * iv.pi * j / n
        c, s = iv.cos(angle), iv.sin(angle)
        lo_c, hi_c = float(c.a), float(c.b)
        lo_s, hi_s = float(s.a), float(s.b)
        z = complex((lo_c + hi_c) / 2, (lo_s + hi_s) / 2)
        radius = max(radius, abs(hi_c-lo_c)+abs(hi_s-lo_s)+8*U)
        out.append(z)
    return np.asarray(out), _up(radius)


def bounded_fft(data, radius=0.0, *, inverse=False, work=None):
    """Three-dimensional DFT and a uniform complex absolute-error bound."""
    z = np.asarray(data, dtype=np.complex128).copy()
    if z.ndim != 3 or len(set(z.shape)) != 1:
        raise ValueError("A cubic mesh is required")
    n = z.shape[0]
    if n < 2 or n & (n-1) or not np.all(np.isfinite(z)):
        raise ValueError("The mesh must be finite and a power of two")
    bits = int(math.log2(n))
    order = np.array([int(f"{j:0{bits}b}"[::-1], 2) for j in range(n)])
    r = _up(radius,work)
    for axis in range(3):
        z = np.take(z, order, axis=axis)
        z = np.moveaxis(z, axis, -1)
        for width in (2**j for j in range(1, bits+1)):
            view = z.reshape(*z.shape[:-1], n//width, width)
            first = view[..., :width//2].copy()
            second = view[..., width//2:].copy()
            tw, tr = _twiddles(width, inverse)
            magnitude = _magnitude(z,work)
            product = second * tw
            view[..., :width//2] = first + product
            view[..., width//2:] = first - product
            # Complex multiplication and addition need at most eight real
            # operations.  The factor 32 covers their Euclidean norm bound
            # and the outward rounding of this scalar recurrence.
            r = _add(_add(_mul(_add(2,tr,work),r,work),
                         _mul(_add(tr,32*U,work),magnitude,work),work),32*ETA,work)
            butterflies=z.size//2
            _record(work,real_multiplications=4*butterflies,
                    real_additions=6*butterflies,complex_butterflies=butterflies)
        z = np.moveaxis(z, -1, axis)
    if inverse:
        z /= n**3  # exact binary scaling
        r = _div(r,n**3,work)
        _record(work,binary_scalings=2*z.size)
    if not np.all(np.isfinite(z)) or not math.isfinite(r):
        raise ArithmeticError("Nonfinite Fourier enclosure")
    return z, r


def _multipliers(setup, n):
    """Full Hermitian multiplier of the real half-spectrum PM operator."""
    k2 = np.asarray(setup["k2"], dtype=float)
    reverse = (-np.arange(n)) % n
    result = []
    for axis in range(3):
        shape = [1, 1, 1]
        shape[axis] = -1
        half = 1j*np.asarray(setup[("kx", "ky", "kz")[axis]], dtype=float).reshape(shape)/k2
        half[0,0,0] = 0
        full = np.zeros((n,n,n), complex)
        full[:,:,:n//2+1] = half
        for j in range(1, n//2):
            full[:,:,n-j] = np.conj(half[:,:,j][np.ix_(reverse,reverse)])
        for j in (0,n//2):
            plane = full[:,:,j]
            full[:,:,j] = (plane + np.conj(plane[np.ix_(reverse,reverse)]))/2
        result.append(full)
    return result


def _weights(positions, config,work=None):
    h = config.box_size_mpch/config.mesh_size
    # The validation lattice uses an exact binary mesh spacing.  This
    # makes cell identification exact away from an interval boundary.
    if math.frexp(h)[0] != 0.5:
        raise ValueError("The certified mesh spacing must be a power of two")
    coord = np.asarray(positions,dtype=float)/h
    base = np.floor(coord).astype(np.int64)
    frac = coord-base
    one_minus=1-frac
    _record(work,binary_scalings=3*len(positions),real_additions=6*len(positions),
            floor_operations=3*len(positions))
    weights = []
    for sx in (0,1):
        for sy in (0,1):
            for sz in (0,1):
                shift = np.array([sx,sy,sz])
                factors = np.where(shift,frac,one_minus)
                w = factors[:,0]*factors[:,1]*factors[:,2]
                index = tuple(((base+shift)%config.mesh_size).T)
                weights.append((index,w))
                _record(work,real_multiplications=2*len(positions),
                        integer_additions=3*len(positions),integer_remainders=3*len(positions))
    return weights


def certified_force(positions, config, multipliers,work=None, *, spectrum_capture=None):
    """CIC force with explicit deposition, Fourier and sampling roundoff."""
    n, count = config.mesh_size, len(positions)
    weights = _weights(positions, config,work)
    mass = np.zeros((n,n,n),float)
    for index,w in weights:
        np.add.at(mass,index,w)
    mean = count/n**3
    if math.frexp(mean)[0] != 0.5:
        raise ValueError("The particle-per-cell normalization must be binary")
    density = mass/mean-1
    _record(work,real_additions=8*count+n**3,binary_scalings=n**3)
    # 8 N additions is a safe bound for every accumulator.  The individual
    # positive CIC products are bounded by one; no cancellation is used.
    numerator=_mul(8*count+16,U,work)
    _record(work,real_additions=1)
    # Express the forward error in terms of the computed positive sum.
    gamma=_div(numerator,_down(1-_mul(2,numerator,work),work),work)
    density_radius=_add(_div(_mul(gamma,max(float(mass.max()),1),work),mean,work),4*U,work)
    _record(work,comparisons=mass.size)
    modes, mr = bounded_fft(density,density_radius,work=work)
    if spectrum_capture is not None:
        spectrum_capture.update(modes=modes,radius=mr)
    mode_max=_magnitude(modes,work)
    force, error, grids = [], 0.0, []
    for multiplier in multipliers:
        # The stored float frequency arrays are exact target coefficients;
        # the multiplier division contributes its own rounding enclosure.
        multmax = _magnitude(multiplier,work)
        cr = _add(_mul(multmax,mr,work),_mul(_mul(16*U,multmax,work),mode_max,work),work)
        field, fr = bounded_fft(modes*multiplier, cr, inverse=True,work=work)
        _record(work,real_multiplications=4*n**3,real_additions=2*n**3)
        field = field.real
        sample = np.zeros(count,float)
        for index,w in weights:
            sample += w*field[index]
        _record(work,real_multiplications=8*count,real_additions=8*count)
        sr=_add(fr,_mul(128*U,max(float(np.max(np.abs(field))),1),work),work)
        _record(work,absolute_values=field.size,comparisons=field.size-1)
        error = max(error,sr)
        force.append(sample)
        grids.append(field)
    return np.stack(force,axis=1), _up(error), grids


def kernel_bounds(multipliers, h,work=None):
    """CIC current-divergence bound; no factor of particle count occurs."""
    slopes, errors = [], []
    for multiplier in multipliers:
        multiplier_error = _up(16*U*float(np.max(np.abs(multiplier))))
        kernel, error = bounded_fft(multiplier,multiplier_error,inverse=True,work=work)
        kernel = kernel.real
        norm = sum(float(np.abs(kernel-np.roll(kernel,1,axis=j)).sum()) for j in range(3))/h
        norm += 6*kernel.size*error/h
        slopes.append(_up(norm*(1+32*kernel.size*U)))
        errors.append(error)
        _record(work,real_additions=6*kernel.size-1,absolute_values=3*kernel.size,
                real_divisions=2,real_multiplications=4)
    return max(slopes), max(errors)


def initial_sheet_error(initial,config,params,a0,weight):
    """Quantization error for the actual generated Lagrangian cold sheet.

    The continuous initial map is X(q)=q+s I_h D(q), with uniform q,
    and P(q)=a0^2 E(a0) f(a0) s I_h D(q).  The arrays and scale are those
    used to generate the stored particle state, without a second draw.
    """
    fields=initial["_initial_sheet"]["displacement_grid"]
    scale=abs(initial["_initial_sheet"]["displacement_scale"])
    h=config.box_size_mpch/config.mesh_size
    derivative=max(sum(float(np.max(np.abs(np.asarray(field,float)-np.roll(field,1,axis=d))))
                       for d in range(3))/h for field in fields)
    derivative=_up(scale*derivative*(1+64*U))
    amplitude=_up(scale*max(float(np.max(np.abs(field))) for field in fields))
    momentum_factor=a0*a0*expansion(a0,params)*initial["_initial_sheet"]["growth_rate"]
    # Float32 CIC sampling, scaling, and addition to the initial lattice.
    u32=2**-24
    positional_roundoff=_up(32*u32*amplitude+2*u32*config.box_size_mpch)
    momentum_roundoff=_up(32*u32*amplitude*abs(momentum_factor)+1024*U)
    half_cell=config.box_size_mpch/(2*config.particles_per_dim)
    radius=_up(max(half_cell*(1+derivative)+positional_roundoff,
                      weight*(half_cell*abs(momentum_factor)*derivative+momentum_roundoff))*(1+64*U))
    return {"initial_map":"X(q)=q+s I_h D(q), P(q)=a0^2 E(a0) f(a0) s I_h D(q), uniform q",
            "source":"the same displacement arrays, amplitude normalization and cell-centered initial lattice",
            "displacement_Lipschitz_upper":derivative,
            "initial_weighted_W_infinity_upper":radius,
            "position_roundoff_upper":positional_roundoff,
            "quantization_half_cell":half_cell,
            "same_initial_tube_admissible":False,
            "spatial_force_consistency_required":True}


def certify_pm_trajectory(config, params, *, a0=0.2, a1=0.4, tube=0.1, return_state=False):
    """Enclose the actual kick--drift--kick trajectory in max-particle norm.

    Piecewise linear canonical coordinates provide the reconstruction.  Its
    defect is bounded against the nonlinear CIC equation on each interval.
    The returned radius covers both positions and weighted canonical momenta.
    """
    if config.dynamics != "cosmological":
        raise ValueError("The trajectory certificate uses the cosmological clock")
    if config.particles_per_dim != config.mesh_size:
        raise ValueError("The local occupancy proof requires one initial particle per mesh cell")
    if not a1>a0>0 or not tube>0:
        raise ValueError("Invalid interval or continuation tube")
    if not math.isclose(1/(1+config.z_initial),a0,rel_tol=0,abs_tol=4*U*a0):
        raise ValueError("The initial state and certificate clocks must agree")
    work=ArithmeticWork()
    setup = _spectral_setup(config)
    setup["params"] = params
    setup["capture_initial_sheet"] = True
    initial = _initial_state(config,setup,np)
    x = np.asarray(initial["positions"],float)
    p = a0**3*expansion(a0,params)*np.asarray(initial["velocities"],float)
    state_digest=hashlib.sha256(x.tobytes()+p.tobytes()).hexdigest()
    operator_digest=hashlib.sha256(b"".join(np.asarray(setup[name]).tobytes()
                          for name in ("kx","ky","kz","k2"))).hexdigest()
    n = config.mesh_size
    h = config.box_size_mpch/n
    cell = np.floor(x/h).astype(int)
    if len(np.unique(cell,axis=0)) != len(x):
        raise ValueError("Initial occupancy does not satisfy the cell proof")
    centres = (cell+0.5)*h
    multipliers = _multipliers(setup,n)
    kernel_slope, kernel_roundoff = kernel_bounds(multipliers,h,work)
    omega = params.omega_b+params.omega_sigma
    alpha = lambda a: 1/(a**3*expansion(a,params))
    beta = lambda a: 1.5*omega/(a**2*expansion(a,params))
    step = (a1-a0)/config.steps_per_interval
    # Both coefficients decrease for positive matter and vacuum density.
    if omega <= 0 or params.omega_lambda < 0:
        raise ValueError("The coefficient bound requires positive matter and nonnegative vacuum")
    weight = math.sqrt(alpha(a0)/(beta(a0)*kernel_slope*1.5))
    sheet=initial_sheet_error(initial,config,params,a0,weight)
    sheet["same_initial_tube_admissible"]=sheet["initial_weighted_W_infinity_upper"]<tube
    particle_ids=initial["ids"]
    error, max_error, max_defect, max_lipschitz = 0.,0.,0.,0.
    max_force_roundoff, max_current = 0.,0.
    initial_spectrum={}
    g0, gerr0, grids0 = certified_force(x,config,multipliers,work,spectrum_capture=initial_spectrum)
    from cosmic_cg.pm_spatial_certificate import initial_force_certificate
    spatial=initial_force_certificate(initial["_initial_sheet"],config,multipliers,initial_spectrum,
                                      bands=tuple(b for b in (1,2) if b<n//2))
    if not spatial["all_checks_passed"]:
        raise ArithmeticError("Nonfinite initial spatial force enclosure")
    del initial,initial_spectrum
    production_force=_particle_acceleration(x,config,setup)
    coefficient0=coefficient_intervals(a0,params,work)
    rows = []
    for j in range(config.steps_per_interval):
        a, end = a0+j*step,a0+(j+1)*step
        dt=end-a
        middle=(a+end)/2
        _,k1=drift_kick_integrals(a,middle,params)
        drift,_=drift_kick_integrals(a,end,params)
        _,k2=drift_kick_integrals(middle,end,params)
        oldx,oldp=x.copy(),p.copy()
        # These are the production force and update operations.
        p=p+k1*production_force
        x=x+drift*p
        production_force=_particle_acceleration(x%config.box_size_mpch,config,setup)
        p=p+k2*production_force
        g1,gerr1,grids1=certified_force(x%config.box_size_mpch,config,multipliers,work)
        excursion=float(np.max(np.abs(x-centres)))
        excursion=_add(max(excursion,float(np.max(np.abs(oldx-centres)))),tube,work)
        excursion=_mul(excursion,1+8*U,work)
        if excursion>=h/2:
            raise ValueError("The continuation tube crosses the occupancy cells")
        # Each face current receives four positive transverse CIC weights.
        q=_add(1,_div(_mul(2,excursion,work),h,work),work)
        current=_mul(q,q,work)
        # Each cell contains one particle whose distance from its centre
        # is at most excursion.  The eight CIC weights give
        # (1-2e/h)^3 <= density <= (1+2e/h)^3 throughout the tube.
        # The sampling derivative is bounded by the same kernel gradient.
        density_contrast=_up(_mul(current,q,work)-1,work)
        force_lip=_mul(kernel_slope,_add(current,density_contrast,work),work)
        dx=(x-oldx)/dt
        dp=(p-oldp)/dt
        coefficient1=coefficient_intervals(end,params,work)
        alphas=(coefficient1[0][0],coefficient0[0][1])
        betas=(coefficient1[1][0],coefficient0[1][1])
        # Linear reconstruction x(theta),p(theta).  Convexity gives the
        # endpoint maximum for alpha,p separately; monotonic coefficients
        # cover the complete time interval.
        r_x=max(float(np.max(np.abs(dx-coeff*momentum)))
                for coeff in alphas for momentum in (oldp,p))
        distance=_mul(float(np.max(np.abs(x-oldx))),1+8*U,work)
        r_p=min(
            _add(max(float(np.max(np.abs(dp-coeff*g0))) for coeff in betas),
                 _mul(betas[1],_add(_mul(force_lip,distance,work),gerr0,work),work),work),
            _add(max(float(np.max(np.abs(dp-coeff*g1))) for coeff in betas),
                 _mul(betas[1],_add(_mul(force_lip,distance,work),gerr1,work),work),work))
        # Each endpoint expression uses fewer than eight rounded real
        # operations.  This absolute allowance covers cancellation in
        # coordinate differences and the two reconstructed derivatives.
        rounding=1024*U*max(1,float(np.max(np.abs(x)))/dt,
                            weight*max(float(np.max(np.abs(p))),float(np.max(np.abs(oldp))))/dt,
                            alphas[1]*max(float(np.max(np.abs(p))),float(np.max(np.abs(oldp)))),
                            weight*betas[1]*max(float(np.max(np.abs(g0))),float(np.max(np.abs(g1)))))
        defect=_add(_mul(max(r_x,_mul(weight,r_p,work)),1+1024*U,work),rounding,work)
        ell=max(_div(alphas[1],weight,work),_mul(_mul(weight,betas[1],work),force_lip,work))
        growth=exponential_upper(_mul(ell,dt,work),work)
        error=_mul(growth,_add(error,_mul(dt,defect,work),work),work)
        if error>=tube:
            raise ValueError(f"Trajectory enclosure leaves tube at step {j}: {error}")
        rows.append({"a":end,"joint_error_upper":error,"defect_upper":defect,"lipschitz_upper":ell})
        max_error,max_defect,max_lipschitz=max(max_error,error),max(max_defect,defect),max(max_lipschitz,ell)
        max_force_roundoff=max(max_force_roundoff,gerr0,gerr1)
        max_current=max(max_current,current)
        g0,gerr0,grids0=g1,gerr1,grids1
        coefficient0=coefficient1
        # The state arrays have 3 N real components.  These are the
        # reconstruction differences, eight endpoint residuals and norm
        # scans; scalar recurrence arithmetic is counted where executed.
        _record(work,real_additions=39*len(x),real_multiplications=24*len(x),
                real_divisions=6*len(x),absolute_values=54*len(x),
                comparisons=54*len(x)-18)
    result={
        "target":"finite periodic CIC particle equations with the stored Fourier multipliers",
        "state_norm":"max(max_i |x_i|_infinity, momentum_weight * max_i |p_i|_infinity)",
        "interval":[a0,a1],"particles":len(x),"mesh":[n,n,n],"steps":config.steps_per_interval,
        "momentum_weight":weight,"continuation_tube":tube,"uniform_joint_error_upper":max_error,
        "relative_position_to_cell_upper":max_error/h,"max_defect_upper":max_defect,
        "max_lipschitz_upper":max_lipschitz,"kernel_current_slope_upper":kernel_slope,
        "kernel_roundoff_upper":kernel_roundoff,"force_roundoff_upper":max_force_roundoff,
        "face_current_upper":max_current,"spatial_continuum_certified":False,
        "original_continuous_initial_sheet":sheet,
        "initial_Newtonian_force_observables":spatial,
        "initial_data":"exact stored binary particle coordinates and canonical momenta",
        "initial_state_sha256":state_digest,
        "finite_operator_data_sha256":operator_digest,
        "work":{"production_force_evaluations":config.steps_per_interval+1,
                "certificate_force_evaluations":config.steps_per_interval+1,
                "certificate_depositions":config.steps_per_interval+1,
                "certificate_3d_complex_ffts":4*(config.steps_per_interval+1)+3,
                "certificate_force_sample_components":3*(config.steps_per_interval+1)*len(x),
                "retained_trajectory_states":2,"certificate_steps":config.steps_per_interval,
                "explicit_arithmetic_counts":work.counts,
                "remaining_scalar_arithmetic_upper":2048*(config.steps_per_interval+1),
                "one_time_multiplier_arithmetic_upper":64*n**3,
                "one_time_initial_sheet_arithmetic_upper":32*n**3,
                "count_scope":"certificate arithmetic; vector indexing/copies and production integrator counted separately",
                "twiddle_preprocessing":{"interval_precision_decimal_digits":40,
                    "trigonometric_enclosures":4*(n-1),"reused_for_all_steps":True}},
        "rows":rows,"all_checks_passed":True,
    }
    if return_state:
        result["_final_state"]={"positions":x%config.box_size_mpch,
                    "velocities":p/(a1**3*expansion(a1,params)),"ids":particle_ids}
    return result
