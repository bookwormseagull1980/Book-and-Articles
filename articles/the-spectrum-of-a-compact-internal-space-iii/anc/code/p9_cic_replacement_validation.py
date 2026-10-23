"""Exact CIC identities and the saved 4096-particle trajectory.

The default path requires the published cache. Algebra checks use exact
fractions independently of the production deposition routine.
"""
from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys

import numpy as np

SOURCE = Path(__file__).resolve().parent
P9 = SOURCE.parents[1]
WORK = P9 / "anc" / "results" / "pm_continuum_current"
sys.path.insert(0, str(P9 / "anc" / "code"))
from cosmic_cg.params import load_model_params
from cosmic_cg.pm_nbody import (PMNBodyConfig, _spectral_setup, _initial_state,
    _particle_acceleration)
from cosmic_cg.cosmological_dynamics import expansion, drift_kick_integrals
from cosmic_cg.pm_certificate import (ArithmeticWork, bounded_fft,
    certified_force, _multipliers, _weights, U)

SUBSETS = tuple(tuple(d for d in range(3) if mask & (1 << d)) for mask in range(8))
SIGNS = tuple(itertools.product((0, 1), repeat=3))


def up(value):
    return float(np.nextafter(float(value), math.inf))


def positive_sum(values):
    # Directed at every addition; no libm/compensated-sum accuracy assumption.
    total=0.
    for value in values:
        if value<0:raise ValueError('The positive sum bound requires nonnegative terms')
        total=up(total+float(value))
    return total


def exact_weights(point, n):
    """Exact periodic CIC weights in mesh units, arbitrary signed positions."""
    base = tuple(math.floor(x) for x in point)
    frac = tuple(point[d] - base[d] for d in range(3))
    return [(tuple((base[d] + shift[d]) % n for d in range(3)),
             math.prod(frac[d] if shift[d] else 1-frac[d] for d in range(3)))
            for shift in SIGNS]


def exact_sample(grid, point):
    return sum(w * grid[index] for index, w in exact_weights(point, grid.shape[0]))


def exact_cell_moments(points, masses):
    base = tuple(math.floor(x) for x in points[0])
    assert all(tuple(math.floor(x) for x in p) == base for p in points)
    moments = []
    for subset in SUBSETS:
        moments.append(sum(m * math.prod(p[d]-base[d]-Q(1, 2) for d in subset)
                           for p, m in zip(points, masses)))
    corners = []
    for shift in SIGNS:
        corners.append(sum(Q(1, 2**(3-len(subset))) *
            math.prod(2*shift[d]-1 for d in subset) * moments[s]
            for s, subset in enumerate(SUBSETS)))
    return corners


def crossing_weight(q, displacement, h=1):
    """Sum 1-t over grid-plane crossings, including a departing t=0 face.

    The derivative at q is chosen from the half-open floor(q/h) cell.
    A path directed left from an exact grid face has its first event at zero.
    Endpoint events have zero weight and can be included.
    """
    if displacement == 0:
        return q * 0
    end = q + displacement
    if displacement > 0:
        indices = range(math.floor(q/h)+1, math.floor(end/h)+1)
    else:
        indices = range(math.floor(q/h), math.ceil(end/h)-1, -1)
    return sum((1-(j*h-q)/displacement for j in indices), q*0)


def exact_grid_constants(grid):
    first = [max(abs(x) for x in (grid-np.roll(grid, 1, d)).flat) for d in range(3)]
    mixed = {(d, e): max(abs(x) for x in
        (grid-np.roll(grid,1,d)-np.roll(grid,1,e)+np.roll(np.roll(grid,1,d),1,e)).flat)
        for d in range(3) for e in range(d+1,3)}
    jumps = [max(abs(x) for x in
        (np.roll(grid,-1,d)-2*grid+np.roll(grid,1,d)).flat) for d in range(3)]
    return first, mixed, jumps


def exact_replacement_bound(grid, points, masses, centre):
    first, mixed, jumps = exact_grid_constants(grid)
    residual = [[p[d]-centre[d] for d in range(3)] for p in points]
    bound = sum(first[d]*abs(sum(m*r[d] for m,r in zip(masses,residual)))
                for d in range(3))
    bound += sum(mixed[d,e]*sum(m*abs(r[d]*r[e]) for m,r in zip(masses,residual))
                 for d,e in mixed)
    bound += sum(jumps[d]*sum(m*abs(r[d])*crossing_weight(centre[d],r[d])
                              for m,r in zip(masses,residual)) for d in range(3))
    return bound


def algebra_checks():
    rng = np.random.default_rng(814317)
    for _ in range(30):
        points = [tuple(Q(int(v),32) for v in row) for row in rng.integers(1,32,(9,3))]
        masses = [Q(int(v),8) for v in rng.integers(1,8,9)]
        direct = [sum(m * math.prod(p[d] if s[d] else 1-p[d] for d in range(3))
                      for p,m in zip(points,masses)) for s in SIGNS]
        assert exact_cell_moments(points,masses) == direct
    for _ in range(40):
        grid = np.array([Q(int(v),8) for v in rng.integers(-16,17,64)],object).reshape(4,4,4)
        points = [tuple(Q(int(v),8) for v in row) for row in rng.integers(-12,29,(7,3))]
        masses = [Q(1,7)]*7
        centre = tuple(sum(m*p[d] for m,p in zip(masses,points)) for d in range(3))
        difference = abs(sum(m*exact_sample(grid,p) for m,p in zip(masses,points))-exact_sample(grid,centre))
        assert difference <= exact_replacement_bound(grid,points,masses,centre)
    grid=np.zeros((4,4,4),dtype=object)
    for d in range(4):
        grid[d,:,:]=Q(d%2)
    centre=(Q(0),Q(1,3),Q(2,3))
    points=[(Q(-1,16),centre[1],centre[2]),(Q(1,16),centre[1],centre[2])]
    actual=abs(sum(Q(1,2)*exact_sample(grid,p) for p in points)-exact_sample(grid,centre))
    bound=exact_replacement_bound(grid,points,[Q(1,2)]*2,centre)
    assert actual == bound == Q(1,16)
    return {"cell_moment_identities_exact":30,"arbitrary_grid_crossing_tests_exact":40,
            "grid_face_kink_sharp_example":{"actual":str(actual),"bound":str(bound)},
            "scope":"independent exact algebra tests; not a smaller dynamical target"}


def original_trajectory(*, require_cache=True):
    params=load_model_params()
    cfg=PMNBodyConfig(dynamics="cosmological",particles_per_dim=16,mesh_size=16,
        box_size_mpch=32.,z_initial=4.,steps_per_interval=64,
        initial_displacement_fraction=.002,mode_bins=8,max_halos_per_snapshot=16)
    formal=json.loads((P9/'anc/results/cosmological_active_current/results.json').read_text(encoding='utf-8'))
    target=formal['PM']['trajectory_certificate']
    cache=WORK/'p9_actual_trajectory_20261004.npz'
    if cache.exists():
        saved=np.load(cache,allow_pickle=False)
        if str(saved['initial_sha']) != target['initial_state_sha256'] or str(saved['operator_sha']) != target['finite_operator_data_sha256']:
            raise RuntimeError('Cached trajectory target differs from formal certificate')
        if hashlib.sha256(saved['x'][0].tobytes()+saved['p'][0].tobytes()).hexdigest()!=target['initial_state_sha256']:
            raise RuntimeError('Cached initial arrays differ from their recorded digest')
        if saved['x'].shape!=(65,4096,3) or saved['p'].shape!=saved['x'].shape or saved['g'].shape!=saved['x'].shape:
            raise RuntimeError('The cache is not the specified actual trajectory')
        if not all(np.all(np.isfinite(saved[key])) for key in ('x','p','g','D','a','s')):
            raise RuntimeError('Nonfinite cached trajectory data')
        return cfg,params,saved['a'],saved['x'],saved['p'],saved['g'],saved['D'],float(saved['s']),target,False
    if require_cache:
        raise FileNotFoundError('This validation requires the saved real trajectory; regeneration is not authorized')
    setup=_spectral_setup(cfg);setup['params']=params;setup['capture_initial_sheet']=True
    initial=_initial_state(cfg,setup,np)
    x=np.asarray(initial['positions'],float)
    p=.2**3*expansion(.2,params)*np.asarray(initial['velocities'],float)
    initial_sha=hashlib.sha256(x.tobytes()+p.tobytes()).hexdigest()
    operator_sha=hashlib.sha256(b''.join(np.asarray(setup[name]).tobytes() for name in ('kx','ky','kz','k2'))).hexdigest()
    assert initial_sha==target['initial_state_sha256']
    assert operator_sha==target['finite_operator_data_sha256']
    times=.2+np.arange(65)*(.2/64)
    force=_particle_acceleration(x,cfg,setup)
    xs=[x.copy()];ps=[p.copy()];gs=[force.copy()]
    for a,end in zip(times[:-1],times[1:]):
        middle=(a+end)/2
        _,k1=drift_kick_integrals(a,middle,params)
        drift,_=drift_kick_integrals(a,end,params)
        _,k2=drift_kick_integrals(middle,end,params)
        p=p+k1*force;x=x+drift*p
        force=_particle_acceleration(x%cfg.box_size_mpch,cfg,setup)
        p=p+k2*force
        xs.append(x.copy());ps.append(p.copy());gs.append(force.copy())
    D=np.stack(initial['_initial_sheet']['displacement_grid'])
    s=initial['_initial_sheet']['displacement_scale']
    np.savez_compressed(cache,a=times,x=np.asarray(xs),p=np.asarray(ps),g=np.asarray(gs),D=D,s=s,
                        initial_sha=initial_sha,operator_sha=operator_sha)
    return cfg,params,times,np.asarray(xs),np.asarray(ps),np.asarray(gs),D,s,target,True


def actual_kernel_constants(multipliers,h):
    work=ArithmeticWork();kernels=[];radius=0.
    for multiplier in multipliers:
        mult_error=up(16*U*np.max(np.abs(multiplier)))
        arr,rad=bounded_fft(multiplier,mult_error,inverse=True,work=work)
        kernels.append(arr.real);radius=max(radius,rad)
    kernels=np.asarray(kernels)
    # Every difference has at most four terms. Explicit binary arithmetic
    # plus the common kernel ball encloses the exact stored-grid operator.
    maxk=up(float(np.max(np.abs(kernels)))+radius)
    allowance=up(32*U*max(maxk,1.))
    first=[up((float(np.max(np.abs(kernels-np.roll(kernels,1,d+1))))+2*radius+allowance)/h) for d in range(3)]
    mixed={}
    for d in range(3):
        for e in range(d+1,3):
            difference=kernels-np.roll(kernels,1,d+1)-np.roll(kernels,1,e+1)+np.roll(np.roll(kernels,1,d+1),1,e+1)
            mixed[d,e]=up((float(np.max(np.abs(difference)))+4*radius+allowance)/(h*h))
    jumps=[up((float(np.max(np.abs(np.roll(kernels,-1,d+1)-2*kernels+np.roll(kernels,1,d+1))))+4*radius+allowance)/h) for d in range(3)]
    return first,mixed,jumps,{"complex_ffts":3,"explicit_arithmetic":work.counts,
        "kernel_difference_real_additions":54*len(kernels[0].flat),
        "kernel_difference_binary_scalings":9*len(kernels[0].flat),
        "kernel_radius":radius}


def sample64(grids,x,cfg):
    return np.stack([sum(w*field[index] for index,w in _weights(x,cfg)) for field in grids],axis=1)


def actual_monopole_tests(cfg,x,multipliers,constants):
    first,mixed,jumps=constants
    n=cfg.mesh_size;h=cfg.box_size_mpch/n
    lattice=np.indices((n,n,n)).reshape(3,-1).T
    basework=ArithmeticWork()
    base_force,base_error,_=certified_force(x%cfg.box_size_mpch,cfg,multipliers,basework)
    rows=[]
    for width in (2,4):
        labels=((lattice//width)*np.array([(n//width)**2,n//width,1])).sum(axis=1)
        replacement=np.empty_like(x);bound=0.;bound_first=0.;bound_mixed=0.;bound_cross=0.
        for label in np.unique(labels):
            ids=np.flatnonzero(labels==label);points=x[ids]
            q=np.array([math.fsum(points[:,d])/len(ids) for d in range(3)])
            replacement[ids]=q
            r=points-q
            # Exact binary inputs are enclosed by a subtraction allowance;
            # crossing factors are evaluated on exact Fractions below.
            for d in range(3):
                m1=sum((Q(float(value))-Q(float(q[d])) for value in points[:,d]),Q(0))
                bound_first=up(bound_first+up(first[d]*up(float(abs(m1)))))
            for d,e in mixed:
                m2=sum((abs((Q(float(point[d]))-Q(float(q[d])))*(Q(float(point[e]))-Q(float(q[e])))) for point in points),Q(0))
                bound_mixed=up(bound_mixed+up(mixed[d,e]*up(float(m2))))
            for d in range(3):
                qc=Q(float(q[d]));hh=Q(float(h))
                event=sum((abs(Q(float(point[d]))-qc)*crossing_weight(qc,Q(float(point[d]))-qc,hh) for point in points),Q(0))
                bound_cross=up(bound_cross+up(jumps[d]*up(float(event))))
        bound=positive_sum((bound_first,bound_mixed,bound_cross))
        work=ArithmeticWork()
        _,candidate_error,grids=certified_force(replacement%cfg.box_size_mpch,cfg,multipliers,work)
        sample=sample64(grids,x%cfg.box_size_mpch,cfg)
        sample_allowance=up(128*U*max(1.,float(np.max(np.abs(grids)))))
        measured=float(np.max(np.abs(sample-base_force)))
        radius=positive_sum((base_error,candidate_error,sample_allowance,32*U*max(1.,measured)))
        lower=max(0.,float(np.nextafter(measured-radius,-math.inf)));upper=up(measured+radius)
        assert upper <= bound
        rows.append({"group_width_in_original_cells":width,"monopoles":int(len(np.unique(labels))),
            "original_particles":len(x),"grid_force_samples_at_original_positions":[lower,upper],
            "rigorous_global_bound":bound,"first_moment_term":bound_first,
            "mixed_curvature_term":bound_mixed,"grid_crossing_term":bound_cross,
            "force_work":work.counts,
            "accepted_at_formal_trajectory_tolerance":False})
    return rows,basework.counts


def grid_lipschitz_matrix(field,h):
    matrix=np.zeros((3,3))
    for b in range(3):
        absolute=up(float(np.max(np.abs(field[b]))))
        for d in range(3):
            difference=field[b]-np.roll(field[b],-1,d)
            matrix[b,d]=up((float(np.max(np.abs(difference)))+up(4*U*absolute))/h)
    return matrix


def continuous_reconstruction_metrics(times,x,p,D,s,cfg,target):
    n=cfg.mesh_size;h=cfg.box_size_mpch/n
    centres=((np.indices((n,n,n)).reshape(3,-1).T+.5)*h)
    # Xhat(q,a)=q+psi0(q)+I_centred[x(a)-x(0)](q). At a0 it
    # agrees exactly with the original continuous sheet. Both meshes have
    # spacing h; their faces differ by h/2 and share a finite overlay.
    initial=np.asarray(D,float)*s
    assert float(np.float32(s))==s
    initial_matrix=grid_lipschitz_matrix(initial,h)
    rows=[]
    for a,positions,momenta in zip(times,x,p):
        inc=(positions-x[0]).reshape(n,n,n,3).transpose(3,0,1,2)
        # The original cell-centred coordinates stay positive and differ
        # by less than a factor two: Sterbenz makes x(a)-x(0) exact.
        assert np.all((positions>=x[0]/2)&(positions<=2*x[0]))
        matrix=np.nextafter(initial_matrix+grid_lipschitz_matrix(inc,h),np.inf)
        lip=max(positive_sum(row) for row in matrix)
        col=max(positive_sum(column) for column in matrix.T)
        lip2=up(math.sqrt(up(lip*col)))
        excursion=float(np.max(np.abs(positions-centres)))
        if excursion>=h/2:
            raise ValueError('The saved reconstruction leaves its assigned occupancy cells')
        rows.append({"a":float(a),"reconstructed_displacement_Lip_infinity_upper":lip,
                     "reconstructed_displacement_Lip_Euclidean_upper":lip2,
                     "cell_occupancy":len(positions),
                     "max_particle_displacement_from_original_centre":excursion})
    # This bound concerns only the saved reconstruction; a C0 state tube
    # does not imply derivative/Jacobian control for the Newtonian solution.
    return {"construction":"q + original trilinear psi0(q) + centred-grid interpolation of particle displacement increments",
            "initial_continuous_reconstruction_error":0,
            "uniform_reconstruction_Lip_upper":max(r['reconstructed_displacement_Lip_infinity_upper'] for r in rows),
            "uniform_reconstruction_Lip_Euclidean_upper":max(r['reconstructed_displacement_Lip_Euclidean_upper'] for r in rows),
            "uniform_empirical_W_infinity_lower":cfg.box_size_mpch/(2*len(x[0])**(1/3)),
            "rows":rows,"Newtonian_evolved_trajectory_certified":False,
            "missing":["time-dependent Newtonian force residual including unretained modes",
                       "stability in a regularity space controlling density/force and actual flow Jacobian",
                       "continuous momentum reconstruction and time-integrated defect enclosure"]}


def main():
    tests=algebra_checks()
    cfg,params,times,x,p,g,D,s,target,regenerated=original_trajectory()
    setup=_spectral_setup(cfg);multipliers=_multipliers(setup,cfg.mesh_size)
    first,mixed,jumps,kwork=actual_kernel_constants(multipliers,2.)
    mono,baseline_work=actual_monopole_tests(cfg,x[-1],multipliers,(first,mixed,jumps))
    reconstruction=continuous_reconstruction_metrics(times,x,p,D,s,cfg,target)
    count=len(x[0]);occupied=[r['cell_occupancy'] for r in reconstruction['rows']]
    result={"upstream_run_id":params.source_run_id,"target_initial_sha256":target['initial_state_sha256'],
        "target_operator_sha256":target['finite_operator_data_sha256'],"interval":[.2,.4],
        "original_particles":count,"original_mesh":[16]*3,"original_steps":64,
        "regenerated_once_this_execution":regenerated,
        "trajectory_cache":"p9_actual_trajectory_20261004.npz",
        "algebra_checks":tests,"kernel_derivative_bounds":{"first":first,"mixed":{str(k):v for k,v in mixed.items()},"face_jumps":jumps},
        "same_target_monopole_validation":mono,"continuous_reconstruction":reconstruction,
        "exact_cell_aggregation":{"minimum_occupied_cells":min(occupied),"maximum_occupied_cells":max(occupied),
            "particle_reduction_for_this_target":0,"mass_and_source_force_preserved_by_8_moments":True,
            "direct_unit_mass_deposition_work_per_state":{"multiplications":16*count,"accumulation_additions":8*count},
            "cell_moment_work_per_state":{"multiplications":4*count,"moment_additions":8*count,
                "Walsh_additions":24*max(occupied),"corner_accumulation_additions":8*max(occupied),
                "binary_scalings":8*max(occupied)},
            "break_even_ignoring_grouping_counting_binary_scalings_as_one":"N_particles > (10/3) * occupied_cells",
            "advancement_closed_by_these_moments_alone":False},
        "validation_work":{"production_force_evaluations_if_cache_built":65 if regenerated else 0,
            "additional_enclosed_force_evaluations":3,"additional_complex_ffts":15,
            "kernel_work":kwork,"baseline_force_work":baseline_work,
            "scope":"one-time candidate validation; exact-Fraction shape budgets separately identifiable in source"},
        "all_checks_passed":True}
    output=WORK/'p9_cic_replacement_validation_20261004.json'
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({"output":str(output),"all_checks_passed":True,
        "max_reconstruction_Lip":reconstruction['uniform_reconstruction_Lip_upper'],
        "actual_target_monopole_errors":[r['grid_force_samples_at_original_positions'] for r in mono],
        "rigorous_bounds":[r['rigorous_global_bound'] for r in mono],
        "occupied_cells": [min(occupied),max(occupied)]},indent=2))


if __name__=='__main__':
    main()
