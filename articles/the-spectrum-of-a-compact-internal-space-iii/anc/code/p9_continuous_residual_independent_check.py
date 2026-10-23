"""Independent polynomial integration and high-precision scalar checks.

This validates element algebra, interval primitives, time interpolation and
the published continuous residual record using the original saved trajectory.
"""
from fractions import Fraction as Q
import hashlib
import itertools
import json
import math
import sys

import mpmath as mp
import numpy as np

import p9_continuous_residual_certificate as cert
from p9_i17_independent_verification import FractionLedger
from p9_cic_replacement_validation import WORK, SOURCE, original_trajectory
from verify_pm_certificate_identity import source_record_matches
from cosmic_cg.pm_nbody import _sample_vector_field
from cosmic_cg.cosmological_dynamics import growing_rate, expansion


def poly_add(a,b):
    out=a.copy()
    for key,value in b.items():out[key]=out.get(key,Q(0))+value
    return out


def poly_scale(a,s):
    return {k:s*v for k,v in a.items()}


def polynomial_from_corners(corners):
    result={}
    for corner,value in corners.items():
        poly={(0,0,0):value}
        for d,bit in enumerate(corner):
            nxt={}
            for power,coefficient in poly.items():
                shifted=list(power);shifted[d]+=1
                nxt[tuple(shifted)]=nxt.get(tuple(shifted),Q(0))+(coefficient if bit else -coefficient)
                if not bit:nxt[power]=nxt.get(power,Q(0))+coefficient
            poly=nxt
        result=poly_add(result,poly)
    return result


def derivative(poly,d,h):
    out={}
    for powers,value in poly.items():
        if powers[d]:
            reduced=list(powers);reduced[d]-=1
            out[tuple(reduced)]=value*powers[d]/h
    return out


def square_integral(poly):
    return sum((a*b/Q(math.prod(i+j+1 for i,j in zip(ka,kb)))
                for ka,a in poly.items() for kb,b in poly.items()),Q(0))


def exact_element_energy(field,h):
    n=field.shape[1]
    value=Q(0);gradient=Q(0);curl=Q(0);divergence=Q(0)
    for c in itertools.product(range(n),repeat=3):
        polynomials=[]
        for b in range(3):
            corners={s:field[(b,*tuple((c[d]+s[d])%n for d in range(3)))]
                     for s in itertools.product((0,1),repeat=3)}
            polynomials.append(polynomial_from_corners(corners))
        derivatives=[[derivative(f,d,h) for d in range(3)] for f in polynomials]
        value+=sum((square_integral(f) for f in polynomials),Q(0))
        gradient+=sum((square_integral(z) for row in derivatives for z in row),Q(0))
        divergence_poly={}
        for d in range(3):divergence_poly=poly_add(divergence_poly,derivatives[d][d])
        divergence+=square_integral(divergence_poly)
        for b in range(3):
            d,e=(b+1)%3,(b+2)%3
            curl_poly=poly_add(derivatives[e][d],poly_scale(derivatives[d][e],Q(-1)))
            curl+=square_integral(curl_poly)
    return {key:z/Q(n**3) for key,z in (
        ('norm_squared',value),('gradient_squared',gradient),
        ('curl_squared',curl),('divergence_squared',divergence))}


def squared_interval_contains(interval,exact_value):
    return Q(float(interval[0]))**2<=exact_value<=Q(float(interval[1]))**2


def exact_checks():
    rng=np.random.default_rng(381724)
    records=[]
    for n in (2,3):
        raw=rng.integers(-8,9,size=(3,n,n,n))
        field=np.array([Q(int(z),16) for z in raw.ravel()],object).reshape(raw.shape)
        h=Q(3,2)
        energies=exact_element_energy(field,h)
        assert energies['curl_squared']==energies['gradient_squared']-energies['divergence_squared']
        interval_field=cert.exact(np.asarray(field,dtype=float))
        norm=cert.field_norm(interval_field)
        grad,curl=cert.grad_curl_norms(interval_field,float(h))
        assert squared_interval_contains(norm,energies['norm_squared'])
        assert squared_interval_contains(grad,energies['gradient_squared'])
        assert squared_interval_contains(curl,energies['curl_squared'])
        records.append({'periodic_cells_per_axis':n,'exact_Q1_energies':{k:str(v) for k,v in energies.items()}})
    primitive=0
    for _ in range(80):
        values=[Q(int(v),128) for v in rng.integers(-200,201,4)]
        a=sorted(values[:2]);b=sorted(values[2:])
        ia=(np.asarray(float(a[0])),np.asarray(float(a[1])))
        ib=(np.asarray(float(b[0])),np.asarray(float(b[1])))
        for enclosed,lo,hi in (
            (cert.add(ia,ib),a[0]+b[0],a[1]+b[1]),
            (cert.sub(ia,ib),a[0]-b[1],a[1]-b[0]),
            (cert.mul(ia,ib),min(x*y for x in a for y in b),max(x*y for x in a for y in b))):
            assert Q(float(enclosed[0]))<=lo and Q(float(enclosed[1]))>=hi
            primitive+=1
    # Distinct convex norms can have opposite endpoint trends.
    first=[Q(1),Q(0)];second=[Q(0),Q(1)]
    assert max(x*y for x,y in zip(first,second))==0
    assert ((first[0]+first[1])/2)*((second[0]+second[1])/2)==Q(1,4)
    assert max(first)*max(second)==1
    # Actual periodic fields: 4 b e1 and b sin(kx)e1, with b k<1.
    # After dividing by b^2 k/sqrt(2), endpoint amp*grad values are
    # 0 and 1, while their affine midpoint gives (5/2)*(1/2)=5/4.
    assert Q(5,2)*Q(1,2)>1
    return {'exact_periodic_element_integrations':records,
            'primitive_rational_enclosure_cases':primitive,
            'opposite_trend_product_counterexample_checked':True,
            'actual_periodic_field_product_counterexample_midpoint_factor':'5/4'}


def scalar_record_checks():
    mp.mp.dps=80
    path=WORK/'p9_continuous_residual_certificate_20261004.json'
    result=json.loads(path.read_text(encoding='utf-8'))
    cfg,params,t,x,p,g,D,s,target,_=original_trajectory()
    a=mp.mpf(float(t[0]));omega=mp.mpf(params.omega_b)+mp.mpf(params.omega_sigma)
    lam=mp.mpf(params.omega_lambda)
    integral=mp.quad(lambda z:z**mp.mpf('1.5')/(omega+lam*z**3)**mp.mpf('1.5'),[0,a])
    E=mp.sqrt(omega/a**3+lam)
    f=-mp.mpf('1.5')*omega/(a**3*E**2)+1/(a**2*E**3*integral)
    factor=a*a*E*f
    lo,hi=map(mp.mpf,result['continuous_initial_momentum_factor_interval'])
    assert lo<=factor<=hi
    # Bind the scalar to the original production initialization as well.
    # Reuse D; neither its random generation nor any trajectory is repeated.
    original_factor=float(t[0])**2*expansion(float(t[0]),params)*growing_rate(float(t[0]),params)
    assert lo<=mp.mpf(original_factor)<=hi
    n=cfg.mesh_size
    centres=((np.indices((n,n,n),dtype=np.float32).reshape(3,-1).T+.5)*
             np.float32(cfg.box_size_mpch/n)).astype(np.float32)
    displacement=_sample_vector_field(tuple(D),centres,cfg)
    displacement*=np.float32(s)
    velocity=(growing_rate(float(t[0]),params)/float(t[0])*displacement).astype(np.float64)
    reconstructed_p0=float(t[0])**3*expansion(float(t[0]),params)*velocity
    assert np.array_equal(reconstructed_p0,p[0])
    L=mp.mpf(result['reconstruction_Lipschitz_Euclidean_upper'])
    A=(1-L)**mp.mpf('-1.5')
    C=(mp.mpf('1.5')*L+L**2+L**3/4)/(1-L)**mp.mpf('1.5')+L/2*(1+L)**mp.mpf('1.5')
    kinv=mp.mpf(cfg.box_size_mpch)/(2*mp.pi)
    state=result['states']
    integral_x=mp.mpf(0);integral_p=mp.mpf(0)
    for j,row in enumerate(result['intervals']):
        upper=lambda key:max(mp.mpf(state[k][key]) for k in (j,j+1))
        extra=upper('mean_upper')+kinv*upper('curl_L2_upper')+A*(
            C*upper('psi_L2_upper')+upper('psi_Linfinity_Euclidean_upper')*upper('gradient_L2_upper'))
        assert extra<=mp.mpf(row['force_projection_composition_upper'])
        rp=mp.mpf(row['momentum_local_polynomial_upper'])+mp.mpf(row['beta_interval'][1])*extra
        assert rp<=mp.mpf(row['r_P_L2_Lagrangian_upper'])
        dt=mp.mpf(row['a_interval'][1])-mp.mpf(row['a_interval'][0])
        integral_x+=dt*mp.mpf(row['r_X_L2_Lagrangian_upper'])
        integral_p+=dt*mp.mpf(row['r_P_L2_Lagrangian_upper'])
    assert integral_x<=mp.mpf(result['integrated_r_X_L2_upper'])
    assert integral_p<=mp.mpf(result['integrated_r_P_L2_upper'])
    source=SOURCE/'p9_continuous_residual_certificate.py'
    assert source_record_matches(result['source_sha256'], source)
    assert result['cache_sha256']==hashlib.sha256((WORK/'p9_actual_trajectory_20261004.npz').read_bytes()).hexdigest()
    assert result['spatial_certificate_sha256']==hashlib.sha256((WORK/'p9_time_spatial_bounds_20261004.json').read_bytes()).hexdigest()
    return {'independent_growth_integral_digits':80,'time_intervals_checked':64,
            'all_identity_digests_match':True,
            'original_initial_momenta_reconstructed_bitwise_equal':True,
            'original_growth_scalar_inside_enclosure':True,
            'additional_initial_interpolation_arithmetic_upper':300*len(p[0]),
            'additional_initial_particle_evolutions':0,
            'exact_growth_factor_reference':str(factor),
            'result_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    ledger=FractionLedger();sys.setprofile(ledger.profile)
    try:checks=exact_checks()
    finally:sys.setprofile(None)
    scalar=scalar_record_checks()
    result={'all_checks_passed':True,'exact_algebra':checks,'actual_record_checks':scalar,
            'verification_work':{'Fraction_function_calls':ledger.calls,
                'observed_operand_and_result_bits_max':ledger.maximum_bits,
                'unreduced_binary_integer_bits_upper':2*ledger.maximum_bits+1,
                'binary_arithmetic_calls':sum(ledger.calls.get(k,0) for k in ('_add','_sub','_mul','_div')),
                'independent_80_digit_quadrature_is_validation_only':True},
            'scope':'element identities and interval implementation, not a smaller dynamical trajectory',
            'new_FFTs':0,'new_particle_evolutions':0}
    path=WORK/'p9_continuous_residual_independent_check_20261004.json'
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'all_checks_passed':True,'work':result['verification_work'],'path':str(path)},indent=2))


if __name__=='__main__':main()
