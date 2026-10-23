"""Independent CIC algebra and numerical-envelope verification.

Fraction-operation calls and observed operand/result bit lengths are recorded,
separately from the spatial certificate's production arithmetic.
"""
import hashlib
import itertools
import json
import math
import sys
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
import mpmath as mp
from p9_cic_replacement_validation import SOURCE

from p9_cic_replacement_validation import (
    WORK, algebra_checks, exact_cell_moments, exact_sample,
    exact_replacement_bound, SIGNS, SUBSETS, original_trajectory, continuous_reconstruction_metrics,
)
from p9_continuous_force_bound import (
    spline_norm_upper, spectral_factors, interval_constant,
)
from p9_time_spatial_bounds import projected_balls
from cosmic_cg.pm_spatial_certificate import ComplexBall


class FractionLedger:
    def __init__(self):
        self.calls={};self.maximum_bits=0
    def inspect(self,value):
        if isinstance(value,Q):
            self.maximum_bits=max(self.maximum_bits,abs(value.numerator).bit_length(),value.denominator.bit_length())
        elif isinstance(value,(tuple,list)):
            for item in value:self.inspect(item)
    def profile(self,frame,event,arg):
        if frame.f_globals.get('__name__')!='fractions':return
        if event=='call':
            name=frame.f_code.co_name
            self.calls[name]=self.calls.get(name,0)+1
            for value in frame.f_locals.values():self.inspect(value)
        elif event=='return':self.inspect(arg)


def walsh_corners(points,masses):
    base=tuple(math.floor(v) for v in points[0])
    moments=[sum(m*math.prod(p[d]-base[d]-Q(1,2) for d in S)
                 for p,m in zip(points,masses)) for S in SUBSETS]
    transformed=[Q(2**len(S),8)*moments[index] for index,S in enumerate(SUBSETS)]
    additions=0
    for d in range(3):
        for index in range(8):
            if index&(1<<d):continue
            a,b=transformed[index],transformed[index|(1<<d)]
            transformed[index]=a+b;transformed[index|(1<<d)]=a-b;additions+=2
    return [transformed[sum((1-epsilon[d])<<d for d in range(3))] for epsilon in SIGNS],additions


def moment_and_crossing_checks():
    rng=np.random.default_rng(78129);walsh=0;cross=0
    for origin in ((0,0,0),(-3,2,4)):
        points=[tuple(Q(origin[d])+Q(int(v),16) for d,v in enumerate(row))
                for row in rng.integers(0,16,(11,3))]
        masses=[Q(int(v),8) for v in rng.integers(1,9,11)]
        corner,count=walsh_corners(points,masses)
        assert corner==exact_cell_moments(points,masses) and count==24
        walsh+=1
    grid=np.array([Q(int(v),16) for v in rng.integers(-16,17,64)],object).reshape(4,4,4)
    for centre in ((Q(0),)*3,(Q(-4),Q(2),Q(4)),(Q(1,3),Q(-1,4),Q(7,8))):
        for offsets in (((-1,-1,-1),(1,1,1)),((-17,9,-13),(5,-7,11)),((0,-1,2),(0,2,-1))):
            points=[tuple(centre[d]+Q(r[d],8) for d in range(3)) for r in offsets]
            masses=[Q(1,3),Q(2,3)]
            actual=abs(sum(m*exact_sample(grid,p) for p,m in zip(points,masses))-exact_sample(grid,centre))
            assert actual<=exact_replacement_bound(grid,points,masses,centre)
            cross+=1
    return {'fast_Walsh_exact_checks':walsh,'additions_per_Walsh':24,
            'additional_corner_tangent_winding_and_noncentroid_checks':cross}


def independent_mass_matrix():
    """Assemble cell basis integrals directly, independently of the stencil."""
    rng=np.random.default_rng(99);n=4;rows=[]
    for case in range(3):
        values=np.array([Q(int(v),16) for v in rng.integers(-16,17,n**3)],object).reshape(n,n,n)
        integral=Q(0)
        for cell in itertools.product(range(n),repeat=3):
            vertices=[values[tuple((cell[d]+shift[d])%n for d in range(3))] for shift in SIGNS]
            for i,eps in enumerate(SIGNS):
                for j,eta in enumerate(SIGNS):
                    integral+=vertices[i]*vertices[j]*math.prod(Q(1,3) if eps[d]==eta[d] else Q(1,6) for d in range(3))
        integral/=n**3
        averaged=values.copy()
        for d in range(3):
            averaged=(np.roll(averaged,1,d)+4*averaged+np.roll(averaged,-1,d))/6
        stencil=sum((a*b for a,b in zip(values.flat,averaged.flat)),Q(0))/n**3
        assert integral==stencil
        arr=np.asarray(values,float)[None,:,:,:]
        upper=spline_norm_upper(arr)
        with mp.workdps(80):
            exact=mp.sqrt(mp.mpf(integral.numerator)/integral.denominator)
            assert exact<=mp.mpf(upper)
        rows.append({'case':case,'exact_squared_norm':str(integral),'floating_norm_upper':upper})
    return rows


def current_matrix_checks():
    """Exact row and column sums for the CIC face-current matrix."""
    rng=np.random.default_rng(28);n=4;e=Q(1,8)
    base=list(itertools.product(range(n),repeat=3))
    points=[tuple(Q(c[d])+Q(1,2)+Q(int(v),32) for d,v in enumerate(row))
            for c,row in zip(base,rng.integers(-4,5,(n**3,3)))]
    C=(1+2*e)**2;rows=[]
    for d in range(3):
        matrix=np.full((n**3,n**3),Q(0),dtype=object)
        transverse=[v for v in range(3) if v!=d]
        for j,(c,p) in enumerate(zip(base,points)):
            frac=[p[b]-c[b] for b in range(3)]
            for bits in itertools.product((0,1),repeat=2):
                index=list(c);weight=Q(1)
                for axis,bit in zip(transverse,bits):
                    index[axis]=(index[axis]+bit)%n
                    weight*=frac[axis] if bit else 1-frac[axis]
                matrix[np.ravel_multi_index(tuple(index),(n,n,n)),j]+=weight
        row_sums=[sum(row,Q(0)) for row in matrix]
        col_sums=[sum(col,Q(0)) for col in matrix.T]
        assert max(row_sums)<=C and all(value==1 for value in col_sums)
        rows.append({'axis':d,'max_row_sum':str(max(row_sums)),
                     'column_sum':'1','analytic_row_bound':str(C)})
    return rows


def independent_shifted_spline_fourier():
    """Closed affine-cell integrals at 80 digits, including aliases/zero mode."""
    n=4;L=8.;h=L/n;values=[Q(1,8),Q(-1,4),Q(3,16),Q(1,16)]
    factors=spectral_factors(n,7);checks=0
    with mp.workdps(80):
        v=[mp.mpf(q.numerator)/q.denominator for q in values]
        for shift in (0.,h/2):
            for m in range(-7,8):
                if m==0:
                    exact=sum(v)/n;formula=exact
                else:
                    k=2*mp.pi*m/L;z=-1j*k*h
                    I0=mp.expm1(z)/z;I1=(mp.exp(z)*(z-1)+1)/(z*z)
                    exact=sum(mp.exp(-1j*k*(shift+j*h))*(v[j]*I0+(v[(j+1)%n]-v[j])*I1) for j in range(n))/n
                    discrete=sum(v[j]*mp.exp(-2j*mp.pi*m*j/n) for j in range(n))/n
                    formula=discrete*(mp.sin(mp.pi*m/n)/(mp.pi*m/n))**2*mp.exp(-1j*k*shift)
                assert abs(exact-formula)<mp.mpf('1e-70')
                checks+=1
        for mode in ((1,2,-3),(-4,1,0),(1,1,1)):
            data=[ComplexBall(complex(.003,-.007),0.),ComplexBall(complex(-.011,.004),0.),ComplexBall(complex(.006,.002),0.)]
            result=projected_balls(mode,data)
            m2=sum(m*m for m in mode)
            for b,ball in enumerate(result):
                exact=sum(mp.mpf(mode[b]*mode[d])/m2*mp.mpc(data[d].centre) for d in range(3))
                assert abs(exact-mp.mpc(ball.centre))<=mp.mpf(ball.radius)
    return {'affine_cell_Fourier_checks':checks,'projection_mode_checks':3,
            'covers':'node and half-cell origins; zero mode; aliases beyond Nyquist'}


def current_time_record_checks():
    cfg,params,times,x,p,g,D,s,target,regenerated=original_trajectory()
    assert not regenerated
    result=json.loads((WORK/'p9_time_spatial_bounds_20261004.json').read_text(encoding='utf-8'))
    assert result['initial_state_sha256']==target['initial_state_sha256']
    assert result['operator_sha256']==target['finite_operator_data_sha256']
    assert len(result['intervals'])==64 and len(result['anchors'])==9
    geometry=continuous_reconstruction_metrics(times,x,p,D,s,cfg,target)
    assert geometry['uniform_reconstruction_Lip_Euclidean_upper']==result['uniform_continuous_reconstruction_Lip_upper']
    assert all(row['cell_occupancy']==4096 for row in geometry['rows'])
    with mp.workdps(80):
        for j,row in enumerate(result['intervals']):
            assert row['a_interval']==[float(times[j]),float(times[j+1])]
            assert times[j+1]>times[j]
            anchor=result['anchors'][str(row['anchor_index'])]
            d=mp.mpf(row['particle_RMS_distance_upper'])
            exact_hi=mp.mpf(anchor['error_interval'][1])+mp.mpf(result['combined_Lipschitz_upper'])*d
            exact_lo=mp.mpf(anchor['Newton_force_norm_lower'])-mp.mpf(result['Newton_L2_source_Lipschitz_upper'])*d
            assert mp.mpf(row['all_frequency_reconstruction_force_error_upper'])>=exact_hi
            assert mp.mpf(row['Newton_force_norm_lower'])<=max(0,exact_lo)
    return {'continuous_time_segments_verified':64,'re_evolved_particle_steps':0,
            'cache_sha256':hashlib.sha256((WORK/'p9_actual_trajectory_20261004.npz').read_bytes()).hexdigest(),
            'initial_state_sha256':target['initial_state_sha256'],
            'operator_sha256':target['finite_operator_data_sha256']}


def main():
    ledger=FractionLedger();sys.setprofile(ledger.profile)
    try:
        rational={'original_independent_checks':algebra_checks(),
                  'additional_moments_and_crossings':moment_and_crossing_checks(),
                  'mass_matrix':independent_mass_matrix(),
                  'face_current':current_matrix_checks()}
    finally:sys.setprofile(None)
    result={'scope':'independent finite algebra, whole-frequency interpolation and interval arithmetic checks',
        'rational_checks':rational,'Fourier_checks':independent_shifted_spline_fourier(),
        'real_trajectory_record_checks':current_time_record_checks(),
        'verification_work':{'Fraction_function_calls':ledger.calls,
            'observed_Fraction_operand_and_result_bits_max':ledger.maximum_bits,
            'binary_Fraction_arithmetic_calls':sum(ledger.calls.get(name,0) for name in ('_add','_sub','_mul','_div')),
            'unreduced_binary_integer_bits_upper':2*ledger.maximum_bits+1,
            'binary_arithmetic_gcd_calls_upper':2*sum(ledger.calls.get(name,0) for name in ('_add','_sub','_mul','_div')),
            'scope':'one-time verification cost, separate from spatial-certificate production counts; Python control and library internals not a claimed exact bit total'},
        'source_files_sha256':{name:hashlib.sha256((SOURCE/name).read_bytes()).hexdigest() for name in (
            'p9_cic_replacement_validation.py','p9_continuous_force_bound.py',
            'p9_time_spatial_bounds.py','p9_i17_independent_verification.py')},
        'spatial_result_sha256':hashlib.sha256((WORK/'p9_time_spatial_bounds_20261004.json').read_bytes()).hexdigest(),
        'all_checks_passed':True}
    output=WORK/'p9_i17_independent_verification_20261004.json'
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'output':str(output),'all_checks_passed':True,
                     'Fraction_bits_max':ledger.maximum_bits},indent=2))


if __name__=='__main__':main()
