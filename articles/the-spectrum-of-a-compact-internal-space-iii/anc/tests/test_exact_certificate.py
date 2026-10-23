"""Independent polynomial identities and complete rounding/basis transitions."""
from fractions import Fraction as Q
import math
from pathlib import Path
import random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'code'))

from cosmic_cg.active_spectral import matter_rhs, _certificate_step
from cosmic_cg.exact_certificate import ExactStepCertificate, compile_polynomial, rhs_numerators


def integer_rhs(state, modes, outputs):
    ratios=[float(x).as_integer_ratio() for x in state]
    D=max(d for _,d in ratios)
    y=[n*(D//d) for n,d in ratios]
    polynomial=compile_polynomial(modes,outputs)
    numerators,_=rhs_numerators(y,D,polynomial)
    return [Q(n,polynomial.coefficient_denominator*D*D) for n in numerators]


def test_two_mode_convolution_has_distinct_density_and_symmetric_velocity_products():
    a,c,b,d=map(Q,[1,-2,3,-4])
    state=list(map(float,(a,c,b,d)))
    expected=[-b-a*d/2+c*b/2, -2*d-a*b, -3*(a*d+c*b)/2, -2*c*d,
              -b/2-3*a/2+b*d/2, -d/2-3*c/4-b*b/2, -3*b*d/2, -d*d]
    assert integer_rhs(state,2,4)==expected


def test_exact_output_support_for_signed_states_and_odd_denominators():
    rng=random.Random(8412)
    for modes in (1,3,5):
        state=[rng.uniform(-.03,.03) for _ in range(2*modes)]
        expected,_=matter_rhs([Q.from_float(x) for x in state],modes,2*modes)
        assert integer_rhs(state,modes,2*modes)==expected


def test_step_certificate_retains_subnormals_and_outward_rounding():
    tiny=math.ulp(0.)
    for state in ([tiny,-tiny,0.,tiny], [1e-320,-1e-315,2.**-500,-2.**-450],
                  [.011,-.003,-.012,.004]):
        candidate=[math.nextafter(x,math.inf) for x in state]
        h,radius=Q(1,10000),Q(3,100)
        exact=ExactStepCertificate(8,h,radius)
        upper,integrated,tube,pairs=exact.evaluate(state,candidate,2,1e-12)
        reference=_certificate_step(state,candidate,2,8,h,1e-12,radius)
        assert (upper,float(integrated/h),tube,pairs)==reference
        assert exact.counters['maximum_binary_denominator_bits'] >= max(
            x.as_integer_ratio()[1].bit_length() for x in state)
        assert exact.counters['maximum_direct_integer_bits_upper'] >= integrated.denominator.bit_length()
        assert exact.counters['maximum_scalar_integer_bits_upper'] >= Q.from_float(upper).denominator.bit_length()


def test_rejected_candidate_and_promotion_do_not_corrupt_cached_state():
    h,radius=Q(1,10000),Q(3,100)
    exact=ExactStepCertificate(8,h,radius)
    state,candidate=[.01,-.01],[.010001,-.010001]
    exact.evaluate(state,candidate,1,1e-8)
    # Reject the candidate and add a zero inside each coordinate block.
    promoted=[state[0],0.,state[1],0.]
    trial=[candidate[0],1e-8,candidate[1],-2e-8]
    upper,defect,tube,pairs=exact.evaluate(promoted,trial,2,1e-8)
    assert (upper,float(defect/h),tube,pairs)==_certificate_step(promoted,trial,2,8,h,1e-8,radius)
    # Accept this trial. Only equality with the actual next left state permits reuse.
    following=[x+1e-9 for x in trial]
    result=exact.evaluate(trial,following,2,upper)
    expected=_certificate_step(trial,following,2,8,h,upper,radius)
    assert (result[0],float(result[1]/h),result[2],result[3])==expected
    assert exact.counters['reused_coordinate_conversions']>=6
    assert exact.counters['plan_builds']==2


def test_work_ledger_counts_zero_coefficient_products_and_cold_construction():
    exact=ExactStepCertificate(2,Q(1,10000),Q(3,100))
    exact.evaluate([0.,0.],[0.,0.],1,0.)
    polynomial=compile_polynomial(1,2)
    counts=exact.counters.copy()
    assert counts['integer_quadratic_products']==0
    # Zero RHS data skip coordinate products, not coefficient products/row factors.
    assert counts['integer_multiplications']==4+polynomial.linear_terms+polynomial.unique_quadratic_terms+12+3
    assert counts['integer_additions']==6+polynomial.linear_terms+polynomial.unique_quadratic_terms+4
    assert counts['integer_exact_divisions']==4
    assert counts['scalar_fraction_operations']==5
    assert counts['scalar_gcd_calls_upper']>=13
    assert counts['build_fraction_operations']>0
    assert counts['build_pair_slots']==6
    exact.evaluate([0.,0.],[0.,0.],1,0.)
    assert exact.counters['build_fraction_operations']==counts['build_fraction_operations']
    assert exact.counters['integer_multiplications']==2*counts['integer_multiplications']
