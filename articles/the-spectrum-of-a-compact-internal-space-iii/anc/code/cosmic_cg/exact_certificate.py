"""Exact certificates for the stored finite Galerkin trajectory.

Stored IEEE coordinates share a binary denominator D. The fixed polynomial
coefficients share a denominator C, so the complete RHS and its omitted-mode
defect can be evaluated with integers. This changes the representation of
the arithmetic, not the floating-point trajectory or the error estimate.

For a step h and displacement Delta, the original piecewise-linear defect
bound is exactly

    (E + ||Delta-h F(y)||_1 + (L*h/2)||Delta||_1) / (1-L*h).

The recurrence is rounded upwards at every step, and the continuation tube
is checked in exact rational arithmetic. All nonlinear output modes up to
min(2*M,N) are included. No sampled reference solution is needed.

The work ledger separates construction, direct integer arithmetic, and scalar
Fraction arithmetic. Counts concern the written arithmetic operations, including
zero coefficient products, rather than machine instructions. Loop/index/cache
control, storage and ledger maintenance are excluded. GCD counts and integer
bit lengths are conservative bounds, not instrumented operation counts. With
B the relevant phase's reported bit bound, a schoolbook implementation admits
O(A*B + (M+D+F+K)*B**2 + G*B**3) bit work: A integer additions/subtractions,
sign/absolute-value operations and dyadic conversions; M multiplications,
D divisions, F Fraction operations/comparisons/conversions,
K LCM calls and G the GCD-call upper bound. This deliberately conservative
symbolic bound is not a wall-clock model and is not added to floating work.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
import math


@dataclass(frozen=True)
class GalerkinPolynomial:
    modes: int
    outputs: int
    coefficient_denominator: int
    linear: tuple
    quadratic: tuple
    logical_pairs: int
    unique_quadratic_terms: int
    linear_terms: int
    coefficient_bits: int
    construction_work: dict


def compile_polynomial(modes, outputs):
    """Collect the exact density/velocity polynomial on the stated support.

    Velocity monomials b_j*b_k are commutative. Density coordinates a_j and
    velocity coordinates b_j have disjoint indices, so a_j*b_k is never
    confused with a_k*b_j. Only exactly equal monomials are collected.
    """
    if type(modes) is not int or type(outputs) is not int or min(modes, outputs) < 1:
        raise ValueError("positive integer mode counts required")
    linear = [dict() for _ in range(2*outputs)]
    quadratic = [dict() for _ in range(2*outputs)]
    pairs = 0
    def add(row, key, coefficient):
        row[key] = row.get(key, Q(0)) + coefficient
    for n in range(1, outputs+1):
        if n <= modes:
            add(linear[n-1], modes+n-1, Q(-n))
            add(linear[outputs+n-1], modes+n-1, Q(-1, 2))
            add(linear[outputs+n-1], n-1, Q(-3, 2*n))
        for p in range(-modes, modes+1):
            q = n-p
            if p == 0 or q == 0 or abs(q) > modes:
                continue
            sp, sq = (1 if p > 0 else -1), (1 if q > 0 else -1)
            a, b, other_b = abs(p)-1, modes+abs(p)-1, modes+abs(q)-1
            add(quadratic[n-1], (a, other_b), Q(-n*sq, 2))
            add(quadratic[outputs+n-1], tuple(sorted((b, other_b))), Q(-n*sp*sq, 4))
            pairs += 1
    denominator = math.lcm(*(c.denominator for rows in (linear, quadratic)
                             for row in rows for c in row.values()))
    def integer_rows(rows):
        return tuple(tuple((index, int(c*denominator)) for index, c in sorted(row.items()) if c)
                     for row in rows)
    linear_int, quadratic_int = integer_rows(linear), integer_rows(quadratic)
    linear_terms = sum(map(len, linear_int))
    quadratic_terms = sum(map(len, quadratic_int))
    entries = sum(len(row) for rows in (linear, quadratic) for row in rows)
    additions = 3*min(modes, outputs)+2*pairs
    terms = linear_terms+quadratic_terms
    coefficient_bits = max(abs(c).bit_length() for rows in (linear_int, quadratic_int)
                           for row in rows for _, c in row)
    # Q(0), the default argument of row.get, is constructed on every add.
    # Each Fraction binary operation needs at most two GCDs; each explicit
    # constructor at most one. A variadic LCM needs at most one per entry.
    work = dict(build_pair_slots=outputs*(2*modes+1),
        build_fraction_operations=additions+terms,
        build_explicit_fraction_constructors=2*additions,
        build_lcm_calls=1, build_gcd_calls_upper=4*additions+2*terms+entries,
        build_integer_multiplications_upper=3*pairs+min(modes, outputs)+entries,
        build_integer_divisions_upper=terms+entries,
        build_integer_sign_changes=2*pairs+min(modes, outputs),
        build_max_integer_bits_upper=4*(denominator.bit_length()+
            max(modes, outputs).bit_length()+additions.bit_length()+2))
    return GalerkinPolynomial(modes, outputs, denominator, linear_int,
        quadratic_int, pairs, quadratic_terms, linear_terms, coefficient_bits, work)


def _ratios(values):
    key = tuple(float(x) for x in values)
    if not key or len(key) % 2 or not all(math.isfinite(x) for x in key):
        raise ValueError("finite density/velocity coordinate blocks required")
    return key, tuple(x.as_integer_ratio() for x in key)


def _embed(values, old_modes, new_modes, zero):
    return (tuple(values[:old_modes]) + (zero,)*(new_modes-old_modes) +
            tuple(values[old_modes:]) + (zero,)*(new_modes-old_modes))


def rhs_numerators(integers, denominator, polynomial):
    """Return exact RHS numerators over C*D^2 and actual quadratic products."""
    if len(integers) != 2*polynomial.modes or denominator < 1:
        raise ValueError("state and polynomial dimensions disagree")
    result, products = [], {}
    actual_multiplications = 0
    for linear, quadratic in zip(polynomial.linear, polynomial.quadratic):
        value = denominator*sum(c*integers[index] for index, c in linear)
        for pair, coefficient in quadratic:
            if pair not in products:
                a, b = pair
                if integers[a] == 0 or integers[b] == 0:
                    products[pair] = 0
                else:
                    products[pair] = integers[a]*integers[b]
                    actual_multiplications += 1
            value += coefficient*products[pair]
        result.append(value)
    return result, actual_multiplications


def _round_up(value):
    rounded = float(value)
    rounded_q = Q.from_float(rounded)
    conversions = 1
    if rounded_q < value:
        rounded = math.nextafter(rounded, math.inf)
        rounded_q = Q.from_float(rounded)
        conversions += 1
    return rounded, rounded_q, conversions


class ExactStepCertificate:
    """Per-run exact arithmetic with state-keyed caches and cold construction.

    A rejected candidate is reused only if its coordinates really become a
    subsequent left state. Zero-padding promotions preserve both coordinate
    blocks separately. The cache contains no reference trajectory or fitted
    defect estimate. Every new retained support builds its own polynomial
    inside the first evaluation; there is no shared warm-run cache.
    """
    def __init__(self, full_modes, step, radius):
        if type(full_modes) is not int or full_modes < 1 or not isinstance(step, Q) or step <= 0:
            raise ValueError("integer full carrier and positive rational step required")
        if not isinstance(radius, Q) or radius <= 0:
            raise ValueError("positive rational radius required")
        self.full_modes, self.step, self.radius = full_modes, step, radius
        self.lip = Q(full_modes) + Q(1, 2) + 2*full_modes*radius
        if self.lip*step >= 1:
            raise ValueError("step exceeds the rational amplification bound")
        self.amplification = 1/(1-self.lip*step)
        self.linear_remainder_factor = self.lip*step/2
        self._left = self._right = None
        self._error_key = self._error_value = None
        self._polynomials = {}
        self._constant_bits = max(abs(x.numerator).bit_length() for x in
            (self.amplification, self.linear_remainder_factor, radius, step))
        self._constant_bits = max(self._constant_bits, *(x.denominator.bit_length()
            for x in (self.amplification, self.linear_remainder_factor, radius, step)))
        self._step_bits = (abs(step.numerator).bit_length(), step.denominator.bit_length())
        setup_bits = max(full_modes.bit_length(), abs(step.numerator).bit_length(),
            step.denominator.bit_length(), abs(radius.numerator).bit_length(), radius.denominator.bit_length())
        self.counters = dict(plan_builds=0, polynomial_quadratic_terms=0,
            dyadic_coordinate_conversions=0, reused_coordinate_conversions=0,
            integer_quadratic_products=0, certificate_calls=0, logical_pairs=0,
            setup_fraction_operations=9, setup_fraction_comparisons=3,
            setup_explicit_fraction_constructors=2, setup_integer_multiplications=1,
            setup_gcd_calls_upper=20, setup_max_integer_bits_upper=8*setup_bits+16,
            build_pair_slots=0, build_fraction_operations=0,
            build_explicit_fraction_constructors=0, build_lcm_calls=0,
            build_gcd_calls_upper=0, build_integer_multiplications_upper=0,
            build_integer_divisions_upper=0, build_integer_sign_changes=0,
            build_max_integer_bits_upper=0,
            integer_additions=0, integer_subtractions=0, integer_multiplications=0,
            integer_exact_divisions=0, integer_absolute_values=0,
            scalar_fraction_operations=0, scalar_fraction_comparisons=0,
            scalar_explicit_fraction_constructors=0, scalar_dyadic_conversions=0,
            scalar_fraction_to_float_conversions=0, scalar_gcd_calls_upper=0,
            maximum_binary_denominator_bits=0, maximum_coordinate_numerator_bits_upper=0,
            maximum_direct_integer_bits_upper=0, maximum_scalar_integer_bits_upper=0)

    def _prepare_left(self, values):
        key = tuple(float(x) for x in values)
        modes = len(key)//2
        for cached in (self._right, self._left):
            if cached is None:
                continue
            old_key, ratios = cached
            if key == old_key:
                self.counters['reused_coordinate_conversions'] += len(key)
                return cached
            old_modes = len(old_key)//2
            if old_modes < modes and _embed(old_key, old_modes, modes, 0.0) == key:
                self.counters['reused_coordinate_conversions'] += len(old_key)
                return key, _embed(ratios, old_modes, modes, (0, 1))
        prepared = _ratios(key)
        self.counters['dyadic_coordinate_conversions'] += len(key)
        return prepared

    def evaluate(self, state, candidate, modes, error):
        """Return upper error, integrated defect, exact tube result, pair count.

        The last count is the number of convolution-pair contributions in the
        original polynomial. The actual integer quadratic multiplications are
        recorded separately in counters; these are different work measures.
        """
        if type(modes) is not int or not 1 <= modes <= self.full_modes:
            raise ValueError("active modes must lie in the finite reference carrier")
        if len(state) != 2*modes or len(candidate) != 2*modes:
            raise ValueError("state dimension does not match the active basis")
        if not isinstance(error, float) or not math.isfinite(error) or error < 0:
            raise ValueError("nonnegative finite stored IEEE error bound required")
        left, right = self._prepare_left(state), _ratios(candidate)
        self.counters['dyadic_coordinate_conversions'] += len(candidate)
        D = max(d for _, d in (*left[1], *right[1]))
        Y = [n*(D//d) for n, d in left[1]]
        next_Y = [n*(D//d) for n, d in right[1]]
        delta = [b-a for a, b in zip(Y, next_Y)]
        norm_numerator = max(sum(abs(x) for x in Y), sum(abs(x) for x in next_Y))
        outputs = min(2*modes, self.full_modes)
        if modes not in self._polynomials:
            self._polynomials[modes] = compile_polynomial(modes, outputs)
            self.counters['plan_builds'] += 1
            self.counters['polynomial_quadratic_terms'] += self._polynomials[modes].unique_quadratic_terms
            for key, value in self._polynomials[modes].construction_work.items():
                self.counters[key] = (max(self.counters[key], value) if key.endswith('bits_upper')
                                      else self.counters[key]+value)
        polynomial = self._polynomials[modes]
        numerators, multiplications = rhs_numerators(Y, D, polynomial)
        self.counters['integer_quadratic_products'] += multiplications
        self.counters['certificate_calls'] += 1
        self.counters['logical_pairs'] += polynomial.logical_pairs
        embedded_delta = _embed(delta, modes, outputs, 0)
        a, b = self.step.numerator, self.step.denominator
        C = polynomial.coefficient_denominator
        integral_scale = b*C*D
        defect_numerator = sum(abs(integral_scale*change-a*F)
                               for change, F in zip(embedded_delta, numerators))
        integrated_defect = Q(defect_numerator, integral_scale*D)
        displacement = Q(sum(abs(x) for x in delta), D)
        error_conversion = int(error != self._error_key)
        if error_conversion:
            self._error_key, self._error_value = error, Q.from_float(error)
        bound = self.amplification*(self._error_value + integrated_defect +
                                    self.linear_remainder_factor*displacement)
        norm = Q(norm_numerator, D)
        tube = norm+bound < self.radius
        upper, rounded_q, rounding_conversions = _round_up(bound)
        # Fixed-path accounting includes products with zero and the D factor
        # for every RHS row. The common b*C*D factor is computed only once.
        linear_terms, quadratic_terms = polynomial.linear_terms, polynomial.unique_quadratic_terms
        self.counters['integer_additions'] += 6*modes+linear_terms+quadratic_terms+2*outputs
        self.counters['integer_subtractions'] += 2*modes+2*outputs
        self.counters['integer_multiplications'] += 4*modes+linear_terms+quadratic_terms+multiplications+6*outputs+3
        self.counters['integer_exact_divisions'] += 4*modes
        self.counters['integer_absolute_values'] += 6*modes+2*outputs
        self.counters['scalar_fraction_operations'] += 5
        self.counters['scalar_fraction_comparisons'] += 2
        self.counters['scalar_explicit_fraction_constructors'] += 3
        self.counters['scalar_dyadic_conversions'] += error_conversion+rounding_conversions
        self.counters['scalar_fraction_to_float_conversions'] += 1
        self.counters['scalar_gcd_calls_upper'] += 13+error_conversion+rounding_conversions
        # O(1) bit-length checks on already-computed scalar quantities replace
        # per-product instrumentation. The coordinate l1 norm bounds each Y_i.
        d_bits, y_bits = D.bit_length(), norm_numerator.bit_length()
        c_bits, coefficient_bits = C.bit_length(), polynomial.coefficient_bits
        rhs_bits = max(d_bits+y_bits+coefficient_bits, 2*y_bits+coefficient_bits)
        rhs_bits += (linear_terms+quadratic_terms).bit_length()+3
        a_bits, b_bits = self._step_bits
        defect_bits = max(b_bits+c_bits+d_bits+y_bits+1, a_bits+rhs_bits)
        defect_bits += (2*outputs).bit_length()+2
        direct_bits = max(defect_bits, b_bits+c_bits+2*d_bits, rhs_bits)
        scalar_bits = max(self._constant_bits, *(max(abs(x.numerator).bit_length(),
            x.denominator.bit_length()) for x in
            (integrated_defect, displacement, self._error_value, norm, rounded_q)))
        for key, value in (
            ('maximum_binary_denominator_bits', d_bits),
            ('maximum_coordinate_numerator_bits_upper', y_bits),
            ('maximum_direct_integer_bits_upper', direct_bits),
            ('maximum_scalar_integer_bits_upper', 8*scalar_bits+16)):
            self.counters[key] = max(self.counters[key], value)
        self._left, self._right = left, right
        # Carry the upward-rounded dyadic error, preventing denominator growth
        # while preserving exactly the same accepted-step certification rule.
        self._error_key, self._error_value = upper, rounded_q
        return upper, integrated_defect, bool(tube), polynomial.logical_pairs
