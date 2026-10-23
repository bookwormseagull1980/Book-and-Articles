"""Sustained active-only evolution of the finite matter system.

The reference is the alias-free, matter-dominated Galerkin model from
Appendix B, with N positive Fourier modes. Only M<=N modes are advanced.
Every accepted step has a rationally evaluated defect enclosure for its
actual floating-point, piecewise-linear trajectory. Production computes no
full reference trajectory. Independent validation is an explicit separate
entry point, with its work reported separately.
"""
from fractions import Fraction as Q
import math
from .exact_certificate import ExactStepCertificate


def matter_rhs(state, modes, output_modes=None):
    """Cosine density and sine velocity coefficients; no FFT aliasing."""
    output_modes = modes if output_modes is None else output_modes
    zero = state[0] * 0
    result = [zero] * (2*output_modes)
    pairs = 0
    for n in range(1, output_modes+1):
        if n <= modes:
            result[n-1] = -n*state[modes+n-1]
            result[output_modes+n-1] = -state[modes+n-1]/2 - 3*state[n-1]/(2*n)
        for p in range(-modes, modes+1):
            q = n-p
            if p == 0 or q == 0 or abs(q) > modes:
                continue
            sign_p, sign_q = (1 if p > 0 else -1), (1 if q > 0 else -1)
            ap, bp = state[abs(p)-1], state[modes+abs(p)-1]
            bq = state[modes+abs(q)-1]
            result[n-1] -= n*sign_q*ap*bq/2
            result[output_modes+n-1] -= n*sign_p*sign_q*bp*bq/4
            pairs += 1
    return result, pairs


def _embed(state, modes, target):
    zero = state[0]*0
    return list(state[:modes]) + [zero]*(target-modes) + list(state[modes:]) + [zero]*(target-modes)


def _up(value):
    """Round a nonnegative exact rational to an enclosing IEEE float."""
    result = float(value)
    if Q.from_float(result) < value:
        result = math.nextafter(result, math.inf)
    return result


def _certificate_step(state, candidate, modes, full_modes, step, error, radius):
    # Every actual stored float is an exact dyadic rational in this check.
    y = [Q.from_float(float(v)) for v in state]
    next_y = [Q.from_float(float(v)) for v in candidate]
    evaluated = min(2*modes, full_modes)
    exact_rhs, pairs = matter_rhs(y, modes, evaluated)
    slope = [(right-left)/step for left, right in zip(y, next_y)]
    slope_embedded = _embed(slope, modes, evaluated)
    defect = sum(abs(slope_embedded[i]-exact_rhs[i]) for i in range(2*evaluated))
    speed = sum(abs(v) for v in slope)
    lip = Q(full_modes) + Q(1, 2) + 2*full_modes*radius
    if lip*step >= 1:
        raise ValueError("Step exceeds the rational amplification bound")
    bound = (Q.from_float(error) + step*defect + lip*speed*step**2/2)/(1-lip*step)
    tube = max(sum(abs(v) for v in y), sum(abs(v) for v in next_y)) + bound < radius
    return _up(bound), float(defect), bool(tube), pairs



def initial_state(full_modes):
    """Shared initial datum for the declared finite reference problem."""
    initial = [0.0]*(2*full_modes)
    initial[0], initial[full_modes] = 0.01, -0.01
    initial[full_modes-1], initial[-1] = 1e-8, -1e-8/full_modes
    return initial


def float_advance_work(modes, pairs):
    """Actual scalar float arithmetic of matter_rhs followed by Euler.

    Loop indices, integer coefficients, comparisons and allocation are not
    float operations. Ordered convolution pairs are a separate logical count.
    """
    return {"multiplications": 1+4*modes+4*pairs,
            "additions_subtractions": 3*modes+2*pairs,
            "divisions": 2*modes+2*pairs, "unary_negations": modes,
            "rhs_coordinates": 2*modes,
            "examined_pair_indices": modes*(2*modes+1),
            "valid_ordered_convolution_pairs": pairs}


def _add_work(total, work):
    for key, value in work.items():
        total[key] = total.get(key, 0)+value


class ReserveBudget:
    """A float selection policy; the independent exact certificate decides safety.

    For q=(1-Lh)^(-1), K steps and initial error E0, the ideal schedule
    A[j+1]=q*A[j]+b, b=(eps-q**K*E0)*(q-1)/(q**K-1), ends at eps.
    Floating-point roundoff here cannot weaken certification: every accepted
    step must separately satisfy the exact outward bound <= eps and the tube.
    """
    def __init__(self, amplification, error, tolerance, steps):
        self.q = float(amplification)
        power = self.q**steps
        if not math.isfinite(power) or power*error >= tolerance:
            raise ValueError("The initial truncation error exhausts this selection budget")
        self.injection = (tolerance-power*error)*(self.q-1)/(power-1)
        self.allowance, self.tolerance, self.steps, self.index = error, tolerance, steps, 0

    def advance(self):
        self.index += 1
        if self.index == self.steps:
            self.allowance = self.tolerance
        else:
            self.allowance = min(self.tolerance, self.q*self.allowance+self.injection)
        return self.allowance

    def work(self):
        return {"scope": "floating-point selection policy, not a second certificate",
                "float_power_evaluations": 1,
                "multiplications": 3+max(0, self.index-1),
                "additions_subtractions": 3+max(0, self.index-1),
                "divisions": 1, "allowance_comparisons": max(0, self.index-1),
                "recurrence_updates": max(0, self.index-1)}


def sustained_active_run(*, full_modes=8, duration=Q(1, 5), step=Q(1, 640),
                         tolerance=1e-4, full=False, validation_history=False):
    """Advance and rigorously certify the specified finite Fourier model.

    Production defaults to active coordinates only and does not compute a
    full trajectory, DOP reference, or full-dimensional history. Set full=True
    for the same certified method on the complete finite carrier. History is
    optional validation data and never available to the selector.
    """
    if type(full_modes) is not int or not isinstance(duration, Q) or not isinstance(step, Q):
        raise TypeError("Exact certification requires an integer mode count and Fraction duration and step")
    if full_modes < 3 or duration <= 0 or step <= 0 or duration/step != int(duration/step):
        raise ValueError("Invalid finite carrier or time grid")
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("Tolerance must be positive")
    radius, count = Q(3, 100), int(duration/step)
    modes = full_modes if full else 1
    state = initial_state(full_modes) if full else [0.01, -0.01]
    error = 0.0 if full else _up(abs(Q.from_float(1e-8))+abs(Q.from_float(-1e-8/full_modes)))
    certificate = ExactStepCertificate(full_modes, step, radius)
    budget = None if full else ReserveBudget(certificate.amplification, error, tolerance, count)
    h = float(step)
    records, events = [], []
    accepted_work, rejected_work = {}, {}
    candidate_attempts = rejected_attempts = certificate_coordinates = 0
    record_interval = max(1, count//20)
    history = [_embed(state, modes, full_modes)] if validation_history else None
    error_bounds = [error] if validation_history else None
    max_error, min_margin = error, math.inf
    for n in range(count):
        allowance = tolerance if full else budget.advance()
        while True:
            rhs, pairs = matter_rhs(state, modes)
            candidate = [value+h*deriv for value, deriv in zip(state, rhs)]
            candidate_attempts += 1
            work = float_advance_work(modes, pairs)
            next_error, integrated_defect, tube, _ = certificate.evaluate(state, candidate, modes, error)
            certificate_coordinates += 2*min(2*modes, full_modes)
            if tube and next_error <= tolerance and next_error <= allowance:
                _add_work(accepted_work, work)
                break
            _add_work(rejected_work, work)
            rejected_attempts += 1
            if modes == full_modes:
                raise RuntimeError("The requested tolerance needs a smaller step or a different selection budget")
            previous = modes
            state = _embed(state, modes, modes+1)
            modes += 1
            events.append({"step": n, "time": float(n*step), "from_modes": previous,
                           "to_modes": modes, "reason": "next-step error reserve", "reconstruction_jump": 0.0})
        state, error = candidate, next_error
        state_norm = sum(abs(v) for v in state)
        max_error = max(max_error, error)
        min_margin = min(min_margin, float(radius)-state_norm-error)
        if validation_history:
            history.append(_embed(state, modes, full_modes))
            error_bounds.append(error)
        if (n+1) % record_interval == 0 or n+1 == count:
            records.append({"step": n+1, "time": float((n+1)*step), "active_modes": modes,
                            "error_upper": error, "defect_norm": float(integrated_defect/step),
                            "state_norm": state_norm})
    total_work = dict(accepted_work)
    _add_work(total_work, rejected_work)
    checks = {"uniform_bound_below_tolerance": max_error <= tolerance,
              "strict_tube": min_margin > 0, "positive_density": float(radius) < 1}
    if not all(checks.values()):
        raise RuntimeError(f"Certified evolution check failed: {checks}")
    record_work = {"records": len(records), "events": len(events),
                   "state_norm_absolute_values": accepted_work.get("rhs_coordinates", 0),
                   "state_norm_additions": accepted_work.get("rhs_coordinates", 0),
                   "diagnostic_fraction_divisions": len(records),
                   "diagnostic_time_fraction_multiplications": len(records)+len(events),
                   "tube_margin_float_subtractions": 2*count,
                   "bound_and_margin_comparisons": 2*count,
                   "validation_state_scalars_stored": 2*full_modes*(count+1) if validation_history else 0,
                   "validation_bound_scalars_stored": count+1 if validation_history else 0,
                   "final_active_state_scalars": len(state)}
    result = {"scope": "alias-free finite matter-dominated Fourier model; conditional matter evolution",
              "initial_density": "0.01 cos(x)+1e-8 cos(Nx)",
              "initial_velocity": "-0.01 sin(x)-1e-8 sin(Nx)/N",
              "full_positive_modes": full_modes, "steps": count, "duration_log_a": float(duration),
              "step": str(step), "tolerance": tolerance, "error_norm": "joint Fourier l1",
              "certificate_arithmetic": "exact integer polynomial of stored IEEE coordinates; exact rational scalar bound rounded upward every step",
              "uniform_error_upper": max_error, "minimum_tube_margin": min_margin,
              "density_lower_bound": 1-float(radius), "events": events, "records": records,
              "final_active_modes": modes, "certified_full_carrier": full,
              "selection_policy": None if full else "A[j+1]=q*A[j]+b; floating policy; each accepted step independently certified",
              "work": {"scope": "typed arithmetic and storage counts; categories are not added with arbitrary weights",
                       "advancement_float": total_work, "accepted_advancement_float": accepted_work,
                       "rejected_advancement_float": rejected_work,
                       "candidate_attempts": candidate_attempts, "rejected_attempts": rejected_attempts,
                       "certificate": dict(certificate.counters),
                       "certificate_rhs_coordinates": certificate_coordinates,
                       "selection": {} if full else budget.work(),
                       "input_setup": {"initial_state_float_divisions": 1,
                           "initial_error_fraction_from_float": 0 if full else 3,
                           "initial_error_fraction_additions": 0 if full else 1,
                           "initial_error_fraction_absolute_values": 0 if full else 2,
                           "initial_error_fraction_to_float": 0 if full else 1,
                           "initial_error_fraction_comparisons": 0 if full else 1,
                           "initial_error_nextafter_upper": 0 if full else 1},
                       "acceptance_comparisons_upper": 3*candidate_attempts,
                       "records_and_storage": record_work},
              "checks": checks, "all_checks_passed": all(checks.values())}
    if validation_history:
        result["_validation_history"] = history
        result["_validation_bounds"] = error_bounds
    return result


def verify_spectral_pair(*, full_modes=8, duration=Q(1, 5), step=Q(1, 640), tolerance=1e-4):
    """Certify active/full paths separately and share one optional numerical check.

    Both methods have a uniform exact certificate at the same grid and target.
    The DOP trajectory is computed after both runs, used only for validation,
    evaluated once, and never fed into selection or production.
    """
    import numpy as np
    from scipy.integrate import solve_ivp
    active = sustained_active_run(full_modes=full_modes, duration=duration, step=step,
                                   tolerance=tolerance, validation_history=True)
    full = sustained_active_run(full_modes=full_modes, duration=duration, step=step,
                                 tolerance=tolerance, full=True, validation_history=True)
    reference_work = {}
    calls = 0
    def reference_rhs(time, state):
        nonlocal calls
        result, pairs = matter_rhs(state, full_modes)
        work = float_advance_work(full_modes, pairs)
        work["multiplications"] -= 2*full_modes
        work["additions_subtractions"] -= 2*full_modes
        _add_work(reference_work, work)
        calls += 1
        return result
    reference = solve_ivp(reference_rhs, (0.0, float(duration)), initial_state(full_modes),
                          method="DOP853", rtol=1e-12, atol=1e-14, dense_output=True)
    if not reference.success:
        raise RuntimeError(reference.message)
    points = active["steps"]+1
    times = np.arange(points)*float(step)
    sampled_reference = reference.sol(times).T
    for run in (active, full):
        states = np.asarray(run.pop("_validation_history"))
        bounds = np.asarray(run.pop("_validation_bounds"))
        errors = np.sum(np.abs(states-sampled_reference), axis=1)
        positive = bounds > 0
        fraction = float(np.max(errors[positive]/bounds[positive]))
        run["max_sampled_reference_error"] = float(np.max(errors))
        run["max_sampled_error_fraction_of_bound"] = fraction
        run["checks"]["sampled_reference_agrees"] = fraction <= 1.00000001 and bool(np.all(errors[~positive] == 0))
        run["all_checks_passed"] = all(run["checks"].values())
        if not run["all_checks_passed"]:
            raise RuntimeError(f"Independent sampled check failed: {run['checks']}")
    active["full_certified_reference"] = full
    active["verification_work"] = {
        "scope": "optional shared independent check, separate from both required certified evolutions",
        "DOP853_rhs_evaluations": calls, "DOP853_rhs_float": reference_work,
        "solver_internal_arithmetic": "not included in RHS counts; no total solver operation count claimed",
        "dense_reference_evaluations": 1, "dense_reference_sample_points": points,
        "dense_reference_scalar_values": 2*full_modes*points,
        "comparison_subtractions": 4*full_modes*points, "comparison_absolute_values": 4*full_modes*points,
        "comparison_additions": 2*points*(2*full_modes-1),
        "stored_trajectory_scalars": 4*full_modes*points}
    return active
