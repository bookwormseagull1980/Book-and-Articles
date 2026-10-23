"""Independent checks of the finite particle trajectory enclosure."""
import unittest
from dataclasses import replace
from unittest.mock import patch

import mpmath as mp
import numpy as np

from cosmic_cg.params import load_model_params
from cosmic_cg.pm_nbody import (
    PMNBodyConfig, _spectral_setup, _initial_state, _advance_to,
    _particle_acceleration,
)
from cosmic_cg.pm_certificate import (
    ArithmeticWork, bounded_fft, certified_force, kernel_bounds,
    coefficient_intervals, exponential_upper, certify_pm_trajectory,
    _multipliers,
)


class PMArithmeticTests(unittest.TestCase):
    def test_fourier_ball_contains_independent_high_precision_dft(self):
        rng=np.random.default_rng(817)
        data=rng.integers(-32,33,size=(4,4,4))/64
        work=ArithmeticWork()
        result,radius=bounded_fft(data,work=work)
        with mp.workdps(70):
            roots=[mp.exp(-2j*mp.pi*j/4) for j in range(4)]
            for k in np.ndindex(data.shape):
                exact=sum(mp.mpf(float(data[x]))*roots[sum(a*b for a,b in zip(k,x))%4]
                          for x in np.ndindex(data.shape))
                self.assertLessEqual(abs(mp.mpc(complex(result[k]))-exact),radius)
        self.assertEqual(work.counts["complex_butterflies"],3*64*2//2)
        self.assertEqual(work.counts["real_multiplications"],4*192+2*6)

    def test_positive_coefficient_intervals_against_high_precision(self):
        params=load_model_params()
        with mp.workdps(70):
            omega=mp.mpf(params.omega_b)+mp.mpf(params.omega_sigma)
            lam=mp.mpf(params.omega_lambda)
            for a in (0.01,0.2,0.2137,0.4,1.0,3.0):
                aa=mp.mpf(a)
                alpha=1/mp.sqrt(omega*aa**3+lam*aa**6)
                beta=mp.mpf('1.5')*omega*aa*alpha
                bounds=coefficient_intervals(a,params)
                for exact,(lo,hi) in zip((alpha,beta),bounds):
                    self.assertLessEqual(mp.mpf(lo),exact)
                    self.assertGreaterEqual(mp.mpf(hi),exact)

    def test_taylor_growth_is_an_upper_bound(self):
        with mp.workdps(70):
            for z in (0.,0.001,0.05,0.9):
                self.assertGreaterEqual(mp.mpf(exponential_upper(z)),mp.exp(mp.mpf(z)))
        with self.assertRaises(ValueError):
            exponential_upper(1.)

    def test_cic_ball_contains_independent_real_half_spectrum_force(self):
        config=PMNBodyConfig(particles_per_dim=4,mesh_size=4,box_size_mpch=8)
        setup=_spectral_setup(config)
        rng=np.random.default_rng(517)
        points=(np.indices((4,4,4)).reshape(3,-1).T+0.5)*2
        points+=rng.integers(-8,9,size=points.shape)/128
        computed,error,_=certified_force(points,config,_multipliers(setup,4))
        # These dyadic CIC weights and their small sums are exact in binary64.
        density=np.zeros((4,4,4))
        records=[]
        for point in points:
            base=np.floor(point/2).astype(int)
            frac=point/2-base
            local=[]
            for shift in np.ndindex((2,2,2)):
                index=tuple((base+shift)%4)
                weight=np.prod([frac[d] if shift[d] else 1-frac[d] for d in range(3)])
                density[index]+=weight
                local.append((index,weight))
            records.append(local)
        modes=np.fft.rfftn(density-1)
        reference=[]
        for d,name in enumerate(("kx","ky","kz")):
            shape=[1,1,1]
            shape[d]=-1
            field=np.fft.irfftn(1j*np.asarray(setup[name],float).reshape(shape)*modes/
                               np.asarray(setup["k2"],float),s=(4,4,4),axes=(0,1,2))
            reference.append([sum(weight*field[index] for index,weight in record) for record in records])
        np.testing.assert_array_less(np.abs(computed-np.asarray(reference).T),error)

    def test_current_kernel_bound_controls_finite_displacements(self):
        config=PMNBodyConfig(particles_per_dim=4,mesh_size=4,box_size_mpch=8)
        multipliers=_multipliers(_spectral_setup(config),4)
        slope,_=kernel_bounds(multipliers,2)
        centres=(np.indices((4,4,4)).reshape(3,-1).T+0.5)*2
        rng=np.random.default_rng(77)
        excursion=0.1
        c=(1+excursion)**2
        d=(1+excursion)**3-1
        for _ in range(12):
            x=centres+rng.uniform(-excursion,excursion,centres.shape)
            y=centres+rng.uniform(-excursion,excursion,centres.shape)
            gx,ex,_=certified_force(x,config,multipliers)
            gy,ey,_=certified_force(y,config,multipliers)
            self.assertLessEqual(float(np.max(np.abs(gx-gy))),
                                 slope*(c+d)*float(np.max(np.abs(x-y)))+ex+ey)


class PMTrajectoryTests(unittest.TestCase):
    def config(self,steps=16):
        return PMNBodyConfig(particles_per_dim=4,mesh_size=4,box_size_mpch=8,
                            z_initial=4,steps_per_interval=steps,
                            initial_displacement_fraction=.002)

    def test_certificate_reuses_the_same_production_trajectory(self):
        config=self.config()
        params=load_model_params()
        setup=_spectral_setup(config)
        setup["params"]=params
        ordinary=_initial_state(config,setup,np)
        with patch("cosmic_cg.pm_nbody._particle_acceleration",wraps=_particle_acceleration) as force:
            _advance_to(ordinary,.2,.4,config,setup)
            self.assertEqual(force.call_count,config.steps_per_interval+1)
        result=certify_pm_trajectory(config,params,return_state=True)
        certified=result.pop("_final_state")
        np.testing.assert_array_equal(ordinary["positions"],certified["positions"])
        np.testing.assert_array_equal(ordinary["velocities"],certified["velocities"])
        self.assertLess(result["uniform_joint_error_upper"],result["continuation_tube"])
        self.assertFalse(result["spatial_continuum_certified"])
        self.assertEqual(result["work"]["certificate_force_evaluations"],17)

    def test_wrong_clock_and_failed_tube_are_rejected(self):
        params=load_model_params()
        with self.assertRaisesRegex(ValueError,"clocks"):
            certify_pm_trajectory(replace(self.config(),z_initial=9),params)
        with self.assertRaisesRegex(ValueError,"leaves tube"):
            certify_pm_trajectory(self.config(),params,tube=1e-12)

    def test_force_update_does_not_transform_unused_potential(self):
        config=self.config()
        setup=_spectral_setup(config)
        setup["params"]=load_model_params()
        state=_initial_state(config,setup,np)
        with patch("numpy.fft.irfftn",wraps=np.fft.irfftn) as inverse:
            _particle_acceleration(state["positions"],config,setup)
            self.assertEqual(inverse.call_count,3)


if __name__=="__main__":
    unittest.main()
