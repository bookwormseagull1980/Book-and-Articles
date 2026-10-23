"""Independent continuum-integral checks of the initial force enclosure."""
import unittest
from itertools import product
import mpmath as mp
import numpy as np

from cosmic_cg.pm_nbody import PMNBodyConfig,_spectral_setup
from cosmic_cg.pm_certificate import _multipliers,certified_force
from cosmic_cg.pm_spatial_certificate import initial_force_certificate


class InitialNewtonForceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config=PMNBodyConfig(particles_per_dim=4,mesh_size=4,box_size_mpch=8)
        cls.setup=_spectral_setup(cls.config)
        cls.multipliers=_multipliers(cls.setup,4)

    def evaluate(self,amplitude):
        grid=np.zeros((3,4,4,4),np.float32)
        grid[0]=np.asarray([0,amplitude,0,-amplitude],np.float32)[:,None,None]
        centres=(np.indices((4,4,4)).reshape(3,-1).T+0.5)*2
        cell=np.floor(centres/2).astype(int)
        displacement=(np.asarray(grid[0,:,0,0],float)+np.roll(np.asarray(grid[0,:,0,0],float),-1))/2
        positions=centres.copy()
        positions[:,0]+=displacement[cell[:,0]]
        positions=positions.astype(np.float32).astype(float)
        capture={}
        certified_force(positions,self.config,self.multipliers,spectrum_capture=capture)
        result=initial_force_certificate({"displacement_grid":grid,"displacement_scale":1.},
                                         self.config,self.multipliers,capture,bands=(1,))
        return result,grid

    def test_exact_cell_integrals_are_inside_density_and_force_enclosures(self):
        # A piecewise affine cold sheet permits closed-form cell integrals,
        # independent of the displacement expansion and Gauss quadrature.
        for amplitude in (0.01,0.125):
            result,grid=self.evaluate(amplitude)
            with mp.workdps(80):
                values=[mp.mpf(float(v)) for v in grid[0,:,0,0]]
                for row in result["mode_enclosures"]:
                    mx,my,mz=row["mode"]
                    exact=mp.mpc(0)
                    if my==0 and mz==0:
                        k=2*mp.pi*mx/8
                        for j in range(4):
                            origin=2*j+values[j]
                            slope=2+values[(j+1)%4]-values[j]
                            exact+=mp.exp(-1j*k*origin)*mp.expm1(-1j*k*slope)/(-1j*k*slope)/4
                    ball=row["Newton_density_coefficient"]
                    centre=mp.mpc(ball["real"],ball["imaginary"])
                    self.assertLessEqual(abs(exact-centre),mp.mpf(ball["radius"]))
                    knorm=2*mp.pi*mp.sqrt(mx*mx+my*my+mz*mz)/8
                    norm=abs(exact)/knorm
                    lo,hi=row["Newton_force_norm"]
                    self.assertLessEqual(mp.mpf(lo),norm)
                    self.assertGreaterEqual(mp.mpf(hi),norm)
            self.assertTrue(result["all_checks_passed"])

    def test_uniform_background_cancels_and_zero_reference_has_no_relative_error(self):
        result,_=self.evaluate(0.)
        self.assertTrue(result["all_checks_passed"])
        self.assertIsNone(result["bands"][0]["relative_L2_force_error"])
        for row in result["mode_enclosures"]:
            ball=row["Newton_density_coefficient"]
            self.assertLessEqual(abs(complex(ball["real"],ball["imaginary"])),ball["radius"])

    def test_all_bands_share_transforms_and_do_not_integrate_a_new_trajectory(self):
        result,_=self.evaluate(.01)
        self.assertEqual(result["work"]["additional_3d_complex_FFTs"],51)
        self.assertEqual(result["work"]["reused_density_transforms"],1)
        self.assertEqual(result["work"]["new_particle_trajectories"],0)
        self.assertFalse(result["all_mode_spatial_error_certified"])
        self.assertFalse(result["evolved_Newtonian_trajectory_certified"])

    def test_three_dimensional_cross_terms_against_unexpanded_cell_integral(self):
        wave=np.asarray([0,1,0,-1],np.float32)
        grid=np.empty((3,4,4,4),np.float32)
        grid[0]=.04*wave[:,None,None]*wave[None,:,None]
        grid[1]=.03*wave[None,None,:]
        grid[2]=.02*wave[:,None,None]*wave[None,None,:]
        base=np.indices((4,4,4),dtype=float)
        centre=np.zeros_like(grid,dtype=float)
        for shift in product((0,1),repeat=3):
            centre+=np.roll(grid,tuple(-s for s in shift),axis=(1,2,3))/8
        positions=((base+.5)*2+centre).reshape(3,-1).T
        capture={}
        certified_force(positions,self.config,self.multipliers,spectrum_capture=capture)
        result=initial_force_certificate({"displacement_grid":grid,"displacement_scale":1.},
                                         self.config,self.multipliers,capture,bands=(1,))
        selected=((1,1,1),(1,0,1),(-1,1,0))
        integral={mode:0j for mode in selected}
        nodes,weights=np.polynomial.legendre.leggauss(12)
        nodes=(nodes+1)/2
        weights=weights/2
        for ids in product(range(12),repeat=3):
            t=[nodes[j] for j in ids]
            displacement=np.zeros_like(grid,dtype=float)
            for shift in product((0,1),repeat=3):
                weight=np.prod([t[d] if shift[d] else 1-t[d] for d in range(3)])
                displacement+=weight*np.roll(grid,tuple(-s for s in shift),axis=(1,2,3))
            positions=(base+np.asarray(t)[:,None,None,None])*2+displacement
            quadrature_weight=np.prod([weights[j] for j in ids])
            for mode in selected:
                phase=2*np.pi*np.sum(np.asarray(mode)[:,None,None,None]*positions,axis=0)/8
                integral[mode]+=quadrature_weight*np.mean(np.exp(-1j*phase))
        for row in result["mode_enclosures"]:
            mode=tuple(row["mode"])
            if mode in integral:
                ball=row["Newton_density_coefficient"]
                self.assertLessEqual(abs(integral[mode]-complex(ball["real"],ball["imaginary"])),ball["radius"])


if __name__=="__main__":
    unittest.main()
