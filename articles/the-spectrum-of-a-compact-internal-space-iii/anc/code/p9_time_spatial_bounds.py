"""Whole-interval spatial-force bounds using sparse actual-trajectory anchors.

The target is the explicitly reconstructed continuous state. No extra particle
trajectory is evolved and no Newtonian time-evolution claim is made.
"""
import itertools
import json
import math
import numpy as np
from mpmath import iv

from p9_cic_replacement_validation import (original_trajectory,
    continuous_reconstruction_metrics, up, positive_sum, WORK,
    _spectral_setup, _multipliers, _weights, bounded_fft, ArithmeticWork, U)
from p9_continuous_force_bound import (spline_norm_upper, spectral_factors,
    multiply_ball, vector_norm_interval, nonlinear_coefficient, interval_constant,
    ComplexBall, ball_add, ball_mul, _add, _mul, _div, _down,
    positive_square_root_upper, positive_square_root_lower, lower_difference,
    tail_norm_upper, sum_of_squared_norms_bounds)


def density_modes(x,cfg,work):
    n=cfg.mesh_size;N=len(x);G=n**3
    mass=np.zeros((n,n,n))
    for index,weight in _weights(x,cfg,work):
        np.add.at(mass,index,weight)
    mean=N/G
    assert math.frexp(mean)[0]==.5
    numerator=up((8*N+16)*U)
    gamma=up(numerator/float(np.nextafter(1-2*numerator,-math.inf)))
    radius=up(up(gamma*max(float(mass.max()),1)/mean)+4*U)
    modes,mr=bounded_fft(mass/mean-1,radius,work=work)
    work.counts['real_additions']=work.counts.get('real_additions',0)+8*N+G
    work.counts['binary_scalings']=work.counts.get('binary_scalings',0)+G
    return modes/G,up(mr/G)


def field_modes(delta,dr,multipliers):
    modes=[];radii=[]
    for multiplier in multipliers:
        modes.append(delta*multiplier)
        mm=up(float(np.max(np.abs(multiplier.real)+np.abs(multiplier.imag))))
        mx=up(float(np.max(np.abs(delta.real)+np.abs(delta.imag))))
        radii.append(_add(_mul(mm,dr),_add(_mul(32*U,_mul(mm,mx)),32*np.nextafter(0.,1.))))
    return np.asarray(modes),radii


def mass_eigenvalues(n):
    iv.dps=40
    one=[]
    for r in range(n):
        value=(2+iv.cos(2*iv.pi*r/n))/3
        one.append(up(float(value.b)))
    one=np.asarray(one)
    return np.nextafter(np.nextafter(one[:,None,None]*one[None,:,None],np.inf)*one[None,None,:],np.inf)


def fourier_spline_norm(modes,radii,mass):
    terms=[]
    for values,radius in zip(modes,radii):
        absolute=np.nextafter(np.sqrt(np.nextafter(
            np.nextafter(values.real*values.real,np.inf)+np.nextafter(values.imag*values.imag,np.inf),np.inf)),np.inf)
        mag=np.nextafter(absolute+radius,np.inf)
        sq=np.nextafter(mag*mag,np.inf)
        terms.extend(np.nextafter(mass*sq,np.inf).ravel())
    return positive_square_root_upper(positive_sum(terms))


def projected_balls(mode,values,work=None):
    """Rational longitudinal projector, including coefficient division error."""
    m2=sum(d*d for d in mode);result=[]
    for b in range(3):
        value=ComplexBall(0j,0.)
        for d in range(3):
            numerator=mode[b]*mode[d]
            coefficient=numerator/m2
            radius=0. if numerator==0 else _add(_mul(2*U,abs(coefficient),work),np.nextafter(0.,1.),work)
            if work is not None:work.add(real_divisions=1)
            value=ball_add(value,ball_mul(ComplexBall(complex(coefficient),radius),values[d],work),work)
        result.append(value)
    return result


def one_anchor(source,increment,source_modes,source_radii,force_modes,force_radii,factors,M,mass,L):
    work=ArithmeticWork();G=source.shape[-1]**3
    modes=[];radii=[]
    for field in increment:
        value,rad=bounded_fft(field,work=work);modes.append(value/G);radii.append(up(rad/G))
    modes=np.asarray(modes)
    psi_norm=_add(spline_norm_upper(source,work=work),spline_norm_upper(increment,work=work),work)
    force_norm=fourier_spline_norm(force_modes,force_radii,mass)
    nonlinear=_mul(nonlinear_coefficient(L),psi_norm,work)
    errors=[];psis=[];forces=[];newtons=[]
    n=source.shape[-1]
    for m in itertools.product(range(-M,M+1),repeat=3):
        if m==(0,0,0):continue
        w,wr=1.,0.;phase,pr=1.+0j,0.
        for component in m:
            f,fr,z,zr=factors[component]
            w,wr=multiply_ball(w,wr,f,fr,work);phase,pr=multiply_ball(phase,pr,z,zr,work)
        index=tuple(component%n for component in m)
        psi=[];force=[]
        for b in range(3):
            shifted=ball_mul(ComplexBall(complex(modes[(b,*index)]),radii[b]),ComplexBall(phase,pr),work)
            total=ball_add(ComplexBall(complex(source_modes[(b,*index)]),source_radii[b]),shifted,work)
            psi.append(ball_mul(total,ComplexBall(w,wr),work))
            force.append(ball_mul(ComplexBall(complex(force_modes[(b,*index)]),force_radii[b]),ComplexBall(w,wr),work))
        proj=projected_balls(m,psi,work)
        difference=[ball_add(proj[b],ComplexBall(-force[b].centre,force[b].radius),work) for b in range(3)]
        for vec,acc in ((difference,errors),(psi,psis),(force,forces),(proj,newtons)):
            acc.append(vector_norm_interval([z.centre for z in vec],[z.radius for z in vec],work))
    error_lo,error_hi=sum_of_squared_norms_bounds(errors,work)
    psi_lo,_=sum_of_squared_norms_bounds(psis,work)
    force_lo,_=sum_of_squared_norms_bounds(forces,work)
    newton_lo,_=sum_of_squared_norms_bounds(newtons,work)
    psi_tail=tail_norm_upper(psi_norm,max(0.,_down(psi_lo*psi_lo,work)),work)
    force_tail=tail_norm_upper(force_norm,max(0.,_down(force_lo*force_lo,work)),work)
    tails=_add(psi_tail,force_tail,work)
    ehi=_add(positive_square_root_upper(_add(_mul(error_hi,error_hi,work),_mul(tails,tails,work),work),work),nonlinear,work)
    elo=lower_difference(error_lo,nonlinear,work)
    normlo=lower_difference(newton_lo,nonlinear,work)
    return {'error_interval':[elo,ehi],'Newton_force_norm_lower':normlo,
        'psi_norm_upper':psi_norm,'nonlinear_remainder_upper':nonlinear,
        'displacement_tail_upper':psi_tail,'CIC_force_tail_upper':force_tail,
        'relative_error_upper':None if not normlo else _div(ehi,normlo,work)},work.counts


def source_operator_bound(multipliers,h,n,e):
    iv.dps=40
    difference=[]
    for r in range(n):
        value=abs(2*iv.sin(iv.pi*r/n)/h)
        difference.append(up(float(value.b)))
    difference=np.asarray(difference)
    d2=np.nextafter(difference**2,np.inf)
    dnorm=np.nextafter(np.sqrt(np.nextafter(np.nextafter(d2[:,None,None]+d2[None,:,None],np.inf)+d2[None,None,:],np.inf)),np.inf)
    h2=np.zeros((n,n,n))
    for value in multipliers:
        absolute=np.nextafter(np.abs(value)*(1+32*U),np.inf)
        h2=np.nextafter(h2+np.nextafter(absolute**2,np.inf),np.inf)
    singular=np.nextafter(np.sqrt(h2),np.inf)*dnorm
    symbol=up(float(np.max(singular))*(1+32*U))
    transverse=_add(1,_div(_mul(2,e),h))
    return symbol,_mul(symbol,transverse)


def main():
    cfg,params,times,x,p,g,D,s,target,_=original_trajectory()
    reconstruction=continuous_reconstruction_metrics(times,x,p,D,s,cfg,target)
    n=cfg.mesh_size;G=n**3;h=cfg.box_size_mpch/n
    source=np.asarray(D,float)*s;setup=_spectral_setup(cfg);multipliers=_multipliers(setup,n)
    work=ArithmeticWork();sm=[];sr=[]
    for field in source:
        mode,rad=bounded_fft(field,work=work);sm.append(mode/G);sr.append(up(rad/G))
    sm=np.asarray(sm);M=4;factors=spectral_factors(n,M);mass=mass_eigenvalues(n)
    L=reconstruction['uniform_reconstruction_Lip_Euclidean_upper']
    e=max(row['max_particle_displacement_from_original_centre'] for row in reconstruction['rows'])
    e=up(e*(1+16*U))
    symbol,cic_lip=source_operator_bound(multipliers,h,n,e)
    iv.dps=40
    newton_lip=up(float((1/iv.sqrt((1-iv.mpf(L))**3)).b))
    combined=up(newton_lip+cic_lip)
    anchors={};anchor_work={};anchor_fields={}
    for index in range(0,65,8):
        dm,dr=density_modes(x[index]%cfg.box_size_mpch,cfg,work)
        fm,fr=field_modes(dm,dr,multipliers)
        inc=(x[index]-x[0]).reshape(n,n,n,3).transpose(3,0,1,2)
        row,local_work=one_anchor(source,inc,sm,sr,fm,fr,factors,M,mass,L)
        row['a']=float(times[index]);anchors[index]=row;anchor_work[index]=local_work
        anchor_fields[index]=(fm,fr)
    distances={}
    def distance(first,second):
        key=tuple(sorted((first,second)))
        if key not in distances:
            delta=x[first]-x[second]
            # Every subtraction satisfies Sterbenz. Positive products and
            # summation use a conventional upper rounding allowance.
            assert np.all((x[first]>=x[second]/2)&(x[first]<=2*x[second]))
            sq=np.nextafter(delta*delta,np.inf)
            distances[key]=up(math.sqrt(up(positive_sum(sq.ravel())/len(x[0]))))
        return distances[key]
    rows=[]
    for j in range(64):
        candidates={8*(j//8),8*((j+8)//8)}
        choices=[]
        for k in sorted(candidates):
            d=max(distance(j,k),distance(j+1,k))
            choices.append((up(anchors[k]['error_interval'][1]+up(combined*d)),k,d))
        error,k,d=min(choices)
        normlo=lower_difference(anchors[k]['Newton_force_norm_lower'],_mul(newton_lip,d))
        rows.append({'a_interval':[float(times[j]),float(times[j+1])],'anchor_index':k,
            'particle_RMS_distance_upper':d,'all_frequency_reconstruction_force_error_upper':error,
            'Newton_force_norm_lower':normlo,'relative_error_upper':None if not normlo else up(error/normlo)})
    source_checks=[]
    for k in range(0,64,8):
        f0,r0=anchor_fields[k];f1,r1=anchor_fields[k+8]
        dr=[up(r0[b]+r1[b]+32*U*max(1.,float(np.max(np.abs(f0[b]))),float(np.max(np.abs(f1[b]))))) for b in range(3)]
        field_change=fourier_spline_norm(f1-f0,dr,mass)
        displacement=distance(k,k+8)
        assert field_change<=cic_lip*displacement
        source_checks.append({'a_interval':[float(times[k]),float(times[k+8])],
            'force_L2_change_upper':field_change,'Lipschitz_prediction_upper':up(cic_lip*displacement),
            'force_change_per_particle_RMS_upper':up(field_change/displacement)})
    result={'scope':'Eulerian L2 spatial force discrepancy on explicit continuous reconstruction, uniformly in a; not the exact Newtonian trajectory',
        'initial_state_sha256':target['initial_state_sha256'],'operator_sha256':target['finite_operator_data_sha256'],
        'upstream_run_id':params.source_run_id,'time_interval':[.2,.4],'anchor_indices':list(anchors),
        'modes_per_anchor':(2*M+1)**3-1,'anchors':anchors,'intervals':rows,
        'uniform_absolute_force_error_upper':max(row['all_frequency_reconstruction_force_error_upper'] for row in rows),
        'uniform_relative_force_error_upper':max(row['relative_error_upper'] for row in rows),
        'grid_source_symbol_norm_upper':symbol,'CIC_L2_source_Lipschitz_upper':cic_lip,
        'Newton_L2_source_Lipschitz_upper':newton_lip,'combined_Lipschitz_upper':combined,
        'uniform_continuous_reconstruction_Lip_upper':L,
        'reconstructed_density_interval':[
            float(np.nextafter(float((1/(1+iv.mpf(L))**3).a),-math.inf)),
            up(float((1/(1-iv.mpf(L))**3).b))],
        'Newtonian_evolved_trajectory_certified':False,
        'independent_actual_anchor_source_checks':source_checks,
        'work':{'additional_3d_complex_ffts':3+4*len(anchors),'additional_particle_evolutions':0,
            'additional_density_depositions':len(anchors),'sparse_anchor_count':len(anchors),
            'mode_observations':len(anchors)*((2*M+1)**3-1),'RMS_distance_evaluations':len(distances),
            'RMS_real_arithmetic_upper':len(distances)*(9*len(x[0])+8),
            'RMS_directed_roundings_upper':len(distances)*(6*len(x[0])+8),
            'core_fft_and_density_arithmetic':work.counts,'anchor_arithmetic':anchor_work,
            'remaining_scalar_and_spline_arithmetic_upper':1024*len(anchors)*((2*M+1)**3-1)+256*G*len(anchors)+256*G,
            'reconstruction_geometry_arithmetic_upper':24*G*len(times)+128*len(times),
            'interval_trigonometric_calls_upper':3*(2*M+1)+2*n,
            'FFT_root_trigonometric_enclosures':4*(n-1),
            'interval_square_root_calls_upper':2*len(anchors)+1,
            'interval_precision_decimal_digits':40,
            'future_reuse':'Capture production density Fourier arrays at these anchors; only three displacement-increment FFTs per anchor remain.'},
        'all_checks_passed':all(math.isfinite(row['all_frequency_reconstruction_force_error_upper']) for row in rows)}
    output=WORK/'p9_time_spatial_bounds_20261004.json'
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'output':str(output),'absolute_error_upper':result['uniform_absolute_force_error_upper'],
        'relative_error_upper':result['uniform_relative_force_error_upper'],
        'operator_norms':[symbol,cic_lip,newton_lip],'complex_ffts':result['work']['additional_3d_complex_ffts'],
        'RMS_evaluations':len(distances)},indent=2))


if __name__=='__main__':
    main()
