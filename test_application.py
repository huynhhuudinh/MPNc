"""Independent checks of the changed observation design and likelihood."""
import unittest
from itertools import product
import numpy as np
from scipy.special import logsumexp
from scipy.stats import multinomial, poisson
from run_application import (ApplicationLikelihood, unequal_probabilities, load_data,
                             working_start, composition_gof, ROOT)
from softclosure import Likelihood, Scenario, generate, seeded_rng, working_fit
from rda_reader import read_cswa


def reference_site(y, mu, r, dep, model):
    """Original pass-allocation sum, without using the totals factorization."""
    pi,a,_=unequal_probabilities(r)
    visits=y[~np.isnan(y).all(axis=1)].astype(int)
    if model=='MPN0':
        return float(sum(poisson.logpmf(v,mu*pi).sum() for v in visits))
    lam=mu*dep if model=='MPNc' else (mu/dep if model=='CRTE' else mu)
    result=0.
    for k in range(int(poisson.isf(1e-14,lam))+15):
        contribution=poisson.pmf(k,lam)
        for v in visits:
            if model=='MPNc':
                conditional=0.
                for z in product(*(range(int(c)+1) for c in v)):
                    if sum(z)<=k:
                        conditional+=multinomial.pmf([*z,k-sum(z)],k,[*pi,1-a])*np.prod(
                            poisson.pmf(v-np.array(z),(1-dep)*mu*pi))
            else:
                pp=dep*pi if model=='CRTE' else pi
                conditional=multinomial.pmf([*v,k-int(v.sum())],k,[*pp,1-pp.sum()]) if k>=v.sum() else 0.
            contribution*=conditional
        result+=contribution
    return np.log(result)


class ApplicationTests(unittest.TestCase):
    def setUp(self):
        self.y=np.array([[[1.,0.,0.],[0.,1.,0.],[np.nan]*3],
                         [[0.,0.,0.],[1.,0.,1.],[0.,1.,0.]]])
        self.area=np.array([.4,.8])

    def test_unequal_intervals_from_ten_minute_probabilities(self):
        r=.27
        minute=r*(1-r)**np.arange(10)
        pi,a,q=unequal_probabilities(r)
        np.testing.assert_allclose(pi,[minute[:2].sum(),minute[2:5].sum(),minute[5:].sum()],atol=1e-15)
        self.assertAlmostEqual(a,1-(1-r)**10,places=14)
        self.assertAlmostEqual(q.sum(),1.,places=14)

    def test_data_axes_area_join_and_missingness(self):
        y,area,ids,audit=load_data()
        self.assertEqual((audit['sites'],audit['observed_visits'],audit['complete_sites']),(43,128,42))
        np.testing.assert_array_equal(y[0],[[1,0,0],[2,0,0],[1,0,1]])
        self.assertEqual(audit['missing_visits'],[dict(point_id='darty',visit=1)])
        self.assertAlmostEqual(area[ids.index('zeala')],.78)
        self.assertAlmostEqual(area[ids.index('zealb')],.50)
        self.assertAlmostEqual(area[ids.index('zeal4a')],.20)
        self.assertEqual(audit['spot_count'],74)
        self.assertAlmostEqual(audit['spot_area_ha'],66.69)

    def test_all_four_against_original_multinomial_likelihood(self):
        like=ApplicationLikelihood(self.y,self.area)
        for model in ['MPN0','MPN1','MPNc','CRTE']:
            expected=[reference_site(y,2.3*area,.19,.43,model) for y,area in zip(self.y,self.area)]
            np.testing.assert_allclose(like.site_loglik(model,2.3,.19,.43),expected,atol=2e-11,rtol=0)

    def test_missing_visit_is_omitted_not_zero_count(self):
        like=ApplicationLikelihood(self.y[:1],self.area[:1])
        subset=ApplicationLikelihood(self.y[:1,:2],self.area[:1])
        zero=ApplicationLikelihood(np.nan_to_num(self.y[:1]),self.area[:1])
        for model in ['MPN0','MPN1','MPNc','CRTE']:
            actual=like.site_loglik(model,2.,.25,.5)
            np.testing.assert_allclose(actual,subset.site_loglik(model,2.,.25,.5),atol=1e-12)
            self.assertGreater(float(actual[0]-zero.site_loglik(model,2.,.25,.5)[0]),0.)

    def test_homogeneous_equal_intervals_match_simulation_implementation(self):
        y=np.nan_to_num(self.y).astype(int)
        app=ApplicationLikelihood(y,np.ones(2),durations=[1,1,1])
        sim=Likelihood(y)
        for model in ['MPN0','MPN1','MPNc','CRTE']:
            np.testing.assert_allclose(app.site_loglik(model,2.1,.23,.61),sim.site_loglik(model,2.1,.23,.61),atol=1e-12)

    def test_two_visits_bivariate_poisson_with_area_offsets(self):
        y=self.y[:,:2]
        like=ApplicationLikelihood(y,self.area)
        pi,a,q=unequal_probabilities(.17)
        expected=[]
        for yy,area in zip(y,self.area):
            d=yy.sum(axis=1).astype(int);mu=1.7*area
            comp=sum(multinomial.logpmf(v,int(v.sum()),q) for v in yy)
            shared=.6*mu*a*a;exclusive=mu*a-shared
            z=np.arange(min(d)+1)
            expected.append(comp+logsumexp(poisson.logpmf(z,shared)+poisson.logpmf(d[0]-z,exclusive)+poisson.logpmf(d[1]-z,exclusive)))
        for model in ['MPNc','CRTE']:
            np.testing.assert_allclose(like.site_loglik(model,1.7,.17,.6),expected,atol=2e-11,rtol=0)

    def test_reject_partial_missing_and_bad_exposure(self):
        partial=self.y.copy();partial[0,0,0]=np.nan
        with self.assertRaises(ValueError):ApplicationLikelihood(partial,self.area)
        with self.assertRaises(ValueError):ApplicationLikelihood(self.y,[0.,.8])

    def test_independent_rda_decoder_matches_normalized_data(self):
        counts,dimnames,covs,spot=read_cswa(ROOT/'data'/'raw'/'cswa.rda')
        saved=np.load(ROOT/'data'/'analysis_data.npz')
        self.assertEqual(list(saved['point_id']),dimnames[0])
        np.testing.assert_array_equal(np.isnan(saved['y']),np.isnan(counts.transpose(0,2,1)))
        np.testing.assert_array_equal(np.nan_to_num(saved['y']),np.nan_to_num(counts.transpose(0,2,1)))
        np.testing.assert_allclose(saved['area_ha'],covs.loc[dimnames[0],'plotArea'].to_numpy(float))
        self.assertEqual(int(spot.CSWA.sum()),74)

    def test_composition_goodness_of_fit(self):
        y,area,_,_=load_data()
        gof=composition_gof(ApplicationLikelihood(y,area),bootstrap=199)
        self.assertEqual(gof['observed'],[51,14,10])
        self.assertAlmostEqual(gof['pearson_x2'],5.444,places=3)
        self.assertAlmostEqual(gof['likelihood_ratio_g2'],5.681,places=3)
        self.assertAlmostEqual(sum(gof['expected']),75.,places=10)
        # Counts proportional to a constant-hazard composition fit almost exactly.
        q=unequal_probabilities(.3)[2]
        close=ApplicationLikelihood(np.array([[np.round(q*100000)]]),[1.])
        result=composition_gof(close,bootstrap=9)
        self.assertLess(result['likelihood_ratio_g2'],1e-3)
        self.assertAlmostEqual(result['r_composition_mle'],.3,places=4)

    def test_simulation_working_profile_is_independence_mle(self):
        from scipy.optimize import minimize
        from scipy.special import expit
        s=Scenario('core','MPNc',200,3,2.,.25,.5)
        like=Likelihood(generate(s,seeded_rng(1,s.id,0)))
        w=working_fit(like)
        f=lambda x:-like.site_loglik('MPN0',np.exp(x[0]),expit(x[1])).sum()
        o=minimize(f,[0.,0.],method='Nelder-Mead',options=dict(xatol=1e-10,fatol=1e-12,maxiter=5000))
        self.assertLessEqual(f(w['x']),o.fun+1e-9)
        np.testing.assert_allclose([w['mu'],w['r']],[np.exp(o.x[0]),expit(o.x[1])],rtol=1e-6)


if __name__=='__main__':
    unittest.main(verbosity=2)
