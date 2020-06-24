mata: mata set matastrict on

mata:
real matrix crt_crpb  ///
	(real scalar m, real scalar trans,  ///
	real scalar rho, real scalar tau,  ///
	real scalar ord, real scalar transrecr,  ///
	real scalar step)
	
/*  Designing a cluster randomised trial with continuous    */
/*  recruitment and prospective baseline.                   */

/*  This Mata function returns a matrix of coordinates for  */
/*  plotting a graph of variance of treatment effect        */
/*  estimator against cross-over time.                      */

/*  Arguments:                                              */
/*  m is the number of eligible participants presenting at  */
/*  a cluster during the recruitment period; trans is the   */
/*  proportion of the recruitment period taken up by the    */
/*  transition period; rho is the ICC for two participants  */
/*  sampled from the same cluster at the same time; tau is  */
/*  the decay in the ICC over the recruitment period; ord   */
/*  is the degree or order of the polynomial time effect;   */
/*  transrecr is a flag to indicate whether recruitment     */
/*  continues in the control arm during transition (any     */
/*  non-zero value) or is suspended in the control arm      */
/*  during transition (zero value); step is a flag to       */
/*  indicate whether a step-change at cross-over should be  */
/*  included with the underlying polynomial time effect.    */
/*  To specify a step-change time effect with no other      */
/*  time trend, ord=0 and step=1.                           */
	
{
	real matrix varplot, seq
	real scalar j
	
	trans = floor(m*trans)
	varplot = (((0::(m-1)) :/ m), J(m, 1, .))
	
	for (j=trans; j<m; j++) {
	    if (transrecr==0) {
			seq = (J(1, j-trans, 0),  ///
				J(1, trans, .), J(1, m-j, 1))  ///
					\ (J(1, j-trans, 0),  ///
					J(1, trans, .), J(1, m-j, 0))
		}
		else {
			seq = (J(1, j-trans, 0),  ///
				J(1, trans, .), J(1, m-j, 1))  ///
					\ J(1,m,0)
		}
		varplot[1+j,2] = glsvar_polystept(seq, rho, tau,  ///
			ord, j*(step != 0))
	}
	
	return(varplot)
}
end
		
mata:
real scalar glsvar_polystept(real matrix seq,  ///
		real scalar rho, real scalar tau,  ///
		real scalar ord, real scalar jstep)
{
	real matrix x, y, x1, y1, ccom, design, cincom
	real colvector t
	real scalar nper, nseq, i, j, var
	
	nper = cols(seq)
	nseq = rows(seq)
	
	y = J(2+ord+(jstep > 0), 0, .)
	x = J(0, 2+ord+(jstep > 0), .)

	t = (1 :: nper) :/ nper :- 0.5
	ccom = rho :* (J(nper, nper, tau) :^ abs  ///
		(J(nper, 1, t') - J(1, nper, t)))  ///
		+ (1-rho) :* I(nper)
		
	for (i=1; i<=nseq; i++) {
		design = (seq[i,.]', J(nper, 1, 1))
		for (j=1; j<=ord; j++) {
			design = (design, (t :^ j))
		}
		if (jstep > 0) {
			design = (design, (J(jstep, 1, 0) \  ///
				J(nper-jstep, 1, 1)))
		}
		x1 = select(design, seq[i,.]' :!= .)
		cincom = select(  ///
			select(ccom, seq[i,.]' :!= .),  ///
			seq[i,.] :!= .)
		y1 = x1' * invsym(cincom)
		x = x \ x1
		y = y , y1
	}
	
	var = invsym(y*x)[1,1]
	return(var)
}
end


/*  For example, here's how to generate (x,y) coordinates   */
/*  for the curve plotted top left in Figure 4 of the       */
/*  article -- m=100, no transition period, rho=0.05,       */
/*  tau=0.5, linear effect of time (no discontinuity):      */

mata:
crt_crpb(100,0,0.05,0.5,1,0,0)
end

