*********************************************
* Thirteen new DID estimators 
* 1. Simulated data  
*********************************************
// Generate a complete panel of 400 units observed in 15 periods

	clear all
	timer clear
	set seed 10
	global T = 15
	global I = 400

	set obs `=$I*$T'
	gen i = int((_n-1)/$T )+1 		// unit id
	gen t = mod((_n-1),$T )+1		// calendar period
	tsset i t

	
// Randomly generate treatment rollout periods uniformly across Ei=10..16
// (note that periods t>=16 would not be useful since all units are treated by then)

	*Period when unit is first treated
	gen Ei = ceil(runiform()*7)+$T -6 if t==1
	bys i (t): replace Ei = Ei[1]
	
	*Relative time, i.e. number of periods since treated (could be missing if never treated)
	gen K = t-Ei
	
	*Treatment indicator
	gen D = K>=0 & Ei!=.

	
// Generate the outcome with parallel trends and heterogeneous treatment effects

	*Heterogeneous treatment effects (in this case vary over calendar periods)
	gen tau = cond(D==1, (t-12.5), 0)
	
	*Error term
	gen eps = rnormal()
	
	*Outcome (FEs play no role since all methods control for them)
	gen Y = i + 3*t + tau*D + eps

		
*********************************************
* 2. Estimators
*********************************************		

**************************************************		
// （1）did_imputation of Borusyak et al. (2021)
**************************************************
	*Estimation
	did_imputation Y i t Ei, allhorizons pretrends(5)
	// Y:	outcome variable
	// i:	unit id variable
	// t:	time period variable
	// Ei:	variable for unit-specific treatment date (never-treated: Ei == missing)

	// allhorizons: include all non-negative horizons available
	// pretrends(): number of pre-treatment coefficients to be estimated (with too many pre-trend coefficients, the power of the joint test will be lower.)
	// standard errors are clustered at unit level by default

	*Plotting
	event_plot, default_look graph_opt(xtitle("Periods since the event") ytitle("Average causal effect") ///
		title("Borusyak et al. (2021) imputation estimator") xlabel(-5(1)5) name(BJS, replace)) together

	*Storing estimates for later
	estimates store bjs
	
**********************************************************************	
// （2）did_multiplegt of de Chaisemartin and D'Haultfoeuille (2020)
**********************************************************************

	*Estimation
	did_multiplegt Y i t D, robust_dynamic dynamic(5) placebo(5) longdiff_placebo breps(100) cluster(i)
	// Y:	outcome variable
	// i:	unit id variable
	// t:	time period variable
	// D:	treatment variable
	
	// robust_dynamic: uses estimator from dCdH (2021) on DID with intertemporal effects
	// dynamic(): number of dynamic treatment effects to be estimated (can only be used with robust_dynamic)
	// placebo(): number of placebo estimates to be estimated (When the longdiff_placebo and dynamic options are requested, the number of placebos requested cannot be larger than the number of dynamic effects)
	// longdiff_placebo: estimates placebo effects using long differences (comparable to dynamic TE estimates)
	// breps(): number of bootstrap iterations for computation of standard errors
	// cluster(i): computes standard errors using block bootstrap at level specified
	
	// Note that, according to the help file, by default "placebos are
	// first-difference estimators, while dynamic effects are long-difference
	// estimators, so they are not really comparable." Thus, we should not plot
	// and compare them in the same graph, if the "longdiff_placebo" option is
	// not specified.

	*Plotting
	event_plot e(estimates)#e(variances), default_look graph_opt(xtitle("Periods since the event") ///
		ytitle("Average causal effect") title("de Chaisemartin and D'Haultfoeuille (2020)") xlabel(-5(1)5) ///
		name(dCdH, replace)) stub_lag(Effect_#) stub_lead(Placebo_#) together
	// bmat#vmat: name of point estimate matrix and name of variance-covariance matrix
	// stub_lag((prefix#postfix): name of lag coefficients in estimation output
	// stub_lead(prefix#postfix): name of lead coefficients in estimation output
	// graph_opt(): twoway options for graph overall
	// together: show leads and lags as one line

	*Storing estimates for later
	matrix dcdh_b = e(estimates)
	matrix dcdh_v = e(variances)

***********************************************************************************************
// （3）csdid of Callaway and Sant'Anna (2020) (v1.6 written by Fernando Rios-Avila @friosavila)
***********************************************************************************************
	*Preparation
	gen gvar = cond(Ei>15, 0, Ei) // group variable as required for the csdid command

	*Estimation
	csdid Y, ivar(i) time(t) gvar(gvar) agg(event)
	// Y: 		outcome variable
	// ivar():	unit id variable
	// time():	time period variable
	// gvar():	variable for unit-specific treatment date (never treated: gvar == 0)
	// (defines "group" in CS jargon)
	
	// agg(): aggregation to use
	// wboot: Wild Bootstrap standard errors (default is asymptotic normal)
	// cluster(): should in principle be possible but throws an error (for me at least)
	// by default uses never treated units as control group (could specify "notyet")
	
	// Note that this command is work in progress. As such, it may subject to
	// ongoing changes. For example, the wboot option currently seems to
	// throw an error if specified. Further, the confidence intervals are not
	// yet correct (not uniform as in CS).
	
	// Also, note that Nick Huntington-Klein provides a Stata package that
	// acts as a wrapper for the "did" package by CS in R (via rcall, i.e. 
	// need to have R installed). It is available on his github page:
	// https://github.com/NickCH-K/did

	*Plotting
	event_plot e(b)#e(V), default_look graph_opt(xtitle("Periods since the event") ///
		ytitle("Average causal effect") xlabel(-14(1)5) title("Callaway and Sant'Anna (2020)") name(CS, replace)) ///
		stub_lag(T+#) stub_lead(T-#) together

	*Storing estimates for later
	matrix cs_b = e(b)
	matrix cs_v = e(V)

**********************************************************************
	// （4）eventstudyinteract of Sun and Abraham (2020)
**********************************************************************
	*Preparation
	sum Ei
	gen lastcohort = Ei==r(max) // dummy for the latest- or never-treated cohort
	forvalues l = 0/5 {
		gen L`l'event = K==`l'
	}
	forvalues l = 1/14 {
		gen F`l'event = K==-`l'
	}
	drop F1event // normalize K=-1 (and also K=-15) to zero

	*Estimation
	eventstudyinteract Y L*event F*event, vce(cluster i) absorb(i t) cohort(Ei) control_cohort(lastcohort)
	// Y: outcome variable
	// L*event: lags to include
	// F*event: leads to include
	// vce(): options for variance-covariance matrix (cluster SE)
	// absorb(): absorb unit and time fixed effects
	// cohort(): variable for unit-specific treatment date (never-treated: Ei == missing)
	// control_cohort(): indicator variable for control cohort (either latest-treated or never-treated units)

	*Plotting
	event_plot e(b_iw)#e(V_iw), default_look graph_opt(xtitle("Periods since the event") ///
		ytitle("Average causal effect") xlabel(-14(1)5) title("Sun and Abraham (2020)") name(SA, replace)) ///
		stub_lag(L#event) stub_lead(F#event) together

	*Storing estimates for later
	matrix sa_b = e(b_iw)
	matrix sa_v = e(V_iw)
	
**********************************************************************
// （5）did2s of Gardner (2021)
**********************************************************************
	* Estimation
	did2s Y, first_stage(i.i i.t) second_stage(F*event L*event) treatment(D) cluster(i)
	// Y: outcome variable
	// first_stage(): fixed effects used to estimate counterfactual Y_it(0). This should be everything besides treatment
	// second_stage(): treatment such as dummy variable or event-study leads and lags
	// cluster(): variable to cluster on 

	
	*Plotting
	event_plot, default_look stub_lag(L#event) stub_lead(F#event) together ///
		graph_opt(xtitle("Periods since the event") ytitle("Average causal effect") xlabel(-14(1)5) ///
		title("Gardner (2021)") name(DID2S, replace))
		
	*Saving estimates for later
	matrix did2s_b = e(b)
	matrix did2s_v = e(V)

**********************************************************************	
// （6）stackedev of Cengiz et al. (2019) 
**********************************************************************
	*Create a variable equal to the time that the unit first received treatment. It must be missing for never treated units.
	gen treat_year=.
	replace treat_year=Ei if Ei!=16

	*Create never treated indicator that equals one if a unit never received treatment and zero if it did.
	gen no_treat= (Ei==16)

	*Preparation
	cap drop F*event L*event
	sum Ei
	forvalues l = 0/5 {
		gen L`l'event = K==`l'
		replace L`l'event = 0 if no_treat==1
	}
	forvalues l = 1/14 {
		gen F`l'event = K==-`l'
		replace F`l'event = 0 if no_treat==1
	}
	drop F1event // normalize K=-1 (and also K=-15) to zero

	* Run stackedev		
	preserve
	stackedev Y F*event L*event, cohort(treat_year) time(t) never_treat(no_treat) unit_fe(i) clust_unit(i) 
	restore
	// Y: outcome variable
	// L*event: lags to include
	// F*event: leads to include
	// time(): numerical variable equal to time
	// never_treat(): binary indicator that equals one if a unit never received treatment and zero if it did
	// unit_fe(): variable indicating unit fixed effects
	// clust_unit(): variable indicating the unit by which to cluster variances
	

	*Plotting
	event_plot e(b)#e(V), default_look graph_opt(xtitle("Periods since the event") ///
		ytitle("Average causal effect") xlabel(-14(1)5) title("Cengiz et al. (2019)") name(CDLZ, replace)) ///
		stub_lag(L#event) stub_lead(F#event) together	
		
	*Saving estimates for later
	matrix stackedev_b = e(b)
	matrix stackedev_v = e(V)
	
**********************************************************************	
// （7）TWFE OLS estimation
**********************************************************************
	*Estimation
	reghdfe Y F*event L*event, absorb(i t) vce(cluster i)

	*Plotting
	event_plot, default_look stub_lag(L#event) stub_lead(F#event) together ///
		graph_opt(xtitle("Days since the event") ytitle("OLS coefficients") xlabel(-14(1)5) ///
		title("OLS") name(OLS, replace))

	*Saving estimates for later
	estimates store ols

***************************************************************************************
// （8）Pedro H.C. Sant'Anna , Jun Zhao (2020). Doubly robust difference-in-differences 
// estimators, Journal of Econometrics.(Written by Fernando Rios-Avila, Pedro H.C. Sant'Anna, 
// Asjad Naqvi )
***************************************************************************************
// drdid
*     drdid Y, ivar(i) time(t) tr(D) all



// The general syntax of the command is as follows
// drdid depvar [indepvar] [if] [in] [iw], [ivar(varname)] time(varname) tr(varname) [estimator] [wboot rseed(#) cluster(varname) gmm stub(str) replace]

// Here an explanation of all the pieces:

// depvar : is your dependent variable or outcome of interest
// indepvar: are your independent variables, you may or may not have variables here. These variables will included in the outcome regression specification and the propensity score estimation.
// ivar: is a variable that identifies the panel ID. If you drop this, the command will use repeated crossection estimators instead. If included, it will estimate the panel estimators, after checking you indeed have a panel setup.
// time: identifies the time variable (for example year). It does not matter if the periods are contiguous or not. However, its important that you have only 2 values in time for the working sample. The earlier period will be used as pre, whereas the later period will be used as post.
// tr: is the treatment variable. It does not matter what values you use, as long as there are only two values in the used sample. Observations with lower values are the control group, whereas observations with the higher values are the treated groups.
// estimator is used to indicate which estimator you want to use. Below the list of all that is available:
// drimp Estimates the DR improved estimator. If you add rc1 it provides you with the alternative estimator (that is not locally efficient)
// dripw Estimates the DR IPW estimator. You can also use rc1 to provide the alternative (not locally efficient) estimator.
// reg Estimates the Outcome regression estimator.
// stdipw Estimates the Standard IPW estimator.
// ipw Estimates the estimator similar to Abadies (2005)
// ipwra Estimates the IPW regression augmented estimator. This is estimated through Stata command teffects ipwra. Not available for Repeated Crossection
// all Provides all estimators. Perhaps you may find it useful for robustness. (I did it for replication). This results cannot be used to make comparison across estimators.
// wboot: When used, the command will estimate the 1-step wild-bootstrap (or multiplier bootstrap), with 999 repetions. It currently implements Mammen(1993) approach. This does not work after all or after `ipwra'.
// cluster: When used, the command will estimate the clustered standard errors based. It can be used in combination with wboot.
// gmm: When used, requests the estimation of SE using gmm. You can combine this with other gmm SE estimators. See help gmm for other details. This estimator is available only for panel estimators.
// rseed(#): Used for replication purposes. Sets the Seed for random numbers.
// stub: When used, the recentered influence function will be stored in the dataset with the name stubatt. If the variable already exists, you can request to overwrite it using replace. This option does not work when one requests all estimators, or when ipwra estimator is requested.


*****************************************************************************************
// （9）Eva Dettmann, Alexander Giebler, Antje Weyh (2020). Flexpaneldid: A Stata Toolbox for Causal Analysis with Varying Treatment Time and Duration. IWH Discussion Papers No. 3/2020
*******************************************************************************************
// flexpaneldid
// Syntax

//    flexpaneldid depvar, id(varname) treatment(varname) time(varname) prepdataset(string)
//                         statmatching(con(varlist) cat(varlist)) [ties radius(float)]
//                         cemmatching(varname1 [(cutpoints1)] [varname2 [(cutpoints2)]...] [k2k])
//                         outcometimerelstart(integer) outcometimerelend(integer)
//                         [outcomedev(integer [integer]) test didmodel outcomemissing]

//    options                         Description
//    id(varname)                   name of the panel id variable
//      treatment(varname)            name of the treatment variable
//      time(varname)                 name of the time variable
//      prepdataset(filename)         filename where preprocessed dataset is stored

//    Distance metric
//      statmatching()                use statistical distance function for matching
//      or                            
//      cemmatching()                 use Coarsened Exact Matching cem

//    Outcome development
//      outcometimerelstart()         end of observed outcome development related to treatment start
//      or                            
//      outcometimerelend()           end of observed outcome development related to treatment end

//      outcomedev(int [int])         pre treatment outcome level[development] is included in the matching
//      test                          quality tests after matching
//      didmodel                      robustness checks on standard two-way DID model
//     outcomemissing                disables a check if depvar is observable for the defined outcome development period


 * flexpaneldid Y, id(i) treatment(D) time(t)






**********************************************************************
// （10）Simon Freyaldenhoven, Christian Hansen, Jesse M. Shapiro (2019). Pre-event Trends in the Panel Event-Study Design. American Economic Review.	
**********************************************************************
// xtevent
    xtevent Y, policyvar(D) panelvar(i) timevar(t) window(4) plot

// Syntax


// xtevent depvar [indepvars] [if] [in] , policyvar(varname) panelvar(varname) timevar(varname) [options]

//    options                    Description
    
//    Main
//    policyvar(varname)        policy variable
//      panelvar(varname)         variable that identifies the panels
//      timevar(varname)          variable that identifies the time periods
//      window(numlist)           # of periods in the estimation window
//      pre(integer)              # of periods with anticipation effects
//      post(integer)             # of periods with policy effects
//      overidpre(integer)        # of periods to test pretrends
//      overidpost(integer)       # of periods to test effects leveling off
//{synopt:norm(integer) event-time coefficient to normalize to 0{p_end}
//      proxy(varname)            proxy for the confound
//      proxyiv(string)           instruments for the proxy variable
//      nofe                      omit panel fixed effects
//      note                      omit time fixed effects
//      nostaggered               calculate endpoints without staggered adoption
//      static                    estimate static model
//      trend(#1)                 extrapolate linear trend from time period #1 before treatment
//      savek(stub)               save time-to-event, event-time and trend variables
//      kvars(stub)               use previously generated even-time variables
//      reghdfe                   use reghdfe for estimation
//      addabsorb(varlist)        absorb additional variables in reghdfe
//      plot                      display plot. See xteventplot.
 
//      additional_options        additional options to be passed to estimation command

	
**********************************************************************	
// （11）Damian Clarke, Kathya Tapia (2020). Implementing the Panel Event Study.
**********************************************************************
//eventdd
// Syntax
//        eventdd depvar [indepvars] [if] [in] [weight], timevar(varname) [options]
    	
	eventdd Y i.t,timevar(K) method(fe, cluster(i)) balanced graph_op(ytitle("Y"))

	
	
	
**********************************************************************	
// （12）Jonathan Roth , Pedro H.C. Sant'Anna  (2021). Efficient Estimation for Staggered Rollout Designs
**********************************************************************
// staggered_stata
// the efficient estimator under random timing proposed in Roth & Sant'Anna (2021), via the staggered command
// To install the staggered_stata package, the user first needs to install the github and RCall packages. 
// The following R packages are required by rcall. rcall attempts to detect R Statistical Software on your system automatically and install the dependency R packages. If the installation fails, read the rcall help file and install the dependencies manually.

	
 *   staggered, y("Y") g("gvar") t("t") i("i") estimand("simple")
	
	
	
	
	
**********************************************************************	
// （13）Jeff Wooldridge  (2021). Two-Way Fixed Effects, the Two-Way Mundlak Regression, and Difference-in-Differences Estimators
**********************************************************************	
//Syntax： 	jwdid  lemp , i(countyreal) t(year) gvar(first_treat)

   jwdid Y，i(i) t(t) gvar(Ei)
   
   
   
   // Construct vector of true average treatment effects by number of periods since treatment

	matrix btrue = J(1,6,.)
	matrix colnames btrue = tau0 tau1 tau2 tau3 tau4 tau5
	qui forvalues h = 0/5 {
		sum tau if K==`h'
		matrix btrue[1,`h'+1]=r(mean)
	}


// Combine all plots using the stored estimates (5 leads and lags around event)

event_plot btrue# bjs dcdh_b#dcdh_v cs_b#cs_v sa_b#sa_v did2s_b#did2s_v stackedev_b#stackedev_v ols, ///
	stub_lag(tau# tau# Effect_# T+# L#event L#event L#event L#event) stub_lead(pre# pre# Placebo_# T-# F#event F#event F#event F#event) ///
	plottype(scatter) ciplottype(rcap) ///
	together perturb(-0.325(0.1)0.325) trimlead(5) noautolegend ///
	graph_opt(title("Event study estimators in a simulated panel (400 units, 15 periods)", size(med)) ///
		xtitle("Periods since the event", size(small)) ytitle("Average causal effect", size(small)) xlabel(-5(1)5)  ///
		legend(order(1 "True value" 2 "Borusyak et al." 4 "de Chaisemartin-D'Haultfoeuille" ///
				6 "Callaway-Sant'Anna" 8 "Sun-Abraham" 10 "Gardner" 12 "Cengiz et al." 14 "TWFE OLS") rows(2) position(6) region(style(none))) ///
	/// the following lines replace default_look with something more elaborate
		xline(-0.5, lcolor(gs8) lpattern(dash)) yline(0, lcolor(gs8)) graphregion(color(white)) bgcolor(white) ylabel(, angle(horizontal)) ///
	) ///
	lag_opt1(msymbol(+) color(black)) lag_ci_opt1(color(black)) ///
	lag_opt2(msymbol(O) color(cranberry)) lag_ci_opt2(color(cranberry)) ///
	lag_opt3(msymbol(Dh) color(navy)) lag_ci_opt3(color(navy)) ///
	lag_opt4(msymbol(Th) color(forest_green)) lag_ci_opt4(color(forest_green)) ///
	lag_opt5(msymbol(Sh) color(dkorange)) lag_ci_opt5(color(dkorange)) ///
	lag_opt6(msymbol(Th) color(blue)) lag_ci_opt6(color(blue)) ///
	lag_opt7(msymbol(Dh) color(red)) lag_ci_opt7(color(red)) ///
	lag_opt8(msymbol(Oh) color(purple)) lag_ci_opt8(color(purple))
graph export "$output/seven_estimators_example_5t.png", replace

// Combine all plots using the stored estimates (5 leads and lags around event)

event_plot btrue# bjs dcdh_b#dcdh_v cs_b#cs_v sa_b#sa_v did2s_b#did2s_v stackedev_b#stackedev_v ols, ///
	stub_lag(tau# tau# Effect_# T+# L#event L#event L#event L#event) stub_lead(pre# pre# Placebo_# T-# F#event F#event F#event F#event) ///
	plottype(scatter) ciplottype(rcap) ///
	together perturb(-0.325(0.1)0.325) trimlead(14) noautolegend ///
	graph_opt(title("Event study estimators in a simulated panel (400 units, 15 periods)", size(med)) ///
		xtitle("Periods since the event", size(small)) ytitle("Average causal effect", size(small)) xlabel(-14(1)5)  ///
		legend(order(1 "True value" 2 "Borusyak et al." 4 "de Chaisemartin-D'Haultfoeuille" ///
				6 "Callaway-Sant'Anna" 8 "Sun-Abraham" 10 "Gardner" 12 "Cengiz et al." 14 "TWFE OLS") rows(2) position(6) region(style(none))) ///
	/// the following lines replace default_look with something more elaborate
		xline(-0.5, lcolor(gs8) lpattern(dash)) yline(0, lcolor(gs8)) graphregion(color(white)) bgcolor(white) ylabel(, angle(horizontal)) ///
	) ///
	lag_opt1(msymbol(+) color(black)) lag_ci_opt1(color(black)) ///
	lag_opt2(msymbol(O) color(cranberry)) lag_ci_opt2(color(cranberry)) ///
	lag_opt3(msymbol(Dh) color(navy)) lag_ci_opt3(color(navy)) ///
	lag_opt4(msymbol(Th) color(forest_green)) lag_ci_opt4(color(forest_green)) ///
	lag_opt5(msymbol(Sh) color(dkorange)) lag_ci_opt5(color(dkorange)) ///
	lag_opt6(msymbol(Th) color(blue)) lag_ci_opt6(color(blue)) ///
	lag_opt7(msymbol(Dh) color(red)) lag_ci_opt7(color(red)) ///
	lag_opt8(msymbol(Oh) color(purple)) lag_ci_opt8(color(purple))
graph export "$output/seven_estimators_example_allt.png", replace


