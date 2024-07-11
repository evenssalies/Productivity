/******************************************************************************************/
/*	Formating in progress																																	*/
/*  Estimation de la productivité à la Levinsohn et Petrin (2003) 												*/
/*	Evens Salies, OFCE (Sciences Po)																											*/
/*	v1: obtenir les mêmes résultats que levpet et prodest dans le cas de la VA, 03/2023.	*/
/*	v2: acélération de ka routine, 2-3 fois plus rapide que levpet et prodest, 04/2023.		*/
/*	v3: comparer à un algorithme de descente https://www.geeksforgeeks.org/how-						*/
/*	to-implement-a-gradient-descent-in-python-to-find-a-local-minimum/										*/
/*																																												*/	
/*	Données reprises de la commande prodest (Rovigatti et Mollisi, 2018, Theory and				*/
/*	practice of total-factor productivity estimation: The control function approach 			*/
/*	using Stata, The Stata Journal, 18(3), pp. 618-662)                                   */
/******************************************************************************************/

cd					"C:\Users\evens\Documents"
cls
macro drop all

/*	To test my routine I replicate LP results obtained from Rovigati and Molisi's (2018) sample */
use					"chile_analysis.dta", clear

/*	Keep relevant variables only, rename ... */
order				*, alpha
drop				CIIU3 ele2 ele3 ele4 exit k2* k3* k4 k_ele* fuel // water
rename			(ANIO CIIU2 ID    ele    go   inv k)(year a100  siren l_ener l_ca l_i l_k)
rename			(skilled   unskilled   va water )(l_skilled l_unskilled l_v l_water)
label variable 	year ""
label variable 	a100 ""
label variable 	siren ""
order			year a100 siren, last
xtset			siren year

generate		l_hours=l_skilled+l_unskilled
clonevar		l_m=l_ener	/* Proxy, notée l_m, même si c'est l'énergie (élec.) */
des

/*	La routine commence ici !!! */

/*	Variables utilisées : siren, year, l_v, l_k, l_hours, l_m 	*/
/* STAGE 1.	beta_l */
/* 		Steps 1-4. 	Régressions polynomiales sur la proxy (m) et k.
				Estimation de phi(m,k), polynôme en 
				m, k, m^2, k^2, mk, m^3, k^3, m^2k, mk^2.
				Remarque : dans le sous-échantillon des y observé, pas
						   seulement des l_hours observés */
display         c(rmsg_time)
forvalues		i=0(1)3 {
 local			k=3-`i'
 forvalues		j=0(1)`k' {
  generate		mk`i'`j'=l_m^`i'*l_k^`j' if l_v!=.
 }
}

drop			mk00
*order			mk10 mk01 mk20 mk02 mk11 mk30 mk03 mk21 mk12, after(l_m)

foreach			x of varlist l_v l_hours { /* et l_ener si CA */
 regress 		`x' mk* if l_v!=.
 predict 		`x'_res if l_v!=., residuals
}

* 		Step 5.	 No-constant regression
regress 		l_v_res l_hours_res, nocons /// /* et l_ener_res si CA */
				vce(bootstrap, reps(3))
/*		Récupère beta_l */
scalar 			beta_l=_b[l_hours_res]	/* scalar beta_e = _b[l_ener_res] si CA */

drop			l_v_res l_hours_res 

/*	STAGE 2. beta_k, omega */
/* 		Step 1.	phi_hat_t */
generate		phi=l_v-beta_l*l_hours /* et -beta_e*l_ener si CA */
regress			phi mk*
predict			phi_hat if l_v!=., xb
drop			mk* 

/*		Step 2.		beta_k et processus de Markov par MCNL */
generate		vphihat1=L1.phi_hat
generate		vphihat2=vphihat1^2
generate		vphihat3=vphihat1^3
generate		vk1=L1.l_k
generate		vk2=vk1^2
generate		vk3=vk1^3
generate		vphihatk1=2*vphihat1*vk1
generate		vphihatk2=3*vphihat2*vk1
generate		vphihatk3=3*vphihat1*vk2
generate		nottouse=(vphihat1==.|vphihat2==.|vphihat3==.|vk1==.|vk2==.|vk3==.| ///
				vphihatk1==.|vphihatk2==.|vphihatk3==.)
generate		newsiren=siren
*xtset			newsiren year	/* Si clustered bootstrap, sinon, pas la peine */
nl 				(phi={b0}+{b1}*vphihat1+{b2}*vphihat2+{b3}*vphihat3 ///
					+{bk}*l_k-{bk}*{b1}*vk1-{bk}*{b2}*vphihatk1-{bk}*{b3}*vphihatk2 ///
					+{bk}^2*{b2}*vk2+{bk}^2*{b3}*vphihatk3-{bk}^3*{b3}*vk3) ///
				if nottouse==0, initial(bk 0) delta(1e-5)
				// vce(bootstrap, reps(3))
				// cluster(siren) idcluster(newsiren)
				// vce(hc2)
local			prodest_evens=c(rmsg_time)

/*	Productivité solution */
/*		Récupère beta_k   */
matrix define	betaomega=e(b)
scalar			beta_k=betaomega[1,5]

/*		Récupère xi (l'erreur du processus de Markov) */
quietly {
 xtset			siren year
 generate		omega=phi_hat-beta_k*l_k 		/* -`beta_m'*l_m si CA */
 generate		omega_1=L1.omega
 generate		omega_2=omega_1*L1.omega
 generate		omega_3=omega_2*L1.omega
 regress		omega omega_*, noheader notable
 predict		omegaxi, residuals
}
drop			omega_1-omega_3

/* 		levpet estimation de la productivité + eta */
generate		myomega_levpet=phi-beta_k*l_k
/* 		prodest estimation de la productivité */
generate		myomega_prodest=phi_hat-beta_k*l_k

/* La routine s'arrête ici !!! */

/**************************************************/
/*	Comparaison des temps de calcul :			  */
/*	Productivités estimées avec levpet et prodest */
/*		Remarque : prodest estime 2 productivités */
/*				   à la LP (avec et sans eta) 	  */
/**************************************************/

display			c(rmsg_time)
levpet			l_v, free(l_hours) proxy(l_m) capital(l_k) valueadded ///
				reps(3)
local			levpet=c(rmsg_time)
predict			levpet_omega, omega
replace			levpet_omega=log(levpet_omega)
display			c(rmsg_time)
prodest			l_v, free(l_hours) proxy(l_m) state(l_k) valueadded ///
				method(lp) poly(3) fsresiduals(prodest_resid1) ///
				tolerance(1e-5) reps(3)
local			prodest=c(rmsg_time)
set rmsg off
display		_newline(20)
display		"me (" `prodest_evens' "), levpet (" `levpet' "), prodest (" `prodest' ")"
display		_newline(20)

s

/*************************************************************************************/
/* Superposition des densités de la productivité de ma routine avec celle de prodest */
/*************************************************************************************/

/*		prodest (productivité à la LP nette de l'erreur idiosyncratique) */ 
predict			prodest_omega, omega
/* 		levpet (productivité à la LP) */
predict			prodest_resid2, residuals

/* 	Vérifie dans un nuage les estimations des productivités de prodest et la mienne */
capture drop	omegadiff omegakeep fractile
generate		phidiff=phi-phi_hat
generate		omegadiff=levpet_omega-prodest_omega
gsort			-omegadiff
generate		fractile=100*_n/_N

/* 		Sur les densités empiriques */
*histogram		myomega_prodest
*graph export 	"tfp_lp_omega.png", as(png) replace
replace			myomega_prodest=log(myomega_prodest)
replace			prodest_omega=log(prodest_omega)

set seed		21041971
generate		normal=rnormal()
summarize		myomega_prodest
generate		myomega_prodestcr=(myomega_prodest-r(mean))/r(sd)

summarize		prodest_omega
generate		prodest_omegacr=(prodest_omega-r(mean))/r(sd)

generate		prod_hours=l_v-l_hours
summarize		prod_hours
generate		prod_hourscr=(prod_hours-r(mean))/r(sd)

graph twoway 	(histogram normal if (normal>-3)&(normal<3), ///
					gap(10) color(green*.2) fcolor(none)) ///
					(kdensity myomega_prodestcr if ///
						(myomega_prodestcr>-3)&(myomega_prodestcr<3), ///
				kernel(epanechnikov) bwidth(0.5) lcolor(red%75) ///
					lwidth(medthin)) ///
				(kdensity prod_hourscr if (prod_hourscr>-3)&(prod_hourscr<3), ///
					kernel(epanechnikov) bwidth(0.5) lcolor(blue) ///
					lwidth(medthin)) ///
			(kdensity prodest_omegacr if (prodest_omegacr>-3)&(prodest_omegacr<3), ///
					kernel(epanechnikov) bwidth(0.5) lcolor(yellow%75) ///
					lwidth(medthin)), ///
				xtitle("Productivité standardisée") ytitle("Densité") ///
				xscale(noline titlegap(3)) yscale(titlegap(3)) ///
				legend(label(1 "N(0,1)") ///
				       label(2 "P.G.F. (evens)") ///
   				       label(3 "P. du travail") ///
				       label(4 "P.G.F. (prodest)") ///
				region(lstyle(none)) rows(2)) ///
				scheme(s1color) ///
				plotregion(fcolor(white)) graphregion(fcolor(white))
graph export 	"tfp_lp.png", as(png) replace
save			"file.dta", replace
