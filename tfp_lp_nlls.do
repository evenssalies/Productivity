/*	LP, Moindres carrés non-linéaires

		Evens Salies, 04/2023			*/

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

/* 		prodest estimation de la productivité */
generate		myomega_prodest=phi_hat-beta_k*l_k
