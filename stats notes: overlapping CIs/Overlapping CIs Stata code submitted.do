 ***** PART A
 ***** Analysis based on the separate CIs for groups A and B given in Figure 1 
 
* calculate ratio of odds ratios (ROR) for A compared to B
 disp "ROR = " 1.30/0.98
 
 
* calculate the standard errors (SE) of each group's treatment effect on lnOR sclae
disp "SE group A = " (ln(1.52)- ln(1.11))/(2*1.96)	
* SE group A = .08019141	
 
disp "SE group B = " (ln(1.23)- ln(0.78))/(2*1.96)
* SE group B = .11619274

* calculate the lower and upper bound of the CI for the ROR
disp "lower ROR = " exp(ln(1.30/0.98) - (1.96*sqrt((0.08019141^2) +(0.11619274^2 ))))
 
disp "upper ROR = " exp(ln(1.30/0.98) + (1.96*sqrt((0.08019141^2) +(0.11619274^2 ))))
  
* open dataset containing group specific ORs and the ROR, with their CIs
use "dataset 1.dta", replace
  
* generate data on lnOR scale
gen lnOR = ln(OR)
gen lowerlnOR = ln(lowerOR)
gen upperlnOR = ln(upperOR)
  
* generate forest plot (Figure 1)
metan lnOR lowerlnOR upperlnOR, eform  forest(lcols(Group )) astext(70) effect(odds ratio) nooverall xlab(0.75 1 1.5 2) nowt  xtitle(Odds Ratio, size(small))
  

***** PART B
***** Meta-Analysis of Figure 2
  
* open meta-analysis data 
use "dataset 2.dta", replace

* generate Figure 2  
metan lnoe selnoe ,  re(reml, hksj) forest(lcols(study ) effect(Odds Ratio) xlab(0.5 1 2) boxsc(50) xtitle(Odds Ratio, size(small)) rfopts(sepline))  eform  astext(50) nohet nowarning 
		
		
		