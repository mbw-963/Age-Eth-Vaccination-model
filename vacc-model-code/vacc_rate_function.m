function r=vacc_rate_function(nu,f,e,S,alpha_max)

sus_elig=e.*f.*S;
alpha=min([alpha_max,nu/sum(sus_elig)]);

r=alpha*e.*f;
