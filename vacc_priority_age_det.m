function dydt = vacc_priority_age_det(t,y,m,N,sz,b,a,nu,alpha_max,C)

sus=y(1:sz,:);
vacc=y(sz+1:2*sz,:);
z=y(2*sz+1:3*sz,:);
expsd=y(3*sz+1:4*sz,:);
v_expsd=y(4*sz+1:5*sz,:);
infd=y(5*sz+1:6*sz,:);
v_infd=y(6*sz+1:7*sz,:);
rcvd=y(7*sz+1:8*sz,:);
v_rcvd=y(8*sz+1:9*sz,:);

r_zero=1.3;
gamma=0.2;
beta=gamma*r_zero/max(eigs(C.*N'));
s=1;
ve_i=0;
ve_t=0;

indicators=vacc_indicator_function(m,z,N,sus,b,a);
f=indicators(:,1);
e=indicators(:,2);

r=vacc_rate_function(nu,f,e,sus,alpha_max);

dsdt=-beta*C*(infd+(1-ve_t)*v_infd).*sus-r.*sus;
dvdt=-beta*C*(1-ve_i)*(infd+(1-ve_t)*v_infd).*vacc+r.*sus;
dzdt=r.*(N-z);

dedt=beta*C*(infd+(1-ve_t)*v_infd).*sus-s*expsd;
dvedt=beta*C*(1-ve_i)*(infd+(1-ve_t)*v_infd).*vacc-s*v_expsd;

didt=s*expsd-gamma*infd;
dvidt=s*v_expsd-gamma*v_infd;

drdt=gamma*infd;
dvrdt=gamma*v_infd;

dydt=[dsdt;dvdt;dzdt;dedt;dvedt;didt;dvidt;drdt;dvrdt];