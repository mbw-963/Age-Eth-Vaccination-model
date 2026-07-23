clear

NZ_conmat=readmatrix('NZ_conmat.csv');
popdata=readtable('ERP_2025_base_2023.csv');
N=table2array(popdata(:,2:5));
age_vec=table2cell(popdata(:,1));
N_age=sum(N,2);
C=0.5*(NZ_conmat+(NZ_conmat.*N_age)'./N_age);

timespan=1:365;
N_t_col=reshape(N',[],1);
age_sz=length(N_age);
sz=length(N_t_col);
init_vacc_rate=0;
N_vacc=N_t_col*init_vacc_rate;
N_unvacc=N_t_col-N_vacc;

init_vals=readmatrix('init_vals.csv');

m=init_vals(:,1);
b=5300;
priority=init_vals(:,2);
nu=0;
alpha_max=1;

epsilon=0;
F=[1;1.19;1.87;0.85];
eth_sz=length(F);
M=matrix_extension_v2(N,F,epsilon,C);

infd0=ones(sz,1).*N_t_col*1e-4;
v_infd0=init_vals(:,4);
sus0=N_unvacc-infd0;
vacc0=N_vacc-v_infd0;
init_sus_ages=sum(reshape((sus0+vacc0),eth_sz,age_sz));
init_sus_eth=sum(reshape((sus0+vacc0),eth_sz,age_sz),2)';

z0=N_vacc;
expsd0=init_vals(:,5);
v_expsd0=init_vals(:,6);
rcvd0=init_vals(:,7);
v_rcvd0=init_vals(:,8);

[t,y]=ode45(@(t,y)vacc_priority_age_det(t,y,m,N_t_col,sz,b,priority,nu,alpha_max,M),timespan,[sus0;vacc0;z0;expsd0;v_expsd0;infd0;v_infd0;rcvd0;v_rcvd0]);

y1=reshape(y,size(y,1),eth_sz,[]);
timeseries_ages=squeeze(sum(y1,2));

y2=reshape(y,size(y,1),eth_sz,[]);
y2=reshape(permute(y2,[1,3,2]),size(y,1),[]);
y2=reshape(y2,size(y,1),age_sz,[]);
y2=reshape(squeeze(sum(y2,2)),size(y,1),9,[]);
timeseries_eth=reshape(permute(y2,[1,3,2]),size(y,1),[]);

total_infd_ages=init_sus_ages-timeseries_ages(end,1:age_sz)+...
    timeseries_ages(end,(age_sz+1):2*age_sz); % calculates total number of infected in each age group during the epidemic
total_prop_infd_ages=total_infd_ages./init_sus_ages;

infd_sus_matrix_ages=[(timeseries_ages(end,1:age_sz)+...
    timeseries_ages(end,(age_sz+1):2*age_sz));total_infd_ages];

total_infd_eth=init_sus_eth-timeseries_eth(end,1:eth_sz)+...
    timeseries_eth(end,(eth_sz+1):2*eth_sz); % calculates total number of infected in each ethnic group during the epidemic
total_prop_infd_eth=total_infd_eth./init_sus_eth; 

infd_sus_matrix_eth=[(timeseries_eth(end,1:eth_sz)+...
    timeseries_eth(end,(eth_sz+1):2*eth_sz));total_infd_eth];

total_infd=N_t_col'-(y(end,1:sz)+y(end,(sz+1):2*sz)); % total number of recovered in each age+ethnic group
total_prop_infd=total_infd./N_t_col';
total_prop_infd_mat=reshape(total_prop_infd,eth_sz,age_sz)';

age_cmap=colormap(slanCM('viridis',age_sz));
eth_cmap=colormap(slanCM('guppy',eth_sz));

%% age timeseries


figure(1)
plot(t,timeseries_ages(:,1:age_sz),'LineWidth',1);
hold on
plot(t,timeseries_ages(:,(age_sz+1):2*age_sz),'LineWidth',1,'LineStyle','--');
xlabel('Time (in days)')
ylabel('Number of susceptible people')
colororder(age_cmap)
xlim tight
%legend(age_vec)
hold off

%x1=["Child" "Adult" "Elderly"];
%bar(x1,[(sus0-y(end,1:3)')./sus0,(vacc0-y(end,4:6)')./vacc0])
%ylabel('Attack Rate')

figure(2)
plot(t,timeseries_ages(:,(3*age_sz+1):4*age_sz),'LineWidth',1);
hold on
plot(t,timeseries_ages(:,(4*age_sz+1):5*age_sz),'LineWidth',1,'LineStyle','--');
xlabel('Time (in days)')
ylabel('Number of exposed people')
colororder(age_cmap)
xlim tight
%legend(age_vec)
hold off

figure(3)
plot(t,timeseries_ages(:,(5*age_sz+1):6*age_sz),'LineWidth',1);
hold on
plot(t,timeseries_ages(:,(6*age_sz+1):7*age_sz),'LineWidth',1,'LineStyle','--');
xlabel('Time (in days)')
ylabel('Number of infected people')
colororder(age_cmap)
xlim tight
%legend(age_vec)
hold off

figure(4)
plot(t,timeseries_ages(:,(7*age_sz+1):8*age_sz),'LineWidth',1);
hold on
plot(t,timeseries_ages(:,(8*age_sz+1):9*age_sz),'LineWidth',1,'LineStyle','--');
xlabel('Time (in days)')
ylabel('Number of recovered people')
colororder(age_cmap)
xlim tight
hold off
%}
%% ethnicity timeseries

figure(5)
plot(t,timeseries_eth(:,1:eth_sz),'LineWidth',1);
hold on
plot(t,timeseries_eth(:,(eth_sz+1):2*eth_sz),'LineWidth',1,'LineStyle','--');
xlabel('Time (in days)')
ylabel('Number of susceptible people')
colororder(eth_cmap)
xlim tight
%legend({'NZEO','Māori','Pasifika','Asian'})
hold off

figure(6)
plot(t,timeseries_eth(:,(3*eth_sz+1):4*eth_sz),'LineWidth',1);
hold on
plot(t,timeseries_eth(:,(4*eth_sz+1):5*eth_sz),'LineWidth',1,'LineStyle','--');
xlabel('Time (in days)')
ylabel('Number of exposed people')
colororder(eth_cmap)
xlim tight
%legend({'NZEO','Māori','Pasifika','Asian'})
hold off

figure(7)
plot(t,timeseries_eth(:,(5*eth_sz+1):6*eth_sz),'LineWidth',1);
hold on
plot(t,timeseries_eth(:,(6*eth_sz+1):7*eth_sz),'LineWidth',1,'LineStyle','--');
xlabel('Time (in days)')
ylabel('Number of infectious people')
colororder(eth_cmap)
xlim tight
%legend({'NZEO','Māori','Pasifika','Asian'})
hold off


figure(8)
plot(t,timeseries_eth(:,(7*eth_sz+1):8*eth_sz),'LineWidth',1);
hold on
plot(t,timeseries_eth(:,(8*eth_sz+1):9*eth_sz),'LineWidth',1,'LineStyle','--');
xlabel('Time (in days)')
ylabel('Number of recovered people')
colororder(eth_cmap)
xlim tight
hold off

%% Attack rates

figure(9)
x1=age_vec;
bar(x1,total_prop_infd_ages)
ylabel('Attack Rate')

figure(10)
x2=["NZE/Other" "Māori" "Pasifika" "Asian"];
bar(x2,total_prop_infd_eth)
ylabel('Attack Rate')

figure(11)
plot(reshape(total_prop_infd,eth_sz,age_sz),'Linewidth',1,'Marker','o')
ylabel('Attack Rate')
xticks(1:eth_sz);
xticklabels({'NZE/Other','Māori','Pasifika','Asian'})
%legend(age_vec)
colororder(age_cmap)

figure(12)
plot(reshape(total_prop_infd,eth_sz,age_sz)','Linewidth',1,'Marker','o')
ylabel('Attack Rate')
xticks(1:age_sz);
xticklabels(age_vec)
legend('NZE/Other','Māori','Pasifika','Asian')
colororder(eth_cmap)

figure(13)
bar(["Susceptible/Vaccinated","Recovered"],infd_sus_matrix_ages,'stacked')
ylabel('Number of people')
legend(age_vec)
colororder(age_cmap)

figure(14)
bar(["Susceptible/Vaccinated","Recovered"],infd_sus_matrix_eth,'stacked')
ylabel('Number of people')
legend('NZE/Other','Māori','Pasifika','Asian')
colororder(eth_cmap)

figure(15)
xvalues={'NZE/Other','Māori','Pasifika','Asian'};
yvalues=age_vec;
h=heatmap(xvalues,yvalues,total_prop_infd_mat,'CellLabelColor','none');
h.XLabel='Ethnicity';
h.YLabel='Age';
h.Colormap=parula;
clim([0 1]);

disp(total_prop_infd_mat)
