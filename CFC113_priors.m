function [longProdsamps,shortProdsamps,RFlongBankLT,DElong, DEshort,Years] = CFC113_priors(N)

y1 = 1955; % Year Start Date
yswitch = 1988;
y2 = 2018; % Year End Date
Sim_yrs = y1:y2;


Prod_reported = NaN(size(Sim_yrs)); % Production is in Tonnes
load('WMO2002.mat') % Use until 1988 - WMO data is the adjusted AFEAS data
ytmp1 = find(WMO2002year==y1); 
ytmp2 = find(WMO2002year==yswitch);
Prod_reported(1:yswitch+1-y1) = Production(ytmp1:ytmp2);

load('a5_na5.mat') % Use Article 5 and non A5 prod total from 1989 onwards
ytmp1 = find(a5_na5(:,1) == yswitch+1);
ytmp2 = find(Sim_yrs == a5_na5(end,1));
Prod_reported(1989-y1+1:ytmp2) = sum(a5_na5(ytmp1:end,[2:3]),2);
Prod_reported(ytmp2:end) = sum(a5_na5(end,[2:3]),2);

str = strcat('AFEAS_cfc113production.mat');
load(str,'longbank','shortbank','yr'); % in thousands of tonnes

ind = find(yr == y1); 
add_yrs = [yr(end)+1:Sim_yrs(end)]'; 

tmp_s = shortbank(ind:end) - shortbank(ind-1:end-1);
tmp_l = longbank(ind:end) - longbank(ind-1:end-1);

tmp_s = [tmp_s; tmp_s(end)*ones(size(add_yrs))];
tmp_l = [tmp_l; tmp_l(end)*ones(size(add_yrs))];
tmp_total = tmp_s + tmp_l;

prod_s = (tmp_s./(tmp_total)).*Prod_reported'; 
prod_l = (tmp_l./(tmp_total)).*Prod_reported'; 

CFC113ReportedProd.BankedProd = Prod_reported;
CFC113ReportedProd.yrs = y1:y2;
%save('Prior_Production.mat','CFC113ReportedProd','-append')



% creating production distribution
for ii = 1:2
    
    switch ii
        case 1 
            Prod_tmp = prod_s';
        case 2
            Prod_tmp = prod_l';    
    end
    
   % creating production distribution
    Nyrs = length(y1:yswitch);
    rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
    rho_tmp = repmat(rho,Nyrs,Nyrs,1);
    exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
    Rhom = rho_tmp.^exp_val;
    clear rho_tmp exp_val
    Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
    clear Rhom
    Um = normcdf(Zm,0,1);
    clear Zm
    Prod_samps(:,1:Nyrs) = (repmat(0.3*Prod_tmp(1:Nyrs),N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.70*repmat(Prod_tmp(1:Nyrs),N,1);
    clear Um

    Nyrs2 = length(Sim_yrs) - Nyrs; 
    rho_tmp = repmat(rho,Nyrs2,Nyrs2,1);

    exp_val = repmat(abs(repmat([1:Nyrs2],Nyrs2,1)-repmat([1:Nyrs2]',1,Nyrs2)),1,1,N);
    Rhom = rho_tmp.^exp_val;
    clear rho_tmp exp_val
    Zm = mvnrnd(zeros(N,Nyrs2), Rhom, N);
    clear Rhom
    Um = normcdf(Zm,0,1);
    clear Zm
    Prod_samps(:,Nyrs+1:length(Sim_yrs)) = (repmat(0.3*Prod_tmp(Nyrs+1:end),N,1)).*logninv(Um,zeros(N,Nyrs2),0.5*ones(N,Nyrs2))+0.70*repmat(Prod_tmp(Nyrs+1:end),N,1);
    clear Um

    switch ii
        case 1 
            shortProdsamps = Prod_samps;
        case 2         
            longProdsamps = Prod_samps;   
    end
end

DEshort = betarnd(12,12,1,N);
RFshort = 1-DEshort;

m = 0.02; v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
DElong= lognrnd(mu,sigma,1,N); %0.02*DEu(2,:);

m = 20; v = (0.5*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFlongBankLT = lognrnd(mu,sigma,1,N); 


Years = 1955:2018;
%save('CFC113Priors.mat','Years','longProdsamps','shortProdsamps','RFlongBankLT','DElong','DEshort','RFshort','Prod_reported','RFfeedstock', 'Feedstocksamps');

