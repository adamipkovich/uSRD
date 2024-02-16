clear all
close all
clc

%% Loading data
PerCapitaEmissions = readtable('Data\PerCapitaEmission.csv');
historicalemissions = readtable('Data\HistoricalEmission.csv');

poptab = readtable("Data\Pop2018", 'Sheet', "Data");
GDPtab = readtable("Data\GDP2018", 'Sheet', "Data");
RuralPoptab = readtable("Data\RuralPop2018", 'Sheet', "Data");
GDPGrowthtab = readtable("Data\GDPGrowth2018", 'Sheet', "Data");
UrbanGrowthtab = readtable("Data\UrbanGrowth2018", 'Sheet', "Data");

%% Renaming the variables
variables = table2array(unique(historicalemissions(:, 3)))
vars = {'AGR';'BLD'; 'BNK'; 'ELH'; 'ENG';'FUE'; 'IND'; 'LUCF'; 'MAN'; 'OFC'; 'SEC'; 'TOT'; 'TOTI'; 'TRP'; 'WAS'}
for j = 1:length(vars)
    index = strcmp(table2array(historicalemissions(:, 3)),variables(j));
    historicalemissions(index, 3) =  vars(j);
end
%% Pretreat data for SRD - Emissions
data=[historicalemissions(:, 1), historicalemissions(:, 3:4), historicalemissions(:, 6)]
data(1,:) = []
data = renamevars(data,["Var1", "Var3", "Var4" "Var6"],["CountryName", "Sector", "Gas", "Emission"]);

%data.Gas = cellfun(@(x) char(x),data.CountryName,'UniformOutput',false);
data.Gas = cellfun(@(x) strrep(x,"CH4","CH_4"), data.Gas,'UniformOutput',false);
data.Gas = cellfun(@(x) strrep(x,"N2O","N_2O"), data.Gas,'UniformOutput',false);
data.Gas = cellfun(@(x) strrep(x,"CO2","CO_2"), data.Gas,'UniformOutput',false);
data.Gas = cellfun(@(x) strrep(x,"F-Gas","FGas"), data.Gas,'UniformOutput',false);

for i=1:height(data)
    data.Sector(i) = {char(append(string(data.Sector(i)), 'e' + string(data.Gas(i))))};
end
data(:, 3) = [];
data= unstack(data,'Emission', 'Sector');

%% Pretreat data for SRD - Economic

economictab = [poptab(:, 1), poptab(:, 63)]
economictab= outerjoin(economictab, [GDPtab(:, 1), GDPtab(:, 63)],'Keys','CountryName', 'MergeKeys',true);
economictab= outerjoin(economictab, [RuralPoptab(:, 1), RuralPoptab(:, 63)],'Keys','CountryName', 'MergeKeys',true);
economictab= outerjoin(economictab, [GDPGrowthtab(:, 1), GDPGrowthtab(:, 63)],'Keys','CountryName', 'MergeKeys',true);
economictab= outerjoin(economictab, [UrbanGrowthtab(:, 1), UrbanGrowthtab(:, 63)],'Keys','CountryName', 'MergeKeys',true);
economictab = renamevars(economictab,["x2018_economictab","x2018_right","x2018_economictab_1", "x2018_right_1", "x2018"],["Pop","GDP","RuralPop", "GDPGrowth", "UrbanGrowth"]);
%% Join CAIT  and World Bank

climate_data= outerjoin(data, economictab, 'Keys','CountryName', 'MergeKeys',true);

%% Remove problematic rows and columns
toDelete = find(isnan(climate_data.Pop));
climate_data(toDelete,:) = [];
toDelete = find(isnan(climate_data.GDP));
climate_data(toDelete,:) = [];

%Remove columns with too many NaN

index = []
for i=2:(width(climate_data))
    if(sum(isnan(climate_data{: ,i})) > 120)
       index = [index, i] 
    end
end
climate_data(:, index) = [];
  

% Remove rows that have too much NaN

index = []
for i=1:height(climate_data)
    if(sum(isnan(climate_data{i,2:end-5})) > 10)
       index = [index, i] 
    end
end
climate_data(index,:) = [];

% Remove micronations with less than 500000 population
index = []
for i=1:height(climate_data)
    if(climate_data.Pop(i) < 500000)
       index = [index, i] 
    end
end
climate_data(index,:) = [];
variables = char(climate_data.Properties.VariableNames(2:end));

% Removing TOTI
climate_data.TOTIeAllGHG = [];
climate_data.TOTIeCO_2 = [];
climate_data.TOTIeN_2O = [];
climate_data.TOTIeCH_4 = [];
climate_data.TOTIeFGas = [];

climate_data.BLDeAllGHG = [];
climate_data.BNKeAllGHG = [];
climate_data.ELHeAllGHG = [];
climate_data.MANeAllGHG = [];
climate_data.TRPeAllGHG = [];

% Clearing variables

clear economictab
clear GDPGrowthtab
clear GDPtab
clear historicalemissions
clear poptab
clear RuralPoptab
clear UrbanGrowthtab
clear data

writetable(climate_data, 'climate_data.xlsx');

PerCapitaEmissions(1, :) = [];
PerCapitaEmissions(end-2: end, :) = [];
PerCapitaEmissions(end-2, :) = [];
PerCapitaEmissions(153, :) = [];
PerCapitaEmissions(:, 2:30) = [];
PerCapitaEmissions = renamevars(PerCapitaEmissions,["Var1", "Var31"],["CountryName","2018"]);

writetable(PerCapitaEmissions, 'mapData.xlsx');
clear PerCapitaEmissions
%% SRD
figure(1)
num = [table2array(climate_data(:, 2:36))./(climate_data.Pop), table2array(climate_data(:, 38))./(climate_data.Pop), table2array(climate_data(:, 39:end))];
gmax = num(:, 27);

u = [num(:, 1:26), num(:, 28:end)];
[N,n] = size(u);
R = tiedrank(u);


%gmax = num(:, 36)./num(:, 46);
[nrow,ncol]=size(R);
%max srd
    if rem(nrow,2)==1
        k=(nrow-1)/2;
        m=2*k*(k+1);
    else
        k=nrow/2;
        m=2*k^2;
    end
%the best "virtiual" method  / ideal objective / ideal ranking  
nrk=tiedrank(gmax, 'omitnan'); 
names= [string(climate_data.Properties.VariableNames(2:36)), string(climate_data.Properties.VariableNames(38:end))];
names(27) = [];
% Calculate the SRD
srd=sum(abs(R-repmat(nrk,1,n)),1, 'omitnan')/m*100;
[srdi,si]=sort(srd);
nsrdi=names(si); 
srdi;
% Simulation of the probabolity distribution 
nSim=1e4;
S=[];
io=[1:nrow];
for i=1:nSim
       S(i)=sum(abs(io-randperm(nrow)))/m*100;
end
%yyaxis left
for i=1:2
subplot(2, 1, i)
yyaxis left
[prob,srdc]= histcounts(S,unique(S),'Normalization','cdf');
plot([srdc],[0 prob] * 100)
xlabel('SRD')
ylabel('P(SRD) [%]')
hold on 
% Calculate the probabilities of the orderings (based on cum prob)
psrdi = interp1([0 srdc],[ 0 0 prob], srdi);
psrdi = psrdi +1;
yyaxis right
axv=axis;
ylabel('SRD')
%axis([0 2*max(srdi) 0 1.5*max(psrdi)])
for j=1:n 
an=text(srdi(j),srdi(j),nsrdi(j),'fontsize',16);
set(an,'Rotation',45);
if (i == 1 && srdi(j) > 30 && srdi(j) < 50)
delete(an)
end
if (i == 2 && (srdi(j) < 30 ||  srdi(j) > 50))
delete(an)
end
line([srdi(j) srdi(j) ],[0 srdi(j)], 'LineWidth', 2.5);
hold on 
end
set(gca,'FontSize',18)
end
srdmed = median(S);
xx1 =srdc(abs(prob-0.05) == min(abs(prob-0.05)))
xx19 = srdc(abs(prob-0.95) == min(abs(prob-0.95)))
subplot(2, 1, 1)
line([srdmed srdmed],[0 100], 'Color','black','LineStyle','--')
line([xx1, xx1],[0 100], 'Color','black','LineStyle','-.')
line([xx19, xx19],[0 100], 'Color','black','LineStyle','-.')
subplot(2,1,2)
line([srdmed srdmed],[0 100], 'Color','black','LineStyle','--')
line([xx1,xx1],[0 100], 'Color','black','LineStyle','-.')
line([xx19, xx19],[0 100], 'Color','black','LineStyle','-.')


natrank = []
nat = tiedrank(gmax)
[natsrdi, nnatsrdi] = sort(nat);
natrank = [natrank, climate_data.CountryName(nnatsrdi(1:nrow))];

natrank(nnatsrdi == 145) = []
%% Percentile transformation
num = [num(1:144, :); num(146:end, :)] 
nationnames = [climate_data.CountryName(1:144); climate_data.CountryName(146:end)]

[N,n] = size(num);
p = []

[up, yp] = ecdf(num(:, 1));
for i = 1:n
    [v_f,v_x, ia,ic] = homemade_ecdf(num(:, i)')
    p = [p, v_f(ic)']; 
end

%% Derringer desirability - ABC analysis limits

figure(100)
subplot(2, 1, 1)
[val, ind] = sort(num(:, 27))
plot(val*1000^2, p(ind, 27), '-')
xlim([0 30])
ylim([0 1])
ylabel("Percentile")
xlabel("Emission per capita [t/cap]")
set(gca,'FontSize',18)
subplot(2, 1, 2)
lwl = val(round(N*0.5), 1);
upl = val(round(N*0.8), 1);
u = derringer(val, upl,  lwl, 1);
plot(val*1000^2, u, '-')
xlim([0 30])
ylim([0 1])
ylabel("Percentile")
xlabel("Emission per capita [t/cap]")
set(gca,'FontSize',18)

%% Derringer desirability for all indicators

refemission = sum(num(:, 1:35), 1, 'omitnan')./sum(num(:, 27), 1, 'omitnan');
refemission(12) = [];
w = [refemission, ones(1, 4)];
u = [];
A = 0.8;
B = 0.5;
s = 1;

%l = derringer(sortedperCap, sortedperCap(round(N*A), 1)*refemission(i), sortedperCap(round(N*B), 1)* refemission(i), s);
for i=1:11
[sortedperCap, sortindex] = sort(num(:, i));
l = derringer(num(:, i), sortedperCap(round(N*A), 1), sortedperCap(round(N*B), 1), s);
u =[u l]; 
end
for i=13:35
[sortedperCap, sortindex] = sort(num(:, i));
l = derringer(num(:, i), sortedperCap(round(N*A), 1), sortedperCap(round(N*B), 1), s);
u =[u l]; 
end

[sortedperCap, sortindex] = sort(1-normalize(num(:, 36), 'range'));
l = derringer(sortedperCap, sortedperCap(round(N*A), 1),  sortedperCap(round(N*B), 1), s);
u =[u l(sortindex)]; 
[sortedperCap, sortindex] = sort(num(:, 37));
l = derringer(sortedperCap, sortedperCap(round(N*A), 1),  sortedperCap(round(N*B), 1), s);
u =[u l(sortindex)]; 
[sortedperCap, sortindex] = sort(num(:, 38));
l = derringer(sortedperCap, sortedperCap(round(N*A), 1),  sortedperCap(round(N*B), 1), s);
u =[u l(sortindex)]; 
[sortedperCap, sortindex] = sort(num(:, 39));
l = derringer(sortedperCap, sortedperCap(round(N*A), 1),  sortedperCap(round(N*B), 1), s);
u =[u l(sortindex)]; 


%% Derringer SRD
figure(4)

R = tiedrank(u);
gmax = sum(u.*w, 2, 'omitnan');


[nrow,ncol]=size(R);
    if rem(nrow,2)==1
        k=(nrow-1)/2;
        m=2*k*(k+1);
    else
        k=nrow/2;
        m=2*k^2;
    end
%the best "virtiual" method  / ideal objective / ideal ranking  
nrk=tiedrank(gmax, 'omitnan'); 
names = [string(climate_data.Properties.VariableNames(2:12)), string(climate_data.Properties.VariableNames(14:36)), string(climate_data.Properties.VariableNames(38:end))];
% Calculate the SRD
srd=sum(abs(R-repmat(nrk,1,ncol)),1, 'omitnan')/m*100;
[srdi,si]=sort(srd);
nsrdi=names(si); 
srdi;
% Simulation of the probabolity distribution 
nSim=1e4;
S=[];
io=[1:nrow];
for i=1:nSim
       S(i)=sum(abs(io-randperm(nrow)))/m*100;
end
%yyaxis left
for i=1:2
subplot(2, 1, i)
yyaxis left
[prob,srdc]= histcounts(S,unique(S),'Normalization','cdf');
plot([srdc],[0 prob] * 100)
xlabel('SRD')
ylabel('P(SRD) [%]')
hold on 
% Calculate the probabilities of the orderings (based on cum prob)
psrdi = interp1([0 srdc],[ 0 0 prob], srdi);
psrdi = psrdi +1;
yyaxis right
axv=axis;
ylabel('SRD')
%axis([0 2*max(srdi) 0 1.5*max(psrdi)])
for j=1:ncol 
an=text(srdi(j),srdi(j),nsrdi(j),'fontsize',16);
set(an,'Rotation',45);
if (i == 2 && srdi(j) > 10 && srdi(j) < 50)
delete(an)
end
if (i == 1 && (srdi(j) >= 50 ))
delete(an)
end
line([srdi(j) srdi(j) ],[0 srdi(j)], 'LineWidth', 2.5);
hold on 
end
set(gca,'FontSize',18)
end
srdmed = median(S);
xx1 =srdc(abs(prob-0.05) == min(abs(prob-0.05)));
xx19 = srdc(abs(prob-0.95) == min(abs(prob-0.95)));
subplot(2, 1, 1)
line([srdmed(1) srdmed(1)],[0 100], 'Color','black','LineStyle','--');
line([xx1, xx1],[0 100], 'Color','black','LineStyle','-.');
line([xx19, xx19],[0 100], 'Color','black','LineStyle','-.');
subplot(2,1,2)
line([srdmed(1) srdmed(1)],[0 100], 'Color','black','LineStyle','--');
line([xx1,xx1],[0 100], 'Color','black','LineStyle','-.');
line([xx19, xx19],[0 100], 'Color','black','LineStyle','-.');
%% Ranking of nations based on Derringer
pc = 1 - normalize(num(:, end-3), 'range');
[sortedperCap, sortindex] = sort(pc);
l = derringer(pc, sortedperCap(round(N*A), 1),  sortedperCap(round(N*B), 1), s);
u(:, end-3) = l; 


pc = normalize(num(:, end-2), 'range')
[sortedperCap, sortindex] = sort(pc);
l = derringer(pc, sortedperCap(round(N*A), 1),  sortedperCap(round(N*B), 1), s);
u(:, end-2) = l; 

pc = num(:, end-1)
pc (pc < 0) = 0
pc = 1 - normalize( pc, 'range')
[sortedperCap, ~] = sort(pc);
l = derringer(pc, sortedperCap(round(N*A), 1),  sortedperCap(round(N*B), 1), s);
u(:, end-1) = l; 

pc = num(:, end)
pc (pc < 0) = 0
pc = 1 - normalize(pc, 'range')
[sortedperCap, sortindex] = sort(pc);
l = derringer(pc, sortedperCap(round(N*A), 1),  sortedperCap(round(N*B), 1), s);
u(:, end) = l; 


ru = [u(:, 27:30), u(:, end-3:end)] 
wu = [w(27:30), w(end-3:end)]
%w(26:31) = []
gmax = sum(ru.*wu, 2, 'omitnan');

[sSUD, sSUDi] = sort(gmax);
sru = round(ru(sSUDi, :), 4)
nSUD = round(sSUD./sum(wu)*100, 4) %Normalize with the maximum desirability
nSUDName = nationnames

%% Mitigation assessment
mit_vars = climate_data.Properties.VariableNames
n_vars = []
for i =1:length(mit_vars)
    var_str = char(mit_vars(i));
    
    if length(strfind(var_str,'AllGHG')) ~= 0 | (length(strfind(var_str,'CO_2')) ~= 0 &  length(strfind(var_str,'BLD')) ~= 0 | length(strfind(var_str,'ELH')) ~= 0 | length(strfind(var_str,'BNK')) ~= 0 |length(strfind(var_str,'MAN')) ~= 0 |length(strfind(var_str,'TRP')) ~= 0)
        n_vars = [n_vars, i];
    elseif length(strfind(var_str,'FUEeCH_4')) ~= 0 
        n_vars = [n_vars, i];
    end

end

mt_names = mit_vars(n_vars)
mit_dat = climate_data(:, n_vars)
mdr = tiedrank(table2array(mit_dat))

gmax = table2array(mit_dat(:, end-2));
nrk=tiedrank(gmax, 'omitnan'); 

mit_srd=sum(abs(mdr-repmat(nrk,1,size(mit_dat, 2))),1, 'omitnan')/m*100;
[msrdi,msi]=sort(mit_srd);
mnsrdi=mt_names(msi); 

[mit_act, mit_act_labels, ma_table] = xlsread("mitigation_actions.xlsx")
mit_act(:, 1:4) = []
mit_act(isnan(mit_act))=0


%D = pdist(mit_act_ord, 'jaccard');
%leafOrder = optimalleaforder(tree,D)

figure(3031)
tree = linkage(mit_act, 'complete', 'jaccard')
[~,~,dind] = dendrogram(tree,0)

figure(3032)
tree2 = linkage(mit_act', 'single', 'jaccard')
[~,~,dind2] = dendrogram(tree2)

mit_act_ord = mit_act(:, msi)
mit_act_lord = mit_act_labels(2:end, 4)
%mit_act_lord = mit_act_lord(dind, :)
figure(3030)
h = heatmap(mit_act_ord(dind, :), 'CellLabelColor','none') %, "ColumnLabels"
h.XDisplayLabels = mt_names(msi)
h.YDisplayLabels = mit_act_lord(dind)  %mit_act_lord
h.ColorbarVisible = 'off'


figure(3033)
h = heatmap(mit_act(dind,dind2), 'CellLabelColor','none') %, "ColumnLabels"
h.XDisplayLabels = mt_names(dind2)
%h.XDisplayLabels = mt_names(msi)
h.YDisplayLabels = mit_act_lord(dind)
h.ColorbarVisible = 'off'

%% Mitigation actions to sectors to countries .....

%Get utility, co2, ch4, n2o
avar_names = names(:, 1:end-4) 
n_vars = []
for i =1:length(avar_names)
    var_str = char(avar_names(i));  
    if length(strfind(var_str,'AllGHG')) == 0
        n_vars = [n_vars, i];
    end
end
avar_names(n_vars)

%%
% matrix of the utilities for *all* variables, even those that were not
% available
% 146 countries, and 3*12 criteria
U = zeros(146, 36) %% countries - sectors+gases
S = zeros(36,12) %% gasses - sectors
d_w = zeros(36, 1)
% generate labels
ds_lab = {'AGR';'BLD'; 'BNK'; 'ELH'; 'ENG';'FUE'; 'IND'; 'LUCF'; 'MAN'; 'OFC'; 'TRP'; 'WAS'}
ds_gases = {'CO_2', 'CH_4', 'N_2O'}
nat_mit_labels = []
for i =1:length(ds_lab)
   for j = 1:length(ds_gases)
        
        nat_mit_labels =[nat_mit_labels,  (string(ds_lab(i)) + "e" + string(ds_gases(j)))]
        idx = []
        for k = 1:length(names(1, 1:end-4))
            idx = strfind(names(k), nat_mit_labels(end));
            if length(idx) ~= 0
                break;
            end
        end
        if length(idx) ~= 0
            U(:,(i-1)*3+j) = u(:, k); 
            d_w((i-1)*3+j, 1) = w(:, k);
        else
            U(:,(i-1)*3+j) = 1;
            d_w((i-1)*3+j, 1) = 1;
        end
        
        S((i-1)*3+j,i) = 1; 
        
   end   
end

U(isnan(U)) = 1

%%

M = mit_act' %% sectors - mitigation actions
M(11, :) = []

mit_labels = mit_act_labels(2:end, 4)

%mitigation action decision support matrix:
DSM =((repmat(d_w', size(U, 1), 1).*(1-(U)))*S/3)*(M./repmat(sum(M),size(M,1),1));
%heatmap(DSM)
cg = clustergram((DSM -0.04), 'Linkage','average', 'RowLabels', nationnames, 'ColumnLabels', mit_labels ,  'RowPdist', 'spearman', 'ColumnPdist', 'spearman'); %, 'Dendrogram',0.025
cg.Colormap = redbluecmap;

[vd, didx] =  min(DSM, [], 2);
[dum1, dum2] = unique(didx)


rDSM = DSM([sSUDi(1:10);flip(sSUDi((end-9):end))], :) %mit_labels
cg = clustergram((rDSM' -0.04)', 'Linkage','average', 'ColumnLabels',string(dind) , 'RowLabels',   nationnames([sSUDi(1:10);flip(sSUDi((end-9):end))]),  'RowPdist', 'spearman', 'ColumnPdist', 'spearman', 'DisplayRatio', 0.1); %, 'Dendrogram',0.025
cg.Colormap = redbluecmap;
figureHandle = gcf;
%# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',16,'fontWeight','bold')
% Clustergram of the ranking



mit_labels()

for i =1:length(unique(didx))
    dummy = unique(didx);
    DSM(didx == dummy(i) ,dummy(i) ) = 1;
end

[uvd, udidx] =  min(DSM, [], 2);
mit_labels(unique(udidx))


