%% Oxygen utilisation rate
%July 13th, 2021
%Code by Olivier Sulpis

clear all

%% Load data files

%load raw 2016 GLODAP dataset (Olsen et al., 2016, ESSD)
load('glodap_reduced1.mat');
load('glodap_reduced2.mat');

%G2rho=gsw_rho(G2salinity,G2temperature,G2pressure);

%load Biomes from Fay and McKinley (2014 ESSD) and their 2D and 3D mapped versions
load('G2biomes2016_new.mat')
load('GLODAP_map_biomes_new.mat')
%load('GLODAP_map_biomes_3d_new.mat')

%load chlorophyll and compute Ez
%Ez=depth at which the downward PAR irradiance falls to 1% of its
%sub-surface value (Morel et al., RSoE2007)
%Chl data from [Copernicus marine service,
%GLOBAL-ANALYSIS-FORECAST-BIO-001-028-MONTHLY, surface data averaged between May 2020 and May 2021] 
load('copernicus_chl.mat')

%load seawater age from TTD analysis (from Emil Jeansson)
load('GLODAPv2_TTDage.mat');

%load seawater age from 14C (from Gebbie and Huybers 2012, JPO, shared by Adam Subhas)
load('GHage.mat');
%GH_tmean=ncread('/Users/oliviersulpis/Library/CloudStorage/GoogleDrive-olivier.sulpis@gmail.com/My Drive/Research/Dark Respiration/Datasets/Age_2deg_GH2012.nc','tmean');
%GH_lat=ncread('/Users/oliviersulpis/Library/CloudStorage/GoogleDrive-olivier.sulpis@gmail.com/My Drive/Research/Dark Respiration/Datasets/Age_2deg_GH2012.nc','latitude');
%GH_depth=ncread('/Users/oliviersulpis/Library/CloudStorage/GoogleDrive-olivier.sulpis@gmail.com/My Drive/Research/Dark Respiration/Datasets/Age_2deg_GH2012.nc','depth');
%GH_lon=ncread('/Users/oliviersulpis/Library/CloudStorage/GoogleDrive-olivier.sulpis@gmail.com/My Drive/Research/Dark Respiration/Datasets/Age_2deg_GH2012.nc','longitude');

%make longitudes monotonic
for i=1:size(G2longitude)
    if G2longitude(i)<0
       G2longitude_monotonic(i,1)=G2longitude(i)+360;
    else
        G2longitude_monotonic(i,1)=G2longitude(i);
    end
end

%load oxygen data including preformed concentrations (from Mark Holzer)
%load('O2_obs_pre_GD24L_CTR.mat');
%TOU=O2_pre_glodap_3d-O2_obs_glodap_3d;
%G2_TOU(:,1)=interp3(xt,yt,zt,TOU,G2longitude_monotonic,G2latitude,G2depth,'linear');
%AGE(:,1)=interp3(xt,yt,zt,Gamma_3d,G2longitude_monotonic,G2latitude,G2depth,'linear');
%load('O2_obs_pre_GD24L_KiLOW.mat');
%TOU=O2_pre_glodap_3d-O2_obs_glodap_3d;
%G2_TOU(:,2)=interp3(xt,yt,zt,TOU,G2longitude_monotonic,G2latitude,G2depth,'linear');
%AGE(:,2)=interp3(xt,yt,zt,Gamma_3d,G2longitude_monotonic,G2latitude,G2depth,'linear');
%load('O2_obs_pre_GD24L_KvHIGH.mat');
%TOU=O2_pre_glodap_3d-O2_obs_glodap_3d;
%G2_TOU(:,3)=interp3(xt,yt,zt,TOU,G2longitude_monotonic,G2latitude,G2depth,'linear');
%AGE(:,3)=interp3(xt,yt,zt,Gamma_3d,G2longitude_monotonic,G2latitude,G2depth,'linear');
%load('O2_obs_pre_GD24L_KvHigh_KiLOW.mat');
%TOU=O2_pre_glodap_3d-O2_obs_glodap_3d;
%G2_TOU(:,4)=interp3(xt,yt,zt,TOU,G2longitude_monotonic,G2latitude,G2depth,'linear');
%AGE(:,4)=interp3(xt,yt,zt,Gamma_3d,G2longitude_monotonic,G2latitude,G2depth,'linear');
%load('O2_obs_pre_GD48L_CTR.mat');
%TOU=O2_pre_glodap_3d-O2_obs_glodap_3d;
%G2_TOU(:,5)=interp3(xt,yt,zt,TOU,G2longitude_monotonic,G2latitude,G2depth,'linear');
%AGE(:,5)=interp3(xt,yt,zt,Gamma_3d,G2longitude_monotonic,G2latitude,G2depth,'linear');
%load('O2_obs_pre_GD48L_Kvhalf.mat');
%TOU=O2_pre_glodap_3d-O2_obs_glodap_3d;
%G2_TOU(:,6)=interp3(xt,yt,zt,TOU,G2longitude_monotonic,G2latitude,G2depth,'linear');
%AGE(:,6)=interp3(xt,yt,zt,Gamma_3d,G2longitude_monotonic,G2latitude,G2depth,'linear');
%load('O2_obs_pre_GD48L_Kvx2.mat');
%TOU=O2_pre_glodap_3d-O2_obs_glodap_3d;
%G2_TOU(:,7)=interp3(xt,yt,zt,TOU,G2longitude_monotonic,G2latitude,G2depth,'linear');
%AGE(:,7)=interp3(xt,yt,zt,Gamma_3d,G2longitude_monotonic,G2latitude,G2depth,'linear');

%load oxygen data including preformed concentrations (from Mark Holzer)
load('O2_obs_pre_GD.mat');
load('O2_obs_pre_GD24L_CTR.mat','Gamma_3d');


%% Find Ez for each biome

Ez=10.^(1.524-0.436.*log10(chl)-0.0145.*log10(chl).^2+0.0186.*log10(chl).^3); %in meters

for i=1:1440
    if chl_lon(i)<0
       chl_lon(i,1)=chl_lon(i,1)+360;
    else
        chl_lon(i,1)=chl_lon(i,1);
    end
end

for i=1:1440
    chl_lat(:,i)=chl_lat(:,1);
end
chl_lat=chl_lat';

for i=1:681
    chl_lon(:,i)=chl_lon(:,1);
end

%assign a biome number to each sediment-trap data: interpolation to the nearest neighbour
chl_biomes=interp2(lon,lat,squeeze(G2map_Biomes(:,:,1)),chl_lon,chl_lat,'nearest');

for i=1:10
    Ez_biome(i)=nanmean(nanmean(Ez(chl_biomes==i)));
end

%% Compute seawater age

%assign an 14C age to each GLODAP sample: linear interpolation
C14_age=interp3(GH_lat,GH_depth,GH_lon,GH_tmean,G2latitude,G2depth,G2longitude_monotonic,'linear');

%G2_age is the age that will be used in the analysis. 
% First, we state that the age of each GLODAP sample is that from the TTD analysis.
TTD_age=mageF12;
G2_age=TTD_age;

 %Then we create a transition zone in each the used age comes from a mixture of TTD and 14C
flagage=ones(size(G2_age,1),1);
flagage(isnan(G2_age))=NaN;
for i=1:size(TTD_age)
    if TTD_age(i) >= 200 && TTD_age(i) <= 300  %transition period between 200 and 300 years
        transition_coef(i) = (TTD_age(i)-200)./100;  %linear transition
        G2_age(i) = (1-transition_coef(i))*TTD_age(i) + transition_coef(i)*C14_age(i);
        flagage(i)=1+transition_coef(i);
    elseif TTD_age(i) > 300 %Finally, above 300 years, we use 14C ages only 
        G2_age(i)=C14_age(i);
        flagage(i)=2;
    end 
end

%% Compute preformed oxygen

%Compute total oxygen utilization using preformed oxygen provided by Mark Holzer 
TOU=O2_pre_glodap_3d-O2_obs_glodap_3d;
G2_TOU=interp3(xt,yt,zt,TOU,G2longitude_monotonic,G2latitude,G2depth,'linear');
G2_OCIMage=interp3(xt,yt,zt,Gamma_3d,G2longitude_monotonic,G2latitude,G2depth,'linear');


%% Compute CaCO3 dissolution rates

%"1" to have bins centered around certain sigma0 values, "0" to have bins with constant nb of data pts
binningmethod = 1;

%find what the biome numbers are 
G2_b=unique(G2biomes(~isnan(G2biomes)))';

%Set up MonteCarlo
MC=20; %number of Monte Carlo simulations
minimum_sigma_increment=0.0005;
maximum_sigma_increment=0.01;
O2_uncertainty=0.1; % 10% standard relative uncertainty for GLODAP oxygen
age_uncertainty=0.2; % 20%, standard relative uncertainty for CFC and 14C ages
minimum_sigma0=20;
maximum_sigma0=28.78;
max_sigma0bins=[minimum_sigma0:minimum_sigma_increment:maximum_sigma0]; 

%Results will be stored here: 1st dimension is Monte Carlo simulation #, 2nd is biome #, 3rd dimension in density bin# 
O2_rate=NaN(MC,size(G2_b,2),size(max_sigma0bins,2)); 
O2_rate_std=NaN(MC,size(G2_b,2),size(max_sigma0bins,2)); 
O2_rate_depth=NaN(MC,size(G2_b,2),size(max_sigma0bins,2));
O2_rate_maxdepth=NaN(MC,size(G2_b,2),size(max_sigma0bins,2));
O2_rate_mindepth=NaN(MC,size(G2_b,2),size(max_sigma0bins,2));
O2_rate_std_depth=NaN(MC,size(G2_b,2),size(max_sigma0bins,2));
O2_rate_r2=NaN(MC,size(G2_b,2),size(max_sigma0bins,2));
Nbofpoints=NaN(MC,size(G2_b,2),size(max_sigma0bins,2));
O2_range=NaN(MC,size(G2_b,2),size(max_sigma0bins,2));
Age_range=NaN(MC,size(G2_b,2),size(max_sigma0bins,2));
Age_source=NaN(MC,size(G2_b,2),size(max_sigma0bins,2));
p_value=NaN(MC,size(G2_b,2),size(max_sigma0bins,2));

%FIRST FOR-LOOP Monte Carlo for-loop
for h=1:MC 
h
    %O2, age, and sigma increment will randomly change at each MC simulation
%    erraou=G2aou.*2.*O2_uncertainty.*rand(size(G2_age,1),1);
    err_TOU=G2_TOU.*2.*O2_uncertainty.*rand(size(G2_age,1),1);
%    err_TOU=G2doc./1.027.*2.*O2_uncertainty.*rand(size(G2_age,1),1);
    errO2=G2oxygen.*2.*O2_uncertainty.*rand(size(G2_age,1),1);
    errage=G2_age.*2.*age_uncertainty.*rand(size(G2_age,1),1);
    sigma_increment=(maximum_sigma_increment-minimum_sigma_increment)*rand(1)+minimum_sigma_increment;
    %"random" values
    rG2_age=G2_age-G2_age.*age_uncertainty+errage;    
    rG2_age(rG2_age<0)=NaN;    
    rO2=G2_TOU-G2_TOU.*O2_uncertainty+err_TOU;
%    rO2=G2doc./1.027-G2doc./1.023.*O2_uncertainty+err_TOU;

        %THIRD FOR-LOOP biome for-loop
        for i=1:size(G2_b,2) 

            bsigma0=G2gamma(G2biomes==G2_b(i));  %densities in a given biome
            bO2=rO2(G2biomes==G2_b(i));                    %oxygen in a given biome
            bage=rG2_age(G2biomes==G2_b(i));            %seawater age in a given biome
            bflagage=flagage(G2biomes==G2_b(i));       %age source in a given biome, '1' is pure CFC-based, '2' is pure 14C-based
            bdepth=G2depth(G2biomes==G2_b(i));        %depths in a given biome
    
            [~,a]=sort(bsigma0,'descend');  %a is a dummy vector that sorts data according to descending sigma0
            bO2=bO2(a,:);                            %oxygen is sorted according to decreasing density
            bage=bage(a,:);                          %sw age is sorted according to decreasing density
            bflagage=bflagage(a,:);              %sw age source is sorted according to decreasing density
            bsigma0=bsigma0(a,:);               %sigma0 is sorted according to decreasing density
            bdepth=bdepth(a,:);                   %depth is sorted according to decreasing density
            clear a
    
            a=~isnan(bO2) & ~isnan(bage) & ~isnan(bsigma0) & ~isnan(bdepth) & bdepth>100 & bage>0; 
            %a is a dummy vector that selects samples for which we have Alkstar, age, sigma, depth, positive age and below 100 m depth
            bO2=bO2(a==1);
            bage=bage(a==1);
            bflagage=bflagage(a==1);
            bsigma0=bsigma0(a==1);
            bdepth=bdepth(a==1);
    
            sigma0bins=[minimum_sigma0:sigma_increment:maximum_sigma0]; 
    
            bin_nb=1;  %bin_nb is the number of the bin. We start at bin_nb = 1
        
            %FOURTH FOR-LOOP density bin for-loop
            for j=1:size(sigma0bins,2) 
        
                bsigma0min=sigma0bins(j)-sigma_increment/2;                                              %minimum sigma0 of the density bin
                bsigma0max=sigma0bins(j)+sigma_increment/2;                                            %maximum sigma0 of the density bin
                bin_age=bage(bsigma0<=bsigma0max & bsigma0>=bsigma0min);             %seawater ages within bin
                bin_flagage=bflagage(bsigma0<=bsigma0max & bsigma0>=bsigma0min); %seawater age source within bin
                bin_O2=bO2(bsigma0<=bsigma0max & bsigma0>=bsigma0min);                %oxygen within bin
                bin_sigma0=bsigma0(bsigma0<=bsigma0max & bsigma0>=bsigma0min);  %denities within bin
                bin_depth=bdepth(bsigma0<=bsigma0max & bsigma0>=bsigma0min);      %depths within bin
                                
                [~,outliers]=rmoutliers(bin_age);
                bin_age(outliers==1)=NaN;
                %[~,outliers]=rmoutliers(bin_O2);
                %bin_O2(outliers==1)=NaN;  
                
                bin_nbsp=sum(~isnan(bin_age));                                                                            %number of data points within bin
                
                if bin_nbsp>2 && (max(bin_age)-min(bin_age))>0.1*nanmean(bin_age) && (max(bin_O2)-min(bin_O2))>0.01*nanmean(bin_O2) && (max(bin_depth)-min(bin_depth))<(0.3864*mean(bin_depth)+261.36)
%                if (max(bin_age)-min(bin_age))>0.1*nanmean(bin_age) && (max(bin_O2)-min(bin_O2))>0.01*nanmean(bin_O2) && (max(bin_depth)-min(bin_depth))<(0.3864*mean(bin_depth)+261.36)

                %[b,b_int,~,~,stats]=regress(bin_O2,[ones(size(bin_age)) bin_age],0.1);
                mdl=fitlm(bin_age,bin_O2,'linear','RobustOpts','talwar','Weights',1./bin_age);
                %change 'talwar' to 'off' to remove robust
                %robust linear regression to reduce effects of outliers https://nl.mathworks.com/help/stats/robust-regression-reduce-outlier-effects.html
                
                b(1)=table2array(mdl.Coefficients(1,1)); %intercept
                b(2)=table2array(mdl.Coefficients(2,1)); %slope
                b(3)=table2array(mdl.Coefficients(1,2)); %standard error in the intercept
                b(4)=table2array(mdl.Coefficients(2,2)); %standard error in the slope
                pearson=cell2mat(struct2cell(mdl.Rsquared)); %R squared of the linear fit - cell format
                stats(1)=pearson(1,1); %R squared of the linear fit - array format                
                stats(3)=table2array(mdl.Coefficients(2,4)); %p value associated with the slope 
                
                    if stats(3)<1
                        
                        %plot randomly selected O2 versus age plots
                        %random_figure=rand();
                        %if random_figure<0.01
                        %    figure
                        %    plot([round(min(bin_age)) round(max(bin_age))],[round(min(bin_age)) round(max(bin_age))].*b(2)+b(1),'k-','LineWidth',1.5) %plot the linear fit
                        %    hold on
                        %    plot([round(min(bin_age)) round(max(bin_age))],[round(min(bin_age)) round(max(bin_age))].*(b(2)-b(4))+b(1)-b(3),'k-','LineWidth',1) %linear fit lower bound
                        %    plot([round(min(bin_age)) round(max(bin_age))],[round(min(bin_age)) round(max(bin_age))].*(b(2)+b(4))+b(1)+b(3),'k-','LineWidth',1) %linear fit upper bound
                        %    scatter(bin_age,bin_O2,30,bin_depth,'filled','MarkerEdgeColor',[0 0 0]) %scatter plot of O2 versus age, depth in color
                        %    cmocean('deep')
                        %    colorbar
                        %end
                        
                    %binned_fits is the array where the results are saved
                    binned_fits(i,bin_nb,1)=b(1);                           %intercept of regression oxygen versus age (umol / kg / a)
                    binned_fits(i,bin_nb,2)=b(2);                           %slope of regression oxygen versus age (umol / kg / a)
                    binned_fits(i,bin_nb,3)=stats(1);                     %r squared of linear regression 
                    binned_fits(i,bin_nb,4)=mean(bin_depth);      %mean depth of samples within a bin
                    binned_fits(i,bin_nb,5)=std(bin_depth);          %standard deviation of depths within a bin
                    binned_fits(i,bin_nb,6)=min(bin_depth);         %minimum depth within a bin
                    binned_fits(i,bin_nb,7)=max(bin_depth);        %maximum depth within a bin
                    binned_fits(i,bin_nb,8)=mean(bin_sigma0);   %mean sigma0 within a bin
                    binned_fits(i,bin_nb,9)=std(bin_sigma0);       %standard deviation of sigma0 within a bin
                    binned_fits(i,bin_nb,10)=min(bin_sigma0);    %minimum sigma0 within a bin
                    binned_fits(i,bin_nb,11)=max(bin_sigma0);   %maximum sigma0 within a bin
                    binned_fits(i,bin_nb,12)=bin_nbsp;                 %number of data points within a bin
                    binned_fits(i,bin_nb,13)=mean(bin_O2);        %mean oxygen 
                    binned_fits(i,bin_nb,14)=std(bin_O2);            %standard deviation of oxygen
                    binned_fits(i,bin_nb,15)=min(bin_O2);           %minimum oxygen
                    binned_fits(i,bin_nb,16)=max(bin_O2);          %maximum oxygen
                    binned_fits(i,bin_nb,17)=mean(bin_age);       %mean age
                    binned_fits(i,bin_nb,18)=std(bin_age);           %standard deviation of age
                    binned_fits(i,bin_nb,19)=min(bin_age);          %minimum age
                    binned_fits(i,bin_nb,20)=max(bin_age);         %maximum age 
                    binned_fits(i,bin_nb,21)=nanmean(bin_flagage);
                    binned_fits(i,bin_nb,22)=stats(3);
                    binned_fits(i,bin_nb,23)=b(4);         %standard error of the slope
                    else
                    binned_fits(i,bin_nb,1:23)=NaN;
                    end
                                
                else
                binned_fits(i,bin_nb,1:23)=NaN;
                end
                                
                bin_nb=bin_nb+1;
                            
            end     %FOURTH FOR-LOOP density bin for-loop
            
%[nanmax(binned_fits(i,:,4));nanmax(binned_fits(i,:,7))]

        end    %THIRD FOR-LOOP biome for-loop

    binned_fits(binned_fits==0)=NaN; 
  
    %important results are stored separately
    %1st dimension is Monte Carlo #, #2nd dimension is biome #, 3rd dimension is bin #
    O2_rate(h,:,1:size(binned_fits,2))=binned_fits(:,:,2);                                        %oxygen respiration rate in [umol/kg/a]
    O2_rate_std(h,:,1:size(binned_fits,2))=binned_fits(:,:,23);                               %oxygen respiration rate standard error in [umol/kg/a]
    O2_rate_depth(h,:,1:size(binned_fits,2))=binned_fits(:,:,4);                            %mean depth for oxygen respiration rate
    O2_rate_maxdepth(h,:,1:size(binned_fits,2))=binned_fits(:,:,7);                     %max depth for oxygen respiration rate
    O2_rate_mindepth(h,:,1:size(binned_fits,2))=binned_fits(:,:,6);                      %min depth for oxygen respiration rate
    O2_rate_std_depth(h,:,1:size(binned_fits,2))=binned_fits(:,:,5);                     %std dev of depth for oxygen respiration rate
    O2_rate_r2(h,:,1:size(binned_fits,2))=binned_fits(:,:,3);                                  %R2 of oxygen versus age fit
    Nbofpoints(h,:,1:size(binned_fits,2))=binned_fits(:,:,12);                               %nb of points for bin
    O2_range(h,:,1:size(binned_fits,2))=binned_fits(:,:,16)-binned_fits(:,:,15);     %oxygen range for bin
    Age_range(h,:,1:size(binned_fits,2))=binned_fits(:,:,20)-binned_fits(:,:,19);   %Age range for bin
    Age_source(h,:,1:size(binned_fits,2))=binned_fits(:,:,21);                               %Age source,  '1' is pure CFC-based, '2' is pure 14C-based
    p_value(h,:,1:size(binned_fits,2))=binned_fits(:,:,22);                                      %p value for slope of O2 versus age linear fit 

    sigma_range(h)=sigma_increment;%sigma increment used for a given Monte Carlo simulation

 end    %FIRST FOR-LOOP Monte Carlo for-loop

%O2_rate=-O2_rate;
%O2_rate(O2_rate<=0)=NaN;
 
 %% Plot discrete dissolution rates

for i=1:10
       a=O2_rate;
       a(Nbofpoints<3 | p_value>0.05)=NaN;
       b=O2_rate_depth;
       b(Nbofpoints<3 | p_value>0.05)=NaN;
       c=p_value;
       c(Nbofpoints<3 | p_value>0.05)=NaN;
       figure
%    scatter(squeeze(nanmean(O2_rate(:,i,:))),squeeze(nanmean(O2_rate_depth(:,i,:))),[],squeeze(nanmean(O2_rate_std(:,i,:)))./squeeze(nanmean(O2_rate(:,i,:))))
       scatter(squeeze(nanmean(a(:,i,:))),squeeze(nanmean(b(:,i,:))),5,squeeze(nanmean(c(:,i,:))));
%    hold on
%    scatter(squeeze(nanmean(O2_rate(:,i,:)))-squeeze(nanmean(O2_rate_std(:,i,:))),squeeze(nanmean(O2_rate_depth(:,i,:))),[],squeeze(nanmean(p_value(:,i,:))))
%    scatter(squeeze(nanmean(O2_rate(:,i,:)))+squeeze(nanmean(O2_rate_std(:,i,:))),squeeze(nanmean(O2_rate_depth(:,i,:))),[],squeeze(nanmean(p_value(:,i,:))))
    colorbar
end


