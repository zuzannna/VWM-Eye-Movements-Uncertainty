sub=4;
n=3; %run number
data_=[];


for i=n:10;
file=sprintf('%d_%d.csv',sub,n);
[subject,run1,trial,accuracy,delay,disc_size,reward_at_stake,target_value,error1,memorandum_x,memorandum_y,actual_location_x,actual_location_y]=importfile_uncertainty(file);
data=[subject target_value delay disc_size error1 memorandum_x memorandum_y actual_location_x actual_location_y accuracy reward_at_stake];
data_=[data_;data];
n=n+1;
end

% prompt = 'Do you want to save the data? 0/1 ';
% z = input(prompt);
%if z==1 

save(sprintf('data_%d.mat',sub),'data_')


%%
data=group_data;
%some cleaning of outliers
too_far=find(data(:,5)>=10); %when the eye landed (final saccade)
%more than 11 dva from the target (wrong quadrant / not following the instructions
data(too_far,:)=NaN;
too_big_disc=find(data(:,4)>9);
data(too_big_disc,:)=NaN;
too_small_disc=find(data(:,4)<0);
data(too_small_disc,:)=NaN;

subjects=[4 5 6 7 9];

prompt = 'Do you want to normalize the data? 0/1 ';
z = input(prompt);
if z==1; data=uncertainty_priority_normalize(data); end

plot(data(:,5),data(:,4),'ko');
xlabel('error');
ylabel('disk size');
xlim([-6 6])
ylim([-6 6])



for h=1:length(subjects)
sub_data=choose_subject(data,subjects(h),1);

high=choose_priority(sub_data,0.6,2);
med=choose_priority(sub_data,0.3,2);
low=choose_priority(sub_data,0.1,2);

figure()
subplot(1,4,1)
plot(high(:,5),high(:,4),'ro');
xlabel('error');
ylabel('disk size');
xlim([-6 6])
ylim([-6 6])
title(sprintf('subject %d',subjects(h)));

subplot(1,4,2)
plot(med(:,5),med(:,4),'bo');
xlabel('error');
ylabel('disk size');
xlim([-6 6])
ylim([-6 6])


subplot(1,4,3)
plot(low(:,5),low(:,4),'ko');
xlabel('error');
ylabel('disk size');
xlim([-6 6])
ylim([-6 6])


subplot(1,4,4)
[f,xi] = ksdensity(high(4));

ksdensity(high(:,4));hold on;
ksdensity(med(:,4));
ksdensity(low(:,4));
legend('0.6','0.3','0.1')
title('disc size')
end

%% uncertainty vs error within bins

%choosing the bins
high=choose_priority(data,0.6,2);
med=choose_priority(data,0.3,2);
low=choose_priority(data,0.1,2);

edges=[ 0    1.5000    3.0000    4.5000    6.0000];
bin_centers=edges(2:end)-(edges(2)-edges(1))/2;

figure()
subplot(1,3,1)
hist(high(:,4),10);
xlabel('disc size')
title('high priority')
subplot(1,3,2)
hist(med(:,4),10);
xlabel('disc size')
title('medium priority')
subplot(1,3,3)
hist(low(:,4),10);
xlabel('disc size')
title('low priority')

figure()
subplot(1,3,1)
hist(high(:,5),10);
xlabel('error')
title('high priority')
subplot(1,3,2)
hist(med(:,5),10);
xlabel('error')
title('medium priority')
subplot(1,3,3)
hist(low(:,5),10);
xlabel('error')
title('low priority')
 

    bin1count=1;
    bin2count=1;
    bin3count=1;
    bin4count=1;

for i=1:length(high)

    if high(i,4)<edges(2)
        bin1(bin1count)=high(i,5);
        bin1count=bin1count+1;
        
    elseif high(i,4)>edges(2) && high(i,4)<edges(3)
        bin2(bin2count)=high(i,5);
        bin2count=bin2count+1;
        
    elseif high(i,4)>edges(3) && high(i,4)<edges(4)
        bin3(bin3count)=high(i,5);
        bin3count=bin3count+1;
        
    elseif high(i,4)>edges(4)
        bin4(bin4count)=high(i,5);
        bin4count=bin4count+1;

    end
end

figure()
discs_high=[variance_2DGaussian(bin1') variance_2DGaussian(bin2') variance_2DGaussian(bin3') variance_2DGaussian(bin4')];

plot([bin_centers(1) bin_centers(2) bin_centers(3) bin_centers(4)],discs_high,'r--o');
title('high')
xlabel('uncertainty (binned)')
ylabel('error variance')
xlim([0 6])
ylim([0 6])


    bin1count=1;
    bin2count=1;
    bin3count=1;
    bin4count=1; %add one more bin

for i=1:length(med)

    if med(i,4)<edges(2)
        binm1(bin1count)=med(i,5);
        bin1count=bin1count+1;
        
    elseif med(i,4)>edges(2) && med(i,4)<edges(3)
        binm2(bin2count)=med(i,5);
        bin2count=bin2count+1;
        
    elseif med(i,4)>edges(3) && med(i,4)<edges(4)
        binm3(bin3count)=med(i,5);
        bin3count=bin3count+1;
        
    elseif med(i,4)>edges(4)
        binm4(bin4count)=med(i,5);
        bin4count=bin4count+1;

    end
end

discs_med=[variance_2DGaussian(binm1') variance_2DGaussian(binm2') variance_2DGaussian(binm3') variance_2DGaussian(binm4')];
%subplot(1,3,2)
figure();plot([bin_centers(1) bin_centers(2) bin_centers(3) bin_centers(4)],discs_med,'b--o'); 
title('medium')
xlabel('uncertainty (binned)')
ylabel('error variance')
xlim([0 6])
ylim([0 6])


    bin1count=1;
    bin2count=1;
    bin3count=1;
    bin4count=1; %add one more bin
    

for i=1:length(low)

    if low(i,4)<edges(2)
        binl1(bin1count)=low(i,5);
        bin1count=bin1count+1;
        
    elseif low(i,4)>edges(2) && low(i,4)<edges(3)
        binl2(bin2count)=low(i,5);
        bin2count=bin2count+1;
        
    elseif low(i,4)>edges(3) && low(i,4)<edges(4)
        binl3(bin3count)=low(i,5);
        bin3count=bin3count+1;
        
    elseif low(i,4)>edges(4)
        binl4(bin4count)=low(i,5);
        bin4count=bin4count+1;

    end
end

discs_low=[variance_2DGaussian(binl1') variance_2DGaussian(binl2') variance_2DGaussian(binl3') variance_2DGaussian(binl4')];


%subplot(1,3,3)
figure();plot([bin_centers(1) bin_centers(2) bin_centers(3) bin_centers(4)],discs_low,'k--o');
title('low')
xlabel('uncertainty (binned)')
ylabel('error variance')
xlim([0 6])
ylim([0 6])


% subplot(1,2,2)
% discs_low_mean=[mean_2DGaussian(binl1') mean_2DGaussian(binl2') mean_2DGaussian(binl3') mean_2DGaussian(binl4')];
% discs_med_mean=[mean_2DGaussian(binm1') mean_2DGaussian(binm2') mean_2DGaussian(binm3') mean_2DGaussian(binm4')];
% discs_high_mean=[mean_2DGaussian(bin1') mean_2DGaussian(bin2') mean_2DGaussian(bin3') mean_2DGaussian(bin4')];

% plot([bin_centers(1) bin_centers(2) bin_centers(3) bin_centers(4)],discs_high_mean,'r--o'); hold on;
% plot([bin_centers(1) bin_centers(2) bin_centers(3) bin_centers(4)],discs_med_mean,'b--o');
% plot([bin_centers(1) bin_centers(2) bin_centers(3) bin_centers(4)],discs_low_mean,'k--o');
% legend('high', 'medium', 'low')
% xlabel('uncertainty (binned)')
% ylabel('mean error')

%%
% ===== BINNED DISK SIZE AND ERROR SD PLOT =====
nQuants = 4;
nCond=3;
% D=data(:,4);
% error_sacc=(data(:,5));



sd_errorsacc = nan(nCond,nQuants);
med_disksize = nan(nCond,nQuants);

%for icond = 1:nCond;
  
  % pair disk size and saccade error data.
  dataMat = high(:,4);
 
  [dataMat,k] = sort(dataMat,1);
  dataMat(:,2)=high(k,5);
  
  nTrials=length(dataMat);
  
  quantileEnds = linspace(0,nTrials,nQuants+1);
  figure()
  icond=1;
  for iquant = 1:nQuants; % for each bin (by quantile)...
      
      % get appropriate quantile from dataMat
      quantMat = dataMat(quantileEnds(iquant)+1:quantileEnds(iquant+1),:); 
      
      % calculate SD of saccade error for this quantile
      sd_errorsacc(icond,iquant) = nanstd(quantMat(:,2));
      
      % median of disksizes
      med_disksize(icond,iquant) = nanmedian(quantMat(:,1));
      plot(med_disksize(icond,:),sd_errorsacc(icond,:),'ro--');
  hold on; 
  end
  
    % pair disk size and saccade error data.
  dataMat = med(:,4);
 
  [dataMat,k] = sort(dataMat,1);
  dataMat(:,2)=med(k,5);
  
  nTrials=length(dataMat);
  quantileEnds = linspace(0,nTrials,nQuants+1);
  icond=2;
  for iquant = 1:nQuants; % for each bin (by quantile)...
      
      % get appropriate quantile from dataMat
      quantMat = dataMat(quantileEnds(iquant)+1:quantileEnds(iquant+1),:); 
      
      % calculate SD of saccade error for this quantile
      sd_errorsacc(icond,iquant) = nanstd(quantMat(:,2));
      
      % median of disksizes
      med_disksize(icond,iquant) = nanmedian(quantMat(:,1));
      plot(med_disksize(icond,:),sd_errorsacc(icond,:),'bo--');

  end
  
    % pair disk size and saccade error data.
  dataMat = low(:,4);
 
  [dataMat,k] = sort(dataMat,1);
  dataMat(:,2)=low(k,5);
 
  nTrials=length(dataMat);
  quantileEnds = linspace(0,nTrials,nQuants+1);
  
  icond=3;
  for iquant = 1:nQuants; % for each bin (by quantile)...
      
      % get appropriate quantile from dataMat
      quantMat = dataMat(quantileEnds(iquant)+1:quantileEnds(iquant+1),:); 
      
      % calculate SD of saccade error for this quantile
      sd_errorsacc(icond,iquant) = nanstd(quantMat(:,2));
      
      % median of disksizes
      med_disksize(icond,iquant) = nanmedian(quantMat(:,1));
      plot(med_disksize(icond,:),sd_errorsacc(icond,:),'ko--');

  end  
  
%end

%defaultplot;
xlabel('disk size')
ylabel('SD of saccade errors')
legend('0.6','0.3','0.1')









