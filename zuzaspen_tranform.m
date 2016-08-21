
[group_data(:,14), group_data(:,15)]=zuzaspenthingtoday(group_data(:,13));
errr=sacc_error(group_data(:,6), group_data(:,7),group_data(:,14),group_data(:,15));
sum((errr-group_data(:,5))>0.0001)

%% look at rows that corrective and final saccades in different quadrants (and nan rows)
idx = find(sign(group_data(:,6)) ~= sign(group_data(:,14)));
idx = unique([idx; find(sign(group_data(:,7)) ~= sign(group_data(:,15)))]);

data_exclude = group_data(idx,:);
hist(data_exclude(:,2)) % this plot shows that there are higher lapse rates for 0.1 priority than for the other two.
% something potentially of interest.

%% throwing out nans and trials with diff quadrants (see above)
group_data(idx,:) = [];

%%
[theta_f,r_f]=cart2pol(abs(group_data(:,6)),abs(group_data(:,7)));
[theta_c,r_c]=cart2pol(abs(group_data(:,14)),abs(group_data(:,15)));

theta_f_deg=theta_f*(180/pi);
theta_c_deg=theta_c*(180/pi);

error_degrees=theta_f_deg-theta_c_deg;

group_data(:,16)=error_degrees;
group_data(:,17)=r_f;


%% BIASES (FOR EACH SUBJECT INDIVIDUALLY)

prompt = 'Do you want to normalize the data? 0/1 ';
z = input(prompt);
if z==1; data=uncertainty_priority_normalize(data); end

prompt = 'Do you want to save the figures? 0/1 ';
c = input(prompt);

prompt = 'Do you want to exclude outliers? 0/1 ';
d = input(prompt);

subjects=unique(group_data(:,1));
data=group_data;

%too_big_disc=find(data(:,4)>10); % we asked participants to draw the disc size of
% the screen if something was screwed up with final saccade location
%data(too_big_disc,:)=NaN;
if d==1
too_small_disc=find(data(:,4)<0); % if the size of the disc was smaller than 0
%too_close=find(data(:,17)<=6);
too_far=find(data(:,17)>=14); %when the eye landed (final saccade)

data(too_far,:)=NaN;
data(too_close,:)=NaN;
data(too_small_disc,:)=NaN;

end
%too_big_disc=find(data(:,4)>9);
%data(too_big_disc,:)=NaN;

data(:,18)=zeros;

% slice the data into 4 within the quadrant & plot radial and circular
% error for each slice & subject separately
for i=1:length(data)
    
    if 1<theta_f_deg(i) & theta_f_deg(i)<22.5
        data(i,18)=1;
    elseif 22.6<theta_f_deg(i) & theta_f_deg(i)<45
        data(i,18)=2;
    elseif 46<theta_f_deg(i) & theta_f_deg(i)<67.5
        data(i,18)=3;
    elseif 67.6<theta_f_deg(i) & theta_f_deg(i)<90
        data(i,18)=4;
    end
    
end

for r=1:length(subjects)
    
    % ============= radius ===============
    %plot average radius difference (r_f-10) as a function of priority
    % we want to see no change in bias with priority
    
    sub_data=choose_subject(data,subjects(r),1);
    high=choose_priority(sub_data,0.6,2);
    med=choose_priority(sub_data,0.3,2);
    low=choose_priority(sub_data,0.1,2);
    
    % average radial error as a function of priority
    
    r_f_mean(r,1)=nanmean(low(:,17))-10;
    r_f_mean(r,2)=nanmean(med(:,17))-10;
    r_f_mean(r,3)=nanmean(high(:,17))-10;
    
    r_f_mean_abs(r,1)=abs(nanmean(low(:,17))-10);
    r_f_mean_abs(r,2)=abs(nanmean(med(:,17))-10);
    r_f_mean_abs(r,3)=abs(nanmean(high(:,17))-10);
    SEMs_radial_mean(r,1:3)=[(std(low(:,17)))/sqrt(length(low(:,17))) (std(med(:,17)))/sqrt(length(med(:,17))) (std(high(:,17)))/sqrt(length(high(:,17)))];
    
    
    
    if r==length(subjects)
        
        figure()
        errorbar(mean(r_f_mean_abs), std(r_f_mean_abs)/sqrt(length(r_f_mean_abs)),'rx') %errrorbar function doesn't accept
        ylabel('mean radial error (r-10), absolute value')
        xlabel('priority')
        title('magnitude of mean radial error as a function of priority (all subjects)')
    end
    
    % variance of radial error as a function of priority
    
    r_f_var(r,1)=nanstd(low(:,17));
    r_f_var(r,2)=nanstd(med(:,17));
    r_f_var(r,3)=nanstd(high(:,17));
    
%     figure()
%         plot([0.5 3.5],[0 0],'k-')
%         hold on;
%         errorbar([r_f_mean(r,1) r_f_mean(r,2) r_f_mean(r,3)], [r_f_var(r,1)/sqrt(length(low(:,17))) r_f_var(r,1)/sqrt(length(med(:,17))) r_f_var(r,1)/sqrt(length(high(:,17)))],'rx') %errrorbar function doesn't accept
%         ylabel('mean radial error (r-10), absolute value')
%         xlabel('priority')
%         titlee=sprintf('subject %d: mean radial error as a function priority',subjects(r));
%         title(titlee)
%         if c==1
% 
%         file=sprintf('radial_error_%d',subjects(r));
%         print(file,'-dpng')
%         end
    
    
    if r==length(subjects)
        figure()
        bar(r_f_var)
        xlabel('subject and priority')
        ylabel('std of radial error')
        legend('0.1','0.3','0.6')
        title('std of radial error as a function of priority (all subjects)')
        if c==1
        
        file=sprintf('radial_error_std_allsubs');
        print(file,'-dpng')
        end
       
        figure()
        bar(r_f_mean)
        xlabel('subject and priority')
        ylabel('mean radial error')
        legend('0.1','0.3','0.6')
        title('mean radial error as a function of priority (all subjects)')
        if c==1
        
        file=sprintf('radial_error_mean_allsubs');
        print(file,'-dpng')
        end
 
        
        figure()
        plot(mean(r_f_mean),'.','Markersize',20)
        hold on;
        plot([0.5 3.5],[0 0],'k-')
        errorbar(mean(r_f_mean), std(r_f_mean)/sqrt(length(r_f_mean)),'rx')
        xlabel('priority')
        ylabel('mean radial error')
        title('mean radial error as a function of priority (all subjects combined)')
        if c==1
        
        file=sprintf('radial_error_mean_allsubs_combined');
        print(file,'-dpng')
        end
        
        figure()
        plot(mean(r_f_var),'.','Markersize',20)
        hold on;
        errorbar(mean(r_f_var), std(r_f_var)/sqrt(length(r_f_var)),'rx')
        xlabel('priority')
        ylabel('std of radial error')
        title('std of radial error as a function of priority (all subjects combined)')
        if c==1
        
        file=sprintf('radial_error_mean_allsubs_combined');
        print(file,'-dpng')
        end
    end
    
    %==========circle size==========
    circle_mean(r,1)=nanmean(low(:,4));
    circle_mean(r,2)=nanmean(med(:,4));
    circle_mean(r,3)=nanmean(high(:,4));
    
    circle_var(r,1)=nanstd(low(:,4));
    circle_var(r,2)=nanstd(med(:,4));
    circle_var(r,3)=nanstd(high(:,4));
    
      if r==length(subjects)
        figure()
        plot(mean(circle_mean),'.','Markersize',20)
        hold on;
        errorbar(mean(circle_mean), std(circle_mean)/sqrt(length(circle_mean)),'rx')
        xlabel('priority')
        ylabel('mean circle size')
        title('mean circle size as a function of priority (all subjects combined)')
        if c==1
        
        file=sprintf('circle_size_allsubs_combined');
        print(file,'-dpng')
        end
      end   
    
    % ========== theta ============
    %separate the data based on target location (quadrant into four slices)
    
    theta_mean(r,1)=nanmean(abs((low(:,16))));
    theta_mean(r,2)=nanmean(abs((med(:,16))));
    theta_mean(r,3)=nanmean(abs((high(:,16))));
    
    theta_var(r,1)=nanstd(abs((low(:,16))));
    theta_var(r,2)=nanstd(abs((med(:,16))));
    theta_var(r,3)=nanstd(abs((high(:,16))));
    
    slices_m(1,1,r)=nanmean(abs(low(find(low(:,18)==1),16)));
    slices_v(1,1,r)=nanstd(low(find(low(:,18)==1),16));
    
    slices_m(1,2,r)=nanmean(abs(med(find(med(:,18)==1),16)));
    slices_v(1,2,r)=nanstd(med(find(med(:,18)==1),16));
    
    slices_m(1,3,r)=nanmean(abs(high(find(high(:,18)==1),16)));
    slices_v(1,3,r)=nanstd(high(find(high(:,18)==1),16));
    
    slices_m(2,1,r)=nanmean(abs(low(find(low(:,18)==2),16)));
    slices_v(2,1,r)=nanstd(low(find(low(:,18)==2),16));
    
    slices_m(2,2,r)=nanmean(abs(med(find(med(:,18)==2),16)));
    slices_v(2,2,r)=nanstd(med(find(med(:,18)==2),16));
    
    slices_m(2,3,r)=nanmean(abs(high(find(high(:,18)==2),16)));
    slices_v(2,3,r)=nanstd(high(find(high(:,18)==2),16));
    
    slices_m(3,1,r)=nanmean(abs(low(find(low(:,18)==3),16)));
    slices_v(3,1,r)=nanstd(low(find(low(:,18)==3),16));
    
    slices_m(3,2,r)=nanmean(abs(med(find(med(:,18)==3),16)));
    slices_v(3,2,r)=nanstd(med(find(med(:,18)==3),16));
    
    slices_m(3,3,r)=nanmean(abs(high(find(high(:,18)==3),16)));
    slices_v(3,3,r)=nanstd(high(find(high(:,18)==3),16));
    
    slices_m(4,1,r)=nanmean(abs(low(find(low(:,18)==4),16)));
    slices_v(4,1,r)=nanstd(low(find(low(:,18)==4),16));
    
    slices_m(4,2,r)=nanmean(abs(med(find(med(:,18)==4),16)));
    slices_v(4,2,r)=nanstd(med(find(med(:,18)==4),16));
    
    slices_m(4,3,r)=nanmean(abs(high(find(high(:,18)==4),16)));
    slices_v(4,3,r)=nanstd(high(find(high(:,18)==4),16));
    
    
    
    slices_m_nabs(1,1,r)=nanmean(low(find(low(:,18)==1),16));
    
    slices_m_nabs(1,2,r)=nanmean(med(find(med(:,18)==1),16));

    slices_m_nabs(1,3,r)=nanmean(high(find(high(:,18)==1),16));
    
    slices_m_nabs(2,1,r)=nanmean(low(find(low(:,18)==2),16));
    
    slices_m_nabs(2,2,r)=nanmean(med(find(med(:,18)==2),16));
    
    slices_m_nabs(2,3,r)=nanmean(high(find(high(:,18)==2),16));
    
    slices_m_nabs(3,1,r)=nanmean(low(find(low(:,18)==3),16));
     
    slices_m_nabs(3,2,r)=nanmean(med(find(med(:,18)==3),16));
    
    slices_m_nabs(3,3,r)=nanmean(high(find(high(:,18)==3),16));
    
    slices_m_nabs(4,1,r)=nanmean(low(find(low(:,18)==4),16));
    
    slices_m_nabs(4,2,r)=nanmean(med(find(med(:,18)==4),16));
    
    slices_m_nabs(4,3,r)=nanmean(high(find(high(:,18)==4),16));
     
    % plot average theta difference for each slice as a function of priority
    
%     figure()
%     bar(slices_m_nabs(:,:,r))
%     xlabel('slice and priority')
%     ylabel('mean of theta')
%     legend('0.6','0.3','0.1')
%     titlee=sprintf('subject %d: mean theta error as a function of within quadrant location and priority',subjects(r));
%     title(titlee)
    
    
    % variance of error as a functin of priority
    
%     figure()
%     bar(slices_v(:,:,r))
%     xlabel('slice and priority')
%     ylabel('std of theta')
%     legend('0.6','0.3','0.1')
%     titlee=sprintf('subject %d: variance of theta error as a function of within quadrant location and priority',subjects(r));
%     title(titlee)
    
    if r==length(subjects)
        figure()
        bar(mean(slices_m_nabs,3))
        hold on;
        errorbar(mean(slices_m_nabs,3), std(slices_m_nabs,0,3)/sqrt(length(slices_m_nabs)),'rx') %errrorbar function doesn't accept
        % multiple rows so it'll look weird, gotta fix it in illustrator
        % and move the bars horizontally
        xlabel('slice and priority')
        legend('0.1','0.3','0.6')
        ylabel('mean of theta error')
        title('mean theta error as a function of within quadrant location and priority (all subjects)')
        
        figure()
        plot(mean(theta_mean),'.','Markersize',20)
        hold on;
        %plot([0.5 3.5],[0 0],'k-')
        errorbar(mean(theta_mean), std(theta_mean)/sqrt(length(theta_mean)),'rx')
        xlabel('priority')
        ylabel('saccade error (theta, in degrees)')
        title('saccade error as a function of priority (all subjects combined)')
        
        figure()
        plot(mean(theta_var),'.','Markersize',40)
        %hold on;
        %plot([0.5 3.5],[0 0],'k-')
        %errorbar(mean(theta_var),'rx')
        xlabel('priority')
        ylabel('std of saccade error (theta, in degrees)')
        title('std of saccade error as a function of priority (all subjects combined)')
        
    end
    
    
end

% correlation between error and confidence (quartile binning thing)
% ========== maybe add for the whole group too ========


for r=1:length(subjects)
    
    sub_data=choose_subject(data,subjects(r),1);
    high=choose_priority(sub_data,0.6,2);
    med=choose_priority(sub_data,0.3,2);
    low=choose_priority(sub_data,0.1,2);
    
    figure()
    subplot(1,3,1)
    hist(high(:,4),10);
    xlabel('disc size')
    title('high priority')
    ylim([0 100])
    subplot(1,3,2)
    hist(med(:,4),10);
    xlabel('disc size')
    title('medium priority')
    ylim([0 100])
    subplot(1,3,3)
    hist(low(:,4),10);
    xlabel('disc size')
    title('low priority')
    ylim([0 100])
    if c==1;
        file=sprintf('discsize_%d',subjects(h));
        print(file,'-dpng')
    end
    
    
    figure()
    subplot(1,3,1)
    hist(high(:,16),10);
    ylim([0 100])
    xlabel('error')
    title('high priority')
    subplot(1,3,2)
    hist(med(:,16),10);
    xlabel('error')
    title('medium priority')
    ylim([0 100])
    subplot(1,3,3)
    hist(low(:,16),10);
    xlabel('error')
    title('low priority')
    ylim([0 100])
    if c==1;
        file=sprintf('error_%d',subjects(r));
        print(file,'-dpng')
    end
    
    
    nQuants = 4;
    nCond=3;
    % D=data(:,4);
    % error_sacc=(data(:,16));
    
    sd_errorsacc = nan(nCond,nQuants);
    med_disksize = nan(nCond,nQuants);
    
    %for icond = 1:nCond;
    
    % pair disk size and saccade error data.
    dataMat = high(:,4);
    
    [dataMat,k] = sort(dataMat,1);
    dataMat(:,2)=abs(high(k,16));
    
    nTrials=length(dataMat);
    
    quantileEnds = linspace(0,nTrials,nQuants+1);
    
    icond=1;
    for iquant = 1:nQuants; % for each bin (by quantile)...
        
        % get appropriate quantile from dataMat
        quantMat = dataMat(quantileEnds(iquant)+1:quantileEnds(iquant+1),:);
        
        % calculate SD of saccade error for this quantile
        mean_errorsacc1(icond,iquant) = nanmean(quantMat(:,2));
        
        % median of disksizes
        mean_disksize1(icond,iquant) = nanmean(quantMat(:,1));
        %std_quant_disc1(icond,iquant)=nanstd(quantMat(:,1));
        %std_quant_errorsacc1(icond,iquant)=nanstd(quantMat(:,2));
        SEMs_disksize1(icond,iquant)=(nanstd(quantMat(:,1)))/(sqrt(length(quantMat(:,1))));
        SEMs_errorsacc1(icond,iquant)=(nanstd(quantMat(:,2)))/(sqrt(length(quantMat(:,2))));
        %plot(mean_disksize1(icond,:),mean_errorsacc1(icond,:),'ro--');
        
        %hold on;
    end
    
    figure()
    errorbar(mean_disksize1,mean_errorsacc1, SEMs_disksize1,'ro--')
    hold on;
    
    % pair disk size and saccade error data.
    dataMat = med(:,4);
    
    [dataMat,k] = sort(dataMat,1);
    dataMat(:,2)=abs(med(k,16));
    
    nTrials=length(dataMat);
    quantileEnds = linspace(0,nTrials,nQuants+1);
    icond=2;
    
    for iquant = 1:nQuants; % for each bin (by quantile)...
        
        % get appropriate quantile from dataMat
        quantMat = dataMat(quantileEnds(iquant)+1:quantileEnds(iquant+1),:);
        
        % calculate SD of saccade error for this quantile
        mean_errorsacc2(1,iquant) = nanmean(quantMat(:,2));
        
        % median of disksizes
        mean_disksize2(1,iquant) = nanmean(quantMat(:,1));
        SEMs_disksize2(1,iquant)=(nanstd(quantMat(:,1)))/(sqrt(length(quantMat(:,1))));
        SEMs_errorsacc2(1,iquant)=(nanstd(quantMat(:,2)))/(sqrt(length(quantMat(:,2))));
        
        
        %plot(mean_disksize2(icond,:),mean_errorsacc2(icond,:),'bo--');
        
    end
    %figure()
    errorbar(mean_disksize2,mean_errorsacc2, SEMs_disksize2,'bo--')
    
    % pair disk size and saccade error data.
    dataMat = low(:,4);
    
    [dataMat,k] = sort(dataMat,1);
    dataMat(:,2)=abs(low(k,16));
    
    nTrials=length(dataMat);
    quantileEnds = linspace(0,nTrials,nQuants+1);
    
    icond=3;
    for iquant = 1:nQuants; % for each bin (by quantile)...
        
        % get appropriate quantile from dataMat
        quantMat = dataMat(quantileEnds(iquant)+1:quantileEnds(iquant+1),:);
        
        % calculate SD of saccade error for this quantile
        mean_errorsacc3(1,iquant) = nanmean(quantMat(:,2));
        
        % median of disksizes
        mean_disksize3(1,iquant) = nanmean(quantMat(:,1));
        SEMs_disksize3(1,iquant)=(nanstd(quantMat(:,1)))/(sqrt(length(quantMat(:,1))));
        SEMs_errorsacc3(1,iquant)=(nanstd(quantMat(:,2)))/(sqrt(length(quantMat(:,2))));
        
        %plot(mean_disksize3(icond,:),mean_errorsacc3(icond,:),'ko--');
        
    end
    
    %end
    
    errorbar(mean_disksize3,mean_errorsacc3, SEMs_disksize3,'ko--')
    
    
    %defaultplot;
    xlabel('Mean circle size')
    ylabel('Mean saccade error (theta, in degrees)')
    legend('0.6','0.3','0.1')
    
    %if c==1;
        file=sprintf('error_discsize_quantiles_%d',subjects(r));
        print(file,'-dpng')
    %end
    
       
end
