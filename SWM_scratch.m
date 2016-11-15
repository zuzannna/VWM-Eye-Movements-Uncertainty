%% understanding Jbar, tau, and sigma

sigma = 5;
Jbar = 16;
tau = .0001;

nSamps = 1e5;
samps = sqrt(1./gamrnd(Jbar/tau,tau,1,nSamps));

hist(1./samps.^2)
hist(samps)

figure; hist(samps)
mean(samps)

%% ==========================================================
%             STUFF WITH ACTUAL DATA!
% ===========================================================
% 08.25.2016

% load in data file
load('group_data.mat')

% rename useful columns of group_data.mat
data_subj = group_data(:,1);
data_priority = group_data(:,2);
data_discsize = group_data(:,4);
response_x = group_data(:,6);
response_y = group_data(:,7);
target_x = 10*cosd(group_data(:,13)); % x location
target_y = 10*sind(group_data(:,13)); % y location

% % delete nans (and potentially also trials where targt and response are not in same quadrant)
% idx = find(sign(response_x) ~= sign(target_x));
% idx = unique([idx; find(sign(response_y) ~= sign(target_y))]);
idx = find(isnan(response_x) | isnan(response_y) | isnan(target_x) | isnan(target_y)  );
response_x(idx) = [];
response_y(idx) = [];
target_x(idx) = [];
target_y(idx) = [];
data_priority(idx) = [];
data_subj(idx) = [];
data_discsize(idx) = [];

% % project everything onto first quadrant
% sign_x = sign(target_x);
% sign_y = sign(target_y);
% response_x = response_x.*sign_x;
% response_y = response_y.*sign_y;
% target_x = target_x.*sign_x;
% target_y = target_y.*sign_y;

% priority and subject numbers
priorityVec = unique(data_priority); % priority condition
nPriorities = length(priorityVec);
titleVec = {'low','med','high'};
subjVec = unique(data_subj);
nSubj = length(subjVec);

% changing to polar coordinates
[target_theta, target_rho] = cart2pol(target_x,target_y);
[response_theta, response_rho] = cart2pol(response_x,response_y);

% errors
error_x = response_x - target_x;
error_y = response_y - target_y;
error_euclid = sqrt(error_x.^2 + error_y.^2);
error_rho = response_rho - target_rho;
error_theta = response_theta - target_theta;



%% plot raw data (targets and responses)
figure;
plot(target_x,target_y,'ro')
hold on; 
plot(response_x,response_y,'k.')

%% plotting mean and covariance plot of each of 8 target locations
% best if you project onto one quadrant (see section above)

for isubj = 1:nSubj
    subjnum = subjVec(isubj);
    
    for ipriority = 1:nPriorities;
        priority = priorityVec(ipriority);
        
        figure;
        targets_x = unique(target_x);
        for itarget = 1:length(targets_x);
            currtarget_y = target_y(find(target_x == targets_x(itarget),1,'first')); % value of corresponding y_cordingate
            idxx = target_x == targets_x(itarget); % trials which have the correct target
            idxxx = idxx & (data_priority == priority) & (data_subj == subjnum); % correct target, priority condition, and subject number
            
            hold on;
            plot(targets_x(itarget),currtarget_y,'k.','MarkerSize',20)
            plot(response_x(idxxx),response_y(idxxx),'.')
            
            % get and plot cov circle for each target
            responses = [response_x(idxxx) response_y(idxxx)];
            covMat = cov(responses);
            circle = bsxfun(@plus,mean(responses),[cos(linspace(0,2*pi,50))' sin(linspace(0,2*pi,50))']*covMat);
            plot(circle(:,1),circle(:,2),'Color',aspencolors('coral')); % plot covariance circle
            
        end
        defaultplot;
        title(['subj' num2str(subjnum) ' ' titleVec{ipriority} ' priority'])
        axis equal
%         savefig(gcf,['subj' num2str(subjnum) '_' titleVec{ipriority} ])
%         savefig(gcf,['subj' num2str(subjnum) '_' titleVec{ipriority} '.png' ])
    end
end


%% main effects of error and disc size from priority

discsizeMat = nan(nSubj,nPriorities);
euclideanerrorMat = discsizeMat;
rhoerrorMat = discsizeMat; thetaerrorMat = discsizeMat;
for isubj = 1:nSubj;
    subjnum = subjVec(isubj);
    
    for ipriority = 1:nPriorities;
        priority = priorityVec(ipriority);
        
        idx = (data_subj == subjnum) & (data_priority == priority); % indices of trials for current subject
        
        % euclidean error
        euclideanerrorMat(isubj,ipriority) = mean(error_euclid(idx));
        euclideanvarMat(isubj,ipriority) = var(error_euclid(idx));
        
        % mean other errors
        rhoerrorMat(isubj,ipriority) = mean(error_rho(idx));
        thetaerrorMat(isubj,ipriority) = mean(error_theta(idx));
        
        % mean disc size
        discsizeMat(isubj,ipriority) = mean(data_discsize(idx));
    end
    
end

% mean euclidean error main effect
figure;
meaneucliderror = mean(euclideanerrorMat);
semeucliderror = std(euclideanerrorMat)/sqrt(nSubj);
errorbar(meaneucliderror,semeucliderror,'k','LineStyle','none')
set(gca,'XTick',[1 2 3],'XTickLabel',[0.1 0.3 0.6])
xlabel('priority'); title('mean euclidean error'); defaultplot;

% var euclidean error main effect
figure;
meanvareucliderror = mean(euclideanvarMat);
semvareucliderror = std(euclideanvarMat)/sqrt(nSubj);
errorbar(meanvareucliderror,semvareucliderror,'k','LineStyle','none')
set(gca,'XTick',[1 2 3],'XTickLabel',[0.1 0.3 0.6])
xlabel('priority'); title('euclidean error variance'); defaultplot;

% % mean disc size
% figure;
% meandiscsize = mean(discsizeMat);
% semdiscsize = std(discsizeMat)/sqrt(nSubj);
% errorbar(meandiscsize,semdiscsize,'k','LineStyle','none')
% set(gca,'XTick',[1 2 3],'XTickLabel',[0.1 0.3 0.6])
% xlabel('priority'); title('discsize'); defaultplot;

%% plotting bias as a function of priority across or within each subject

aveplot = 1;
indvlplot = 0;

biasVec_theta = nan(nSubj,nPriorities);
varVec_theta = nan(nSubj,nPriorities);
biasVec_rho = nan(nSubj,nPriorities);
varVec_rho = nan(nSubj,nPriorities);
for isubj = 1:nSubj % for each subject
    subjnum = subjVec(isubj);
    
    for ipriority = 1:nPriorities % for each priority
        priority = priorityVec(ipriority);
        idx = (data_subj == subjnum) & (data_priority == priority); % indices of trials for current subject
        
        % calculating bias and variance of rho and theta
        biasVec_rho(isubj,ipriority) = mean(error_rho(idx));
        varVec_rho(isubj,ipriority) = var(error_rho(idx));
        biasVec_theta(isubj,ipriority) = mean(error_theta(idx));
        varVec_theta(isubj,ipriority) = var(error_theta(idx));
    end
    
    if (indvlplot)
        figure;
        subplot(2,2,1); 
        plot(biasVec_rho(isubj,:),'ko')
        title('bias rho')
        defaultplot;
        
        subplot(2,2,2); 
        plot(biasVec_theta(isubj,:),'ko')
        title('bias theta')
        defaultplot;
        
        subplot(2,2,3); 
        plot(varVec_rho(isubj,:),'ko')
        title('variance rho')
        defaultplot;
        
        subplot(2,2,4); 
        plot(varVec_theta(isubj,:),'ko')
        title('variance theta')
        defaultplot;
        
        pause;
    end
end

if (aveplot)
    figure;
    subplot(2,2,1);
    mean_bias_rho = mean(biasVec_rho);
    sem_bias_rho = std(biasVec_rho)/sqrt(size(biasVec_rho,1));
    errorbar(mean_bias_rho,sem_bias_rho,'k')
    set(gca,'XTick',[1 2 3],'XTickLabel',[0.1 0.3 0.6])
    xlabel('priority'); title('bias rho'); defaultplot;
    
    subplot(2,2,2);
    mean_bias_theta = mean(biasVec_theta);
    sem_bias_theta = std(biasVec_theta)/sqrt(size(biasVec_theta,1));
    errorbar(mean_bias_theta,sem_bias_theta,'k')
    set(gca,'XTick',[1 2 3],'XTickLabel',[0.1 0.3 0.6])
    xlabel('priority'); title('bias theta'); defaultplot;
    
    subplot(2,2,3);
    mean_var_rho = mean(varVec_rho);
    sem_var_rho = std(varVec_rho)/sqrt(size(varVec_rho,1));
    errorbar(mean_var_rho,sem_var_rho,'k')
    set(gca,'XTick',[1 2 3],'XTickLabel',[0.1 0.3 0.6])
    xlabel('priority'); title('variance rho'); defaultplot;
    
    subplot(2,2,4);
    mean_var_theta = mean(varVec_theta);
    sem_var_theta = std(varVec_theta)/sqrt(size(varVec_theta,1));
    errorbar(mean_var_theta,sem_var_theta,'k')
    set(gca,'XTick',[1 2 3],'XTickLabel',[0.1 0.3 0.6])
    xlabel('priority'); title('variance theta'); defaultplot;
end

%% scatterplot of disc size vs error (euclidean or circular tangential)

plot(data_discsize,error_euclid,'.')

%% plot of variances for different priorities (with cov matrix)

% project everything onto first quadrant (may or may not have already done
% this)
sign_x = sign(target_x);
sign_y = sign(target_y);
response_x = response_x.*sign_x;
response_y = response_y.*sign_y;
target_x = target_x.*sign_x;
target_y = target_y.*sign_y;

% changing to polar coordinates
[target_theta, target_rho] = cart2pol(target_x,target_y);
[response_theta, response_rho] = cart2pol(response_x,response_y);

% collapse them all onto the same spot
error_x = response_x - target_x;
error_y = response_y - target_y;

colorVec = {'r','b','k'};
figure; 
for ipriority = 1:nPriorities
    priority = priorityVec(ipriority);
    idx = (data_priority == priority); % indices of trials for current priority
    
    % plot target and responses
    subplot(2,2,ipriority); hold on
    plot(error_x(idx),error_y(idx),'.','Color',[0.7 0.7 0.7]) % errors
    plot(0,0,'k.','MarkerSize',14) % target
    title(num2str(priority))
    defaultplot
    
    % plot covariance matrix
    covMat = cov([error_x(idx) error_y(idx)]);
    circle = bsxfun(@plus,mean([error_x(idx) error_y(idx)]),[cos(linspace(0,2*pi,50))' sin(linspace(0,2*pi,50))']*covMat);
    plot(circle(:,1),circle(:,2),colorVec{ipriority}); % plot covariance circle
    
    % plot just covariance matrices on top of one another
    subplot(2,2,4); hold on
    plot(circle(:,1),circle(:,2),colorVec{ipriority}); % plot covariance circle
    axis([-10 10 -20 20])
    defaultplot
end


%% plotting quantile bins of error with discsize

nQuantiles = 4;
meandiscsizeVec = nan(nSubj,nPriorities,nQuantiles);
meaneuclideanerrorVec = meandiscsizeVec; 
meanangularerrorVec = meandiscsizeVec;

for isubj = 1:nSubj
    subjnum = subjVec(isubj);
    
    for ipriority = 1:nPriorities
        priority = priorityVec(ipriority);
        
        idx = (data_subj == subjnum) & (data_priority == priority); % correct subject number and priority
        
        % get the relevant informtion for current subject and priority
        discsize = data_discsize(idx);
        euclideanerror = error_euclid(idx);
        angularerror = error_theta(idx);
        
        % sort by increasing discsize
        [discsize, I] = sort(discsize);
        euclideanerror = euclideanerror(I);
        angularerror = angularerror(I);
        
        % chop up into quantiles
        ndatapoints = length(discsize);
        idx_quantile = [0 round([1:nQuantiles].*(ndatapoints/nQuantiles))];

        for iquantile = 1:nQuantiles
            meandiscsizeVec(isubj,ipriority,iquantile) = mean(discsize(idx_quantile(iquantile)+1:idx_quantile(iquantile+1)));
            meaneuclideanerrorVec(isubj,ipriority,iquantile) = mean(euclideanerror(idx_quantile(iquantile)+1:idx_quantile(iquantile+1)));
            meanangularerrorVec(isubj,ipriority,iquantile) = mean(angularerror(idx_quantile(iquantile)+1:idx_quantile(iquantile+1)));
        end
        
    end
end

% means and SEMs across subjects
mean_discsize = squeeze(mean(meandiscsizeVec));
mean_euclideanerror = squeeze(mean(meaneuclideanerrorVec));
mean_angularerror = squeeze(mean(meanangularerrorVec));
sem_discsize = squeeze(std(meandiscsizeVec))/sqrt(nSubj);
sem_euclideanerror = squeeze(std(meaneuclideanerrorVec))/sqrt(nSubj);
sem_angularerror = squeeze(std(meanangularerrorVec))/sqrt(nSubj);

% plot euclidean error
colorVec = {'r','b','k'};
figure; hold on;
for ipriority = 1:nPriorities;
    errorbar(mean_euclideanerror(ipriority,:),sem_euclideanerror(ipriority,:),colorVec{ipriority});
end
defaultplot;
ylabel('euclidean error')
xlabel('disc size quantile')


% plot
figure; hold on;
for ipriority = 1:nPriorities;
    errorbar(mean_angularerror(ipriority,:),sem_angularerror(ipriority,:),colorVec{ipriority});
end
defaultplot;
ylabel('angular error')
xlabel('disc size quantile')
