% load in data file
load('group_data.mat')

%% calculate target locations

% calculate x and y locations
target_x = 10*cosd(group_data(:,13));
target_y = 10*sind(group_data(:,13));

% rename useful columns of group_data.mat
response_x = group_data(:,6);
response_y = group_data(:,7);

% delete lapses and nans

idx = find(sign(response_x) ~= sign(target_x));
idx = unique([idx; find(sign(response_y) ~= sign(target_y))]);
response_x(idx) = [];
response_y(idx) = [];
target_x(idx) = [];
target_y(idx) = [];

% put them all over to one side
response_x = abs(response_x);
response_y = abs(response_y);
target_x = abs(target_x);
target_y = abs(target_y);

% plot error around every target
figure;
plot(target_x,target_y,'k.')

%% calculating which repsonse goes with which target
data_priority = group_data(:,2);
data_priority(idx) = [];
priorityVec = unique(data_priority); % priority conditio
titleVec = {'low','med','high'};
data_subj = group_data(:,1);
data_subj(idx) = [];
subjVec = unique(data_subj);

for isubj = 1:length(subjVec)
    subjnum = subjVec(isubj);
    
    for ipriority = 1:length(priorityVec);
        priority = priorityVec(ipriority);
        
        figure;
        targets_x = unique(target_x);
        for itarget = 1:length(targets_x);
            currtarget_y = target_y(find(target_x == targets_x(itarget),1,'first')); % value of corresponding y_cordingate
            idxx = target_x == targets_x(itarget);
            idxxx = idxx & (data_priority == priority) & (data_subj == subjnum);
            
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
        title(['subj' num2str(subjnum) titleVec{ipriority} 'priority'])
        axis equal
        savefig(gcf,['subj' num2str(subjnum) '_' titleVec{ipriority} ])
%         savefig(gcf,['subj' num2str(subjnum) '_' titleVec{ipriority} '.png' ])
    end
end
