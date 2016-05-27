function priority_data=choose_priority(data, priority, column)

priority_data_ind=find(data(:,column)==priority);
priority_data=data(priority_data_ind,:);
%putvar sub_data
end