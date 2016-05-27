function data_group_standardized=priority_normalize(data)
subjects=unique(data(:,1));
data_group_standardized=NaN(size(data));
%normalizing
for h=1:length(subjects)
    
    sub_data_ind=find(data(:,1)==subjects(h));
    sub_data=choose_subject(data,subjects(h),1);
    
    miu=nanmean(sub_data);
    sigma=nanstd(sub_data);
    
    for j=3:size(sub_data,2) %don't standardize subject number or target value
        if j==4 %don't standardize delay periods
            j=j+1;
        else
            standardized_column=(sub_data(:,j)-miu(j))/sigma(j);
            data_group_standardized(sub_data_ind,j)=standardized_column;
        end
        
        
    end
end

data_group_standardized(:,[1 2 4 9 10 11 12 13 14])=data(:,[1 2 4 9 10 11 12 13 14]);
