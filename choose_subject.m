function sub_data=choose_subject(data, subject_number, column)

sub_data_ind=find(data(:,column)==subject_number);
sub_data=data(sub_data_ind,:);
%putvar sub_data
end