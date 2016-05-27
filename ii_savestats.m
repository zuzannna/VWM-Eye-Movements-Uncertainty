function ii_savestats(subject, run)
sub=num2str(subject);
run1=num2str(run);

cd ~
cd ..
cd ..

cd data/cueprob/EyeTracking_Probability/v2/processed/
save('ii_stats','ii_stats');
workspace=sprintf('%d_%d',subject,run);
save(sprintf('%d_%d',subject,run))
