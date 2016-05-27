function fitted_mean=mean_2DGaussian(data)
%this function returns fitted variance of a 2D Gaussian 

fit=fitgmdist(data,1);
fitted_mean=fit.mu;
end