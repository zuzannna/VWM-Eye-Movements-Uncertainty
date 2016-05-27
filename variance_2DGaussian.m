function fitted_variance=variance_2DGaussian(data)
%this function returns fitted variance of a 2D Gaussian 

fit=fitgmdist(data,1);
fitted_variance=fit.Sigma;
end