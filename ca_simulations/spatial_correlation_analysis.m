function [ width ] = spatial_correlation_analysis( for_corr1, start_step,range)

for_corr1=for_corr1(:,:,start_step:size(for_corr1,3));
ICS2DCorr = corrfunc(for_corr1);
ics_mean=mean(ICS2DCorr,3);
%figure
%plot(ics_mean(49,:))
ICS2DCorrCrop = autocrop(ics_mean,24);
warning('off','all');
a  = gaussfit(ICS2DCorrCrop,'2d',1,'n');
warning('on','all');
sigma_x=mean(a(:,2))/sqrt(2);
sigma_y=mean(a(:,3))/sqrt(2);
sigma_average=(sigma_x+sigma_y)/2;
width_px=sigma_average*2*sqrt(2*log(2));
width=width_px/2;
fprintf('width (px): %4.2f \n', width_px);
fprintf('width (um): %4.2f \n', width);
fprintf('a2: %4.2f \n', mean(a(:,2)));
fprintf('a3: %4.2f \n', mean(a(:,3)));
            