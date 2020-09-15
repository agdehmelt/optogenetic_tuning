clear;

%options

show_traces = false;
sim_time=8000;

%select what kind of data should be generated
%please note that only one repeat of stochastic
%simulations will be generated and that the
%parameter ranges are smaller than in the
%publication to facilitate rapid data generation
list = {'Figure2J','Figure3B','Figure4A','Figure4B','Figure4C','Figure4D','Figure4FG'};
[indx,tf] = listdlg('PromptString','Generate simulations related to:','SelectionMode','single','ListString',list);

settings=char(list(indx))

%settings='Figure2J';
%settings='Figure3B';
%settings='Figure4A';
%settings='Figure4B';
%settings='Figure4C';
%settings='Figure4D';
%settings='Figure4FG';

%default settings

%range and interval for varying Gt
Gt_range=(0.1:0.1:1);


%range and interval for varying sigma_M
sigma_M_range=[0.000];

titration = 0;

%settings to generate data related to Figures

switch settings
    case {'Figure2J'}
        show_plot='AMP_vs_Gt';
        Gt_range=(0.1:0.1:1);
        sigma_M_range=[0.000];
    case {'Figure3B'}
        show_plot='Activity_vs_Time_Titration';
        Gt_range=[0.08];
        sigma_M_range=[0.000];
        titration = 1/5000;
        sim_time=3000;
    case {'Figure4A'}
        show_plot='Activity_vs_Time';
        Gt_range=[0.1];
        sigma_M_range=[0.001];
    case {'Figure4B'}
        show_plot='Activity_vs_Time';
        Gt_range=[0.3];
        sigma_M_range=[0.001];
    case {'Figure4C'}
        show_plot='Activity_vs_Time';
        Gt_range=[1.3];
        sigma_M_range=[0.001];
    case {'Figure4D'}
        show_plot='CV_vs_Gt';
        Gt_range=(0.1:0.1:1);
        sigma_M_range=[0.001];
    case {'Figure4FG'}
        show_plot='Noise_vs_Maxamp';
        Gt_range=[0.1 0.3 1.3];
        sigma_M_range=(0.0000:0.0005:0.0025);
    otherwise
        show_plot='AMP_vs_Gt';
        Gt_range=(0.1:0.1:1);
        sigma_M_range=[0.000];
end

%initialize arrays to save simulation results
data_amp = zeros([1 1]);
data_diff1 = zeros([1 1]);
data_diff2 = zeros([1 1]);
data_diff_p1 = zeros([1 1]);
data_diff_p2 = zeros([1 1]);
data_freq = zeros([1 1]);
data_osc_amp = zeros([1 1]);
data_X0 = zeros([1 1]);
data_Y0 = zeros([1 1]);
data_conv = zeros([1 1]);
data_CV = zeros([1 1]);
data_interval = zeros([1 1]);
data_Gt = zeros([1 1]);
data_sigma_M = zeros([1 1]);

%initialize variables
x_count=0;
y_count=0;


for sigma_M=sigma_M_range
    x_count=x_count+1;    
    for Gt=Gt_range
        y_count=y_count+1;
        
        %first find the steady state or limit cycle
        
        %initial conditions for simulation
        X0=Gt*0.9;
        Y0=0.0425;
        Z0=0.0425;
        
        %perform inital simuilation
        [T,Y,G] = ode_euler_noise_model2([X0 Y0 Z0 sigma_M Gt 0 6000 0]);
 
        %set new initial conditions based on initial simulation
        X0=Y(end,1);
        Y0=Y(end,2);
        Z0=Y(end,3);
        X0s=X0;
        Y0s=Y0;
        Z0s=Z0;

        %converge to steady state
        
        for conv=1:1:10
            [T,Y,G] = ode_euler_noise_model2([X0 Y0 Z0 sigma_M Gt 0 600*conv 0]);
                
            X0=Y(end,1);
            Y0=Y(end,2);
            Z0=Y(end,3);
            diff1 = abs(Y(end,1)-Y(end-1,1));
            diff2 = abs(Y(end,2)-Y(end-1,2));
            if(diff1<0.00000000000001 && diff2<0.00000000000001) break;
            end
        end
        diff_p1=NaN;
        diff_p2=NaN;
        clear pks_Y;
        %if not converged, converge to limit cycle with constant amplitude
        if(diff1>=0.000001 || diff2>=0.000001) 
            diff1=NaN;
            diff2=NaN;
            X0=X0s;
            Y0=Y0s;
            for conv=1:1:10
                [T,Y,G] = ode_euler_noise_model2([X0 Y0 Z0 sigma_M Gt 0 6000*conv 0]);
                
                [pks_X,locs_X] = findpeaks(Y(:,1));
                if size(pks_X,1)>3                
                    diff_p1=abs(pks_X(end-1)-pks_X(end-2));
                end
            
                [pks_Y,locs_Y] = findpeaks(Y(:,2));
                if size(pks_Y,1)>3               
                    diff_p2=abs(pks_Y(end-1)-pks_Y(end-2));
                end

                if(diff_p1<0.001 && diff_p2<0.001)                    
                    break;
                end
            end
            if(diff_p1>=0.001 || diff_p2>=0.001) conv=0;
            end
        end
        Z=Y;
        ZT=T;
        
        %perform the simulation

        [T,Y,G] = ode_euler_noise_model2([X0 Y0 Z0 sigma_M Gt 0 sim_time titration]);

        Y=Y(1:end-1,:);

        T=T(1:end-1,:);

        if (show_traces == true)        
            figID = figure;
            plot(T,Y(:,1),'-',T,Y(:,2),'-.')
        end
        
        %analyze simulation

        average_g=mean(Y(:,1));
        average_m=mean(Y(:,2));

        data_average_g(x_count,y_count)=average_g;
        data_average_m(x_count,y_count)=average_m;

        %set variables
        Mt=1.24;

        Ga=average_g;
        Gi=Gt-Ga;

        Ma=average_m;
        Mi=Mt-Ma;
        
        %detect peaks

        min_of_trace=min(Y(:,3));
        max_of_trace=max(Y(:,3));

        stdev_of_trace=std(Y(:,3));


        prominence_threshold=0.1*(max_of_trace-min_of_trace);

        [pks,peaks, w, p] = findpeaks(Y(:,3),'MinPeakProminence',prominence_threshold);

        [pks_nth,peaks_nth, w_nth, p_nth] = findpeaks(Y(:,3));

        %analyze peak amplitude and peak width

        peak_counter=0;
        amp_counter=0;
        width_counter=0;
        amp=zeros(1);
        width=zeros(1);
        peak_distance=zeros(1);
        if (size(peaks,1) >3);
            for i = 1:(size(peaks,1)-1)
                peak_counter=peak_counter+1;
                peak_distance(peak_counter) = T(peaks(i+1))-T(peaks(i));

                amp_counter=amp_counter+1;
                amp(peak_counter)=p(i);

                width_counter=width_counter+1;
                width(peak_counter)=w(i);
            end 
        end

        average_amp=mean(amp);
        average_width=mean(width);
        average_dist=mean(peak_distance);

        data_average_amp(x_count,y_count)=average_amp;
        data_max_amp(x_count,y_count)=max(amp);
        data_average_width(x_count,y_count)=average_width;
        data_average_dist(x_count,y_count)=average_dist;

        %compute coefficient of variation
        
        CV = nanstd(peak_distance)/nanmean(peak_distance);
        data_CV(x_count,y_count)=CV;

        interval = nanmean(peak_distance);
        data_interval(x_count,y_count)=interval;
        %figure;
        %histogram(peak_distance)

        % log data

        ap_amplitude=mean(max(Y(:,3)))-mean(min(Y(:,3)));
        data_amp(x_count,y_count)=ap_amplitude;
        data_diff1(x_count,y_count)=diff1;
        data_diff2(x_count,y_count)=diff2;
        data_diff_p1(x_count,y_count)=diff_p1;
        data_diff_p2(x_count,y_count)=diff_p2;
        data_X0(x_count,y_count)=X0;
        data_Y0(x_count,y_count)=Y0;
        data_conv(x_count,y_count)=conv;
        data_Gt(x_count,y_count)=Gt;
        data_sigma_M(x_count,y_count)=sigma_M;
        fprintf('sigma_M: %4.2f Gt:%4.2f conv:%i amp:%4.2f CV:%4.4f \n',sigma_M,Gt,conv,ap_amplitude,CV);

    end
    y_count=0;
end
figure;
if strcmp(show_plot,'AMP_vs_Gt')  
    plot(data_Gt,data_average_amp)
    xlabel('Gt');
    ylabel('Average Rho Amplitude');
    title(['data related to ' settings]);
end
if strcmp(show_plot,'CV_vs_Gt')  
    plot(data_Gt,data_CV)
    xlabel('Gt');
    ylabel('Coefficient of Variation');
    title(['data related to ' settings]);
end
if strcmp(show_plot,'Activity_vs_Time')  
    plot(T,Y(:,1))
    xlabel('Time');
    ylabel('Ga');
    title(['data related to ' settings]);
    figure;
    plot(T,Y(:,2))
    xlabel('Time');
    ylabel('Ma');
    title(['data related to ' settings]);
    figure;
    plot(T,Y(:,3))
    xlabel('Time');
    ylabel('Ra');
    title(['data related to ' settings]);
end
if strcmp(show_plot,'Activity_vs_Time_Titration')  
    plot(T,G(:,1))
    xlabel('Time');
    ylabel('Gt');
    title(['data related to ' settings]);
    figure;
    plot(T,Y(:,3))
    xlabel('Time');
    ylabel('Ra');
    title(['data related to ' settings]);
end
if strcmp(show_plot,'Noise_vs_Maxamp')  
    plot(data_sigma_M,data_max_amp)
    xlabel('Sigma_M');
    ylabel('Maximal Rho amplitude');
    legend('Gt=0.1','Gt=0.3','Gt=1.3')
    title(['data related to Figure4F']);
    figure;
    plot(data_sigma_M,data_CV)
    xlabel('Sigma_M');
    ylabel('Coefficient of Variation');
    legend('Gt=0.1','Gt=0.3','Gt=1.3')
    title(['data related to Figure4G']);
end
