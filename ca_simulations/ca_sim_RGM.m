clearvars -except settings

delete(findall(0,'Type','figure'))


%general parameters for CA simulation
Param.Stacks=1;
Param.StochVal=0.5;
Param.OutFreq=1;
frame_start=1;
CA_size=100;
noise=0.01; %initial homogenous noise between individual cells of the CA

%run options
save_1=false;
save_2=false;
save_3=false;
save_4=false;
save_5=true;
save_6=false;

save_tiffs=true;

sim_num=1;

% %point perturbation
length_factor=8;

% %quick sim
% length_factor=1;
% CA_size=50;

Param.Reps=2250*length_factor;
save_start=1;
save_end=Param.Reps;

%simulation step to start the spatio-temporal correlation analysis
start_step=1000*length_factor;

%fit with Metropolis-Hastings Monte Carlo Chain 190515

%Michaelis constants
Km0_GEF_act_Rho=2.422982;
Km1_act_Rho=2.422982;
Km2_inact_Rho=0.074542;
Km3=1;
Km4=1;
Km5_Rho_act_Myo=0.013965;
Km6_inact_Myo=0.785576;
Km7_act_Myo=1;
Km8_Myo_inact_Rho=Km2_inact_Rho;
Hill_Rho_act_Myo=1.0;

%rate constants
k0_GEF_act_Rho=3.884716*Km0_GEF_act_Rho;
k1_act_Rho=0;
k2_inact_Rho=2.035871*Km2_inact_Rho;
k3_Rho_act_GEF=1.194935;
k4_Myo_inact_GEF=3.982494;
k5_Rho_act_Myo=0.4174425*Km5_Rho_act_Myo;
k6_inact_Myo=0.005086*Km6_inact_Myo;
k7_act_Myo=0.0;
k2prime_Myo_inact_Rho=0;

%multiply all rate constants with factor for timestep 0.25s
all_rates=0.25;

%rate constants with factor
k0_GEF_act_Rho=k0_GEF_act_Rho*all_rates;
k1_act_Rho=k1_act_Rho*all_rates;
k2_inact_Rho=k2_inact_Rho*all_rates;
k3_Rho_act_GEF=k3_Rho_act_GEF*all_rates;
k4_Myo_inact_GEF=k4_Myo_inact_GEF*all_rates;
k5_Rho_act_Myo=k5_Rho_act_Myo*all_rates;
k6_inact_Myo=k6_inact_Myo*all_rates;
k7_act_Myo=k7_act_Myo*all_rates;
k2prime_Myo_inact_Rho=k2prime_Myo_inact_Rho*all_rates;

%concentrations in million molecules per cell
myo_conc=1.24;
rho_conc=0.443;
gef_titration=0.0000;
gef_conc=0.06;
%gef_conc=2.0;

%diffusion factors for the CA based on fitting
%simulations to known diffusion coefficients

rho_diff=1.95;
myo_diff=1.2;
gef_diff=rho_diff;

rho_cyto_diff=8.78;
myo_cyto_diff=3.01;
gef_cyto_diff=rho_cyto_diff;

%sigmas for noise input
krand1=0.0;
krand2=0.001;
krand3=0.0;

switch settings
    case {'scan_gefpm'}
        %gefpm
        save_tiffs=false;
        x_range = [0.615384615 0.641025641	0.687179487	0.738461538	1	1.282051282	1.794871795	2.307692308];
        y_range = [1 2 3 4];
        vary_parameter = 1; 
        do_spatial_correlation_analysis=true;
    case {'scan_gefcyto'}
        %gefcyto
        save_tiffs=false;
        x_range = [0.375854214	0.512528474	0.683371298	1	1.412300683	1.958997722	2.801822323];
        y_range = [1 2 3 4];
        vary_parameter = 2; 
        do_spatial_correlation_analysis=true;
    case {'scan_myopm'}
        %myopm
        save_tiffs=false;
        x_range = [0.841666667	0.916666667	0.966666667	1	1.083333333	1.175833333	1.5];
        y_range = [1 2 3 4];
        vary_parameter = 3; 
        do_spatial_correlation_analysis=true;
    case {'scan_myocyto'}
        %myocyto
        save_tiffs=false;
        x_range = [0.465116279	0.574750831	0.754152824	1	1.328903654	1.827242525	2.491694352];
        y_range = [1 2 3 4];
        vary_parameter = 4; 
        do_spatial_correlation_analysis=true;
    case {'scan_rhopm'}
        %rhopm
        save_tiffs=false;
        x_range = [0.615384615 0.641025641	0.687179487	0.738461538	1	1.282051282	1.794871795	2.307692308];
        y_range = [1 2 3 4];
        vary_parameter = 5; 
        do_spatial_correlation_analysis=true;
    case {'scan_rhocyto'}
        %rhocyto
        save_tiffs=false;
        x_range = [0.375854214	0.512528474	0.683371298	1	1.412300683	1.958997722	2.801822323];
        y_range = [1 2 3 4];
        vary_parameter = 6; 
        do_spatial_correlation_analysis=true;
    case {'single_low_gt'}
        %single
        save_tiffs=true;
        vary_parameter = 0; 
        x_range = [1];
        y_range = [1];
        gef_conc=0.06;
        do_spatial_correlation_analysis=false;
    case {'single_high_gt'}
        %single
        save_tiffs=true;
        vary_parameter = 0; 
        x_range = [1];
        y_range = [1];
        gef_conc=2.0;
        do_spatial_correlation_analysis=false;
    otherwise
        error('Settings not found! Please start with launch_CA.m');
end

x_count=0;
y_count=0;
for y=y_range
    y_count=y_count+1;
    x_count=0;
    for x=x_range

        x_count=x_count+1;

        %assign reaction rate parameters into struct for simulation
        Param.IntPar = [k0_GEF_act_Rho k1_act_Rho k2_inact_Rho k3_Rho_act_GEF k4_Myo_inact_GEF k5_Rho_act_Myo k6_inact_Myo k7_act_Myo k2prime_Myo_inact_Rho 0.0 ...
        Km0_GEF_act_Rho Km1_act_Rho Km2_inact_Rho 1.0 1.0 Km5_Rho_act_Myo Km6_inact_Myo Km7_act_Myo Km8_Myo_inact_Rho ...
        Hill_Rho_act_Myo rho_conc gef_titration 0 krand1 krand2 krand3 ...
        ];

        %assign diffusion parameters into struct for simulation

        %                     1          2            3            4            5            6           
        DiffRadName = {'DiffGEFPM' 'DiffGEFCyto' 'DiffMyoPM' 'DiffMyoCyto' 'DiffRhoPM' 'DiffRhoCyto' };

        Param.DiffRad = [gef_diff gef_cyto_diff myo_diff myo_cyto_diff rho_diff rho_cyto_diff];

        %readable names for parameters
        %                     1             2              3                4                 5                    6               7              8                   9                10          11               12              13               14                   15              16                 17              18               19                        20                 21               22                   23       24       25       26      
        IntParName = {'k0-GEF_act_Rho', 'k1-act_Rho', 'k2-inact_Rho', 'k3-Rho_act_GEF', 'k4-Myo_inact_GEF', 'k5-Rho_act_Myo', 'k6-inact_Myo', 'k7-act_Myo', 'k2prime-Myo_inact_Rho', 'empty', 'Km0-GEF_act_Rho', 'Km1-act_Rho', 'Km2-inact_Rho', 'Km3-Rho_act_GEF', 'Km4-inact_GEF', 'Km5-Rho_act_Myo', 'Km6-inact_Myo', 'Km7-act_Myo', 'Km2prime-Myo_inact_Rho', 'Hill-Rho_act_Myo' 'Rho_total_conc', 'GEF_titration_rate', 'not_used''krand1','krand2','krand3'};

        %assign general simulation parameters for the dynamic
        %simulation

        %6PDE calculation
        Param.Array=[CA_size CA_size gef_conc myo_conc 10 noise rho_conc];

        %multiply with factor from x_range to scan effect of
        %varying one parameter
        if(vary_parameter ~= 0)
            Param.DiffRad(vary_parameter)=Param.DiffRad(vary_parameter)*x;
        end

        %perform CA simulation
        calc2D_RGM;
        output_folder='';
        if save_tiffs==true;
            saveTiffs;
        end
        %extract Rho activity data from result matrix
        Rho_a=single(squeeze(O(:,:,1,5,:)));
        if do_spatial_correlation_analysis==true;
            [half_size] = spatial_correlation_analysis(Rho_a, start_step,49);
            data_spatial_width_Rho(x_count,y_count)=half_size;
        end

    end
end

fprintf('Parameters:\n');

for pval=1:1:25
    fprintf('%s=%4.6f\n', char(IntParName(pval)), Param.IntPar(pval));
end

for pval=1:1:6
    fprintf('%s=%4.6f\n', char(DiffRadName(pval)), Param.DiffRad(pval));
end
if do_spatial_correlation_analysis==true
    dlmwrite(['data_spatial_width_' settings '.csv'],data_spatial_width_Rho);
    fprintf(['saved csv file with results: data_spatial_width_' settings '.csv']);
end