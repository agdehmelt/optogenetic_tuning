function [t,dy, g] = ode_euler_noise_model2(y)
Sx    = y(4); % pass parameter in y vector 
F    = y(5); % pass parameter in y vector
titration = y(8);

%options

adjust_time_step=true;

%constants

%fit with Metropolis-Hastings Monte Carlo Chain 050219
Km0=2.422982;
Km1=2.422982;
Km2=0.074542;
Km3=1;
Km4=1;
Km5=0.013965;
Km6=0.785576;
Km7=1;

k0=3.884716*Km0;
k1=0;
k2=2.035871*Km2;
k3=1.194935;
k4=3.982494;
k5=0.4174425*Km5;
k6=0.005086*Km6;
k7=0.0;

Rt=0.443;
Mt=1.24;

sigma_G = 0.000;
sigma_M = Sx;
Gt=F;

%simulation variables

T = y(7);

N = T*20;
h = T/N;                   % time step

var=sqrt(h);

t = zeros(N,1);            % prepare a place to store times
g = zeros(N,1);   

% initial time
t(1) = 0; 
g(1) = Gt; 

% initial values
Ga=y(1);
Ma=y(2);
Ra=y(3);

dy = zeros(3,1); 

dy(1,1)=Ga;
dy(1,2)=Ma;
dy(1,3)=Ra;


for i=1:N                % take N steps

    %time  
    t(i+1) = t(i) + h;
    g(i) = Gt; 
    
    Gt = Gt + titration * h;
    
        %variables
        Gi=Gt-Ga;
        Mi=Mt-Ma;
        Ri=Rt-Ra;
        
        %equations

       
        dy(i+1,1) = Ga+(k3*Gi*Ra-k4*Ga*Ma)*h;

        dy(i+1,2) = Ma+(((k5*Ra*Mi)/(Km5+Mi))+((k7*Mi)/(Km7+Mi))-((k6*Ma)/(Km6+Ma)))*h;
  
        dy(i+1,3) = Ra+(((k0*Ri*Ga)/(Km1+Ri))-((k2*Ra)/(Km2+Ra)))*h;
  
                
        %randn_G=randn;
        randn_M=randn;
        
        %apply additive noise
        dy(i+1,2) =  dy(i+1,2) + sigma_M*randn_M*var;
        
        %detect zero crossing in Ga
        if(dy(i+1,1)<0)
            dy(i+1,1)=0;
             %fprintf('0 crossing Ga  \n');
        end

        %detect zero crossing in Ma
        if(dy(i+1,2)<0)        
            %fprintf('Warning: Ma=%4.2f  \n',dy(i+1,2));
            if(adjust_time_step == true)
                factor=2;
                while dy(i+1,2)<0
                    %adjust time step
                    h_a=h/factor;
                    %recalculate with adapted timestep
                    
                    Gt = Gt - titration * h;
                    Gt = Gt + titration * h_a;
                    
                    Gi=Gt-Ga;
                    Mi=Mt-Ma;
                    Ri=Rt-Ra;

                    dy(i+1,1) = Ga+(k3*Gi*Ra-k4*Ga*Ma)*h_a;

                    dy(i+1,2) = Ma+(((k5*Ra*Mi)/(Km5+Mi))+((k7*Mi)/(Km7+Mi))-((k6*Ma)/(Km6+Ma)))*h_a;

                    dy(i+1,3) = Ra+((k0*Ri*Ga)/(Km1+Ri))-((k2*Ra)/(Km2+Ra))*h_a;
                    
                    %pick new random numbers
                    
                    randn_G=randn;
                    randn_M=randn;
        
                    %apply additive noise with adjusted timestep

                    dy(i+1,1) =  dy(i+1,1) + sigma_G*randn_G*sqrt(h_a);
                    dy(i+1,2) =  dy(i+1,2) + sigma_M*randn_M*sqrt(h_a);
                    
                    factor=factor*2;
                end
                t(i+1) = t(i) + h_a;
                %fprintf('Adjusted timestep: h=%f  \n',h_a);
                %fprintf('Value with adjusted timestep: Ma=%12.8f  \n',dy(i+1,2));
            else
                dy(i+1,2)=0;
            end
            
       
        end
        
        Ga=dy(i+1,1);
        Ma=dy(i+1,2);
        Ra=dy(i+1,3);
end
end
