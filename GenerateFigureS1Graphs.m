function GenerateFigureS1Graphs

% Function to examine signalling directionality in planar polarity.

% Biochemical reactions:
% 1. Fmi + Fmi -> Fmi:Fmi
% 2. Fz + Fmi:Fmi -> Fz-Fmi:Fmi
% 3. Vang + Fmi:Fmi -> Vang-Fmi:Fmi
% 4. Fz-Fmi:Fmi + Vang -> Fz-Fmi:Fmi-Vang
% 5. Vang-Fmi:Fmi + Fz -> Fz-Fmi:Fmi-Vang
% 6. Vang-Fmi:Fmi + Vang -> Vang-Fmi:Fmi-Vang
% 7. Fz-Fmi:Fmi + Fz -> Fz-Fmi:Fmi-Fz

% All complexes containing four proteins are stored in left hand edge of
% cell. Complex formation is simulated in a 1D row of cells with periodic
% boundary conditions.

% Tidy up
clear all
close all

%% Conditions for figure panels

% Select a clone condition:
% Uncomment appropriate line to generate loss-of-function clones
clone = ' fmi'; % Used to generate Fig. S1D
% clone = ' none'; % Used to generate Fig. S1E
% clone = ' fz'; % Used to generate Figs. S1F

% Select appropriate Vmax,F and Vmax,V values
% both are equal to 7 in FigS1D
% both are equal to 38 in FigS1E upper and 32 in FigS1E lower
% both are equal to 20 in FigS1F upper and 10 in FigS1F lower
Vmaxf = 7; % Vmax,F for destabilising feedback
Vmaxv = 7; % Vmax,V for destabilising feedback

% Select flag for stabilising or destabilising feedback
% use destabilising only for Fig.S1E,F and Fig.S1D lower
% use stabilising for Fig.S1D upper
destabilising = 1;
stabilising = 0;

% Select flag for whether to apply an initial bias in Fz localisation
% set to 1 for Fig.S1D and F
% set to 0 for Fig.S1E
apply_bias = 1;

%% Parameters

num_cells = 30; % number of num_cells

% Mass action rate constants for direct monodirectional signalling
k1 = 1; k2 = 1; k3 = 1; k4 = 1; k5 = 1; k6 = 1; k7 = 1;
v1 = 5; v2 = 0.1; v3 = 0.1; v4 = 0.01; v5 = 0.1; v6 = 0.1; v7 = 0.1;

% Feedback parameters for Hill functions
thresh = 0.5; %referred to as K in paper
expon = 2; %referred to as w in paper

% Diffusion rates
coeff = 0.03;
delta_x = 5;
diff_coeff = coeff/(delta_x^2);
D1 = diff_coeff; %Fz
D2 = diff_coeff; %Vang
D3 = diff_coeff; %Fmi

%% Simulation settings

tstart = 0;  % start time
tend = 100000; % end time
nrep = 200;  % number of solution reports

%% Initial conditions

% Initial protein concentrations
start = 1; % initial amount on each side of each cell

init_Fmi_L = start*ones(num_cells, 1);
init_Fmi_R = start*ones(num_cells, 1);
init_Fz_L = start*ones(num_cells, 1);
init_Fz_R = start*ones(num_cells, 1);
init_Vang_L = start*ones(num_cells, 1);
init_Vang_R = start*ones(num_cells, 1);
init_Fmi_Fmi = zeros(num_cells, 1);
init_Fmi_Fmi_Fz_L = zeros(num_cells, 1);
init_Fz_Fmi_Fmi_R = zeros(num_cells, 1);
init_Fmi_Fmi_Vang_L = zeros(num_cells, 1);
init_Vang_Fmi_Fmi_R = zeros(num_cells, 1);
init_Fz_Fmi_Fmi_Vang = zeros(num_cells, 1);
init_Vang_Fmi_Fmi_Fz = zeros(num_cells, 1);
init_Vang_Fmi_Fmi_Vang = zeros(num_cells, 1);
init_Fz_Fmi_Fmi_Fz = zeros(num_cells, 1);

totF = init_Fz_L(1,1)+ init_Fz_R(1,1);

% Apply small initial distal bias of Fz
if apply_bias
    bias = init_Fz_L(1,1) * 0.001; %referred to as b in paper
    init_Fz_R(:) = (0.5 * totF) + bias;
    init_Fz_L(:) = (0.5 * totF) - bias;
end

% Implement clones

if strcmp(clone, ' fz')
    init_Fz_L(6:10) = 0;
    init_Fz_R(6:10) = 0;
end

if strcmp(clone, ' Vang')
    init_Vang_L(6:10) = 0;
    init_Vang_R(6:10) = 0;
end

if strcmp(clone,' fmi')
    init_Fmi_L(6:10) = 0;
    init_Fmi_R(6:10) = 0;
end

if strcmp(clone,' Vangfz')
    init_Fz_L(6:10) = 0;
    init_Fz_R(6:10) = 0;
    init_Vang_L(6:10) = 0;
    init_Vang_R(6:10) = 0;
end

%% Run simulation

init = [init_Fmi_L; init_Fmi_R; init_Fz_L; init_Fz_R; init_Vang_L; ...
    init_Vang_R; init_Fmi_Fmi; init_Fmi_Fmi_Fz_L; init_Fz_Fmi_Fmi_R; ...
    init_Fmi_Fmi_Vang_L; init_Vang_Fmi_Fmi_R; init_Fz_Fmi_Fmi_Vang; ...
    init_Vang_Fmi_Fmi_Fz; init_Vang_Fmi_Fmi_Vang; init_Fz_Fmi_Fmi_Fz];

options = odeset('MaxStep', 1e2, 'RelTol', 1e-8, 'AbsTol', 1e-8);
sol = ode45(@odes, [tstart, tend], init, options, num_cells);
tout = linspace(tstart, tend, nrep);
stout = deval(sol, tout);

%% Find solutions and plot

time = tout;
sol_Fmi_L = stout(1:num_cells,:);
sol_Fmi_R = stout(num_cells+1:2*num_cells,:);
sol_Fz_L = stout(2*num_cells+1:3*num_cells,:);
sol_Fz_R = stout(3*num_cells+1:4*num_cells,:);
sol_Vang_L = stout(4*num_cells+1:5*num_cells,:);
sol_Vang_R = stout(5*num_cells+1:6*num_cells,:);
sol_Fmi_Fmi = stout(6*num_cells+1:7*num_cells,:);
sol_Fmi_Fmi_Fz_L = stout(7*num_cells+1:8*num_cells,:);
sol_Fz_Fmi_Fmi_R = stout(8*num_cells+1:9*num_cells,:);
sol_Fmi_Fmi_Vang_L = stout(9*num_cells+1:10*num_cells,:);
sol_Vang_Fmi_Fmi_R = stout(10*num_cells+1:11*num_cells,:);
sol_Fz_Fmi_Fmi_Vang = stout(11*num_cells+1:12*num_cells,:);
sol_Vang_Fmi_Fmi_Fz = stout(12*num_cells+1:13*num_cells,:);
sol_Vang_Fmi_Fmi_Vang = stout(13*num_cells+1:14*num_cells,:);
sol_Fz_Fmi_Fmi_Fz = stout(14*num_cells+1:15*num_cells,:);

% Calculate bound quantities
sol_bound_Fz_L = sol_Fmi_Fmi_Fz_L + sol_Fz_Fmi_Fmi_Fz + sol_Vang_Fmi_Fmi_Fz;
sol_bound_Fz_R = sol_Fz_Fmi_Fmi_R + circshift(sol_Fz_Fmi_Fmi_Fz,[-1,0]) + circshift(sol_Fz_Fmi_Fmi_Vang,[-1,0]) ;
sol_bound_Vang_L = sol_Fmi_Fmi_Vang_L + sol_Fz_Fmi_Fmi_Vang + sol_Vang_Fmi_Fmi_Vang;
sol_bound_Vang_R = sol_Vang_Fmi_Fmi_R + circshift(sol_Vang_Fmi_Fmi_Fz,[-1,0]) + circshift(sol_Vang_Fmi_Fmi_Vang,[-1,0]);
sol_bound_Fz = sol_bound_Fz_L + sol_bound_Fz_R;
sol_bound_Vang = sol_bound_Vang_L + sol_bound_Vang_R;
sol_bound_Fmi_L = sol_Fmi_Fmi + sol_Fmi_Fmi_Fz_L + circshift(sol_Fz_Fmi_Fmi_R,[1,0]) + sol_Fmi_Fmi_Vang_L + circshift(sol_Vang_Fmi_Fmi_R,[1,0]) + sol_Vang_Fmi_Fmi_Fz + sol_Fz_Fmi_Fmi_Vang + sol_Fz_Fmi_Fmi_Fz + sol_Vang_Fmi_Fmi_Vang;
sol_bound_Fmi_R =  circshift(sol_bound_Fmi_L,[-1,0]);

% Plot solutions over time to check for steady state
figure,
plot(tout,sol_bound_Fz_L,tout,sol_bound_Fz_R, 'LineWidth', 2)
set(gca,'fontsize',10);
xlabel('simulated time')
ylabel('Bound Fz')

% Plot bound amounts of Fz/Vang in each cell
figure,
cell_plotter_LR(num_cells, sol_bound_Fz_L(:,end), sol_bound_Fz_R(:,end), time(:,end), 'Fz')
set(gca,'fontsize',10);
xlabel('cell number')
ylabel('Bound Fz')

figure,
cell_plotter_LR(num_cells, sol_bound_Vang_L(:,end), sol_bound_Vang_R(:,end), time(:,end), 'Vang')
set(gca,'fontsize',10);
xlabel('cell number')
ylabel('Bound Vang')

%% Subfunctions
    function f = hill(x, V, t)
        % Increasing Hill function with min of 1 and max of V
        f = 1+((V-1).*(x.^expon))./(t^expon + x.^expon);
    end

    function dy_dt = odes(~,y,num_cells)
        
        % Split initial vector y into different species
        Fmi_L = y(1:num_cells,1);
        Fmi_R = y(num_cells+1:2*num_cells,1);
        Fz_L = y(2*num_cells+1:3*num_cells,1);
        Fz_R = y(3*num_cells+1:4*num_cells,1);
        Vang_L = y(4*num_cells+1:5*num_cells,1);
        Vang_R = y(5*num_cells+1:6*num_cells,1);
        Fmi_Fmi = y(6*num_cells+1:7*num_cells,1);
        Fmi_Fmi_Fz_L = y(7*num_cells+1:8*num_cells,1);
        Fz_Fmi_Fmi_R = y(8*num_cells+1:9*num_cells,1);
        Fmi_Fmi_Vang_L = y(9*num_cells+1:10*num_cells,1);
        Vang_Fmi_Fmi_R = y(10*num_cells+1:11*num_cells,1);
        Fz_Fmi_Fmi_Vang = y(11*num_cells+1:12*num_cells,1);
        Vang_Fmi_Fmi_Fz = y(12*num_cells+1:13*num_cells,1);
        Vang_Fmi_Fmi_Vang = y(13*num_cells+1:14*num_cells,1);
        Fz_Fmi_Fmi_Fz = y(14*num_cells+1:15*num_cells,1);
        
        % Store bound quantities
        bound_Fz_L = Fmi_Fmi_Fz_L + Fz_Fmi_Fmi_Fz + Vang_Fmi_Fmi_Fz;
        bound_Fz_R = Fz_Fmi_Fmi_R + circshift(Fz_Fmi_Fmi_Fz,[-1,0]) + circshift(Fz_Fmi_Fmi_Vang,[-1,0]);  %% both should shift -1
        bound_Vang_L = Fmi_Fmi_Vang_L + Fz_Fmi_Fmi_Vang + Vang_Fmi_Fmi_Vang;
        bound_Vang_R = Vang_Fmi_Fmi_R + circshift(Vang_Fmi_Fmi_Fz,[-1,0]) + circshift(Vang_Fmi_Fmi_Vang,[-1,0]);
        
        % Compute feedback multipliers
        V1 = v1;
        V2_L = v2;
        V2_R = v2;
        V3_L = v3;
        V3_R = v3;
        V4_L = v4;
        V4_R = v4;
        V5_L = v5;
        V5_R = v5;
        V6_L = v6;
        V6_R = v6;
        V7_L = v7;
        V7_R = v7;
        
        if destabilising
            V1 = v1;
            V2_L = V2_L*hill(bound_Vang_L,Vmaxv, thresh);
            V2_R = V2_R*hill(bound_Vang_R,Vmaxv, thresh);
            V3_L = V3_L*hill(bound_Fz_L,Vmaxf, thresh);
            V3_R = V3_R*hill(bound_Fz_R,Vmaxf, thresh);
            V4_L = V4_L*hill(bound_Fz_L,Vmaxf, thresh);
            V4_R = V4_R*hill(circshift(bound_Fz_R,[1,0]),Vmaxf, thresh);
            V5_L = V5_L*hill(bound_Vang_L,Vmaxv, thresh);
            V5_R = V5_R*hill(circshift(bound_Vang_R,[1,0]),Vmaxv, thresh);
            V6_L = V6_L*hill(bound_Fz_L,Vmaxf, thresh);
            V6_R = V6_R*hill(circshift(bound_Fz_R,[1,0]),Vmaxf, thresh);
            V7_L = V7_L*hill(bound_Vang_L,Vmaxv, thresh);
            V7_R = V7_R*hill(circshift(bound_Vang_R,[1,0]),Vmaxv, thresh);
        end;
        
        if stabilising
            K1 = k1;
            K2_L = k2*hill(bound_Fz_L,Vmaxf, thresh);
            K2_R = k2*hill(bound_Fz_R,Vmaxf, thresh);
            K3_L = k3*hill(bound_Vang_L,Vmaxv, thresh);
            K3_R = k3*hill(bound_Vang_R,Vmaxv, thresh);
            K4_L = k4*hill(bound_Vang_L,Vmaxv, thresh);
            K4_R = k4*hill(circshift(bound_Vang_R,[1,0]),Vmaxv, thresh);
            K5_L = k5*hill(bound_Fz_L,Vmaxf, thresh);
            K5_R = k5*hill(circshift(bound_Fz_R,[1,0]),Vmaxf, thresh);
            K6_L = k6*hill(bound_Vang_L,Vmaxspos, thresh);
            K6_R = k6*hill(circshift(bound_Vang_R,[1,0]),Vmaxv, thresh);
            K7_L = k7*hill(bound_Fz_L,Vmaxf, thresh);
            K7_R = k7*hill(circshift(bound_Fz_R,[1,0]),Vmaxf, thresh);
        else
            K1 = k1*ones(num_cells,1);
            K2_L = k2*ones(num_cells,1);
            K2_R = k2*ones(num_cells,1);
            K3_L = k3*ones(num_cells,1);
            K3_R = k3*ones(num_cells,1);
            K4_L = k4*ones(num_cells,1);
            K4_R = k4*ones(num_cells,1);
            K5_L = k5*ones(num_cells,1);
            K5_R = k5*ones(num_cells,1);
            K6_L = k6*ones(num_cells,1);
            K6_R = k6*ones(num_cells,1);
            K7_L = k7*ones(num_cells,1);
            K7_R = k7*ones(num_cells,1);
        end
        
        % Calculate reaction rates
        rate_1 = K1.*Fmi_L.*circshift(Fmi_R,[1,0]) - V1*Fmi_Fmi;
        rate_2L = K2_L.*Fz_L.*Fmi_Fmi - V2_L.*Fmi_Fmi_Fz_L;
        rate_2R = K2_R.*Fz_R.*circshift(Fmi_Fmi,[-1,0]) - V2_R.*Fz_Fmi_Fmi_R;
        rate_3L = K3_L.*Vang_L.*Fmi_Fmi - V3_L.*Fmi_Fmi_Vang_L;
        rate_3R = K3_R.*Vang_R.*circshift(Fmi_Fmi,[-1,0]) - V3_R.*Vang_Fmi_Fmi_R;
        rate_4L = K4_L.*Vang_L.*circshift(Fz_Fmi_Fmi_R,[1,0]) - V4_L.*Fz_Fmi_Fmi_Vang;
        rate_4R = K4_R.*circshift(Vang_R,[1,0]).*Fmi_Fmi_Fz_L - V4_R.*Vang_Fmi_Fmi_Fz;
        rate_5L = K5_L.*Fz_L.*circshift(Vang_Fmi_Fmi_R,[1,0]) - V5_L.*Vang_Fmi_Fmi_Fz;
        rate_5R = K5_R.*circshift(Fz_R,[1,0]).*Fmi_Fmi_Vang_L - V5_R.*Fz_Fmi_Fmi_Vang;
        rate_6L = K6_L.*Vang_L.*circshift(Vang_Fmi_Fmi_R,[1,0]) - V6_L.*Vang_Fmi_Fmi_Vang;
        rate_6R = K6_R.*circshift(Vang_R,[1,0]).*Fmi_Fmi_Vang_L - V6_R.*Vang_Fmi_Fmi_Vang;
        rate_7L = K7_L.*Fz_L.*circshift(Fz_Fmi_Fmi_R,[1,0]) - V7_L.*Fz_Fmi_Fmi_Fz;
        rate_7R = K7_R.*circshift(Fz_R,[1,0]).*Fmi_Fmi_Fz_L - V7_R.*Fz_Fmi_Fmi_Fz;
        rate_8 = D1*(Fz_L - Fz_R);
        rate_9 = D2*(Vang_L - Vang_R);
        rate_10 = D3*(Fmi_L - Fmi_R);
        
        d_Fmi_L_dt = -rate_1 -rate_10;
        d_Fmi_R_dt = -circshift(rate_1,[-1,0]) +rate_10;
        d_Fz_L_dt = -rate_2L -rate_5L -rate_7L -rate_8;
        d_Fz_R_dt = -rate_2R -circshift(rate_5R,[-1,0]) -circshift(rate_7R,[-1,0]) +rate_8;
        d_Vang_L_dt = -rate_3L -rate_4L -rate_6L -rate_9;
        d_Vang_R_dt = -rate_3R  -circshift(rate_4R,[-1,0]) -circshift(rate_6R,[-1,0]) +rate_9;
        d_Fmi_Fmi_dt = rate_1 -rate_2L -circshift(rate_2R,[1,0]) -rate_3L -circshift(rate_3R,[1,0]);
        d_Fmi_Fmi_Fz_L_dt = rate_2L -rate_4R -rate_7R;
        d_Fz_Fmi_Fmi_R_dt = rate_2R -circshift(rate_4L,[-1,0]) -circshift(rate_7L,[-1,0]);
        d_Fmi_Fmi_Vang_L_dt = rate_3L -rate_5R - rate_6R;
        d_Vang_Fmi_Fmi_R_dt = rate_3R -circshift(rate_5L,[-1,0]) -circshift(rate_6L,[-1,0]);
        d_Fz_Fmi_Fmi_Vang_dt = rate_4L +rate_5R;
        d_Vang_Fmi_Fmi_Fz_dt = rate_4R +rate_5L;
        d_Vang_Fmi_Fmi_Vang_dt = rate_6L +rate_6R;
        d_Fz_Fmi_Fmi_Fz_dt = rate_7L +rate_7R;
        
        dy_dt = [d_Fmi_L_dt; d_Fmi_R_dt; d_Fz_L_dt; d_Fz_R_dt; ...
            d_Vang_L_dt; d_Vang_R_dt; d_Fmi_Fmi_dt; ...
            d_Fmi_Fmi_Fz_L_dt; d_Fz_Fmi_Fmi_R_dt; ...
            d_Fmi_Fmi_Vang_L_dt; d_Vang_Fmi_Fmi_R_dt; ...
            d_Fz_Fmi_Fmi_Vang_dt; d_Vang_Fmi_Fmi_Fz_dt; ...
            d_Vang_Fmi_Fmi_Vang_dt; d_Fz_Fmi_Fmi_Fz_dt];
    end

    function cell_plotter_LR(num_cells, array_L, array_R, t, name)
        % Plots a 1D row of cells, bar heights demonstrate amount of
        % material in each cell edge
        
        % Use RGB for Fz or Vang
        if strcmp(name,'Fz')
            col = [0.0118, 0.7725, 0.0588];
        elseif strcmp(name,'Vang')
            col = [0.9843, 0.5059, 0.0745];
        end
        
        % Initalise the coords, BL=bottom left, TR=top right etc
        BLy = zeros(1, num_cells);
        BRy = zeros(1, num_cells);
        BLx = zeros(1, num_cells);
        BRx = zeros(1, num_cells);
        
        TLy = zeros(1, num_cells);
        TRy = zeros(1, num_cells);
        TLx = zeros(1, num_cells);
        TRx = zeros(1, num_cells);
        
        % Set x coords based on cell number
        for i=1:num_cells
            BLx(1,i)=i-0.45;
            TLx(1,i)=i-0.45;
            
            BRx(1,i)=i+0.45;
            TRx(1,i)=i+0.45;
        end
        
        % Set y coords for top of bar based on simulation
        TLy = array_L';
        TRy = array_R';
        
        % Concatanate coords for use in patch function
        X = [BLx; TLx; TRx; BRx];
        Y = [BLy; TLy; TRy; BRy];
        
        p = patch(X,Y,[0,0,0]);
        clear cdata
        set(p,'FaceColor', col)
        axis([0 num_cells+1 0 1]);
    end

end
