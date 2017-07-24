%% Numerical simulation code
% Code designed to estimate distributions for G1 and G2 growth as a
% function of arbitrary inputs.  Separate inputs for both cell cycle phases
% with only the outputs recorded:
%
% Pre: (birth size, Start size, growth pre-Start, duration pre-Start)
% Post: (Start size, division size, growth post-Start, duration post-Start)

close all
clear all

load('20170316 ACT1pr-mCi Fit Parameters Mother G2 2D Lambda.mat')

%% Simulation Parameters

% This will determine the size of the starting population (total number of
% cells will =cell_num*cell_bin_num arranged between the size range
% indicated for each cell cycle stage

cell_num = 1000;
cell_bin_num = 1000;
g1_bin_range = [1e4 5e4];

%% Cell Parameters
G1_lambda = coeffvalues(min_fit1_act);

bud_prob_beta_1 = G1_lambda(1);      %per minute;
bud_prob_M0_1 = -G1_lambda(2)/bud_prob_beta_1;      %fl

G2_lambda_2d = exp_2d_coeff;

constant_growth_rate = all_fit_act(1); %min_a_act; %per min
constant_daughter_growth_rate = log(daughter_growth_fit_act(1)*3 + 1)/3;
bud_mass_correction = d_b_diff;

poststart_G1_timer = PSG1_timer_fit;
SG2M_timer = SG2M_timer_fit;

%% G1 Simulation

% Simulation Input

n_cells = cell_num * cell_bin_num;
birth_size_list = ...
    linspace(g1_bin_range(1),g1_bin_range(2),cell_bin_num + 1);
birth_size_list = birth_size_list(1:(length(birth_size_list)-1)) + ...
    diff(birth_size_list);

cell_v = zeros(n_cells,1); %sample from log-normal volume dist
for i = 1:cell_bin_num
    cell_v((cell_num*(i-1) + 1):(cell_num*i),1) = ...
        birth_size_list(i)*ones(cell_num,1);
end
cell_cycle = zeros(n_cells,1);  %begin with G1 cells; G1 = 0, G2 = 1
cell_bud = zeros(n_cells,1);  %begin with no buds
G2_counter = zeros(n_cells,1);  %cells count down G2
mother_daughter = zeros(n_cells,1); %begin with all daughters; daughter = 0, mother = 1
mother_G1_counter = zeros(n_cells,1); %all cells start in G1 as daughters
post_start_G1 = zeros(n_cells,1);
post_start_counter = zeros(n_cells,1);
start_size = zeros(n_cells,1);
cell_list_G1 = 1:n_cells;

bud_prob = zeros(length(cell_v(:,1)),1);

% Simulation

size_at_birth = cell_v;
size_at_start = zeros(n_cells,1);
G1_duration = zeros(n_cells,1);
G1_growth = zeros(n_cells,1);

current_time = 0;

while n_cells - (sum(size_at_start > 0)) > 0
    
    current_time = current_time + 1;
    
    ['Current Time: ' num2str(current_time)]
    ['Cells Remaining: ' num2str(n_cells - (sum(size_at_start > 0)))]
    
    i = current_time;
    
    %Begin simulation
    bud_prob = bud_prob_beta_1.*(cell_v(:,i) - bud_prob_M0_1); %piecewise linear fits from data
    bud_prob(bud_prob < 0) = 0;
    bud_prob(bud_prob > 1) = 1;
%     bud_prob = .0183*ones(length(cell_v(:,i)),1); % For constant k
    bud_prob = 1 - exp(-bud_prob*(1));
    
    bud = zeros(length(bud_prob),1);
    for j = 1:length(bud_prob)
        bud(j) = randsample([0 1], 1, true, [1-bud_prob(j), bud_prob(j)]);
    end
    bud = transpose(bud);
    
    for k = cell_list_G1   %for all cells
        growth_size = cell_v(k,i) + cell_bud(k,i);
        dv = (exp(constant_daughter_growth_rate) - 1) * growth_size;
        if dv < 0
            dv = 0;
        end
        if cell_cycle(k,i) == 0 %if cells are in G1
            if bud(k) == 0 %if cell does not bud this time
                cell_v(k,i+1) = cell_v(k,i) + dv;
                cell_cycle(k,i+1) = 0; %stays in G1
                cell_bud(k,i+1) = 0;
                mother_daughter(k,i+1) = mother_daughter(k,i);
            elseif bud(k) == 1 %if cell passes start
                cell_v(k,i+1) = cell_v(k,i) + dv;
                cell_cycle(k,i+1) = 0; %stays in G1
                cell_bud(k,i+1) = 0;
                mother_daughter(k,i+1) = mother_daughter(k,i);
                
                size_at_start(k) = cell_v(k,i);
                G1_duration(k) = current_time;
                G1_growth(k) = cell_v(k,i) - size_at_birth(k);
                
                cell_list_G1(cell_list_G1 == k) = [];
                
            end
        end
    end
end

%% G2 Simulation

% Simulation Input

n_cells = cell_num * cell_bin_num;

cell_v = size_at_start; %sample from log-normal volume dist

cell_cycle = zeros(n_cells,1);  %begin with G1 cells; G1 = 0, G2 = 1
cell_bud = zeros(n_cells,1);  %begin with no buds
G2_counter = zeros(n_cells,1);  %cells count down G2
mother_daughter = zeros(n_cells,1); %begin with all daughters; daughter = 0, mother = 1
mother_G1_counter = polyval(poststart_G1_timer,cell_v); %all cells start in G1 as daughters
post_start_G1 = ones(n_cells,1);
start_size = cell_v(:,1);
post_start_counter = round(polyval(poststart_G1_timer,start_size));
post_start_counter(post_start_counter < 0) = 0;
G2_length = polyval(SG2M_timer,start_size);

current_time = 0;
cell_list_G2 = 1:n_cells;
size_at_div = zeros(n_cells,1);
growth_post_start = zeros(n_cells,1);
duration_post_start = zeros(n_cells,1);
duration_post_start_G1 = zeros(n_cells,1);
mother_bud_mass_defect = zeros(n_cells,1);

% Simulation

while isempty(cell_list_G2) == 0
    
    current_time = current_time + 1;
    i = current_time;
    
    ['Current Time: ' num2str(current_time)]
    ['Cells Remaining: ' num2str(length(cell_list_G2))]
    
    for k = cell_list_G2
        
        growth_size = cell_v(k,i) + cell_bud(k,i);
        dv = (exp(constant_daughter_growth_rate) - 1) * growth_size;
        if dv < 0
            dv = 0;
        end
        
        if cell_cycle(k,i) == 0
            if post_start_G1(k,i) == 1
                if post_start_counter(k,i) > 0
                    cell_v(k,i+1) = cell_v(k,i) + dv;
                    cell_cycle(k,i+1) = 0; %stays in G1
                    cell_bud(k,i+1) = 0;
                    post_start_G1(k,i+1) = 1;
                    post_start_counter(k,i+1) = post_start_counter(k,i) - 1;
                elseif post_start_counter(k,i) == 0
                    mother_bud_mass_defect(k) = ...
                        bud_mass_correction / G2_length(k);
                    cell_v(k,i+1) = cell_v(k,i) + mother_bud_mass_defect(k);
                    cell_bud(k,i+1) = cell_bud(k,i) + ...
                        dv - mother_bud_mass_defect(k); %grow bud
                    cell_cycle(k,i+1) = 1; %change to G2
                    post_start_G1(k,i+1) = 0;
                    post_start_counter(k,i+1) = 0;
                    
                    duration_post_start_G1(k) = current_time;
                end
            end
            
        elseif cell_cycle(k,i) == 1 %if cells in G2
            
            div_prob = G2_lambda_2d(1)*cell_bud(k,i) + ...
                G2_lambda_2d(2)*start_size(k) + G2_lambda_2d(3);
            
            div_prob = 1 - exp(-div_prob*(1));
            
            if div_prob < 0
                div_prob = 0;
            elseif div_prob > 1
                div_prob = 1;
            end
            
            G2_counter(k,i) = ...
                randsample([0 1], 1, true, [div_prob, 1 - div_prob]);
            
            if G2_counter(k,i) > 0 %growth in G2
                cell_v(k,i+1) = cell_v(k,i) + mother_bud_mass_defect(k); %bud grows, mother does not
                cell_bud(k,i+1) = cell_bud(k,i) + ...
                    dv - mother_bud_mass_defect(k); %grow bud
                cell_cycle(k,i+1) = 1; %Stay in G2
                post_start_G1(k,i+1) = 0;
                post_start_counter(k,i+1) = 0;
            elseif G2_counter(k,i) == 0
                
                size_at_div(k) = cell_v(k,i) + cell_bud(k,i);
                growth_post_start(k) = size_at_div(k) - start_size(k);
                duration_post_start(k) = current_time;
                
                cell_list_G2(cell_list_G2 == k) = [];
                
            end
        end
    end
end

total_growth = G1_growth + growth_post_start;

%% Saving Module

save('20170403 Numerical Analaysis Adder.mat',...
    'size_at_div','growth_post_start','duration_post_start',...
    'start_size','size_at_birth','size_at_start','G1_duration',...
    'G1_growth','total_growth')

%% Analysis and Plotting

% G1 Growth Characteristics

n_bins = 50;
n_bins_data = 25;

load('ActmCitr_20161220.mat')

[y_total_growth, x_birth_size] = KS_bindata_mean_20141016(...
    size_at_birth,total_growth,n_bins);

[y_total_growth_err, x_birth_size_err] = KS_bindata_std_20141016(...
    size_at_birth,total_growth,n_bins);

[y_total_growth_measured, x_birth_size_measured] = ...
    KS_bindata_mean_20141016(volumebirth_act,...
    deltavolume_fullcellcycle_act,n_bins_data);

[y_total_growth_measured_err, x_birth_size_measured_err] = ...
    KS_bindata_std_20141016(volumebirth_act,...
    deltavolume_fullcellcycle_act,n_bins_data);

figure('Name','Growth in G1')
hold on
shadedErrorBar(x_birth_size,y_total_growth,y_total_growth_err,'r',.5)
shadedErrorBar(x_birth_size_measured,y_total_growth_measured,...
    y_total_growth_measured_err,'b',.5)
xlabel('Size at Birth (AU)')
ylabel('Growth During G1 (AU)')
axis([0 inf 0 inf])
hold off
