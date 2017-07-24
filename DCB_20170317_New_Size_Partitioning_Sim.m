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

%Define simulation length parameters; time/iteration = 1 minute

start_phase_length = 2000;      %number of iterations in simulation
control_pop_size = 1;           %0 = like flask growth, 1 = like chemostat
generate_random_pop_start = 1;  %0 = random seed using parameters below; 1 = load seed pop; 2 = load and seed with contraints
start_phase_pop_cap = 10000; %maximum number of cells during steady state phase
sim_stage = 2;  %Seeding (1) or running (2) simulation; (0) does not save
max_seed_size = 100;
min_seed_size = 40;
n_cells = 10;          %starting number of cells

pop_mu_i = 3e4;          %femtoliters
pop_sigma_i = .17 * pop_mu_i;   %femptoliters

save(['Starting Parameters Temp.mat'], 'start_phase_length',...
    'control_pop_size','generate_random_pop_start','sim_stage',...
    'max_seed_size','min_seed_size','n_cells');

%% Cell Parameters

% Daughters

G1_lambda = coeffvalues(min_fit1_act);

bud_prob_beta_1 = G1_lambda(1);      %per minute;
bud_prob_M0_1 = -G1_lambda(2)/bud_prob_beta_1;      %fl

G2_lambda_2d = exp_2d_coeff;

constant_growth_rate = log(all_fit_act(1)*3 + 1)/3; %min_a_act; %per min
constant_daughter_growth_rate = log(daughter_growth_fit_act(1)*3 + 1)/3;
bud_mass_correction = d_b_diff;

poststart_G1_timer = PSG1_timer_fit;
SG2M_timer = SG2M_timer_fit;

% Mothers

mother_constant_growth_rate = constant_growth_rate;
mother_G1_growth_rate = mother_g1_gr_fit;
mother_lambda = exp_2d_coeff_mother;

mother_G2_timer = mother_G2_dur;
mother_G1_timer = mother_G1_fit;
mother_bud_mass_correction = m_b_diff;

% Define protein production parameters
p1_syn_G2 = 0;              %AU/min constant portion of synthesis rate
p1_syn_rate_G2 = .5;          %AU/(min*fl) scaling portion of synthesis rate
p1_syn_G1 = 0;              %AU/min constant portion of synthesis rate
p1_syn_rate_G1 = 0;          %AU/(min*fl) scaling portion of synthesis rate
p2_syn_G2 = .5;              %AU/min constant portion of synthesis rate
p2_syn_rate_G2 = 0;          %AU/(min*fl) scaling portion of synthesis rate
p2_syn_G1 = 0;              %AU/min constant portion of synthesis rate
p2_syn_rate_G1 = 0;          %AU/(min*fl) scaling portion of synthesis rate
partition1 = 0;          %partitioning behavior, 0 means proportional, 1 fixed
partition2 = 0;          %partitioning behavior, 0 means proportional, 1 fixed

%%
%Build starting population

if generate_random_pop_start == 0
    
    mu = log((pop_mu_i^2)/sqrt(pop_sigma_i+pop_mu_i^2));    %values for log-normal
    sigma = sqrt(log(pop_sigma_i/(pop_mu_i^2)+1));  %values for log-normal
    cell_v = lognrnd(mu,sigma,n_cells,1);    %sample from log-normal volume dist
    cell_p1 = zeros(n_cells,1);  %begin with empty cells
    cell_p2 = zeros(n_cells,1);  %begin with empty cells
    cell_cycle = zeros(n_cells,1);  %begin with G1 cells; G1 = 0, G2 = 1
    cell_bud = zeros(n_cells,1);  %begin with no buds
    G2_counter = zeros(n_cells,1);  %cells count down G2
    mother_daughter = zeros(n_cells,1); %begin with all daughters; daughter = 0, mother = 1
    mother_G1_counter = zeros(n_cells,1); %all cells start in G1 as daughters
    post_start_G1 = zeros(n_cells,1);
    post_start_counter = zeros(n_cells,1);
    start_size = zeros(n_cells,1);
    G2_length_daughter = []; %time in G2 (if buds)
    mother_counter = [];
    
elseif generate_random_pop_start == 1
    
    load('Seed Data Size Partitioning Sim.mat');
    load('Starting Parameters Temp.mat');
    delete('Starting Parameters Temp.mat');
    seed_time = length(cell_v(1,:));
    choose_cells = randsample(length(cell_v(:,seed_time)),n_cells); %select cells from previous population
    cell_v = cell_v(choose_cells,seed_time);
    cell_p1 = cell_p1(choose_cells,seed_time);  %begin with empty cells
    cell_p2 = cell_p2(choose_cells,seed_time);  %begin with empty cells
    cell_cycle = cell_cycle(choose_cells,seed_time);  %begin with G1 cells; G1 = 0, G2 = 1
    cell_bud = cell_bud(choose_cells,seed_time);  %begin with no buds
    G2_counter = G2_counter(choose_cells,seed_time);  %cells count down G2
    mother_daughter = mother_daughter(choose_cells,seed_time); %begin with all daughters; daughter = 0, mother = 1
    mother_G1_counter = mother_G1_counter(choose_cells,seed_time); %all cells start in G1 as daughters
    post_start_G1 = post_start_G1(choose_cells,seed_time);
    post_start_counter = post_start_counter(choose_cells,seed_time);
    mother_bud_mass_defect = mother_bud_mass_defect(choose_cells);
    start_size = start_size(choose_cells);
    mother_counter = mother_counter(choose_cells);
    G2_length_mother = G2_length_mother(choose_cells);
    
end

%% Simulation

%Grow cells
%Cells first choose between bud, divide, grow; cells whose divide counters
%are greater than 0 will grow bud, else, divide; cells who are in G1 and
%don't get a bud draw will grow mother; cells with bud draws will initiate
%bud and grow bud

%Diagnostic Parameters
pop_vol_mean = zeros(start_phase_length,1);
pop_vol_sigma = zeros(start_phase_length,1);
pop_p1 = zeros(start_phase_length,1);
pop_p1_sigma = zeros(start_phase_length,1);
pop_p2 = zeros(start_phase_length,1);
pop_p2_sigma = zeros(start_phase_length,1);
cell_num = zeros(start_phase_length,1);
mother_fraction = zeros(start_phase_length,1);
pop_cell_cycle = zeros(start_phase_length,1);

bud_prob = [];  %chance of budding (volume dependent)
bud = [];   %does this cell bud (if in G1)?

for i = 1:start_phase_length
    
    ['Current time: ' num2str(i)]
    
    %Update diagnostics
    pop_vol_mean(i) = mean(cell_v(:,i));
    pop_vol_sigma(i) = std(cell_v(:,i));
    pop_p1(i) = mean(cell_p1(:,i));
    pop_p1_sigma(i) = std(cell_p1(:,i));
    pop_p2(i) = mean(cell_p2(:,i));
    pop_p2_sigma(i) = std(cell_p2(:,i));
    cell_num(i) = length(cell_v(:,i));
    mother_fraction(i) = mean(mother_daughter(:,i));
    pop_cell_cycle(i) = mean(cell_cycle(:,i));
    
    cell_num(i)
    
    %Begin simulation - define global event stats
    bud_prob = bud_prob_beta_1.*(cell_v(:,i) - bud_prob_M0_1); %piecewise linear fits from data
    bud_prob(bud_prob < 0) = 0;
    bud_prob = 1 - exp(-bud_prob*(1));
    
    bud = zeros(length(bud_prob),1);
    bud = transpose(bud);
    
    for k = 1:length(cell_v(:,i))   %for all cells
        growth_size = cell_v(k,i) + cell_bud(k,i);
        
        if mother_daughter(k,i) == 0 %if cells are newborn
            dv = (exp(constant_growth_rate) - 1) * growth_size;
            dv(dv < 0) = 0;
            if cell_cycle(k,i) == 0 %if cells are in G1
                if post_start_G1(k,i) == 0 %if cells are pre-Start
                    
                    bud(k) = randsample([0 1], 1, true, [1-bud_prob(k), ...
                        bud_prob(k)]);
                    
                    if bud(k) == 0
                        cell_v(k,i+1) = cell_v(k,i) + dv;
                        cell_cycle(k,i+1) = 0; %stays in G1
                        cell_bud(k,i+1) = 0;
                        cell_p1(k,i+1) = cell_p1(k,i) + (cell_v(k,i)*p1_syn_rate_G1 + p1_syn_G1);
                        cell_p2(k,i+1) = cell_p2(k,i) + (cell_v(k,i)*p2_syn_rate_G1 + p2_syn_G1);
                        mother_daughter(k,i+1) = mother_daughter(k,i);
                        G2_counter(k,i+1) = 0;
                        mother_G1_counter(k,i+1) = 0;
                        post_start_G1(k,i+1) = 0;
                        post_start_counter(k,i+1) = 0;
                        mother_counter(k) = 0;
                    elseif bud(k) == 1
                        cell_v(k,i+1) = cell_v(k,i) + dv;
                        cell_cycle(k,i+1) = 0; %stays in G1
                        cell_bud(k,i+1) = 0;
                        cell_p1(k,i+1) = cell_p1(k,i) + (cell_v(k,i)*p1_syn_rate_G1 + p1_syn_G1);
                        cell_p2(k,i+1) = cell_p2(k,i) + (cell_v(k,i)*p2_syn_rate_G1 + p2_syn_G1);
                        mother_daughter(k,i+1) = mother_daughter(k,i);
                        G2_counter(k,i+1) = 1;
                        mother_G1_counter(k,i+1) = 0;
                        post_start_G1(k,i+1) = 1;
                        post_start_counter(k,i+1) = ...
                            round(polyval(poststart_G1_timer,cell_v(k,i)));
                        if post_start_counter(k,i+1) < 0;
                            post_start_counter(k,i+1) = 0;
                        end
                        start_size(k) = cell_v(k,i);
                        G2_length(k) = polyval(SG2M_timer,start_size(k));
                        if G2_length(k) < 0;
                            G2_length(k) = 0;
                        end
                        mother_bud_mass_defect(k) = bud_mass_correction / G2_length(k);
                        mother_counter(k) = 0;
                    end
                elseif post_start_G1(k,i) == 1 %if cells are post-Start
                    if post_start_counter(k,i) > 0
                        cell_v(k,i+1) = cell_v(k,i) + dv;
                        cell_cycle(k,i+1) = 0; %stays in G1
                        cell_bud(k,i+1) = 0;
                        cell_p1(k,i+1) = cell_p1(k,i) + (cell_v(k,i)*p1_syn_rate_G1 + p1_syn_G1);
                        cell_p2(k,i+1) = cell_p2(k,i) + (cell_v(k,i)*p2_syn_rate_G1 + p2_syn_G1);
                        mother_daughter(k,i+1) = mother_daughter(k,i);
                        G2_counter(k,i+1) = G2_counter(k,i);
                        mother_G1_counter(k,i+1) = 0;
                        post_start_G1(k,i+1) = 1;
                        post_start_counter(k,i+1) = post_start_counter(k,i) - 1;
                        mother_counter(k) = 0;
                    elseif post_start_counter(k,i) == 0
                        cell_v(k,i+1) = cell_v(k,i); %bud grows, mother does not
                        cell_cycle(k,i+1) = 1; %change to G2
                        G2_counter(k,i+1) = G2_counter(k,i); %randomly select length of G2
                        cell_bud(k,i+1) = dv;
                        cell_p1(k,i+1) = cell_p1(k,i) + (cell_v(k,i)*p1_syn_rate_G2 + p1_syn_G2);
                        cell_p2(k,i+1) = cell_p2(k,i) + (cell_v(k,i)*p2_syn_rate_G2 + p2_syn_G2);
                        mother_daughter(k,i+1) = mother_daughter(k,i);
                        mother_G1_counter(k,i+1) = 0;
                        post_start_G1(k,i+1) = 0;
                        post_start_counter(k,i+1) = 0;
                        mother_counter(k) = 0;
                    end
                end
            elseif  cell_cycle(k,i) == 1 %if cells are in G2
                div_prob = G2_lambda_2d(1)*cell_bud(k,i) + ...
                    G2_lambda_2d(2)*start_size(k) + G2_lambda_2d(3);
                
                if div_prob < 0
                    div_prob = 0;
                end
                
                div_prob = 1 - exp(-div_prob*(1));
                
                G2_counter(k,i) = ...
                    randsample([0 1], 1, true, [div_prob, 1 - div_prob]);
                
                if G2_counter(k,i) > 0 %growth in G2
                    cell_v(k,i+1) = cell_v(k,i) + mother_bud_mass_defect(k); %bud grows, mother does not
                    cell_bud(k,i+1) = cell_bud(k,i) + ...
                        dv - mother_bud_mass_defect(k); %grow bud
                    cell_p1(k,i+1) = cell_p1(k,i) + ...
                        ((cell_v(k,i)+cell_bud(k,i))*p1_syn_rate_G2 + p1_syn_G2); %synthesize protein
                    cell_p2(k,i+1) = cell_p2(k,i) + ...
                        ((cell_v(k,i)+cell_bud(k,i))*p2_syn_rate_G2 + p2_syn_G2); %synthesize protein
                    G2_counter(k,i+1) = G2_counter(k,i) - 1; %shorten remaining G2
                    cell_cycle(k,i+1) = 1; %Stay in G2
                    mother_daughter(k,i+1) = mother_daughter(k,i);
                    mother_G1_counter(k,i+1) = 0;
                    post_start_G1(k,i+1) = 0;
                    post_start_counter(k,i+1) = 0;
                    mother_counter(k) = mother_counter(k) + 1;
                elseif G2_counter(k,i) == 0
                    new_cell_num = length(cell_v(:,i))+1;
                    cell_v(new_cell_num,i+1) = cell_bud(k,i); %Make new daughter
                    cell_v(k,i+1) = cell_v(k,i); %No growth during cytokinesis
                    cell_bud(new_cell_num,i+1) = 0; %No bud on daughter
                    cell_cycle(new_cell_num,i+1) = 0; %New cell in G1
                    cell_bud(k,i+1) = 0; %Seperate bud from mother
                    cell_cycle(k,i+1) = 0; %Return mother to G1
                    mother_daughter(k,i+1) = 1; %Not first time mother
                    mother_daughter(new_cell_num,i+1) = 0; %New never budded daughter
                    mother_G1_counter(k,i+1) = ...
                        round(polyval(mother_G1_timer,cell_v(k,i)));
                    if mother_G1_counter(k,i+1) < 0
                        mother_G1_counter(k,i+1) = 0;
                    end
                    mother_G1_counter(new_cell_num,i+1) = 0;
                    G2_counter(k,i+1) = 0;
                    G2_counter(new_cell_num,i+1) = 0;
                    post_start_G1(k,i+1) = 0;
                    post_start_counter(k,i+1) = 0;
                    post_start_G1(new_cell_num,i+1) = 0;
                    post_start_counter(new_cell_num,i+1) = 0;
                    mother_bud_mass_defect(new_cell_num) = 0;
                    mother_bud_mass_defect(k) = 0;
                    start_size(k) = cell_v(k,i+1);
                    start_size(new_cell_num) = 0;
                    mother_counter(k) = 0;
                    mother_counter(new_cell_num) = 0;
                    G2_length_mother(k) = ...
                        polyval(mother_G2_timer,cell_v(k,i));
                    G2_length_mother(new_cell_num) = 0;
                    cell_p1(k,i+1) = (1 - partition1) * ...
                        cell_p1(k,i) * (cell_v(k,i+1) / ...
                        (cell_v(k,i+1) + cell_bud(k,i))) + ...
                        partition1 * cell_p1(k,i) * .5; %allot protein to mother
                    cell_p1(new_cell_num,i+1) = (1 - partition1) * ...
                        cell_p1(k,i) * (cell_bud(k,i) / (cell_v(k,i+1) + ...
                        cell_bud(k,i))) + partition1 * ...
                        cell_p1(k,i) * .5; %allot protein to daughter
                    cell_p2(k,i+1) = (1 - partition2) * ...
                        cell_p2(k,i) * (cell_v(k,i+1) / ...
                        (cell_v(k,i+1) + cell_bud(k,i))) + ...
                        partition2 * cell_p2(k,i) * .5; %allot protein to mother
                    cell_p2(new_cell_num,i+1) = (1 - partition2) * ...
                        cell_p2(k,i) * (cell_bud(k,i) / (cell_v(k,i+1) + ...
                        cell_bud(k,i))) + partition2 * ...
                        cell_p2(k,i) * .5; %allot protein to daughter
                    M_D_parent_at_birth(new_cell_num) = mother_daughter(k,i);
                    
                end
                
            end
        elseif mother_daughter(k,i) == 1 % if cells are 2nd gen
            growth_size = cell_v(k,i) + cell_bud(k,i);
            if cell_cycle(k,i) == 0
                dv = ((mother_g1_gr_fit(1)*growth_size + ...
                    mother_g1_gr_fit(2))*exp(mother_g1_gr_fit(1)) - ...
                    mother_g1_gr_fit(2))/mother_g1_gr_fit(1) - growth_size;
                dv(dv < 0) = 0;
                if mother_G1_counter(k,i) > 0 %if cell does not bud this time
                    cell_v(k,i+1) = cell_v(k,i) + dv;
                    cell_cycle(k,i+1) = 0; %stays in G1
                    cell_bud(k,i+1) = 0;
                    cell_p1(k,i+1) = cell_p1(k,i) + (cell_v(k,i)*p1_syn_rate_G1 + p1_syn_G1);
                    cell_p2(k,i+1) = cell_p2(k,i) + (cell_v(k,i)*p2_syn_rate_G1 + p2_syn_G1);
                    G2_counter(k,i+1) = 0;
                    mother_daughter(k,i+1) = mother_daughter(k,i);
                    mother_G1_counter(k,i+1) = mother_G1_counter(k,i) - 1;
                    post_start_G1(k,i+1) = 0;
                    post_start_counter(k,i+1) = 0;
                    mother_counter(k) = mother_counter(k) + 1;
                    
                elseif mother_G1_counter(k,i) == 0 %if cell buds
                    cell_v(k,i+1) = cell_v(k,i); %bud grows, mother does not
                    cell_cycle(k,i+1) = 1; %change to G2
                    G2_counter(k,i+1) = 1; %randomly select length of G2
                    cell_bud(k,i+1) = dv;
                    cell_p1(k,i+1) = cell_p1(k,i) + (cell_v(k,i)*p1_syn_rate_G2 + p1_syn_G2);
                    cell_p2(k,i+1) = cell_p2(k,i) + (cell_v(k,i)*p2_syn_rate_G2 + p2_syn_G2);
                    mother_daughter(k,i+1) = mother_daughter(k,i);
                    mother_G1_counter(k,i+1) = 0;
                    post_start_G1(k,i+1) = 0;
                    post_start_counter(k,i+1) = 0;
                    mother_bud_mass_defect(k) = mother_bud_mass_correction / G2_length_mother(k);
                    mother_counter(k) = mother_counter(k) + 1;
                    
                end
                
            elseif cell_cycle(k,i) == 1
                
                dv = (exp(constant_growth_rate) - 1) * growth_size;
                dv(dv < 0) = 0;
                
                div_prob = mother_lambda(1)*cell_bud(k,i) + ...
                    mother_lambda(2)*start_size(k) + mother_lambda(3);
%                 div_prob = mother_lambda(1) * mother_counter(k) ...
%                     + mother_lambda(2); 
%                 div_prob = mother_lambda(1) * cell_bud(k,i) ...
%                     + mother_lambda(2);                
                div_prob = 1 - exp(-div_prob*(1));
                div_prob(div_prob < 0) = 0;
                
                G2_counter(k,i) = ...
                    randsample([0 1], 1, true, [div_prob, 1 - div_prob]);
                
                if G2_counter(k,i) > 0 %growth in G2
                    cell_v(k,i+1) = cell_v(k,i) + mother_bud_mass_defect(k); %bud grows, mother does not
                    cell_bud(k,i+1) = cell_bud(k,i) + ...
                        dv - mother_bud_mass_defect(k); %grow bud
                    cell_p1(k,i+1) = cell_p1(k,i) + ...
                        ((cell_v(k,i)+cell_bud(k,i))*p1_syn_rate_G2 + p1_syn_G2); %synthesize protein
                    cell_p2(k,i+1) = cell_p2(k,i) + ...
                        ((cell_v(k,i)+cell_bud(k,i))*p2_syn_rate_G2 + p2_syn_G2); %synthesize protein
                    G2_counter(k,i+1) = G2_counter(k,i); %shorten remaining G2
                    cell_cycle(k,i+1) = 1; %Stay in G2
                    mother_daughter(k,i+1) = mother_daughter(k,i);
                    mother_G1_counter(k,i+1) = 0;
                    post_start_G1(k,i+1) = 0;
                    post_start_counter(k,i+1) = 0;
                    mother_counter(k) = mother_counter(k) + 1;
                    
                elseif G2_counter(k,i) == 0
                    new_cell_num = length(cell_v(:,i))+1;
                    cell_v(new_cell_num,i+1) = cell_bud(k,i); %Make new daughter
                    cell_v(k,i+1) = cell_v(k,i); %No growth during cytokinesis
                    cell_bud(new_cell_num,i+1) = 0; %No bud on daughter
                    cell_cycle(new_cell_num,i+1) = 0; %New cell in G1
                    cell_bud(k,i+1) = 0; %Seperate bud from mother
                    cell_cycle(k,i+1) = 0; %Return mother to G1
                    mother_daughter(k,i+1) = 1; %Not first time mother
                    mother_daughter(new_cell_num,i+1) = 0; %New never budded daughter
                    mother_G1_counter(k,i+1) = ...
                        round(polyval(mother_G1_timer,cell_v(k,i)));
                    if mother_G1_counter(k,i+1) < 0
                        mother_G1_counter(k,i+1) = 0;
                    end
                    mother_G1_counter(new_cell_num,i+1) = 0;
                    G2_counter(k,i+1) = 0;
                    G2_counter(new_cell_num,i+1) = 0;
                    post_start_G1(k,i+1) = 0;
                    post_start_counter(k,i+1) = 0;
                    post_start_G1(new_cell_num,i+1) = 0;
                    post_start_counter(new_cell_num,i+1) = 0;
                    mother_bud_mass_defect(new_cell_num) = 0;
                    mother_bud_mass_defect(k) = 0;
                    start_size(k) = cell_v(k,i+1);
                    start_size(new_cell_num) = 0;
                    mother_counter(k) = 0;
                    mother_counter(new_cell_num) = 0;
                    G2_length_mother(k) = ...
                        polyval(mother_G2_timer,cell_v(k,i));
                    if G2_length_mother(k) < 0;
                            G2_length_mother(k) = 0;
                        end
                    G2_length_mother(new_cell_num) = 0;
                    
                    cell_p1(k,i+1) = (1 - partition1) * ...
                        cell_p1(k,i) * (cell_v(k,i+1) / ...
                        (cell_v(k,i+1) + cell_bud(k,i))) + ...
                        partition1 * cell_p1(k,i) * .5; %allot protein to mother
                    cell_p1(new_cell_num,i+1) = (1 - partition1) * ...
                        cell_p1(k,i) * (cell_bud(k,i) / (cell_v(k,i+1) + ...
                        cell_bud(k,i))) + partition1 * ...
                        cell_p1(k,i) * .5; %allot protein to daughter
                    cell_p2(k,i+1) = (1 - partition2) * ...
                        cell_p2(k,i) * (cell_v(k,i+1) / ...
                        (cell_v(k,i+1) + cell_bud(k,i))) + ...
                        partition2 * cell_p2(k,i) * .5; %allot protein to mother
                    cell_p2(new_cell_num,i+1) = (1 - partition2) * ...
                        cell_p2(k,i) * (cell_bud(k,i) / (cell_v(k,i+1) + ...
                        cell_bud(k,i))) + partition2 * ...
                        cell_p2(k,i) * .5; %allot protein to daughter
                    M_D_parent_at_birth(new_cell_num) = mother_daughter(k,i);
                    
                end
            end
        end
    end
    
    if control_pop_size == 1
        select_cells = [];
        if length(cell_v(:,i+1)) > start_phase_pop_cap %caps population by randomly removing excess cells
            select_cells = randsample(length(cell_v(:,i+1)),start_phase_pop_cap);
            cell_v = cell_v(select_cells,:);
            cell_p1 = cell_p1(select_cells,:);
            cell_p2 = cell_p2(select_cells,:);
            cell_cycle = cell_cycle(select_cells,:);
            cell_bud = cell_bud(select_cells,:);
            G2_counter = G2_counter(select_cells,:);
            mother_G1_counter = mother_G1_counter(select_cells,:);
            mother_daughter = mother_daughter(select_cells,:);
            post_start_G1 = post_start_G1(select_cells,:);
            post_start_counter = post_start_counter(select_cells,:);
            mother_bud_mass_defect(select_cells);
            M_D_parent_at_birth = M_D_parent_at_birth(select_cells);
            start_size = start_size(select_cells);
            mother_counter = mother_counter(select_cells);
            G2_length_mother = G2_length_mother(select_cells);
            
        end
    end
    
    if length(find(cell_v(:,i+1) == 0)) > 0 % Negative growth is an error condition
        ['Negative Growth']
        return
    end
end


if sim_stage == 1
    save(['Seed Data Size Partitioning Sim.mat'],...
        'cell_bud', 'cell_cycle', 'cell_p1', 'cell_p2', ...
        'cell_v', 'G2_counter', 'mother_daughter', ...
        'mother_G1_counter', 'post_start_counter', 'post_start_G1',...
        'mother_bud_mass_defect','start_size','G2_length_mother',...
        'mother_counter');
end

%%
%Investigative Plot Data Processing

steady_state_time = 1000; %Time to start investigating cells

cell_cycle_change = [];
START_change = [];
birth_time = [];
first_daughter_bud = [];
first_daughter_cytokinesis = [];
G1_growth = [];
size_at_birth = [];
size_at_bud =  [];
cell_cycle_growth = [];
p1_birth_amount = [];
p2_birth_amount = [];
p1_over_cell_cycle = [];
p2_over_cell_cycle = [];
p1_at_bud = [];
p2_at_bud = [];
G1_amount_volume = [];
G1_amount_p1 = [];
G1_amount_p2 = [];
first_daughter_start = [];
bud_size_at_cytokinesis = [];
volume_at_START = [];
n_data_cells = 0;
cell_cycle_length = [];

for cell = 1:length(cell_v(:,1))
    cell_cycle_change(cell,:) = diff(cell_cycle(cell,:));
    START_change(cell,:) = diff(post_start_G1(cell,:));
    if isempty(find(cell_cycle_change(cell, :) == 1, 1)) == 0
        if isempty(find(cell_cycle_change(cell, :) == -1, 1)) == 0
            if isempty(find(START_change(cell,:) == 1, 1)) == 0
                n_data_cells = n_data_cells + 1;
                birth_time(n_data_cells) = find(cell_v(cell,:) > 0, 1);
                first_daughter_bud(n_data_cells) = find(cell_cycle_change(cell, :) == 1, 1);
                first_daughter_cytokinesis(n_data_cells) = find(cell_cycle_change(cell, :) == -1, 1);
                first_daughter_start(n_data_cells) = find(START_change(cell,:) == 1, 1);
                cell_cycle_length(n_data_cells) = first_daughter_cytokinesis(n_data_cells) - birth_time(n_data_cells);
                G1_growth(n_data_cells) = cell_v(cell,first_daughter_bud(n_data_cells)) - ...
                    cell_v(cell,birth_time(n_data_cells));
                size_at_birth(n_data_cells) = cell_v(cell,birth_time(n_data_cells));
                size_at_bud(n_data_cells) = cell_v(cell,first_daughter_bud(n_data_cells));
                cell_cycle_growth(n_data_cells) = cell_v(cell,first_daughter_cytokinesis(n_data_cells)) - ...
                    cell_v(cell,birth_time(n_data_cells)) + cell_bud(cell,first_daughter_cytokinesis(n_data_cells));
                p1_birth_amount(n_data_cells) = cell_p1(cell, birth_time(n_data_cells));
                p2_birth_amount(n_data_cells) = cell_p2(cell, birth_time(n_data_cells));
                p1_at_bud(n_data_cells) = cell_p1(cell, first_daughter_bud(n_data_cells));
                p2_at_bud(n_data_cells) = cell_p2(cell, first_daughter_bud(n_data_cells));
                p1_over_cell_cycle(n_data_cells) = cell_p1(cell,first_daughter_cytokinesis(n_data_cells)) - ...
                    cell_p1(cell,birth_time(n_data_cells));
                p2_over_cell_cycle(n_data_cells) = cell_p2(cell,first_daughter_cytokinesis(n_data_cells)) - ...
                    cell_p2(cell,birth_time(n_data_cells));
                G1_amount_volume = [G1_amount_volume, ...
                    cell_v(cell,birth_time(n_data_cells):first_daughter_bud(n_data_cells))];
                G1_amount_p1 = [G1_amount_p1, cell_p1(cell, ...
                    birth_time(n_data_cells):first_daughter_bud(n_data_cells))];
                G1_amount_p2 = [G1_amount_p2, cell_p2(cell, ...
                    birth_time(n_data_cells):first_daughter_bud(n_data_cells))];
                bud_size_at_cytokinesis(n_data_cells) = cell_bud(cell, ...
                    first_daughter_cytokinesis(n_data_cells) - 1);
                volume_at_START(n_data_cells) = cell_v(cell,first_daughter_start(n_data_cells));
                M_D_at_birth_processed(n_data_cells) = M_D_parent_at_birth(cell);
            end
        end
    end
end

first_gen_mother_birth = [];
first_gen_mother_bud = [];
first_gen_mother_cytokinesis = [];
mother_birth_size = [];
mother_adder = [];
n_mother_cells = 0;
for cell = 1:length(cell_v(:,1))
    if length(find(cell_cycle_change(cell, :) == 1, 2)) > 1
        if length(find(cell_cycle_change(cell, :) == -1, 2)) > 1
            n_mother_cells = n_mother_cells + 1;
            m_birth = find(cell_cycle_change(cell, :) == -1, 1);
            first_gen_mother_birth(n_mother_cells) = ...
                m_birth(1) + 1;
            m_bud = find(cell_cycle_change(cell, :) == 1, 2);
            first_gen_mother_bud(n_mother_cells) = ...
                m_bud(2);
            m_cytokinesis = find(cell_cycle_change(cell, :) == -1, 2);
            first_gen_mother_cytokinesis(n_mother_cells) = ...
                m_cytokinesis(2);
            mother_birth_size(n_mother_cells) = ...
                cell_v(cell,first_gen_mother_birth(n_mother_cells)) + ...
                cell_bud(cell,first_gen_mother_birth(n_mother_cells));
            mother_cytokinesis_size(n_mother_cells) = ...
                cell_v(cell,first_gen_mother_cytokinesis(n_mother_cells)) + ...
                cell_bud(cell,first_gen_mother_cytokinesis(n_mother_cells));
            mother_adder(n_mother_cells) = ...
                mother_cytokinesis_size(n_mother_cells) - ...
                mother_birth_size(n_mother_cells);
            mother_bud_size(n_mother_cells) = ...
                cell_v(cell,first_gen_mother_bud(n_mother_cells)) + ...
                cell_bud(cell,first_gen_mother_bud(n_mother_cells));
        end
    end
end

cells_after_steady_state = find(birth_time > steady_state_time);
mothers_after_steady_state = find(first_gen_mother_birth > steady_state_time);
first_gen_mother_birth = first_gen_mother_birth(mothers_after_steady_state);
first_gen_mother_bud = first_gen_mother_bud(mothers_after_steady_state);
first_gen_mother_cytokinesis = first_gen_mother_cytokinesis(mothers_after_steady_state);
birth_time = birth_time(cells_after_steady_state);
first_daughter_bud = first_daughter_bud(cells_after_steady_state);
first_daughter_cytokinesis = first_daughter_cytokinesis(cells_after_steady_state);
G1_growth = G1_growth(cells_after_steady_state);
size_at_birth = size_at_birth(cells_after_steady_state);
size_at_bud =  size_at_bud(cells_after_steady_state);
cell_cycle_growth = cell_cycle_growth(cells_after_steady_state);
p1_birth_amount = p1_birth_amount(cells_after_steady_state);
p2_birth_amount = p2_birth_amount(cells_after_steady_state);
p1_over_cell_cycle = p1_over_cell_cycle(cells_after_steady_state);
p2_over_cell_cycle = p2_over_cell_cycle(cells_after_steady_state);
p1_at_bud = p1_at_bud(cells_after_steady_state);
p2_at_bud = p2_at_bud(cells_after_steady_state);
bud_size_at_cytokinesis = bud_size_at_cytokinesis(cells_after_steady_state);
volume_at_START = volume_at_START(cells_after_steady_state);
first_daughter_start = first_daughter_start(cells_after_steady_state);
M_D_at_birth_processed = M_D_at_birth_processed(cells_after_steady_state);


%Real data for comparison
load('ActmCitr_20161220.mat')
load('ActmCitr_fulldistribution_20161220.mat')
load('ActmCitr_mothers_20170117')

counter = 0;
G1_length_microscope = [];
for i = 1:length(act_cwhi5_bp(:,3))
    if act_cwhi5_bp(i,3) == 1
        G1_length_microscope = [G1_length_microscope, counter];
        counter = 0;
    elseif act_cwhi5_bp(i,3) == 0
        counter = counter + 3;
    end
end

if sim_stage == 2
    save('20170326 Size-Partitioning Simulation Results.mat')
end

%% Diagnostic Plots

%population volume over time
figure('Name','Volume')
hold on
shadedErrorBar(1:start_phase_length,pop_vol_mean,pop_vol_sigma)
hold off

%average cell cycle stage over time
figure('Name','Cell Cycle Synchrony')
hold on
plot(1:start_phase_length,pop_cell_cycle)
hold off

%% Analysis

%population v microscopy full distribution
figure('Name','Size Distribution - All Cells, Microscope')
hold on
ecdf(cell_v(:,length(cell_v(1,:))) + cell_bud(:,length(cell_v(1,:))), ...
    'bounds','on')
ecdf(allvolumes_act,'bounds','on')
axis([0 inf 0 inf])
xlabel('Volume (fl)')
ylabel('Cumulative Probability')
legend('Simulated','','','','Measured','','','')
hold off

%birth size real v sim, microscope
figure('Name', 'Size at Birth, Microscope')
hold on
ecdf(size_at_birth,'bounds','on')
ecdf(volumebirth_act,'bounds','on')
legend('Simulated','','','','Measured','','','')
hold off

%size @ START real v sim, microscope
figure('Name', 'Size at Start, Microscope')
hold on
ecdf(volume_at_START,'bounds','on')
ecdf(volumeSTART_act,'bounds','on')
legend('Simulated','','','','Measured','','','')
hold off

%size @ cytokinesis real v sim, microscope
figure('Name', 'Size at Cytokinesis, Microscope')
hold on
ecdf(size_at_birth + cell_cycle_growth,'bounds','on')
ecdf(volumebirth_act + transpose(deltavolume_fullcellcycle_act),'bounds','on')
legend('Simulated','','','','Measured','','','')
hold off

%Growth in cell cycle v size at birth
figure('Name','Adder Model Binned')
hold on
% scatter(size_at_birth, cell_cycle_growth)
[inc_mean_y,inc_mean_x] = KS_bindata_mean_20141016(size_at_birth(size_at_birth < 32000), ...
    cell_cycle_growth(size_at_birth < 32000), 15);
[inc_err_y,inc_err_x] = KS_bindata_std_20141016(size_at_birth(size_at_birth < 32000), ...
    cell_cycle_growth(size_at_birth < 32000), 15);
% errorbar(inc_mean_x, inc_mean_y, inc_err_y)
[inc_mean_y_real,inc_mean_x_real] = KS_bindata_mean_20141016(volumebirth_act(volumebirth_act < 32000), ...
    deltavolume_fullcellcycle_act(volumebirth_act < 32000), 15);
[inc_err_y_real,inc_err_x_real] = KS_bindata_std_20141016(volumebirth_act(volumebirth_act < 32000), ...
    deltavolume_fullcellcycle_act(volumebirth_act < 32000), 15);
shadedErrorBar(inc_mean_x, inc_mean_y, inc_err_y, 'r',.5)
shadedErrorBar(inc_mean_x_real, inc_mean_y_real, inc_err_y_real, 'b',.5)
axis([0 inf 0 inf])
xlabel('Size at Birth')
ylabel('Growth over the Cell Cycle')
legend('Simulated','','','','Measured','','','','Location','northwest')
hold off

%More size control
figure('Name','G1 Size Control')
hold on
% scatter(size_at_birth, size_at_bud)
axis([0 40000 0 50000])
[bud_y,birth_x] = KS_bindata_mean_20141016(size_at_birth, ...
    volume_at_START, 25);
[bud_err_y,birth_err_x] = KS_bindata_std_20141016(size_at_birth, ...
    volume_at_START, 25);
[bud_y_real,birth_x_real] = KS_bindata_mean_20141016(volumebirth_act, ...
    volumeSTART_act, 25);
[bud_err_y_real,birth_err_x_real] = KS_bindata_std_20141016(volumebirth_act, ...
    volumeSTART_act, 25);
shadedErrorBar(birth_x, bud_y, bud_err_y, 'r',.5)
shadedErrorBar(birth_x_real, bud_y_real, bud_err_y_real, 'b',.5)
xlabel('Size at Birth')
ylabel('Size at Start')
legend('Simulated','','','','Measured','','','')
hold off

figure('Name','Added Post-Start')
hold on
axis([0 50000 0 60000])
[bud_y,birth_x] = KS_bindata_mean_20141016(volume_at_START, ...
    size_at_birth + cell_cycle_growth - volume_at_START, 25);
[bud_err_y,birth_err_x] = KS_bindata_std_20141016(volume_at_START, ...
    size_at_birth + cell_cycle_growth - volume_at_START, 25);
[bud_y_real,birth_x_real] = KS_bindata_mean_20141016(volumeSTART_act, ...
    deltavolumepostSTART_act, 25);
[bud_err_y_real,birth_err_x_real] = KS_bindata_std_20141016(volumeSTART_act, ...
    deltavolumepostSTART_act, 25);
shadedErrorBar(birth_x, bud_y, bud_err_y, 'r',.5)
shadedErrorBar(birth_x_real, bud_y_real, bud_err_y_real, 'b',.5)
xlabel('Size at Start')
ylabel('Amount Added SG2M')
legend({'Simulated','','','','Measured','','',''},'Location','northwest')
hold off

figure('Name','Growth During G1')
hold on
axis([0 40000 0 20000])
[bud_y,birth_x] = KS_bindata_mean_20141016(size_at_birth, ...
    volume_at_START - size_at_birth, 25);
[bud_err_y,birth_err_x] = KS_bindata_std_20141016(size_at_birth, ...
    volume_at_START - size_at_birth, 25);
[bud_y_real,birth_x_real] = KS_bindata_mean_20141016(volumebirth_act, ...
    volumeSTART_act - volumebirth_act', 25);
[bud_err_y_real,birth_err_x_real] = KS_bindata_std_20141016(volumebirth_act, ...
    volumeSTART_act - volumebirth_act', 25);
% scatter(volumebirth_act, volumeSTART_act - volumebirth_act')
shadedErrorBar(birth_x, bud_y, bud_err_y, 'r',.5)
shadedErrorBar(birth_x_real, bud_y_real, bud_err_y_real, 'b',.5)
xlabel('Size at Start')
ylabel('Amount Added G1')
legend('Simulated','','','','Measured','','','')
hold off

% Mother Growth Check
figure('Name', 'Mother Growth Check')
hold on
[inc_mean_y,inc_mean_x] = KS_bindata_mean_20141016(mother_birth_size,...
    mother_cytokinesis_size - mother_birth_size, 15);
[inc_err_y,inc_err_x] = KS_bindata_std_20141016(mother_birth_size,...
    mother_cytokinesis_size - mother_birth_size, 15);
[inc_mean_y_real,inc_mean_x_real] = KS_bindata_mean_20141016(...
    mothersize_beginG1_act, ...
    mothersize_totalsizecytokinesis_act - mothersize_beginG1_act, 10);
[inc_err_y_real,inc_err_x_real] = KS_bindata_std_20141016(...
    mothersize_beginG1_act, ...
    mothersize_totalsizecytokinesis_act - mothersize_beginG1_act, 10);
shadedErrorBar(inc_mean_x, inc_mean_y, inc_err_y, 'r',.5)
shadedErrorBar(inc_mean_x_real, inc_mean_y_real, inc_err_y_real, 'b',.5)
axis([0 inf 0 inf])
xlabel('Size at Birth')
ylabel('Growth over the Cell Cycle')
legend('Simulated','','','','Measured','','','')
hold off

% Mother Growth Check
figure('Name', 'Mother G2 Growth Check')
hold on
[inc_mean_y,inc_mean_x] = KS_bindata_mean_20141016(mother_bud_size,...
    mother_cytokinesis_size - mother_bud_size, 15);
[inc_err_y,inc_err_x] = KS_bindata_std_20141016(mother_bud_size,...
    mother_cytokinesis_size - mother_bud_size, 15);
[inc_mean_y_real,inc_mean_x_real] = KS_bindata_mean_20141016(...
    mothersize_budemerge_act, ...
    mothersize_totalsizecytokinesis_act - mothersize_budemerge_act, 15);
[inc_err_y_real,inc_err_x_real] = KS_bindata_std_20141016(...
    mothersize_budemerge_act, ...
    mothersize_totalsizecytokinesis_act - mothersize_budemerge_act, 15);
shadedErrorBar(inc_mean_x, inc_mean_y, inc_err_y, 'r',.5)
shadedErrorBar(inc_mean_x_real, inc_mean_y_real, inc_err_y_real, 'b',.5)
axis([0 inf 0 inf])
xlabel('Size at Bud')
ylabel('Growth G2')
legend('Simulated','','','','Measured','','','')
hold off

% No Correlation Plot
figure('Name','Post v Pre Start growth')
hold on
[inc_mean_y,inc_mean_x] = KS_bindata_mean_20141016(volume_at_START - ...
    size_at_birth,...
    size_at_birth + cell_cycle_growth - volume_at_START, 15);
[inc_err_y,inc_err_x] = KS_bindata_std_20141016(volume_at_START - ...
    size_at_birth,...
    size_at_birth + cell_cycle_growth - volume_at_START, 15);
[inc_mean_y_real,inc_mean_x_real] = KS_bindata_mean_20141016(...
    volumeSTART_act - volumebirth_act', ...
    deltavolumepostSTART_act, 15);
[inc_err_y_real,inc_err_x_real] = KS_bindata_std_20141016(...
    volumeSTART_act - volumebirth_act', ...
    deltavolumepostSTART_act, 15);
shadedErrorBar(inc_mean_x, inc_mean_y, inc_err_y, 'r',.5)
shadedErrorBar(inc_mean_x_real, inc_mean_y_real, inc_err_y_real, 'b',.5)
axis([-inf inf 0 inf])
xlabel('G1 Growth')
ylabel('G2 Growth')
legend('Simulated','','','','Measured','','','')
hold off