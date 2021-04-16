% script to collate the abs error mat files from server

load('matrix_dk68_empirical_89_subs_MEG.mat')

edge_matrix_values = zeros(9, 3, 3);
mean_global_values = zeros(9, 3, 3);
Rspatime_values = zeros(9, 3, 3);
std_global_values = zeros(9, 3, 3);
esp_values = zeros(9, 3, 3);

for i = 1:9
    for j = 7:9
        load(sprintf('plan_b_matrix_MEG_values_supercrit_89_beta_%03d_G_%03d.mat',i, j)) % copy from other file
        edge_matrix_values(i, j-6, :) = edge_matrix_std_av_sim';
        mean_global_values(i, j-6, :) = mean_global_av_sim';
        Rspatime_values(i, j-6, :) = Rspatime_av_sim';
        std_global_values(i, j-6, :) = std_global_av_sim';
        esp_values(i, j-6, :) = EdgeSpaTimePredictability_av_sim';
    end
end


% need to find the working point for the boxplots and then define the
% surrogate

figure(1)

boxplot([edge_matrix_std_av(1:3)' squeeze(edge_matrix_values(9, 2, :))]);

figure(2)

boxplot([EdgeSpaTimePredictability_av(1:3) squeeze(esp_values(9, 1, :))]);


edge_matrix_value = mean(edge_matrix_values, 3);
mean_global_value = mean(mean_global_values, 3);
Rspatime_value = mean(Rspatime_values, 3);
std_global_value = mean(std_global_values, 3);
esp_value = mean(esp_values, 3);

edge_matrix_error = abs(edge_matrix_value - mean(edge_matrix_std_av));


G_range=linspace(0.0125, 0.015, 3);%linspace(0.005, 0.015, 9);
beta_range = linspace(1, 1.2, 9);%linspace(1, 1.2, 9); % loop over different beta values 1-3 and optimise

% figure(6)
% 
% cdata = Rspatime_value;
% xvalues = G_range;
% yvalues = beta_range;
% h = heatmap(xvalues,yvalues,cdata);
% 
% h.caxis([0, 0.3]);
% h.Title = 'Std of LKOP (simulation)';
% h.XLabel = 'G';
% h.YLabel = '\beta';
% 
% 
% figure(7)
% 
% 
% cdata = mean_global_value;
% xvalues = G_range;
% yvalues = beta_range;
% h = heatmap(xvalues,yvalues,cdata);
% 
% h.caxis([0, 1]);
% h.Title = 'Mean of GKOP (simulation)';
% h.XLabel = 'G';
% h.YLabel = '\beta';
% 
% 
% 
% figure(8)
% 
% 
% cdata = std_global_value;
% xvalues = G_range;
% yvalues = beta_range;
% h = heatmap(xvalues,yvalues,cdata);
% 
% h.caxis([0, 0.3]);
% h.Title = 'Std of GKOP (simulation)';
% h.XLabel = 'G';
% h.YLabel = '\beta';


figure(3)

cdata = edge_matrix_value;
xvalues = G_range;
yvalues = beta_range;
h = heatmap(xvalues,yvalues,cdata);

%h.caxis([0, 0.3]);
h.Title = 'Std of Edge Measure Matrix (simulation)';
h.XLabel = 'G';
h.YLabel = '\beta';

figure(4)

cdata = esp_value;
xvalues = G_range;
yvalues = beta_range;
h = heatmap(xvalues,yvalues,cdata);

% h.caxis([0, 0.3]);
h.Title = 'ESP (simulation)';
h.XLabel = 'G';
h.YLabel = '\beta';


figure(5)

cdata = edge_matrix_error;
xvalues = G_range;
yvalues = beta_range;
h = heatmap(xvalues,yvalues,cdata);

%h.caxis([0, 0.2]);
h.Title = 'Error in Std of Edge Measure Matrix';
h.XLabel = 'G';
h.YLabel = '\beta';