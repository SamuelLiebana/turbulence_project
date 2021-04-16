% script to collate the abs error mat files from server

edge_matrix_value = zeros(11, 11);
mean_global_value = zeros(11, 11);
Rspatime_value = zeros(11, 11);
std_global_value = zeros(11, 11);
esp_value = zeros(11, 11);

for i = 1:11
    load(sprintf('values_supercrit_1000_beta_%03d.mat',i)) % copy from other file
    edge_matrix_value(i, :) = edge_matrix_std_av_sim';
    mean_global_value(i, :) = mean_global_av_sim';
    Rspatime_value(i, :) = Rspatime_av_sim';
    std_global_value(i, :) = std_global_av_sim';
    esp_value(i, :) = EdgeSpaTimePredictability_av_sim';
end

G_range=0:0.05:0.5;
beta_range = 1:0.2:3; % loop over different beta values 1-3 and optimise

figure(6)

cdata = Rspatime_value;
xvalues = G_range;
yvalues = beta_range;
h = heatmap(xvalues,yvalues,cdata);

h.caxis([0, 0.3]);
h.Title = 'Std of LKOP (simulation)';
h.XLabel = 'G';
h.YLabel = '\beta';


figure(7)


cdata = mean_global_value;
xvalues = G_range;
yvalues = beta_range;
h = heatmap(xvalues,yvalues,cdata);

h.caxis([0, 1]);
h.Title = 'Mean of GKOP (simulation)';
h.XLabel = 'G';
h.YLabel = '\beta';



figure(8)


cdata = std_global_value;
xvalues = G_range;
yvalues = beta_range;
h = heatmap(xvalues,yvalues,cdata);

h.caxis([0, 0.3]);
h.Title = 'Std of GKOP (simulation)';
h.XLabel = 'G';
h.YLabel = '\beta';


figure(9)

cdata = edge_matrix_value;
xvalues = G_range;
yvalues = beta_range;
h = heatmap(xvalues,yvalues,cdata);

h.caxis([0, 0.3]);
h.Title = 'Std of Edge Measure Matrix (simulation)';
h.XLabel = 'G';
h.YLabel = '\beta';

figure(10)

cdata = esp_value;
xvalues = G_range;
yvalues = beta_range;
h = heatmap(xvalues,yvalues,cdata);

h.caxis([0, 1]);
h.Title = 'ESP (simulation)';
h.XLabel = 'G';
h.YLabel = '\beta';