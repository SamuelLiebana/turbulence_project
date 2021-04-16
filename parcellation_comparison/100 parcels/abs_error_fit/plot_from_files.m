% script to collate the abs error mat files from server

edge_matrix_error = zeros(11, 11);
mean_global_error = zeros(11, 11);
Rspatime_error = zeros(11, 11);
std_global_error = zeros(11, 11);
esp_error = zeros(11, 11);

for i = 1:11
    load(sprintf('abs_errors_supercrit_100_beta_%03d.mat',i)) % copy from other file
    edge_matrix_error(i, :) = edge_matrix_std_diff';
    mean_global_error(i, :) = mean_global_diff';
    Rspatime_error(i, :) = Rspatime_diff';
    std_global_error(i, :) = std_global_diff';
    esp_error(i, :) = EdgeSpaTimePredictability_diff';
end

G_range=0:0.05:0.5;
beta_range = 1:0.2:3; % loop over different beta values 1-3 and optimise

figure(1)

cdata = Rspatime_error;
xvalues = G_range;
yvalues = beta_range;
h = heatmap(xvalues,yvalues,cdata);

h.caxis([0, 0.2]);
h.Title = 'Error in Std of LKOP';
h.XLabel = 'G';
h.YLabel = '\beta';


figure(2)


cdata = mean_global_error;
xvalues = G_range;
yvalues = beta_range;
h = heatmap(xvalues,yvalues,cdata);

h.caxis([0, 0.7]);
h.Title = 'Error in Mean of GKOP';
h.XLabel = 'G';
h.YLabel = '\beta';



figure(3)


cdata = std_global_error;
xvalues = G_range;
yvalues = beta_range;
h = heatmap(xvalues,yvalues,cdata);

h.caxis([0, 0.2]);
h.Title = 'Error in Std of GKOP';
h.XLabel = 'G';
h.YLabel = '\beta';


figure(4)

cdata = edge_matrix_error;
xvalues = G_range;
yvalues = beta_range;
h = heatmap(xvalues,yvalues,cdata);

h.caxis([0, 0.2]);
h.Title = 'Error in Std of Edge Measure Matrix';
h.XLabel = 'G';
h.YLabel = '\beta';

figure(5)

cdata = esp_error;
xvalues = G_range;
yvalues = beta_range;
h = heatmap(xvalues,yvalues,cdata);

h.Title = 'Error in ESP';
h.XLabel = 'G';
h.YLabel = '\beta';
