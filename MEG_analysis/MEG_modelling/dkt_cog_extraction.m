function cog = dkt_cog_extraction(cog_file, n_parcels)
% 	COG_EXTRACTION  extract cogs from a Thomas Yeo Lab text file.
%   cog = cog_extraction(filename, number of parcels) extracts n_parcels x 3 matrix of cogs.

    input_table = readtable(cog_file);
    
    cog = table2array(input_table);
    
    filename = strcat('cog_dkt', num2str(n_parcels), '.mat');
    
    save(filename, 'cog');
    
end

