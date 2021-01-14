function cog = cog_extraction(cog_file, n_parcels)

    input_table = readtable(cog_file)
    
    cog = table2array(input_table(:,end-2:end))
    
    filename = strcat('cog_schaefer', num2str(n_parcels), '.mat')
    
    save(filename, 'cog') 
    
end

