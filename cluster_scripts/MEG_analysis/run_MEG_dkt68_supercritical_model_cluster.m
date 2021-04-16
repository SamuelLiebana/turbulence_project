% get a handle on the cluster
c = parcluster
% days-hh:mm:s, e.g. 2h = 02:00:0 | 2d12h = 2-12:00:0
c.AdditionalProperties.WallTime = '2-23:00:0'
% memory usage, this will per node for parallel / multiple cpu runs
c.AdditionalProperties.MemUsage = '2g'
% partition on the cluster, see: https://docs.csc.fi/#computing/running/batch-job-partitions/
c.AdditionalProperties.QueueName = 'large' 
% this is the account name for HPC Big Data Neuroscience Platform (international)
c.AdditionalProperties.AccountName = 'project_2001640'
% email to receive notification to, default is begin and end i.e. "ALL"
c.AdditionalProperties.EmailAddress = 'sc.liebanagarcia@gmail.com'
% Check configured values
c.AdditionalProperties
% save profile for the submitted scripts
c.saveProfile

%long range

for i=1:9
    for o = 1:9
        j = c.batch(@dkt68_supercritical_model_cluster, 0, {i, o}, 'pool', 2, 'CurrentFolder', '.', 'AutoAddClientPath',false)
    end
end;