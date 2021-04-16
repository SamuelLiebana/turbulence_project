% THIS IS a model for Configuring Jobs that are submitted from your Matlab
% the settings are reset by running the command "ConfigCluster", after
% which you will receive a prompt to input your CSC user name
% you may want to do this before starting new runs of analyses to empty the job history 
% Please always consider the parameters below carefully
% INSTRUCTIONS ONLINE: https://docs.csc.fi/#apps/matlab/
%
% TIP: you can connect (ssh) to Puhti to see your home directory and also
% to monitor the queue online (this seems to update slowly to Matlab)
% after you have the ssh ready type "squeue -l -u $USER" to see the queue
% to see the available licenses (for workers) type the following
% "scontrol show lic=mdcs"
% 
% Jetro J. Tuulari 01 / 2020, jetro.tuulari@utu.fi

% get a handle on the cluster
c = parcluster
% days-hh:mm:s, e.g. 2h = 02:00:0 | 2d12h = 2-12:00:0
c.AdditionalProperties.WallTime = '02:00:0'
% memory usage, this will per node for parallel / multiple cpu runs
c.AdditionalProperties.MemUsage = '2g'
% partition on the cluster, see: https://docs.csc.fi/#computing/running/batch-job-partitions/
c.AdditionalProperties.QueueName = 'small' 
% this is the account name for HPC Big Data Neuroscience Platform (international)
c.AdditionalProperties.AccountName = 'project_2001640'
% email to receive notification to, default is begin and end i.e. "ALL"
c.AdditionalProperties.EmailAddress = 'jetijeitii@gmail.com'
% Check configured values
c.AdditionalProperties
% save profile for the submitted scripts
c.saveProfile
% to check currently running and finished jobs
jobs = c.Jobs