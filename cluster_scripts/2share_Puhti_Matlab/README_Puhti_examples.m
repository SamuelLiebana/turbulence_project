%% ---- Copy paste the following to the Matlab command line prompt stepwise
% This is really just checking that all settings are in check
%
% Jetro J. Tuulari 01 / 2020, jetro.tuulari@utu.fi

%%
% to send a task that prints the working directory (Puhti server)
j = batch(c, @pwd, 1, {}, 'CurrentFolder', '.', 'AutoAddClientPath', false)

% OPTIONAL: wait for the results and fetch the outputs
j.wait

% to check currently running and finished jobs
jobs = c.Jobs

% >> % Get a handle to the job with sequence number 1
j1 = c.Jobs(1)
% Fetch results for job number 1
fetchOutputs(j1)

% The output should read (not surprisingly)
% 
% ans =
% 
%   1×1 cell array
% 
%     {'/users/tuulari'}

%%
% If there is an issue try finding the debug log
j.Parent.getDebugLog(j.Tasks(1))

%%
% TESTING A PARALLEL JOB

% Submitting a parallel job to 8 cores, taking N + 1 cores
% submits a function called "parallel_example.m"
% this should be located in the tutorial folder
j = batch(c, @parallel_example, 1, {}, 'pool', 8, 'CurrentFolder','.', 'AutoAddClientPath',false)

% to check currently running and finished jobs
jobs = c.Jobs

% Once the job is finished, (job number 2)
% >> % Get a handle to the job with sequence number 2
j2 = c.Jobs(2)

% Fetch results for job number 2
fetchOutputs(j2)

% An example output
% 
% ans =
% 
%   1×1 cell array
% 
%     {[5.8952]}

%%
% If there is an issue, get the debug log .. 
j = c.findJob('ID', 2);
j.Parent.getDebugLog(j)

%% REAL LIFE EXAMPLE USING HEDONIA LAB CODE (partial LEiDA)
%
% % Parallel run for partial LEiDA across 3 groups, using 12 workers
% % (K means clustering uses parallel processing in the model script)
% % Key here is that we are reading the V1_all_new.mat file that is stored in the
% % cluster in tuulari's home directory, and the output will be written there as well 
%
% use scp to move the file (rsync would be used for folders)
%
% MODIFY THE CurrentFolder to the place where V1_all_new.mat has been
% copied
j = c.batch(@Kmeans_fromLEiDA_function, 1, {}, 'CurrentFolder', '/users/tuulari/2020-01-Kmeans-tryout/', 'AutoAttachFiles',false, 'AutoAddClientPath', false, 'pool', 12)

%%
% to check currently running and finished jobs
jobs = c.Jobs

%%
% OPTIONAL: check what the matlab print looks like
j.diary

% Once the job is finished, (job number 3)
% >> % Get a handle to the job with ID 3
j3 = c.Jobs(3)

% Fetch results for job number 3
fetchOutputs(j3)
%% END