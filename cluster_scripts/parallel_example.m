function t = parallel_example
t0 = tic;
parfor idx = 1:16
    A(idx) = idx;
    pause(2)
end
t = toc(t0);
