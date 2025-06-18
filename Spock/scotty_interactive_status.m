% ---- 1. Quick status line (RUNNING / PENDING / etc.) ------------
rez.squeue = ssh2_command_scotty(s_conn, 'squeue -j 1114359 --Format=jobid,state,timeleft,nodelist');

% ---- 2. Full job description (node, memory, submit dir, etc.) ---
rez.scontrol = ssh2_command_scotty(s_conn, 'scontrol show job 1114359');

% ---- 3. Historical record after it finishes ---------------------
rez.sacct = ssh2_command_scotty(s_conn, ...
    'sacct -j 1114359 --format=JobID,State,Elapsed,ExitCode');
