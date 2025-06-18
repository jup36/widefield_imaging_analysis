%% update_ganymed_ssh2.m
% -------------------------------------------------------------
% Purpose:  Replace any old Ganymed-SSH2 library with build 262,
%           add it permanently via classpath.txt, and verify it’s
%           on the dynamic class path for the current session.
%
% How to use:
%   • Edit `newJar` if you put the JAR somewhere else.
%   • Run the script once in MATLAB.
%   • Restart MATLAB so the static class path is refreshed.
%
% (c) 2025  J.Park — feel free to modify / redistribute
% -------------------------------------------------------------

%% 1.  Path to the newest Ganymed-SSH2 JAR
newJar = '/Users/jp3025/matlab/jars/ganymed-ssh2-262.jar';

%% 2.  Remove any *older* Ganymed JARs that were loaded dynamically
jcp = javaclasspath('-dynamic');
bad = contains(jcp, 'ganymed-ssh2');
if any(bad)
    fprintf('Removing outdated JAR(s):\n');
    disp(jcp(bad)');
    cellfun(@javarmpath, jcp(bad));
else
    fprintf('No older Ganymed JARs found on dynamic path.\n');
end

%% 3.  Add the new JAR for *this* MATLAB session
javaaddpath(newJar);
fprintf('Added new JAR to dynamic path:\n  %s\n', newJar);

%% 4.  Show dynamic class path for confirmation
disp('Current dynamic class path now contains:');
disp(javaclasspath('-dynamic')');

%% 5.  Write/overwrite classpath.txt in prefdir for *future* sessions
cpFile = fullfile(prefdir, 'classpath.txt');
fid = fopen(cpFile, 'w');
fprintf(fid, '%s\n', newJar);
fclose(fid);
fprintf('classpath.txt written to:\n  %s\n', cpFile);

%% 6.  Friendly reminder
fprintf([ ...
    '\nAll done!  Please restart MATLAB so the new JAR is loaded\n', ...
    'as part of the *static* class path. After restart, run your\n', ...
    'ssh2_smoketest to confirm the Scotty handshake.\n']);
