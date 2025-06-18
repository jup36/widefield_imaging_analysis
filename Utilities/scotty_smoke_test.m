% ---------------- ssh2_smoketest.m -------------------------------
% Purpose: 1) open an SSH2 session to Scotty
%          2) run a trivial command ('hostname')
%          3) print the result and close the session
%
% Prerequisites:
%   • ganymed-ssh2.jar is on the Java path  ➜  javaaddpath('path/to/ganymed-ssh2.jar')
%   • ssh2_config / ssh2_command / ssh2_close are in the MATLAB path

try
    % --- 1. Gather credentials ---
    user = input('Scotty username (e.g. puID): ', 's');
    pass = passcode();                    % your helper for Duo / password

    % --- 2. Configure simple connection ---
    conn = ssh2_config('scotty.pni.princeton.edu', user, pass);

    % --- 3. Issue a trivial command ---
    [conn, reply] = ssh2_command(conn, 'hostname', 1);  % 1 prints to screen

    % --- 4. Close the session cleanly ---
    conn = ssh2_close(conn);

    fprintf('\n✅ SSH2 smoke test succeeded.  Reply was: %s\n', strjoin(reply));

catch ME
    fprintf('\n❌ SSH2 smoke test failed:\n%s\n', getReport(ME, 'extended'));
    fprintf('• If the error mentions authentication or Duo, try key-pair login.\n');
    fprintf('• If the error mentions “no matching key exchange / cipher,”\n');
    fprintf('  update to the latest ganymed-ssh2.jar (v262) or use JSch.\n');
end
% -----------------------------------------------------------------
