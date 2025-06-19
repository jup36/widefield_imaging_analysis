function out = ssh2_command_scotty(varargin)
% ssh2_command_scotty  – lightweight wrapper around system ssh/scp for Scotty
%
%  out = ssh2_command_scotty(s_conn, remoteCmd)
%      Executes remoteCmd on Scotty and returns a struct:
%         .status          – system() exit code
%         .command_result  – cell array of each output line
%
%  s_conn = ssh2_command_scotty('connect', userhost, keyFile)
%      Prompts for password (if keyFile omitted) and returns a struct
%      that can be reused for subsequent calls.
%
%  Notes
%  • keyFile is optional and should be a full path to your private key
%    (e.g. '~/.ssh/id_ed25519').  If omitted, ssh will fall back to
%    agent / password / Duo as usual.
%
%  J.Park — 2025-06-18

% ---------------------------------------------------------------------
if nargin > 0 && ischar(varargin{1}) && strcmpi(varargin{1},'connect')
    % interactive setup
    userhost = varargin{2};
    if nargin >= 3
        keyFile = varargin{3};
        keyOpt  = ['-i ', keyFile];
    else
        keyOpt  = '';
    end
    
    % test login with a harmless command
    fprintf('Connecting to %s …\n', userhost);
    cmd = sprintf('ssh %s %s "hostname"', keyOpt, userhost);
    status = system(cmd);
    
    if status ~= 0
        error('ssh2_command_scotty:connectFail', ...
              'Could not log in to %s – check key or Duo.', userhost);
    end
    
    s_conn.username  = regexp(userhost,'^[^@]+','match','once');
    s_conn.userhost  = userhost;
    s_conn.keyOpt    = keyOpt;
    s_conn.ready     = true;
    
    out = s_conn;
    return
end
% ---------------------------------------------------------------------
% normal command mode
if nargin ~= 2
    error('ssh2_command_scotty:badInput', ...
          'Usage: ssh2_command_scotty(s_conn, remoteString)');
end
s_conn     = varargin{1};
remoteCmd  = varargin{2};

sshCmd = sprintf('ssh %s %s "%s"', ...
                 s_conn.keyOpt, s_conn.userhost, remoteCmd);

[status, raw] = system(sshCmd);
out.status = status;
out.command_result = splitlines(string(strtrim(raw)));
end
