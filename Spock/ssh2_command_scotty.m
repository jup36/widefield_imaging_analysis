function out = ssh2_command_scotty(varargin)
% ssh2_command_scotty  – lightweight wrapper around system ssh for Scotty
%
% Usage:
%   s_conn = ssh2_command_scotty('connect','scotty');
%   s_conn = ssh2_command_scotty('connect','jp3025@scotty.pni.princeton.edu');
%   out    = ssh2_command_scotty(s_conn, 'hostname');
%
% Optional:
%   s_conn = ssh2_command_scotty('connect','scotty', keyFile);
%
% J.Park / revised Windows-friendly version

% -----------------------------
% Defaults
% -----------------------------
defaultUser = 'jp3025';
defaultHost = 'scotty.pni.princeton.edu';

sshExe = i_get_ssh_exe();

% ============================================================
% Connect mode
% ============================================================
if nargin > 0 && ischar(varargin{1}) && strcmpi(varargin{1}, 'connect')

    if nargin < 2 || isempty(varargin{2})
        error('ssh2_command_scotty:badInput', ...
              'Usage: ssh2_command_scotty(''connect'', userhost, optionalKeyFile)');
    end

    userhost = varargin{2};

    % Resolve short names into explicit user@host.
    userhost = i_resolve_userhost(userhost, defaultUser, defaultHost);

    % Optional key file
    if nargin >= 3 && ~isempty(varargin{3})
        keyFile = i_expand_user_path(varargin{3});
        keyOpt  = sprintf('-i %s', i_quote_path(keyFile));
    else
        keyOpt = '';
    end

    % SSH options.
    % accept-new helps first-time Windows connection avoid hanging at:
    % "Are you sure you want to continue connecting?"
    sshOpts = [ ...
        '-o ConnectTimeout=15 ', ...
        '-o ServerAliveInterval=30 ', ...
        '-o ServerAliveCountMax=2 ', ...
        '-o StrictHostKeyChecking=accept-new ' ...
        ];

    fprintf('Connecting to %s ...\n', userhost);

    cmd = sprintf('%s %s %s %s "hostname"', ...
                  sshExe, sshOpts, keyOpt, userhost);

    fprintf('\nRunning:\n%s\n\n', cmd);

    % Use -echo so password/Duo/passphrase prompts are visible.
    status = system(cmd, '-echo');

    if status ~= 0
        fprintf('\nSSH command failed with status %d.\n', status);
        fprintf('Command was:\n%s\n\n', cmd);

        error('ssh2_command_scotty:connectFail', ...
              ['Could not log in to %s.\n' ...
               'Check VPN, username, key/passphrase, Duo prompt, and host-key confirmation.'], ...
               userhost);
    end

    s_conn.username = regexp(userhost, '^[^@]+', 'match', 'once');
    s_conn.userhost = userhost;
    s_conn.keyOpt   = keyOpt;
    s_conn.sshExe   = sshExe;
    s_conn.sshOpts  = sshOpts;
    s_conn.ready    = true;

    out = s_conn;
    return
end

% ============================================================
% Normal command mode
% ============================================================
if nargin ~= 2
    error('ssh2_command_scotty:badInput', ...
          'Usage: ssh2_command_scotty(s_conn, remoteString)');
end

s_conn    = varargin{1};
remoteCmd = varargin{2};

if ~isstruct(s_conn) || ~isfield(s_conn, 'ready') || ~s_conn.ready
    error('ssh2_command_scotty:badConnection', ...
          'First input must be a valid s_conn returned by connect mode.');
end

if isfield(s_conn, 'sshExe') && ~isempty(s_conn.sshExe)
    sshExe = s_conn.sshExe;
else
    sshExe = i_get_ssh_exe();
end

if isfield(s_conn, 'sshOpts') && ~isempty(s_conn.sshOpts)
    sshOpts = s_conn.sshOpts;
else
    sshOpts = '-o ConnectTimeout=15 -o ServerAliveInterval=30 -o ServerAliveCountMax=2 ';
end

if isfield(s_conn, 'keyOpt') && ~isempty(s_conn.keyOpt)
    keyOpt = s_conn.keyOpt;
else
    keyOpt = '';
end

userhost = s_conn.userhost;

% Escape double quotes in remote command
remoteCmdEsc = strrep(remoteCmd, '"', '\"');

sshCmd = sprintf('%s %s %s %s "%s"', ...
                 sshExe, sshOpts, keyOpt, userhost, remoteCmdEsc);

[status, raw] = system(sshCmd);

out.status = status;
out.raw = raw;
out.sshCmd = sshCmd;
out.command_result = splitlines(string(strtrim(raw)));

end

% ============================================================
% Helper functions
% ============================================================

function sshExe = i_get_ssh_exe()
% Return platform-specific ssh executable.

if ispc
    winSSH = fullfile(getenv('WINDIR'), 'System32', 'OpenSSH', 'ssh.exe');

    if exist(winSSH, 'file')
        sshExe = i_quote_path(winSSH);
    else
        sshExe = 'ssh';
    end
else
    sshExe = 'ssh';
end

end

function userhost = i_resolve_userhost(userhost, defaultUser, defaultHost)
% Resolve short Scotty aliases into explicit user@host.

userhost = char(userhost);

if strcmpi(userhost, 'scotty')
    userhost = sprintf('%s@%s', defaultUser, defaultHost);
    return
end

if strcmpi(userhost, 'scotty.pni.princeton.edu')
    userhost = sprintf('%s@%s', defaultUser, defaultHost);
    return
end

% If user already supplied user@host, keep it.
if contains(userhost, '@')
    return
end

% If they supplied some other host without user, attach default user.
userhost = sprintf('%s@%s', defaultUser, userhost);

end

function p = i_expand_user_path(p)
% Expand ~ to user home directory.

p = char(p);

if isempty(p)
    return
end

if startsWith(p, '~')
    if ispc
        homeDir = getenv('USERPROFILE');
    else
        homeDir = getenv('HOME');
    end

    if strcmp(p, '~')
        p = homeDir;
    elseif startsWith(p, '~/') || startsWith(p, '~\')
        p = fullfile(homeDir, p(3:end));
    end
end

end

function q = i_quote_path(p)
% Quote path safely.

p = char(p);

if startsWith(p, '"') && endsWith(p, '"')
    q = p;
else
    q = ['"', p, '"'];
end

end