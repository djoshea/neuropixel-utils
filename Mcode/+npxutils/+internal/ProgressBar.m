classdef ProgressBar < handle
    % Stripped down version for Neuropixel which does not support parallel mode
    % Generates a progress bar which uses ANSI color codes to colorize output
    % in a color terminal. The background color steadily advances along from
    % left to right beneath a custom message and progress percentage string.
    %
    % Usage:
    %   pbar = ProgressBar('Message goes here', nThingsToProcess);
    %   for i = 1:nThingsToProcess
    %       pbar.update(i, [optional message update]);
    %       ...
    %   end
    %   pbar.finish([optional final message]);
    %
    % Daniel J. O'Shea, 2019
    
    properties(SetAccess=protected)
        message
        n = 0; % last value
        N
        cols
        trueColor = false;
        firstUpdate
        timeStart
        
        lastCalled
        
        usingTerminal
        
        fnamePrefix
        
        %         objWorker
        
        nCompleteByWorker
        
        lastNBoxes = 0;
        lastNSpaces = 0;
        
        trueCmap;
        
        textprogressStrCR = '';
    end
    
    properties(Constant)
        minInterval = 0.1; % seconds
    end
    
    methods
        function pbar = ProgressBar(N, message, varargin)
            if nargin >= 2
                pbar.message = sprintf(message, varargin{:});
            else
                pbar.message = '';
            end
            
            if nargin >= 1
                pbar.N = N;
            else
                pbar.N = 1;
            end
            
            pbar.usingTerminal = ~usejava('desktop');
            
            [~, pbar.cols] = npxutils.internal.ProgressBar.getTerminalSize();
            pbar.trueColor = ~isempty(getenv('ITERM_PROFILE')) && false;
            
            if pbar.trueColor
                hsv = ones(pbar.cols, 3);
                hsv(:, 1) = 0.5;
                hsv(:, 2) = 0.6;
                b = 0.5 * (1+sin((1:pbar.cols) / 8));
                x = 0.3;
                b = b*x + 0.95-x;
                hsv(:, 3) = b;
                pbar.trueCmap = round(256*hsv2rgb(hsv));
            end
            
            pbar.firstUpdate = true;
            pbar.timeStart = clock;
            pbar.lastNBoxes = 0;
            pbar.lastNSpaces = 0;
            pbar.update(0);
        end
        
        function increment(this, varargin)
            if nargin > 1 && isnumeric(varargin{1})
                by = varargin{1};
                varargin = varargin(2:end);
            else
                by = 1;
            end
            this.update(this.n+by, varargin{:});
        end
        
        function update(this, n, message, varargin)
            if nargin > 2
                this.message = sprintf(message, varargin{:});
            else
            end
            this.n = n; % store last update
            
            % don't run too often
            if isempty(this.lastCalled)
                this.lastCalled = clock;
            elseif etime(clock, this.lastCalled) < this.minInterval
                return;
            end
            
            this.lastCalled = clock;
            
            if this.N > 0
                numWidth = ceil(log10(this.N));
            else
                numWidth = 1;
            end
            
            if n < 0
                n = 0;
            end
            if isempty(this.N) || this.N == 1
                ratio = n;
                if ratio < 0
                    ratio = 0;
                end
                if ratio > 1
                    ratio = 1;
                end
                percentage = ratio * 100;
                progStr = sprintf('[ %5.1f%% ]', percentage);
            else
                ratio = (n-1)/this.N;
                percentage = min(max(ratio*100, 0), 100);
                progStr = sprintf('%*d / %*d [ %5.1f%% ]', numWidth, n, numWidth, this.N, percentage);
            end
            
            progLen = length(progStr);
            
            if ratio < 0
                ratio = 0;
            end
            if ratio > 1
                ratio = 1;
            end
            
            if length(this.message) + progLen + 3 > this.cols
                message = [this.message(1:(this.cols - progLen - 6)), '...'];
            else
                message = this.message;
            end
            
            gap = this.cols - 1 - (length(message)+1) - progLen;
            spaces = repmat(' ', 1, gap);
            if this.usingTerminal
                str = [message spaces progStr];
            else
                str = [message spaces blanks(numel(progStr))];
            end
            
            % separate into colored portion of bar and non-colored portion of bar
            ind = min(length(str), ceil(ratio*this.cols));
            preStr = str(1:ind);
            postStr = str(ind+1:end);
            
            % try using 24 color
            if this.trueColor
                newPreStr = '';
                for i = 1:numel(preStr)
                    color = this.trueCmap(i, :);
                    %                     color = pbar.trueCmap(mod(i-1, size(pbar.trueCmap, 1))+1, :);
                    newPreStr = [newPreStr, sprintf('\x1b[48;2;%d;%d;%dm%s' , color, preStr(i))];
                end
                preStr = newPreStr;
            end
            
            if this.usingTerminal
                if this.firstUpdate
                    fprintf(' '); % don't delete whole line on first update
                end
                if this.trueColor
                    fprintf('\b\r%s\033[49;37m%s\033[0m ', preStr, postStr);
                else
                    fprintf('\b\r\033[1;44;37m%s\033[49;37m%s\033[0m  ', preStr, postStr);
                end
            else
                this.textprogressbar(ratio);
            end
            this.firstUpdate = false;
            
            %str = sprintf('\b\r\033[1;44;37m %s\033[49;37m%s\033[0m ', preStr, postStr);
            %disp(str);
            
        end
        
        function finish(this, message, varargin)
            % if message is provided (also in printf format), the message
            % will be displayed. Otherwise, the progress bar will disappear
            % and output will resume on the same line.
            
            if this.usingTerminal
                fprintf('\033[1Acll\033[2K\r');
            else
                this.textprogressbar(1);
                fprintf('\n');
            end
            if nargin > 1
                this.message = sprintf(message, varargin{:});
                fprintf('%s\n', this.message);
            end
            
        end
        
        function textprogressbar(this, c)
            %
            % Original Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
            % Version: 1.0
            % Changes tracker:  29.06.2010  - First version
            %
            % Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/
            %
            % Modified by Dan O'Shea
            
            %% Initialization
            
            strPercentageLength = 9;   %   Length of percentage string (must be >5)
            strDotsMaximum      = 15;   %   The total number of dots in a progress bar
            
            if this.firstUpdate
                fprintf('%s : ', this.message);
                this.textprogressStrCR = -1;
            end
            
            c = round(c*100, 1);
            
            percentageOut = [num2str(c) '%%'];
            percentageOut = [repmat(' ',1,strPercentageLength-length(percentageOut)-1) percentageOut ' '];
            nDots = floor(c/100*strDotsMaximum);
            dotOut = ['[' repmat('_',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']   '];
            strOut = [percentageOut dotOut];
            
            % Print it on the screen
            if this.textprogressStrCR == -1
                % Don't do carriage return during first run
                fprintf(strOut);
            else
                % Do it during all the other runs
                fprintf([this.textprogressStrCR strOut]);
            end
            
            % Update carriage return
            this.textprogressStrCR = repmat('\b',1,length(strOut)-1);
            
        end
    end
    
    methods (Static)
        function [rows, cols] = getTerminalSize()
            usingTerminal = ~usejava('desktop');
            
            % use sensible defaults
            rows = 24;
            cols = 80;
            
            if (ismac || isunix) && usingTerminal
                % actual terminal: get terminal width using tput
                cmd = 'tput lines';
                [~, r] = unix(cmd);
                num = sscanf(r, '%d');
                if ~isempty(num)
                    rows = num(end);
                end
                
                cmd = 'tput cols';
                [~, r] = unix(cmd);
                num = sscanf(r, '%d');
                if ~isempty(num)
                    cols = num(end);
                end
                
            elseif ~usingTerminal %#ok<*PROP
                % matlab command window size
                try
                    jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
                    cmdWin = jDesktop.getClient('Command Window');
                    
                    jTextArea = cmdWin.getComponent(0).getViewport.getComponent(0);
                    height = jTextArea.getHeight();
                    width = jTextArea.getParent.getWidth();
                    font = jTextArea.getParent().getFont();
                    metrics = cmdWin.getFontMetrics(font);
                    charWidth = metrics.charWidth('A');
                    charHeight = metrics.getHeight();
                    
                    rows = floor(height/charHeight);
                    cols = floor(width/charWidth);
                catch
                end
            end
        end
        
    end
end
