function mdfile = latex2markdown(filename,options)
% Convert Live Script LaTeX output to Markdown
%
% mdfile = npxutils.internal.livescript2markdown.latex2markdown(filename, options)
%
% Converts a LaTeX-published Matlab LiveScript file to Markdown. Produces a .md
% file, and if there are any images, those are saved in a *_images folder next to
% the output .md file.
%
% Filename is the path to a LaTeX file produced by publishing a Matlab Live Script
% to LaTeX format. It may be either an .mlx or a .tex file. If it is .mlx, it is assumed
% that the derived .tex file is sitting right next to it.
%
% Options:
%   outputfilename
%   format - 'quiita' or 'github'*
%   png2jpeg - (logical, false*)
%   tableMaxWidth - (numeric, 20*)
%
% By default, outputs to a .md file next to the input file.

%  Copyright 2020 The MathWorks, Inc.

arguments
    filename (1,1) string
    options.outputfilename char = filename
    options.format char {mustBeMember(options.format,{'qiita','github'})} = 'github'
    options.png2jpeg logical = false
    options.tableMaxWidth (1,1) double = 20
end

% Latex filename
[filepath,name,ext] = fileparts(filename);

if filepath == ""
    filepath = pwd;
end
    
if ext == "" % if without extention, add .tex
    latexfile = fullfile(filepath, name + ".tex");
else %
    latexfile = filename;
end

% check if latexfile exists
assert(exist(latexfile,'file')>0, ...
    latexfile + " does not exist. " ...
    + "If you haven't generated latex file from a live script please do so.");

% Import latex file and have it as a string variable
% If the mlx contains double byte charactors, it needs to use fgets.
% Otherwise the following does the work.
% str = string(fileread(latexfile));
fid = fopen(latexfile,'r','n','UTF-8');
str = string;
tmp = fgets(fid);
while tmp > 0
    str = append(str,string(tmp));
    tmp = fgets(fid);
end
fclose(fid);

% Extract body from latex
str = extractBetween(str,"\begin{document}","\end{document}");

% Delete table of contents: 目次は現時点で削除（TODO）
% ex: \label{H_D152BAC0}
str = regexprep(str,"\\matlabtableofcontents{([^{}]+)}", "");
str = regexprep(str,"\\label{[a-zA-Z_0-9]+}","");

%% Devide the body into each environment
% Preprocess 1:
% Add 'newline' to the end of the following.
% \end{lstlisting}, \end{verbatim}, \end{matlabcode}, \end{matlaboutput},\end{center}
% \end{matlabtableoutput}, \end{matlabsymbolicoutput}  \vspace{1em}
% 要素ごとに分割しやすいように \end の後に改行が無い場合は１つ追加
str = replace(str,"\end{lstlisting}"+newline,"\end{lstlisting}"+newline+newline);
str = replace(str,"\end{verbatim}"+newline,"\end{verbatim}"+newline+newline);
str = replace(str,"\end{matlabcode}"+newline,"\end{matlabcode}"+newline+newline);
str = replace(str,"\end{matlaboutput}"+newline,"\end{matlaboutput}"+newline+newline);
str = replace(str,"\end{matlabtableoutput}"+newline,"\end{matlabtableoutput}"+newline+newline);
str = replace(str,"\end{matlabsymbolicoutput}"+newline,"\end{matlabsymbolicoutput}"+newline+newline);
str = replace(str,"\end{center}"+newline,"\end{center}"+newline+newline);
str = replace(str,"\vspace{1em}"+newline,"\vspace{1em}"+newline+newline);

% Preprocess 2:
% Replace more than three \n to \n\n.
% 3行以上の空白は2行にしておく
str = regexprep(str,'\n{3,}','\n\n');
% Devide them into parts by '\n\n'
% 空白行で要素に分割
str = strsplit(str,'\n\n')';

% Preprocess 3:
% The following environments will be merge into one string
% for the ease of process.
% MATLABコードが複数行にわたるとうまく処理できないので 対応する \end が見つかるまで連結処理
% \begin{verbatim}
% \begin{lstlisting}
% \begin{matlabcode}
% \begin{matlaboutput}
% \begin{matlabtableoutput}
% \begin{matlabsymbolicoutput}
str = mergeSameEnvironments(str,"lstlisting");
str = mergeSameEnvironments(str,"verbatim");
str = mergeSameEnvironments(str,"matlabcode");
str = mergeSameEnvironments(str,"matlaboutput");
str = mergeSameEnvironments(str,"matlabtableoutput");
str = mergeSameEnvironments(str,"matlabsymbolicoutput");

%% Let's convert latex to markdown
% 1: Process parts that require literal output.
[str, idxLiteral] = processLiteralOutput(str);

% 2: Process that other parts
str2md = str(~idxLiteral);
str2md = processDocumentOutput(str2md,options.tableMaxWidth);

% Equations (数式部分)
str2md = processEquations(str2md, options.format);

% includegraphics (画像部分)
str2md = processincludegraphics(str2md, options.format, options.png2jpeg, name, filepath);

% Apply vertical/horizontal space
% markdown: two spaces for linebreak
% latex: \vspace{1em}
% latex: \hskip1em
str2md = regexprep(str2md,"\\vspace{1em}","  ");
str2md = regexprep(str2md,"\\hskip1em","  ");
str(~idxLiteral) = str2md;

%% Done! Merge them together
strmarkdown = join(str,newline);

%% File outputファイル出力
mdfile = options.outputfilename + ".md";
fileID = fopen(mdfile,'w');
fprintf(fileID,'%s\n',strmarkdown);
fclose(fileID);
