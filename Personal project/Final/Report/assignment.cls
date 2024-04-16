% assignment.cls
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{assignment}[2024/01/19 v1.0 ME621 Assignments class]

\LoadClass{article}

% Page layout
\RequirePackage[a4paper, margin=20mm]{geometry}

% Colors and graphics
\RequirePackage[dvipsnames,table]{xcolor}
\RequirePackage{graphicx}
\RequirePackage{subcaption}

% Language settings
\RequirePackage[english]{babel}
\RequirePackage[acronym]{glossaries}

% Hyperlinks
\RequirePackage[hidelinks]{hyperref}

% Mathematical typesetting
\RequirePackage{amsmath}

% Code listings
\RequirePackage{listings}
\RequirePackage{xcolor}

\definecolor{codegreen}{rgb}{0,0.6,0}
\definecolor{codegray}{rgb}{0.5,0.5,0.5}
\definecolor{codepurple}{rgb}{0.58,0,0.82}
\definecolor{backcolour}{rgb}{0.95,0.95,0.92}

\definecolor{dkgreen}{rgb}{0,0.6,0}
\definecolor{gray}{rgb}{0.5,0.5,0.5}
\definecolor{mauve}{rgb}{0.58,0,0.82}
\definecolor{gray}{rgb}{0.4,0.4,0.4}
\definecolor{darkblue}{rgb}{0.0,0.0,0.6}
\definecolor{lightblue}{rgb}{0.0,0.0,0.9}
\definecolor{cyan}{rgb}{0.0,0.6,0.6}
\definecolor{darkred}{rgb}{0.6,0.0,0.0}

\lstset{
    basicstyle        = \sffamily\footnotesize,
    keywordstyle      = \bfseries,
    commentstyle      = \rmfamily\itshape,
    stringstyle       = \ttfamily,
    backgroundcolor   = \color{white},
    numberstyle       = \tiny\color{codegray},
    flexiblecolumns,
    frame             = single,
    xleftmargin       = 15pt,
    xrightmargin      = 5pt,
    breakatwhitespace = false,
    breaklines        = true,
    basewidth         = 0.5em,
    captionpos        = b,
    columns           = fixed,
    keepspaces        = true,
    numbers           = left,
    numbersep         = 5pt,
    showspaces        = false,
    showstringspaces  = false,
    showtabs          = false,
    tabsize           = 4
}

\lstdefinestyle{mathematica}{
    commentstyle = \color[HTML]{338AAF}\ttfamily,
    stringstyle  = \color{codepurple},
}

\lstdefinestyle{C}{
    language     = C,
    keywordstyle = \color{green!40!black},
    commentstyle = \color{purple!40!black},
    stringstyle  = \color{brown},
    escapeinside = {(*@}{@*)},
}


% \lstset{
%   basicstyle=\ttfamily\footnotesize,
%   columns=fullflexible,
%   showstringspaces=false,
%   numbers=left,                   % where to put the line-numbers
%   numberstyle=\tiny\color{gray},  % the style that is used for the line-numbers
%   stepnumber=1,
%   numbersep=5pt,                  % how far the line-numbers are from the code
%   backgroundcolor=\color{white},      % choose the background color. You must add \usepackage{color}
%   showspaces=false,               % show spaces adding particular underscores
%   showstringspaces=false,         % underline spaces within strings
%   showtabs=false,                 % show tabs within strings adding particular underscores
%   frame=none,                   % adds a frame around the code
%   rulecolor=\color{black},        % if not set, the frame-color may be changed on line-breaks within not-black text (e.g. commens (green here))
%   tabsize=2,                      % sets default tabsize to 2 spaces
%   captionpos=b,                   % sets the caption-position to bottom
%   breaklines=true,                % sets automatic line breaking
%   breakatwhitespace=false,        % sets if automatic breaks should only happen at whitespace
%   title=\lstname,                   % show the filename of files included with \lstinputlisting;
%                                   % also try caption instead of title
%   commentstyle=\color{gray}\upshape
% }


\lstdefinestyle{XML}
{
    basicstyle=\ttfamily\footnotesize,
    morestring=[s][\color{mauve}]{"}{"},
    morestring=[s][\color{black}]{>}{<},
    morecomment=[s]{<?}{?>},
    morecomment=[s][\color{dkgreen}]{<!--}{-->},
    stringstyle=\color{black},
    identifierstyle=\color{lightblue},
    keywordstyle=\color{red},
    morekeywords={xmlns,xsi,noNamespaceSchemaLocation,type,id,x,y,source,target,version,tool,transRef,roleRef,objective,eventually}% list your attributes here
}

\endinput
