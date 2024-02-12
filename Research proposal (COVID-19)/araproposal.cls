% Version 1. Zhixuan Yang (zhixuan.s.yang@gmail.com). 12 Nov 2023.
%
% This work may be distributed and/or modified under the
% conditions of the LaTeX Project Public License, either version 1.3
% of this license or (at your option) any later version.
% The latest version of this license is in
%   https://www.latex-project.org/lppl.txt
% and version 1.3c or later is part of all distributions of LaTeX
% version 2008 or later.
%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{araproposal}[2023/11/06 Unofficial class for ARA proposal]

\DeclareOption{titlesec}{\newcommand{\@loadtitlesec}{}}
\DeclareOption{noparskip}{\newcommand{\@noparskip}{}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax

\LoadClass[11pt]{article}

\RequirePackage[letterpaper, margin=1in]{geometry}
\RequirePackage{charter}
\RequirePackage[numbib]{tocbibind}
\RequirePackage{hyperref}
\RequirePackage{amsmath}
\RequirePackage{amssymb}
\RequirePackage{microtype}
\RequirePackage{enumitem}
\RequirePackage{amsthm}

\ifcsname @loadtitlesec\endcsname%
\RequirePackage{titlesec}
\titlespacing\section{0pt}{10pt plus 5pt minus 5pt}{6pt plus 3pt minus 3pt}
\titlespacing\subsection{0pt}{10pt plus 5pt minus 5pt}{6pt plus 3pt minus 3pt}
\fi

\ifcsname @noparskip\endcsname%
\else
\RequirePackage[parfill]{parskip}
\fi


\newcommand{\piauthor}[1]{\gdef\@piauthor{#1}}%
\newcommand{\@piauthor}{\@latex@warning@no@line{No \noexpand\piauthor given}}

%TODO: I don't know how to implement a _list_ of things in LaTeX so at the moment
% there are just two hard-coded slots for co-PIs.
\newcommand{\copiauthorA}[1]{\gdef\@copiauthorA{#1}}%
\newcommand{\@copiauthorA}{}
\copiauthorA{}
\newcommand{\copiauthorB}[1]{\gdef\@copiauthorB{#1}}%
\newcommand{\@copiauthorB}{}
\copiauthorB{}

\newcommand{\cash}[1]{\gdef\@cash{#1}}%
\newcommand{\@cash}{}
\cash{}

\newcommand{\awscredits}[1]{\gdef\@awscredits{#1}}%
\newcommand{\@awscredits}{}
\awscredits{}

\newcommand{\contact}[1]{\gdef\@contact{#1}}%
\newcommand{\@contact}{}
\contact{}

\renewcommand\@maketitle{%
  \let \footnote \thanks
    {\noindent\LARGE\bfseries \@title \par}%
    {\large
      \vspace{.5em}%
      \begin{flushleft}%
        \textbf{PI}: \@piauthor\\%
        \ifx\@copiauthorA\empty\else\textbf{Co-PI}: \@copiauthorA\\\fi%
        \ifx\@copiauthorB\empty\else\textbf{Co-PI}: \@copiauthorB\\\fi%
        \ifx\@cash\empty\else\textbf{Cash funding needed}: \@cash\\\fi%
        \ifx\@awscredits\empty\else\textbf{Promotional Credits needed}: \@awscredits\\\fi%
        \ifx\@contact\empty\else\textbf{Contact}: \@contact\\\fi%
        \ifx\@date\empty\else\textbf{Date}: \@date\fi%
      \end{flushleft}\par}%
  \par}

\renewcommand\abstract[1]{\noindent\textbf{Abstract} #1}

\newcommand\keywords[1]{\noindent\textbf{Keywords} #1}

\newcommand\email[1]{\href{mailto:#1}{#1}}

\theoremstyle{plain}
\newtheorem{question}{Q}
\newenvironment{answer}{\noindent}{}

