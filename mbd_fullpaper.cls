%----------------------------------------------------------------------------------------------------
%
% Class file for the extended abstract of the 
% ECCOMAS on Multibody Dynamics
%
%----------------------------------------------------------------------------------------------------

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{fullpaper_MBD_2017_LaTeX}[2017/03/28 MBD fullpaper class]  
\ProcessOptions\relax

\DeclareOption{10pt}{\OptionNotUsed}  
\DeclareOption{11pt}{\OptionNotUsed}  
\DeclareOption{twocolumn}{\OptionNotUsed}  
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}  
  
\ProcessOptions\relax  
  
\LoadClass[10pt,a4paper]{article}  
  
\RequirePackage{ucs} %Extended UTF-8 input encoding support for LaTeX - CTAN
\RequirePackage[utf8x]{inputenc}
\RequirePackage{mathptmx} % Use Times as default text font, and provide maths support. This package defines Adobe Times Roman (or equivalent) as default text font
\RequirePackage[T1]{fontenc}

\RequirePackage{graphicx}
\RequirePackage{nth}
\RequirePackage{verbatim}
\RequirePackage{tabularx}
\RequirePackage{amsmath,amsfonts,amssymb,amsthm,amsbsy,amsmath}
\RequirePackage{latexsym}
\RequirePackage{bm}
\RequirePackage{a4}
\RequirePackage{fancyhdr}
\RequirePackage{engord}
\RequirePackage{epsfig,vmargin,booktabs}


\setpapersize{A4}  % needs vmargin package
%\setmarginsrb{25mm}{20mm}{25mm}{20mm}{15mm}{7.6mm}{0mm}{12.5mm}
\setmarginsrb{22mm}{5mm}{22mm}{20mm}{15mm}{0mm}{0mm}{12.5mm}
%\setmarginsrb{leftmargin}{topmargin}{rightmargin}{bottommargin}{headheight}{headsep}{footheight}{footskip}
\parindent=0em 
\parskip=6pt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SECTIONS 
%\renewcommand\section{\@startsection {section}{1}{\z@}{-3.5ex plus -1ex minus-.2ex}{0.2ex}{\normalsize\bf}}
\def\@seccntDot{.}
\def\@seccntformat#1{\csname the#1\endcsname\@seccntDot\hskip 0.3em}
\renewcommand\section{\@startsection {section}{1}{\z@}{-1.0ex plus -1ex minus-.2ex}{0.1ex}{\hskip8mm\large\bf}}
\renewcommand\subsection{\@startsection {subsection}{1}{\z@}{-1.0ex plus -1ex minus-.2ex}{0.1ex}{\hskip8mm\large\bf}}
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KEYWORDS 
\newcommand{\keywords}{\medskip\noindent \textbf{Keywords: }}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcommand{\ts}{\textsuperscript}
%%%%%%%%%%%%%%%%%%%%

\AtBeginDocument{%    
  \pagestyle{empty}
  \thispagestyle{fancy}
  \fancyhf{} % clear all header and footer fields
  \renewcommand{\headrulewidth}{0pt}
  \rhead{
    \scriptsize{
      ECCOMAS Thematic Conference on Multibody Dynamics\\
      June 19-22, 2017, Prague, Czech Republic\\
      \vspace{2mm}}
	}  
}  
  
\endinput


%%% Additional settings like font sizes etc.
\def\figurename{Figure}
\def\tablename{Table}

%Command             10pt    11pt    12pt
%\tiny               5       6       6
%\scriptsize         7       8       8
%\footnotesize       8       9       10
%\small              9       10      10.95
%\normalsize         10      10.95   12
%\large              12      12      14.4
%\Large              14.4    14.4    17.28
%\LARGE              17.28   17.28   20.74
%\huge               20.74   20.74   24.88
%\Huge               24.88   24.88   24.88



