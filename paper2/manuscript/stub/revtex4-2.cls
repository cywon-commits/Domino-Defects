% PREVIEW STUB ONLY - not REVTeX. Mimics the commands used in main.tex so the draft can be compiled
% where REVTeX is not installed. Compile with the real revtex4-2 for the actual manuscript.
\NeedsTeXFormat{LaTeX2e}\ProvidesClass{revtex4-2}[2026/09/29 preview stub]
\DeclareOption*{}\ProcessOptions\relax
\LoadClass[10pt,twocolumn]{article}
\RequirePackage[margin=1.8cm,columnsep=0.7cm]{geometry}\RequirePackage{array}
\newcommand{\affiliation}[1]{\gdef\@affil{#1}}\newcommand{\email}[1]{\gdef\@email{#1}}
\def\@affil{}\def\@email{}
\newsavebox\@absbox
\renewenvironment{abstract}{\global\setbox\@absbox\vbox\bgroup\hsize=0.82\textwidth\small\noindent}{\egroup}
\renewcommand{\maketitle}{\twocolumn[{\centering{\Large\bfseries\@title\par}\vskip 1em{\large\@author\par}{\small\itshape\@affil\par}{\small\ttfamily\@email\par}{\small\@date\par}\vskip 1em\box\@absbox\vskip 1.5em}]}
\newenvironment{ruledtabular}{\begin{center}\small}{\end{center}}
