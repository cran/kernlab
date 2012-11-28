\def\fileversion{1.0}
\def\filename{A}
\def\filedate{2004/10/08}
%%
%%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{A}[\filedate\space\fileversion\space A class ]

%% options
\LoadClass[10pt,a4paper,twoside]{article}
\newif\if@notitle
\@notitlefalse
\DeclareOption{notitle}{\@notitletrue}
\ProcessOptions

%% required packages
\RequirePackage{graphicx,a4wide,color,hyperref,ae,fancyvrb,thumbpdf}
\RequirePackage[T1]{fontenc}
\usepackage[authoryear,round,longnamesfirst]{natbib}
\bibpunct{(}{)}{;}{a}{}{,}
\bibliographystyle{jss}

%% paragraphs
\setlength{\parskip}{0.7ex plus0.1ex minus0.1ex}
\setlength{\parindent}{0em}

%% commands
\let\code=\texttt
\let\proglang=\textsf
\newcommand{\pkg}[1]{{\normalfont\fontseries{b}\selectfont #1}}
\newcommand{\email}[1]{\href{mailto:#1}{\normalfont\texttt{#1}}}
\newcommand{\E}{\mathsf{E}}
\newcommand{\VAR}{\mathsf{VAR}}
\newcommand{\COV}{\mathsf{COV}}
\newcommand{\Prob}{\mathsf{P}}

%% for all publications
\newcommand{\Plaintitle}[1]{\def\@Plaintitle{#1}}
\newcommand{\Shorttitle}[1]{\def\@Shorttitle{#1}}
\newcommand{\Plainauthor}[1]{\def\@Plainauthor{#1}}
\newcommand{\Keywords}[1]{\def\@Keywords{#1}}
\newcommand{\Plainkeywords}[1]{\def\@Plainkeywords{#1}}
\newcommand{\Abstract}[1]{\def\@Abstract{#1}}

%% defaults
\author{Firstname Lastname\\Affiliation}
\title{Title}
\Abstract{---!!!---an abstract is required---!!!---}
\Plainauthor{\@author}
\Plaintitle{\@title}
\Shorttitle{\@title}
\Keywords{---!!!---at least one keyword is required---!!!---}
\Plainkeywords{\@Keywords}

%% Sweave(-like)
%\DefineVerbatimEnvironment{Sinput}{Verbatim}{fontshape=sl}
%\DefineVerbatimEnvironment{Soutput}{Verbatim}{}
%\DefineVerbatimEnvironment{Scode}{Verbatim}{fontshape=sl}
%\newenvironment{Schunk}{}{}
\DefineVerbatimEnvironment{Code}{Verbatim}{}
\DefineVerbatimEnvironment{CodeInput}{Verbatim}{fontshape=sl}
\DefineVerbatimEnvironment{CodeOutput}{Verbatim}{}
\newenvironment{CodeChunk}{}{}
\setkeys{Gin}{width=0.8\textwidth}

%% new \maketitle
\def\maketitle{
 \begingroup
   \def\thefootnote{\fnsymbol{footnote}}
   \def\@makefnmark{\hbox to 0pt{$^{\@thefnmark}$\hss}}
   \long\def\@makefntext##1{\parindent 1em\noindent
			    \hbox to1.8em{\hss $\m@th ^{\@thefnmark}$}##1}
   \@maketitle \@thanks
 \endgroup
 \setcounter{footnote}{0}
 \thispagestyle{empty}
 \markboth{\centerline{\@Shorttitle}}{\centerline{\@Plainauthor}}
 \pagestyle{myheadings}

 \let\maketitle\relax \let\@maketitle\relax
 \gdef\@thanks{}\gdef\@author{}\gdef\@title{}\let\thanks\relax
}

\def\@maketitle{\vbox{\hsize\textwidth \linewidth\hsize 
 {\centering
 {\LARGE\bf \@title\par}
   \def\And{\end{tabular}\hfil\linebreak[0]\hfil
	    \begin{tabular}[t]{c}\large\bf\rule{\z@}{24pt}\ignorespaces}%
    \begin{tabular}[t]{c}\large\bf\rule{\z@}{24pt}\@author\end{tabular}%
 \vskip 0.3in minus 0.1in
 \hrule
 \begin{abstract}
 \@Abstract
 \end{abstract}}
 \textit{Keywords}:~\@Keywords.
 \vskip 0.1in minus 0.05in
 \hrule
 \vskip 0.2in minus 0.1in
}}


%% sections, subsections, and subsubsections
\newlength{\preXLskip}
\newlength{\preLskip}
\newlength{\preMskip}
\newlength{\preSskip}
\newlength{\postMskip}
\newlength{\postSskip}
\setlength{\preXLskip}{1.8\baselineskip plus 0.5ex minus 0ex}
\setlength{\preLskip}{1.5\baselineskip plus 0.3ex minus 0ex}
\setlength{\preMskip}{1\baselineskip plus 0.2ex minus 0ex}
\setlength{\preSskip}{.8\baselineskip plus 0.2ex minus 0ex}
\setlength{\postMskip}{.5\baselineskip plus 0ex minus 0.1ex}
\setlength{\postSskip}{.3\baselineskip plus 0ex minus 0.1ex}

\newcommand{\jsssec}[2][default]{\vskip \preXLskip%
  \pdfbookmark[1]{#1}{Section.\thesection.#1}%
  \refstepcounter{section}%
  \centerline{\textbf{\Large \thesection. #2}} \nopagebreak
  \vskip \postMskip \nopagebreak}
\newcommand{\jsssecnn}[1]{\vskip \preXLskip%
  \centerline{\textbf{\Large #1}} \nopagebreak
  \vskip \postMskip \nopagebreak}

\newcommand{\jsssubsec}[2][default]{\vskip \preMskip%
  \pdfbookmark[2]{#1}{Subsection.\thesubsection.#1}%
  \refstepcounter{subsection}%
  \textbf{\large \thesubsection. #2} \nopagebreak
  \vskip \postSskip \nopagebreak}
\newcommand{\jsssubsecnn}[1]{\vskip \preMskip%
  \textbf{\large #1} \nopagebreak
  \vskip \postSskip \nopagebreak}

\newcommand{\jsssubsubsec}[2][default]{\vskip \preSskip%
  \pdfbookmark[3]{#1}{Subsubsection.\thesubsubsection.#1}%
  \refstepcounter{subsubsection}%
  {\large \textit{#2}} \nopagebreak
  \vskip \postSskip \nopagebreak}
\newcommand{\jsssubsubsecnn}[1]{\vskip \preSskip%
  {\textit{\large #1}} \nopagebreak
  \vskip \postSskip \nopagebreak}

\newcommand{\jsssimplesec}[2][default]{\vskip \preLskip%
%%  \pdfbookmark[1]{#1}{Section.\thesection.#1}%
  \refstepcounter{section}%
  \textbf{\large #1} \nopagebreak
  \vskip \postSskip \nopagebreak}
\newcommand{\jsssimplesecnn}[1]{\vskip \preLskip%
  \textbf{\large #1} \nopagebreak
  \vskip \postSskip \nopagebreak}

\renewcommand{\section}{\secdef \jsssec \jsssecnn}
\renewcommand{\subsection}{\secdef \jsssubsec \jsssubsecnn}
\renewcommand{\subsubsection}{\secdef \jsssubsubsec \jsssubsubsecnn}

%% colors
\definecolor{Red}{rgb}{0.7,0,0}
\definecolor{Blue}{rgb}{0,0,0.8}
\hypersetup{%
  hyperindex = {true},
  colorlinks = {true},
  linktocpage = {true},
  plainpages = {false},
  linkcolor = {Blue},
  citecolor = {Blue},
  urlcolor = {Red},
  pdfstartview = {Fit},
  pdfpagemode = {UseOutlines},
  pdfview = {XYZ null null null}
}

\AtBeginDocument{
  \hypersetup{%
    pdfauthor = {\@Plainauthor},
    pdftitle = {\@Plaintitle},
    pdfkeywords = {\@Plainkeywords}
  }
}
\if@notitle
  %% \AtBeginDocument{\maketitle}
\else
  \AtBeginDocument{\maketitle}
\fi
