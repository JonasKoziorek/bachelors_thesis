\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{cwru_thesis}
              [2023/02/10 v0.1 CWRU thesis class]

% Load memoir
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions\relax
\LoadClass[oneside,12pt]{memoir}

\RequirePackage{etoolbox}

% Overleaf 28 Jan: Page geometry
\RequirePackage[letterpaper,top=1in,bottom=1in,left=1.5in,right=1.5in,headsep=\dimexpr 0.25in-1em,marginparwidth=1.25in]{geometry}

% Page style
\makeoddhead{headings}{}{}{\thepage}
\makeevenhead{headings}{}{}{\thepage}

\RequirePackage{graphicx}
\RequirePackage{titlecaps}

%% Overleaf 24 Jan: Formatting of ToC and lists
\setcounter{tocdepth}{2}
\let\listoffigurestables\listoffigures
\appto{\tableofcontents}{\clearpage}
\appto{\listoffigurestables}{\clearpage}
\appto{\listoffigures}{\clearpage}
\appto{\listoftables}{\clearpage}
\renewcommand{\contentsname}{Table of Contents}
\renewcommand{\listfigurename}{List of Figures}
\renewcommand{\listtablename}{List of Tables}


%\addtocontents{toc}{\protect\vskip 35pt}
%\addtocontents{toc}{\protect\SingleSpacing}
%\addtocontents{lof}{\protect\vskip 35pt}
%% In case single spacing is desired later
% \addtocontents{lof}{\protect\SingleSpacing}
% \addtocontents{lot}{\protect\SingleSpacing}
\addtocontents{lof}{{\protect\itshape Number}\protect\hfill{\protect\itshape Page}\protect\par}
\addtocontents{lot}{{\protect\itshape Number}\protect\hfill{\protect\itshape Page}\protect\par}

\setlength{\cftbeforechapterskip}{0pt}

%% For formatting chapter numbers as Roman but
%% not touch appendix alphabets
\RequirePackage{xstring}
\newcommand*{\IsInteger}[3]{%
    \IfStrEq{#1}{ }{#3}{%
        \IfInteger{#1}{#2}{#3}%
    }%
}%
\newcommand{\formatchapternumber}[1]{%
  \IsInteger{#1}{\uppercase\expandafter{\romannumeral #1\relax}}{#1}%
}
\renewcommand{\chapternumberline}[1]{\chaptername\space\formatchapternumber{#1}:\space}
\patchcmd{\cftchapterleader}{\bfseries}{\normalfont}{}{}
\renewcommand{\cftchapterfont}{\normalfont}
\renewcommand{\cftchapterpagefont}{\normalfont}
\renewcommand{\cftchapterdotsep}{\cftdotsep}
\setlength{\cftsectionnumwidth}{2em}

% Merged list of figures and tables -- not needed
% \def\table{\def\figurename{Table}\figure}
% \let\endtable\endfigure

\setlength{\cftbeforefigureskip}{0pt}
\renewcommand{\cftfiguredotsep}{\cftdotsep}
\settowidth\cftfigurenumwidth{\itshape Number}
\renewcommand\cftfigurepresnum{\hfill}
\renewcommand\cftfigureaftersnum{\hskip1em}

\setlength{\cftbeforetableskip}{0pt}
\renewcommand{\cfttabledotsep}{\cftdotsep}
\settowidth\cfttablenumwidth{\itshape Number}
\renewcommand\cfttablepresnum{\hfill}
\renewcommand\cfttableaftersnum{\hskip1em}

\apptocmd{\appendix}{%
  \addtocontents{toc}{%
    \protect\renewcommand{\protect\chaptername}{Appendix}
  }%
}{}{}

\RequirePackage{letterspace} % not so heavy...

% Chapter headings
\renewcommand{\insertchapterspace}{}
\setlength{\midchapskip}{\onelineskip}
\setlength{\afterchapskip}{1.5\onelineskip}
\renewcommand{\chapternamenum}{\space\space}
\renewcommand{\chapterheadstart}{\clearforchapter}
\renewcommand{\chapnamefont}{\SingleSpacing\centering\normalfont\large\itshape\textls[250]}
\renewcommand{\chapnumfont}{\normalfont\large\itshape}
\renewcommand{\chaptitlefont}{\SingleSpacing\centering\normalfont\large}
\renewcommand{\printchaptertitle}[1]{\chaptitlefont\MakeUppercase{#1}}

% Add a command to create an extra chapter with TOC entry.
\newcommand{\extrachapter}[2][\@empty]{%
  \chapter*{#2}%
  \addcontentsline{toc}{chapter}{\ifx#1\@empty{#2}\else{#1}\fi}}

% Section headings
\usepackage{titlesec}
\titleformat*{\section}{\normalsize\bfseries}
\titleformat*{\subsection}{\normalsize\itshape}

\aliaspagestyle{chapter}{headings}

%
% Define spacing for use on titlepage
%
\def\titleskip{\vskip 4\bigskipamount}      %% Whitespace after title
\def\authorskip{\vskip 2\bigskipamount}     %% Whitespace after author

%
% Additional titlepage definitions
%
\newcommand{\degreeaward}{\gdef\@degreeaward}
\degreeaward{Doctor of Philosophy}                  % Default is Ph.D.

\newcommand{\department}{\gdef\@department}
\degreeaward{Mechanical and Aerospace Engineering}  %Default is Mechanical and Aerospace Engineering

\newcommand{\university}{\gdef\@university}
\university{Case Western Reserve University}     % Default is CWRU

\newcommand{\unilogo}{\gdef\@unilogo}
\unilogo{cwru_logo.eps}                                  % CWRU logo

\newcommand{\address}{\gdef\@address}
\address{Cleveland, Ohio}                      % Default is CWRU

\newcommand{\committeeChair}{\gdef\@committeeChair}
\newcommand{\committeeOne}{\gdef\@committeeOne}
\newcommand{\committeeTwo}{\gdef\@committeeTwo}
\newcommand{\committeeThree}{\gdef\@committeeThree}
\newcommand{\numberOfMembers}{\gdef\@numberOfMembers}
\numberOfMembers{3}

\newcommand{\copyyear}{\gdef\@copyyear}
\newcommand{\defenddate}{\gdef\@defenddate}
\newcommand{\defendmonth}{\gdef\@defendmonth}

%Make Committee Page
\newcommand{\makecommittee}{
\begin{center}
    \textbf{\MakeUppercase{Case Western Reserve University}}\\
    \textbf{\MakeUppercase{School of Graduate Studies}}
    \authorskip
    We hereby approve the thesis/dissertation of \par
    \textbf{\@author} \par
    candidate for the degree of \textbf{\@degreeaward}\footnote{We also certify that written approval has been obtained for any proprietary material contained therein.}. \authorskip
    Committee Chair \\
    \textbf{\@committeeChair} \authorskip
    Committee Member \\
    \textbf{\@committeeOne} \authorskip
    Committee Member \\
    \textbf{\@committeeTwo} \authorskip
    Committee Member \\
    \textbf{\@committeeThree} \authorskip
    Date of Defense \\
    \textbf{\@defenddate}
\end{center}
}

\renewcommand{\maketitle}[1][]{%
 \begin{titlingpage}
 \large
   \OnehalfSpacing
   \let\footnotesize\small \let\footnoterule\relax \setcounter{page}{1}
   \vskip-\headsep
   \begin{center}
     {\SingleSpacing\huge\MakeUppercase\@title\par}
     \titleskip
     by
     \titleskip
     {\Large \MakeUppercase\@author}\par
     \authorskip
     Submitted in partial fulfillment of the requirements for the degree of\\
     \@degreeaward
     \authorskip
     \@department\par
     \ifstrequal{#1}{logo}{
       \vspace{\stretch{1}}
       \includegraphics[width=3in]{\@unilogo} \\
       \vspace{\stretch{1}}
     }
     {\vfill}
      \MakeUppercase{\@university}\\
     \authorskip
     \@defendmonth\\
   \end{center}
   \vspace*{.5in}
 \end{titlingpage}
 \pagenumbering{roman}
 \setcounter{page}{2}
 \makecommittee
}               %% End of \maketitle


% The abstract and acknowledgments are treated as unnumbered chapters
% with Table of Contents entries.
\renewenvironment{abstract}
  {\extrachapter[Abstract]{ }\begin{center}
  \authorskip
  \chaptitlefont \@title \authorskip Abstract \authorskip
      by \authorskip
      \MakeUppercase\@author
  \end{center}}
  {\clearpage}
\newenvironment{dedication}
  {\chapter*{Dedication}}
  {\clearpage}
\newenvironment{preface}
  {\extrachapter{Preface}}
  {\clearpage}
\newenvironment{acknowledgements}
  {\extrachapter{Acknowledgements}}
  {\clearpage}

%Double space text
\linespread{1.6}

\nobibintoc

\endinput
