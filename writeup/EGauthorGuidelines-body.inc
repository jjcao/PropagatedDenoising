% ---------------------------------------------------------------------
% EG author guidelines plus sample file for EG publication using LaTeX2e input
% D.Fellner, v1.17, Sep 23, 2010


\title[EG \LaTeX\ Author Guidelines]%
      {Intrinsic filtering for 2D manifold surface}

% for anonymous conference submission please enter your SUBMISSION ID
% instead of the author's name (and leave the affiliation blank) !!
\author[D. Fellner \& S. Behnke]
       {D.\,W. Fellner\thanks{Chairman Eurographics Publications Board}$^{1,2}$
        and S. Behnke$^{2}$
%        S. Spencer$^2$\thanks{Chairman Siggraph Publications Board}
        \\
% For Computer Graphics Forum: Please use the abbreviation of your first name.
         $^1$TU Darmstadt \& Fraunhofer IGD, Germany\\
         $^2$Institut f{\"u}r ComputerGraphik \& Wissensvisualisierung, TU Graz, Austria
%        $^2$ Another Department to illustrate the use in papers from authors
%             with different affiliations
       }

% ------------------------------------------------------------------------

% if the Editors-in-Chief have given you the data, you may uncomment
% the following five lines and insert it here
%
% \volume{27}   % the volume in which the issue will be published;
% \issue{1}     % the issue number of the publication
% \pStartPage{1}      % set starting page


%-------------------------------------------------------------------------
\begin{document}


\maketitle

\begin{abstract}
    The bilateral filter is a classical filtering operator with the goal of preserving edge structure while smoothing signals.
    The reason it works is mainly on the spatial and range distances are enough close, but it ignores the local relations between signals.
    In this paper, we introduce a intrinsic filter method for 2D manifold surfaces.
    This approach builds the connections between desired signal and its neighbors.% shown in figure~\ref{Fig:relation}.
    Therefore, it has a outstanding power on filtering neighbor signal while preserving the important feature.
    A novel filtering model based on our intrinsic filter is proposed with application in mesh denoising.
    As the traditional mesh filtering, our framework also has two-stage process:
    first, we apply geodesic path to build the face normal connections, then filter the face normals basing two different accumulative distance weights.
    Afterwards, the vertex positions are updated according to the filtered face normals.
    For accelerating our intrinsic algorithm, we also introduce a simple, fast and effective pattern to compensate the time deficiency in using geodesic path.
    It performs well on a wide variety of meshes and is competitive with other state-of-the-art methods.

\begin{classification} % according to http://www.acm.org/about/class/1998
\CCScat{Computer Graphics}{I.3.3}{Picture/Image Generation}{Line and curve generation}
\end{classification}

\end{abstract}


\input{introduction}
\input{related}
\input{algorithm}
\input{results}
%\input{conclusion}

%-------------------------------------------------------------------------

%\bibliographystyle{eg-alpha}
\bibliographystyle{eg-alpha-doi}

\bibliography{egbibsample}

%-------------------------------------------------------------------------

\end{document}
