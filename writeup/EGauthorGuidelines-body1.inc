% ---------------------------------------------------------------------
% EG author guidelines plus sample file for EG publication using LaTeX2e input
% D.Fellner, v1.17, Sep 23, 2010


\title[Mesh denoising via generalized intrinsic filter on 2d manifold]%
    {Mesh denoising via generalized intrinsic filter on 2d manifold}
    %{Mesh denoising via intrinsic normal filter}

% for anonymous conference submission please enter your SUBMISSION ID
% instead of the author's name (and leave the affiliation blank) !!

%\author[B. Liu \& J. J. Cao \&N. Ma \& W. M. Wang \& L. G. Liu \& X. P. Liu]
       %{B. Liu$^{1}$  W. M. Wang$^{1}$  J. J. Cao$^{1}$  N. Ma$^{1}$  L. G. Liu$^{2}$ and
        %X. P. Liu\thanks{xpliu@dlut.edu.cn}$^{1}$\\
         %$^1$Dalian University of Technology, China\\
         %$^2$University of Science and Technology of China, China
       %}

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
    %The bilateral filter is a classical filtering operator with the goal of preserving edge structure while smoothing signals.
    %The reason it works is mainly on the spatial and range distances are enough close, but it ignores the local relations between signals.
    %In this paper, we introduce a intrinsic filter method for 2D manifold surfaces.
    %This approach builds the connections between desired signal and its neighbors.% shown in figure~\ref{Fig:relation}.
    %Therefore, it has a outstanding power on filtering neighbor signal while preserving the important feature.
    %A novel filtering model based on our intrinsic filter is proposed with application in mesh denoising.
    %As the traditional mesh filtering, our framework also has two-stage process:
    %first, we apply geodesic path to build the face normal connections, then filter the face normals basing two different accumulative distance weights.
    %Afterwards, the vertex positions are updated according to the filtered face normals.
    %For accelerating our intrinsic algorithm, we also introduce a simple, fast and effective pattern to compensate the time deficiency in using geodesic path.
    %It performs well on a wide variety of meshes and is competitive with other state-of-the-art methods.
    
%Most mesh denoising methods are weighted average normal filters and their performance depends on the weight design. 
%When computing the weight between current face and one of its neighbour, existing methods only consider properties of the two faces, such as positions and normals. 
%However if the two faces belong to different regions which are separated by a narrow feature edge, 
%even when the two properties are close, it is improper to assign a large weight which will damage the feature.

%In this paper, we present an intrinsic mesh normal filter. 
%It estimates the weight between the current face and its neighbour based on the integral of two kinds of normal differences along the geodesic path connecting their barycenters. 
%Therefore, features are better preserved when removing noises. 
%It is also a very generalized model and many classic filters are just special cases of it. 
%Furthermore, by projection on the tangent plane, we introduce a simple pattern to estimate the weight which are faster and more effective than using geodesic path. 
%We apply the filter to mesh denoising and experiments illustrate the efficacy of our method comparing with other state-of-the-art methods.

Weighted average is the most common strategy used in various digital signal filters and their performance depend on the weight design. 
When computing the weight between the current point and one of its neighbours, existing methods consider only properties of the two points, such as positions and signals. 
However if the two points belong to different regions which are separated by a narrow feature edge, even when the two properties are close, 
It is improper to assign a large weight which will damage the feature.

In this paper, we present an intrinsic filter on 2d manifold. 
It estimates the weight between the current point and its neighbour based on the integral of some properties along the geodesic path connecting them. 
Therefore, features are better preserved when removing noises. 
It is also a very generalized model and many classic filters are just special cases of it. 
Furthermore, by projections on the tangent plane, we introduce a simple pattern to estimate the weight which are faster and more effective than using geodesic path. 
Finally, we apply the filter to mesh denoising and experiments illustrate the efficacy of our method comparing with other state-of-the-art methods.

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
