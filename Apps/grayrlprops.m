function stats = grayrlprops(varargin)

%GRAYCOPROPS Properties of gray-level run-length matrix.
%  -------------------------------------------
%  STATS = GRAYCOPROPS(GLRLM,PROPERTIES) Each element in  GLRLM, (r,c),
%   is the probability occurrence of pixel having gray level values r, run-length c in the image.
%   GRAYCOPROPS is to calculate PROPERTIES.
%  Current supported statistics include:
%  -------------------------------------------
%   Short Run Low Gray-Level Emphasis (SRLGE)
%   Short Run High Gray-Level Emphasis (SRHGE)
%   Long Run Low Gray-Level Emphasis (LRLGE)
%   Long Run High Gray-Level Emphasis (LRHGE)
%  --------------------------------------------


% Check GLRLM
[GLRLM numGLRLM] = ParseInputs(varargin{:});

% Initialize output stats structure.
% 11 statistics for each GLRLM
numStats = 4;

% % count number of GLRLM
% numGLRLM = length(GLRLM);

% Initialization default 4*11 matrix
stats = zeros(numGLRLM,numStats);

for p = 1 : numGLRLM
    %N-D indexing not allowed for sparse.

    if numGLRLM ~= 1
        % transfer to double matrix
        tGLRLM = GLRLM{p};
    else
        tGLRLM = GLRLM;
    end
    %     if numGLRLM ~= 1
    %         % transfer to double matrix
    %         tGLRLM = normalizeGLRL(GLRLM{p});
    %     else
    %         tGLRLM = normalizeGLRL(GLRLM);
    %     end
    % Get row and column subscripts of GLRLM.  These subscripts correspond to the
    % pixel values in the GLRLM.
    s = size(tGLRLM);
    % colum indicator
    c_vector =1:s(1);
    % row indicator
    r_vector =1:s(2);
    % matrix element indicator
    % Matrix form col and row: using meshgrid, you should transpose before using
    [c_matrix,r_matrix] = meshgrid(c_vector,r_vector);

    % Total number of runs
    N_runs = sum(sum(tGLRLM));

    % total number of elements
    N_tGLRLM = s(1)*s(2);

    %--------------------Prepare four matrix for speedup--------------
    % 1.Gray Level Run-Length Pixel Number Matrix
    %     p_p = calculate_p_p(tGLRLM,c_matrix');

    % 2.Gray-Level Run-Number Vector
    %   This vector represents the sum distribution of the number of runs
    %   with gray level i.
    p_g = sum(tGLRLM);

    % 3.Run-Length Run-Number Vector
    %   This vector represents the sum distribution of the number of runs
    %   with run length j.
    p_r = sum(tGLRLM,2)';

    % 4.Gray-Level Run-Length-One Vector
    %
    % p_o = tGLRLM(:,1); % Not used yet
    % ----------------------End four matrix---------------------------
    %
    %------------------------Statistics-------------------------------

    % 8. Short Run Low Gray-Level Emphasis (SRLGE)
    SGLGE =calculate_SGLGE(tGLRLM,r_matrix',c_matrix',N_runs);
    % 9. Short Run High Gray-Level Emphasis (SRHGE)
    SRHGE =calculate_SRHGE(tGLRLM,r_matrix',c_matrix',N_runs);
    % 10. Long Run Low Gray-Level Emphasis (LRLGE)
    LRLGE =calculate_LRLGE(tGLRLM,r_matrix',c_matrix',N_runs);
    % 11.Long Run High Gray-Level Emphasis (LRHGE
    LRHGE =calculate_LRHGE(tGLRLM,r_matrix',c_matrix',N_runs);
    %----------------insert statistics----------------------------
    stats(p,:)=[ SGLGE SRHGE LRLGE  LRHGE ];
end % end all run length matrixs


%---------------------------------
function SGLGE =calculate_SGLGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Short Run Low Gray-Level Emphasis (SRLGE):

term = tGLRLM./((r_matrix.*c_matrix).^2);
SGLGE= sum(sum(term))./N_runs;

%------------------------------------
function  SRHGE =calculate_SRHGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Short Run High Gray-Level Emphasis (SRHGE):
%
term  = tGLRLM.*(r_matrix.^2)./(c_matrix.^2);
SRHGE = sum(sum(term))/N_runs;
%------------------------------------
function   LRLGE =calculate_LRLGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Long Run Low Gray-Level Emphasis (LRLGE):
%
term  = tGLRLM.*(c_matrix.^2)./(r_matrix.^2);
LRLGE = sum(sum(term))/N_runs;
%---------------------------------------
function  LRHGE =calculate_LRHGE(tGLRLM,r_matrix,c_matrix,N_runs)
% Long Run High Gray-Level Emphasis (LRHGE):

term  = tGLRLM.*(c_matrix.^2).*(r_matrix.^2);
LRHGE = sum(sum(term))/N_runs;

%-----------------------------------------------------------------------------
function [glrlm num_glrlm] = ParseInputs(varargin)
% check stability of inputs
%
% first receive all inputs
glrlm = varargin{:};
% get numbers total
num_glrlm=length(glrlm);
% then for each element, check its stability
for i=1:num_glrlm
    % The 'nonnan' and 'finite' attributes are not added to iptcheckinput because the
    % 'integer' attribute takes care of these requirements.
    % iptcheckinput(glrlm,{'cell'},{'real','nonnegative','integer'}, ...
    % mfilename,'GLRLM',1);
    iptcheckinput(glrlm{i},{'logical','numeric'},{'real','nonnegative','integer'},...
        mfilename,'GLRLM',1);
    % Cast GLRLM to double to avoid truncation by data type. Note that GLRLM is not an
    % image.
    if ~isa(glrlm,'double')
        glrlm{i}= double(glrlm{i});
    end
end