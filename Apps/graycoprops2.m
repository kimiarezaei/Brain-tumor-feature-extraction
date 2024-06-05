function stats = graycoprops2(varargin)
%GRAYCOPROPS Properties of gray-level co-occurrence matrix.  
%   STATS = GRAYCOPROPS(GLCM,PROPERTIES) normalizes the gray-level
%   co-occurrence matrix (GLCM) so that the sum of its elements is one. Each
%   element in the normalized GLCM, (r,c), is the joint probability occurrence
%   of pixel pairs with a defined spatial relationship having gray level
%   values r and c in the image. GRAYCOPROPS uses the normalized GLCM to
%   calculate PROPERTIES.
%
%   GLCM can be an m x n x p array of valid gray-level co-occurrence
%   matrices. Each gray-level co-occurrence matrix is normalized so that its
%   sum is one.
%
%   PROPERTIES can be a comma-separated list of strings, a cell array
%   containing strings, the string 'all', or a space separated string. They
%   can be abbreviated, and case does not matter.
%
%   Properties include:
%  
%   'Contrast'      the intensity contrast between a pixel and its neighbor 
%                   over the whole image. Range = [0 (size(GLCM,1)-1)^2]. 
%                   Contrast is 0 for a constant image.
%
%   'Correlation'   statistical measure of how correlated a pixel is to its 
%                   neighbor over the whole image. Range = [-1 1]. 
%                   Correlation is 1 or -1 for a perfectly positively or
%                   negatively correlated image. Correlation is NaN for a 
%                   constant image.
%
%   'Energy'        summation of squared elements in the GLCM. Range = [0 1].
%                   Energy is 1 for a constant image.
%  
%   'Homogeneity'   closeness of the distribution of elements in the GLCM to
%                   the GLCM diagonal. Range = [0 1]. Homogeneity is 1 for
%                   a diagonal GLCM.
%  


allStats = {'variance','entropy','Energy','IDM'};
  
[glcm, requestedStats] = ParseInputs(allStats, varargin{:});

% Initialize output stats structure.
numStats = length(requestedStats);
numGLCM = size(glcm,3);
empties = repmat({zeros(1,numGLCM)},[numStats 1]);
stats = cell2struct(empties,requestedStats,1);

for p = 1 : numGLCM
  
  if numGLCM ~= 1 %N-D indexing not allowed for sparse. 
    tGLCM = normalizeGLCM(glcm(:,:,p));
  else 
    tGLCM = normalizeGLCM(glcm);
  end
  
  % Get row and column subscripts of GLCM.  These subscripts correspond to the
  % pixel values in the GLCM.
  s = size(tGLCM);
  [c,r] = meshgrid(1:s(1),1:s(2));
  r = r(:);
  c = c(:);

  % Calculate fields of output stats structure.
  for k = 1:numStats
    name = requestedStats{k};  
    switch name
     case 'variance'
      stats.(name)(p) = calculatevariance(tGLCM,r,c);
      
     case 'entropy'
      stats.(name)(p) = calculateentropy(tGLCM);
      
     case 'Energy'
      stats.(name)(p) = calculateEnergy(tGLCM);
      
     case 'IDM'
      stats.(name)(p) = calculateIDM(tGLCM,r,c);
    end
  end

end


%-----------------------------------------------------------------------------
function glcm = normalizeGLCM(glcm)
  
% Normalize glcm so that sum(glcm(:)) is one.
if any(glcm(:))
  glcm = glcm ./ sum(glcm(:));
end
  
%-----------------------------------------------------------------------------
function ent = calculateentropy(glcm)
term1=glcm(:);
term2=log(glcm(:));
term = term1 .* term2;
ent =-sum(term);

%-----------------------------------------------------------------------------
function S = stdIndex(index,glcm,m)

term1 = (index - m).^2 .* glcm(:);
S = sqrt(sum(term1));

%-----------------------------------------------------------------------------
function M = meanIndex(index,glcm)

M = index .* glcm(:);
M = sum(M);

%----------------------------------------------------------------------------
function v = calculatevariance(glcm,r,c)
% Reference: Haralick RM, Shapiro LG. Computer and Robot Vision: Vol. 1,
% Addison-Wesley, 1992, p. 460.  
mr = meanIndex(r,glcm);
term1=(r-mr).^2;
term2=glcm(:);
term = term1 .* term2;
v = sum(term);

%-------------------------------------------------------------
function E = calculateEnergy(glcm)
% Reference: Haralick RM, Shapiro LG. Computer and Robot Vision: Vol. 1,
% Addison-Wesley, 1992, p. 460.  
  
foo = glcm.^2;
E = sum(foo(:));

%-----------------------------------------------------------------------------
function IDM = calculateIDM(glcm,r,c)
term1=1+(r-c).^2;
term2=1/term1;
term3=glcm(:);
term=term2 * term3;
IDM=sum(term);



%-----------------------------------------------------------------------------
function [glcm,reqStats] = ParseInputs(allstats,varargin)
  
numstats = length(allstats);
iptchecknargin(1,numstats+1,nargin,mfilename);

reqStats = '';
glcm = varargin{1};

% The 'nonnan' and 'finite' attributes are not added to iptcheckinput because the
% 'integer' attribute takes care of these requirements.
iptcheckinput(glcm,{'logical','numeric'},{'real','nonnegative','integer'}, ...
              mfilename,'GLCM',1);

if ndims(glcm) > 3
  eid = sprintf('Images:%s:invalidSizeForGLCM',mfilename);
  msg = 'GLCM must either be an m x n or an m x n x p array.';
  error(eid,'%s',msg);
end

% Cast GLCM to double to avoid truncation by data type. Note that GLCM is not an
% image.
if ~isa(glcm,'double')
  glcm = double(glcm);
end

list = varargin(2:end);

if isempty(list)
  % GRAYCOPROPS(GLCM) or GRAYCOPROPS(GLCM,PROPERTIES) where PROPERTIES is empty.
  reqStats = allstats;
else
  if iscell(list{1}) || numel(list) == 1
    % GRAYCOPROPS(GLCM,{...})
    list = list{1};
  end

  if ischar(list)
    %GRAYCOPROPS(GLCM,SPACE-SEPARATED STRING)
    list = strread(list, '%s');
  end

  anyprop = allstats;
  anyprop{end+1} = 'all';
  
    match = iptcheckstrs(list{k}, anyprop, mfilename, 'PROPERTIES', k+1);
  for k = 1 : length(list)
    if strcmp(match,'all')
      reqStats = allstats;
      break;
    end
    reqStats{k} = match;
  end
  
end

% Make sure that reqStats are in alphabetical order.
reqStats = sort(reqStats);

if isempty(reqStats)
  eid = sprintf('Images:%s:internalError',mfilename);
  msg = 'Internal error: requestedStats has not been assigned.';
  error(eid,'%s',msg);
end
