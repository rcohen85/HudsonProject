function varargout = corrplot(varargin)
%CORRPLOT Plot variable correlations
%
% Syntax:
%
%   corrplot(X)
%   corrplot(Tbl)
%   [R,PValue,H] = corrplot(...)
%   [...] = corrplot(...,param,val,...)
%   [...] = corrplot(ax,...)
%
% Description:
%
%   CORRPLOT creates a numVars-by-numVars grid of scatter plots between
%   pairs of variables. Histograms of individual variables appear along the
%   diagonal. The scatter plots display a least-squares reference line with
%   slope equal to a displayed correlation coefficient.
%
% Input Arguments:
%
%   X - Time series data, specified as a numObs-by-numVars numeric matrix.
%
%   Tbl - Time series data, specified as a table or timetable. Specify
%       series for X using the 'DataVariables' parameter.
%
%   ax - Axes object in which to plot. If unspecified, CORRPLOT plots to
%        the current axes (gca). CORRPLOT does not support uiaxes targets.
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME        VALUE
%
%   'Type'      String or character vector indicating the type of
%               correlation coefficient to compute. Values are:
%
%               'Pearson'  Pearson's linear correlation coefficient
%
%               'Kendall'  Kendall's rank correlation coefficient (tau)
%
%               'Spearman' Spearman's rank correlation coefficient (rho)
%
%               The default is 'Pearson'.
%
%   'Rows'      String or character vector indicating how to treat NaN
%               values in the data. Values are:
%
%               'all'       Use all rows, regardless of NaNs
%
%               'complete'  Use only rows with no NaNs
%
%               'pairwise'  Use rows with no NaNs in column i or j to
%                           compute R(i,j)
%
%               The default is 'pairwise'.
%
%   'Tail'      String or character vector indicating the alternative
%               hypothesis used to compute the PValue output. Values are:
%
%               'both'	Ha: Correlation is not zero
%
%               'right'	Ha: Correlation is greater than zero
%
%               'left'  Ha: Correlation is less than zero
%
%               The default is 'both'.
%
%	'VarNames' 	String vector or cell vector of character vectors of
%               length numVars of unique variable names to be used in the
%               plots. Names are truncated to the first five characters.
%               The default for matrix input X is {'var1','var2',...}. The
%               default for tabular input Tbl is Tbl.Properties.VarNames.
%
%   'TestR'     String or character vector indicating whether or not to
%               test for significant correlations and highlight them in
%               red. Values are 'off' and 'on'. The default is 'off'.
%
%   'Alpha'     Scalar level for tests of correlation significance. Values
%               must be between 0 and 1. The default value is 0.05.
%
%   'DataVariables' Variables in Tbl to use for X, specified as names in
%               Tbl.Properties.VariableNames. Variable names are cell
%               vectors of character vectors, string vectors, integer
%               vectors or logical vectors. The default is all variables
%               in Tbl.
%
% Output Arguments:
%
%	R - Correlations of the displayed plots. R is a numVars-by-numVars
%       matrix or table, corresponding to the input type.
%
%	PValue - P-values of the displayed plots, used to test the hypothesis
%       of no correlation against the alternative of a nonzero correlation
%       (see Notes). PValue is a numVars-by-numVars matrix or table,
%       corresponding to the input type.
%
%   H - Handles to the displayed graphics objects. H is a numVars-by-numVars
%       matrix or table, corresponding to the input type.
%
% Notes:
%
%   o P-values for Pearson's correlation are computed by transforming the
%     correlation to create a t statistic with numObs-2 degrees of freedom.
%     The transformation is exact when X is normal. P-values for Kendall's
%     and Spearman's rank correlations are computed using either the exact
%     permutation distributions (for small sample sizes), or large-sample
%     approximations. P-values for two-tailed tests are computed by
%     doubling the more significant of the two one-tailed p-values.
%
%   o Using the 'pairwise' option for the 'rows' parameter may return a
%     correlation matrix that is not positive definite. The 'complete'
%     option always returns a positive definite matrix, but in general the
%     estimates are based on fewer observations.
%
% Example:
%
%   load Data_Canada
%   corrplot(DataTable)
%
% See also COLLINTEST, CORR.

% Copyright 2023 The MathWorks, Inc.

% Preprocess varargin for target axes:

try
    
    [ax,args] = internal.econ.axesparser(varargin{:});
    
catch ME
    
    throw(ME)
    
end

% This function produces a single plot:

if ~isempty(ax) && ~isscalar(ax)
    
    error(message('econ:internal:econ:axesparser:InvalidParent'));
    
end

% Input arguments:

Data = args{1};
args = args(2:end);

isTabular = istable(Data) || istimetable(Data);

% Parse inputs and set defaults:

parseObj = inputParser;
parseObj.addRequired('Data',...
                     @(x)validateattributes(x,{'double','table','timetable'},{'nonempty','2d'}));
parseObj.addParameter('Type','Pearson',@typeCheck);
parseObj.addParameter('Rows','pairwise',@rowsCheck);
parseObj.addParameter('Tail','both',@tailCheck);
parseObj.addParameter('VarNames',{},@varNamesCheck);
parseObj.addParameter('TestR','off',@testRCheck);
parseObj.addParameter('Alpha',0.05,@alphaCheck);
parseObj.addParameter('DataVariables',[],...
                      @(x)validateattributes(x,{'double','logical','cell','string'},{'vector'}));

try
    
  parseObj.parse(Data,args{:});
  
catch ME
    
  throwAsCaller(ME)
  
end

Data = parseObj.Results.Data;
corrType = validatestring(parseObj.Results.Type,{'pearson','kendall','spearman'});
whichRows = validatestring(parseObj.Results.Rows,{'all','complete','pairwise'});
tail = validatestring(parseObj.Results.Tail,{'both','right','left'});

varNames = parseObj.Results.VarNames;

if length(varNames) ~= length(unique(varNames))

    error(message('econ:corrplot:NonUniqueVariableNames'))

end

varNames = cellstr(varNames);

testR = validatestring(parseObj.Results.TestR,{'off','on'});
testRFlag = strcmpi(testR,'on');

alpha = parseObj.Results.Alpha;
varSpec = parseObj.Results.DataVariables;

% Select X with 'DataVariables':

if isnumeric(Data)
    
   	X = Data;
    
    if isvector(Data)
        
        error(message('econ:corrplot:DataIsVector'))
        
    end
    
    if ~isempty(varSpec)
        
        warning(message('econ:corrplot:DataVariablesUnused'))
        
    end
    
else % Tabular data

    if ~isempty(varSpec)

        try

            X = Data(:,varSpec);

        catch ME

            throwAsCaller(ME)

        end

    else

        X = Data; % Default

    end
    
    try

        internal.econ.TableAndTimeTableUtilities.isTabularFormatValid(X,'X')
        internal.econ.TableAndTimeTableUtilities.isTabularDataSinglePath(X,'X')

    catch ME

        throwAsCaller(ME)

    end

    XVarNames = X.Properties.VariableNames;
    X = table2array(X);
    X = double(X);

end

numVars = size(X,2);

% Create variable names:

if isempty(varNames)
    
    if isTabular

    	varNames = XVarNames;

    else

        varNames = strcat({'var'},num2str((1:numVars)','%-u'));
        
    end

else

    if length(varNames) < numVars

        error(message('econ:corrplot:VarNamesTooFew'))

    elseif length(varNames) > numVars

        error(message('econ:corrplot:VarNamesTooMany'))
        
    end

end

% Truncate variable names to first five characters for display:

varNames = cellfun(@(s)[s,'         '],varNames,'UniformOutput',false);
varNames = cellfun(@(s)s(1:6),varNames,'UniformOutput',false);

% Compute plot information:

[R,PValue] = corr(X,'type',corrType,'rows',whichRows,'tail',tail);

Mu = mean(X,"omitnan");
Sigma = std(X,"omitnan");

if any(Sigma < eps)
    
    error(message('econ:corrplot:NoVariationInData'))
    
end

Z = bsxfun(@minus,X,Mu);
Z = bsxfun(@rdivide,Z,Sigma);
ZLims = [min(Z(:),[],"omitnan"),max(Z(:),[],"omitnan")];

% Plot matrix:

ax = newplot(ax);

[hS,Ax,bigAx,hH] = plotmatrix(ax,X,'bo');

set(hS,'MarkerSize',2)
set(hH,'FaceColor','b','FaceAlpha',0.6,'EdgeColor','c')

for i = 1:numVars
        
    for j = 1:numVars
        
        hFig = ancestor(bigAx,'figure');
        
        set(hFig,'CurrentAxes',Ax(i,j))
        set(Ax(i,j),'XLim',Mu(j)+(1.1)*ZLims*Sigma(j),...
                    'YLim',Mu(i)+(1.1)*ZLims*Sigma(i))
        axis(Ax(i,j),'normal')   
            
        if i ~= j
            
            hls = lsline(Ax(i,j));
            set(hls,'Color','m','Tag','lsLines');
            
            if testRFlag && (PValue(i,j) < alpha)
                
                corrColor = 'r';
                
            else
                
                corrColor = 'k';
                
            end

            text(Ax(i,j),0.1,0.9,...
                 num2str(R(i,j),'%3.2f'),...
                 'Units','normalized',...
                 'FontWeight','bold',...
                 'Color',corrColor,...
                 'Tag','corrCoefs')
            
        end
        
    end
    
end

% Set axes properties:

ax.Tag = 'CorrPlot';

Xlabels = gobjects(numVars,1);
Ylabels = gobjects(numVars,1);

for i = 1:numVars
    
    Xlabels(i) = xlabel(Ax(numVars,i),varNames{i});
    Ylabels(i) = ylabel(Ax(i,1),varNames{i});

end

set(get(bigAx,'Title'),'String','{\bf Correlation Matrix}')

% Restore current axes to bigAx:

set(hFig,'CurrentAxes',bigAx)

% Return plot object:

H = hS;
H(logical(eye(numVars))) = hH;

% Create tabular output:

if isTabular
    
    R = array2table(R,'RowNames',varNames,'VariableNames',varNames);
    PValue = array2table(PValue,'RowNames',varNames,'VariableNames',varNames);
    H = array2table(H,'RowNames',varNames,'VariableNames',varNames);
    
end

% Suppress assignment to ans:

nargoutchk(0,3);

if nargout > 0
    
    varargout = {R,PValue,H};
    
end

%-------------------------------------------------------------------------
% Check value of 'type' parameter
function OK = typeCheck(corrType)

if ~isvector(corrType)

    error(message('econ:corrplot:CorrTypeNonVector'))

elseif isnumeric(corrType)

    error(message('econ:corrplot:CorrTypeNumeric'))

end

OK = true;

%-------------------------------------------------------------------------
% Check value of 'rows' parameter
function OK = rowsCheck(whichRows)

if ~isvector(whichRows)

    error(message('econ:corrplot:RowsParamNonVector'))

elseif isnumeric(whichRows)

    error(message('econ:corrplot:RowsParamNumeric'))

end

OK = true;

%-------------------------------------------------------------------------
% Check value of 'tail' parameter
function OK = tailCheck(tail)

if ~isvector(tail)

    error(message('econ:corrplot:TailParamNonVector'))

elseif isnumeric(tail)

    error(message('econ:corrplot:TailParamNumeric'))

end

OK = true;

%-------------------------------------------------------------------------
% Check value of 'varNames' parameter
function OK = varNamesCheck(varNames)
    
if ischar(varNames) || ~isvector(varNames) || ...
   isnumeric(varNames) || (iscell(varNames) && any(cellfun(@isnumeric,varNames)))

    error(message('econ:corrplot:VarNamesInvalid'))

end

OK = true;

%-------------------------------------------------------------------------
% Check value of 'testR' parameter
function OK = testRCheck(testR)

if ~isvector(testR)

    error(message('econ:corrplot:testRNonVector'))

elseif isnumeric(testR)

    error(message('econ:corrplot:testRNumeric'))

end

OK = true;

%-------------------------------------------------------------------------
% Check value of 'alpha' parameter
function OK = alphaCheck(alpha)
    
if ~isnumeric(alpha)

    error(message('econ:corrplot:AlphaNonNumeric'))

elseif ~isscalar(alpha)

    error(message('econ:corrplot:AlphaNonScalar'))

elseif alpha < 0 || alpha > 1

    error(message('econ:corrplot:AlphaOutOfRange'))

end

OK = true;
