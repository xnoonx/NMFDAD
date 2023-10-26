function [wbest,hbest,normbest,i] = cnmf(a,k,num,varargin)
%     cnmf--Constrained Non-negative Matrix Factorization.
%     This function is based on MATLAB NNMF function. Constraints are added
%     to deal with spectra and eluting peaks. For more infomation please
%     refer to our paper:
%     Xu N, Hu M, Li X et al. Resolving ultraviolet-visible spectral for
%     complex dissolved mixtures of multitudinous organic matters in aerosol.
%     Analytical Chemistry.

%     [W,H] = cnmf(A,K,N) factors the N-by-M matrix A into non-negative
%     factors W (N-by-K) and H (K-by-M). N represents the number of spliced
%     samples.
%     [W,H,D] = cnmf(...) also returns D, the root mean square residual.
%     [W,H,D,i] = cnmf(...) also returns i, the final iteration times.%
%     [...] = cnmf(A,K,'PARAM1',val1,'PARAM2',val2,...) specifies one or
%     more of the following parameter name/value pairs:
%
%     Parameter    Value
%     'algorithm'  Either 'als' (default) to use an alternating least
%     squares algorithm, or 'mult' to use a multiplicative
%     update algorithm.
%     'w0'         An N-by-K matrix to be used as the initial value for W.
%     'h0'         A K-by-M matrix to be used as the initial value for H.
%     'replicates' The number of times to repeat the factorization, using
%     new random starting values for W and H, except at the
%     first replication if w0 and h0 are given (default 1).
%     This tends to be most beneficial with the 'mult'
%     algorithm.
%     'options'  An options structure as created by the STATSET
%     function.  cnmf uses the following fields:
%
%     'Display'   Level of display output.  Choices are 'off'
%     (the default), 'final', and 'iter'.
%     'MaxIter'   Maximum number of steps allowed. The default
%     is 100.  Unlike in optimization settings,
%     reaching MaxIter is regarded as convergence.
%     'TolFun'    Positive number giving the termination tolerance
%     for the criterion.  The default is 1e-4.
%     'TolX'      Positive number giving the convergence threshold
%     for relative change in the elements of W and H.
%     The default is 1e-4.
%     'UseParallel'
%     'UseSubstreams'
%     'Streams'   These fields specify whether to perform multiple
%     replicates in parallel, and how to use random
%     numbers when generating the starting points for
%     the replicates. For information on these fields
%     see PARALLELSTATS.

% Check required arguments
if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

narginchk(3,Inf);
[n,m] = size(a);
if ~isscalar(k) || ~isnumeric(k) || k<1 || k>min(m,n) || k~=round(k)
    error(message('stats:nnmf:BadK'));
end

% Process optional arguments
pnames = {'algorithm' 'w0' 'h0' 'replicates' 'options'};
dflts  = {'als'       []   []   1            []       };
[alg,w0,h0,tries,options] = ...
    internal.stats.parseArgs(pnames,dflts,varargin{:});

% Check optional arguments
alg = internal.stats.getParamVal(alg,{'mult' 'als'},'ALGORITHM');
ismult = strncmp('mult',alg,numel(alg));
checkmatrices(a,w0,h0,k);
if ~isscalar(tries) || ~isnumeric(tries) || tries<1 || tries~=round(tries)
    error(message('stats:nnmf:BadReplicates'));
end

defaultopt = statset('nnmf');
tolx = statget(options,'TolX',defaultopt,'fast');
tolfun = statget(options,'TolFun',defaultopt,'fast');
maxiter = statget(options,'MaxIter',defaultopt,'fast');
%maxiter=1000;
dispopt = statget(options,'Display',defaultopt,'fast');

[~,dispnum] = internal.stats.getParamVal(dispopt, {'off','notify','final','iter'},'Display');
dispnum = dispnum - 1;

[useParallel, RNGscheme, poolsz] = ...
    internal.stats.parallel.processParallelAndStreamOptions(options,true);

usePool = useParallel && poolsz>0;

% Special case, if K is full rank we know the answer
if isempty(w0) && isempty(h0)
    if k==m
        w0 = a;
        h0 = eye(k);
    elseif k==n
        w0 = eye(k);
        h0 = a;
    end
end


% Define the function that will perform one iteration of the
% loop inside smartFor
loopbody = @loopBody;

% Suppress undesired warnings.
if usePool
    % On workers and client
    pctRunOnAll internal.stats.parallel.muteParallelStore('rankDeficientMatrix', ...
        warning('off','MATLAB:rankDeficientMatrix') );
else
    % On client
    ws = warning('off','MATLAB:rankDeficientMatrix');
end

% Prepare for in-progress
if dispnum > 1 % 'iter' or 'final'
    if usePool
        % If we are running on a parallel pool, each worker will generate
        % a separate periodic report.  Before starting the loop, we
        % seed the parallel pool so that each worker will have an
        % identifying label (eg, index) for its report.
        internal.stats.parallel.distributeToPool( ...
            'workerID', num2cell(1:poolsz) );
        
        % Periodic reports behave differently in parallel than they do
        % in serial computation (which is the baseline).
        % We advise the user of the difference.
        warning(message('stats:nnmf:displayParallel2'));
        
        % Leave formatted by \t UI strings untranslated. 8/17/2011
        fprintf('    worker\t      rep\t   iteration\t     rms resid\t     |delta x|\n' );
    else
        if useParallel
            warning(message('stats:nnmf:displayParallel'));
        end
        fprintf('    rep\t   iteration\t   rms resid\t  |delta x|\n');
    end
end

try
    whbest = internal.stats.parallel.smartForReduce(...
        tries, loopbody, useParallel, RNGscheme, 'argmin');
catch ME
    % Revert warning setting for rankDeficientMatrix to value prior to nnmf.
    if usePool
        % On workers and on client
        pctRunOnAll warning(internal.stats.parallel.statParallelStore('rankDeficientMatrix').state,'MATLAB:rankDeficientMatrix');
    else
        % On client
        warning(ws);
    end
    rethrow(ME);
end

normbest = whbest{1};
wbest = whbest{3};
hbest = whbest{4};
i=whbest{2};
% whbest{2} contains the iteration chosen for the best factorization,
% but it has no meaning except as a "reproducible" tie-breaker, and
% is not supplied as a return value.

if dispnum > 1   % 'final' or 'iter'
    fprintf('%s\n',getString(message('stats:nnmf:FinalRMSResidual',sprintf('%g',normbest))));
end

% Revert warning setting for rankDeficientMatrix to value prior to nnmf.
if usePool
    % On workers and on client
    pctRunOnAll warning(internal.stats.parallel.statParallelStore('rankDeficientMatrix').state,'MATLAB:rankDeficientMatrix');
else
    % On client
    warning(ws);
end

if normbest==Inf
    error(message('stats:nnmf:NoSolution'))
end

% Put the outputs in a standard form - first normalize the max of w to 1

wlen = max(wbest);
wbest = wbest./(wlen+eps(wlen));
hbest = hbest.*(wlen+eps(wlen))';

% ---- Nested functions ----

    function cellout = loopBody(iter,S)
        if isempty(S)
            S = RandStream.getGlobalStream;
        end
        
        % whtry is a "temporary variable" and hence needs to be
        % reinitialized at start of each loop.
        whtry = cell(4,1);
        
        % Get random starting values if required
        if( ~isempty(w0) && iter ==1 )
            whtry{3} = w0;
        else
            whtry{3} = rand(S,n,k);
        end
        if( ~isempty(h0) && iter ==1 )
            whtry{4} = h0;
        else
            whtry{4} = rand(S,k,m);
        end
        
        % Perform a factorization
        [whtry{3},whtry{4},whtry{1},whtry{2}] = ...
            nnmf1(a,k,num,whtry{3},whtry{4},ismult,maxiter,tolfun,tolx,...
            dispnum,iter,usePool);
        %whtry{2} = iter;
        
        cellout = whtry;
    end

end  % of nnmf

% -------------------
function [w,h,dnorm,j] = nnmf1(a,k,num,w0,h0,ismult,maxiter,tolfun,tolx,...
    dispnum,repnum,usePool)
% Single non-negative matrix factorization
nm = numel(a);
sqrteps = sqrt(eps);

% Display progress.  For parallel computing, the replicate number will be
% displayed under the worker performing the replicate.
if dispnum>1 % 'final' or 'iter'
    if usePool
        labindx = internal.stats.parallel.workerGetValue('workerID');
        dispfmt = '%8d\t%8d\t%8d\t%14g\t%14g\n';
    else
        dispfmt = '%7d\t%8d\t%12g\t%12g\n';
    end
end
[~,m]=size(a);

for j=1:maxiter
    
    if ismult
        % Multiplicative update formula
        numer = a*h0';
        w = max(0,w0 .* (numer ./ (w0*(h0*h0') + eps(numer))));
        numer = w'*a;
        h = max(0,h0 .* (numer ./ ((w'*w)*h0 + eps(numer))));
        
    else
        % Alternating least squares
        w = max(0, a/h0);
        h = max(0, w\a);
    end
    
    % Add unimodal constraints for elution peaks in each sample
    mn=m/num;
    for ii=1:num
        hn=h(:,(ii-1)*mn+1:ii*mn);
        for kk = 1:k
            [M,I]=max(hn(kk,:));
            FWHM=sum(hn(kk,:)>M/2);
            for r=I-20:-1:1
                if hn(kk,r)<M/1.5
                    if hn(kk,r)-hn(kk,r+5)>0
                        if max(hn(kk,1:r))-hn(kk,r)>M/10  || max(abs(hn(kk,max(r-round(FWHM),1):r)-hn(kk,r)))<M/50
                            hn(kk,1:r)=0;
                            break
                        end
                    end
                end
            end
            for r=I+20:mn
                if hn(kk,r)<M/1.5
                    if  hn(kk,r)-hn(kk,r-5)>0
                        if max(hn(kk,r:mn))-hn(kk,r)>M/10  || max(abs(hn(kk,r:min(r+round(FWHM),mn))-hn(kk,r)))<M/50
                            hn(kk,r:mn)=0;
                            break
                        end
                    end
                end
            end
        end
        h(:,(ii-1)*mn+1:ii*mn)=hn;
    end
    
    % Get norm of difference and max change in factors
    d = a - w*h;
    dnorm = sqrt(sum(sum(d.^2))/nm);
    dw = max(max(abs(w-w0) / (sqrteps+max(max(abs(w0))))));
    dh = max(max(abs(h-h0) / (sqrteps+max(max(abs(h0))))));
    delta = max(dw,dh);
    
    % Check for convergence
    if j>2
        if delta <= tolx
            break;
        elseif dnorm0-dnorm <= tolfun*max(1,dnorm0)
            break;
        elseif j==maxiter
            break
        end
    end
    
    if dispnum>2 % 'iter'
        if usePool
            fprintf(dispfmt,labindx,repnum,j,dnorm,delta);
        else
            fprintf(dispfmt,repnum,j,dnorm,delta);
        end
    end
    
    % Remember previous iteration results
    w0 =w;
    h0 = h;
    dnorm0 = dnorm;
end

% Smooth data for each sample
for ii=1:num
    h(:,(ii-1)*mn+1:ii*mn)= smoothdata(h(:,(ii-1)*mn+1:ii*mn),2,'lowess',20);
end

if dispnum>1   % 'final' or 'iter'
    if usePool
        fprintf(dispfmt,labindx,repnum,j,dnorm,delta);
    else
        fprintf(dispfmt,repnum,j,dnorm,delta);
    end
end

end

% ---------------------------
function checkmatrices(a,w,h,k)
% check for non-negative matrices of the proper size

if ~ismatrix(a) || ~isnumeric(a) || ~isreal(a) || any(any(~isfinite(a)))
    error(message('stats:nnmf:BadA'))
end
[n,m] = size(a);
if ~isempty(w)
    if ~ismatrix(w) || ~isnumeric(w)|| ~isreal(w) || any(any(w<0)) || any(any(~isfinite(w)))
        error(message('stats:nnmf:BadWNegativeValues'))
    elseif ~isequal(size(w),[n k])
        error(message('stats:nnmf:BadWSizeIsWrong', sprintf( '%d', n ), sprintf( '%d', k )));
    end
end
if ~isempty(h)
    if ~ismatrix(h) || ~isnumeric(h)|| ~isreal(h) || any(any(h<0)) || any(any(~isfinite(h)))
        error(message('stats:nnmf:BadHNegativeValues'))
    elseif ~isequal(size(h),[k m])
        error(message('stats:nnmf:BadHSizeIsWrong', sprintf( '%d', k ), sprintf( '%d', m )));
    end
end
end % checkmatrices

