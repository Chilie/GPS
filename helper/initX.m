% Initialize x0 using the specified initialization method
function x0 = initX(A, At, b0, n, opts)
switch lower(opts.initMethod)
    case {'truncatedspectral','truncated'}
        x0 = initSpectral(A,At,b0,n,true,true,opts.verbose);  % isTruncated, isScaled?
    case 'spectral'
        x0 = initSpectral(A,At,b0,n,false,true,opts.verbose); % isTruncated, isScaled?
    case {'amplitudespectral','amplitude'}
        x0 = initAmplitude(A,At,b0,n,opts.verbose);
    case {'weightedspectral','weighted'}
        x0 = initWeighted(A,At,b0,n,opts.verbose);
    case {'orthogonalspectral','orthogonal'}
        x0 = initOrthogonal(A,At,b0,n,opts.verbose);
    case {'optimal','optimalspectral'}
        x0 = initOptimalSpectral(A,At,b0,n,true,opts.verbose);
    case 'angle'
        assert(isfield(opts,'xt'),'The true solution, opts.xt, must be specified to use the angle initializer.')
        assert(isfield(opts,'initAngle'),'An angle, opts.initAngle, must be specified (in radians) to use the angle initializer.')
        x0 = initAngle(opts.xt, opts.initAngle);
    case 'custom'
        x0 = opts.customx0;
    otherwise
        error('Unknown initialization method "%s"', opts.initMethod);
end
end