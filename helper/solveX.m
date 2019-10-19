function [sol, outs] = solveX(A, At, b0, x0, opts)
switch lower(opts.algorithm)
    case 'custom'
        [sol, outs] = optsCustomAlgorithm(A, At, b0, x0, opts);
    case 'amplitudeflow'
        [sol, outs] = solveAmplitudeFlow(A, At, b0, x0, opts);
    case 'coordinatedescent'
        [sol, outs] = solveCoordinateDescent(A, At, b0, x0, opts);
    case 'fienup'
        [sol, outs] = solveFienup(A, At, b0, x0, opts);
    case 'gerchbergsaxton'
        [sol, outs] = solveGerchbergSaxton(A, At, b0, x0, opts);
    case 'kaczmarz'
        [sol, outs] = solveKaczmarzSimple(A, At, b0, x0, opts);
    case 'phasemax'
        [sol, outs] = solvePhaseMax(A, At, b0, x0, opts);
    case 'phaselamp'
        [sol, outs] = solvePhaseLamp(A, At, b0, x0, opts);
    case 'phaselift'
        [sol, outs] = solvePhaseLift(A, At, b0, x0, opts);
    case 'raf'
        [sol, outs] = solveRAF(A, At, b0, x0, opts);
    case 'rwf'
        [sol, outs] = solveRWF(A, At, b0, x0, opts);
    case 'sketchycgm'
        [sol, outs] = solveSketchyCGM(A, At, b0, x0, opts);
    case 'taf'
        [sol, outs] = solveTAF(A, At, b0, x0, opts);
    case 'twf'
        [sol, outs] = solveTWF(A, At, b0, x0, opts);
    case {'wirtflow','wf'}
        [sol, outs] = solveWirtFlow(A, At, b0, x0, opts);
    otherwise
        error('Unknown algorithm "%s"', opts.algorithm);
end
end