function alphainterp = ControlInterp(t,AlphaList,tsearch)

% alphainterp = interp1(t,AlphaList,tsearch,'previous');
alphainterp = interp1(t,AlphaList,tsearch,'next');
% alphainterp = interp1(t,AlphaList,tsearch);
% alphainterp = spline(t,AlphaList,tsearch);
end