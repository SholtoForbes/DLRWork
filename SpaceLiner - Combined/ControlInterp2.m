function alphainterp = ControlInterp2(t,AlphaList,tsearch)

alphainterp = interp1(t,AlphaList,tsearch,'previous'); 
% alphainterp = interp1(t,AlphaList,tsearch,'next');
% alphainterp = spline(t,AlphaList,tsearch);
end