function check = isNearValue(A,B)
% The ability to compare numeric values to within a tolerance is not available as a built-in function in MATLAB 8.0 (R2012b). As a workaround to avoid issues resulting from roundoff error, you can compare the absolute difference of the operands to a tolerance. Instead of:
% A==B
% you can use the following:
% abs(A-B) < 1e4*eps(min(abs(A),abs(B)))
% There are also a couple of entries for this on the MATLAB Central file exchange. These files are not produced by or supported by MathWorks, so any questions should be directed to the authors, however other customers have found these to be useful.
check = abs(A-B) < 1e4*eps(min(abs(A),abs(B)));

end