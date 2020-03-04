function out = reshape4montage(in)

if length(size(in)) ~= 3
    error('This function is specfically for reshaping 3D arrays for montage function.')
else
    out = reshape(in, size(in,1), size(in,2), 1, size(in,3));
end