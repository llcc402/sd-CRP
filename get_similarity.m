% This function returns a symmetric matrix that the (i,j)-th element is the
% distance between data(i,:) and data(j,:). 
function S = get_similarity(data, a)
if nargin < 2
    a = 1;
end
mydist = dist(data');
S = exp(-mydist / a);