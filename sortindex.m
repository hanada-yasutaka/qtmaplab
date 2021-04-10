function index = sortindex(targetmat, refmat)
arguments
    targetmat (:, :) {mustBeNumeric}
    refmat    (:, :) {mustBeNumeric}
end
if size(targetmat) ~= size(refmat)
    error("size must be same");
end

mat = conj(targetmat).' * refmat;
mat2 = conj(mat) .* mat;

dim = length(targetmat);
qnumber = 1:dim;
index = zeros(1,dim);
for i=1:dim
    [~, ind] = max(mat2(:,i));
    index(i) = qnumber(ind);
    %index = [index, qnumber(ind)];
    qnumber(ind) = [];
    mat2(ind,:) = [];
end

if length(unique(index)) ~= dim
    error("index is not unique");
end
end
