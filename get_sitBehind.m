function sitBehind = get_sitBehind(pos, c)
sitBehind = [];
x = pos;
while ~isempty(x)
    sitBehind = union(sitBehind, x);
    x = find(ismember(c, x));
    x = setdiff(x, sitBehind);
end