function p = make_association_lines(a, b, Apts,Bpts)
p = [];
for i=1:length(Apts)
    if Bpts(i)~=0
        p = [p a(:, Apts(i)) b(:,Bpts(i)) [NaN;NaN]];
    end
end