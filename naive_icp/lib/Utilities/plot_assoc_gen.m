function plot_assoc_gen(scA, scB, Apts, Bpts)
pa = make_association_lines(scA, scB, Apts, Bpts);
pt = make_association_lines(scA, scB, Apts, Bpts);
%plot(scA(1,:),scA(2,:),'r+', scB(1,:),scB(2,:),'b+');
if ~isempty(pa)
    plot(pa(1,:),pa(2,:),'g');
end
if ~isempty(pt)
    plot(pt(1,:),pt(2,:),'m');
end