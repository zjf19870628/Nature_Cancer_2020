function MakeFatTriangle_ALL (Data, Colors, FileOut)
%%
File=readtable('MEL.txt','Delimiter','\t')
Data=table2array(File(:,2:size(File,2)))
patients = table2array(File(:,1))';
pi3kpatients = {'Pat_01','Pat_07','Pat_08','Pat_11','Pat_16','Pat_24','Pat_29','Pat_30','Pat_44','Pat_48','Pat_58','Pat_60','Pat_63','Pat_70','Pat_74'};
BrafNraspatients = {'Pat_01','Pat_04','Pat_05','Pat_06','Pat_07','Pat_08','Pat_14','Pat_15','Pat_16','Pat_22','Pat_24','Pat_26','Pat_30','Pat_31','Pat_32','Pat_40','Pat_41','Pat_44','Pat_45','Pat_46','Pat_48','Pat_58','Pat_59','Pat_60','Pat_63','Pat_70','Pat_73','Pat_74','Pat_76'};
MapKpatients = {'Pat_04','Pat_05','Pat_15','Pat_31','Pat_32','Pat_41','Pat_45','Pat_63','Pat_76'};

[sharedVals,pi3kIntoA] = intersect(patients,pi3kpatients);
[sharedVals,brafIntoA] = intersect(patients,BrafNraspatients);
[sharedVals,mapkIntoA] = intersect(patients,MapKpatients);
[sharedVals,mapkpi3kIntoA] = intersect(pi3kIntoA, mapkIntoA);
Colors=repmat('k',size(patients,2),1);
Colors(pi3kIntoA)='r';
Colors(mapkIntoA)='b';
Colors(mapkpi3kIntoA)='g';
% Colors(brafIntoA)='g';

dx = Data(:,1)';
rel = Data(:,2)';
com = Data(:,3)';

fig = figure;
a=0;
[x,y,z] = sphere(20);
ll = (x >= a & y >= a & z >= a);
hh = surf(x.*ll,y.*ll,z.*ll, 'FaceColor','w', 'FaceLighting', 'phong', 'LineStyle', '-', 'LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
view([1,1,1]);
xlim([0 1])
axis off;
axis equal;
hold on;
%%
t = [0: 0.0001 : pi/2];
a = sin(t);
b = cos(t);
plot3(a,b,0*t,'k','LineWidth',3);
plot3(0*t,a,b,'k','LineWidth',3);
plot3(b,0*t,a,'k','LineWidth',3);
hold on;
%%
a=[com; dx; rel];
b=1./(com.^2 + dx.^2 + rel.^2)'*ones(1,3);
x2s = a.*sqrt(b');

[~, ls] = sort(sum(a), 'descend');

%%
for ii=1:length(ls)
    i= ls(ii);
    Msize=log2(com(i)+rel(i)+dx(i))+1;
    ss = 0.2;
    plot3(ss+x2s(2,i),ss+x2s(3,i),ss+x2s(1,i),'o','LineWidth',1.5,'MarkerEdgeColor', [0.5 0.5 0.5],'MarkerFaceColor', Colors(i) ,'MarkerSize',Msize,'DisplayName',patients{i});
    %plot3(ss+x2s(2,i),ss+x2s(3,i),ss+x2s(1,i),'ko','LineWidth',1.5, 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor','r', 'MarkerSize',Msize);
end
h1 = plot(NaN,NaN,'ro');set(h1,'MarkerFaceColor','r');
h2 = plot(NaN,NaN,'bo');set(h2,'MarkerFaceColor','b');
h3 = plot(NaN,NaN,'go');set(h3,'MarkerFaceColor','g');
legend([h1,h2,h3],'PI3K','MAPK','PI3K+MAPK');
% h4 = plot(NaN,NaN,'go');set(h4,'MarkerFaceColor','g');
% legend(h4,'BRAF+NRAS');
%%
% if (FileOut)
%     saveas(gca, sprintf ('%s.eps', FileOut), 'eps2c');
%     print(fig, '-dtiff', '-r600', sprintf ('%s.tif', FileOut)); 
% end
