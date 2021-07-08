clc;
clear;
%% Read data from the Input file
fileID=fopen('Input_File.txt','r');
fgetl(fileID);
formatSpec='%s%f';
A=textscan(fileID,formatSpec,[1 2]);
global N_Types
N_Types=A{2};
for i=1:N_Types
    formatSpec='%s%f';
    AA=textscan(fileID,formatSpec,[1 2]);
    formatSpec='%s';
    B=textscan(fileID,formatSpec,[1 1]);
    formatSpec='%f';
    C=fscanf(fileID,formatSpec,[1 Inf]);
    fgetl(fileID);
    EdgesType{AA{2}}=C;
    for j=1:size(C,2)
        textscan(fileID,'%s',[1 1]);
        global ConnectorID
        if AA{2}==C(j)
            ConnectorID{AA{2}}{j}(1,:)=fscanf(fileID,formatSpec,[1 Inf]);
            Xo=textscan(fileID,'%s',[1 1]);
            Xoo=cell2mat(Xo{1});
            Yo=fscanf(fileID,formatSpec,[1 Inf]);
            ConnectorID{AA{2}}{j}(2,:)=[str2num(Xoo(2)) Yo];
            textscan(fileID,'%s',[1 1]);
        else
            ConnectorID{AA{2}}{j}=fscanf(fileID,formatSpec,[1 Inf]);
        end
    end
    formatSpec='%s%f';
    AAA=textscan(fileID,formatSpec,[1 2]);
    textscan(fileID,'%s',[1 1]);
    for j=1:AAA{2}
        X=fscanf(fileID,'%f',[1 3]);
        global ConnectorCord
        ConnectorCord{AA{2}}{j}=X;
    end   
    BBB=textscan(fileID,formatSpec,[1 2]);
    QQ=textscan(fileID,'%s',[1 1]);
    for j=1:BBB{2}
        ff=textscan(fileID,'%s',[1 1]);
        global AtomsType
        AtomsType{AA{2}}{j}=ff;
        y=fscanf(fileID,'%f',[1 3]);
        global AtomsCord
        AtomsCord{AA{2}}{j}=y;
    end      
    
end    
fgetl(fileID);    
fgetl(fileID);
fgetl(fileID);
formatSpec='%s%f';
Q=textscan(fileID,formatSpec,[1 2]);
global N_Vertices
N_Vertices=Q{2};
fgetl(fileID);
fgetl(fileID);
for i=1:N_Vertices
    T=fscanf(fileID,'%f',[1 Inf]);
    global Type
    Type(T(1))=T(2);
    fscanf(fileID,'%s',[1 1]);
     for j=1:T(3)
    TT=fscanf(fileID,'%f',[1 Inf]);
    global Net_Edges
    Net_Edges{i}{j}(1:2)=TT;
    CI=fscanf(fileID,'%s',[1 1]);
    indCI=[];
    counter=1;
    for k=1:size(CI,2)
        if CI(k)=='(' | CI(k)==')' | CI(k)==','
            indCI(counter)=k;
            counter=counter+1;
        end
    end
    CI(indCI)=[];
    if size(CI,2)==3
        Net_Edges{i}{j}(3:5)=[str2num(CI(1)) str2num(CI(2)) str2num(CI(3))];
    else
        FC=find(CI=='-');
        if FC==1
           Net_Edges{i}{j}(3:5)=[-str2num(CI(2)) str2num(CI(3)) str2num(CI(4))];
        elseif FC==2
           Net_Edges{i}{j}(3:5)=[str2num(CI(1)) -str2num(CI(3)) str2num(CI(4))];
        elseif FC==3
           Net_Edges{i}{j}(3:5)=[str2num(CI(1)) str2num(CI(2)) -str2num(CI(4))];
        end
    end

     end
end
fgetl(fileID);
fgetl(fileID);
for i=1:N_Vertices
    Initial_Coordinates(i,:)=fscanf(fileID,'%f',[1 6]);
end
fgetl(fileID);
fgetl(fileID);

clear C D

formatSpec='%f%f%s%s';
C=textscan(fileID,formatSpec,[3 4]);
formatSpec='%f';
D=textscan(fileID,formatSpec);

xlo=C{1}(1);
ylo=C{1}(2);
zlo=C{1}(3);
xhi=C{2}(1);
yhi=C{2}(2);
zhi=C{2}(3);
xy=D{1}(1);
xz=D{1}(2);
yz=D{1}(3);
InitialBox=[xhi-xlo xy xz;0 yhi-ylo yz;0 0 zhi-zlo];
Origin=[xlo ylo zlo];

clear A AA AAA B C CI counter D FC fileID formatSpec i indCI j k Q T TT X ans
%% Maincode
 XYZ0=Initial_Coordinates(:,1:3);
 Angle0=Initial_Coordinates(:,4:6);
 Box0=InitialBox;

XYZ=XYZ0;%%+rand(size(XYZ0));
Angle=Angle0;
%Angle=2*pi*rand(size(Angle0));
Box=Box0;
L1=sqrt(Box(1,1)^2+Box(2,1)^2+Box(3,1)^2);
L2=sqrt(Box(1,2)^2+Box(2,2)^2+Box(3,2)^2);
L3=sqrt(Box(1,3)^2+Box(2,3)^2+Box(3,3)^2);
L_Box=[L1 L2 L3];
alpha=acos(dot(Box(:,2),Box(:,3))/(L2*L3));
beta=acos(dot(Box(:,1),Box(:,3))/(L1*L3));
gamma=acos(dot(Box(:,2),Box(:,1))/(L2*L1));
Angle_Box=[alpha beta gamma];

Vol=L_Box(1)*L_Box(2)*L_Box(3)*sqrt(1-(cos(Angle_Box(1)))^2 -(cos(Angle_Box(2)))^2 -(cos(Angle_Box(3)))^2 +2*cos(Angle_Box(1))*cos(Angle_Box(2))*cos(Angle_Box(3)))
Cart_to_frac=[1/L_Box(1) -(cot(Angle_Box(3)))/L_Box(1) L_Box(2)*L_Box(3)*((cos(Angle_Box(1))*cos(Angle_Box(3))-cos(Angle_Box(2)))/(Vol*sin(Angle_Box(3))));...
    0 1/(L_Box(2)*sin(Angle_Box(3))) L_Box(1)*L_Box(3)*((cos(Angle_Box(2))*cos(Angle_Box(3))-cos(Angle_Box(1)))/(Vol*sin(Angle_Box(3))));...
    0 0 (L_Box(1)*L_Box(2)*sin(Angle_Box(3)))/Vol];


%Angle_Box=pi*rand(1,3)/2;


Vol=L_Box(1)*L_Box(2)*L_Box(3)*sqrt(1-(cos(Angle_Box(1)))^2 -(cos(Angle_Box(2)))^2 -(cos(Angle_Box(3)))^2 +2*cos(Angle_Box(1))*cos(Angle_Box(2))*cos(Angle_Box(3)))
frac_to_Cart=[L_Box(1) L_Box(2)*cos(Angle_Box(3)) L_Box(3)*cos(Angle_Box(2));...
    0 L_Box(2)*sin(Angle_Box(3)) L_Box(3)*((cos(Angle_Box(1))-cos(Angle_Box(2))*cos(Angle_Box(3)))/(sin(Angle_Box(3))));...
    0 0 Vol/(L_Box(1)*L_Box(2)*sin(Angle_Box(3)))];

KisiXYZ=(Cart_to_frac*XYZ')';


Sigma=0;
dSigmadKisiXYZ=zeros(size(KisiXYZ));
dSigmadAngle=zeros(size(Angle));
dSigmadLBox=zeros(size(L_Box));
dSigmadAngleBox=zeros(size(Angle_Box));...
Box=frac_to_Cart;



XYZ_Updated=(frac_to_Cart*KisiXYZ')';
Box_Updated=Box;

figure(1)
scatter3(XYZ_Updated(:,1),XYZ_Updated(:,2),XYZ_Updated(:,3),'filled','r')
hold on
PPLL=[Origin;...
    Origin+Box_Updated(:,1)';...
    Origin+Box_Updated(:,1)'+Box_Updated(:,2)';...
    Origin+Box_Updated(:,2)';...
    Origin;...
    Origin+Box_Updated(:,3)';...
    Origin+Box_Updated(:,1)'+Box_Updated(:,3)';...
    Origin+Box_Updated(:,1)';...
    Origin+Box_Updated(:,1)'+Box_Updated(:,3)';...
    Origin+Box_Updated(:,1)'+Box_Updated(:,2)'+Box_Updated(:,3)';...
    Origin+Box_Updated(:,1)'+Box_Updated(:,2)';...
    Origin+Box_Updated(:,1)'+Box_Updated(:,2)'+Box_Updated(:,3)';...
    Origin+Box_Updated(:,2)'+Box_Updated(:,3)';...
    Origin+Box_Updated(:,2)';...
    Origin+Box_Updated(:,2)'+Box_Updated(:,3)';...
    Origin+Box_Updated(:,3)';...
    Origin];

 plot3(PPLL(:,1),PPLL(:,2),PPLL(:,3),'b')
 


OptimizationFunctionParameters=[KisiXYZ;Angle;L_Box;Angle_Box];
 
tic
[OptimizationFunctionValue,dOptimizationFunction] = OptimizationFunction(OptimizationFunctionParameters);
toc  

% 
tic
options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true);
x0 =OptimizationFunctionParameters;
fun = @OptimizationFunction;
[x,fval] = fminunc(fun,x0,options)
toc

KisiXYZ=x(1:N_Vertices,:);
Angle=x(N_Vertices+1:2*N_Vertices,:);
L_Box=x(2*N_Vertices+1,:);
Angle_Box=x(2*N_Vertices+2,:);
 
 
Vol=L_Box(1)*L_Box(2)*L_Box(3)*sqrt(1-(cos(Angle_Box(1)))^2 -(cos(Angle_Box(2)))^2 -(cos(Angle_Box(3)))^2 +2*cos(Angle_Box(1))*cos(Angle_Box(2))*cos(Angle_Box(3)))

Cart_to_frac=[1/L_Box(1) -(cot(Angle_Box(3)))/L_Box(1) L_Box(2)*L_Box(3)*((cos(Angle_Box(1))*cos(Angle_Box(3))-cos(Angle_Box(2)))/(Vol*sin(Angle_Box(3))));...
    0 1/(L_Box(2)*sin(Angle_Box(3))) L_Box(1)*L_Box(3)*((cos(Angle_Box(2))*cos(Angle_Box(3))-cos(Angle_Box(1)))/(Vol*sin(Angle_Box(3))));...
    0 0 (L_Box(1)*L_Box(2)*sin(Angle_Box(3)))/Vol];

frac_to_Cart=[L_Box(1) L_Box(2)*cos(Angle_Box(3)) L_Box(3)*cos(Angle_Box(2));...
    0 L_Box(2)*sin(Angle_Box(3)) L_Box(3)*((cos(Angle_Box(1))-cos(Angle_Box(2))*cos(Angle_Box(3)))/(sin(Angle_Box(3))));...
    0 0 Vol/(L_Box(1)*L_Box(2)*sin(Angle_Box(3)))];

XYZ=(frac_to_Cart*KisiXYZ')';
Box=frac_to_Cart;


XYZ_Updated=XYZ;
Angle_Updated=Angle;
Box_Updated=Box;



for i=1:N_Vertices
    T=Type(i);
    RCenter=XYZ_Updated(i,:);
    C=ConnectorCord{T};
    RR=[];
    for j=1:size(C,2)
        RR=[RR;C{j}];
    end
    QQ=(rotz((180/pi)*Angle_Updated(i,1))*roty((180/pi)*Angle_Updated(i,2))*rotx((180/pi)*Angle_Updated(i,3))*RR')';
    AAQQ=RCenter+QQ;
    %scatter3(AAQQ(:,1),AAQQ(:,2),AAQQ(:,3),'k')
    FFDD{i}=((Box_Updated^-1)*(RCenter'+QQ'))';
end

maxX=-Inf;
maxY=-Inf;
maxZ=-Inf;
minX=Inf;
minY=Inf;
minZ=Inf;

for i=1:size(FFDD,2)
    if max(FFDD{i}(:,1))>maxX
        maxX=max(FFDD{i}(:,1));
    end
    if max(FFDD{i}(:,2))>maxY
        maxY=max(FFDD{i}(:,2));
    end    
    if max(FFDD{i}(:,3))>maxZ
        maxZ=max(FFDD{i}(:,3));
    end
    
    if min(FFDD{i}(:,1))<minX
        minX=min(FFDD{i}(:,1));
    end
    if min(FFDD{i}(:,2))<minY
        minY=min(FFDD{i}(:,2));
    end    
    if min(FFDD{i}(:,3))<minZ
        minZ=min(FFDD{i}(:,3));
    end
end
figure(2)
Origin=(Box_Updated*[minX; minY; minZ])';
scatter3(XYZ_Updated(:,1),XYZ_Updated(:,2),XYZ_Updated(:,3),'filled','r')
hold on
PPLL=[Origin;...
    Origin+Box_Updated(:,1)';...
    Origin+Box_Updated(:,1)'+Box_Updated(:,2)';...
    Origin+Box_Updated(:,2)';...
    Origin;...
    Origin+Box_Updated(:,3)';...
    Origin+Box_Updated(:,1)'+Box_Updated(:,3)';...
    Origin+Box_Updated(:,1)';...
    Origin+Box_Updated(:,1)'+Box_Updated(:,3)';...
    Origin+Box_Updated(:,1)'+Box_Updated(:,2)'+Box_Updated(:,3)';...
    Origin+Box_Updated(:,1)'+Box_Updated(:,2)';...
    Origin+Box_Updated(:,1)'+Box_Updated(:,2)'+Box_Updated(:,3)';...
    Origin+Box_Updated(:,2)'+Box_Updated(:,3)';...
    Origin+Box_Updated(:,2)';...
    Origin+Box_Updated(:,2)'+Box_Updated(:,3)';...
    Origin+Box_Updated(:,3)';...
    Origin];

 plot3(PPLL(:,1),PPLL(:,2),PPLL(:,3),'b')
 for i=1:size(FFDD,2)
     AFG=(Box*FFDD{i}')';
     hold on
     scatter3(AFG(:,1),AFG(:,2),AFG(:,3),'*','k')
 end
 
AllAtoms=[]; 
AllTypes=[];
for i=1:N_Vertices
    T=Type(i);
    QCenter=XYZ_Updated(i,:);
    C=AtomsCord{T};
    D=AtomsType{T};
    RR=[];
    for j=1:size(C,2)
        RR=[RR;C{j}];
    end
    QQ=(rotz((180/pi)*Angle_Updated(i,1))*roty((180/pi)*Angle_Updated(i,2))*rotx((180/pi)*Angle_Updated(i,3))*RR')';
    AllAtoms=[AllAtoms;QCenter+QQ];
    for j=1:size(D,2)
        AllTypes=[AllTypes;D{j}{1}];
    end
end 
 
X=unique(AllTypes);
Y=[1:size(X,1)]';
for i=1:size(AllTypes,1)
    TTI=AllTypes{i};
    for j=1:size(X,1)
        if X{j}==TTI
            res=j;
        end
    end
    AllTypesq(i,:)=res;
end
 outputfile = fopen('Atoms_Coordinates.data','w');
 fprintf(outputfile,'%s\n\n','Input_File');
 fprintf(outputfile,'%d atoms\n\n%d atom types\n\n%f %f xlo xhi\n%f %f ylo yhi \n%f %f zlo zhi\n%f %f %f xy xz yz\n\n',...
                     size(AllAtoms,1),size(X,1),Origin(1),Origin(1)+Box_Updated(1,1),Origin(2),Origin(2)+Box_Updated(2,2),Origin(3),Origin(3)+Box_Updated(3,3),Box_Updated(1,2),Box_Updated(1,3),Box_Updated(2,3));
 fprintf(outputfile,'Masses\n\n');
 for i=1:size(X,1)
     fprintf(outputfile,'%s %f\n',X{i},1);
 end
 fprintf(outputfile,'\nAtoms\n\n');
 for i=1:size(AllAtoms,1)
     fprintf(outputfile,'%d %d %f %f %f\n',i,AllTypesq(i),AllAtoms(i,1),AllAtoms(i,2),AllAtoms(i,3));
 end
     
 fclose(outputfile);
