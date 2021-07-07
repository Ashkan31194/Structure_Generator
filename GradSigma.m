function [Sigma,dSigmadKisiXYZ,dSigmadAngle,dSigmadLBox,dSigmadAngleBox] = GradSigma(KisiXYZ,Angle,L_Box,Angle_Box)


Vol=L_Box(1)*L_Box(2)*L_Box(3)*sqrt(1-(cos(Angle_Box(1)))^2 -(cos(Angle_Box(2)))^2 -(cos(Angle_Box(3)))^2 +2*cos(Angle_Box(1))*cos(Angle_Box(2))*cos(Angle_Box(3)));
frac_to_Cart=[L_Box(1) L_Box(2)*cos(Angle_Box(3)) L_Box(3)*cos(Angle_Box(2));...
    0 L_Box(2)*sin(Angle_Box(3)) L_Box(3)*((cos(Angle_Box(1))-cos(Angle_Box(2))*cos(Angle_Box(3)))/(sin(Angle_Box(3))));...
    0 0 Vol/(L_Box(1)*L_Box(2)*sin(Angle_Box(3)))];
Box=frac_to_Cart;

Sigma=0;
dSigmadKisiXYZ=zeros(size(KisiXYZ));
dSigmadAngle=zeros(size(Angle));
dSigmadLBox=zeros(size(L_Box));
dSigmadAngleBox=zeros(size(Angle_Box));

global  ConnectorID ConnectorCord N_Vertices Type Net_Edges

for i=1:N_Vertices
     gI=Net_Edges{i};
     I1=i;
     for j=1:size(Net_Edges{i},2)
         I2=gI{j}(1);
         cncI1=gI{j}(2);
         gJ=Net_Edges{I2};
         for k=1:size(gJ,2)
             indind(k)=gJ{k}(1);
         end
         indcnc=find(indind==I1);
         cncI2=gJ{indcnc}(2);
         clear indind indcnc
         ehtaI2=gI{j}(3:5);
         typeI1=Type(I1);
         typeI2=Type(I2);
         KisiXYZI1=KisiXYZ(I1,:);
         KisiXYZI2=KisiXYZ(I2,:);
         
         if cncI1>0 && cncI2>0
         rk0I1=ConnectorCord{typeI1}(ConnectorID{typeI1}{cncI1}(1,:));
         rk0I2=ConnectorCord{typeI2}(ConnectorID{typeI2}{cncI2}(1,:));
         elseif cncI1>0 && cncI2<0
         rk0I1=ConnectorCord{typeI1}(ConnectorID{typeI1}{cncI1}(1,:));
         rk0I2=ConnectorCord{typeI2}(ConnectorID{typeI2}{-cncI2}(2,:));
         elseif cncI1<0 && cncI2>0
         rk0I1=ConnectorCord{typeI1}(ConnectorID{typeI1}{-cncI1}(2,:));
         rk0I2=ConnectorCord{typeI2}(ConnectorID{typeI2}{cncI2}(1,:));         
         elseif cncI1<0 && cncI2<0
         rk0I1=ConnectorCord{typeI1}(ConnectorID{typeI1}{-cncI1}(2,:));
         rk0I2=ConnectorCord{typeI2}(ConnectorID{typeI2}{-cncI2}(2,:)); 
         end
         
         
         rI1=(Box*KisiXYZI1');
         rI2=(Box*KisiXYZI2');
         for l=1:size(rk0I1,2)
             rkI1=rI1+rotz((180/pi)*Angle(I1,1))*roty((180/pi)*Angle(I1,2))*rotx((180/pi)*Angle(I1,3))*rk0I1{l}';
             rkI2=rI2+rotz((180/pi)*Angle(I2,1))*roty((180/pi)*Angle(I2,2))*rotx((180/pi)*Angle(I2,3))*rk0I2{l}';
             Sigma=Sigma+0.5*(rkI1-rkI2-Box*ehtaI2')'*(rkI1-rkI2-Box*ehtaI2');
             
             KDI1XYZ=zeros(size(KisiXYZ));
             KDI1XYZ(I1,:)=[1 1 1];
             KDI2XYZ=zeros(size(KisiXYZ));
             KDI2XYZ(I2,:)=[1 1 1];
             First_term=Box'*(rkI1-rkI2-Box*ehtaI2');
             FFTT=repmat(First_term',size(KisiXYZ,1),1);
             dSigmadKisiXYZ=dSigmadKisiXYZ+FFTT.*(KDI1XYZ-KDI2XYZ);
             
             
             First_TermAI1=(rkI1-rkI2-Box*ehtaI2')*rk0I1{l};
             First_TermAI2=(rkI1-rkI2-Box*ehtaI2')*rk0I2{l};
             Second_TermAI1=[-sin(Angle(I1,1)) -cos(Angle(I1,1)) 0;cos(Angle(I1,1)) -sin(Angle(I1,1)) 0;0 0 0]*roty((180/pi)*Angle(I1,2))*rotx((180/pi)*Angle(I1,3));
             Second_TermAI2=[-sin(Angle(I2,1)) -cos(Angle(I2,1)) 0;cos(Angle(I2,1)) -sin(Angle(I2,1)) 0;0 0 0]*roty((180/pi)*Angle(I2,2))*rotx((180/pi)*Angle(I2,3));
             FTerm1A=trace(First_TermAI1'*Second_TermAI1);
             FTerm2A=trace(First_TermAI2'*Second_TermAI2);
             Second_TermBI1=rotz((180/pi)*Angle(I1,1))*[-sin(Angle(I1,2)) 0 cos(Angle(I1,2));0 0 0;-cos(Angle(I1,2)) 0 -sin(Angle(I1,2))]*rotx((180/pi)*Angle(I1,3));
             Second_TermBI2=rotz((180/pi)*Angle(I2,1))*[-sin(Angle(I2,2)) 0 cos(Angle(I2,2));0 0 0;-cos(Angle(I2,2)) 0 -sin(Angle(I2,2))]*rotx((180/pi)*Angle(I2,3));
             FTerm1B=trace(First_TermAI1'*Second_TermBI1);
             FTerm2B=trace(First_TermAI2'*Second_TermBI2);
             Second_TermCI1=rotz((180/pi)*Angle(I1,1))*roty((180/pi)*Angle(I1,2))*[0 0 0;0 -sin(Angle(I1,3)) -cos(Angle(I1,3));0 cos(Angle(I1,3)) -sin(Angle(I1,3))];
             Second_TermCI2=rotz((180/pi)*Angle(I2,1))*roty((180/pi)*Angle(I2,2))*[0 0 0;0 -sin(Angle(I2,3)) -cos(Angle(I2,3));0 cos(Angle(I2,3)) -sin(Angle(I2,3))];
             FTerm1C=trace(First_TermAI1'*Second_TermCI1);
             FTerm2C=trace(First_TermAI2'*Second_TermCI2);
             KDI1Angle=zeros(size(Angle));
             KDI1Angle(I1,:)=[1 1 1];
             KDI2Angle=zeros(size(Angle));
             KDI2Angle(I2,:)=[1 1 1];
             dSigmadAngle=dSigmadAngle+repmat([FTerm1A FTerm1B FTerm1C],size(Angle,1),1).*KDI1Angle-repmat([FTerm2A FTerm2B FTerm2C],size(Angle,1),1).*KDI2Angle;
             
             dSdB=(rkI1-rkI2-Box*ehtaI2')*(KisiXYZI1-KisiXYZI2-ehtaI2);
             dBdL1=[1 0 0;0 0 0;0 0 0];
             dBdL2=[0 cos(Angle_Box(3)) 0;0 sin(Angle_Box(3)) 0; 0 0 0];
             dBdL3=[0 0 cos(Angle_Box(2));0 0 (cos(Angle_Box(1))-cos(Angle_Box(2))*cos(Angle_Box(3)))/sin(Angle_Box(3));0 0 Vol/((L_Box(1)*L_Box(2)*L_Box(3))*sin(Angle_Box(3)))];
             dBdalpha=[0 0 0;0 0 -L_Box(3)*(sin(Angle_Box(1))/sin(Angle_Box(3)));0 0 L_Box(1)*L_Box(2)*(L_Box(3))^2 *((0.5*sin(2*Angle_Box(1))-cos(Angle_Box(3))*cos(Angle_Box(2))*sin(Angle_Box(1)))/(Vol*sin(Angle_Box(3))))];
             dBdbeta=[0 0 -L_Box(3)*sin(Angle_Box(2));0 0 L_Box(3)*sin(Angle_Box(2))*cot(Angle_Box(3));0 0 L_Box(1)*L_Box(2)*(L_Box(3))^2 *((0.5*sin(2*Angle_Box(2))-cos(Angle_Box(3))*cos(Angle_Box(1))*sin(Angle_Box(2)))/(Vol*sin(Angle_Box(3))))];
             dBdgamma=[0 -L_Box(2)*sin(Angle_Box(3)) 0;0 L_Box(2)*cos(Angle_Box(3)) L_Box(3)*(cos(Angle_Box(2))-((cot(Angle_Box(3))/sin(Angle_Box(3)))*(cos(Angle_Box(1))-cos(Angle_Box(2))*cos(Angle_Box(3)))));0 0 L_Box(1)*L_Box(2)*(L_Box(3))^2 *((cos(Angle_Box(3))-cos(Angle_Box(1))*cos(Angle_Box(2)))/Vol)-(1/(L_Box(1)*L_Box(2)))*((cos(Angle_Box(3))*Vol)/(sin(Angle_Box(3)))^2)];
             dSigmadLBox=dSigmadLBox+[trace(dSdB'*dBdL1) trace(dSdB'*dBdL2) trace(dSdB'*dBdL3)];
             dSigmadAngleBox=dSigmadAngleBox+[trace(dSdB'*dBdalpha) trace(dSdB'*dBdbeta) trace(dSdB'*dBdgamma)];
             

         end
     end
 end
end

