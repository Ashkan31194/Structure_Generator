function [OptimizationFunctionValue,dOptimizationFunction] = OptimizationFunction(OptimizationFunctionParameters)

N=(size(OptimizationFunctionParameters,1)-2)/2;
KisiXYZ=OptimizationFunctionParameters(1:N,:);
Angle=OptimizationFunctionParameters(N+1:2*N,:);
L_Box=OptimizationFunctionParameters(2*N+1,:);
Angle_Box=OptimizationFunctionParameters(2*N+2,:);
[Sigma,dSigmadKisiXYZ,dSigmadAngle,dSigmadLBox,dSigmadAngleBox] = GradSigma(KisiXYZ,Angle,L_Box,Angle_Box);
OptimizationFunctionValue=Sigma;
dOptimizationFunction=[dSigmadKisiXYZ;dSigmadAngle;dSigmadLBox;dSigmadAngleBox];

end

