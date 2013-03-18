clc;clear all;close all;
C=Euler3212C([pi/4,pi/6,pi/9]);
%% convert C to ehat, phi
Phi = acos(.5*(C(1,1)+C(2,2)+C(3,3)-1));

ehat=1/(2*sin(Phi)).*[C(2,3)-C(3,2);
                      C(3,1)-C(1,3);
                      C(1,2)-C(2,1)];
etilde=[0       -ehat(3)     ehat(2);
        ehat(3)     0       -ehat(1);
        -ehat(2) ehat(1)        0       ];


%% Output
fprintf('Numerically verifying C=expm(phi*etilde)\n\n')

C1=expm(-Phi*etilde);
maxdif=max(max(C1-C));
fprintf('Maximum difference was:  %3.5f \n\n\n',maxdif)

fprintf('Now numerically verifying EQ3.82\n\n')
C2=eye(3,3)*cos(Phi)-sin(Phi)*etilde + (1-cos(Phi))*ehat*ehat';

maxdif=max(max(C2-C));
fprintf('Maximum difference was:  %3.5f \n\n\n',maxdif)

