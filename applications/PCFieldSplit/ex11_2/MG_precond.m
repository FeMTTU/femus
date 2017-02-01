clear all
clc
load 'preconditioner.txt'
pre = preconditioner;
pre = pre';
% pre = inv(pre);
eigrealP = real(eig(pre));
eigimagP = imag(eig(pre));
figure
plot( eigrealP , eigimagP, 'k+')

figure
spy(pre)
return
load 'matrix.txt'
aij = spconvert(matrix);
aij = full(aij);
realJ = real(eig(pre*aij));
imagJ = imag(eig(pre*aij));
figure
plot(realJ, imagJ, 'b*')

% figure
% spy(aij)
 
 
 
% % sizeT = 289;
% % sizeU = 289;
% % sizeV = 289;
% % sizeP = 192;
% % total_size = sizeT + sizeU + sizeV +sizeP;
% % Prematrix = zeros(total_size, total_size);
% % for j =1:total_size
% % 
% %     Umatrix = zeros(sizeU,1);
% %     for j2 = 1:sizeU
% %        Umatrix(j2,1) = pre(j,j2);
% %     end
% %     
% %     Vmatrix = zeros(sizeV,1);
% %     for j3 = 1:sizeV
% %        Vmatrix(j3,1) = pre(j,j3);
% %     end
% %     
% %     Pmatrix = zeros(sizeP,1);
% %     for j4 = 1:sizeP
% %        Pmatrix(j4,1) = pre(j,j4);
% %     end
% %     
% %         Pmatrix = zeros(sizeP,1);
% %     for j4 = 1:sizeP
% %        Pmatrix(j4,1) = pre(j,j4);
% %     end
% %     
% %     Prematrix(1:sizeU,j)= Umatrix;
% %     Prematrix(1:sizeV,j)= Vmatrix;
% %     Prematrix(1:sizeP,j)= Pmatrix;
% %     Prematrix(1:sizeT,j)= Tmatrix;
% % end
% 
% load 'matrix.txt'
% 
%  aij = spconvert(matrix);
% spy(aij)
%  aij = full(aij);
 
%  
% %  plot(eigrealA,eigimagA,'*')
% % hold on
% % bij = zeros(total_size, total_size);
% % for j =1:total_size
% %     Tmatrix = zeros(sizeT,1);
% %     for j1 = 1:sizeT
% %        Tmatrix(j1,1) = aij(j,j1);
% %     end
% %     
% %     Umatrix = zeros(sizeU,1);
% %     for j2 = 1:sizeU
% %        Umatrix(j2,1) = aij(j,j2);
% %     end
% %     
% %     Vmatrix = zeros(sizeV,1);
% %     for j3 = 1:sizeV
% %        Vmatrix(j3,1) = aij(j,j3);
% %     end
% %     
% %     Pmatrix = zeros(sizeP,1);
% %     for j4 = 1:sizeP
% %        Pmatrix(j4,1) = aij(j,j4);
% %     end
% %     
% %     bij(1:sizeU,j)= Umatrix;
% %     bij(1:sizeV,j)= Vmatrix;
% %     bij(1:sizeP,j)= Pmatrix;
% %     bij(1:sizeT,j)= Tmatrix;
% % end
% % 
% % % realJ = real(eig(bij*Prematrix));
% % % imagJ = imag(eig(bij*Prematrix));
% % % plot(realJ, imagJ, '*')
% % 
% % realK = real(eig(bij));
% % imagK = imag(eig(bij));
% % plot(realK, imagK, 'b*')
