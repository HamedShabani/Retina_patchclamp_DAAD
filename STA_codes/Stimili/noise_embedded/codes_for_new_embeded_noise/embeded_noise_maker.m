load('D:\Hame2\Scripts\enoise_embeded\sin_noise_embeded.mat')

sig=noiseembeded';
sig3=[];
for i=1:length(sig)
    sig2=[sig(i)+400 1000];
    sig3=[sig3;sig2; 0.000000 39000.000000];
end
  
T = table(sig3,  'VariableNames', { 'sig3'} )
% Write data to text file
writetable(T, 'MyFile.txt')
dlmwrite('Myfile.txt',sig3,'roffset',0,'coffset',1,'-append','delimiter',' ','newline','pc','precision','%.6f')