
% Example for using the function replacewinbyACG and replacewinbyCG.
% Note that replacewinbyACG performs much faster than replacewinbyGC as
% Computation of the Confined Gaussian window requires finding the
% lowest eigenvector of a large matrix.
% Compute Hann-window
% w = hamming(128);
w = blackmanharris(128); 
print('Hamming')
% w = exp(zeros(128,1));
% w = tukeywin(128,0.5);
% Compute Confined Gaussian Window (CGW) and Approximate Confined Gaussian
% window (ACGW) with same RMS temporal width as Hann-window w.
[wACG,alf] = replacewinbyACG(w);
alf
[wCG, alpha] =  replacewinbyCG(w);
alpha
% Display all three windows for comparison
wvtool(wCG, wACG,w)
