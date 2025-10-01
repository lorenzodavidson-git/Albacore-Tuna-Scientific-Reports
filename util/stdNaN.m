% function [Mmatrix, Num] = stdNaN(matrix,n)
%
% Take the STD for the index N in MATRIX	, where
% MATRIX can contain NaN values. These will be treated
% as giving zero contribution to the mean.
% OPT is optional. If assigned instead of computing the mean
% it will just compute the sum.
% NUM returns the numbers of elements found in the mean (or sum
% if OPT is defined).

function [Mmatrix, Num] = stdNaN(matrix,n,varargin)

[MyMean, Num] = meanNaN(matrix,n);
isize=size(matrix);
isize2=isize(n);
isize(:)=1; isize(n)=isize2;


%matrix=abs(matrix-repmat(MyMean,isize));
matrix=(matrix-repmat(MyMean,isize)).^2;

Num=ones(size(matrix));

inans=find(isnan(matrix));
matrix(inans)=0;
Num(inans)=0;

Mmatrix=sum(matrix,n);
Num=sum(Num,n);

inans=find(Num == 0);

%   Mmatrix =  Mmatrix./Num;
   Mmatrix =  sqrt(Mmatrix./(Num-1));
   Mmatrix(inans) = NaN;
   Num(inans) = NaN;

