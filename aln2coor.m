function [coor,dist]=aln2coor(msamat,qgr) % take character matrix of a sequence alignment and compute numerical coordinates

if (nargin<2) ; qgr=0; end % whether to use grantham metric


% factor definitions from Atchley et al. 2005

if (~qgr)
% alphabet
 abet='ACDEFGHIKLMNPQRSTVWY'; % natural alphabet
 f1=[ -0.591 -1.343 1.050 1.357 -1.006 -0.384  0.336 -1.239 1.831 -1.019 -0.663 0.945 0.189 0.931 1.538 -0.228 -0.032 -1.337 -0.595 0.260 ];
 f2=[ -1.302  0.465 0.302 -1.453 -0.590 1.652 -0.417 -0.547 -0.561 -0.987 -1.524 0.828 2.081 -0.179 -0.055 1.399 0.326 -0.279 0.009 0.830 ];
 f3=[ -0.733 -0.862 -3.656 1.477 1.891 1.330 -1.673 2.131 0.533 -1.505 2.219 1.299 -1.628 -3.005 1.502 -4.760 2.213 -0.544 0.672 3.097 ];
 f4=[ 1.570 -1.020 -0.259 0.113 -0.397 1.045 -1.474 0.393 -0.277 1.266 -1.005 -0.169 0.421 -0.503 0.440 0.670 0.908 1.242 -2.128 0.838 ];
 f5=[ -0.146 -0.255 -3.242 -0.837 0.412 2.064 -0.078 0.816 1.648 -0.912 1.212 0.933 -1.392 -1.853 2.897 -2.647 1.313 -1.262 -0.184 1.512 ];
 F=[f1; f2; f3; f4; f5 ]';
else
% alphabet : easier to maintain order in the grantham 74 paper
 abet='SRLPTAVGIFYCHQNKDEMW'; % natural alphabet
 a=1.833 ; b=0.1018 ; g=0.000399 ;
 f1=sqrt(a)*[1.42 0.65 0 0.39 0.71 0 0 0.74 0 0 0.2 2.75 0.58 0.89 1.33 0.33 1.38 0.92 0 0.13];
 f2=sqrt(b)*[ 9.2 10.5 4.9 8.0 8.6 8.1 5.9 9.0 5.2 5.2 6.2 5.5 10.4 10.5 11.6 11.3 13.0 12.3 5.7 5.4 ];
 f3=sqrt(g)*[ 32 124 111 32.5 61 31 84 3 111 132 136 55 96 85 56 119 54 83 105 170];
 F=[f1; f2; f3]';
end

% function to map letters to factors
fall=@(x) F(find(abet==x,1),:);

% extend alphabet :
aext='BJZ-';
abet=[abet aext];
% B
F(find(abet=='B'),:) = 0.5 * ( fall('D') + fall('N') ) ; 
% J
F(find(abet=='J'),:) = 0.5 * ( fall('I') + fall('L') ) ; 
% Z
F(find(abet=='Z'),:) = 0.5 * ( fall('E') + fall('Q') ) ; 
% - (blank/deletion)
F(find(abet=='-'),:) = 0 ; % not clear what to do with missing letters (i.e. insertions/deletions)
%
% redefine fi's :
f1=F(:,1) ;
f2=F(:,2) ;
f3=F(:,3) ;
if ~qgr
 f4=F(:,4) ;
 f5=F(:,5) ;
end

msanum=zeros(size(msamat));
for i = 1:numel(abet)
 msanum(find(msamat==abet(i)))=i;
end

msaf1=f1(msanum);
msaf2=f2(msanum);
msaf3=f3(msanum);
if ~qgr
 msaf4=f4(msanum);
 msaf5=f5(msanum);
end

nseq=size(msamat,1) ; % one sequence per row
% concatenate alignment variables, so that the "coordinates" appear together :
if (qgr)
 coor=reshape( [ msaf1 ; msaf2 ; msaf3 ], nseq, [] );
else
 coor=reshape( [ msaf1 ; msaf2 ; msaf3 ; msaf4 ; msaf5 ], nseq, [] );
end

% decide whether to compute pairwise distances :
if (nargin > 1) % yes
 [nseq,nres]=size(msamat);
 nvar=size(coor,1)/nres;
 dist=zeros(1,nseq*(nseq-1)/2);
 ind=1;
 for i=1:nseq-1
  di=nseq-i ;
  dist(ind:ind+di-1)=sqrt(sum(bsxfun(@minus,coor(i+1:end,:),coor(i,:)).^2,2));
  ind=ind+di;
 end
end

