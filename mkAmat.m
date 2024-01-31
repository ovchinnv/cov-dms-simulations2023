% create matrix for linear noninteracting per-aa model
function Atr=mkAmat(kd,xwseq,xaseq,xbseq,xeseq,inds,offset,qgrantham)
 Atr=zeros(numel(inds),numel(xwseq));
 if(qgrantham)
  ndim=3;
 else
  ndim=5;
 end % will not worry about other encodings, yet
%
 for i=1:size(Atr,1) ;
%
  if (mod(i,100)==0) ; fprintf('processing %d of %d\n',i,numel(inds)) ; end;
%
  ii=inds(i);
  if(kd.strain(ii)=='W')
   Atr(i,:)=xwseq;
  elseif(kd.strain(ii)=='A')
   Atr(i,:)=xaseq;
  elseif(kd.strain(ii)=='B')
   Atr(i,:)=xbseq;
  elseif(kd.strain(ii)=='E')
   Atr(i,:)=xeseq;
  end
% apply mutation
  for imut=1:kd.nmut(ii)
   ipos=kd.imuts{ii}(imut)-offset+1; % location of mutation
   xmut=aln2coor(kd.muts{ii}(imut),qgrantham); % coordinate of mutant residue
   Atr(i,1+ndim*(ipos-1):ndim*ipos)=xmut;
  end
 end
end

