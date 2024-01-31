%
qoct=exist('OCTAVE_VERSION');
if (qoct)
  graphics_toolkit('gnuplot');
end
%
qgr=1; % grantham

if ~exist('qread')
 qread=1 ; 
end
if (qread)
 load ./ka-xpr.mat
 qread=0; % do not reread upon reexecution
end

qcorrect=1
% train/test/split :
rseed=pi;
qrnd=1; rng(rseed) ;% set seed for reproducibility
testf=0.50 ;% test fraction (relevant if qrnd=1)
nrep=1 ; % number of samples for computing mean + std

% select data set
dataset='xpr' ; % all strains ; a problem for this could be the different backgrounds
eval(['xp=',dataset,';']);

% 7.26.22 : we did not get rid of stop codons, they are marked as stars in the mutations ; do it here :
stop=cellfun(@(x) any(x=='*'), xp.muts);

tr_selection='xp.mfi(:)>0';
%tr_selection='xp.nmut(:)>0 & xp.mfi(:)>0';
%tr_selection='xp.nmut(:)==1 & xp.mfi(:)>0';
%tr_selection='xp.nmut(:)==2 & xp.mfi(:)>0'; %
%tr_selection='xp.nmut(:)>=3 & xp.mfi(:)>0'; %
%te_selection='xp.nmut(:)==2 & xp.mfi(:)>0'; %
%te_selection='xp.nmut(:)==1 & xp.mfi(:)>0'; %
%te_selection='xp.nmut(:)>=3 & xp.mfi(:)>0'; %
%te_selection='xp.nmut(:)>=2 & xp.mfi(:)>0'; %
te_selection=tr_selection;
%
error_selection='xp.nxp(:)>1';% selection to use to estimate measurement error
%
 cptr=[];
 cstr=[];
 rmstr=[];
 cpte=[];
 cste=[];
 rmste=[];
%
for irep=1:1+(nrep-1)*qrnd

if (qrnd)
 selection=tr_selection;
 te_selection=tr_selection;
 isel = eval ( ['find(~stop & ',selection,')']) ; % evaluate selection specified as text above
 itrain=isel(randperm( numel(isel), round( (min(1,1-testf)*numel(isel)) ))) ;
 if (~exist('rseed'))
  qtest=1 ; qtrain=1 ; qinv=1 ; %need this, otherwise A matrices are incorrect for this (random) strain choice
 end
 itest=setdiff(isel,itrain) ;
% make sure test/train split makes sense
 assert( isempty(intersect(itrain,itest)) )
 assert( isempty(setdiff(union(itrain,itest),isel)) )
else
 itrain=eval ( ['find(~stop & ',tr_selection,')']) ; ;
 itest=eval ( ['find(~stop & ',te_selection,')']) ; ;
end
%
ierror=eval ( ['find(~stop & ',error_selection,')']) ; 

% convert all seqs to coords
ifirst=331;
ilast=531;
for str='wabe'
 eval(['x',str,'seq=aln2coor(xpr.',str,'seq(ifirst:ilast),qgr);']);
%xwseq=aln2coor(kds.wseq(331:531),qgr);
end
nd=numel(xwseq)/numel(xpr.wseq(ifirst:ilast));

% whether to recompute pseudo inverse
if ~exist('qinv')
 qinv=1 ;
end
%
if ~exist('qtrain')
 qtrain=1 ;
end
%

if (qtrain | qinv) % whether to compute training matrix or inverse (which requires training matrix)

% create coordinate matrix
 if (qtrain)
  fprintf('Computing coordinate matrix Atrain\n');
  Atr=mkAmat(xp,xwseq,xaseq,xbseq,xeseq,itrain,ifirst,qgr) ; % rewrote as a function
 end
 if ( qinv)
  fprintf('Computing weights w = Api mfi via pseudo-inverse of Atr\n');
%  Atri=pinv(Atr);
% compute weights
%  w=Atri * xp.mfi(itrain)(:);
%  w = (Atr'*Atr)\(Atr'*xp.mfi(itrain)(:)) ; % much faster to bypass Atri ; probably because inverse computation is very expensive
%  w = Atr\xp.mfi(itrain)(:) ; % same thing as above -- even simpler notation
  w = Atr\xp.mfi(itrain)';
  qinv=0; % to avoid recompute
 end
%
% qtrain=0;
end

% now we can use these weights to apply the model to various sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST
qsplit=1 - (numel(itest)==numel(itrain) && all(sort(itest)==sort(itrain)) ) ;

if (irep==1)
 wgts=w(:)';
else
 wgts=[wgts;w(:)'];
end

mfitrain=Atr*w;

f=figure(1);set(f,'position',[100 100 1000 400]); clf;
if(qsplit) ; subplot(1,2,1) ; end ;
scatter(xp.mfi(itrain),mfitrain);hold on
xlim([3 12])
ylim([3 12])
plot( get(gca,'xlim'), get(gca,'ylim'), 'k--' )
axis equal
xlabel('experiment (log10[MFI])')
ylabel('model')
cp=corr(xp.mfi(itrain)',mfitrain)
if (qoct)
 cs=spearman(xp.mfi(itrain),mfitrain)
else
 cs=corr(xp.mfi(itrain)',mfitrain,'type','spearman')
end
%rmse= sqrt( mean( (mfitrain(:)-xp.mfi(itrain)(:)).^2 ))
rmse= sqrt( mean( (mfitrain(:)-xp.mfi(itrain)').^2 ))
% approximate chi^2 & Q value :
%1 - approximate average error (var) in experimental data:
e2av=mean(xp.mfie(ierror).^2);
ntrain=numel(itrain);
npar=numel(xwseq);
chi2=(rmse.^2)*ntrain / e2av ;
Q=gammainc(0.5*chi2, 0.5*(ntrain-npar)) % typically get near 1 -- which is "perfect"

cptr=[cptr cp];
cstr=[cstr cs];
rmstr=[rmstr rmse];

title({['Training set (',num2str(100*(1-testf*qrnd)),'%)'], ['Selection: ',tr_selection]})
text( 8.5,5  , [ 'c_P=',num2str( cp ) ])
text( 8.5,4.5, [ 'c_S=',num2str( cs ) ])
text( 8.5,4  , [ 'RMSE=',num2str( rmse ) ])
text(1.5,12.5,'A','fontsize', 16);

 set(gca, 'fontsize', 13);

if (qcorrect)
if (qtrain)
% now compute "enrichment" stats
% (1) for a small total sample, form a pair matrix of signal differnces
mfiexp=xp.mfi(itrain);
dmfiexp=bsxfun(@minus,mfiexp(:),mfiexp(:)');
mfisim=mfitrain; % recall that rmsa correlated negatively with the expression !
dmfisim=bsxfun(@minus,mfisim(:),mfisim(:)');
okmat=sign(dmfiexp.*dmfisim);
nmut=numel(itrain) ;
ok = 0.5*sum(okmat(:)>0) / (nchoosek(nmut,2)) % halve because we are overcounting, i.e. the matrix is antisymmetric, reflecting which mut we chose first
notok= 0.5*sum(okmat(:)<0) / (nchoosek(nmut,2)) % should be 1-ok
end
%
text( 8.5, 3.5, ['Correct(%)=',sprintf('%.1f',ok * 100 ) ])
%
end
qtrain=0 ;
%
box on;
%
if (qsplit)
 fprintf('Computing coordinate matrix Atest\n');
 if (~exist('qtest')) ; qtest=1 ; end
 if(qtest)
  Atst=mkAmat(xp,xwseq,xaseq,xbseq,xeseq,itest,ifirst,qgr) ;
%  qtest=0 ;
 end
 mfitest=Atst*w;
 subplot(1,2,2) ;
 scatter(xp.mfi(itest),mfitest,'r');hold on;
 xlim([3 12])
 ylim([3 12])
 plot( get(gca,'xlim'), get(gca,'ylim'), 'k--' )
 axis equal
 xlabel('experiment (log10[MFI])')
 ylabel('model')
 [cp,pp]=corr(xp.mfi(itest)',mfitest)
 if (qoct)
  [cs,ps]=spearman(xp.mfi(itest),mfitest)
 else
  [cs,ps]=corr(xp.mfi(itest)',mfitest,'type','spearman')
 end
% rmse=sqrt( mean( (mfitest(:)-xp.mfi(itest)(:)).^2 ))
 rmse=sqrt( mean( (mfitest(:)-xp.mfi(itest)').^2 ))
%
 cpte=[cpte cp];
 cste=[cste cs];
 rmste=[rmste rmse];
%
 title({['Test set (', num2str(100*(1-(1-testf)*qrnd)),'%)'],['Selection: ',te_selection]})
 text( 8.5,5  , [ 'c_P=',num2str( cp ) ]);
 text( 8.5,4.5, [ 'c_S=',num2str( cs ) ]);
 text( 8.5,4  , [ 'RMSE=',num2str( rmse ) ])
 text(1.5,12.5,'B','fontsize', 16);

 set(gca, 'fontsize', 13);

 if (qcorrect)
% now compute "enrichment" stats
% (1) for a small total sample, form a pair matrix of signal differences
 if (qtest)
 mfiexp=xp.mfi(itest);
 dmfiexp=bsxfun(@minus,mfiexp(:),mfiexp(:)');
 mfisim=mfitest; % recall that rmsa correlated negatively with the expression !
 dmfisim=bsxfun(@minus,mfisim(:),mfisim(:)');
 okmat=sign(dmfiexp.*dmfisim);
 nmut=numel(itest) ;
 ok = 0.5*sum(okmat(:)>0) / (nchoosek(nmut,2)) % halve because we are overcounting, i.e. the matrix is antisymmetric, reflecting which mut we chose first
     notok= 0.5*sum(okmat(:)<0) / (nchoosek(nmut,2)) % should be 1-ok
%
 end
 text( 8.5, 3.5, ['Correct(%)=',sprintf('%.1f',ok * 100 ) ])
 end
 box on ;
 qtest=0;

end

end

wm=mean(wgts,1);

figure(2);
wd=reshape(wm,nd,[]);
if (qoct)
 wn=norm(wd,'cols');
else
 wn = sqrt(sum (wd.^2, 1)) ;
end
plot(ifirst:ilast,wn);

figure(1) ;
if (strcmp(tr_selection, te_selection))
 fname=[dataset,'-',tr_selection];
else
 fname=[dataset,'-',tr_selection,'_',te_selection];
end
%set(gcf, 'paperpositionmode','auto')
%return
%print(gcf, '-depsc2', [fname,'.eps']);
%print(gcf, '-dtiff', [fname,'.tif']);

