qoct=exist('OCTAVE_VERSION');
if (qoct)
  graphics_toolkit('gnuplot');
end
%
qgr=1; % grantham encoding

if ~exist('qread')
 qread=1 ;
end
if (qread)
 load ./ka-xpr.mat
 qread=0; % do not reread upon reexecution
end

% select data set
dataset='xpr' ; % everything
eval(['xp=',dataset,';']);

% 7.26.22 : we did not get rid of stop codons, they are marked as stars in the mutations ; do it here :
stop=cellfun(@(x) any(x=='*'), xp.muts);

qsave=1 ;% whether to save data to a file (possibly overwriting an existing file)
% train/test/split :
qrnd=1;

%testf=0.5 ;% test fraction (relevant if qsplit=1)
for testf=[0.1:0.1:0.9]

testf

nrep=3 ; % number of samples for computing mean + std ; 10 in the paper, reduced here for faster speed

%cut='>7.5';
%cut='<7.5';
cut='>0';

%tr_selection=['xp.mfi(:)',cut]; % all mutations, as in the paper
tr_selection=['xp.nmut(:)<=1 & xp.mfi(:)',cut]; % single-mutant set to run faster

%te_selection='xp.nmut(:)<=1 & xp.mfi(:)>7.5'; % low xp
%tr_selection='xp.nmut(:)>1 & xp.mfi(:)>7.5'; % multiple mutants

%te_selection='xp.nmut(:)<=1'; % low xp

%xp=xwuhan ; % wuhan only
% select indices
%selection=' xp.strain(:)=='''W''' & xp.nmut(:)<=1 '               ;% e.g. wuhan strain with no more than one mutation
%selection=' xp.strain(:)=='''W''' & xp.nmut(:)<=1 & xp.nxp(:)>1 ' ;% e.g. wuhan strain with no more than one mutation & 2 or more meas
%selection=' xp.strain(:)=='''W''' & xp.nmut(:)<=1 & xp.mfi(:)>6 ' ;% e.g. wuhan strain with no more than one mutation & high xp

% all strains
%isel=(1:numel(xpr.nxp)); % lower corr, as expected

%isel=find( xp.nmut(:)<=1 & xp.mfi(:)<7.5 ) ; % low xp


%selection='xp.strain(:)=='''W''' & xp.nmut(:)<=1 & xp.mfi(:)<8 ' ;% e.g. wuhan strain with no more than one mutation & low xp
%selection='xp.strain(:)=='''W''' & xp.nmut(:)<=1 & xp.mfie(:)<0.5 ' ;% e.g. wuhan strain with no more than one mutation & low error
%selection='xp.strain(:)=='''W''' & xp.nmut(:)<=1 & xp.mfie(:)<0.5 & xp.mfi(:)>6 ' ;% you get the idea

%
 cptr=[];
 cstr=[];
 rmstr=[];
 oktr=[];
 cpte=[];
 cste=[];
 rmste=[];
 okte=[];
%
for irep=1:1+(nrep-1)*qrnd

if (qrnd)
 selection=tr_selection;
 te_selection=tr_selection ;
 isel = eval ( ['find(~stop & ',selection,')']) ; % evaluate selection specified as text above
 itrain=isel(randperm( numel(isel), round( (min(1,1-testf)*numel(isel)) ))) ;
 qtest=1 ; qtrain=1 ; qinv=1 ; %need this, otherwise A mtrices are incorrect for this (random) strain choice
 itest=setdiff(isel,itrain) ;
% make sure test/train split makes sense
 assert( isempty(intersect(itrain,itest)) )
 assert( isempty(setdiff(union(itrain,itest),isel)) )
else
 itrain=eval ( ['find(~stop & ',tr_selection,')']) ; ;
 itest=eval ( ['find(~stop & ',te_selection,')']) ; ;
end

% convert all seqs to coords
% RBD sequence limits
ifirst=331;
ilast=531;
for str='wabe'
 eval(['x',str,'seq=aln2coor(xpr.',str,'seq(ifirst:ilast),qgr);']);
%xwseq=aln2coor(xpr.wseq(331:531),qgr);
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
 fprintf('Computing coordinate matrix Atrain\n');
 Atr=mkAmat(xp,xwseq,xaseq,xbseq,xeseq,itrain,ifirst,qgr) ; % rewrote as a function
 fprintf('Computing pseudo-inverse of Atr\n');
 if ( qinv)
  Atri=pinv(Atr);
% compute weights
  fprintf('Computing weights from w = Api mfi \n');
  w=Atri * xp.mfi(itrain)';
  qinv=0; % to avoid recompute
 end
%
 qtrain=0;
end

% now we can use these weights to apply the model to various sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST
qsplit=1 - (numel(itest)==numel(itrain) && all(sort(itest)==sort(itrain)) ) ;

if (irep==1)
 wgts=w(:)';
 itrains=itrain(:)';
else
 wgts=[wgts;w(:)'];
 itrains=[itrains;itrain(:)'];
end

mfitrain=Atr*w;

f=figure(1);set(f,'position',[100 100 1000 400]); clf;
if(qsplit) ; subplot(1,2,1) ; end ;
scatter(xp.mfi(itrain),mfitrain);hold on
xlim([3 12])
ylim([3 12])
plot( get(gca,'xlim'), get(gca,'ylim'), 'k--' )
axis equal
xlabel('experiment')
ylabel('model')
[cp,pp]=corr(xp.mfi(itrain)',mfitrain)
if (qoct)
 [cs,ps]=spearman(xp.mfi(itrain),mfitrain)
else
 [cs,ps]=corr(xp.mfi(itrain)',mfitrain,'type','spearman')
end
rmse=sqrt( mean( (mfitrain(:)-xp.mfi(itrain)').^2 ))

cptr=[cptr cp];
cstr=[cstr cs];
rmstr=[rmstr rmse];

title({['Training set (',num2str(100*(1-testf)),'%)'], ['Selection: ',tr_selection]})
text( 8,6  , [ 'c_P=',num2str( cp ) ])
text( 8,5.5, [ 'c_S=',num2str( cs ) ])
text( 8,5  , [ 'RMSE=',num2str( rmse ) ])
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
oktr=[oktr ok];
%
text( 8, 4.5, ['Correct(%)=',sprintf('%.1f',ok * 100 ) ])
%
box on;
%
if (qsplit)
 fprintf('Computing coordinate matrix Atest\n');
 if (~exist('qtest')) ; qtest=1 ; end
 if(qtest)
  Atst=mkAmat(xp,xwseq,xaseq,xbseq,xeseq,itest,ifirst,qgr) ;
  qtest=0 ;
 end
 mfitest=Atst*w;
 subplot(1,2,2) ;
 scatter(xp.mfi(itest),mfitest,'r');hold on;
 xlim([3 12])
 ylim([3 12])
 plot( get(gca,'xlim'), get(gca,'ylim'), 'k--' )
 axis equal
 xlabel('experiment')
 ylabel('model')
 cp=corr(xp.mfi(itest)',mfitest)
 if (qoct)
  cs=spearman(xp.mfi(itest),mfitest)
 else
  cs=corr(xp.mfi(itest)',mfitest,'type','spearman')
 end
 rmse=sqrt( mean( (mfitest(:)-xp.mfi(itest)').^2 ))
%
 cpte=[cpte cp];
 cste=[cste cs];
 rmste=[rmste rmse];
%
 title({['Test set (', num2str(100*testf),'%)'],['Selection: ',te_selection]})
 text( 8,6  , [ 'c_P=',num2str( cp ) ])
 text( 8,5.5, [ 'c_S=',num2str( cs ) ])
 text( 8,5  , [ 'RMSE=',num2str( rmse ) ])
%
% now compute "enrichment" stats
% (1) for a small total sample, form a pair matrix of signal differnces
 mfiexp=xp.mfi(itest);
 dmfiexp=bsxfun(@minus,mfiexp(:),mfiexp(:)');
 mfisim=mfitest; % recall that rmsa correlated negatively with the expression !
 dmfisim=bsxfun(@minus,mfisim(:),mfisim(:)');
 okmat=sign(dmfiexp.*dmfisim);
 nmut=numel(itest) ;
 ok = 0.5*sum(okmat(:)>0) / (nchoosek(nmut,2)) % halve because we are overcounting, i.e. the matrix is antisymmetric, reflecting which mut we chose first
 notok= 0.5*sum(okmat(:)<0) / (nchoosek(nmut,2)) % should be 1-ok
 okte=[okte ok];
%
 text( 8, 4.5, ['Correct(%)=',sprintf('%.1f',ok * 100 ) ])

 box on ;

end

wm=mean(wgts,1);

fname = [dataset,cut,'-tefr', num2str(testf)];
if (~qgr) ; fname=[fname,'a']; end

if (~exist('qsave'));qsave=1;end
if (qsave)
 save('-mat', [fname,'.mat'], 'testf', 'nrep', 'cptr', 'cstr', 'rmstr', 'cpte', 'cste', 'rmste', 'wgts', 'wm', 'itrains', 'oktr', 'okte')
end

end % testf

end
