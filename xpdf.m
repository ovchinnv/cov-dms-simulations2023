%
qgr=1; % grantham

if ~exist('qread')
 qread=1 ; 
end
if (qread)
 load ./ka-xpr.mat
 qread=0; % do not reread upon reexecution
end

qscatter=0 ; % do not plot scatter plot
% train/test/split :
qrnd=1;
rseed=pi;
qrnd=1; rng(rseed) ;% set seed for reproducibility
testf=0.5 ;% test fraction (relevant if qrnd=1)
nrep=1 ; % number of samples for computing mean + std

% select data set
xp=xpr ; % all strains
% 7.26.22 : we did not get rid of stop codons, they are marked as stars in the mutations ; do it here :
stop=cellfun(@(x) any(x=='*'), xp.muts);

tr_selection='xp.mfi(:)>7.5';
tr_selection='xp.mfi(:)<7.5';
%tr_selection='xp.mfi(:)>0';
tr_selection='xp.nmut(:)>0 & xp.mfi(:)>0';
%tr_selection='xp.nmut(:)<=1 & xp.mfi(:)>0';
%tr_selection='xp.nmut(:)<=1 & xp.mfi(:)>7.5 & xp.mfie(:)>0.25' ; % low exp std (?) not a good test b/c some entries have one meas.
%tr_selection='xp.nmut(:)<=1 & xp.mfi(:)>7.5 & xp.nxp(:)>2' ; % require at least two meas ; no sig diff

%te_selection='xp.nmut(:)<=1 & xp.mfi(:)>7.5'; % low xp
%tr_selection='xp.nmut(:)>1 & xp.mfi(:)>7.5'; % multiple mutants

%te_selection='xp.nmut(:)<=1'; % low xp

%xp=kwuhan ; % wuhan only
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
 itrain=eval ( ['find(~stop & ',tr_selection,')']) ;
 if (~exist('te_selection')); te_selection=tr_selection;end
 itest=eval ( ['find(~stop & ',te_selection,')']) ;
end

% convert all seqs to coords
ifirst=331;
ilast=531;
for str='wabe'
 eval(['x',str,'seq=aln2coor(xpr.',str,'seq(ifirst:ilast),qgr);']);
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
  w=Atri * xp.mfi(itrain(:))';
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
else
 wgts=[wgts;w(:)'];
end

mfitrain=Atr*w;

if (qscatter)

figure(1);clf;
if(qsplit) ; subplot(1,2,1) ; end ;
scatter(xp.mfi(itrain),mfitrain);hold on
xlim([4 11])
ylim([4 11])
plot( get(gca,'xlim'), get(gca,'ylim'), 'k--' )
axis equal
xlabel('experiment')
ylabel('model')
cp=corr(xp.mfi(itrain),mfitrain)
cs=spearman(xp.mfi(itrain),mfitrain)
rmse=log(10)*0.6 * sqrt( mean( (mfitrain(:)-xp.mfi(itrain(:)')).^2 ))

cptr=[cptr cp];
cstr=[cstr cs];
rmstr=[rmstr rmse];

title({['Training set (',num2str(100*(1-testf*qrnd)),'%)'], ['Selection: ',tr_selection]})
text( -3+min(xp.mfi(itrain)), max(xp.mfi(itrain))-1, [ 'c_P=',num2str( cp ) ])
text( -3+min(xp.mfi(itrain)), max(xp.mfi(itrain))-1.25, [ 'c_S=',num2str( cs ) ])
text( -3+min(xp.mfi(itrain)), max(xp.mfi(itrain))-1.5, [ 'RMSE(kcal/mol)=',num2str( rmse ) ])
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
 xlim([4 11])
 ylim([4 11])
 plot( get(gca,'xlim'), get(gca,'ylim'), 'k--' )
 axis equal
 xlabel('experiment')
 ylabel('model')
 cp=corr(xp.mfi(itest),mfitest)
 cs=spearman(xp.mfi(itest),mfitest)
 rmse=log(10)*0.6 * sqrt( mean( (mfitest(:)-xp.mfi(itest(:)')).^2 ))
%
 cpte=[cpte cp];
 cste=[cste cs];
 rmste=[rmste rmse];
%
 title({['Test set (', num2str(100*testf),'%)'],['Selection: ',te_selection]})
 text( -3+min(xp.mfi(itest)), max(xp.mfi(itest))-1, [ 'c_P=',num2str( cp ) ])
 text( -3+min(xp.mfi(itest)), max(xp.mfi(itest))-1.25, [ 'c_S=',num2str( cs ) ])
 text( -3+min(xp.mfi(itest)), max(xp.mfi(itest))-1.5, [ 'RMSE(kcal/mol)=',num2str( rmse ) ])
end

end % qsplit
end % qscatter

wm=mean(wgts,1);

figure(1); % plot pdf here ?
subplot(1,2,1);

xdata=xp.mfi(itrain);
ydata=mfitrain;
xmin=min(xdata);
xmax=max(xdata);
ymin=min(ydata);
ymax=max(ydata);
%ymin=xmin;
%ymax=xmax;
%
dx=0.2 ;
dy=0.2 ;
%
xbin=[xmin:dx:xmax]; % actually, consider these left bin edges
ybin=[ymin:dy:ymax];

xinds = fix ( ( xdata-xmin)/dx ) + 1 ;
yinds = fix ( ( ydata-ymin)/dy ) + 1 ;

pdf = zeros(numel(xbin), numel(ybin) );
for kind=1:numel(xdata)
 pdf(xinds(kind), yinds(kind))=pdf(xinds(kind), yinds(kind)) + 1 ;
end

%pcolor(ybin-0.5*dy, xbin-0.5*dy, pdf) ; shading flat ; 
%mesh(ybin-0.5*dy, xbin-0.5*dy, pdf) ; shading flat ; 
%c=contourf(xbin-0.5*dx, ybin-0.5*dy, pdf', 10) ; shading flat ; 
[cl,ch]=contour(xbin-0.5*dx, ybin-0.5*dy, pdf', 10, 'k') ; % shading flat ; 
clev=ch.LevelList
%shading interp ;
%colorbar ;
%axis equal ; colormap jet ;
%xlim([5 11])
%ylim([5 11])
set(gca, 'tickdir','out') ;  box on ; 

fname=['xppdf-',tr_selection] ;
set(gcf, 'paperpositionmode','auto')
print(gcf, '-depsc2', [fname,'.eps']);
print(gcf, '-dtiff', [fname,'.tif']);

