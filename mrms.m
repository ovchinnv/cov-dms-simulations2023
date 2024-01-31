% RMSF computed using matlab ; add them up to compute overall fluctuation
 qoct=exist('OCTAVE_VERSION');
 if (qoct)
  graphics_toolkit('gnuplot');
 end
%
% qcg2=1 ; % to see data for the subset corresponding to the second coarse-graining ; set here or pass in
%
 qfit=1 ;% whether to plot linear fit
%
% 1 -- read data:
 selection='xp.nmut(:)>0 & xp.mfi(:)>0';
 dcut=0.01;
 fname=['./xprcg-',num2str(dcut),selection,'-vrms-rmsf.mat'];
%
 load(fname);
%
 if (qcg2)
   icgind=load('icgind2.dat'); % indices for second coarse-grained dataset
 end
% residue limits for averaging RMSF curve
%
 ioffset=331 ; % this is the resid of the first residue of the RMD sequence
 ires=337;
 eres=525;
%
 sel='cg' ;
%
% the parameters below specify what to plot ; they can be set here, or passed in from the outside
% mdsamples=1:3 ; tag='ures' ; % unrestrained simulations ...
% mdsamples=4:6 ; tag='res' ; % restrained simulations ...
% if (qcg2)
%  mdsamples=7:9 ;  tag='long' ; % unrestrained long simulations, cg2 subset, only for qcg2
%  mdsamples=10:12 ;tag='cube' ; % unrestrained cube simulations, cg2 subset, only for qcg2
% end
%
 cmd=['xp.rmsf_',sel,'(mdsamples,icgind);'] ;
 rmsf=eval(cmd) ;
%
 qfast=1; % assume a data layout and reshape the data all at once:
%
 nres=numel(rmsf{1,1});
 [nrep,nmut]=size(rmsf);
%
%
 if (qfast) % really, not so much faster !!!
% checked to make sure it's the same operation as the slow method below
  rmsall = reshape(cell2mat(rmsf),nrep,nres,nmut);
 else
  rmsall=zeros(nrep,nres,nmut) ;
  for i=1:numel(icgind)
   for j=1:nrep
    rmsall(j,:,i)=rmsf{j,i};
   end
  end
 end
%
 rmsa=zeros(nrep,nmut);
 rmsa(:,:)=squeeze(mean(rmsall(:,[ires:eres]-ioffset+1,:),2)) ; % average rmsf

 cp=corr(rmsa') % generally low correlation
 if (qoct)
  cs=spearman(rmsa')
 else
  cs=corr(rmsa','type','Spearman') % generally low correlation 
 end
%
 close all ;
 f = figure(1) ; clf ; hold on ; box on ;
%
% put rmsa average at the end ;
 rmsa(nrep+1,:)=mean(rmsa(1:nrep,:),1) ;
 rmsa(nrep+2,:)=max(rmsa(1:nrep,:),[],1) ;
%
 ind=nrep+1; % average pos
% ind=nrep+2 ; %
 scatter ( xp.mfi(icgind), rmsa(ind,:) ) ;
 xlabel('log MFI (experiment)') ;
 ylabel('<RMSF>(Ang)') ;
 set(gca, 'fontsize',16)
%
 [cp,pp]=corr( xp.mfi(icgind(:))', rmsa(ind,:)' )
% cl=corr( xp.mfi(icgind(:))', log(rmsa(ind,:))' )
% ce=corr( exp(xp.mfi(icgind(:)))', rmsa(ind,:)' )
 if (qoct)
  [cs,ps]=spearman( xp.mfi(icgind), rmsa(ind,:) )
% note spearman with log would be the same because log is monotonic and therefore preserves order
 else
  [cs,ps]=corr( xp.mfi(icgind(:))', rmsa(ind,:)' ,'type','Spearman')
 end
%
 if (qfit)
  lfit=polyfit(xp.mfi(icgind), rmsa(ind,:),1) ;
  xx=[min(xp.mfi(icgind)) max(xp.mfi(icgind))];
  plot(xx,lfit(2) + lfit(1)*xx,'k-');
 end
% text( -0+min(xp.mfi(icgind)), 2.5, [ 'c_P=',num2str( cp ) ])
% text( -0+min(xp.mfi(icgind)), 2.6, [ 'c_S=',num2str( cs ) ])
 d= ( max(rmsa(ind,:))-min(rmsa(ind,:)) ) / 10 ;
 text( 9., max(rmsa(ind,:))-0*d, [ 'c_P=',num2str( cp ),'(', sprintf('%.1d',pp),')' ])
 text( 9., max(rmsa(ind,:))-1*d, [ 'c_S=',num2str( cs ),'(', sprintf('%.1d',ps),')' ])
% ylim([0.5 1.5])
% label :
 if ~exist('label') ; label='' ; end % survive in case label is not defined
 text(2.8,max(rmsa(ind,:))+2*d,label,'fontsize',18)
% now compute "enrichment" stats
% (1) for a small total sample, form a pair matrix of signal differnces
 mfiexp=xp.mfi(icgind);
 dmfiexp=bsxfun(@minus,mfiexp',mfiexp);
 mfisim=-rmsa(ind,:); % recall that rmsa correlated negatively with the expression !
 dmfisim=bsxfun(@minus,mfisim',mfisim);
 okmat=sign(dmfiexp.*dmfisim);
 ok = 0.5*sum(okmat(:)>0) / (nchoosek(nmut,2)) % halve because we are overcounting, i.e. the matrix is antisymmetric, reflecting which mut we chose first
 notok= 0.5*sum(okmat(:)<0) / (nchoosek(nmut,2)) % should be 1-ok
%
 text( 9.5, max(rmsa(ind,:))-2*d, ['Correct(%)=',sprintf('%.1f',ok * 100 ) ])
% print fig:
 figure(1);
 set (gcf, 'paperpositionmode','auto');
 print(gcf, '-depsc2', [sel,'-rmsf-mfi-',tag,'.eps'])
%
