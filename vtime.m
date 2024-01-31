% some tests of RMSD's computed from VMD
 qoct=exist('OCTAVE_VERSION');
 if (qoct)
  graphics_toolkit('gnuplot');
 end
%
% qcg2=0 ; % to see data for the subset coresponding to the second coarse-graining
%
 addpath('~/scripts/matlab/anal');
% 1 -- read train data:
 selection='xp.nmut(:)>0 & xp.mfi(:)>0';
 dcut=0.01;
 fname=['./xprcg-',num2str(dcut),selection,'-vrms-rmsf.mat'];
%
 load(fname);
%
 if (qcg2)
   icgind=load('icgind2.dat'); % indices for second coarse-grained dataset
 end
%
 sel='nohr336_526' ; % select residues of the RBD
%
% the parameters below specify what to plot ; they can be set here, or passed in from the outside
% mdsamples=1:3 ; tag='ures' ; % unrestrained simulations ... so-so answer
% mdsamples=4:6 ; tag='res' ; % restrained simulations ...  better answer
% if (qcg2)
%  mdsamples=7:9 ; tag='long' ; % unrestrained long simulations, cg2 subset, only for qcg2
%  mdsamples=10:12 ; tag='cube' ; % unrestrained cube simulations, cg2 subset, only for qcg2
% end
%
 cmd=['xp.rms_',sel,'(mdsamples,icgind);'] ;
 rmsd=eval(cmd) ;
%
 qfast=0; % assume a data layout and reshape the data all at once:
%
 nrms=numel(rmsd{1,1});
 [nrep,nmut]=size(rmsd);
%
 if (qfast) % really, not so much faster !!!
% checked to make sure it's the same operation as the slow method below
  rmsall = reshape(cell2mat(rmsd),nrep,nrms,nmut);
 else
  rmsall=zeros(nrep,nrms,nmut) ;
  for i=1:numel(icgind)
   for j=1:nrep
    rmsall(j,:,i)=rmsd{j,i};
   end
  end
 end
%
 rmsa=zeros(nrep+1,nmut);
%
 irms=0 + 1;
 i=1;
 tinds = irms+1 : nrms ;
 clear crms cmfi srms smfi ;
 for erms = tinds
  rmsa(1:nrep,:)=squeeze(mean(rmsall(:,irms:erms,:),2)) ; % average rmsd
  rmsa(nrep+1,:)=mean(rmsa(1:nrep,:),1) ; % average in last position
  c=corr(rmsa(1:nrep,:)'); % generally low correlation
  crms(i) = ( mean(c(:))*nrep - 1 ) / (nrep-1) ; % average coeff
  cmfi(i) = corr(xp.mfi(icgind)', rmsa(nrep+1,:)');

  if (qoct)
   c=spearman(rmsa(1:nrep,:)'); % generally low correlation
   srms(i) = ( mean(c(:))*nrep - 1 ) / (nrep-1) ; % average coeff
   smfi(i) = spearman(xp.mfi(icgind), rmsa(nrep+1,:));
  else
   c=corr(rmsa(1:nrep,:)','type','spearman'); % generally low correlation
   srms(i) = ( mean(c(:))*nrep - 1 ) / (nrep-1) ; % average coeff
   smfi(i) = corr(xp.mfi(icgind)', rmsa(nrep+1,:)','type','spearman');
  end
  i=i+1 ;
 end
 t = tinds * 100000 * 4 / 1e6 ;
%
 f3 = figure(3) ; clf ; hold on ; box on ;
%
 step=round(length(t)/75) ;
%
 plot(t(1:step:end),crms(1:step:end),'ko-')
 plot(t(1:step:end),cmfi(1:step:end),'rs-')
 plot(t(1:step:end),srms(1:step:end),'bv-')
 plot(t(1:step:end),smfi(1:step:end),'g*-')
%
 xlim([0,t(end)]);
 ylim([-0.6 0.8]);
%
 xlabel('\it t(ns)') ; 
% ylabel('\it correlation') ;
 ylabel('\it c_P, c_S') ;
 set(gca, 'fontsize' , 16 ) ;
%
 text(-0.2*t(end), 0.85, label, 'fontsize', 18)
%
% print(gcf, '-depsc2', 'corr-vs-time.eps')
 print(gcf, '-depsc2', [sel,'-corr-vs-time-',tag,'.eps'])
