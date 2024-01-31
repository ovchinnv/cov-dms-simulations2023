% display model results
qoct=exist('OCTAVE_VERSION');

cut='>0';

testfs=[0.1:0.1:0.9] ;
trainfs=1-testfs;

ind=1;
dataset='xpr'; % everything
if ~exist('flg') ; flg='';end
for testf=testfs
 fname = [dataset,cut,'-tefr', num2str(testf),flg];
 load([fname,'.mat']);

 cptrm(ind)=mean(cptr);
 cptem(ind)=mean(cpte);
 rmstrm(ind)=mean(rmstr);
 rmstem(ind)=mean(rmste);
 cstrm(ind)=mean(cstr);
 cstem(ind)=mean(cste);
 oktrm(ind)=mean(oktr);
 oktem(ind)=mean(okte);

 cptre(ind)=std(cptr);
 cptee(ind)=std(cpte);
 rmstre(ind)=std(rmstr);
 rmstee(ind)=std(rmste);
 cstre(ind)=std(cstr);
 cstee(ind)=std(cste);
 oktre(ind)=std(oktr);
 oktee(ind)=std(okte);

 ind=ind+1;
end

ave=[cptrm;cptem;cstrm;cstem;];
err=2*[cptre;cptee;cstre;cstee;];

ave=[cptrm;cptem;cstrm;cstem;rmstrm;rmstem];
err=2*[cptre;cptee;cstre;cstee;rmstre;rmstee];

ave=[cptrm;cptem;cstrm;cstem;rmstrm;rmstem;oktrm;oktem];
err=2*[cptre;cptee;cstre;cstee;rmstre;rmstee;oktre;oktee];
%bar(trainfs,[cptrm;cptem]')
f=figure(1);
set(f,'position', [100,200, 1000, 400]);
clf
bar(trainfs,ave') ; hold on ; 
%bar(testfs,cptrm)
%
% add error bars :
nbar=size(ave,1);
for i=1:nbar
 if (qoct)
  offset = - ( -0.5*(nbar-1) + i - 1 ) ; % need to invert because plotted from right to left !
  bw=0.072/nbar;
 else
  offset =   ( -0.5*(nbar-1) + i - 1 ) ; % matlab
  bw=0.08/nbar;
 end
 offset*bw
 errorbar(trainfs + offset * bw ,ave(i,:),0*err(i,:),err(i,:),'k.')
end

set(gca, 'xtick', testfs, 'tickdir','out', 'xticklabel',testfs*100)
set(gca, 'fontsize',15)
%set(gca,'yscale','log')
xlabel('% data used for training')
l=legend({'c_P^{TRAIN}','c_P^{TEST}','c_S^{TRAIN}','c_S^{TEST}', 'RMSE^{TRAIN}', 'RMSE^{TEST}','q_{correct}^{TRAIN}','q_{correct}^{TEST}'},'location','northeastoutside')
legend boxoff;
%set(l,'orientation','horizontal')
set(l,'orientation','vertical')
set(l, 'fontsize',14)

ylim([0 2.5])
text(-0.075, 2.7, 'C','fontsize', 15);

fname=[dataset,cut,flg];

set(gcf, 'paperpositionmode','auto')
print(gcf, '-depsc2', [fname,'.eps']);
