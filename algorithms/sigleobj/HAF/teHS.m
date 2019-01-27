[bv,bs]=HS('FDCGA');
% ff = FDCGA;
config.M = bs/sum(bs);
bs = bs/sum(bs);
alGbest( 'HAF1' , 'cec15_f1',bs,config);

[bv,bs]=HS('FD1');
config.M = bs/sum(bs)
bs = cumsum (bs/sum(bs))
alGbest( 'HAF2' , 'cec15_f1',bs,config);

[bv,bs]=HS('FD2');
config.M = bs/sum(bs)
bs = (bs/sum(bs))
alGbest( 'HAF3' , 'cec15_f1',bs,config);
