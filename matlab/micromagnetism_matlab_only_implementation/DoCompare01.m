
'' ;
DemagTensorCompare =  InteractionMatricesCompare.DemagTensor ;
figure ;
plot(InteractionMatrices.DemagTensor.KglobYY{1}(:),'.')
hold on
plot(DemagTensorCompare.KglobYY{1}(:),'or')

sqrt(sum((InteractionMatrices.DemagTensor.KglobXX{1}(:)-DemagTensorCompare.KglobXX{1}(:)).^2)/sum(DemagTensorCompare.KglobXX{1}(:).^2))
sqrt(sum((InteractionMatrices.DemagTensor.KglobXY{1}(:)-DemagTensorCompare.KglobXY{1}(:)).^2)/sum(DemagTensorCompare.KglobXX{1}(:).^2))
sqrt(sum((InteractionMatrices.DemagTensor.KglobXZ{1}(:)-DemagTensorCompare.KglobXZ{1}(:)).^2)/sum(DemagTensorCompare.KglobXZ{1}(:).^2))
sqrt(sum((InteractionMatrices.DemagTensor.KglobYY{1}(:)-DemagTensorCompare.KglobYY{1}(:)).^2)/sum(DemagTensorCompare.KglobXX{1}(:).^2))
sqrt(sum((InteractionMatrices.DemagTensor.KglobYZ{1}(:)-DemagTensorCompare.KglobYZ{1}(:)).^2)/sum(DemagTensorCompare.KglobXX{1}(:).^2))
sqrt(sum((InteractionMatrices.DemagTensor.KglobZZ{1}(:)-DemagTensorCompare.KglobZZ{1}(:)).^2)/sum(DemagTensorCompare.KglobXX{1}(:).^2))


%%
max( abs((InteractionMatrices.DemagTensor.KglobXX{1}(:)-DemagTensorCompare.KglobXX{1}(:)) ))  %  ./(DemagTensorCompare.KglobXX{1}(:).^2))
max( abs((InteractionMatrices.DemagTensor.KglobXY{1}(:)-DemagTensorCompare.KglobXY{1}(:)) )) 
max( abs((InteractionMatrices.DemagTensor.KglobXZ{1}(:)-DemagTensorCompare.KglobXZ{1}(:)) )) 
max( abs((InteractionMatrices.DemagTensor.KglobYY{1}(:)-DemagTensorCompare.KglobYY{1}(:)) )) 
max( abs((InteractionMatrices.DemagTensor.KglobYZ{1}(:)-DemagTensorCompare.KglobYZ{1}(:)) )) 
max( abs((InteractionMatrices.DemagTensor.KglobZZ{1}(:)-DemagTensorCompare.KglobZZ{1}(:)) )) 


figure ;
plot(abs((InteractionMatrices.DemagTensor.KglobXX{1}(:)-DemagTensorCompare.KglobXX{1}(:)) ),'.')

figure
surf(log10(abs(InteractionMatrices.DemagTensor.KglobXX{1}-DemagTensorCompare.KglobXX{1})) ,'linestyle','none')

%%


figure ;
plot(InteractionMatrices.DemagTensor.KglobYY{1}(:)-DemagTensorCompare.KglobYY{1}(:),'.')


isequal(InteractionMatricesCompare.X(:),InteractionMatrices.X(:))
isequal(InteractionMatricesCompare.Y(:),InteractionMatrices.Y(:))
isequal(InteractionMatricesCompare.Z(:),InteractionMatrices.Z(:))

sqrt(sum((InteractionMatrices.X(:)-InteractionMatricesCompare.X(:)).^2)/sum(InteractionMatricesCompare.X(:).^2))

sqrt(sum((InteractionMatrices.Y(:)-InteractionMatricesCompare.Y(:)).^2)/sum(InteractionMatricesCompare.Y(:).^2))

sqrt(sum((InteractionMatrices.Z(:)-InteractionMatricesCompare.Z(:)).^2)/sum(InteractionMatricesCompare.Z(:).^2))

% figure ;
% plot(InteractionMatrices.Y(:),'.')
% hold on
% plot(InteractionMatricesCompare.Y(:),'or')

sqrt(sum((InteractionMatrices.A2(:)-InteractionMatricesCompare.A2(:)).^2)/sum(InteractionMatricesCompare.A2(:).^2))


figure ;
plot(InteractionMatrices.A2(:),'.')
hold on
plot(InteractionMatricesCompare.A2(:),'or')

