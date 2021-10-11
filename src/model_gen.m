function model = model_gen(vs,vp,rho,z,n)

%global model;
model.vsv = [];
model.vpv = [];
model.rhov = [];
model.h = [];

for i =1:length(n);
   model.vsv = [model.vsv vs(i)*ones(1,n(i))];
   model.vpv = [model.vpv vp(i)*ones(1,n(i))];
   model.rhov = [model.rhov rho(i)*ones(1,n(i))];
   model.h = [model.h  z(i)*ones(1,n(i))];
end

model.Nn = length(model.h);
model.hzcum = cumsum(model.h);
