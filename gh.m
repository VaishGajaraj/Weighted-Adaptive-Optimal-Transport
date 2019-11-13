function [grad ,hess ] = gh( a,b,x,mu,y,nu,sigma )
[~,xlength]=size(x);
[~,ylength]=size(y);
[alength,~]=size(a);
[blength,~]=size(b);
AAL=zeros(alength); ABL=zeros(alength,blength); BBL=zeros(blength);
AL=zeros(alength,1); BL=zeros(blength,1);
for j = 1:xlength
    [ XF,AFx,AFy,AAFx,AAFy ] = f( x(:,j),a,sigma );
    [ XG,XGx,XGy,BG,BXG,BBG] = g( XF,b,sigma );
    AL = AL + mu(j)*[AFx AFy]*XG;
    BL = BL + mu(j)*BG;
    AAL = AAL + mu(j)*([AFx AFy]*XGx*AFx' + [AFx AFy]*XGy*AFy' + XG(1)*AAFx + XG(2)*AAFy);
    ABL = ABL +  mu(j)*([AFx AFy]*BXG);
    BBL = BBL  + mu(j)*(BBG);
end
for j = 1:ylength
    [ G,BG,BBG ] = g1( y(:,j),b,sigma );
    BL = BL - nu(j)*exp(G)*BG;
    BBL = BBL - nu(j)*(exp(G)*(BG*BG' + BBG));
end
grad = [AL; BL];
hess = [AAL ABL; ABL' BBL];
end