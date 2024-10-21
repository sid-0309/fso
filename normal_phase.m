function fn = normal_phase(x,y,mu,sigma)
    [X,Y] = meshgrid(x,y);
    Z = [X(:) Y(:)];
    z = mvnpdf(Z,mu,sigma);
    z = reshape(z,length(x),length(y));
    fn = z;
end