function im_r=Interpolation(im_n,c)

    % Delaunay triangulation based interpolation
    % c=0 denotes preserved pixels
    [x,y]=find(c==0);
    [M,N]=size(im_n);
    [x1,y1]=meshgrid(1:M,1:N);
    im_r=griddata(x,y,im_n(find(c==0)),x1,y1);
    im_r=im_r';
    I=find(isnan(im_r)==0);
    I=find(isnan(im_r)==1);
    %J1=max(1,I-1);J2=min(M*N,I+1);
    im_r(I)=128;
end