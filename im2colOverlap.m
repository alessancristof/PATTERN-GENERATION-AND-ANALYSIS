function Ic=im2colOverlap(I,m)

siz=size(I);
ar=floor(siz(1:2)/m)*m;
I1=I(1:ar(1),1:ar(2),:);

I=I(1+floor(m/2):end,1+floor(m/2):end,:);
siz=size(I);
ar=floor(siz(1:2)/m)*m;
I2=I(1:ar(1),1:ar(2),:);


siz=size(I1);
N=(siz(1)*siz(2))/(m^2);
Ic1=zeros([m^2 N size(I1,3)],'like',I);

siz=size(I2);
N=(siz(1)*siz(2))/(m^2);
Ic2=zeros([m^2 N size(I2,3)],'like',I);

for i=1:size(I1,3)
    Ic1(:,:,i) = im2col(I1(:,:,i),[m m],'distinct');
    Ic2(:,:,i) = im2col(I2(:,:,i),[m m],'distinct');
end

Ic=cat(2,Ic1,Ic2);
