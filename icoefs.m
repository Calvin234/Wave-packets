function x = icoefs(a0,an,bn)
    n=length(an);
    fspace=zeros(2*n,1);
    fspace(1)=a0;
    fspace(2:n+1)=an-1i*bn;
    for i = 1:n-1
        fspace(end+1-i)=an(i)+1i*bn(i);
    end
    nx = length(fspace);
    xs = real(ifft((nx/2)*fspace));
    x = ifftshift(xs);
end
  
