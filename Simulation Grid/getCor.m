function [A, Aindex] = getCor(n,L)
    A=zeros(n*n,2);
    Aindex = zeros(n,n);
    delta=1/n;
    for i=0:n-1
        for j=0:n-1
            A(i*n+j+1,2)=(n-1-i)*delta+delta/2;
            A(i*n+j+1,1)=j*delta+delta/2;
            Aindex(i+1,j+1)= n*i + j+1;
        end
    end
    A   =   A*L - L/2;
end