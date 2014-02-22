function c=precondvec(ac,ar,type)
n=length(ac);
c=zeros(n,1);
switch type
    case 's'% strang precondition matrix.
        
    case 'c'% T.Chan precondition matrix.
        c(1)=ac(1);
        for i=2:n
            c(i)=((n-i+1)*ac(i)+(i-1)*ar(n-i+2))/n;
        end
        
    case 't'% superoptimal precondition matrix.
end