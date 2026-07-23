function indicators=vacc_indicator_function(m,cmlt,N,S,b,a)

L=(m-(cmlt./N)).*S;
V=ones(length(L),1);

f=1-min([V,exp(-L/b)],[],2);
L_plus=zeros(length(L),1);

for n=1:length(L)
    L_plus(n)=0;
    for m=1:length(L)
        if a(m) < a(n)
            L_plus(n)=L_plus(n)+L(m);
        end
    end
end

e=min([V,exp(-L_plus/b)],[],2);

indicators=[f,e];
