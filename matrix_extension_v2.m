function M=matrix_extension_v2(N,F,epsilon,mat)

N_age=sum(N,2)';
%C1=A'*A/dot(A,N_age);
M1=mat;

F_bar=sum(N.*F',2);
N_col=reshape(N',[],1);
N_div=N./F_bar;

M_seg=zeros(length(N_col),length(N_col));
M_prop_unit=repmat(F,1,(length(F)));
M_prop_unit=M_prop_unit.*F';
M_prop=repmat(M_prop_unit,length(N_age),length(N_age));
M_prop=M_prop.*N_col';

for a=1:length(N_age)
    for b=1:length(N_age)
        age_terms=(M1(a,b)*N_age(a))/(F_bar(a)*F_bar(b)); 
        M_prop(((a-1)*length(F)+1):a*length(F),((b-1)*length(F)+1):b*length(F))=...
            age_terms*M_prop(((a-1)*length(F)+1):a*length(F),((b-1)*length(F)+1):b*length(F));
    end
end

for c=1:length(N_age)
    for d=1:length(N_age)
        if c==d
            M_seg_diag=sum(N_age(c)*M1(c,:).*(repmat(N_div(c,:),length(N_age),1)'-N_div').*(F./(2*N(c,:)')),2)+(M1(c,d)*N_age(c)/F_bar(c))*F;
        else
            M_seg_diag=M1(c,d)*N_age(c)*(N_div(c,:)+N_div(d,:))'.*(F./(2*N(c,:)'));
        end
        M_seg(length(F)*(c-1)+1:length(F)*c,length(F)*(d-1)+1:length(F)*d)=diag(M_seg_diag);
    end
end

M=(1-epsilon)*M_prop+epsilon*M_seg;
M_total_contact=M.*N_col;

M_eth=zeros(length(F),length(F));

for f=1:length(F)
    for g=1:length(F)
        total=0;
        for h=1:length(N_age)
            for k=1:length(N_age)
                total=total+M_total_contact(f+(h-1)*length(F),g+(k-1)*length(F));
            end
        end
        M_eth(f,g)=total/sum(N(:,f),1);
    end
end
M=M./(N_col');
end
