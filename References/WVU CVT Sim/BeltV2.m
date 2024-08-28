function [Alpha,Beta] = BeltV2(n,Rp,Rs,cc)
%Calculates angle of wrap at primary and secondary through shift
%Belt Wrap
A=Rs>=Rp; %Determine when secondary has larger working radius than primary
A=A(:,any(A,1));
%Preallocation
Alpha=zeros(1,n); %Primary belt wrap through shift (Degrees)
Beta=zeros(1,n); %Secondary belt wrap through shift (Degrees)
for i=1:length(A)
    Alpha(i) = (pi-(2*asin((Rs(i)-Rp(i))/cc)));
    Beta(i)= (pi+(2*asin((Rs(i)-Rp(i))/cc)));
end
for i=(length(A)+1):n
    Alpha(i) = (pi+(2*asin((Rp(i)-Rs(i))/cc)));
    Beta(i) = (pi-(2*asin((Rp(i)-Rs(i))/cc)));
end
end