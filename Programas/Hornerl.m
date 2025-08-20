function Hornerl(a,n,uO,C)
 %Compute point on power basis curve.
 %Input: a,n,uO *1
 %Output: C *1
  C = a(:,n);

  for i=n:-1:1
    C = C*uO + a(:,i);
  endfor
end
