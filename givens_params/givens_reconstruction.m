

function M = givens_reconstruct(phi,theta,t,n)
  q=1;
  w=1;
  M=eye(t);
  for j=1:n
     p=phi(q:q+t-j);
     p=[zeros(j-1,1);p];
     M=M*diag(exp(i*p))
     q=q+t-j+1;
    
    T = theta(w:w-1+ t-j);
    
    L = size(T);
    w = w+t-j;
    for l=1:L
      G = given_matrix(T(l),t-l,t-l+1,t)
      disp("M"),disp(M);
      M=M*G
      
    end
  end
  I_ = eye(n);
  I_ = [I_,zeros(n,t-n)];
  M=M*(I_)';
  
end
function G =  given_matrix(theta,a,b,s)
  G=eye(s);
  G(a,a) = cos(theta);
  G(a,b) = -sin(theta);
  G(b,a) = sin(theta);
  G(b,b) = cos(theta);
end