
t = 4; r = 4;
H=[[-0.200737+0.964821i -0.896303-0.212404i -1.48884-0.652983i -0.636742-0.160972i]
 [0.64064+0.883473i -0.759858+1.67214i -1.33804-0.469741i -1.29916+0.0361926i]
 [1.10155+0.399494i 1.67236+0.868358i -0.614789-0.0356921i 1.01917-1.04693i]
 [-0.421538-0.846632i 0.165559-0.517613i 0.565461+0.72725i 0.760191-0.280465i]];

#H = randn(r,t) + i*randn(r,t);
[U,S,V] = svd(H,'econ');
function [params1,params2,params] = given(V)
  disp("V: "),disp(V);
  D = diag(exp(-1*i*angle(V(1,:))));
  #V = V*D ; #make first row real;
  params1 = [];
  params2 = [];
  params = [];
  #disp("V: "),disp(V);
#  assert(isreal(V(1,:)), "First row must be real");
  [t,n] = size(V);
  for j=1:n
  phi = angle(V(j:t,j));
  params=[params;phi];
  params1=[params1;phi]; 
  #disp("phi"),disp(phi);
  phi=[zeros(j-1,1);phi];
  D_ = diag(exp(-1*i*phi));
  V = D_*V;
  #disp("V: "),disp(j),disp(V);
    for k = 1:t-j
    [c,s] = givens(V(t-k,j),V(t-k+1,j))
    theta = atan2(real(s),real(c))
    params = [params;theta];
    params2=[params2;theta];
    G = eye(t);
    G(t-k,t-k) = c;
    G(t-k,t-k+1) = s;
    G(t-k+1,t-k) = -s;
    G(t-k+1,t-k+1) = c;
    disp("G"),disp(G);
    V =G*V;
    #disp("V: "),disp(t-k+1),disp(V);
    
    end
  end
end


[p1,p2,p] = given(V);

