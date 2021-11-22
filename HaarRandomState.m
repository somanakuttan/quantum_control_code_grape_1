function  v= HaarRandomState (d)
%Returns  a  random  pure  s t a t e  o f  dimension  d
 %pi cked from the Haar measure
 v = randn(d,1)+1i*randn(d,1) ;
 v = v/norm(v) ;

end
