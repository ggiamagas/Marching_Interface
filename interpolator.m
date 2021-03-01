function [P,S,IN,U,V,W]=interpolator(Edge,isolevel,Phi,N,Psi,u,v,w,X,Y,Z)

if Edge==0
      if (abs(isolevel-Phi(1)) < 1e-5)
      P=[X(1),Y(1),Z(1)];
      elseif (abs(isolevel-Phi(2)) < 1e-5)
      P=[X(2),Y(2),Z(2)];
      elseif (abs(Phi(1)-Phi(2)) < 1e-5)
      P=[X(1),Y(1),Z(1)];
      else
      slope=(Phi(2)-Phi(1))/(X(2)-X(1));
      P(1)= X(1) - Phi(1)/slope;
  	  P(2)= Y(1); 
      P(3)= Z(1);
      end
      
      var=(u(2)-u(1))/(X(2)-X(1));
      U=u(1)+var*(P(1)-X(1));
      var=(v(2)-v(1))/(X(2)-X(1));
      V=v(1)+var*(P(1)-X(1));
      var=(w(2)-w(1))/(X(2)-X(1));
      W=w(1)+var*(P(1)-X(1));
      
      grad=(N(:,2)-N(:,1))/(X(2)-X(1));
      IN=N(:,1)+grad*(P(1)-X(1));
      
      ratio=(Psi(2)-Psi(1))/(X(2)-X(1));
      S=Psi(1)+ratio*(P(1)-X(1));
            
elseif Edge ==1
      if (abs(isolevel-Phi(2)) < 1e-5)
      P=[X(2),Y(2),Z(2)];
      elseif (abs(isolevel-Phi(3)) < 1e-5)
      P=[X(3),Y(3),Z(3)];
      elseif (abs(Phi(2)-Phi(3)) < 1e-5)
      P=[X(2),Y(2),Z(2)];
      else
      slope=(Phi(3)-Phi(2))/(Y(3)-Y(2));
      P(1)= X(2);
  	  P(2)= Y(2)-Phi(2)/slope; 
      P(3)= Z(2);
      end
      
      var=(u(3)-u(2))/(Y(3)-Y(2));
      U=u(3)+var*(P(2)-Y(2));
      var=(v(3)-v(2))/(Y(3)-Y(2));
      V=v(3)+var*(P(2)-Y(2));
      var=(w(3)-w(2))/(Y(3)-Y(2));
      W=w(3)+var*(P(2)-Y(2));
      
      grad=(N(:,3)-N(:,2))/(Y(3)-Y(2));
      IN=N(:,2)+grad*(P(2)-Y(2));
      
      ratio=(Psi(3)-Psi(2))/(Y(3)-Y(2));
      S=Psi(2)+ratio*(P(2)-Y(2));
      
elseif Edge ==2
      if (abs(isolevel-Phi(3)) < 1e-5)
      P=[X(3),Y(3),Z(3)];
      elseif (abs(isolevel-Phi(4)) < 1e-5)
      P=[X(4),Y(4),Z(4)];
      elseif (abs(Phi(3)-Phi(4)) < 1e-5)
      P=[X(3),Y(3),Z(3)];
      else
      slope=(Phi(4)-Phi(3))/(X(4)-X(3));
      P(1)= X(3)- Phi(3)/slope;
  	  P(2)= Y(3); 
      P(3)= Z(3);
      end
      
      var=(u(4)-u(3))/(X(4)-X(3));
      U=u(3)+var*(P(1)-X(3));
      var=(v(4)-v(3))/(X(4)-X(3));
      V=v(3)+var*(P(1)-X(3));
      var=(w(4)-w(3))/(X(4)-X(3));
      W=w(3)+var*(P(1)-X(3));
      
      grad=(N(:,4)-N(:,3))/(X(4)-X(3));
      IN=N(:,3)+grad*(P(1)-X(3));
      
      ratio=(Psi(4)-Psi(3))/(X(4)-X(3));
      S=Psi(3)+ratio*(P(1)-X(3));
      
elseif Edge ==3
      if (abs(isolevel-Phi(4)) < 1e-5)
      P=[X(4),Y(4),Z(4)];
      elseif (abs(isolevel-Phi(1)) < 1e-5)
      P=[X(1),Y(1),Z(1)];
      elseif (abs(Phi(4)-Phi(1)) < 1e-5)
      P=[X(4),Y(4),Z(4)];
      else
      slope=(Phi(1)-Phi(4))/(Y(1)-Y(4));
      P(1)= X(4);
  	  P(2)= Y(4)- Phi(4)/slope; 
      P(3)= Z(4);
      end  
      
      var=(u(1)-u(4))/(Y(1)-Y(4));
      U=u(4)+var*(P(2)-Y(4));
      var=(v(1)-v(4))/(Y(1)-Y(4));
      V=v(4)+var*(P(2)-Y(4));
      var=(w(1)-w(4))/(Y(1)-Y(4));
      W=w(4)+var*(P(2)-Y(4));
      
      grad=(N(:,1)-N(:,4))/(Y(1)-Y(4));
      IN=N(:,4)+grad*(P(2)-Y(4));
      
      ratio=(Psi(1)-Psi(4))/(Y(1)-Y(4));
      S=Psi(4)+ratio*(P(2)-Y(4));
      
elseif Edge ==4
      if (abs(isolevel-Phi(5)) < 1e-5)
      P=[X(5),Y(5),Z(5)];
      elseif (abs(isolevel-Phi(6)) < 1e-5)
      P=[X(6),Y(6),Z(6)];
      elseif (abs(Phi(5)-Phi(6)) < 1e-5)
      P=[X(5),Y(5),Z(5)];
      else
      slope=(Phi(6)-Phi(5))/(X(6)-X(5));
      P(1)= X(5)- Phi(5)/slope;
  	  P(2)= Y(5); 
      P(3)= Z(5);
      end  
      
      var=(u(6)-u(5))/(X(6)-X(5));
      U=u(5)+var*(P(1)-X(5));
      var=(v(6)-v(5))/(X(6)-X(5));
      V=v(5)+var*(P(1)-X(5));
      var=(w(6)-w(5))/(X(6)-X(5));
      W=w(5)+var*(P(1)-X(5));
      
      grad=(N(:,6)-N(:,5))/(X(6)-X(5));
      IN=N(:,5)+grad*(P(1)-X(5));
      
      ratio=(Psi(6)-Psi(5))/(X(6)-X(5));
      S=Psi(5)+ratio*(P(1)-X(5));
      
elseif Edge ==5 
      if (abs(isolevel-Phi(6)) < 1e-5)
      P=[X(6),Y(6),Z(6)];
      elseif (abs(isolevel-Phi(7)) < 1e-5)
      P=[X(7),Y(7),Z(7)];
      elseif (abs(Phi(6)-Phi(7)) < 1e-5)
      P=[X(6),Y(6),Z(6)];
      else
      slope=(Phi(7)-Phi(6))/(Y(7)-Y(6));
      P(1)= X(6);
  	  P(2)= Y(6)- Phi(6)/slope; 
      P(3)= Z(6);
      end   
      
      var=(u(7)-u(6))/(Y(7)-Y(6));
      U=u(6)+var*(P(2)-Y(6));
      var=(v(7)-v(6))/(Y(7)-Y(6));
      V=v(6)+var*(P(2)-Y(6));
      var=(w(7)-w(6))/(Y(7)-Y(6));
      W=w(6)+var*(P(2)-Y(6));
      
      grad=(N(:,7)-N(:,6))/(Y(7)-Y(6));
      IN=N(:,6)+grad*(P(2)-Y(6));
      
      ratio=(Psi(7)-Psi(6))/(Y(7)-Y(6));
      S=Psi(6)+ratio*(P(2)-Y(6));
      
elseif Edge ==6  
      if (abs(isolevel-Phi(7)) < 1e-5)
      P=[X(7),Y(7),Z(7)];
      elseif (abs(isolevel-Phi(8)) < 1e-5)
      P=[X(8),Y(8),Z(8)];
      elseif (abs(Phi(7)-Phi(8)) < 1e-5)
      P=[X(7),Y(7),Z(7)];
      else
      slope=(Phi(8)-Phi(7))/(X(8)-X(7));
      P(1)= X(7)- Phi(7)/slope;
  	  P(2)= Y(7); 
      P(3)= Z(7);
      end   
      
      var=(u(8)-u(7))/(X(8)-X(7));
      U=u(7)+var*(P(1)-X(7));
      var=(v(8)-v(7))/(X(8)-X(7));
      V=v(7)+var*(P(1)-X(7));
      var=(w(8)-w(7))/(X(8)-X(7));
      W=w(7)+var*(P(1)-X(7));
      
      grad=(N(:,8)-N(:,7))/(X(8)-X(7));
      IN=N(:,7)+grad*(P(1)-X(7));
      
      ratio=(Psi(8)-Psi(7))/(X(8)-X(7));
      S=Psi(7)+ratio*(P(1)-X(7));
      
elseif Edge ==7 
      if (abs(isolevel-Phi(8)) < 1e-5)
      P=[X(8),Y(8),Z(8)];
      elseif (abs(isolevel-Phi(5)) < 1e-5)
      P=[X(5),Y(5),Z(5)];
      elseif (abs(Phi(8)-Phi(5)) < 1e-5)
      P=[X(8),Y(8),Z(8)];
      else
      slope=(Phi(5)-Phi(8))/(Y(5)-Y(8));
      P(1)= X(8);
  	  P(2)= Y(8)- Phi(8)/slope; 
      P(3)= Z(8);
      end 
      
      var=(u(5)-u(8))/(Y(5)-Y(8));
      U=u(8)+var*(P(2)-Y(8));
      var=(v(5)-v(8))/(Y(5)-Y(8));
      V=v(8)+var*(P(2)-Y(8));
      var=(w(5)-w(8))/(Y(5)-Y(8));
      W=w(8)+var*(P(2)-Y(8));
      
      grad=(N(:,5)-N(:,8))/(Y(5)-Y(8));
      IN=N(:,8)+grad*(P(2)-Y(8));
      
      ratio=(Psi(5)-Psi(8))/(Y(5)-Y(8));
      S=Psi(8)+ratio*(P(2)-Y(8));
      
elseif Edge ==8
      if (abs(isolevel-Phi(1)) < 1e-5)
      P=[X(1),Y(1),Z(1)];
      elseif (abs(isolevel-Phi(5)) < 1e-5)
      P=[X(5),Y(5),Z(5)];
      elseif (abs(Phi(1)-Phi(5)) < 1e-5)
      P=[X(1),Y(1),Z(1)];
      else
      slope=(Phi(5)-Phi(1))/(Z(5)-Z(1));
      P(1)= X(1);
  	  P(2)= Y(1); 
      P(3)= Z(1)- Phi(1)/slope;
      end 
      
      var=(u(5)-u(1))/(Z(5)-Z(1));
      U=u(1)+var*(P(3)-Z(1));
      var=(v(5)-v(1))/(Z(5)-Z(1));
      V=v(1)+var*(P(3)-Z(1));
      var=(w(5)-w(1))/(Z(5)-Z(1));
      W=w(1)+var*(P(3)-Z(1));
      
      grad=(N(:,5)-N(:,1))/(Z(5)-Z(1));
      IN=N(:,1)+grad*(P(3)-Z(1));
      
      ratio=(Psi(5)-Psi(1))/(Z(5)-Z(1));
      S=Psi(1)+ratio*(P(3)-Z(1));
      
elseif Edge ==9
      if (abs(isolevel-Phi(2)) < 1e-5)
      P=[X(2),Y(2),Z(2)];
      elseif (abs(isolevel-Phi(6)) < 1e-5)
      P=[X(6),Y(6),Z(6)];
      elseif (abs(Phi(2)-Phi(6)) < 1e-5)
      P=[X(2),Y(2),Z(2)];
      else
      slope=(Phi(6)-Phi(2))/(Z(6)-Z(2));
      P(1)= X(2);
  	  P(2)= Y(2); 
      P(3)= Z(2)- Phi(2)/slope;
      end 
      
      var=(u(6)-u(2))/(Z(6)-Z(2));
      U=u(2)+var*(P(3)-Z(2));
      var=(v(6)-v(2))/(Z(6)-Z(2));
      V=v(2)+var*(P(3)-Z(2));
      var=(w(6)-w(2))/(Z(6)-Z(2));
      W=w(2)+var*(P(3)-Z(2));
      
      grad=(N(:,6)-N(:,2))/(Z(6)-Z(2));
      IN=N(:,2)+grad*(P(3)-Z(2));
      
      ratio=(Psi(6)-Psi(2))/(Z(6)-Z(2));
      S=Psi(2)+ratio*(P(3)-Z(2));
      
elseif Edge ==10
      if (abs(isolevel-Phi(3)) < 1e-5)
      P=[X(3),Y(3),Z(3)];
      elseif (abs(isolevel-Phi(7)) < 1e-5)
      P=[X(7),Y(7),Z(7)];
      elseif (abs(Phi(3)-Phi(7)) < 1e-5)
      P=[X(3),Y(3),Z(3)];
      else
      slope=(Phi(7)-Phi(3))/(Z(7)-Z(3));
      P(1)= X(3);
  	  P(2)= Y(3); 
      P(3)= Z(3)- Phi(3)/slope;
      end
      
      var=(u(7)-u(3))/(Z(7)-Z(3));
      U=u(3)+var*(P(3)-Z(3));
      var=(v(7)-v(3))/(Z(7)-Z(3));
      V=v(3)+var*(P(3)-Z(3));
      var=(w(7)-w(3))/(Z(7)-Z(3));
      W=w(3)+var*(P(3)-Z(3));
      
      grad=(N(:,7)-N(:,3))/(Z(7)-Z(3));
      IN=N(:,3)+grad*(P(3)-Z(3));
      
      ratio=(Psi(7)-Psi(3))/(Z(7)-Z(3));
      S=Psi(3)+ratio*(P(3)-Z(3));
      
elseif Edge ==11
      if (abs(isolevel-Phi(4)) < 1e-5)
      P=[X(4),Y(4),Z(4)];
      elseif (abs(isolevel-Phi(8)) < 1e-5)
      P=[X(8),Y(8),Z(8)];
      elseif (abs(Phi(4)-Phi(8)) < 1e-5)
      P=[X(4),Y(4),Z(4)];
      else
      slope=(Phi(8)-Phi(4))/(Z(8)-Z(4));
      P(1)= X(4);
  	  P(2)= Y(4); 
      P(3)= Z(4)- Phi(4)/slope;
      end
      
      var=(u(8)-u(4))/(Z(8)-Z(4));
      U=u(4)+var*(P(3)-Z(4));
      var=(v(8)-v(4))/(Z(8)-Z(4));
      V=v(4)+var*(P(3)-Z(4));
      var=(w(8)-w(4))/(Z(8)-Z(4));
      W=w(4)+var*(P(3)-Z(4));
      
      grad=(N(:,8)-N(:,4))/(Z(8)-Z(4));
      IN=N(:,4)+grad*(P(3)-Z(4));
      
      ratio=(Psi(8)-Psi(4))/(Z(8)-Z(4));
      S=Psi(4)+ratio*(P(3)-Z(4));
      
else 
 disp 'Edge is not valid cannot interpolate'
end   

end