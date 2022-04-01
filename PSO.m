function [Best,fo,Gbest,Swarm]=PSO(SS,m,Phi1,Phi2,CI,dmin,dmax) %w)

% PSO estandar con coeficiente de constriccion 
% SS dimension del ejambre  
% m  dimension del espacio de busqueda 
% Phi1, Phi2 factores de aceleracion
% CI cantidad de iteraciones 
% rmin y rmax cotas del espacio de busqueda


fo=zeros(1,CI);

Gbest=fo;

[Swarm]=SwarmGen(SS,m,dmin,dmax);

 X=2/((Phi1+Phi2) - 2 + sqrt((Phi1+Phi2)^2 - 4*(Phi1+Phi2))); % factor de constriccion
 
for i=1:CI % para cada una de las iteraciones 
 
 for j=1:SS % para cada una de las particulas
 
          R1=Phi1*rand;  %  (coefiiente de aceleracion)
    
          R2=Phi2*rand; % (coefiiente de aceleracion)
                              
        if j==Swarm.g && i==1 % para el caso de la mejor particula de enjambre en la primera iteracion
    
          Swarm.v(j,:) =rand(1,m);
                            
          Swarm.r(j,:) = Swarm.r(j,:) + Swarm.v(j,:);     
                               
        else   
                                    
           g=Swarm.g;
         
          % velocidad
           
          Swarm.v(j,:) = X*(Swarm.v(j,:) + R1*(Swarm.rp(j,:) - Swarm.r(j,:)) + R2*(Swarm.rg - Swarm.r(j,:)));
             
          Swarm.r(j,:) = Swarm.r(j,:) + Swarm.v(j,:); % desplazamiento 
          
         % evaluación de las partículas 
          
         % [logIC,logR]=IC(Swarm.r(j,:));        %

         %[coef,gof]=fit(logR,logIC,'poly1');   %
  
         %Swarm.Fo(j)=ackley(Swarm.r(j));  %w(1)*gof.SSE + w(2)*Dist; % falta definir Dist 
          
         Swarm.Fo(j)=alpine(Swarm.r(j,:)); % esta es una función alpine de la base de datos 
         
         % actulización del enjambre 
         
         if Swarm.Fo(j) < Swarm.Fop(j)
            
             Swarm.Fop(j)=Swarm.Fo(j);
             
             Swarm.rp(j,:)=Swarm.r(j,:); 
         end
             
         if Swarm.Fo(j) < Swarm.Fog
                 
                 Swarm.rg=Swarm.r(j,:);
                 
                 Swarm.Fog=Swarm.Fo(j);
                 
                 Swarm.g=j;
             
         end
     

        end                        
 
 end % for 1 para cada particulas

fo(i)=Swarm.Fog;

Gbest(i)=Swarm.g;

end  % 

Best=Swarm.Fog;

end

function [Swarm]=SwarmGen(SS,m,rmin,rmax)

% Crea el enjambre de dimension [SS,m]

% w pesos de la función objetivo

% dimension del enjambre (cantidad de partículas)

% m - cantidad de puntos a ajustar 

Swarm=struct('r' ,zeros(SS,m),'rp' ,zeros(SS,m),'rg',zeros(1,m),'v',randn(SS,m),'Fo',zeros(1,SS),'Fop',zeros(1,SS),'g',0,'Fog',0);

% Swarm: .r posición de la particula , .rp mejor posicion personal, .rg posición mejor del enjambre .v velocidad,
%  .Fo funcion objetivo , .Fop funcion objetivo mejor pocision personal, .Fog mejor partícula del enjambre   

Swarm.r= rmin + (rmax-rmin)*rand(SS,m); % generar las posiciones iniciales r de cada particula

Swarm.rp=Swarm.r; % la mejor posicion  personal(de particula i) rp se iniciliza la con la posicion inicial r

for i=1:SS
    
%[logIC,logR]=IC(Swarm.r(i,:));        %

%[coef,gof]=fit(logR,logIC,'poly1');   %

%Swarm.Fo(i)=ackley(Swarm.r(i)); funcion de prueba  Ackley %w(1)*gof.SSE + w(2)*Dist; % falta definir Dist (distancia entre los puntos)

Swarm.Fo(i)=alpine(Swarm.r(i,:)); % funcion de prueba Alpine

Swarm.Fop(i)=Swarm.Fo(i);       

end

[Swarm.Fog,Swarm.g]=min(Swarm.Fo(:)); 

Swarm.rg=Swarm.r(Swarm.g,:);%

end


function [out]=ackley(x) % test function  1 

% dimension is # of columns of input, x1, x2, ..., xn
 n=length(x(1,:));
 e=exp(1);

 out = (20 + e ...
       -20*exp(-0.2*sqrt((1/n).*sum(x.^2,2))) ...
       -exp((1/n).*sum(cos(2*pi*x),2)));
end

function [out]=alpine(in) % test function 2

 out = sum(abs(in.*sin(in) + 0.1.*in),2);

end
