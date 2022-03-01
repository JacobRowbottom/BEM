%clear
NT=200;
MT=200;
NP=25;

% IPx=0.25;
% IPy=0.25;

IPx=linspace(0,1,NP+2); %gird for square
IPy=linspace(0,1,NP+2);

IPx = IPx(2:end-1);
IPy = IPy(2:end-1);

[IPX,IPY] =meshgrid(IPx,IPy);
 
for ipx=1:NP 
    for ipy=1:NP
        for nn=1:NT
            NN=nn-1;
            for mm=1:MT
                MM=mm-1;
                lambdaNM = (pi^2)*(MM^2 + NN^2);
                unmIP = cos(NN*pi*IPx(ipx))*cos(MM*pi*IPy(ipy));
                unmPS = cos(NN*pi*PSx)*cos(MM*pi*PSy);
                cnm=4;
                if MM==0 && NN==0
                    cnm=1;
                elseif MM==0 && NN>0
                    cnm=2;
                elseif NN==0 && MM>0
                    cnm=2;
                end
                
                GkExact(nn,mm)= (cnm*unmIP*unmPS)/(k^2 - lambdaNM);
                
            end            
        end       
        GkE(ipx,ipy) = (sum(sum(GkExact)));
    end    
end
figure
mesh(IPX,IPY,GkE);