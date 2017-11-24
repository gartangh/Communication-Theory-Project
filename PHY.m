classdef PHY
          
   methods(Static=true)
        
             
        % Functie die een bitstring omzet naar complexe symbolen.
        function a = mapper(bitstring, constellation)
            % INPUT
             % bitstring: vector van ongecodeerde bits met lengte 1xN.
             % constellation : ofwel 'BPSK', '4QAM', '4PSK', '4PAM', etc.
            % OUTPUT
             % a : vector met (complexe) symbolen. 
            
                       
            switch(constellation)
                case 'BPSK'                    
                    % BPSK mapping here.
                    a = (bitstring-1/2)*2; % 0 -> -1 & 1 -> 1.
                    
                case '4QAM'
                    % 4QAM mapping here.
                    dQAM = sqrt(6/(4-1));
                    deltaQAM = dQAM/2;
                    
                    a = zeros(length(bitstring)/2);
                    
                    for j = 1:length(bitstring)/2
                        x = bitstring(2*j-1);
                        y = bitstring(2*j);
                        
                        if (x==0 && y==0)
                            a(j) = sqrt(2)*exp(1*1i*pi/4);
                        elseif (x==0 && y==1)
                            a(j) = sqrt(2)*exp(3*1i*pi/4);
                        elseif (x==1 && y==0)
                            a(j) = sqrt(2)*exp(5*1i*pi/4);
                        elseif(x==1 && y==1)
                            a(j) = sqrt(2)*exp(7*1i*pi/4);
                        else
                            % Incorrect input
                        end
                    end
                    
                case '4PAM'
                    % 4PAM mapping here.
                    dPAM = sqrt(12/(4^2-1));
                    deltaPAM = dPAM/2;
                    
                    a = zeros(length(bitstring)/2);
                    
                    for j = 1:length(bitstring)/2
                        x = bitstring(2*j-1);
                        y = bitstring(2*j);
                        
                        if (x==0 && y==0)
                            a(j) = -3*deltaPAM;
                        elseif (x==0 && y==1)
                            a(j) = -1*deltaPAM;
                        elseif (x==1 && y==0)
                            a(j) = 1*deltaPAM;
                        elseif(x==1 && y==1)
                            a(j) = 3*deltaPAM;
                        else
                            % Incorrect input
                        end
                    end
                               
                otherwise
                    error('Constellation not recognized');
            end
                        
        end
        
        % Functie die complexe symbolen omzet naar bits.
        function bitstring = demapper(a, constellation)
            % INPUT
             % a : vector met (complexe) symbolen.
             % constellation : ofwel 'BPSK', 'QAM', '4PSK', '4PAM', etc.  
            % OUTPUT
             % bitstring : vector met bits horend bij a.
            
            switch(constellation)
                case 'BPSK'                       
                    % BPSK demapping here.
                    bitstring = (a+1)/2;
                    
                case '4QAM'
                    % 4QAM mapping here.
                    dQAM = sqrt(6/(4-1));
                    deltaQAM = dQAM/2;
                    
                    bitstring = zeros(length(a)*2);
                    
                    for j = 1:length(a)
                        if (a(j) == sqrt(2)*exp(1*1i*pi/4))
                            bitstring(2*j-1)    = 0;
                            bitstring(2*j)      = 0;
                        elseif (a(j) == sqrt(2)*exp(3*1i*pi/4))
                            bitstring(2*j-1)    = 0;
                            bitstring(2*j)      = 1;
                        elseif (a(j) == sqrt(2)*exp(5*1i*pi/4))
                            bitstring(2*j-1)    = 1;
                            bitstring(2*j)      = 0;
                        elseif (a(j) == sqrt(2)*exp(7*1i*pi/4))
                            bitstring(2*j-1)    = 1;
                            bitstring(2*j)      = 1;
                        else
                            % Incorrect input
                        end
                    end
                    
                case '4PAM'
                    % 4PAM mapping here.
                    dPAM = sqrt(12/(4^2-1));
                    deltaPAM = dPAM/2;
                    
                    bitstring = zeros(length(a)*2);
                    
                    for j = 1:length(a)
                        if (a(j) == -3*deltaPAM)
                            bitstring(2*j-1)    = 0;
                            bitstring(2*j)      = 0;
                        elseif (a(j) == -1*deltaPAM)
                            bitstring(2*j-1)    = 0;
                            bitstring(2*j)      = 1;
                        elseif (a(j) == 1*deltaPAM)
                            bitstring(2*j-1)    = 1;
                            bitstring(2*j)      = 0;
                        elseif (a(j) == 3*deltaPAM)
                            bitstring(2*j-1)    = 1;
                            bitstring(2*j)      = 1;
                        else
                            % Incorrect input
                        end
                    end
                  
                otherwise
                    error('Constellation not recognized');
            end
        end
        
        % Functie die harde desicie toepast op x.
        function a_estim = hard_decisions(x, constellation)
            % INPUT
             % x : vector met ruizige (complexe) symbolen.
             % constellation: ofwel 'BPSK', '4QAM', '4PSK', '4PAM', etc.   
            % OUTPUT
             % a_estim : vector met geschatte (complexe) symbolen.           
            
            switch(constellation)
                case 'BPSK'                   
                    % BPSK here      
                    
                otherwise
                    error('Constellation not recognized');
                    
            end
                
                
        end
        
        % Functie die de decisie-variabele afleidt uit rdown.
        function u = make_decision_variable(rdown,hch_hat,theta_hat)
            % INPUT
             % rdown : vector met het gedecimeerde ontvangen signaal.
             % hch_hat : schatting van amplitde van het kanaal.
             % theta_hat : schatting van fase van het kanaal.
            % OUTPUT
             % u : vector met decisie-variabele.
            
        end
        
        % Functie die symbolen op een puls zet en dit signaal moduleert op een golf.
        function s = modulate(a,T,Ns,frequency,alpha,Lf)
            % INPUT
             % a : vector van symbolen.
             % T : symboolperiode in seconden.
             % Ns : samples per symbool. 
             % frequency : carrier frequentie in Hz.
             % alpha : roll-off factor.
             % Lf : pulse duur (in aantal symboolperioden).
            % OUTPUT
             % s : vector met gemoduleerde samples.
             
            % Bemonsteringstheorema van Shannon-Nyquist
             % De bemonsteringsfrequentie Ns/T moet minstens tweemaal zo hoog zijn
             % als de hoogste frequentie, aanwezig in het signaal,
             % om het origineel zonder fouten te kunnen reproduceren.
            tijdsspan = -Lf*T:T/Ns:Lf*T;
            x = zeros(length(a)*length(t));
            
            for i = 1:length(a)
                % optelling uit opgave (3) hier gebruiken!
                y = a(i)*pulse(tijdsspan,T,alpha);
                x((i-1)*length(y)+1:i*length(y)) = y;
            end
            
            % Basisbandkanaal op draaggolf met freqentie frequency
            s = zeros(lenght(x));
            s(i) = sqrt(2)*real(x(i)*exp(1i*2*pi*frequency*i*T/length(t)));
            
        end
        
        function rdemod = demodulate(r,T,Ns,frequency,alpha,Lf,theta)
            % Demoduleert het signaal r en voert het matched filter uit. 
            % INPUT
             % r: vector met ontvangen samples.
             % T : symboolperiode in secondend
             % Ns : samples per symbool. 
             % frequency : carrier frequentie in Hz.
             % alpha : roll-off factor.
             % Lf : pulse duur (in aantal symboolperioden).
            % OUTPUT
             % rdemod: vector met gedemoduleerde samples.
                        
        end
        
        function y = pulse(t,T,alpha)
            % Functie die de square root raised cosine pulse maakt
            % INPUT
             % t : samples van tijd.
             % T : tijdsinterval van 1 symbool in seconden.
             % alpha : rolloff factor.
            % OUTPUT
             % y : samples van de pulse
            
            % vb van gebruik:
            % alpha = 0.5;
            % T = 1;
            % t = [-5:0.1:5];
            % s = PHY.pulse(t, T, alpha);
            % plot(t, s)
            % xlabel('tijd t/T');
            % ylabel('y(t)');
            % title(['Square root raised cosine pulse met rollofffactor ' num2str(alpha)]);
                        
            een  = (1-alpha)*sinc(t*(1-alpha)/T);
            twee = (alpha)*cos(pi*(t/T-0.25)).*sinc(alpha*t/T-0.25);
            drie = (alpha)*cos(pi*(t/T+0.25)).*sinc(alpha*t/T+0.25);
            y    = 1/(sqrt(T)) * (een + twee + drie);

        end
        
        % Functie die het gegeven signaal decimeert met factor Ns en het overgansverschijnsel verwijdert
        function rdown=downsample(rdemod,Ns,Lf)
            % INPUT
             % rdemod : vector met Ns samples per symbool.
             % Ns : samples per symbool.
             % Lf : pulse duration (in number of symbol periods)
            % OUTPUT
             % rdown: vector met 1 sample per symbool
        end   
        
        % Functie die het kanaal simuleert
        function r = channel(s, hch, sigma)
            % INPUT
             % s : rijvector met samples van ingang van het kanaal.
             % hch : versterking van het kanaal. 
             % sigma : standaard deviatie van de ruis.
            % OUTPUT
             % r : rijvector met samples van de uitgang van het kanaal.
            
        end
            
    end
    
end