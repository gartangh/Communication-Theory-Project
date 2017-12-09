classdef PHY
          
   methods(Static=true)
         
        % Functie die een bitstring omzet naar complexe symbolen.
        function a = mapper(bitstring, constellation)
            % INPUT
             % bitstring: vector van ongecodeerde bits met lengte 1xN.
             % constellation : ofwel 'BPSK', '4QAM', '4PSK', '4PAM', etc.
            % OUTPUT
             % a : vector met (complexe) symbolen. 
            N = numel(bitstring);
            if bitget(N,1)
                error('Bitstring moet een even lengte hebben');
            end
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
                            a(j) = sqrt(2)*deltaQAM*exp(1*1i*pi/4);
                        elseif (x==0 && y==1)
                            a(j) = sqrt(2)*deltaQAM*exp(3*1i*pi/4);
                        elseif (x==1 && y==1)
                            a(j) = sqrt(2)*deltaQAM*exp(5*1i*pi/4);
                        elseif(x==1 && y==0)
                            a(j) = sqrt(2)*deltaQAM*exp(7*1i*pi/4);
                        else
                            % Incorrect input
                        end
                    end
                    
                case '4PAM'
                    % 4PAM mapping here.
                    dPAM = sqrt(12/(4^2-1));
                    deltaPAM = dPAM/2;
                    a = zeros(1, length(bitstring)/2);
                    
                    for j = 1:length(bitstring)/2
                        % Wordt gelezen van links naar rechts
                        decimal = bi2de([bitstring(2*j) bitstring(2*j-1)]);
                        disp(decimal);
                        switch(decimal)
                            case 0
                                a(j) = -3*deltaPAM;
                            case 1
                                a(j) = -1*deltaPAM;
                            case 3
                                a(j) = 1*deltaPAM;
                            case 2
                                a(j) = 3*deltaPAM;
                            otherwise
                                error('De input is incorrect');
                        end;
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
                        if      (a(j) == sqrt(2)*deltaQAM*exp(1*1i*pi/4))
                            bitstring(2*j-1)    = 0;
                            bitstring(2*j)      = 0;
                        elseif  (a(j) == sqrt(2)*deltaQAM*exp(3*1i*pi/4))
                            bitstring(2*j-1)    = 0;
                            bitstring(2*j)      = 1;
                        elseif  (a(j) == sqrt(2)*deltaQAM*exp(5*1i*pi/4))
                            bitstring(2*j-1)    = 1;
                            bitstring(2*j)      = 1;
                        elseif  (a(j) == sqrt(2)*deltaQAM*exp(7*1i*pi/4))
                            bitstring(2*j-1)    = 1;
                            bitstring(2*j)      = 0;
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
                            bitstring(2*j)      = 1;
                        elseif (a(j) == 3*deltaPAM)
                            bitstring(2*j-1)    = 1;
                            bitstring(2*j)      = 0;
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
            
            a_estim = zeros(length(x)); 
            
            switch(constellation)
                case 'BPSK'                   
                    % BPSK here
                    for i = 1:length(x)
                        if real(x(i)) < 0
                            a_estim(i) = -1;
                        else
                            a_estim(i) = 1;
                        end
                    end
                    
                case '4QAM'
                    % 4QAM here            
                    for i = 1:length(x)
                        if (real(x(i)) >= 0 && imag(x(i)) >= 0)
                            a_estim(i) = sqrt(2)*deltaQAM*exp(1*1i*pi/4);
                        elseif (real(x(i)) < 0 && imag(x(i)) > 0)
                            a_estim(i) = sqrt(2)*deltaQAM*exp(3*1i*pi/4);
                        elseif (real(x(i)) <= 0 && imag(x(i)) <= 0)
                            a_estim(i) = sqrt(2)*deltaQAM*exp(5*1i*pi/4);
                        elseif (real(x(i)) > 0 && imag(x(i)) < 0)
                            a_estim(i) = sqrt(2)*deltaQAM*exp(7*1i*pi/4);
                        else
                            % Incorrect input
                        end
                    end
                    
                case '4PAM'
                    % 4PAM here
                    dPAM = sqrt(12/(4^2-1));
                    deltaPAM = dPAM/2;
                    
                    for i = 1:length(x)
                        a_estim(i) = round(real(x(i))/deltaPAM/2) * 2;
                        if (a_estim(i) >= 0)
                            a_estim(i) = a_estim(1) + 1;
                        else
                            a_estim(i) = a_estim(1) - 1;
                        end
                    end
                
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
            
            u = rdown/hch_hat*exp(1i*theta_hat);
            
            scatter(real(u),imag(u));   
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
            
            t = -Lf*T:T/Ns:Lf*T;
            p = pulse(t,T,alpha);
            % optelling uit opgave (2)
            x = zeros(length(t),1);
            for k = 0:length(a)-1
                x(t/T*Ns+Lf*T) = x(t/T*Ns+Lf*T) + a(k)*p(t/T*Ns+Lf*T - k*T/T*Ns+Lf*T);
            end
            
            % Basisbandkanaal op draaggolf met freqentie frequency
            s(t/T*Ns+Lf*T) = sqrt(2)*real(x(t/T*Ns+Lf*T)*exp(1i*2*pi*frequency*t*T)); 
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
            
            t = -Lf*T:T/Ns:Lf*T;
            p = pulse(t,T,alpha);
            
            r(t/T*Ns+Lf*T) = sqrt(2)*r(t/T*Ns+Lf*T)*exp(-1i*(2*pi*frequency*t*Ns/T + theta));
            
            rdemod = T/Ns * conv(r,p);
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
        
        % Functie die het gegeven signaal decimeert met factor Ns en het overgangsverschijnsel verwijdert
        function rdown=downsample(rdemod,Ns,Lf)
            % INPUT
             % rdemod : vector met Ns samples per symbool.
             % Ns : samples per symbool.
             % Lf : pulse duration (in number of symbol periods)
            % OUTPUT
             % rdown: vector met 1 sample per symbool
            
            % Decimeren met factor Ns
            rdown = decimate(rdemod, Ns);
            
            % Overgangsverschijnsel verwijderen
            rdown = rdown; % Hier nog formule op toepassen!
        end   
        
        % Functie die het kanaal simuleert
        function r = channel(s, hch, sigma)
            % INPUT
             % s : rijvector met samples van ingang van het kanaal.
             % hch : versterking van het kanaal. 
             % sigma : standaard deviatie van de ruis.
            % OUTPUT
             % r : rijvector met samples van de uitgang van het kanaal.
            
            nl = sqrt(sigma) * (1/sqrt(2) * (randn(1,lenght(s)) + 1i * randn(1,lenght(s))));
            r = hch*s + nl;
        end
   end
end