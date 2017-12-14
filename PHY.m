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
                    a = (bitstring-1/2)*2; % 0 -> -1 & 1 -> 1;
                    
                case '4QAM'
                    % 4QAM mapping here.
                    dQAM = sqrt(12/(4-1)); % 2
                    deltaQAM = dQAM/2; % 1
                    
                    a = zeros(1,length(bitstring)/2);
                   
                    for j = 1:length(bitstring)/2
                        % Eerst mappen van 0 -> -1 en 1 -> 1 en deze dan
                        % als coordinaten gebruiken. Door de definitie van
                        % Gray coding komt dit dan neer op telkens de
                        % dichtsbijzijnde afstand (vanaf linksbeneden, met
                        % de klok mee.)
                        a(j) = (bitstring(2*j-1) - 1/2)*2 * deltaQAM + (bitstring(2*j) - 1/2)*2 * deltaQAM * 1i;
                        a(j) = 1/sqrt(2) * a(j);
                    end     
                    
                case '4PAM'
                    % 4PAM mapping here.
                    dPAM = sqrt(12/(4^2-1)); % sqrt(12/15)
                    deltaPAM = dPAM/2; % sqrt(12/15)/2
                    
                    a = zeros(1,length(bitstring)/2);
                    
                    for j = 1:length(bitstring)/2
                        % Wordt gelezen van links naar rechts
                        decimal = bi2de([bitstring(2*j) bitstring(2*j-1)]);
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
                    bitstring = zeros(1,length(a)*2);
                    
                    for j = 1:length(a)
                        bitstring(2*j-1) = (real(a(j))+1)/2;
                        bitstring(2*j) = (imag(a(j))+1)/2;
                    end
                    
                case '4PAM'
                    % 4PAM mapping here.
                    dPAM = sqrt(12/(4^2-1));
                    deltaPAM = dPAM/2;
                    
                    bitstring = zeros(1,length(a)*2);
                    
                    for j = 1:length(a)
                        switch(a(j))
                            case -3*deltaPAM
                                bitstring(2*j-1:2*j) = [0 0];
                            case -1*deltaPAM
                                bitstring(2*j-1:2*j) = [0 1];
                            case 1*deltaPAM
                                bitstring(2*j-1:2*j) = [1 1];
                            case 3*deltaPAM
                                bitstring(2*j-1:2*j) = [1 0];
                            otherwise
                                error('Incorrecte input');
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
            
            a_estim = zeros(1,length(x));
            
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
                            a_estim(i) = 1+j;
                        elseif (real(x(i)) < 0 && imag(x(i)) > 0)
                            a_estim(i) = -1 + j;
                        elseif (real(x(i)) <= 0 && imag(x(i)) <= 0)
                            a_estim(i) = -1-j;
                        elseif (real(x(i)) > 0 && imag(x(i)) < 0)
                            a_estim(i) = 1-j;
                        else
                            error('Incorrecte input');
                        end
                        a_estim(i) = 1/sqrt(2) * a_estim(i);
                    end
                    
                case '4PAM'
                    % 4PAM here
                    dPAM = sqrt(12/(4^2-1));
                    deltaPAM = dPAM/2;
                    
                    for i = 1:length(x)
                        possibilities = [-3*deltaPAM, -deltaPAM, deltaPAM, 3*deltaPAM];
                        [~, index] = min(abs(possibilities-x(i)));
                       a_estim(i) = possibilities(index);
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
             % Ns/T >= 2(f0+B)
            
            t = -Lf*T:T/Ns:Lf*T;
            p = PHY.pulse(t,T,alpha);
            % optelling uit opgave (2)
            x = zeros(1,2*Lf*Ns + (length(a)-1)*Ns + 1); 
            for k = 0:length(a)-1
                x(k*Ns + 1 : k*Ns + 2*Lf*Ns + 1) = x(k*Ns + 1 : k*Ns + 2*Lf*Ns + 1) + a(k+1)*p;
            end
            
            % Basisbandkanaal op draaggolf met freqentie frequency
            s = zeros(1, length(x));
            for k = 1:length(x)
                s(k) = sqrt(2)*real(x(k) *exp(1j*2*pi*frequency*k*T/Ns));
            end
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
            p = PHY.pulse(t,T,alpha);
            r1 = zeros(1, length(r));
            for i = 1: length(r)
                r1(i) = sqrt(2)*r(i)*exp(-1*j*(2*pi*frequency*i*T/Ns + theta));
            end
            
            rdemod = T/Ns * conv(r1,p);
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
            rdemod_trim = rdemod(2*Lf*Ns+1:length(rdemod) - 2*Lf*Ns + 1);
            
            rdown = [];
            i = 1;
            while i < length(rdemod_trim)
                rdown = [rdown rdemod_trim(i)];
                i = i + Ns;
            end
        end   
        
        % Functie die het kanaal simuleert
        function r = channel(s, hch, sigma)
            % INPUT
             % s : rijvector met samples van ingang van het kanaal.
             % hch : versterking van het kanaal. 
             % sigma : standaard deviatie van de ruis.
            % OUTPUT
             % r : rijvector met samples van de uitgang van het kanaal.

            r = hch*s + sqrt(sigma)*randn(1, length(s));
        end
    end
end
