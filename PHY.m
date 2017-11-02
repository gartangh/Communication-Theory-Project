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
                    a = (bitstring-1/2)*2;
                               
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
             % T : symboolperiode in secondend
             % Ns : samples per symbool. 
             % frequency : carrier frequentie in Hz.
             % alpha : roll-off factor.
             % Lf : pulse duur (in aantal symboolperioden).
            % OUTPUT
             % s : vector met gemoduleerde samples.
                        
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
             % T : tijdsinterval van 1 symbool in seconds. 
             % alpha : rolloff factor.
            % OUTPUT
             % y : samples van de pulse
            
            % vb van gebruik:
            % alpha = 0.5;
            % t = [-5:0.1:5];
            % s = PHY.pulse(t, 1, alpha);
            % plot(t, s)
            % xlabel('tijd t/T');
            % ylabel('y(t)');
            % title(['Square root raised cosine pulse met rollofffactor ' num2str(alpha)]);
                        
            een=(1-alpha)*sinc(t*(1-alpha)/T);
            twee=(alpha)*cos(pi*(t/T-0.25)).*sinc(alpha*t/T-0.25);
            drie=(alpha)*cos(pi*(t/T+0.25)).*sinc(alpha*t/T+0.25);
            y=1/(sqrt(T))*(een+twee+drie);

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