classdef Channel_Coding
    
    methods(Static=true)
        
        % Functie die de encoder van de outer code implementeert
        function [bitenc] = Encode_outer(bitstring)
            % INPUT
             % bitstring : vector met ongecodereerde bits
            % OUTPUT
             % bitenc : vector met gecodeerde bits
            
            G = [[1 1 0 0 0 0 0 0 1 0 0 0 0 0];
                 [0 1 0 0 0 1 0 0 1 0 1 0 0 0];
                 [0 0 1 0 0 1 1 0 0 0 0 0 0 0];
                 [0 1 0 1 0 1 1 0 0 0 0 0 0 0];
                 [0 1 0 0 1 1 0 0 0 0 0 0 0 0];
                 [0 0 0 0 0 1 1 0 1 0 0 1 0 0];
                 [0 0 0 0 0 1 0 0 1 0 0 0 1 0];
                 [0 1 0 0 0 1 1 1 1 0 0 0 0 0];
                 [0 1 0 0 0 0 1 0 0 0 0 0 0 1];
                 [0 0 0 0 0 0 1 0 1 1 0 0 0 0];]; % generatormatrix
            
            % Minimale Hamming afstand d_min = 3:
            % Via H^T; minimale set van rijen die de nul-vector uitkomt (c*H^T = 0),
            % waarbij Hamming gewicht = Hamming afstand = aantal enen die set voorstelt
            
            % gegarandeerd foutcorrigerend vermogen t = 1 = floor((d_min - 1)/2)
            % gegarandeerd foutdetecterend vermogen £ = 2 = d_min - 1
            
            bitstring = bitstring(:)'; % bitstring zeker een rij vector 
            N = length(bitstring); 
            N_codewords = ceil(N/10);
            bitstring = [bitstring, zeros(1, N_codewords*10-N)]; % vul aan met nullen als bitstring geen geheel aantal informatiewoorden is.             
            
            bitenc = zeros(1,N_codewords*14);
    
            for i = 1:N_codewords
                b = bitstring(((i-1)*10+1):i*10); 
                bitenc(1, ((i-1)*14+1):i*14) = mod(b * G,2);
            end
        end
        
        % Functie die de decoder van de outer code implementeert
        function [bitsdec,bool_error]  = Decode_outer(bitstring)
            % INPUT
             % bitstring : vector met gecodeerde bits
            % OUTPUT
             % bitsdec : vector met gedecodeerde bits bij volledige foutcorrectie
             % bool_error : 1 als een fout gedetecteerd is bij zuivere foutdetectie, 
             % 0 anders
            
            H = [[1 1 0 1 1 0 0 1 0 0 1 0 0 1];
                 [0 0 1 1 1 1 0 1 0 0 1 1 1 0];
                 [0 0 1 1 0 0 1 1 0 1 0 1 0 1];
                 [1 0 0 0 0 0 0 1 1 1 1 1 1 0]]; % checkmatrix
            H_T = H.'; % getransponeerde checkmatrix
             
            % MANUEEL OPSTELLEN (p.320):
            % H_T als syndroom, rijen aanvullen met 000...0, 100...0, 
            % 0100...0, ..., 0...01, resterende rijen KIEZEN door middel 
            % van exclusieve som van syndroom
            % AUTOMATISCH LATEN GENEREREN MET MATLAB:            
            S = syndtable(H); % syndroomtabel in volgorde van oplopend syndroom
           
            bitstring = bitstring(:)'; % bitstring zeker een rij vector 
            N = length(bitstring); 
            N_codewords = ceil(N/14);
            
            bool_error = zeros(1,N_codewords);
            % foutcorrectie ontvangen codewoord
            for i = 1:N_codewords
                r_i = bitstring(((i-1)*14+1):i*14); % i-de ontvangen codewoord
                s_i = mod(r_i * H_T, 2); % syndroom horend bij i-de woord
                error_vector = S(bi2de(fliplr(s_i)) + 1, :); % geschatte foutvector horend bij syndroom
                bitstring(((i-1)*14+1):i*14) = mod((r_i + error_vector), 2); % i-de woord in bitstring corrigeren a.d.h.v. geschatte foutvector
                
                if (sum(error_vector(:)) ~= 0)
                   bool_error(i) = 1; % als foutpatroon minstens één 1 bevat is er een fout opgetreden
                end
            end
            
            % decoderen ontvangen codewoord
            for i = 1:N_codewords
               r_corr = bitstring(((i-1)*14+1):i*14); % i-de gecorrigeerde ontvangen codewoord
               
               % zelfde permutatie uitvoeren als voor G naar G_sys om
               % systematische deel van ontvangen woord af te splitsen
               bitstring_sys = r_corr(1:10); % gedecodeerde ontvangen informatiewoord vóór bovenvermelde permutaties
               bitstring_sys(2) = r_corr(11);
               bitstring_sys(6) = r_corr(12);
               bitstring_sys(7) = r_corr(13);
               bitstring_sys(9) = r_corr(14);
               
               bitsdec(((i-1)*10+1):i*10) = bitstring_sys; % i-de informatiewoord
            end
        end
        
        % Functie die de encoder van de inner code implementeert
        function bitenc = Encode_inner(bitstring,g_x)
           % INPUT
            % bitstring : vector met ongecodereerde bits
            % g_x : CRC-veelterm
           % OUTPUT
            % bitenc : vector met gecodeerde bits
           
           % bepalen van dimensies
           k = size(bitstring,2);
           n = k + size(g_x, 2) - 1; % size(g_x) = n - k + 1
           
           % berekenen CRC-deel van codewoord
           bitstring = bitstring(:)'; % bitstring zeker een rij vector 
           bitstring = [bitstring, zeros(1,n - k)];
           [~,s] = deconv(bitstring, g_x);
           s = mod(s, 2);
           
           % informatiewoord uitbreiden met CRC-deel
           bitenc = mod(bitstring + s, 2);
        end
        
        % Functie die de decoder van de inner code implementeert
        function [bitsdec,bool_error] = Decode_inner(bitstring,g_x)
            % INPUT
             % bitstring : vector met gecodeerde bits
             % g_x : CRC-veelterm
            % OUTPUT
             % bitsdec : vector met gedecodeerde bits
             % bool_error : 1 als een fout gedetecteerd is bij zuivere foutdetectie, 0 anders
            
            % bepalen van dimensies
            n = size(bitstring, 2);
            k = n - size(g_x, 2) + 1; % size(g_x, 2) = n - k + 1 
            
            bitstring = bitstring(:)'; % bitstring zeker een rij vector 
            
            bitsdec = bitstring(1:k); % informatiewoord van volledige codewoord afsplitsen
            
            % controle door rest van gedecodeerde codewoord te vergelijken
            % met 0
            [~,s_controle] = deconv(bitstring, g_x);
            s_controle = mod(s_controle,2);
            if sum(s_controle) ~= 0
                bool_error = 1;
            else
                bool_error = 0;
            end
                
        end
    end
end
