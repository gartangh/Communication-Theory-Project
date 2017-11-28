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
                 [0 0 0 0 0 0 1 0 1 1 0 0 0 0];]; % vul hier de generatormatrix in
            
            % Minimale Hamming afstand = 3:
            % Via H^T; minimale set van rijen die de nul-vector uitkomt (c*H^T = 0),
            % waarbij Hamming gewicht = Hamming afstand = aantal enen die set voorstelt
            
            % gegarandeerd foutcorrigerend vermogen t = 1 = floor((d_min - 1)/2)
            % gegarandeerd foutdetecterend vermogen £ = 2 = d_min - 1
            
            bitstring = bitstring(:)'; % bitstring zeker een rij vector 
            %bitstring = str2num(bitstring); % bitstring omzetten naar vector
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
            S = syndtable(H) % syndroomtabel
            
            % Manueel opstellen (p.320): H_T als syndroom, rijen aanvullen met
            % 000...0, 100...0, 0100...0, ..., 0...01, resterende rijen
            % KIEZEN door middel van exclusieve som van syndroom
           
            bitstring = bitstring(:)';
            N = length(bitstring);
            N_codewords = ceil(N/14);
            bool_error = zeros(1, N_codewords * 14);
            
            for i = 1:N_codewords
                r = bitstring(((i-1)*14+1):i*14);
                s = mod(r * H_T, 2); % syndroom
                bool_error(((i-1)*14+1):i*14) = S(bi2de(fliplr(s)) + 1, :); % geschatte foutvector
            end
            
            bitsdec = mod(bitstring + bool_error,2);
            
        end
        
        % Functie die de encoder van de inner code implementeert
        function bitenc = Encode_inner(bitstring,g_x)
           % INPUT
            % bitstring : vector met ongecodereerde bits
            % g_x : CRC-veelterm
           % OUTPUT
            % bitenc : vector met gecodeerde bits
          
        end
        
        % Functie die de decoder van de inner code implementeert
        function [bitsdec,bool_error] = Decode_inner(bitstring,g_x)
            % INPUT
             % bitstring : vector met gecodeerde bits
             % g_x : CRC-veelterm
            % OUTPUT
             % bitsdec : vector met gedecodeerde bits
             % bool_error : 1 als een fout gedetecteerd is bij zuivere foutdetectie, 0 anders
        end
        
    end
end
