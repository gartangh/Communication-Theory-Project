classdef Channel_Coding
    
    methods(Static=true)
        
        % Functie die de encoder van de outer code implementeert
        function [bitenc] = Encode_outer(bitstring)
            % INPUT
             % bitstring : vector met ongecodereerde bits
            % OUTPUT
             % bitenc : vector met gecodeerde bits
            G = []; % vul hier de generatormatrix in

            bitstring = bitstring(:)'; % bitstring zeker een rij vector 
            N = length(bitstring);
            N_codewords = ceil(N/10);
            bitstring = [bitstring zeros(1, N_codewords*10-N)]; % vul aan met nullen als bitstring geen geheel aantal informatiewoorden is. 
        end
        
        % Functie die de decoder van de outer code implementeert
        function [bitsdec,bool_error]  = Decode_outer(bitstring)
            % INPUT
             % bitstring : vector met gecodeerde bits
            % OUTPUT
             % bitsdec : vector met gedecodeerde bits bij volledige foutcorrectie
             % bool_error : 1 als een fout gedetecteerd is bij zuivere foutdetectie, 0 anders
            H = [[1 1 0 1 1 0 0 1 0 0 1 0 0 1];[0 0 1 1 1 1 0 1 0 0 1 1 1 0 ];[0 0 1 1 0 0 1 1 0 1 0 1 0 1];[1 0 0 0 0 0 0 1 1 1 1 1 1 0]]; % checkmatrix
            
            bitstring = bitstring(:)';
            N = length(bitstring);
            N_codewords = ceil(N/14);
                        
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
