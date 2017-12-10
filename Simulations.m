classdef Simulations
    properties
        codewords = Calc_codewords();
    end
    methods (Static = true)
        function [codewords] = Calc_codewords
            persistent codewords_temp;
            codewords_temp = zeros(0,14);
            for i=0:1023
                % informatiewoord b (1x10) opstellen
                b = zeros(1,10);
                teller = dec2bin(i);
                for j=(10-length(teller)+1):10
                    b(1, j) = str2num(teller(j-10+length(teller)));
                end
                % codewoord c (1x14) berekenen
                c = Channel_Coding.Encode_outer(b);

                % checken of codewoord al in tabel zit
                [check, ~] = ismember(codewords_temp, c, 'rows');
                if(sum(check) == 0)
                    codewords_temp = [codewords_temp;c];
                end
            end
            codewords = codewords_temp;
        end
        
        % berekeningen voor Channel_Coding
        function [p_c, p_e, p_f, p_m, p_c_aprx, p_e_aprx, p_f_aprx, p_m_aprx] = Calc_prob
            % constanten
            persistent S;
            persistent p;
            S = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0];
                 [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0];
                 [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                 [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];
                 [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0];
                 [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0];];
            p = 0.05;
            
            % alle codewoorden genereren
            persistent codewords;
            
            [codewords, ~] = Simulations.Calc_codewords;
            
            % p_c (kans op correcte decodering)
            p_c = 0;
            for i = 1:size(S,1)
                w = sum(S(i,:));
                p_c = p_c + p^w*(1-p)^(14-w);
            end
            
            % benadering p_c (kans op correcte decodering)
            p_c_aprx = 1 - 90*p^2;
            
            % p_e (kans op fout)
            p_e = 0;
            for j=2:size(codewords,1)
                for k=1:16
                    e = S(k,:);
                    c_plus_e = mod(codewords(j, :) + e,2);
                    w = sum(c_plus_e);
                    p_e = p_e + p^(w)*(1-p)^(14-w);
                end
            end
            
            % benaderde p_e (kans op fout)
            p_e_aprx = 1 - p_c_aprx;
            
            % p_m (kans op niet-detecteerbare fout)
            p_m = 0;
            for j=2:size(codewords,1)
                    w = sum(codewords(j, :));
                    p_m = p_m + p^(w)*(1-p)^(14-w);
            end
            
            % benaderde p_m (kans op niet-detecteerbare fout)
            p_m_aprx = 28*p^3;
            % p_m_aprx = 0;
            % for i=3:12
            %    w = i;
            %    for j = 0:(14-w)
            %        p_m_aprx = p_m_aprx + codewords_w(i)*p^i*(factorial(14-w)/factorial(14-w-j)/factorial(j)*(-p)^j);
            %    end
            % end
                     
            % p_f (kans op decodeerfalen)
            p_f = 0;
            for j=1:size(codewords,1)
                c_plus_e = mod(codewords(j, :) + [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],2);
                w = sum(c_plus_e);
                p_f = p_f + p^(w)*(1-p)^(14-w);
            end
            
            % benaderde p_f (kans op decodeerfalen)
            p_f_aprx = 1 - p_c_aprx - p_m_aprx;
            % p_f_aprx = 0;
            % for i=1:size(codewords,1)
            %     c_plus_e = mod(codewords(1, :) + [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],2);
            %     w = sum(c_plus_e);
            %     for j=0:(14-w)
            %         p_f_aprx = p_f_aprx + p^(w)*(factorial(14-w)/factorial(14-w-j)/factorial(j)*(-p)^j);
            %     end
            % end
            
            % ontwerpcriterium
            f = @Calc_crit;
            x_0 = 0.0034;
            p_min_criterium = vpa(fzero(f,x_0));
            
            % resultaten
            probabiliteiten = [["p_c"; "p_c ben."; "p_e"; "p_e ben."; "p_f"; "p_f ben."; "p_m"; "p_m ben.";], [p_c; p_c_aprx; p_e; p_e_aprx; p_f; p_f_aprx; p_m; p_m_aprx]]        
            criterium_waarde = p_min_criterium
        end
        
        % exacte berekening voor lineaire blokcode met checkmatri
        function p_e = Calc_outer(p)
            persistent S;
            S = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0];
                 [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0];
                 [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                 [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];
                 [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0];
                 [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
                 [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0];];
            
            persistent codewords;
            codewords = Simulations.Calc_codewords;
            
            % kans op een fout
            p_c = 0;
            for i = 1:size(S,1)
                w = sum(S(i,:));
                p_c = p_c + p^w*(1-p)^(14-w);
            end
            p_e = 1-p_c;
        end
        
        % simulatie voor lineaire blokcode met checkmatrix
        function p_e = Simulate_outer(p)
            % INPUT
             % p : kans op bitovergang
            % OUTPUT
             %

            N = 1000; % testgrootte
            
            % random input string
            b = zeros(1,N*10);
            for i=1:N
                b_rand = dec2bin((floor(rand(1).*1023)));
                for j=(10-length(b_rand)+1):10
                    b(1, (i-1)*10+j) = str2num(b_rand(j-10+length(b_rand)));
                end
            end
            
            % encoderen van woord
            c = Channel_Coding.Encode_outer(b);
           
            r = c;
            for j=1:(N*14)
                p_rand = rand(1);
                if(p_rand < p)
                    if(r(1,j) == 0)
                        r(1,j) = 1;
                    else
                        r(1,j) = 0;
                    end
                end
            end
            
            [b_rec, ~] = Channel_Coding.Decode_outer(r);
            
            controle = zeros(1,N);
            
            for k=1:N
                if (isequal(b(10*(k-1)+1:10*(k)),b_rec(10*(k-1)+1:10*(k))))
                    controle(k) = 0;
                else
                    controle(k) = 1;
                end
            end
            
            p_e = sum(controle(1,:))/N;
        end
        
        function Error_outer
            N = 0.5/0.001;
            persistent p_e_ana;
            persistent p_e_sim;
            p_e_ana = zeros(N,1);
            p_e_sim = zeros(N,1);
            p = [0.001:0.001:0.5];
            for i = 1:N
                p_e_ana(i) = vpa(Simulations.Calc_outer(p(i)));
                p_e_sim(i) = vpa(Simulations.Simulate_outer(p(i)));
                %[i, p(i), p_e_ana(i), p_e_sim(i)]
            end
            %[p,p_e_ana,p_e_sim]
            
            hold on
            p_sim = loglog(p,p_e_sim, 'LineWidth', 1, 'Color', 'Red');
            p_ana = loglog(p,p_e_ana, 'LineWidth', 2, 'Color', 'Blue');
            
            
            l = legend([p_ana, p_sim], '$p_{e}$ analytisch', '$p_{e}$ experimenteel', 'Location', 'northeast');
            set(l,'Interpreter','latex')
            
            t = title('Kans op een decodeerfout, $p_{e}$');
            set(t,'Interpreter','latex')
            
            xlabel('p');
            y = ylabel('$p_{e}$');
            set(y,'Interpreter','latex')
            
            axis([0.001 0.5 0 1]);
            hold off;
        end
        
        function Simulate_inner_5()
            
        end
        
        function Simulate_inner_8()
            
        end
    end
end

