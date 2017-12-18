classdef Simulations
    properties
    end
    methods (Static = true)
        % genereren alle codewoorden
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
        
        function Calc_prob
            Simulations.Calc_prob_corr;
            Simulations.Calc_prob_det;
        end
        
        % berekeningen voor Channel_Coding
        function [p_c, p_e, p_f, p_m, p_c_aprx, p_e_aprx, p_f_aprx, p_m_aprx] = Calc_prob_corr
            % syndroomtabel
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
             
            % bitovergangswaarschijnlijkheid
            persistent p;
            p = 0.05;
            
            % array van alle codewoorden
            persistent codewords;
            codewords = Simulations.Calc_codewords;
            
            % p_c (kans op correcte decodering)
            p_c = 0;
            for i = 1:size(S,1)
                w = sum(S(i,:));
                p_c = p_c + p^w * (1-p)^(14-w);
            end
            
            % benadering p_c (kans op correcte decodering)
            p_c_aprx = 1 - 90 * p^2;
            
            % p_e (kans op fout)
            p_e = 0;
            for j=2:size(codewords,1) % alle codewoorden behalve nulwoord
                for k=1:size(S,1) % alle syndromen
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
            %{
            p_m_aprx = 0;
            for i=3:12
                w = i;
                for j = 0:(14-w)
                    p_m_aprx = p_m_aprx + codewords_w(i)*p^i*(factorial(14-w)/factorial(14-w-j)/factorial(j)*(-p)^j);
                end
            end
            %}      
            % p_f (kans op decodeerfalen)
            p_f = 0; % alle foutpatronen staan in tabel dus geen sommatie    
            
            % benaderde p_f (kans op decodeerfalen)
            p_f_aprx = 0;
            
            % ontwerpcriterium
            f = @Calc_crit;
            x_0 = 0.0034;
            p_min_criterium = vpa(fzero(f,x_0));
            
            % resultaten
            probabiliteiten = [["VOLLEDIGE DECODERING";"p_c"; "p_c ben."; "p_e"; "p_e ben."; "p_f"; "p_f ben."; "p_m"; "p_m ben.";], ["";p_c; p_c_aprx; p_e; p_e_aprx; p_f; p_f_aprx; p_m; p_m_aprx]]        
            %criterium_waarde = p_min_criterium
        end
        
        % berekeningen voor Channel_Coding
        function [p_c, p_e, p_f, p_m, p_c_aprx, p_e_aprx, p_f_aprx, p_m_aprx] = Calc_prob_det
            % syndroomtabel
            persistent S_voll;
            S_voll = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
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
            persistent S;
            S = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
             
            % bitovergangswaarschijnlijkheid
            persistent p;
            p = 0.05;
            
            % array van alle codewoorden
            persistent codewords;
            codewords = Simulations.Calc_codewords;
            
            % p_c (kans op correcte decodering)
            p_c = 0;
            for i = 1:1
                w = sum(S(i,:));
                p_c = p_c + p^w * (1-p)^(14-w);
            end
            
            % benadering p_c (kans op correcte decodering)
            p_c_aprx = 1 - 14 * p;
            
            % p_e (kans op fout)
            p_e = 0;
            for j=2:size(codewords,1) % alle codewoorden behalve nulwoord
                for k=1:size(S,1) % alle syndromen
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
            p_f = 1 - p_c - p_m;
            
            % benaderde p_f (kans op decodeerfalen)
            p_f_aprx = 1 - p_c_aprx - p_m_aprx;
            
            % ontwerpcriterium
            f = @Calc_crit;
            x_0 = 0.0034;
            p_min_criterium = vpa(fzero(f,x_0));
            
            % resultaten
            probabiliteiten = [["ZUIVERE FOUTDETECTIE";"p_c"; "p_c ben."; "p_e"; "p_e ben."; "p_f"; "p_f ben."; "p_m"; "p_m ben.";], ["";p_c; p_c_aprx; p_e; p_e_aprx; p_f; p_f_aprx; p_m; p_m_aprx]]        
            %criterium_waarde = p_min_criterium
        end
        
        % exacte berekening voor lineaire blokcode met checkmatrix
        function p_e = Calc_outer(p)
            % INPUT
             % p : kans op bitovergang
            % OUTPUT
             % p_e : kans op decodeerfout
            
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
            p_e = 1 - p_c;
        end
        
        % simulatie voor lineaire blokcode met checkmatrix
        function p_e = Sim_outer(p)
            % INPUT
             % p : kans op bitovergang
            % OUTPUT
             % p_e: kans op decodeerfout

            N = 10000; % testgrootte
            
            % random input string
            b_enc = str2double(num2cell(dec2bin(floor(rand(1,N)*1023),10))); % "kolomvector" van random woorden
            b = reshape(b_enc',[], N*10); % "rijvector" van random woorden
            
            % encoderen van te verzenden woord
            c = Channel_Coding.Encode_outer(b);
           
            % foutpatroon afkomsting van BSC kanaal met 
            % bitovergangwaarschijnlijkheid p
            e = zeros(1, N*14);
            e(randperm(numel(e), round((14*N)*p))) = 1;
            
            % modellering transmissie door modulo-2 som van codewoord met 
            % foutpatroon
            r = mod(c+e,2);
            
            % decoderen van ontvangen woord
            [b_rec, ~] = Channel_Coding.Decode_outer(r);
            b_dec = reshape(b_rec, [10,N])';
            
            % optellen aantal fouten
            amount_errors = sum(sum(mod(b_enc - b_dec,2), 2) > 0);
            p_e = amount_errors/N;
        end
        
        % vergelijken analytische en experimentele kans op decodeerfout
        % voor lineaire blokcode met checkmatrix
        function Error_outer
            p = [0.001:0.001:0.5];
            
            p_e_ana = zeros(size(p,2),1);
            p_e_sim = zeros(size(p,2),1);
            
            for i = 1:size(p,2)
                p(i)
                p_e_ana(i) = Simulations.Calc_outer(p(i));
                p_e_sim(i) = Simulations.Sim_outer(p(i));
            end
            %[["p","p_e_ana", "p_e_sim"];[p,p_e_ana,p_e_sim]]
            
            hold on
            plot_sim = loglog(p,p_e_sim, 'LineWidth', 1, 'Color', 'Red');
            plot_ana = loglog(p,p_e_ana, 'LineWidth', 2, 'Color', 'Blue');
            
            l = legend([plot_ana, plot_sim], '$p_{e}$ analytisch', '$p_{e}$ experimenteel', 'Location', 'northeast');
            set(l,'Interpreter','latex')
            
            t = title('Kans op een decodeerfout, $p_{e}$, in functie van bitovergangsswaarschijnlijkheid, $p$');
            set(t,'Interpreter','latex')
            
            set(gca, 'XScale', 'log');
            set(gca, 'YScale', 'log');
            xlabel('p', 'FontSize',15);
            y = ylabel('$p_{e}$', 'FontSize',15);
            set(y,'Interpreter','latex')
            
            axis([0.001 0.5 0 1]);
            hold off;
        end
        
        % simulatie voor lineaire blokcode met checkmatrix en retrasmissie
        function [p_e_ARQ, deb] = Sim_outer_ARQ(T_max)
            % INPUT
             % T_max : hoeveelste retransmissie wanneer volledige
             % decodering toegepast wordt
            % OUTPUT
             % p_e_ARQ: kans op decodeerfout
             
            p = 0.05; % bitovergangwaarschijnlijkhied
            N = 100000; % testgrootte
            N_retr = N; % aantal woorden in huidige (re)transmissie
            k = 10 * N; % aantal informatiebits (aantal woorden * 10)
            n = 14 * N; % aantal verzonden bits (minstens aantal woorden * 14)
            T = 0; % aantal retransmissies (minstens 1)
           
            b_rec = zeros(0,0);
            b_enc = zeros(0,0);
            b_dec = zeros(0,0);
            
            % random input string
            b_enc_retr = str2double(num2cell(dec2bin(floor(rand(1,N_retr)*1023),10))); % "kolomvector" van random woorden
            
            while (T <= T_max && N_retr > 0)
                b = reshape(b_enc_retr',[], N_retr*10); % "rijvector" van random woorden               
                
                % encoderen van te verzenden woord 
                c = Channel_Coding.Encode_outer(b);

                % foutpatroon afkomsting van BSC kanaal met 
                % bitovergangwaarschijnlijkheid p
                e = zeros(1, N_retr*14);
                e(randperm(numel(e), round((14*N_retr)*p))) = 1;

                % modellering transmissie door modulo-2 som van codewoord met 
                % foutpatroon
                r = mod(c+e,2);
            
                [b_rec, b_err] = Channel_Coding.Decode_outer(r);
                b_rec = reshape(b_rec, [10,N_retr])';
                
                if(T < T_max)
                    b_enc = [b_enc; b_enc_retr((b_err(1,:) == 0),:)];
                    b_dec = [b_dec; b_rec((b_err(1,:) == 0),:)];
                    b_enc_retr = b_rec((b_err(1,:) == 1),:) ;

                    N_retr = sum(b_err(1,:)); % N woorden foutief ontvangen dus retransmissie
                    n = n + 14 * N_retr; % per extra woord dat retransmissie vereist 14 bits extra
                    T = T + 1;
                else
                    b_enc = [b_enc; b_enc_retr];
                    b_dec = [b_dec; b_rec];
                    T = T_max + 1;
                end
            end
            
            % optellen aantal fouten
            amount_errors = sum(sum(mod(b_enc - b_dec,2), 2) > 0);
            p_e_ARQ = amount_errors/N;
            deb = k / n;
        end
        
        % exacte berekening voor lineaire blokcode met checkmatrix en retrasmissie
        function [p_e_ARQ] = Calc_outer_ARQ(T_max)
            % INPUT
             % T_max : hoeveelste retransmissie wanneer volledige
             % decodering toegepast wordt
            % OUTPUT
             % p_e_ARQ: kans op decodeerfout
            
            % probabiliteiten voor p = 0.05
            p_e = 0.15163; % bij volledige decodering
            p_m = 0.0023029; % onafhankelijk van decodering
            p_f = 0.51002; % bij zuiver foutdetectie
            
            p_e_ARQ = (p_m * (1 - (p_f^T_max))/(1 - (p_f))) + (p_f^(T_max)) * p_e;
        end
        
        % vergelijken analytische en experimentele kans op decodeerfout
        % voor lineaire blokcode met checkmatrix en retransmissie
        function Error_outer_ARQ(type)
            p_e_ARQ_ana = zeros(16,1);
            p_e_ARQ_sim = zeros(16,1);
            deb = zeros(16,1);
            T_max = [0:1:15]';
            
            for i=0:15
                p_e_ARQ_ana(i+1) = Simulations.Calc_outer_ARQ(i);
                [p_e_ARQ_sim(i+1),deb(i+1)] = Simulations.Sim_outer_ARQ(i);
            end
            
            % [["T_max", "p_e_ARQ_ana", "p_e_ARQ_sim"];[T_max,p_e_ARQ_ana,(p_e_ARQ_sim)]]
            
            hold on
            if(type == 'p')
                t = title('Kans op een decodeerfout, $p_{e}^{ARQ1}$, in functie van max aantal toegelaten retransmissies, $T_{max}$', 'FontSize', 16);
            
                plot_sim = plot(T_max,p_e_ARQ_sim, 'LineWidth', 1, 'Color', 'Red');
                plot_ana = plot(T_max,p_e_ARQ_ana, 'LineWidth', 2, 'Color', 'Blue');
               
                l = legend([plot_ana, plot_sim], '$p_{e}^{ARQ1}$ analytisch', '$p_{e}^{ARQ1}$ experimenteel', 'Location', 'northeast');
                set(l,'Interpreter','latex')
            
                y = ylabel('$p_{e}^{ARQ1}$', 'FontSize',15);
                axis([0 15 0 0.20]);
            elseif(type == 'd')
                t = title('Informatiedebiet, $\frac{k}{n}$, in functie van max aantal toegelaten retransmissies, $T_{max}$', 'FontSize', 16);
            
                plot_deb = plot(T_max,p_e_ARQ_sim, 'LineWidth', 1, 'Color', 'Red');
               
                l = legend([plot_deb], 'informatiedebiet $\frac{k}{n}$', 'Location', 'northeast');
                set(l,'Interpreter','latex')
            
                y = ylabel('$\frac{k}{n}$', 'FontSize',15);
                %axis([0 15 0 0.20]);
            end
            set(t,'Interpreter','latex');
            
            x = xlabel('$T_{max}$', 'FontSize', 15);
            set(x,'Interpreter','latex');
            
            set(gca, 'YScale', 'log');
            set(y,'Interpreter','latex');
            
            hold off;
        end
        
        function [p_e_ARQ, deb] = Sim_inner_5(T_max)
            % INPUT
             % T_max : hoeveelste retransmissie wanneer volledige
             % decodering toegepast wordt
            % OUTPUT
             % p_e_ARQ: kans op decodeerfout
            
            p = 0.05; % bitovergangwaarschijnlijkhied
            g_x = [1, 1, 0, 1, 0, 1];
            N = 500000; % testgrootte
            N_retr = N; % aantal woorden in huidige (re)transmissie
            k = 5 * N; % aantal informatiebits (aantal woorden * 5)
            n = 14 * N; % aantal verzonden bits (minstens aantal woorden * 14)
            T = 0; % aantal retransmissies (minstens 1)
           
            b_rec = zeros(0,5);
            b_enc = zeros(0,5);
            b_dec = zeros(0,5);
            amount_errors = 0;
            
            % random input string
            b_enc_retr = str2double(num2cell(dec2bin(floor(rand(1,N_retr)*31),5))); % "kolomvector" van random woorden
            
            
            while (T <= T_max && N_retr > 0)
                % inner encoderen van te verzenden woord 
                c_i = cell2mat(arrayfun(@(i) Channel_Coding.Encode_inner(b_enc_retr(i,:), g_x), 1:size(b_enc_retr,1),'un',0)');
                % outer encoderen van te verzenden woord 
                c_o = Channel_Coding.Encode_outer(reshape(c_i',[], N_retr*10));

                % foutpatroon afkomsting van BSC kanaal met 
                % bitovergangwaarschijnlijkheid p
                e = zeros(1, N_retr*14);
                e(randperm(numel(e), round((14*N_retr)*p))) = 1;
                
                % modellering transmissie door modulo-2 som van codewoord met 
                % foutpatroon
                r = mod(c_o+e,2);
            
                [b_rec_o, ~] = Channel_Coding.Decode_outer(r);
                b_rec_o = reshape(b_rec_o, [10,N_retr])';
                
                [b_rec_i, b_err_i] = (arrayfun(@(i) Channel_Coding.Decode_inner(b_rec_o(i,:), g_x), 1:size(b_rec_o,1),'un',0));
                b_rec_i = cell2mat(b_rec_i');
                b_err_i = cell2mat(b_err_i);
                b_enc = [b_enc; b_enc_retr((b_err_i(1,:) == 0),:)];
                b_dec = [b_dec; b_rec_i((b_err_i(1,:) == 0),:)];
                b_enc_retr = b_rec_i((b_err_i(1,:) == 1),:) ;

                if(T < T_max)
                    N_retr = sum(b_err_i(1,:)); % N woorden foutief ontvangen dus retransmissie
                    n = n + 14 * N_retr; % per extra woord dat retransmissie vereist 14 bits extra
                    T = T + 1;
                else
                    amount_errors = amount_errors + sum(b_err_i(1,:));
                    T = T_max + 1;
                end
            end
            
            % optellen aantal fouten
            amount_errors = amount_errors + sum(sum(mod(b_enc - b_dec,2), 1) > 0);
            p_e_ARQ = amount_errors/N;
            deb = k / n;
        end
        
        function Error_inner_5(type)
            % INPUT
             % type ('p' voor p_e_ARQ, 'd' voor debiet)
            
            p_e_ARQ_2a_sim = zeros(7,1);
            p_e_ARQ_2a_ana = zeros(7,1);
            debiet_sim = zeros(7,1);
            T_max = [0:6]';
            for i = 0:size(T_max,1)-1
               [p_e_ARQ_2a_sim(i+1,1), debiet_sim(i+1,1)] = Simulations.Sim_inner_5(i);
               p_e_ARQ_2a_ana(i+1,1) = 0.15163^(i+1);
            end
            
            [["T_max", "p_e_ARQ_ana", "p_e_ARQ_sim", "debiet"]; [T_max,p_e_ARQ_2a_ana,p_e_ARQ_2a_sim, debiet_sim]]
            
            hold on
            
            if(type == 'p')
                t = title('Kans op een decodeerfout, $p_{e}^{ARQ2a}$, in functie van max aantal toegelaten retransmissies, $T_{max}$', 'FontSize', 16);
            
                plot_sim = plot(T_max,p_e_ARQ_2a_sim, 'LineWidth', 1, 'Color', 'Red');                
                plot_ana = plot(T_max,p_e_ARQ_2a_ana, 'LineWidth', 2, 'Color', 'Blue');
            
                l = legend([plot_ana, plot_sim], '$p_{e}^{ARQ2a}$ analytisch', '$p_{e}^{ARQ2a}$ experimenteel', 'Location', 'northeast');
                
                set(gca, 'YScale', 'log');
                y = ylabel('$p_{e}^{ARQ2a}$', 'FontSize',15);
                axis([0 6 0 0.20]);
            elseif(type == 'd')
                t = title('Informatiedebiet, $\frac{k}{n}$, in functie van max aantal toegelaten retransmissies, $T_{max}$', 'FontSize', 16);
            
                plot_deb = plot(T_max,debiet_sim, 'LineWidth', 1, 'Color','Blue');

                l = legend([plot_deb], 'Informatiedebiet $\frac{k}{n}$', 'Location', 'northeast');
                              
                set(gca, 'YScale', 'log')
                y = ylabel('$\frac{k}{n}$', 'FontSize',15);
                axis([0 6 0.28 0.37]);
            end
            set(l,'Interpreter','latex')
            set(t,'Interpreter','latex')
            
            x = xlabel('$T_{max}$', 'FontSize',15);
            set(x,'Interpreter','latex')
            set(y,'Interpreter','latex')
                       
            hold off;
        end
        
        function [p_e_ARQ, deb] = Sim_inner_8(T_max)
            % INPUT
             % T_max : hoeveelste retransmissie wanneer volledige
             % decodering toegepast wordt
            % OUTPUT
             % p_e_ARQ: kans op decodeerfout
            
            p = 0.05; % bitovergangwaarschijnlijkhied
            g_x = [1, 1, 0, 0, 1, 1, 0, 1, 1];
            N = 500000; % testgrootte
            N_retr = N; % aantal woorden in huidige (re)transmissie
            k = 2 * N; % aantal informatiebits (aantal woorden * 2)
            n = 14 * N; % aantal verzonden bits (minstens aantal woorden * 14)
            T = 0; % aantal retransmissies (minstens 1)
           
            b_rec = zeros(0,2);
            b_enc = zeros(0,2);
            b_dec = zeros(0,2);
            amount_errors = 0;
            
            % random input string
            b_enc_retr = str2double(num2cell(dec2bin(floor(rand(1,N_retr)*3),2))); % "kolomvector" van random woorden
            
            
            while (T <= T_max && N_retr > 0)
                % inner encoderen van te verzenden woord 
                c_i = cell2mat(arrayfun(@(i) Channel_Coding.Encode_inner(b_enc_retr(i,:), g_x), 1:size(b_enc_retr,1),'un',0)');
                % outer encoderen van te verzenden woord 
                c_o = Channel_Coding.Encode_outer(reshape(c_i',[], N_retr*10));

                % foutpatroon afkomsting van BSC kanaal met 
                % bitovergangwaarschijnlijkheid p
                e = zeros(1, N_retr*14);
                e(randperm(numel(e), round((14*N_retr)*p))) = 1;

                % modellering transmissie door modulo-2 som van codewoord met 
                % foutpatroon
                r = mod(c_o+e,2);
            
                [b_rec_o, ~] = Channel_Coding.Decode_outer(r);
                b_rec_o = reshape(b_rec_o, [10,N_retr])';
                
                [b_rec_i, b_err_i] = (arrayfun(@(i) Channel_Coding.Decode_inner(b_rec_o(i,:), g_x), 1:size(b_rec_o,1),'un',0));
                b_rec_i = cell2mat(b_rec_i');
                b_err_i = cell2mat(b_err_i);
                b_enc = [b_enc; b_enc_retr((b_err_i(1,:) == 0),:)];
                b_dec = [b_dec; b_rec_i((b_err_i(1,:) == 0),:)];
                b_enc_retr = b_rec_i((b_err_i(1,:) == 1),:) ;

                if(T < T_max)
                    N_retr = sum(b_err_i(1,:)); % N woorden foutief ontvangen dus retransmissie
                    n = n + 14 * N_retr; % per extra woord dat retransmissie vereist 14 bits extra
                    T = T + 1;
                else
                    amount_errors = amount_errors + sum(b_err_i(1,:));
                    T = T_max + 1;
                end
            end
            
            % optellen aantal fouten
            amount_errors = amount_errors + sum(sum(mod(b_enc - b_dec,2), 1) > 0);
            p_e_ARQ = amount_errors/N;
            deb = k / n;
        end
        
        function Error_inner_8(type)
        % INPUT
         % type ('p' voor p_e_ARQ, 'd' voor debiet)
         
            p_e_ARQ_2b_sim = zeros(7,1);
            p_e_ARQ_2b_ana = zeros(7,1);
            debiet_sim = zeros(7,1);
            T_max = [0:6]';
            for i = 0:size(T_max,1)-1
               [p_e_ARQ_2b_sim(i+1,1), debiet_sim(i+1,1)] = Simulations.Sim_inner_8(i);
               p_e_ARQ_2b_ana(i+1,1) = 0.15163^(i+1);
            end
            
            [["T_max", "p_e_ARQ_ana", "p_e_ARQ_sim", "debiet"]; [T_max,p_e_ARQ_2b_ana,p_e_ARQ_2b_sim, debiet_sim]]
            
            hold on
 
            if(type == 'p')
                t = title('Kans op een decodeerfout, $p_{e}^{ARQ2b}$, in functie van max aantal toegelaten retransmissies, $T_{max}$', 'FontSize', 16);
                
                plot_sim = plot(T_max,p_e_ARQ_2b_sim, 'LineWidth', 1, 'Color', 'Red');
                plot_ana = plot(T_max,p_e_ARQ_2b_ana, 'LineWidth', 2, 'Color', 'Blue');
                
                l = legend([plot_ana, plot_sim], '$p_{e}^{ARQ2b}$ analytisch', '$p_{e}^{ARQ2b}$ experimenteel', 'Location', 'northeast');
                
                set(gca, 'YScale', 'log');
                y = ylabel('$p_{e}^{ARQ2b}$', 'FontSize',15);
                axis([0 6 0 0.20]);
            elseif(type == 'd')
                t = title('Informatiedebiet, $\frac{k}{n}$, in functie van max aantal toegelaten retransmissies, $T_{max}$', 'FontSize', 16);
            
                plot_deb = plot(T_max,debiet_sim, 'LineWidth', 1, 'Color','Blue');

                l = legend([plot_deb], 'Informatiedebiet $\frac{k}{n}$', 'Location', 'northeast');
                
                y = ylabel('$\frac{k}{n}$', 'FontSize',15);
                axis([0 6 0.09 0.16]);
            end
            
            set(l,'Interpreter','latex')
            set(t,'Interpreter','latex')
            
            x = xlabel('$T_{max}$', 'FontSize',15);
            set(x,'Interpreter','latex')
            set(y,'Interpreter','latex')
                       
            hold off;
        end
        
        function Error_inner_compare(type)
            % INPUT
             % type ('p' voor p_e_ARQ, 'd' voor debiet)
            
            p_e_ARQ_1 = zeros(7,1);
            p_e_ARQ_2a = zeros(7,1);
            p_e_ARQ_2b = zeros(7,1);
            
            deb_a = zeros(7,1);
            deb_b = zeros(7,1);
            
            T_max = [0:6]';
            
            for i = 0:size(T_max,1)-1
               [p_e_ARQ_2a(i+1,1), deb_a(i+1,1)] = Simulations.Sim_inner_5(i);
               [p_e_ARQ_2b(i+1,1), deb_b(i+1,1)] = Simulations.Sim_inner_8(i);
               p_e_ARQ_1(i+1,1) = 0.15163^(i+1);
            end
            
            [["T_max", "p_e_ARQ_2a", "p_e_ARQ_2b", "deb_a", "deb_b"]; [T_max,p_e_ARQ_2a,p_e_ARQ_2b, deb_a, deb_b]]
            
             hold on
 
            if(type == 'p')
                t = title('Kans op een decodeerfout, $p_{e}^{ARQ2b}$, in functie van max aantal toegelaten retransmissies, $T_{max}$', 'FontSize', 16);
                
                plot_2a = plot(T_max,p_e_ARQ_2a, 'LineWidth', 1, 'Color', 'Blue');
                plot_2b = plot(T_max,p_e_ARQ_2b, 'LineWidth', 1, 'Color', 'Red');
                plot_1 = plot(T_max, p_e_ARQ_1, 'LineWidth', 2, 'Color', 'Green');
                
                l = legend([plot_1,plot_2a, plot_2b], '$p_{e}^{ARQ}$ analytisch','$p_{e}^{ARQ2a}$ experimenteel', '$p_{e}^{ARQ2b}$ experimenteel', 'Location', 'northeast');
                
                set(gca, 'YScale', 'log');
                y = ylabel('$p_{e}^{ARQ}$', 'FontSize',15);
                axis([0 6 0 0.20]);
            elseif(type == 'd')
                t = title('Informatiedebiet, $\frac{k}{n}$, in functie van max aantal toegelaten retransmissies, $T_{max}$', 'FontSize', 16);
            
                plot_deb_a = plot(T_max,deb_a, 'LineWidth', 1, 'Color','Blue');
                plot_deb_b = plot(T_max,deb_b, 'LineWidth', 1, 'Color','Red');

                l = legend([plot_deb_a, plot_deb_b], 'Informatiedebiet $(\frac{k}{n})^{2a}$', 'Informatiedebiet $(\frac{k}{n})^{2b}$', 'Location', 'northeast');
                
                y = ylabel('$\frac{k}{n}$', 'FontSize',15);
                axis([0 6 0.09 0.37]);
            end
            
            set(l,'Interpreter','latex')
            set(t,'Interpreter','latex')
            
            x = xlabel('$T_{max}$', 'FontSize',15);
            set(x,'Interpreter','latex')
            set(y,'Interpreter','latex')
                       
            hold off;
        end
    end
end

