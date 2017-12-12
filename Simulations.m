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
            % p_m_aprx = 0;
            % for i=3:12
            %    w = i;
            %    for j = 0:(14-w)
            %        p_m_aprx = p_m_aprx + codewords_w(i)*p^i*(factorial(14-w)/factorial(14-w-j)/factorial(j)*(-p)^j);
            %    end
            % end
                     
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
            p_e = 1-p_c;
        end
        
        % simulatie voor lineaire blokcode met checkmatrix
        function p_e = Sim_outer(p)
            % INPUT
             % p : kans op bitovergang
            % OUTPUT
             % p_e: kans op decodeerfout

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
           
            % modellering van BSC kanaal met bitovergangwaarschijnlijkheid
            % p
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
        
        % vergelijken analytische en experimentele kans op decodeerfout
        % voor lineaire blokcode met checkmatrix
        function Error_outer
            N = 0.5/0.001;
            persistent p_e_ana;
            persistent p_e_sim;
            p_e_ana = zeros(N,1);
            p_e_sim = zeros(N,1);
            p = [0.001:0.001:0.5];
            for i = 1:N
                p_e_ana(i) = vpa(Simulations.Calc_outer(p(i)));
                p_e_sim(i) = vpa(Simulations.Sim_outer(p(i)));
            end
            %[["p","p_e_ana", "p_e_sim"];[p,p_e_ana,p_e_sim]]
            
            hold on
            plot_sim = loglog(p,p_e_sim, 'LineWidth', 1, 'Color', 'Red');
            plot_ana = loglog(p,p_e_ana, 'LineWidth', 2, 'Color', 'Blue');
            
            l = legend([plot_ana, plot_sim], '$p_{e}$ analytisch', '$p_{e}$ experimenteel', 'Location', 'northeast');
            set(l,'Interpreter','latex')
            
            t = title('Kans op een decodeerfout, $p_{e}$, in functie van bitovergangsswaarschijnlijkheid, $p$');
            set(t,'Interpreter','latex')
            
            xlabel('p', 'FontSize',12);
            y = ylabel('$p_{e}$', 'FontSize',12);
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
            N = 1000; % testgrootte
            N_retr = N; % aantal woorden in huidige (re)transmissie
            k = 10 * N; % aantal informatiebits (aantal woorden * 10)
            n = 14 * N; % aantal verzonden bits (minstens aantal woorden * 14)
            T = 0; % aantal retransmissies (minstens 1)
           
            % random input string
            b = zeros(1,N_retr*10);
            for i=1:N_retr
                b_rand = dec2bin((floor(rand(1).*1023)));
                for j=(10-length(b_rand)+1):10
                    b(1, (i-1)*10+j) = str2num(b_rand(j-10+length(b_rand)));
                end
            end
            
            % bij elke transmissie oorspronkelijke en correcte/niet opgemerkte 
            % foute woorden hierin wegschrijven ter controle
            b_sent = zeros(0,0); % alle verstuurde/geencodeerde woorden
            b_rec = zeros(0,0); % alle ontvangen/gedecodeerde woorden
                
            while (T <= T_max && N_retr > 0)
                % encoderen van woord
                c = Channel_Coding.Encode_outer(b);

                % modellering van BSC kanaal met bitovergangwaarschijnlijkheid
                % p
                r = c;
                for j=1:(N_retr*14)
                    p_rand = rand(1);
                    if(p_rand < p)
                        if(r(1,j) == 0)
                            r(1,j) = 1;
                        else
                            r(1,j) = 0;
                        end
                    end
                end
                
                [b_dec, b_err] = Channel_Coding.Decode_outer(r);
                
                if(T < T_max)
                    b_retr = zeros(0,0);
                    for l=1:N_retr
                        if(b_err(1,l) == 1) % FOUT? Retransmissie
                            b_retr = [b_retr, b(10*(l-1)+1:10*l)];
                        else % GEEN FOUT? Wegschrijven want 'aanvaard' door decoder
                            b_sent = [b_sent, b(10*(l-1)+1:10*l)];
                            b_rec = [b_rec, b_dec(10*(l-1)+1:10*l)];
                        end
                    end
                    b = b_retr;
                    
                    N_retr = sum(b_err(1,:)); % N woorden foutief ontvangen dus retransmissie
                    n = n + 14 * N_retr; % per extra woord dat retransmissie vereist 14 bits extra
                    T = T + 1;
                else
                    for l=1:N_retr
                        b_sent = [b_sent, b(10*(l-1)+1:10*l)];
                        b_rec = [b_rec, b_dec(10*(l-1)+1:10*l)];
                    end
                    T = T_max + 1;
                end
            end
            
            % verstuurde en ontvangen woorden vergelijken
            controle = zeros(1,N_retr); % rij met 1 op i-de plaats als na T_max 
                                        % transmissies i-de woord in b_rec niet
                                        % overeenkomt met i-de woord in
                                        % b_sent
            for l=1:N
                if (isequal(b_sent(10*(l-1)+1:10*(l)),b_rec(10*(l-1)+1:10*(l))))
                    controle(l) = 0;
                else
                    controle(l) = 1;
                end
            end
            
            deb = k / n;
            p_e_ARQ = sum(controle(1,:))/N;
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
            
            p_e_ARQ = vpa(p_m * (1 - (p_f^T_max))/(1 - (p_f))) + (p_f^(T_max)) * p_e;
        end
        
        % vergelijken analytische en experimentele kans op decodeerfout
        % voor lineaire blokcode met checkmatrix en retransmissie
        function Error_outer_ARQ
            p_e_ARQ_ana = zeros(16,1);
            p_e_ARQ_sim = zeros(16,1);
            T_max = [0:1:15]';
            
            for i=0:15
                p_e_ARQ_ana(i+1) = vpa(Simulations.Calc_outer_ARQ(i));
                p_e_ARQ_sim(i+1) = vpa(Simulations.Sim_outer_ARQ(i));
            end
            
            % [["T_max", "p_e_ARQ_ana", "p_e_ARQ_sim"];[T_max,p_e_ARQ_ana,(p_e_ARQ_sim)]]
            
            hold on
            %plot_sim = loglog(T_max,p_e_ARQ_sim, 'LineWidth', 1, 'Color', 'Red');
            %plot_ana = loglog(T_max,p_e_ARQ_ana, 'LineWidth', 2, 'Color', 'Blue');
            plot_sim = plot(T_max,p_e_ARQ_sim, 'LineWidth', 1, 'Color', 'Red');
            plot_ana = plot(T_max,p_e_ARQ_ana, 'LineWidth', 2, 'Color', 'Blue');
            
            l = legend([plot_ana, plot_sim], '$p_{e}^{ARQ1}$ analytisch', '$p_{e}^{ARQ1}$ experimenteel', 'Location', 'northeast');
            set(l,'Interpreter','latex')
            
            t = title('Kans op een decodeerfout, $p_{e}^{ARQ1}$, in functie van max aantal toegelaten retransmissies, $T_{max}$', 'FontSize', 16);
            set(t,'Interpreter','latex')
            
            xlabel('T_{max}', 'FontSize',15);
            y = ylabel('$p_{e}^{ARQ1}$', 'FontSize',15);
            set(y,'Interpreter','latex')
            
            axis([0 15 0 0.20]);
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
            N = 2000; % testgrootte
            N_retr = N; % aantal woorden in huidige (re)transmissie
            k = 5 * N; % aantal informatiebits (aantal woorden * 10)
            n = 0; % aantal verzonden bits (minstens aantal woorden * 14)
            T = 0; % aantal retransmissies (minstens 1)
            controle = zeros(1,N_retr); % rij met 1 op i-de plaats als na T_max 
                                        % transmissies i-de ontvangen woord niet
                                        % overeenkomt met verstuurde woord
                                        
            % bij elke transmissie oorspronkelijke en correcte/niet opgemerkte 
            % foute woorden hierin wegschrijven ter controle
            b_sent = zeros(0,0); % alle verstuurde/geencodeerde woorden
            b_rec = zeros(0,0); % alle ontvangen/gedecodeerde woorden
            
            % random input string
            b = zeros(1,N_retr*5);
            
            for i=1:N_retr
                b_rand = dec2bin((floor(rand(1).*31)));
                for j=(5-length(b_rand)+1):5
                    b(1, (i-1)*5+j) = str2num(b_rand(j-5+length(b_rand)));
                end
            end
                
            while (T <= T_max && N_retr > 0) 
                n = n + 14 * N_retr; % per extra woord dat retransmissie vereist 14 bits extra !!!
                
                % inner encoderen van woord
                c_i = zeros(1,N_retr*10);
                for i = 1:N_retr
                    c_i(10*(i-1)+1:10*i) = Channel_Coding.Encode_inner(b(5*(i-1)+1:5*i), g_x);
                end
                
                % outer encoderen van woord
                c_o = Channel_Coding.Encode_outer(c_i);
                
                % modellering van BSC kanaal met bitovergangwaarschijnlijkheid
                % p
                r = c_o;
                for j=1:(N_retr*14)
                    p_rand = rand(1);
                    if(p_rand < p)
                        if(r(1,j) == 0)
                            r(1,j) = 1;
                        else
                            r(1,j) = 0;
                        end
                    end
                end
                
                b_rec_o = zeros(N_retr*14);
                [b_rec_o, ~] = Channel_Coding.Decode_outer(r);
                
                b_retr = zeros(0,0); % lijst van te retransmitten woorden
                N_retr_prev = N_retr;
                N_retr = 0;
                for i = 1:N_retr_prev
                    [b_rec_i, b_err_i] = Channel_Coding.Decode_inner(b_rec_o(10*(i-1)+1:10*i), g_x); % i-de gedecodeerde woord en foutpatroon
                    if(b_err_i == 1)
                        b_retr = [b_retr, b(5*(i-1)+1:5*i)]; 
                        N_retr = N_retr + 1;
                    else
                        b_sent = [b_sent, b(5*(i-1)+1:5*i)];
                        b_rec = [b_rec, b_rec_i];
                    end
                end
                b = b_retr;
                T = T + 1;
            end
            if(N_retr ~= 0)
                for i = 1:N_retr
                    [b_rec_i, ~] = Channel_Coding.Decode_inner(b_rec_o(10*(i-1)+1:10*i), g_x); % i-de gedecodeerde woord en foutpatroon
                    b_sent = [b_sent, b(5*(i-1)+1:5*i)];
                    b_rec = [b_rec, b_rec_i];
                end
            end
            
            for l=1:N
                if (isequal(b_sent(5*(l-1)+1:5*(l)),b_rec(5*(l-1)+1:5*(l))))
                    controle(l) = 0;
                else
                    controle(l) = 1;
                end
            end
                   
            deb = k / n;
            p_e_ARQ = vpa(sum(controle(1,:))/N);
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
            
            %[["T_max", "p_e_ARQ_ana", "p_e_ARQ_sim", "debiet"]; [T_max,p_e_ARQ_2a_ana,p_e_ARQ_2a_sim, debiet_sim]]
            
            hold on
            
            if(type == 'p')
                t = title('Kans op een decodeerfout, $p_{e}^{ARQ2a}$, in functie van max aantal toegelaten retransmissies, $T_{max}$', 'FontSize', 16);
            
                plot_sim = plot(T_max,p_e_ARQ_2a_sim, 'LineWidth', 1, 'Color', 'Red');                
                plot_ana = plot(T_max,p_e_ARQ_2a_ana, 'LineWidth', 2, 'Color', 'Blue');
            
                l = legend([plot_ana, plot_sim], '$p_{e}^{ARQ2a}$ analytisch', '$p_{e}^{ARQ2a}$ experimenteel', 'Location', 'northeast');
            
                y = ylabel('$p_{e}^{ARQ2a}$', 'FontSize',15);
                axis([0 6 0 0.20]);
            elseif(type == 'd')
                t = title('Informatiedebiet, $\frac{k}{n}$, in functie van max aantal toegelaten retransmissies, $T_{max}$', 'FontSize', 16);
            
                plot_deb = plot(T_max,debiet_sim, 'LineWidth', 1, 'Color','Blue');

                l = legend([plot_deb], 'Informatiedebiet $\frac{k}{n}$', 'Location', 'northeast');
                
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
            N = 2000; % testgrootte
            N_retr = N; % aantal woorden in huidige (re)transmissie
            k = 2 * N; % aantal informatiebits (aantal woorden * 10)
            n = 0; % aantal verzonden bits (minstens aantal woorden * 14)
            T = 0; % aantal retransmissies (minstens 1)
            controle = zeros(1,N_retr); % rij met 1 op i-de plaats als na T_max 
                                        % transmissies i-de ontvangen woord niet
                                        % overeenkomt met verstuurde woord
                                        
            % bij elke transmissie oorspronkelijke en correcte/niet opgemerkte 
            % foute woorden hierin wegschrijven ter controle
            b_sent = zeros(0,0); % alle verstuurde/geencodeerde woorden
            b_rec = zeros(0,0); % alle ontvangen/gedecodeerde woorden
            
            % random input string
            b = zeros(1,N_retr*2);
            
            for i=1:N_retr
                b_rand = dec2bin((floor(rand(1).*3)));
                for j=(2-length(b_rand)+1):2
                    b(1, (i-1)*2+j) = str2num(b_rand(j-2+length(b_rand)));
                end
            end
                
            while (T <= T_max && N_retr > 0) 
                n = n + 14 * N_retr; % per extra woord dat retransmissie vereist 14 bits extra !!!
                
                % inner encoderen van woord
                c_i = zeros(1,N_retr*10);
                for i = 1:N_retr
                    c_i(10*(i-1)+1:10*i) = Channel_Coding.Encode_inner(b(2*(i-1)+1:2*i), g_x);
                end
                
                % outer encoderen van woord
                c_o = Channel_Coding.Encode_outer(c_i);
                
                % modellering van BSC kanaal met bitovergangwaarschijnlijkheid
                % p
                r = c_o;
                for j=1:(N_retr*14)
                    p_rand = rand(1);
                    if(p_rand < p)
                        if(r(1,j) == 0)
                            r(1,j) = 1;
                        else
                            r(1,j) = 0;
                        end
                    end
                end
                
                b_rec_o = zeros(N_retr*14);
                [b_rec_o, ~] = Channel_Coding.Decode_outer(r);
                
                b_retr = zeros(0,0); % lijst van te retransmitten woorden
                N_retr_prev = N_retr;
                N_retr = 0;
                for i = 1:N_retr_prev
                    [b_rec_i, b_err_i] = Channel_Coding.Decode_inner(b_rec_o(10*(i-1)+1:10*i), g_x); % i-de gedecodeerde woord en foutpatroon
                    if(b_err_i == 1)
                        b_retr = [b_retr, b(2*(i-1)+1:2*i)]; 
                        N_retr = N_retr + 1;
                    else
                        b_sent = [b_sent, b(2*(i-1)+1:2*i)];
                        b_rec = [b_rec, b_rec_i];
                    end
                end
                b = b_retr;
                T = T + 1;
            end
            if(N_retr ~= 0)
                for i = 1:N_retr
                    [b_rec_i, ~] = Channel_Coding.Decode_inner(b_rec_o(10*(i-1)+1:10*i), g_x); % i-de gedecodeerde woord en foutpatroon
                    b_sent = [b_sent, b(2*(i-1)+1:2*i)];
                    b_rec = [b_rec, b_rec_i];
                end
            end
            
            for l=1:N
                if (isequal(b_sent(2*(l-1)+1:2*(l)),b_rec(2*(l-1)+1:2*(l))))
                    controle(l) = 0;
                else
                    controle(l) = 1;
                end
            end
                   
            deb = k / n;
            p_e_ARQ = vpa(sum(controle(1,:))/N);
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
            
            %[["T_max", "p_e_ARQ_ana", "p_e_ARQ_sim", "debiet"]; [T_max,p_e_ARQ_2b_ana,p_e_ARQ_2b_sim, debiet_sim]]
            
            hold on
 
            if(type == 'p')
                t = title('Kans op een decodeerfout, $p_{e}^{ARQ2b}$, in functie van max aantal toegelaten retransmissies, $T_{max}$', 'FontSize', 16);
                
                plot_sim = plot(T_max,p_e_ARQ_2b_sim, 'LineWidth', 1, 'Color', 'Red');
                plot_ana = plot(T_max,p_e_ARQ_2b_ana, 'LineWidth', 2, 'Color', 'Blue');
                
                l = legend([plot_ana, plot_sim], '$p_{e}^{ARQ2b}$ analytisch', '$p_{e}^{ARQ2b}$ experimenteel', 'Location', 'northeast');
                
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

