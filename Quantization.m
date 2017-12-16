classdef Quantization
    %klasse om de figuur te kwantiseren;
    
    properties(Constant)
        filename = 'ProjectA.mat'; %vul hier het pad naar de file met de figuur toe. (.mat bestand)
    end
    
    methods(Static=true)
        
        % Functie om distributie en genormalizeerd histogram te plotten
        function plot_distribution()
            % samples bevat de 8-bit pixelwaarden van de originele figuur
            % distr is een anonieme functie die f_U(u) voorstelt.
            [samples, distr] = make_probability_functions(Quantization.filename);

            h = histogram(samples, 'Normalization', 'pdf', 'FaceColor', 'Yellow');
            hold on
            p = plot(0:255, distr(0:255), 'LineWidth', 2, 'Color', 'Blue');

            legend([p, h], 'p.d.f', 'hist', 'Location', 'northeast');
            title('Plot of the distribution and a histogram of the figure');
            xlabel('8-bit pixelvalue');
            ylabel('Probablity');
            axis([0 255 0 0.016]);
            hold off;
        end
        
        % Functie om de optimale uniforme kwantisator te bepalen
        function [Delta_opt,GKD,SQR,entropie,r,q,p]= optimal_linear_quantizer()
            % OUTPUT
              % Delta_opt : optimale stapgrootte
              % GKD : GKD uniforme kwantisator
              % SQR : SQR uniforme kwantisator
              % entropie : entropie van het gekwantiseerde signaal
              % r : kwantisatiedrempels
              % q : kwantisatieniveaus
              % p : relatieve frequentie kwantisatieniveaus
              % baseer je voor alle berekeningen op de anonieme functie distr.
            [~, distr] = make_probability_functions(Quantization.filename);
            x_0 = 127.5;
            M = 8;
            Delta = [5:0.1:30]; % 1 X 251
            
            q = x_0 + ([1:M] - (M+1)/2)' * Delta;
            sigma_e_gr_sqr = sum(arrayfun(@(q_i_j, lo, hi) integral(@(u) (q_i_j-u).^2.*distr(u), lo, hi), q, q-Delta/2, q+Delta/2));
            sigma_e_ol_sqr = arrayfun(@(q1, qM, q1_min, qM_plus)...
                integral(@(u) (q1-u).^2.*distr(u), 0, q1_min) + integral(@(u) (qM-u).^2.*distr(u), qM_plus, 255), ...
                q(1,:), q(M,:), q(1,:)-Delta/2, q(M,:)+Delta/2);
            sigma_e_sqr = sigma_e_gr_sqr + sigma_e_ol_sqr;
             
            [GKD, min_index] = min(sigma_e_sqr(:));
            Delta_opt = Delta(min_index);
            gem=integral(@(u) u.*distr(u), 0, 255);
            var=integral(@(u) distr(u).*(u-gem).^2, 0, 255);
            SQR=10*log10(var/GKD);
            q = q(:, min_index)';
            r = x_0 + (2*[1:M-1] - M)/2 .* Delta_opt;
            r = [0 r 255];
            p = arrayfun(@(i) integral(@(x) distr(x), r(i), r(i+1)), 1:8);
            entropie = -sum(p .* log2(p));
            %{
            subplot(1,2,1)
            gkd_plot = plot(Delta, sigma_e_sqr, 'r', 'LineWidth', 1);
            hold on
            gr_plot = plot(Delta, sigma_e_gr_sqr, 'y--', 'LineWidth', 1);
            ol_plot = plot(Delta, sigma_e_ol_sqr, 'b:', 'LineWidth', 1.5);
            min_marker = plot(Delta_opt, GKD, 'gO', 'LineWidth', 3);
            title('GKD, granulaire bijdrage en overlaadbijdrage');
            xlabel('\Delta');
            ylabel('GKD');
            legend([gkd_plot, gr_plot, ol_plot, min_marker], 'GKD', 'Granulaire bijdrage', 'Overlaadbijdrage', 'Optimale stapgrootte');
            hold off
            
            subplot(1,2,2);
            
            pdf = plot(0:255, distr(0:255), 'LineWidth', 2, 'Color', 'Blue');
            hold on
            q_points = stem(q, distr(q), 'r*', 'LineWidth', 2, 'Color', 'Green');
            r_points = stem(r, distr(r), 'g*', 'LineWidth', 2, 'Color', 'Red');
            xlabel('u');
            ylabel('f_U(u)');
            legend([pdf, q_points, r_points], '$f_{U}(u)$', 'reconstructieniveaus', 'kwantisatiedrempels');
            title('Waarschijnlijkheid van reconstructieniveaus en kwantisatiedrempels');
   
            axis([0 255 0 0.016]);
            hold off
            %}
        end
        
        % Functie om Lloyd-Max kwantisator te bepalen
        function [GKD,SQR,entropie,r,q,cur_p]= Lloyd_max_quantizer()
            % OUTPUT
             % GKD : GKD Lloyd-Max kwantisator
             % SQR : SQR Lloyd-Max kwantisator
             % entropie : entropie van het gekwantiseerde signaal
             % r : kwantisatiedrempels
             % q : kwantisatieniveaus
             % p : relatieve frequentie kwantisatieniveaus
            % baseer je voor alle berekeningen op de anonieme functie distr.
            [~, distr] = make_probability_functions(Quantization.filename);
            M = 8;
            EPS = 1e-9;
            
            r = zeros(1,M+1);
            q = [255/(2*M) : 255/M : 255-255/(2*M)]; % optimale uniforme met bereik 0..255
            l = 1;
            r(1) = 0;
            r(M+1) = 255;
            prev_distortion = 0;
            cur_distortion = 0;
            all_distortions = [];
            all_entropie = [];
            while true
                
                for i = 1:M-1 % update alle kwantisatiedrempels
                    r(i+1) = (q(i) + q(i+1))/2;
                end
                
                for i = 1:M % update alle kwantisatieniveaus
                    teller = integral(@(u) u.*distr(u), r(i), r(i+1));
                    noemer = integral(@(u) distr(u), r(i), r(i+1));
                    if noemer~=0 
                        q(i) = teller/noemer;
                    end
                end
                
                prev_distortion = cur_distortion;
                cur_distortion = 0;
                for i = 1:M 
                    cur_distortion = cur_distortion + integral(@(u) (q(i) - u).^2.*distr(u), r(i), r(i+1));
                end
                all_distortions = [all_distortions cur_distortion];
                
                cur_p = arrayfun(@(i) integral(@(x) distr(x), r(i), r(i+1)), 1:8);

                cur_entropie = -sum((cur_p .* log2(cur_p)));
                all_entropie = [all_entropie cur_entropie];
                
                if (l>1 && (prev_distortion - cur_distortion)/prev_distortion < EPS)
                    break
                end
                l = l + 1;
            end
            
            GKD = cur_distortion;
            gem=integral(@(u) u.*distr(u), 0, 255);
            var=integral(@(u) distr(u).*(u-gem).^2, 0, 255);
            SQR = 10*log10(var/GKD);
            p = cur_p;
            entropie =  cur_entropie;
            
            %{
            d_plot = plot(1:l, all_distortions, 'Color', 'Blue');
            leg = legend([d_plot], 'GKD $\sigma$', 'Location', 'northeast');
            set(leg,'Interpreter','latex')
            t = title('GKD, $\sigma_{e}^{2}$, in functie van het aantal iteraties, $l$.', 'FontSize', 16);
            set(t,'Interpreter','latex');
            x = xlabel('$l$', 'FontSize',15);
            set(x,'Interpreter','latex');
            y = ylabel('$\sigma_{e}^{2}$', 'FontSize',15);
            set(y,'Interpreter','latex')
            axis([0 146 35 74]);
            hold off;
            %}
            
            %{
            hold on
            e_plot = plot(1:l, all_entropie, 'Color', 'Blue');
            ref_plot = plot(1:l, ones(1,l)*3, 'Color', 'Red');
            leg = legend([e_plot ref_plot], 'entropie, $H$', 'referentie entropie, $H_{ref}$', 'Location', 'northeast');
            set(leg,'Interpreter','latex')
            t = title('Entropie, $H$, in functie van het aantal iteraties, $l$.', 'FontSize', 16);
            set(t,'Interpreter','latex');
            x = xlabel('$l$', 'FontSize',15);
            set(x,'Interpreter','latex');
            y = ylabel('$H$', 'FontSize',15);
            set(y,'Interpreter','latex')
            axis([0 146 2 3.2]);
            hold off;
            %}
            
             % PLOT
            f_plot = plot(0:255, distr(0:255), 'Color', 'Blue');
            hold on
            r_plot = stem(r, distr(r), 'Color', 'Red');
            q_plot = stem(q, distr(q), 'Color','Green');
            l = legend([r_plot, q_plot, f_plot], 'kwantisatiedrempels', 'reconstructieniveaus', 'f_{u}(u)', 'Location', 'northeast');
            set(l,'Interpreter','latex')
            t = title('Lloyd\textendash Max kwantisator', 'FontSize',16);
            set(t,'Interpreter','latex')
            
            xlabel('u', 'FontSize',15);
            y = ylabel('$f_{u}(u)$', 'FontSize',15);
            set(y,'Interpreter','latex')
            
            axis([0 255 0 0.016]);
            hold off;
        end
        
        % Functie om de compansie kwantisator te bepalen
        function [GKD,SQR,entropie,r,q,p] = companding_quantizer()
            % OUTPUT
             % GKD : GKD Lloyd-Max kwantisator
             % SQR : SQR Lloyd-Max kwantisator
             % entropie : entropie van het gekwantiseerde signaal
             % r : kwantisatiedrempels
             % q : kwantisatieniveaus
             % p : relatieve frequentie kwantisatieniveaus
            % baseer je voor alle berekeningen op de anonieme functie distr.
            % Om de waarden van inverse van een (stijgende) anonieme functie te bepalen kan je gebruik maken van de functie Quantization.inverse()
            [~, distr] = make_probability_functions(Quantization.filename);
            M = 8;
            g = @(u) (integral(@(x) distr(x), 0, u)) - 1/2;
            
            %g_plot = plot(0:255, arrayfun(g, 0:255));
            %hold on
            
            %l = legend([g_plot], 'g(u)', 'Location', 'northeast');
            %set(l,'Interpreter','latex')
            
            %t = title('Verloop van g(u)', 'FontSize',16);
            %set(t,'Interpreter','latex')
            
            %x = xlabel('$u$', 'FontSize',15);
            %y = ylabel('$g(u)$', 'FontSize',15);
            %set(x,'Interpreter','latex')
            %set(y,'Interpreter','latex')
            
            %axis([0 255 -0.65 0.65]);
            %hold off;
            
            r_i_uni = [-0.5:0.125:0.5];
            q_i_uni = [-0.4375:0.125:0.4375];
            
            r = Quantization.inverse(g,r_i_uni);
            q = Quantization.inverse(g,q_i_uni);
            
            GKD = sum(arrayfun(@(i) integral(@(u) (u - q(i)).^2.*distr(u), r(i), r(i+1)), 1:M));
            gem = integral(@(u) u .* distr(u), 0, 255);
            var = integral(@(u) (u - gem).^2.*distr(u),0,255);
            SQR = 10*log10(var/GKD);
            
            p = arrayfun(@(i) integral(@(x) distr(x), r(i), r(i+1)), 1:M);
            entropie = - sum (p .* log2(p));
            
            % PLOT
            f_plot = plot(0:255, distr(0:255), 'Color', 'Green');
            hold on
            r_plot = stem(r, distr(r), 'Color', 'Red');
            q_plot = stem(q, distr(q), 'Color','Blue');
            
            l = legend([r_plot, q_plot, f_plot], 'reconstructieniveaus', 'kwantisatiedrempels', 'f_{u}(u)', 'Location', 'northeast');
            set(l,'Interpreter','latex')
            
            t = title('Compansiekwantisator', 'FontSize',16);
            set(t,'Interpreter','latex')
            
            xlabel('u', 'FontSize',15);
            y = ylabel('$f_{u}(u)$', 'FontSize',15);
            set(y,'Interpreter','latex')
            
            axis([0 255 0 0.016]);
            hold off;
        end
        
        function [samples_quantized]= quantize(r,q)
            % functie om originele figuur te kwantiseren
            % INPUT            
             % r : kwantisatiedremples
             % q : kwantisatieniveaus
            % OUTPUT
             % samples_quantized : sequentie gekwantiseerd signaal
            [samples, ~] = make_probability_functions(Quantization.filename);
            samples = samples(:)'; % maak rij van samples
            samples_quantized = ones(size(samples));
            
            for i = 1 : length(samples)
                index = 1;
                for j = 1 : length(r)
                   if(r(j) < samples(i))
                       index = j;
                       j = length(r) + 1;
                   end
                end
                samples_quantized(i) = q(index);
            end
        end
        
        %--------- Hierna niets veranderen -------------
              
        % functie om originele en gekwantiseerde figuur te plotten 
        function show_figures(samples_quantized)
            % functie om gekwantiseerde foto te vergelijken met de originele
            % INPUT 
             % samples_quantized : gekwantiseerde samples
            [samples, ~] = make_probability_functions(Quantization.filename);
            if(numel(samples)~=numel(samples_quantized)) % controle
                error('Number of elements in samples_quantized not correct')
            end
            figure('Name','Origineel')
            imshow(uint8(samples));
            figure('Name','Gekwantiseerd')
            imshow(uint8(reshape(samples_quantized,size(samples))));
        end
        
        % functie om inverse te bepalen van stijgende anonieme functie g in de y-waarden in Yvec
        % g is een anonieme functie
        % Yvec bevat de y waarden waarvoor je de inverse wil berekenen
        function Uvec=inverse(g,Yvec)
            Uvec = zeros(size(Yvec));
            eps = 1e-6;
            for i = 1 : length(Yvec)
                y = Yvec(i);
                a = 0; g_a = g(a);
                b = 255; g_b = g(b);
                nit = 1;
                while(abs(g_a-y)>eps && abs(g_b-y)>eps && nit<1000)
                    u_test = (a+b)/2;
                    g_test = g(u_test);
                    if(g_test<=y)
                        a = u_test;
                        g_a = g_test;
                    else
                        b = u_test;
                        g_b = g_test;
                    end
                    nit = nit+1;
                end
                if(abs(g_a-y)<=eps)
                    Uvec(i) = a;
                else
                    Uvec(i) = b;
                end
            end
        end
          
    end
end
    