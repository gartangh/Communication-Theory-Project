function p_e = Calc_crit( p )
    H = [[1 1 0 1 1 0 0 1 0 0 1 0 0 1];
         [0 0 1 1 1 1 0 1 0 0 1 1 1 0];
         [0 0 1 1 0 0 1 1 0 1 0 1 0 1];
         [1 0 0 0 0 0 0 1 1 1 1 1 1 0]];
    H_T = H.';
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
    n = 14;
    G = [[1 1 0 0 0 0 0 0 1 0 0 0 0 0];
         [0 1 0 0 0 1 0 0 1 0 1 0 0 0];
         [0 0 1 0 0 1 1 0 0 0 0 0 0 0];
         [0 1 0 1 0 1 1 0 0 0 0 0 0 0];
         [0 1 0 0 1 1 0 0 0 0 0 0 0 0];
         [0 0 0 0 0 1 1 0 1 0 0 1 0 0];
         [0 0 0 0 0 1 0 0 1 0 0 0 1 0];
         [0 1 0 0 0 1 1 1 1 0 0 0 0 0];
         [0 1 0 0 0 0 1 0 0 0 0 0 0 1];
         [0 0 0 0 0 0 1 0 1 1 0 0 0 0];];

    codewords = zeros(0,14);
    for i=1:1023
        % informatiewoord b (1x10) opstellen
        b = zeros(1,10);
        teller = dec2bin(i);
        for j=(10-length(teller)+1):10
            b(1, j) = str2num(teller(j-10+length(teller)));
        end
        % codewoord c (1x14) berekenen
        c = Channel_Coding.Encode_outer(b);

        % checken of codewoord al in tabel zit
        [check, ~] = ismember(codewords, c, 'rows');
        if(sum(check) == 0)
            codewords = [codewords;c];
        end
    end

    % kans op een fout
    p_e = 0;
    for j=1:size(codewords,1)
        for k=1:16
            S(k,:);
            c_plus_e = mod(codewords(j, :) + S(k,:),2);
            w = sum(c_plus_e);
            p_e = p_e + p^(w)*(1-p)^(14-w);
        end
    end

    p_e = p_e - 0.001; % voor zoeken bovengrens p bij i.f.v. ontwerpcriterium
end

