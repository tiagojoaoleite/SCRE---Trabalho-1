# SCRI 2026/27 \- TP1: Estimação de Estado em Redes de Energia (WLS)

# ==================== MODELO ====================
```matlab
% --- Admitancias serie das linhas (pu) ---
y12a = 2 - 10j;          % linha a entre 1 e 2
y12b = 2 - 10j;          % linha b entre 1 e 2 (paralela)
y13  = 2 - 10j;
y23  = 2 - 10j;


% --- Admitancia shunt do banco de condensadores (barramento 1) ---
Bcap = 2.2;              % shunt = j*Bcap  ->  Qgbat = Bcap*V1^2


% --- Medicoes (pu) ---
V1    = 1.10;
V2    = 1.00;
V3    = 1.07;
P21a  = -0.79;
P21b  =  0.79;
P31   =  1.71;
P32   =  2.34;
Q12a  =  0.97;
Q32   =  1.22;
Qgbat =  2.66;
P1v   =  0;              % medida virtual: barramento 1 sem carga/geracao
Q1v   =  0;              % medida virtual


% --- Variancias dos erros de medida ---
varV1    = 0.0002;
varV2    = 0.0002;
varV3    = 0.0002;
varP21a  = 0.001;
varP21b  = 0.001;
varP31   = 0.001;
varP32   = 0.001;
varQ12a  = 0.005;
varQ32   = 0.005;
varQgbat = 0.005;
varP1v   = 0.00001;
varQ1v   = 0.00001;


% --- Ângulos conhecidos ---
t2 = 0;                  % barramento 2 = slack (referencia). t1 e t3 sao estimados.


% --- Parametros do metodo ---
tol   = 1e-4;            % tolerancia em max|dx|
maxit = 50;
alfa  = 0.05;            % nivel de significancia do teste chi^2
limiar_rN = 3;           % limiar para |rN|


% Estrutura com os parametros da rede
rede.g12a = real(y12a);  rede.b12a = imag(y12a);
rede.g12b = real(y12b);  rede.b12b = imag(y12b);
rede.g13  = real(y13);   rede.b13  = imag(y13);
rede.g23  = real(y23);   rede.b23  = imag(y23);
rede.Bcap = Bcap;        rede.t2   = t2;
```

# =================== MATRIZ Ybus ===================
```matlab
Y11 =  y12a + y12b + y13 + 1j*Bcap;
Y22 =  y12a + y12b + y23;
Y33 =  y13 + y23;
Y12 = -(y12a + y12b);    Y21 = Y12;
Y13 = -y13;              Y31 = Y13;
Y23 = -y23;              Y32 = Y23;
Ybus = [Y11 Y12 Y13;
        Y21 Y22 Y23;
        Y31 Y32 Y33];


BusNames = {'ONE', 'TWO', 'THREE'};
disp_Ybus = array2table(Ybus, 'RowNames', BusNames, 'VariableNames', BusNames)
```

| |ONE|TWO|THREE|
|:--:|:--:|:--:|:--:|
|1 ONE|6.0000 -27.8000i|-4.0000 +20.0000i|-2.0000 +10.0000i|
|2 TWO|-4.0000 +20.0000i|6.0000 -30.0000i|-2.0000 +10.0000i|
|3 THREE|-2.0000 +10.0000i|-2.0000 +10.0000i|4.0000 -20.0000i|

# ================= VETOR DE ESTADO =================
```matlab
% x = [t1, t3, V1, V2, V3]'   (angulos em rad)
nomes_x = {'t1', 't3', 'V1', 'V2', 'V3'};
n = numel(nomes_x);
```

# ================ VETOR DE MEDIDAS ================
```matlab
nomes_z = {'V1','V2','V3','P21a','P21b','P31','P32','Q12a','Q32','Qgbat','P1v','Q1v'};
z   = [V1; V2; V3; P21a; P21b; P31; P32; Q12a; Q32; Qgbat; P1v; Q1v];
var = [varV1; varV2; varV3; varP21a; varP21b; varP31; varP32; ...
       varQ12a; varQ32; varQgbat; varP1v; varQ1v];
virtual = ismember(nomes_z, {'P1v','Q1v'})';
disp(table(z, var, 'RowNames', nomes_z, 'VariableNames', {'z','variancia'}));
```

```matlabTextOutput
               z      variancia
             _____    _________


    V1         1.1     0.0002  
    V2           1     0.0002  
    V3        1.07     0.0002  
    P21a     -0.79      0.001  
    P21b      0.79      0.001  
    P31       1.71      0.001  
    P32       2.34      0.001  
    Q12a      0.97      0.005  
    Q32       1.22      0.005  
    Qgbat     2.66      0.005  
    P1v          0      1e-05  
    Q1v          0      1e-05  
```

```matlab
% Verificacao do Jacobiano analitico (diferencas finitas)
x_teste = [0.05; 0.2; 1.08; 0.97; 1.10];
Hnum = zeros(numel(z), n);  e = 1e-7;
for c = 1:n
    xp = x_teste; xp(c) = xp(c) + e;
    Hnum(:,c) = (h_fun(xp, rede) - h_fun(x_teste, rede)) / e;
end
fprintf('Verificacao do Jacobiano: max|H - Hnum| = %.2e (deve ser ~1e-6)\n\n', ...
        max(max(abs(H_fun(x_teste, rede) - Hnum))));
```

```matlabTextOutput
Verificacao do Jacobiano: max|H - Hnum| = 2.82e-06 (deve ser ~1e-6)
```

# PERGUNTA 1 \- WLS com todas as medidas + 2 medidas virtuais
```matlab
ativas_P1 = true(numel(z), 1);
P1 = gauss_newton(z, var, ativas_P1, nomes_z, nomes_x, rede, tol, maxit, alfa);


disp('Tabela de estado por iteração (angulos em graus):');      disp(P1.T_estado);
```

```matlabTextOutput
Tabela de estado por iteração (angulos em graus):
          It0      It1        It2       It3        It4  
          ___    _______    _______    ______    _______


    t1     0      2.4736     2.5004    2.5041     2.5042
    t3     0      10.334      9.374    9.3781     9.3786
    V1     1      1.0944      1.089     1.089      1.089
    V2     1     0.96359    0.97326    0.9735    0.97351
    V3     1      1.1048     1.0997    1.0995     1.0995
```

```matlab
disp('Tabela de residuos por iteracao (r = z - h(x)):');        disp(P1.T_residuos);
```

```matlabTextOutput
Tabela de residuos por iteracao (r = z - h(x)):
              It0        It1          It2          It3           It4    
             _____    _________    _________    __________    __________


    V1         0.1    0.0055595     0.010973       0.01099      0.010991
    V2           0     0.036408      0.02674        0.0265      0.026493
    V3        0.07     -0.03477    -0.029719      -0.02953     -0.029513
    P21a     -0.79    -0.084638     -0.10428      -0.10395      -0.10393
    P21b      0.79       1.4954       1.4757        1.4761        1.4761
    P31       1.71     0.010912      0.23595       0.23652        0.2365
    P32       2.34     0.083895      0.28997       0.29005       0.29002
    Q12a      0.97     -0.38086     -0.20835       -0.2054       -0.2053
    Q32       1.22     -0.13045     0.035045      0.040036      0.040323
    Qgbat     0.46      0.02484     0.050842      0.050927      0.050932
    P1v          0      0.16636    0.0062759     0.0052511     0.0052518
    Q1v        2.2     -0.39782    -0.003845    0.00013726    0.00013753
```

```matlab
disp('Tabela J(x) por iteracao e teste chi^2:');                disp(P1.T_J);
```

```matlabTextOutput
Tabela J(x) por iteracao e teste chi^2:
               J         k    chi2_crit    max_abs_dx         decisao     
           __________    _    _________    __________    _________________


    It0    4.9425e+05    7     14.067             NaN    {'dados errados'}
    It1         20889    7     14.067         0.18036    {'dados errados'}
    It2        2351.8    7     14.067        0.016753    {'dados errados'}
    It3        2350.1    7     14.067      0.00023929    {'dados errados'}
    It4        2350.1    7     14.067      1.6623e-05    {'dados errados'}
```

```matlab
disp('Tabela covariancia dos residuos (Omega_ii) e rN por iteracao:'); disp(P1.T_rN);
```

```matlabTextOutput
Tabela covariancia dos residuos (Omega_ii) e rN por iteracao:
             Omega_It0     rN_It0     Omega_It1     rN_It1     Omega_It2     rN_It2     Omega_It3     rN_It3     Omega_It4     rN_It4 
             __________    _______    __________    _______    __________    _______    __________    _______    __________    _______


    V1       0.00013923     8.4748     0.0001459    0.46026    0.00014531    0.91025     0.0001453    0.91174     0.0001453    0.91182
    V2       0.00014776          0    0.00014434     3.0305    0.00014478     2.2223    0.00014479     2.2023     0.0001448     2.2017
    V3       0.00014465     5.8203     0.0001508    -2.8314    0.00015047    -2.4227    0.00015046    -2.4075    0.00015045    -2.4061
    P21a     0.00093226    -25.874    0.00093299    -2.7709    0.00093324    -3.4136    0.00093325    -3.4026    0.00093325    -3.4022
    P21b     0.00093226     25.874    0.00093299     48.956    0.00093324     48.307    0.00093325     48.317    0.00093325     48.318
    P31      0.00073237     63.188    0.00069963    0.41255    0.00070371     8.8945    0.00070374     8.9158    0.00070374     8.9151
    P32      0.00040032     116.95    0.00043054     4.0432    0.00042577     14.053    0.00042571     14.058    0.00042571     14.057
    Q12a      0.0044703     14.508     0.0044348     -5.719     0.0044419    -3.1261     0.0044417     -3.082     0.0044417    -3.0805
    Q32      0.00091109     40.418    0.00079481     -4.627    0.00080187     1.2376    0.00080236     1.4134    0.00080239     1.4235
    Qgbat     0.0038236     7.4392     0.0037455    0.40588     0.0037443    0.83089     0.0037442    0.83228     0.0037442    0.83236
    P1v      3.3245e-08          0    3.3317e-08     911.43    3.3211e-08     34.438     3.321e-08     28.815     3.321e-08     28.819
    Q1v      2.9027e-09      40834    2.9612e-09    -7310.6    2.9401e-09    -70.912    2.9396e-09     2.5316    2.9396e-09     2.5366
```

```matlab
disp(table(z, P1.h_final, z - P1.h_final, 'RowNames', nomes_z, 'VariableNames', {'z','z_estimado','residuo'}));
```

```matlabTextOutput
               z      z_estimado      residuo  
             _____    ___________    __________


    V1         1.1          1.089      0.010991
    V2           1        0.97351      0.026493
    V3        1.07         1.0995     -0.029513
    P21a     -0.79       -0.68607      -0.10393
    P21b      0.79       -0.68607        1.4761
    P31       1.71         1.4735        0.2365
    P32       2.34           2.05       0.29002
    Q12a      0.97         1.1753       -0.2053
    Q32       1.22         1.1797      0.040323
    Qgbat     2.66         2.6091      0.050932
    P1v          0     -0.0052518     0.0052518
    Q1v          0    -0.00013753    0.00013753
```

```matlab
% b) teste chi^2
mostra_teste_chi2(P1, alfa);
```

```matlabTextOutput
Convergiu em 4 iteracoes.
Estado: t1 = 2.5042 graus, t3 = 9.3786 graus, V1 = 1.0890, V2 = 0.9735, V3 = 1.0995 pu
J(x) = 2350.12 | k = 7 | chi2(7; 95%) = 14.07 -> existem dados errados (95% de confianca)
```

# PERGUNTA 2 \- Chi^2 + residuos normalizados \-> medida com erro grosseiro
```matlab
% Ordenar as medidas ativas por |rN| (resultado final da Pergunta 1)
idxA = find(ativas_P1);
[~, ord] = sort(abs(P1.rN_final), 'descend');
T_ord = table(abs(P1.rN_final(ord)), virtual(idxA(ord)), ...
              'RowNames', nomes_z(idxA(ord)), 'VariableNames', {'abs_rN','virtual'});
disp('Medidas ordenadas por |rN|:'); disp(T_ord);
```

```matlabTextOutput
Medidas ordenadas por |rN|:
             abs_rN     virtual
             _______    _______


    P21b      48.318     false 
    P1v       28.819     true  
    P32       14.057     false 
    P31       8.9151     false 
    P21a      3.4022     false 
    Q12a      3.0805     false 
    Q1v       2.5366     true  
    V3        2.4061     false 
    V2        2.2017     false 
    Q32       1.4235     false 
    V1       0.91182     false 
    Qgbat    0.83236     false 
```

```matlab
% Candidata a remocao: maior |rN| entre as medidas REAIS (virtuais sao exatas)
rN_reais = abs(P1.rN_final);  rN_reais(virtual(idxA)) = -Inf;
[rmax, q] = max(rN_reais);
cand_P2 = idxA(q);
fprintf('Teste chi^2: J = %.2f >= chi2 = %.2f -> existem dados errados.\n', P1.J_final, P1.chi2_crit);
```

```matlabTextOutput
Teste chi^2: J = 2350.12 >= chi2 = 14.07 -> existem dados errados.
```

```matlab
fprintf('Candidata a remocao: %s (|rN| = %.2f > %.1f)\n', nomes_z{cand_P2}, rmax, limiar_rN);
```

```matlabTextOutput
Candidata a remocao: P21b (|rN| = 48.32 > 3.0)
```

```matlab
% Nova estimacao sem a candidata
ativas_P2 = ativas_P1;  ativas_P2(cand_P2) = false;
P2 = gauss_newton(z, var, ativas_P2, nomes_z, nomes_x, rede, tol, maxit, alfa);


disp('Tabela de estado por iteracao (angulos em graus):');      disp(P2.T_estado);
```

```matlabTextOutput
Tabela de estado por iteracao (angulos em graus):
          It0      It1        It2        It3        It4  
          ___    _______    _______    _______    _______


    t1     0      3.0704     3.0782     3.0812     3.0812
    t3     0      12.093     10.998     10.998     10.998
    V1     1       1.094     1.0911     1.0912     1.0912
    V2     1     0.96136    0.97817    0.97853    0.97853
    V3     1      1.1081      1.099     1.0989     1.0989
```

```matlab
disp('Tabela de residuos por iteracao:');                       disp(P2.T_residuos);
```

```matlabTextOutput
Tabela de residuos por iteracao:
              It0        It1          It2           It3           It4    
             _____    _________    __________    __________    __________


    V1         0.1    0.0059799     0.0089395     0.0088317     0.0088326
    V2           0     0.038645      0.021829      0.021473      0.021473
    V3        0.07    -0.038057     -0.029023      -0.02894     -0.028937
    P21a     -0.79     0.025399    0.00086883     0.0012731     0.0012727
    P31       1.71     -0.25219      0.017513      0.018442      0.018443
    P32       2.34     -0.26405     -0.015935     -0.015628     -0.015626
    Q12a      0.97     -0.38381      -0.16247      -0.15976      -0.15975
    Q32       1.22     -0.19559       0.10455       0.10953       0.10956
    Qgbat     0.46     0.026864      0.041091      0.040574      0.040578
    P1v          0      0.18853     0.0015845    7.5338e-05    7.5835e-05
    Q1v        2.2     -0.45113    -0.0039266    8.7562e-05    8.7137e-05
```

```matlab
disp('Tabela J(x) por iteracao e teste chi^2:');                disp(P2.T_J);
```

```matlabTextOutput
Tabela J(x) por iteracao e teste chi^2:
               J         k    chi2_crit    max_abs_dx         decisao     
           __________    _    _________    __________    _________________


    It0    4.9363e+05    6     12.592             NaN    {'dados errados'}
    It1         24092    6     12.592         0.21106    {'dados errados'}
    It2        17.151    6     12.592        0.019119    {'dados errados'}
    It3        15.304    6     12.592      0.00035504    {'dados errados'}
    It4        15.304    6     12.592       2.597e-06    {'dados errados'}
```

```matlab
disp('Tabela covariancia dos residuos (Omega_ii) e rN por iteracao:'); disp(P2.T_rN);
```

```matlabTextOutput
Tabela covariancia dos residuos (Omega_ii) e rN por iteracao:
             Omega_It0     rN_It0     Omega_It1     rN_It1     Omega_It2      rN_It2     Omega_It3      rN_It3     Omega_It4      rN_It4 
             __________    _______    __________    _______    __________    ________    __________    ________    __________    ________


    V1       0.00013923     8.4748     0.0001459    0.49507    0.00014534     0.74152    0.00014534     0.73258    0.00014534     0.73265
    V2       0.00014775          0    0.00014471     3.2125    0.00014538      1.8104    0.00014539      1.7809    0.00014539      1.7808
    V3       0.00014464     5.8203    0.00015078    -3.0993    0.00015003     -2.3695    0.00015001     -2.3628    0.00015001     -2.3626
    P21a     0.00092734    -25.942    0.00092944     0.8331    0.00092983    0.028493    0.00092985    0.041749    0.00092985    0.041736
    P31      0.00071363     64.012    0.00067814    -9.6845    0.00068308     0.67008    0.00068312     0.70559    0.00068312     0.70563
    P32      0.00035746     123.77    0.00039069    -13.359    0.00038514    -0.81199    0.00038508    -0.79637    0.00038508    -0.79627
    Q12a      0.0044703     14.508     0.0044204    -5.7727     0.0044242     -2.4426     0.0044239      -2.402     0.0044239     -2.4018
    Q32      0.00091107     40.419    0.00079465    -6.9385    0.00081432      3.6636    0.00081478      3.8371    0.00081479      3.8383
    Qgbat     0.0038235     7.4392     0.0037465     0.4389     0.0037403     0.67189       0.00374     0.66345       0.00374     0.66352
    P1v      2.1404e-08          0    2.1839e-08     1275.7    2.1679e-08      10.762    2.1676e-08     0.51171    2.1676e-08     0.51509
    Q1v      2.9027e-09      40834     2.991e-09    -8248.8    2.9476e-09     -72.324     2.947e-09       1.613     2.947e-09      1.6051
```

```matlab
mostra_teste_chi2(P2, alfa);
```

```matlabTextOutput
Convergiu em 4 iteracoes.
Estado: t1 = 3.0812 graus, t3 = 10.9980 graus, V1 = 1.0912, V2 = 0.9785, V3 = 1.0989 pu
J(x) = 15.30 | k = 6 | chi2(6; 95%) = 12.59 -> existem dados errados (95% de confianca)
```

# PERGUNTA 3 \- Remocao sequencial ate obter o melhor subconjunto
```matlab
ativas_P3 = ativas_P2;
P3 = P2;
removidas = {nomes_z{cand_P2}};
while true
    idxA = find(ativas_P3);
    rN_reais = abs(P3.rN_final);  rN_reais(virtual(idxA)) = -Inf;
    [rmax, q] = max(rN_reais);
    if P3.J_final < P3.chi2_crit && rmax < limiar_rN
        fprintf('J = %.2f < chi2 = %.2f e todos |rN| < %.1f -> conjunto final.\n', ...
                P3.J_final, P3.chi2_crit, limiar_rN);
        break
    end
    j = idxA(q);
    fprintf('Remover %s (|rN| = %.2f; J = %.2f vs chi2 = %.2f)\n', ...
            nomes_z{j}, rmax, P3.J_final, P3.chi2_crit);
    if nnz(ativas_P3) - 1 <= n
        warning('Remover %s tornaria o sistema nao observavel (m <= n).', nomes_z{j}); break
    end
    ativas_P3(j) = false;  removidas{end+1} = nomes_z{j}; %#ok<SAGROW>
    P3 = gauss_newton(z, var, ativas_P3, nomes_z, nomes_x, rede, tol, maxit, alfa);
    if rcond(P3.G_final) < 1e-12
        warning('Matriz de ganho singular: estado nao observavel.');
    end
end
```

```matlabTextOutput
Remover Q32 (|rN| = 3.84; J = 15.30 vs chi2 = 12.59)
J = 0.56 < chi2 = 11.07 e todos |rN| < 3.0 -> conjunto final.
```

```matlab
fprintf('Medidas removidas: %s\n', strjoin(removidas, ', '));
```

```matlabTextOutput
Medidas removidas: P21b, Q32
```

```matlab
disp('Tabela de estado por iteracao (angulos em graus):');      disp(P3.T_estado);
```

```matlabTextOutput
Tabela de estado por iteracao (angulos em graus):
          It0      It1       It2       It3       It4  
          ___    _______    ______    ______    ______


    t1     0      3.3455    3.1535    3.1545    3.1545
    t3     0      12.972    11.543    11.547    11.547
    V1     1      1.1051    1.0981     1.098     1.098
    V2     1     0.99605    1.0013    1.0014    1.0014
    V3     1      1.0672    1.0733    1.0735    1.0735
```

```matlab
disp('Tabela de residuos por iteracao:');                       disp(P3.T_residuos);
```

```matlabTextOutput
Tabela de residuos por iteracao:
              It0        It1           It2            It3           It4    
             _____    __________    __________    ___________    __________


    V1         0.1    -0.0051101     0.0019382      0.0019802     0.0019798
    V2           0     0.0039497    -0.0013111      -0.001401     -0.001401
    V3        0.07     0.0027694    -0.0033041     -0.0035122    -0.0035133
    P21a     -0.79      0.065859     0.0052691      0.0052527     0.0052522
    P31       1.71      -0.21467      0.018416       0.017065      0.017066
    P32       2.34      -0.25242     -0.008518       -0.01013     -0.010129
    Q12a      0.97      -0.12552      0.011937       0.013461      0.013456
    Qgbat     0.46     -0.026791     0.0073724      0.0075756     0.0075737
    P1v          0      0.081065    -0.0013873     7.4114e-05    7.4115e-05
    Q1v        2.2       -0.4834    -0.0053685    -2.4834e-05    -2.458e-05
```

```matlab
disp('Tabela J(x) por iteracao e teste chi^2:');                disp(P3.T_J);
```

```matlabTextOutput
Tabela J(x) por iteracao e teste chi^2:
               J         k    chi2_crit    max_abs_dx         decisao      
           __________    _    _________    __________    __________________


    It0    4.9333e+05    5      11.07             NaN    {'dados errados' }
    It1         24143    5      11.07         0.22641    {'dados errados' }
    It2        3.6353    5      11.07        0.024946    {'nao se rejeita'}
    It3       0.56086    5      11.07      0.00020813    {'nao se rejeita'}
    It4       0.56086    5      11.07      1.0976e-06    {'nao se rejeita'}
```

```matlab
disp('Tabela covariancia dos residuos (Omega_ii) e rN por iteracao:'); disp(P3.T_rN);
```

```matlabTextOutput
Tabela covariancia dos residuos (Omega_ii) e rN por iteracao:
             Omega_It0     rN_It0     Omega_It1      rN_It1     Omega_It2      rN_It2     Omega_It3      rN_It3     Omega_It4      rN_It4 
             __________    _______    __________    ________    __________    ________    __________    ________    __________    ________


    V1       0.00013558     8.5882      0.000143    -0.42733     0.0001423     0.16247    0.00014229     0.16601    0.00014229     0.16597
    V2       0.00011199          0    0.00011055     0.37565    0.00011099    -0.12445    0.00011099    -0.13298    0.00011099    -0.13298
    V3        9.512e-05     7.1773    0.00010612     0.26884    0.00010542     -0.3218    0.00010543    -0.34206    0.00010543    -0.34216
    P21a     0.00092732    -25.943    0.00093036      2.1592    0.00093019     0.17276    0.00093021     0.17222    0.00093021     0.17221
    P31      0.00071356     64.015    0.00067756      -8.247    0.00068284     0.70476     0.0006829     0.65303     0.0006829     0.65307
    P32      0.00035729      123.8    0.00038884     -12.801    0.00038316    -0.43516    0.00038309    -0.51757    0.00038309    -0.51751
    Q12a      0.0026774     18.746     0.0023604     -2.5835     0.0023863     0.24435     0.0023863     0.27556     0.0023863     0.27545
    Qgbat     0.0037528      7.509     0.0036524     -0.4433     0.0036531     0.12198      0.003653     0.12534      0.003653     0.12531
    P1v      2.1404e-08          0      2.14e-08      554.16    2.1426e-08     -9.4777    2.1427e-08     0.50631    2.1427e-08     0.50632
    Q1v      2.1986e-09      46919    2.0991e-09      -10551    2.0848e-09     -117.58    2.0849e-09     -0.5439    2.0849e-09    -0.53832
```

```matlab
mostra_teste_chi2(P3, alfa);
```

```matlabTextOutput
Convergiu em 4 iteracoes.
Estado: t1 = 3.1545 graus, t3 = 11.5465 graus, V1 = 1.0980, V2 = 1.0014, V3 = 1.0735 pu
J(x) = 0.56 | k = 5 | chi2(5; 95%) = 11.07 -> nao se pode afirmar que existem dados errados
```

```matlab
% Valores estimados de TODAS as grandezas medidas (incluindo removidas)
estado = repmat({'ativa'}, numel(z), 1);  estado(~ativas_P3) = {'removida'};
disp('P3 - Comparacao com as medidas originais:');
```

```matlabTextOutput
P3 - Comparacao com as medidas originais:
```

```matlab
disp(table(z, P3.h_final, z - P3.h_final, estado, 'RowNames', nomes_z, ...
     'VariableNames', {'z','z_estimado','diferenca','estado'}));
```

```matlabTextOutput
               z      z_estimado     diferenca        estado   
             _____    ___________    __________    ____________


    V1         1.1          1.098     0.0019798    {'ativa'   }
    V2           1         1.0014     -0.001401    {'ativa'   }
    V3        1.07         1.0735    -0.0035133    {'ativa'   }
    P21a     -0.79       -0.79525     0.0052522    {'ativa'   }
    P21b      0.79       -0.79525        1.5853    {'removida'}
    P31       1.71         1.6929      0.017066    {'ativa'   }
    P32       2.34         2.3501     -0.010129    {'ativa'   }
    Q12a      0.97        0.95654      0.013456    {'ativa'   }
    Q32       1.22        0.56133       0.65867    {'removida'}
    Qgbat     2.66         2.6524     0.0075737    {'ativa'   }
    P1v          0    -7.4115e-05    7.4115e-05    {'ativa'   }
    Q1v          0      2.458e-05    -2.458e-05    {'ativa'   }
```

```matlab
% Producao dos geradores (carga no barramento 2 = 6 + j2 pu)
xf = P3.x_final;
Vc = [xf(3)*exp(1j*xf(1));             % V1 /_ t1
      xf(4)*exp(1j*t2);                % V2 /_ t2 = 0
      xf(5)*exp(1j*xf(2))];            % V3 /_ t3
Sinj = Vc .* conj(Ybus * Vc);           % injecoes liquidas
Scarga2 = 6 + 2j;
SG2 = Sinj(2) + Scarga2;
SG3 = Sinj(3);
disp('Producao estimada dos geradores:');
```

```matlabTextOutput
Producao estimada dos geradores:
```

```matlab
disp(table([real(SG2); real(SG3)], [imag(SG2); imag(SG3)], ...
     'RowNames', {'G2 (slack)','G3'}, 'VariableNames', {'P_pu','Q_pu'}));
```

```matlabTextOutput
                   P_pu       Q_pu  
                  ______    ________


    G2 (slack)    2.1568     0.26604
    G3            4.0431    0.080398
```

```matlab
fprintf('Perdas ativas = %.4f pu | S1 (verificacao, ~0) = %.4f %+.4fj\n', ...
        real(sum(Sinj)), real(Sinj(1)), imag(Sinj(1)));
```

```matlabTextOutput
Perdas ativas = 0.1998 pu | S1 (verificacao, ~0) = -0.0001 +0.0000j
```

```matlab
% Resumo para o relatorio
disp('Resumo:');
```

```matlabTextOutput
Resumo:
```

```matlab
disp(table([P1.J_final; P2.J_final; P3.J_final], [P1.k; P2.k; P3.k], ...
     [P1.chi2_crit; P2.chi2_crit; P3.chi2_crit], ...
     [max(abs(P1.rN_final)); max(abs(P2.rN_final)); max(abs(P3.rN_final))], ...
     'RowNames', {'P1_todas','P2','P3_final'}, ...
     'VariableNames', {'J','k','chi2_crit','max_abs_rN'}));
```

```matlabTextOutput
                   J       k    chi2_crit    max_abs_rN
                _______    _    _________    __________


    P1_todas     2350.1    7     14.067        48.318  
    P2           15.304    6     12.592        3.8383  
    P3_final    0.56086    5      11.07       0.65307  
```

```matlab
% Grafico de |rN| nas tres estimacoes
figure('Name','Residuos normalizados');
dados = nan(numel(z), 3);
dados(ativas_P1,1) = abs(P1.rN_final);
dados(ativas_P2,2) = abs(P2.rN_final);
dados(ativas_P3,3) = abs(P3.rN_final);
bar(categorical(nomes_z, nomes_z), dados); hold on;
yline(limiar_rN, 'r--');                       % requer R2018b+
legend('Pergunta 1','Pergunta 2','Pergunta 3','|r^N| = 3'); ylabel('|r_i^N|'); grid on;
set(gca, 'YScale', 'log');
```

![figure_0.png](./scri_trabalho1_media/figure_0.png)

# ============================ FUNCOES ============================
```matlab
function R = gauss_newton(z, var, ativas, nomes_z, nomes_x, rede, tol, maxit, alfa)
% Iteracao de Gauss-Newton:  x(k+1) = x(k) + G^-1 * H' * W * (z - h(x(k)))
    za = z(ativas);    nz = nomes_z(ativas);
    W  = diag(1 ./ var(ativas));
    Rm = diag(var(ativas));
    n  = numel(nomes_x);  m = nnz(ativas);  k = m - n;
    chi2c = chi2crit(k, alfa);


    % Flat start
    flat_zero = [0; 0; 1; 1; 1];          % [t1 t3 V1 V2 V3]
    x = flat_zero;


    X_hist = [];  r_hist = [];  J_hist = [];  dx_hist = [];
    Om_hist = [];  rN_hist = [];


    for it = 0:maxit
        h = h_fun(x, rede);    H = H_fun(x, rede);
        h = h(ativas);         H = H(ativas, :);
        r = za - h;                                   % residuos
        G = H' * W * H;                               % matriz de ganho
        J = sum(r.^2 ./ var(ativas));                 % J(x) = sum r^2/sigma^2
        Omega = Rm - H * (G \ H');                    % covariancia dos residuos
        rN = r ./ sqrt(abs(diag(Omega)));             % residuos normalizados


        X_hist  = [X_hist,  x];                       %#ok<AGROW>
        r_hist  = [r_hist,  r];                       %#ok<AGROW>
        J_hist  = [J_hist,  J];                       %#ok<AGROW>
        Om_hist = [Om_hist, diag(Omega)];             %#ok<AGROW>
        rN_hist = [rN_hist, rN];                      %#ok<AGROW>


        if it > 0 && max(abs(dx)) < tol, break; end   % convergiu
        dx = G \ (H' * W * r);
        x  = x + dx;
        dx_hist = [dx_hist, max(abs(dx))];            %#ok<AGROW>
    end
    nit = size(X_hist, 2);
    col = arrayfun(@(i) sprintf('It%d', i), 0:nit-1, 'UniformOutput', false);


    % Tabela 1: estado por iteracao (angulos em graus)
    Xd = X_hist;  Xd(1:2,:) = rad2deg(Xd(1:2,:));
    R.T_estado = array2table(Xd, 'RowNames', nomes_x, 'VariableNames', col);


    % Tabela 2: residuos por iteracao
    R.T_residuos = array2table(r_hist, 'RowNames', nz, 'VariableNames', col);


    % Tabela 3: J(x) por iteracao, comparacao com chi^2(k, alfa)
    dec = repmat({'nao se rejeita'}, nit, 1);  dec(J_hist >= chi2c) = {'dados errados'};
    R.T_J = table(J_hist', repmat(k, nit, 1), repmat(chi2c, nit, 1), ...
                  [NaN, dx_hist(1:nit-1)]', dec, 'RowNames', col, ...
                  'VariableNames', {'J','k','chi2_crit','max_abs_dx','decisao'});


    % Tabela 4: covariancia dos residuos (Omega_ii) e rN por iteracao
    T = table();
    for i = 1:nit
        T.(sprintf('Omega_It%d', i-1)) = Om_hist(:, i);
        T.(sprintf('rN_It%d',    i-1)) = rN_hist(:, i);
    end
    T.Properties.RowNames = nz;
    R.T_rN = T;


    % Resultados finais
    R.x_final = x;               R.h_final = h_fun(x, rede);
    R.r_final = r_hist(:, end);  R.rN_final = rN_hist(:, end);
    R.J_final = J_hist(end);     R.k = k;   R.chi2_crit = chi2c;
    R.G_final = G;               R.nit = nit - 1;
    R.Omega_final = array2table(Omega, 'RowNames', nz, 'VariableNames', nz);
end


function h = h_fun(x, rede)
% Funcoes de medida - cada grandeza calculada individualmente
    t1 = x(1);  t3 = x(2);  V1 = x(3);  V2 = x(4);  V3 = x(5);  t2 = rede.t2;
    g12a = rede.g12a; b12a = rede.b12a;  g12b = rede.g12b; b12b = rede.b12b;
    g13  = rede.g13;  b13  = rede.b13;   g23  = rede.g23;  b23  = rede.b23;
    Bcap = rede.Bcap;


    t12 = t1 - t2;  c12 = cos(t12);  s12 = sin(t12);
    c21 = c12;      s21 = -s12;                 % t21 = -t12
    t31 = t3 - t1;  c31 = cos(t31);  s31 = sin(t31);
    c13 = c31;      s13 = -s31;                 % t13 = -t31
    t32 = t3 - t2;  c32 = cos(t32);  s32 = sin(t32);


    h = zeros(12, 1);
    h(1)  = V1;                                                   % V1
    h(2)  = V2;                                                   % V2
    h(3)  = V3;                                                   % V3
    h(4)  = V2^2*g12a - V2*V1*(g12a*c21 + b12a*s21);              % P21a
    h(5)  = V2^2*g12b - V2*V1*(g12b*c21 + b12b*s21);              % P21b
    h(6)  = V3^2*g13  - V3*V1*(g13*c31  + b13*s31);               % P31
    h(7)  = V3^2*g23  - V3*V2*(g23*c32  + b23*s32);               % P32
    h(8)  = -V1^2*b12a - V1*V2*(g12a*s12 - b12a*c12);             % Q12a
    h(9)  = -V3^2*b23  - V3*V2*(g23*s32  - b23*c32);              % Q32
    h(10) = Bcap*V1^2;                                            % Qgbat


    % Fluxos que saem do barramento 1 (para as medidas virtuais)
    P12a = V1^2*g12a - V1*V2*(g12a*c12 + b12a*s12);
    P12b = V1^2*g12b - V1*V2*(g12b*c12 + b12b*s12);
    P13  = V1^2*g13  - V1*V3*(g13*c13  + b13*s13);
    Q12a = -V1^2*b12a - V1*V2*(g12a*s12 - b12a*c12);
    Q12b = -V1^2*b12b - V1*V2*(g12b*s12 - b12b*c12);
    Q13  = -V1^2*b13  - V1*V3*(g13*s13  - b13*c13);


    h(11) = P12a + P12b + P13;                                    % P1 (virtual)
    h(12) = Q12a + Q12b + Q13 - Bcap*V1^2;                        % Q1 (virtual)
end


function H = H_fun(x, rede)
% Jacobiano - cada elemento calculado individualmente
% Colunas: 1 = d/dt1 | 2 = d/dt3 | 3 = d/dV1 | 4 = d/dV2 | 5 = d/dV3
    t1 = x(1);  t3 = x(2);  V1 = x(3);  V2 = x(4);  V3 = x(5);  t2 = rede.t2;
    g12a = rede.g12a; b12a = rede.b12a;  g12b = rede.g12b; b12b = rede.b12b;
    g13  = rede.g13;  b13  = rede.b13;   g23  = rede.g23;  b23  = rede.b23;
    Bcap = rede.Bcap;


    t12 = t1 - t2;  c12 = cos(t12);  s12 = sin(t12);
    c21 = c12;      s21 = -s12;
    t31 = t3 - t1;  c31 = cos(t31);  s31 = sin(t31);
    c13 = c31;      s13 = -s31;
    t32 = t3 - t2;  c32 = cos(t32);  s32 = sin(t32);


    H = zeros(12, 5);


    % --- V1, V2, V3 ---
    H(1,3) = 1;
    H(2,4) = 1;
    H(3,5) = 1;


    % --- P21a ---
    H(4,1) = -V2*V1*(g12a*s21 - b12a*c21);
    H(4,2) = 0;
    H(4,3) = -V2*(g12a*c21 + b12a*s21);
    H(4,4) = 2*V2*g12a - V1*(g12a*c21 + b12a*s21);
    H(4,5) = 0;


    % --- P21b ---
    H(5,1) = -V2*V1*(g12b*s21 - b12b*c21);
    H(5,2) = 0;
    H(5,3) = -V2*(g12b*c21 + b12b*s21);
    H(5,4) = 2*V2*g12b - V1*(g12b*c21 + b12b*s21);
    H(5,5) = 0;


    % --- P31 ---
    H(6,1) = -V3*V1*(g13*s31 - b13*c31);
    H(6,2) =  V3*V1*(g13*s31 - b13*c31);
    H(6,3) = -V3*(g13*c31 + b13*s31);
    H(6,4) = 0;
    H(6,5) = 2*V3*g13 - V1*(g13*c31 + b13*s31);


    % --- P32 ---
    H(7,1) = 0;
    H(7,2) = V3*V2*(g23*s32 - b23*c32);
    H(7,3) = 0;
    H(7,4) = -V3*(g23*c32 + b23*s32);
    H(7,5) = 2*V3*g23 - V2*(g23*c32 + b23*s32);


    % --- Q12a ---
    H(8,1) = -V1*V2*(g12a*c12 + b12a*s12);
    H(8,2) = 0;
    H(8,3) = -2*V1*b12a - V2*(g12a*s12 - b12a*c12);
    H(8,4) = -V1*(g12a*s12 - b12a*c12);
    H(8,5) = 0;


    % --- Q32 ---
    H(9,1) = 0;
    H(9,2) = -V3*V2*(g23*c32 + b23*s32);
    H(9,3) = 0;
    H(9,4) = -V3*(g23*s32 - b23*c32);
    H(9,5) = -2*V3*b23 - V2*(g23*s32 - b23*c32);


    % --- Qgbat ---
    H(10,3) = 2*Bcap*V1;


    % --- P1 virtual = P12a + P12b + P13 ---
    H(11,1) = V1*V2*(g12a*s12 - b12a*c12) + V1*V2*(g12b*s12 - b12b*c12) ...
            + V1*V3*(g13*s13 - b13*c13);
    H(11,2) = -V1*V3*(g13*s13 - b13*c13);
    H(11,3) = (2*V1*g12a - V2*(g12a*c12 + b12a*s12)) ...
            + (2*V1*g12b - V2*(g12b*c12 + b12b*s12)) ...
            + (2*V1*g13  - V3*(g13*c13  + b13*s13));
    H(11,4) = -V1*(g12a*c12 + b12a*s12) - V1*(g12b*c12 + b12b*s12);
    H(11,5) = -V1*(g13*c13 + b13*s13);


    % --- Q1 virtual = Q12a + Q12b + Q13 - Bcap*V1^2 ---
    H(12,1) = -V1*V2*(g12a*c12 + b12a*s12) - V1*V2*(g12b*c12 + b12b*s12) ...
            - V1*V3*(g13*c13 + b13*s13);
    H(12,2) =  V1*V3*(g13*c13 + b13*s13);
    H(12,3) = (-2*V1*b12a - V2*(g12a*s12 - b12a*c12)) ...
            + (-2*V1*b12b - V2*(g12b*s12 - b12b*c12)) ...
            + (-2*V1*b13  - V3*(g13*s13  - b13*c13)) - 2*Bcap*V1;
    H(12,4) = -V1*(g12a*s12 - b12a*c12) - V1*(g12b*s12 - b12b*c12);
    H(12,5) = -V1*(g13*s13 - b13*c13);
end


function c = chi2crit(k, alfa)
    if exist('chi2inv', 'file')
        c = chi2inv(1 - alfa, k);
    else                                 % tabela chi^2 a 95% (k = 1..10)
        t = [3.84 5.99 7.81 9.49 11.07 12.59 14.07 15.51 16.92 18.31];
        c = t(k);
    end
end


function mostra_teste_chi2(R, alfa)
    fprintf('Convergiu em %d iteracoes.\n', R.nit);
    fprintf('Estado: t1 = %.4f graus, t3 = %.4f graus, V1 = %.4f, V2 = %.4f, V3 = %.4f pu\n', ...
            rad2deg(R.x_final(1)), rad2deg(R.x_final(2)), R.x_final(3:5));
    fprintf('J(x) = %.2f | k = %d | chi2(%d; %.0f%%) = %.2f -> ', ...
            R.J_final, R.k, R.k, 100*(1-alfa), R.chi2_crit);
    if R.J_final >= R.chi2_crit
        fprintf('existem dados errados (%.0f%% de confianca)\n', 100*(1-alfa));
    else
        fprintf('nao se pode afirmar que existem dados errados\n');
    end
end
```
