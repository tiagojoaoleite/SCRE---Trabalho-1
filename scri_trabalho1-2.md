# SCRI 2026/27 \- TP1: Estimacao de Estado em Redes de Energia (WLS)

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


% --- Angulos ja conhecidos ---
t2 = 0;                  % barramento 2 = slack (referencia). t1 e t3 sao estimados.


% --- Parametros do metodo ---
tol   = 1e-4;            % tolerancia em max|dx|
maxit = 50;
alfa  = 0.05;            % nivel de significancia do teste chi^2
limiar_rN = 3;           % limiar para |rN|


% Estrutura com os parametros da rede (passada as funcoes)
rede.g12a = real(y12a);  rede.b12a = imag(y12a);
rede.g12b = real(y12b);  rede.b12b = imag(y12b);
rede.g13  = real(y13);   rede.b13  = imag(y13);
rede.g23  = real(y23);   rede.b23  = imag(y23);
rede.Bcap = Bcap;        rede.t2   = t2;


```

# ================== MATRIZ Ybus ==================
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

# ================ VETOR DE ESTADO ================
```matlab


% x = [t1, t3, V1, V2, V3]'   (angulos em rad)
StateNames  = {'t1', 't3', 'V1', 'V2', 'V3'};
StateNamesD = {'t1_graus', 't3_graus', 'V1_pu', 'V2_pu', 'V3_pu'};   % para as tabelas
n = numel(StateNames);


```

# ================ VETOR DE MEDIDAS ================
```matlab


MeasNames = {'V1','V2','V3','P21a','P21b','P31','P32','Q12a','Q32','Qgbat','P1v','Q1v'};
z   = [V1; V2; V3; P21a; P21b; P31; P32; Q12a; Q32; Qgbat; P1v; Q1v];
var = [varV1; varV2; varV3; varP21a; varP21b; varP31; varP32; ...
       varQ12a; varQ32; varQgbat; varP1v; varQ1v];
virtual = ismember(MeasNames, {'P1v','Q1v'})';


MeasCols = {'z', 'variancia'};
disp_medidas = array2table([z, var], 'RowNames', MeasNames, 'VariableNames', MeasCols)
```

| |z|variancia|
|:--:|:--:|:--:|
|1 V1|1.1000|2.0000e-04|
|2 V2|1|2.0000e-04|
|3 V3|1.0700|2.0000e-04|
|4 P21a|-0.7900|1.0000e-03|
|5 P21b|0.7900|1.0000e-03|
|6 P31|1.7100|1.0000e-03|
|7 P32|2.3400|1.0000e-03|
|8 Q12a|0.9700|0.0050|
|9 Q32|1.2200|0.0050|
|10 Qgbat|2.6600|0.0050|
|11 P1v|0|1.0000e-05|
|12 Q1v|0|1.0000e-05|

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

```matlab


% Nomes de colunas comuns
JCols    = {'J', 'k', 'chi2_crit', 'max_abs_dx', 'dados_errados'};   % dados_errados: 1 = sim, 0 = nao
TesteCols = {'J', 'k', 'chi2_crit', 'dados_errados', 'iteracoes'};
ativas_P1 = true(numel(z), 1);
P1 = gauss_newton(z, var, ativas_P1, rede, tol, maxit, alfa);
IterNames_P1 = P1.IterNames;
MeasNames_P1 = MeasNames(ativas_P1);


% Tabela de estado por iteracao (angulos em graus)
disp_P1_estado = array2table(P1.X_hist, 'RowNames', StateNamesD, 'VariableNames', IterNames_P1)
```

| |It0|It1|It2|It3|It4|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1 t1_graus|0|2.4736|2.5004|2.5041|2.5042|
|2 t3_graus|0|10.3339|9.3740|9.3781|9.3786|
|3 V1_pu|1|1.0944|1.0890|1.0890|1.0890|
|4 V2_pu|1|0.9636|0.9733|0.9735|0.9735|
|5 V3_pu|1|1.1048|1.0997|1.0995|1.0995|

```matlab


% Tabela de residuos por iteracao (r = z - h(x))
disp_P1_residuos = array2table(P1.r_hist, 'RowNames', MeasNames_P1, 'VariableNames', IterNames_P1)
```

| |It0|It1|It2|It3|It4|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1 V1|0.1000|0.0056|0.0110|0.0110|0.0110|
|2 V2|0|0.0364|0.0267|0.0265|0.0265|
|3 V3|0.0700|-0.0348|-0.0297|-0.0295|-0.0295|
|4 P21a|-0.7900|-0.0846|-0.1043|-0.1039|-0.1039|
|5 P21b|0.7900|1.4954|1.4757|1.4761|1.4761|
|6 P31|1.7100|0.0109|0.2359|0.2365|0.2365|
|7 P32|2.3400|0.0839|0.2900|0.2900|0.2900|
|8 Q12a|0.9700|-0.3809|-0.2083|-0.2054|-0.2053|
|9 Q32|1.2200|-0.1304|0.0350|0.0400|0.0403|
|10 Qgbat|0.4600|0.0248|0.0508|0.0509|0.0509|
|11 P1v|0|0.1664|0.0063|0.0053|0.0053|
|12 Q1v|2.2000|-0.3978|-0.0038|1.3726e-04|1.3753e-04|

```matlab


% Tabela J(x) por iteracao e teste chi^2
disp_P1_J = array2table(P1.J_tab, 'RowNames', IterNames_P1, 'VariableNames', JCols)
```

| |J|k|chi2_crit|max_abs_dx|dados_errados|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1 It0|4.9425e+05|7|14.0671|NaN|1|
|2 It1|2.0889e+04|7|14.0671|0.1804|1|
|3 It2|2.3518e+03|7|14.0671|0.0168|1|
|4 It3|2.3501e+03|7|14.0671|2.3929e-04|1|
|5 It4|2.3501e+03|7|14.0671|1.6623e-05|1|

```matlab


% Tabela covariancia dos residuos (Omega_ii) e rN por iteracao
disp_P1_rN = array2table(P1.OmrN_tab, 'RowNames', MeasNames_P1, 'VariableNames', P1.OmrNCols)
```

| |Omega_It0|rN_It0|Omega_It1|rN_It1|Omega_It2|rN_It2|Omega_It3|rN_It3|Omega_It4|rN_It4|
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
|1 V1|1.3923e-04|8.4748|1.4590e-04|0.4603|1.4531e-04|0.9103|1.4530e-04|0.9117|1.4530e-04|0.9118|
|2 V2|1.4776e-04|0|1.4434e-04|3.0305|1.4478e-04|2.2223|1.4479e-04|2.2023|1.4480e-04|2.2017|
|3 V3|1.4465e-04|5.8203|1.5080e-04|-2.8314|1.5047e-04|-2.4227|1.5046e-04|-2.4075|1.5045e-04|-2.4061|
|4 P21a|9.3226e-04|-25.8737|9.3299e-04|-2.7709|9.3324e-04|-3.4136|9.3325e-04|-3.4026|9.3325e-04|-3.4022|
|5 P21b|9.3226e-04|25.8737|9.3299e-04|48.9563|9.3324e-04|48.3067|9.3325e-04|48.3173|9.3325e-04|48.3177|
|6 P31|7.3237e-04|63.1875|6.9963e-04|0.4125|7.0371e-04|8.8945|7.0374e-04|8.9158|7.0374e-04|8.9151|
|7 P32|4.0032e-04|116.9532|4.3054e-04|4.0432|4.2577e-04|14.0528|4.2571e-04|14.0577|4.2571e-04|14.0566|
|8 Q12a|0.0045|14.5078|0.0044|-5.7190|0.0044|-3.1261|0.0044|-3.0820|0.0044|-3.0805|
|9 Q32|9.1109e-04|40.4184|7.9481e-04|-4.6270|8.0187e-04|1.2376|8.0236e-04|1.4134|8.0239e-04|1.4235|
|10 Qgbat|0.0038|7.4392|0.0037|0.4059|0.0037|0.8309|0.0037|0.8323|0.0037|0.8324|
|11 P1v|3.3245e-08|0|3.3317e-08|911.4281|3.3211e-08|34.4381|3.3210e-08|28.8149|3.3210e-08|28.8188|
|12 Q1v|2.9027e-09|4.0834e+04|2.9612e-09|-7.3106e+03|2.9401e-09|-70.9116|2.9396e-09|2.5316|2.9396e-09|2.5366|

```matlab


% a) Valores estimados das grandezas medidas
EstCols = {'z', 'z_estimado', 'residuo'};
disp_P1_estimados = array2table([z, P1.h_final, z - P1.h_final], ...
                                'RowNames', MeasNames, 'VariableNames', EstCols)
```

| |z|z_estimado|residuo|
|:--:|:--:|:--:|:--:|
|1 V1|1.1000|1.0890|0.0110|
|2 V2|1|0.9735|0.0265|
|3 V3|1.0700|1.0995|-0.0295|
|4 P21a|-0.7900|-0.6861|-0.1039|
|5 P21b|0.7900|-0.6861|1.4761|
|6 P31|1.7100|1.4735|0.2365|
|7 P32|2.3400|2.0500|0.2900|
|8 Q12a|0.9700|1.1753|-0.2053|
|9 Q32|1.2200|1.1797|0.0403|
|10 Qgbat|2.6600|2.6091|0.0509|
|11 P1v|0|-0.0053|0.0053|
|12 Q1v|0|-1.3753e-04|1.3753e-04|

```matlab


% b) Teste chi^2
disp_P1_teste = array2table(P1.teste, 'RowNames', {'Pergunta1'}, 'VariableNames', TesteCols)
```

| |J|k|chi2_crit|dados_errados|iteracoes|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1 Pergunta1|2.3501e+03|7|14.0671|1|4|

```matlab
% Medidas ordenadas por |rN| (resultado final da Pergunta 1)
idxA = find(ativas_P1);
[~, ord] = sort(abs(P1.rN_final), 'descend');
OrdNames = MeasNames(idxA(ord));
OrdCols  = {'abs_rN', 'virtual'};                       % virtual: 1 = sim
disp_P2_ordenadas = array2table([abs(P1.rN_final(ord)), double(virtual(idxA(ord)))], ...
                                'RowNames', OrdNames, 'VariableNames', OrdCols)
```

| |abs_rN|virtual|
|:--:|:--:|:--:|
|1 P21b|48.3177|0|
|2 P1v|28.8188|1|
|3 P32|14.0566|0|
|4 P31|8.9151|0|
|5 P21a|3.4022|0|
|6 Q12a|3.0805|0|
|7 Q1v|2.5366|1|
|8 V3|2.4061|0|
|9 V2|2.2017|0|
|10 Q32|1.4235|0|
|11 V1|0.9118|0|
|12 Qgbat|0.8324|0|

```matlab


% Candidata a remocao: maior |rN| entre as medidas REAIS (virtuais sao exatas)
rN_reais = abs(P1.rN_final);  rN_reais(virtual(idxA)) = -Inf;
[rmax, q] = max(rN_reais);
cand_P2 = idxA(q);
fprintf('Candidata a remocao: %s (|rN| = %.2f > %.1f)\n', MeasNames{cand_P2}, rmax, limiar_rN);
```

```matlabTextOutput
Candidata a remocao: P21b (|rN| = 48.32 > 3.0)
```

```matlab


% Nova estimacao sem a candidata
ativas_P2 = ativas_P1;  ativas_P2(cand_P2) = false;
P2 = gauss_newton(z, var, ativas_P2, rede, tol, maxit, alfa);
IterNames_P2 = P2.IterNames;
MeasNames_P2 = MeasNames(ativas_P2);


% Tabela de estado por iteracao (angulos em graus)
disp_P2_estado = array2table(P2.X_hist, 'RowNames', StateNamesD, 'VariableNames', IterNames_P2)
```

| |It0|It1|It2|It3|It4|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1 t1_graus|0|3.0704|3.0782|3.0812|3.0812|
|2 t3_graus|0|12.0930|10.9976|10.9980|10.9980|
|3 V1_pu|1|1.0940|1.0911|1.0912|1.0912|
|4 V2_pu|1|0.9614|0.9782|0.9785|0.9785|
|5 V3_pu|1|1.1081|1.0990|1.0989|1.0989|

```matlab


% Tabela de residuos por iteracao
disp_P2_residuos = array2table(P2.r_hist, 'RowNames', MeasNames_P2, 'VariableNames', IterNames_P2)
```

| |It0|It1|It2|It3|It4|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1 V1|0.1000|0.0060|0.0089|0.0088|0.0088|
|2 V2|0|0.0386|0.0218|0.0215|0.0215|
|3 V3|0.0700|-0.0381|-0.0290|-0.0289|-0.0289|
|4 P21a|-0.7900|0.0254|8.6883e-04|0.0013|0.0013|
|5 P31|1.7100|-0.2522|0.0175|0.0184|0.0184|
|6 P32|2.3400|-0.2641|-0.0159|-0.0156|-0.0156|
|7 Q12a|0.9700|-0.3838|-0.1625|-0.1598|-0.1598|
|8 Q32|1.2200|-0.1956|0.1045|0.1095|0.1096|
|9 Qgbat|0.4600|0.0269|0.0411|0.0406|0.0406|
|10 P1v|0|0.1885|0.0016|7.5338e-05|7.5835e-05|
|11 Q1v|2.2000|-0.4511|-0.0039|8.7562e-05|8.7137e-05|

```matlab


% Tabela J(x) por iteracao e teste chi^2
disp_P2_J = array2table(P2.J_tab, 'RowNames', IterNames_P2, 'VariableNames', JCols)
```

| |J|k|chi2_crit|max_abs_dx|dados_errados|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1 It0|4.9363e+05|6|12.5916|NaN|1|
|2 It1|2.4092e+04|6|12.5916|0.2111|1|
|3 It2|17.1508|6|12.5916|0.0191|1|
|4 It3|15.3038|6|12.5916|3.5504e-04|1|
|5 It4|15.3038|6|12.5916|2.5970e-06|1|

```matlab


% Tabela covariancia dos residuos (Omega_ii) e rN por iteracao
disp_P2_rN = array2table(P2.OmrN_tab, 'RowNames', MeasNames_P2, 'VariableNames', P2.OmrNCols)
```

| |Omega_It0|rN_It0|Omega_It1|rN_It1|Omega_It2|rN_It2|Omega_It3|rN_It3|Omega_It4|rN_It4|
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
|1 V1|1.3923e-04|8.4748|1.4590e-04|0.4951|1.4534e-04|0.7415|1.4534e-04|0.7326|1.4534e-04|0.7326|
|2 V2|1.4775e-04|0|1.4471e-04|3.2125|1.4538e-04|1.8104|1.4539e-04|1.7809|1.4539e-04|1.7808|
|3 V3|1.4464e-04|5.8203|1.5078e-04|-3.0993|1.5003e-04|-2.3695|1.5001e-04|-2.3628|1.5001e-04|-2.3626|
|4 P21a|9.2734e-04|-25.9423|9.2944e-04|0.8331|9.2983e-04|0.0285|9.2985e-04|0.0417|9.2985e-04|0.0417|
|5 P31|7.1363e-04|64.0115|6.7814e-04|-9.6845|6.8308e-04|0.6701|6.8312e-04|0.7056|6.8312e-04|0.7056|
|6 P32|3.5746e-04|123.7663|3.9069e-04|-13.3590|3.8514e-04|-0.8120|3.8508e-04|-0.7964|3.8508e-04|-0.7963|
|7 Q12a|0.0045|14.5078|0.0044|-5.7727|0.0044|-2.4426|0.0044|-2.4020|0.0044|-2.4018|
|8 Q32|9.1107e-04|40.4187|7.9465e-04|-6.9385|8.1432e-04|3.6636|8.1478e-04|3.8371|8.1479e-04|3.8383|
|9 Qgbat|0.0038|7.4392|0.0037|0.4389|0.0037|0.6719|0.0037|0.6635|0.0037|0.6635|
|10 P1v|2.1404e-08|0|2.1839e-08|1.2757e+03|2.1679e-08|10.7616|2.1676e-08|0.5117|2.1676e-08|0.5151|
|11 Q1v|2.9027e-09|4.0834e+04|2.9910e-09|-8.2488e+03|2.9476e-09|-72.3237|2.9470e-09|1.6130|2.9470e-09|1.6051|

```matlab


% Teste chi^2
disp_P2_teste = array2table(P2.teste, 'RowNames', {'Pergunta2'}, 'VariableNames', TesteCols)
```

| |J|k|chi2_crit|dados_errados|iteracoes|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1 Pergunta2|15.3038|6|12.5916|1|4|

```matlab


ativas_P3 = ativas_P2;
P3 = P2;
removidas = {MeasNames{cand_P2}};
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
            MeasNames{j}, rmax, P3.J_final, P3.chi2_crit);
    if nnz(ativas_P3) - 1 <= n
        warning('Remover %s tornaria o sistema nao observavel (m <= n).', MeasNames{j}); break
    end
    ativas_P3(j) = false;  removidas{end+1} = MeasNames{j}; %#ok<SAGROW>
    P3 = gauss_newton(z, var, ativas_P3, rede, tol, maxit, alfa);
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
IterNames_P3 = P3.IterNames;
MeasNames_P3 = MeasNames(ativas_P3);


% Tabela de estado por iteracao (angulos em graus)
disp_P3_estado = array2table(P3.X_hist, 'RowNames', StateNamesD, 'VariableNames', IterNames_P3)
```

| |It0|It1|It2|It3|It4|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1 t1_graus|0|3.3455|3.1535|3.1545|3.1545|
|2 t3_graus|0|12.9722|11.5429|11.5465|11.5465|
|3 V1_pu|1|1.1051|1.0981|1.0980|1.0980|
|4 V2_pu|1|0.9961|1.0013|1.0014|1.0014|
|5 V3_pu|1|1.0672|1.0733|1.0735|1.0735|

```matlab


% Tabela de residuos por iteracao
disp_P3_residuos = array2table(P3.r_hist, 'RowNames', MeasNames_P3, 'VariableNames', IterNames_P3)
```

| |It0|It1|It2|It3|It4|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1 V1|0.1000|-0.0051|0.0019|0.0020|0.0020|
|2 V2|0|0.0039|-0.0013|-0.0014|-0.0014|
|3 V3|0.0700|0.0028|-0.0033|-0.0035|-0.0035|
|4 P21a|-0.7900|0.0659|0.0053|0.0053|0.0053|
|5 P31|1.7100|-0.2147|0.0184|0.0171|0.0171|
|6 P32|2.3400|-0.2524|-0.0085|-0.0101|-0.0101|
|7 Q12a|0.9700|-0.1255|0.0119|0.0135|0.0135|
|8 Qgbat|0.4600|-0.0268|0.0074|0.0076|0.0076|
|9 P1v|0|0.0811|-0.0014|7.4114e-05|7.4115e-05|
|10 Q1v|2.2000|-0.4834|-0.0054|-2.4834e-05|-2.4580e-05|

```matlab


% Tabela J(x) por iteracao e teste chi^2
disp_P3_J = array2table(P3.J_tab, 'RowNames', IterNames_P3, 'VariableNames', JCols)
```

| |J|k|chi2_crit|max_abs_dx|dados_errados|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1 It0|4.9333e+05|5|11.0705|NaN|1|
|2 It1|2.4143e+04|5|11.0705|0.2264|1|
|3 It2|3.6353|5|11.0705|0.0249|0|
|4 It3|0.5609|5|11.0705|2.0813e-04|0|
|5 It4|0.5609|5|11.0705|1.0976e-06|0|

```matlab


% Tabela covariancia dos residuos (Omega_ii) e rN por iteracao
disp_P3_rN = array2table(P3.OmrN_tab, 'RowNames', MeasNames_P3, 'VariableNames', P3.OmrNCols)
```

| |Omega_It0|rN_It0|Omega_It1|rN_It1|Omega_It2|rN_It2|Omega_It3|rN_It3|Omega_It4|rN_It4|
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
|1 V1|1.3558e-04|8.5882|1.4300e-04|-0.4273|1.4230e-04|0.1625|1.4229e-04|0.1660|1.4229e-04|0.1660|
|2 V2|1.1199e-04|0|1.1055e-04|0.3757|1.1099e-04|-0.1244|1.1099e-04|-0.1330|1.1099e-04|-0.1330|
|3 V3|9.5120e-05|7.1773|1.0612e-04|0.2688|1.0542e-04|-0.3218|1.0543e-04|-0.3421|1.0543e-04|-0.3422|
|4 P21a|9.2732e-04|-25.9425|9.3036e-04|2.1592|9.3019e-04|0.1728|9.3021e-04|0.1722|9.3021e-04|0.1722|
|5 P31|7.1356e-04|64.0149|6.7756e-04|-8.2470|6.8284e-04|0.7048|6.8290e-04|0.6530|6.8290e-04|0.6531|
|6 P32|3.5729e-04|123.7959|3.8884e-04|-12.8011|3.8316e-04|-0.4352|3.8309e-04|-0.5176|3.8309e-04|-0.5175|
|7 Q12a|0.0027|18.7463|0.0024|-2.5835|0.0024|0.2444|0.0024|0.2756|0.0024|0.2755|
|8 Qgbat|0.0038|7.5090|0.0037|-0.4433|0.0037|0.1220|0.0037|0.1253|0.0037|0.1253|
|9 P1v|2.1404e-08|0|2.1400e-08|554.1565|2.1426e-08|-9.4777|2.1427e-08|0.5063|2.1427e-08|0.5063|
|10 Q1v|2.1986e-09|4.6919e+04|2.0991e-09|-1.0551e+04|2.0848e-09|-117.5757|2.0849e-09|-0.5439|2.0849e-09|-0.5383|

```matlab


% Teste chi^2
disp_P3_teste = array2table(P3.teste, 'RowNames', {'Pergunta3'}, 'VariableNames', TesteCols)
```

| |J|k|chi2_crit|dados_errados|iteracoes|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1 Pergunta3|0.5609|5|11.0705|0|4|

```matlab


% Comparacao com as medidas originais (todas, incluindo as removidas)
CompCols = {'z', 'z_estimado', 'diferenca', 'ativa'};        % ativa: 1 = usada, 0 = removida
disp_P3_comparacao = array2table([z, P3.h_final, z - P3.h_final, double(ativas_P3)], ...
                                 'RowNames', MeasNames, 'VariableNames', CompCols)
```

| |z|z_estimado|diferenca|ativa|
|:--:|:--:|:--:|:--:|:--:|
|1 V1|1.1000|1.0980|0.0020|1|
|2 V2|1|1.0014|-0.0014|1|
|3 V3|1.0700|1.0735|-0.0035|1|
|4 P21a|-0.7900|-0.7953|0.0053|1|
|5 P21b|0.7900|-0.7953|1.5853|0|
|6 P31|1.7100|1.6929|0.0171|1|
|7 P32|2.3400|2.3501|-0.0101|1|
|8 Q12a|0.9700|0.9565|0.0135|1|
|9 Q32|1.2200|0.5613|0.6587|0|
|10 Qgbat|2.6600|2.6524|0.0076|1|
|11 P1v|0|-7.4115e-05|7.4115e-05|1|
|12 Q1v|0|2.4580e-05|-2.4580e-05|1|

```matlab


% Producao dos geradores (carga no barramento 2 = 6 + j2 pu)
xf = P3.x_final;
Vc = [xf(3)*exp(1j*xf(1));             % V1 /_ t1
      xf(4)*exp(1j*t2);                % V2 /_ t2 = 0
      xf(5)*exp(1j*xf(2))];            % V3 /_ t3
Sinj = Vc .* conj(Ybus * Vc);          % injecoes liquidas
Scarga2 = 6 + 2j;
SG2 = Sinj(2) + Scarga2;
SG3 = Sinj(3);


GenNames = {'G2_slack', 'G3'};
GenCols  = {'P_pu', 'Q_pu'};
disp_geradores = array2table([real(SG2), imag(SG2); real(SG3), imag(SG3)], ...
                             'RowNames', GenNames, 'VariableNames', GenCols)
```

| |P_pu|Q_pu|
|:--:|:--:|:--:|
|1 G2_slack|2.1568|0.2660|
|2 G3|4.0431|0.0804|

```matlab


InjCols = {'P_inj_pu', 'Q_inj_pu'};
disp_injecoes = array2table([real(Sinj), imag(Sinj)], 'RowNames', BusNames, 'VariableNames', InjCols)
```

| |P_inj_pu|Q_inj_pu|
|:--:|:--:|:--:|
|1 ONE|-7.4115e-05|2.4580e-05|
|2 TWO|-3.8432|-1.7340|
|3 THREE|4.0431|0.0804|

```matlab
fprintf('Perdas ativas = %.4f pu\n', real(sum(Sinj)));
```

```matlabTextOutput
Perdas ativas = 0.1998 pu
```

```matlab


% Resumo para o relatorio
ResumoNames = {'P1_todas', 'P2', 'P3_final'};
ResumoCols  = {'J', 'k', 'chi2_crit', 'max_abs_rN'};
disp_resumo = array2table([P1.J_final, P1.k, P1.chi2_crit, max(abs(P1.rN_final));
                           P2.J_final, P2.k, P2.chi2_crit, max(abs(P2.rN_final));
                           P3.J_final, P3.k, P3.chi2_crit, max(abs(P3.rN_final))], ...
                          'RowNames', ResumoNames, 'VariableNames', ResumoCols)
```

| |J|k|chi2_crit|max_abs_rN|
|:--:|:--:|:--:|:--:|:--:|
|1 P1_todas|2.3501e+03|7|14.0671|48.3177|
|2 P2|15.3038|6|12.5916|3.8383|
|3 P3_final|0.5609|5|11.0705|0.6531|

```matlab


% Grafico de |rN| nas tres estimacoes
img = imread('Screenshot 2026-10-16 at 19.35.08.png');
```


```matlab


% =========================== FUNCOES ===========================
```matlab

function R = gauss_newton(z, var, ativas, rede, tol, maxit, alfa)
% Iteracao de Gauss-Newton:  x(k+1) = x(k) + G^-1 * H' * W * (z - h(x(k)))
% Devolve matrizes numericas; as tabelas sao construidas no script principal.
    za = z(ativas);
    W  = diag(1 ./ var(ativas));
    Rm = diag(var(ativas));
    n  = 5;  m = nnz(ativas);  k = m - n;
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


    % Nomes das colunas de iteracao
    R.IterNames = arrayfun(@(i) sprintf('It%d', i), 0:nit-1, 'UniformOutput', false);


    % Estado por iteracao (angulos em graus)
    Xd = X_hist;  Xd(1:2,:) = rad2deg(Xd(1:2,:));
    R.X_hist = Xd;


    % Residuos por iteracao
    R.r_hist = r_hist;


    % J(x) por iteracao: [J k chi2_crit max_abs_dx dados_errados]
    R.J_tab = [J_hist', repmat(k, nit, 1), repmat(chi2c, nit, 1), ...
               [NaN, dx_hist(1:nit-1)]', double(J_hist' >= chi2c)];


    % Omega_ii e rN por iteracao (colunas intercaladas)
    R.OmrN_tab = zeros(m, 2*nit);
    R.OmrNCols = cell(1, 2*nit);
    for i = 1:nit
        R.OmrN_tab(:, 2*i-1) = Om_hist(:, i);
        R.OmrN_tab(:, 2*i)   = rN_hist(:, i);
        R.OmrNCols{2*i-1} = sprintf('Omega_It%d', i-1);
        R.OmrNCols{2*i}   = sprintf('rN_It%d',    i-1);
    end


    % Resultados finais
    R.x_final  = x;               R.h_final  = h_fun(x, rede);
    R.r_final  = r_hist(:, end);  R.rN_final = rN_hist(:, end);
    R.J_final  = J_hist(end);     R.k = k;   R.chi2_crit = chi2c;
    R.G_final  = G;               R.nit = nit - 1;
    R.Omega_final = Omega;
    R.teste = [R.J_final, k, chi2c, double(R.J_final >= chi2c), R.nit];
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
```
