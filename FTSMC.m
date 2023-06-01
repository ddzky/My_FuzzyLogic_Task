clear, clc;

%Import data
T = readtable('PDRB.xlsx');
pdrb = T.PDRB;
%Mencari maximum dan minimum pdrb
u_min = min(pdrb);
u_max = max(pdrb);
disp("u_min = " + num2str(u_min));
disp("u_max = " + num2str(u_max));
u1 = input("Masukkan u1 : ");
u2 = input("Masukkan u2 : ");

%Membagi himpunan semesta U dengan rumus Sturges
n = round(1 + 3.322*log10(length(pdrb)));
%Menghitung panjang interval
l = ((u_max+u2) - (u_min-u1))/n;

%Membuat partisi dari himpunan semesta
U_par = [];
m = []; %Matriks berisikan nilai tengah dari partisi himpunan semesta
a = u_min-u1;
b = a + l;
for i = 1:n
    if i==1
        U_par = [fixed.Interval(a,b)];
        m = [(a+b)/2];
    else
        a = b;
        b = b+l;
        U_par = [U_par;fixed.Interval(a,b)];
        m = [m;(a+b)/2];
    end
end

%Fuzzifikasi
for i = 1:length(pdrb)
    for j = 1:length(U_par)
        if contains(U_par(j,1),pdrb(i,1))
            if i == 1
                t = [i];
                DataFuzzy = ["A"+ num2str(j)];
            else
                t = [t;i];
                DataFuzzy = [DataFuzzy;"A" + num2str(j)];
            end
        end
    end
end
T1 = renamevars(table(t,pdrb,DataFuzzy),["t","pdrb","DataFuzzy"],["Data Ke-","PDRB","Fuzzy Set"]);

%FLR (Fuzzy Logical Relationship)
A = zeros(n); %Matriks transisi
for i = 1:length(pdrb)-1
    if i == 1
        urutanData = [num2str(T1.("Data Ke-")(i)) + "-" + num2str(T1.("Data Ke-")(i+1))];
        FLR = [T1.("Fuzzy Set")(i) + "-->" + T1.("Fuzzy Set")(i+1)];
        A(str2num(extractAfter(T1.("Fuzzy Set")(i),"A")),str2num(extractAfter(T1.("Fuzzy Set")(i+1),"A"))) = A(str2num(extractAfter(T1.("Fuzzy Set")(i),"A")),str2num(extractAfter(T1.("Fuzzy Set")(i+1),"A"))) + 1;
    else
        urutanData = [urutanData;num2str(T1.("Data Ke-")(i)) + "-" + num2str(T1.("Data Ke-")(i+1))];
        FLR = [FLR;T1.("Fuzzy Set")(i) + "-->" + T1.("Fuzzy Set")(i+1)];
        A(str2num(extractAfter(T1.("Fuzzy Set")(i),"A")),str2num(extractAfter(T1.("Fuzzy Set")(i+1),"A"))) = A(str2num(extractAfter(T1.("Fuzzy Set")(i),"A")),str2num(extractAfter(T1.("Fuzzy Set")(i+1),"A"))) + 1;
    end
end
T2 = renamevars(table(urutanData,FLR),"urutanData","Urutan Data");

%Membuat Tabel FLRG
FLRG = [];
rentangData = [];
for i=1:length(T2.("Urutan Data"))
    if i==1
        FLRG = [T2.FLR(i,1)];
        rentangData = [T2.("Urutan Data")(i,1)];
    else
        lhs1 = extractBefore(T2.FLR(i-1,1),"-->");
        rhs1 = extractAfter(T2.FLR(i-1,1),"-->");

        lhs2 = extractBefore(T2.FLR(i,1),"-->");
        rhs2 = extractAfter(T2.FLR(i,1),"-->");

        if lhs1 == lhs2 && rhs1 == rhs2
            rentangData(length(rentangData),1) = strrep(rentangData(length(rentangData),1),extractAfter(rentangData(length(rentangData),1),"-"),extractAfter(T2.("Urutan Data")(i,1),"-"));
            continue;
        elseif lhs1 == lhs2 && rhs1 ~= rhs2
            rentangData(length(rentangData),1) = strrep(rentangData(length(rentangData),1),extractAfter(rentangData(length(rentangData),1),"-"),extractAfter(T2.("Urutan Data")(i,1),"-"));
            FLRG(length(FLRG),1) = insertAfter(FLRG(length(FLRG),1),strlength(FLRG(length(FLRG),1)),insertBefore(rhs2,1,","));
        elseif lhs1 ~= lhs2
            rentangData = [rentangData;T2.("Urutan Data")(i,1)];
            FLRG = [FLRG;T2.FLR(i,1)];
        end
    end
end

T3 = renamevars(table(rentangData,FLRG),"rentangData","Data Ke-");


%Matriks Probabilitas
P = zeros(n);
for i=1:n
    for j=1:n
        P(i,j) = A(i,j)/sum(A(i,:));
    end
end

%Peramalan Awal F(t)
F=[];
for i=1:length(pdrb)
    if i==1
        F = [0];
    else
        j=1;
        lhs1 = str2num(extractBefore(T2.("Urutan Data")(i-1,1),"-"));
        rhs1 = str2num(extractAfter(T2.("Urutan Data")(i-1,1),"-"));

        %Mencari indeks dari A_i
        lhs3 = str2num(extract(T2.FLR(i-1,1),2));
        rhs3 = str2num(extractAfter(T2.FLR(i-1,1),6));
        %Mencari FLRG
        control = true;
        while control
            lhs2 = str2num(extractBefore(T3.("Data Ke-")(j,1),"-"));
            rhs2 = str2num(extractAfter(T3.("Data Ke-")(j,1),"-"));
            if lhs1 >= lhs2 && rhs1 <= rhs2
                control = false;
            else
                j=j+1;
            end
        end

        %rule 1 One-to-One
        if strlength(T3.FLRG(j,1)) == 7
            F = [F;m(rhs3,1)];
            
        else%rule 2 One-to-Many
            m_sementara = m;
            m_sementara(lhs3,1) = T.PDRB(i-1,1);
            F = [F ; P(lhs3,:)*m_sementara];
            
        end
    end
end
T.("Peramalan Awal F(t)") = F;

%Adjusting data aktual dengan hasil peramalan pertama
D = zeros(length(pdrb),2);
for i=1:length(pdrb)
    if i==1
        continue;
    else
        lhs = str2num(extractAfter(T1.("Fuzzy Set")(i-1,1),"A"));
        rhs = str2num(extractAfter(T1.("Fuzzy Set")(i,1),"A"));
        %Rule 1 dan 2
        if  lhs < rhs && A(lhs,rhs) > 0 && A(rhs,lhs) > 0
            D(i,1) = (l/2);
        elseif lhs > rhs && A(lhs,rhs) > 0 && A(rhs,lhs) > 0
            D(i,1) = -(l/2);
        end
        %Rule 3 dan 4
        if lhs ~= rhs
            D(i,2) = (l/2)*(rhs - lhs);
        end
    end
end
T.Penyesuaian = D(:,1) + D(:,2);
T.("Peramalan Akhir F*(t)") = T.("Peramalan Awal F(t)") + T.Penyesuaian;

for i=1:length(T.PDRB)
    if abs(T.PDRB(i,1) - T.("Peramalan Awal F(t)")(i,1)) < abs(T.PDRB(i,1) - T.("Peramalan Akhir F*(t)")(i,1))
        T.("Peramalan Akhir F*(t)")(i,1) = T.("Peramalan Awal F(t)")(i,1);
    end
end

%Peramalan


%Mape untuk F(t) dan RMSE
Mape1 = mean(abs((T.PDRB-T.("Peramalan Awal F(t)"))./T.PDRB))*100;
rmse1 = sqrt(mean((T.PDRB-T.("Peramalan Awal F(t)")).^2));
disp("Mape untuk F(t) = " + num2str(Mape1));
disp("RMSE untuk F(t) = " + num2str(rmse1));
%Mape untuk F*(t) dan RMSE
Mape2 = mean(abs((T.PDRB-T.("Peramalan Akhir F*(t)"))./T.PDRB))*100;
rmse2 = sqrt(mean((T.PDRB-T.("Peramalan Akhir F*(t)")).^2));
disp("Mape untuk F*(t) = " + num2str(Mape2));
disp("RMSE untuk F*(t) = "+ num2str(rmse2));
x = 1:1:length(pdrb);
plot(x,pdrb,x,F,x,T.("Peramalan Akhir F*(t)"));
legend("Data Actual","Sebelum Penyesuaian","Setelah Penyesuaian");