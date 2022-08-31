clc;
close all;
%%%%%% DECLARTION %%%%%%%%%%%%%%%
mat_bit = zeros(24); %512 QAM space matrix representation in binary
mat_dec = zeros(24); %512 QAM space matrix representation in decimal
phase_mat = zeros(24);

SNR_dB=[30, 10, 5, 0, -3];
SNR=[0,0,0,0,0];
a(5)=0;
No=2e-9;
E_bit(5)=0;
eff = 4.5;
Eb_No_dB =[-10, -5, 0, 5, 10, 15, 20];
Eb_No =[0.1, 0.316227766, 1, 3.16227766, 10, 31.6227766, 100];
for i=1:length(SNR_dB)
    x = SNR_dB(i)/10;
    SNR(i) = 10^(x);
    E_bit(i)= (No*SNR(i))/eff;
    a(i) = sqrt((27*E_bit(i))/1022);
end
    
no =200; %no. of bits
input_stream= randi([0 1],no,9);


for i =1:no           
    A=(input_stream(i,:));
    for ii = 1:length(A)
        if ii == 1      
            C=num2str(A(ii));
        else 
            C = cat(2,C,num2str(A(ii))); 
        end
    end
    y(i,1)= str2num(C);
end

channel_in(1,no)= 0;
for ind=1:1 %5
row=20;
columns=1;
counter=0;

re= -23;
img= -15;
h=1;
for i =1:384           %Grey code for the first 384 number
    A=dec2bin(i-1,9); %Transforms the decimal number into 10-bit binary number
    for ii = 1:length(A) %Every bit in the same number is xor-ed with the previous bit
        if ii == 1      %The MSB is kept the same
            C=A(ii);
        else
       G=xor(str2num(A(ii-1)),str2num(A(ii))); %Result of 2 xor-ed bits
       C = cat(2,C,num2str(G));  %Result of the 10-bit xor-ed put in one variable
        end
    end
        C_bit = str2num(C); %C in bits
        C2 = bin2dec(C); %C in decimal
                
            counter= counter +1; 
            mat_bit(row,columns)=C_bit;
            mat_dec(row,columns)=C2;
            phase_mat(row,columns)= re + 1i*img;
            
            for q=1:no
                if y(q,1) == C_bit
                channel_in(1,h)= (re*a(ind) + 1i* img*a(ind));
                h=h+1;
                end
            end
            
            if mod(counter,16)==0
                columns = columns+1;
                re = re+2;
                    if mod(columns,2)==0
                        row=5;
                        img=15;
                    else
                        row=20;
                        img=-15;
                    end
            else
                 if mod(columns,2)==0
                    row = row+1;
                    img = img-2;
                 else
                     row=row-1;
                     img=img+2;
                 end
            end     
end

row=24;
columns=5;
counter=0;
re= -15;
img= -23;
for i =385:512
    A=dec2bin(i-1,9);
    for ii = 1:length(A)
        if ii == 1
            C=A(ii);
        else
       G=xor(str2num(A(ii-1)),str2num(A(ii))); %#ok<*ST2NM>
       C = cat(2,C,num2str(G));
        end
    end
        C_bit = str2num(C);
        C2 = bin2dec(C);
           counter= counter +1; 
            mat_bit(row,columns)=C_bit;
            mat_dec(row,columns)=C2;
            phase_mat(row,columns)= re + 1i*img;
            
             for q=1:no
                if y(q,1) == C_bit
                channel_in(1,h)= (re*a(ind) + 1i* img*a(ind));
                
                h=h+1;
                end
            end
            
            if mod(counter,8)==0
                columns = columns+1;
                re = re+2;
                    if mod(columns,2)==0
                        row=1;
                        img=23;
                    else
                        row=24;
                        img=-23;
                    end
            else
                 if mod(columns,2)==0
                    row = row+1;
                    img= img-2;
                        if row==5
                            row=21;
                            img=-17;
                        end
                 else
                     row=row-1;
                     img=img+2;
                        if row == 20
                            row=4;
                            img=17;
                        end
                 end
            end     
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scatterplot(channel_in)
grid on;
%%%%%% CHANNEL %%%%%%%%
count = 1;
 for i=1:no
     for j=1:20
noise = sqrt((E_bit(ind)*9)/(2*SNR(ind)))*randn(size(channel_in(i)))+1i*sqrt((E_bit(ind)*9)/(2*SNR(ind)))*randn(size(channel_in(i)));
signal_wn(count) = noise + channel_in(i);
count =count+1;
     end
 end
scatterplot(signal_wn)
grid on;

%%%%%%%%%%%%%%%%%%%%DEMAPPER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

counter2 =1;
min_length=0;
for rt=1:length(Eb_No)
    signal_wn(rt) = channel_in(i) + noise/sqrt(Eb_No(rt)); 
for p = 1:length(signal_wn)
    for f=1:length(channel_in)
        d(f)= signal_wn(p)-channel_in(f);
        if f ==1
            minn=1;
            min_length=abs(d(f));
        else
            if abs(d(f))< min_length
                minn = f;
                min_length=abs(d(f));
            end
        end
    end
    received(p) = channel_in(minn);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:no 
 Rx(m) = received(m*j);
 Rx_norm(m,rt) = rdivide(received(m*j),a(ind));
end

Rx_bin(no,rt)=1;
counter3 = 1;
for k = 1:no
for m=1:24
    for n=1:24
        if phase_mat(n,m) == Rx_norm(k)
           Rx_bin(counter3,rt) = mat_bit(n,m);
           counter3 = counter3+1;
        end
    end
end
end
right(7)=0;

for mx=1:no
    for my=1:length(y)
        if Rx_bin(mx,1) == y(my)
            right(ind) =right(ind)+1;
            break
        end
    end
end

SBER(ind)= right(ind)/no;
end
end
scatterplot(Rx)
grid on

%%%%%%% BER %%%%%%%%
 %i = 1:length(Eb_No_dB)
    Diff(i)=y(i)-Rx_bin(i);
    t(i) = find(Diff);
    %F(i) = size(t);
    %SBER = F/9;
    %end
M=512; %% Symbols number
N=log2(M); % Number of Bits Per Symbol

Eb_No_db=-10:1:20;

% BER THEORETICAL

BER_theoretical = (4/N)*(1-(1/sqrt(M)))*(qfunc(sqrt((3*N)/(M-1))*sqrt(10.^(Eb_No_db/10))));

figure(1);
semilogy(Eb_No_db,BER_theoretical);
title('BER THEORETICAL');
ylabel('BER');
xlabel('Eb/No (dB)');
grid;
