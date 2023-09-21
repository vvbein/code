clc;clear;
hr_zh_data;
Fa=data(:,4:5)'; 
Pst=data(:,6)'; 
Fst=data(:,7)';
PDrg=data(:,13)';
Trg1=data(:,1)'; 
Trg2=data(:,2)'; 
Tst=data(:,3)'; 
Tra1=data(:,8)'; 
Tra2=data(:,9)';
Ta=207; 
Tfg1=data(:,1)';
Tfg2=data(:,2)';
Go=data(:,10)';
Gro=data(:,11)';
Gw=data(:,12)';

Cpc=1.097; 
Cpa=1.409+1.683*0.0363 ;  
Cpo=2.72;
CpL=0.47;
Cpw=2.16;
CO_fg=0.0479; CO2_fg=0.1448; O2_fg=0.0243;
%CO_fg=[0.1062,0.0479]; CO2_fg=[0.1081,0.1448]; O2_fg=[0.0066,0.0243];第一再生器烟气分析不知道，不用了
dense_o=[7,1.87]; %假密度,提升管油气密度
dense_L=7850; %真密度，耐火衬里，基础材料为碳钢
Vr=[33.5,492.2]; %m3
L=[13.2,27.1]; %第一第二反应区
Th=0.12; %衬里厚度, m
r_in=[0.9,2.4];
%delta_H=4.186*(344-408); %油气气相焓值，标准密度dense_20=0.912,K=10.436
beta=CO2_fg./CO_fg; %CO2、CO为再生烟气体积分数
alpha=(8.93-0.425.*(CO2_fg+O2_fg)-0.257.*CO_fg)/(CO2_fg+CO_fg);
voc=22.4/12.*((beta+0.5)./(beta+1)+3.*alpha); %beta为CO2/CO mole fraction ratio，Alpha为H/C mass fraction ratio
HDcb=3.00.*10^4.*((beta+0.303)./(beta+1)+3.53.*alpha);
Gcf=zeros(1,3000);
Cv_new=zeros(1,3000);
X=zeros(2,3000);
Y=zeros(2,3000);
Z=zeros(2,3000);
for turn=0:299
    Qo(turn*10+1)=(Fa(1,turn*10+1).*(0.21-O2_fg).*HDcb)./voc+(Fa(2,turn*10+1).*(0.21-O2_fg).*HDcb)./voc;
    Qq(turn*10+1)=((592.46-2.1.*(Pst(turn*10+1)+1)./0.098+(0.48+0.0045.*(Pst(turn*10+1)+1)./0.098)*200)-166).*Fst(turn*10+1)*4.186*1000;
    Qt(turn*10+1)=0.063*Qo(turn*10+1);
    Qa(turn*10+1)=Fa(1,turn*10+1).*Cpa.*(Tfg1(turn*10+1)-Ta)+Fa(2,turn*10+1).*Cpa.*(Tfg2(turn*10+1)-Ta);
    Ql(turn*10+1)=0.03*Qo(turn*10+1);%热损失
    Qc(turn*10+1)=Qo(turn*10+1)-Qa(turn*10+1)-Qt(turn*10+1)-Qq(turn*10+1)-Ql(turn*10+1);
    Gcs(turn*10+1)=Qc(turn*10+1)./Cpc/(Trg2(turn*10+1)-Tst(turn*10+1));
    Gcf(turn*10+1)=Gcs(turn*10+1); %kg/h
    Cv_new(turn*10+1)=Gcf(turn*10+1)./(PDrg(turn*10+1).^0.5);
  
    for k=turn*10+2:turn*10+10

        X(:,k)=[Tra1(k),Trg2(k)]';
        sum_1(k)=0;
        sum_2(k)=0;
        for i=1:10
            if k-turn*10==i
                break
            end
            sum_1(k)=sum_1(k)+Tra1(k-i);
            sum_2(k)=sum_2(k)+Trg2(k-i);
        end
        Y(:,k)=[sum_1(k)./(i-1),sum_2(k)./(i-1)]';
        Z(:,k)=abs(X(:,k)-Y(:,k));
        if Z(:,k)<=[1.3,1.3]' %阈值
            %稳态ccr计算
            Qo(k)=(Fa(1,k).*(0.21-O2_fg).*HDcb)./voc+(Fa(2,k).*(0.21-O2_fg).*HDcb)./voc; %delta_Hcb燃烧1kg焦炭产生的热量kj/kg焦炭，voc为燃烧1kg焦炭耗氧量Nm3/kg焦炭
            Qq(k)=((592.46-2.1.*(Pst(k)+1)./0.098+(0.48+0.0045.*(Pst(k)+1)./0.098)*200)-166).*Fst(k)*4.186*1000; %Pst外取热器产蒸汽压力PD2150, Fst外取热器产汽量FD2131
            %Qq=Fqw.*Cpw.*(Tw2-Tw1); %Fqw、Cpw取热水流量、热容，Tw2、Tw1取热出入口温度
            Qt(k)=0.063*Qo(k); %Qt焦炭脱附热
            Qa(k)=Fa(1,k).*Cpa.*(Tfg1(k)-Ta)+Fa(2,k).*Cpa.*(Tfg2(k)-Ta); %Ta主风入口温度,Fa主风流量Nm3/h  Cpa为主风热容(乘了dens_air kg/Nm3)kj/Nm3/℃
            Ql(k)=0.03*Qo(k);%热损失
            Qc(k)=Qo(k)-Qa(k)-Qt(k)-Qq(k)-Ql(k); %Qc给催化剂的热量kj/h
            Gcs(k)=Qc(k)./Cpc/(Trg2(k)-Tst(k)); %Trg1、Trg2待生、再生催化剂温度(沉降器密相温度TD101、再生斜立管温度TD2105) Cpc kj/℃/kg
            %更新gama系数
            Gcd_0(k)=Gcs(k);
            Cv(k)=Gcd_0(k)./(PDrg(k).^0.5); %Cv=func(ivp)
            %gama(k)=Cv(k)/ivp(k);
            %gama_new(k)=0.2.*gama_new(k-1)+(1-0.2).*gama(k);
            Cv_new(k)=(1-0.8).*Cv_new(k-1)+0.8.*Cv(k); %一阶低通滤波
        else
            %gama_new(k)=gama_new(k-1);
            Cv_new(k)=Cv_new(k-1);
        end
        %Gcd(k)=gama_new(k).*ivp(k).*(PDrg(k).^0.5);
        Gcd(k)=Cv_new(k).*(PDrg(k).^0.5);
        Gcf(k)=0.2.*Gcf(k-1)+(1-0.2).*Gcd(k);%一阶低通滤波
        Hr(turn*10+1)=(Gcf(turn*10+1).*Cpc.*(Tra1(turn*10+1)-Trg2(turn*10+1))+(Go(turn*10+1)+Gro(turn*10+1)).*Cpo.*(Tra1(turn*10+1)-220)+Gw(turn*10+1).*Cpw.*(Tra1(turn*10+1)-200)+...
            (Tra2(turn*10+1)-Tra1(turn*10+1)).*((Go(turn*10+1)+Gro(turn*10+1)).*Cpo+Gw(turn*10+1).*Cpw+Gcf(turn*10+1).*Cpc))./(Go(turn*10+1)+Gro(turn*10+1))*(-1);
        Hrd(turn*10+1)=Hr(turn*10+1);
        Hr(k)=((dense_o(1).*Vr(1).*Cpo+Gcf(k).*Vr(1).*dense_o(1)./(Go(k)+Gro(k)).*Cpc+2.*dense_L.*pi.*...
            L(1).*(r_in(1)+Th./2).*Th.*CpL*0.7+Gw(k).*dense_o(1).*Vr(1).*Cpw./(Go(k)+Gro(k))).*(Tra1(k)-Tra1(k-1))*60+...
            (dense_o(2).*Vr(2).*Cpo+Gcf(k).*Vr(2).*dense_o(2)./(Go(k)+Gro(k)).*Cpc+2.*dense_L.*pi.* ...
            L(2).*(r_in(2)+Th./2).*Th.*CpL*0.7+Gw(k).*dense_o(2).*Vr(2).*Cpw./(Go(k)+Gro(k))).*(Tra2(k)-Tra2(k-1))*60+...
            +Gcf(k).*Cpc.*(Tra1(k)-Trg2(k))+(Go(k)+Gro(k)).*Cpo.*(Tra1(k)-220)+Gw(k).*Cpw.*(Tra1(k)-200)+...
            (Tra2(k)-Tra1(k)).*((Go(k)+Gro(k)).*Cpo+Gw(k).*Cpw+Gcf(k).*Cpc))./(Go(k)+Gro(k)).*(-1); %积累+延提升管变化=反应热，三项同号
        Hrd(k)=0.2.*Hrd(k-1)+(1-0.2).*((dense_o(1).*Vr(1).*Cpo+Gcf(k).*Vr(1).*dense_o(1)./(Go(k)+Gro(k)).*Cpc+2.*dense_L.*pi.*...
            L(1).*(r_in(1)+Th./2).*Th.*CpL*0.7+Gw(k).*dense_o(1).*Vr(1).*Cpw./(Go(k)+Gro(k))).*(Tra1(k)-Tra1(k-1))*60+...
            (dense_o(2).*Vr(2).*Cpo+Gcf(k).*Vr(2).*dense_o(2)./(Go(k)+Gro(k)).*Cpc+2.*dense_L.*pi.* ...
            L(2).*(r_in(2)+Th./2).*Th.*CpL*0.7+Gw(k).*dense_o(2).*Vr(2).*Cpw./(Go(k)+Gro(k))).*(Tra2(k)-Tra2(k-1))*60+...
            +Gcf(k).*Cpc.*(Tra1(k)-Trg2(k))+(Go(k)+Gro(k)).*Cpo.*(Tra1(k)-220)+Gw(k).*Cpw.*(Tra1(k)-200)+...
            (Tra2(k)-Tra1(k)).*((Go(k)+Gro(k)).*Cpo+Gw(k).*Cpw+Gcf(k).*Cpc))./(Go(k)+Gro(k)).*(-1);
    end
end
Gcf=Gcf';
Hrd=Hrd';
Gcs=Gcs';
