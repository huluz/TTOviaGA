%采用GA算法求解输气管道末段非稳态优化问题

% 代码初始化
clc;clear;

%气体相关参数
mol = 17.44;			%相对分子质量
R = 8314;			%通用气体常数
PreC = 46.59E5;		%气体视临界压力
TempC = 198.19;		%气体视临界温度
Temp = 15 + 273.15;		%输气温度
Den_air = 1.209;		%标况下空气密度，kg/m^3
Temp_sta = 20 + 273.15;	%标况温度
Pre_sta = 101325;		%标况压力
Rel_Den = 0.6;			%相对密度
Den_sta = Rel_Den * Den_air;	%标况气体密度
alpha = Den_sta*Temp_sta/Pre_sta/Temp;	%通用气体状态方程参数
beta = 0.257/PreC - 0.533*TempC/PreC/Temp;	%AGA气体状态方程参数

%管段参数
Len = 100E3;			%管道长度，m
Pin = 4.5E6;			%压缩机入口压力
lamda = 0.008;			%摩阻系数
Din = 0.6096;			%管内径，m
Area = 0.25*pi*Din^2;		%管段横截面积
Pe_min = 4.5E6;		%管段出口允许最低压力
Ps_max = 7.3E6;		%管段进口允许最高压力
Qbasic = 65;			%基础流量，m^3/s
Ff = [0.2; 0.15; 0.1; 0.25; 0.35; 0.58; 1.2; 1.3; 1; 0.97; 0.85; 1.65; 2; 1; 0.8; 0.65; ...
	1.15; 1.9; 2.8; 2.2; 1.2; 0.85; 0.5; 0.35];	%小时不均匀用气系数

%计算参数
Time = 24*3600;		%优化时间段，s
Time_Per_Sec = 3600;		%流量边界条件设定时步，s
dt = 15 * 60;			%时步，s
Time_Secs = Time / Time_Per_Sec;	%时间段数
TimeSteps_Total = Time / dt;	%总时步数
TimeSteps_Per_Sec = Time_Per_Sec / dt;	%每个时间段包含的时步数
dx = 10E3;			%空间步长，m
SpaceSteps = Len / dx;		%空间分段数
creat_transfun_re01(SpaceSteps);		%创建状态转移方程
Qs_Min = 5;			%管段进口流量范围, Nm^3/s
Qs_Max = 185;
options = optimset('Display','off');	%求解非线性方程组参数

%模拟参数
C0 = 0.03848;			%稳态模拟公式参数
%a = [4.4889;0.071;-0.009];	%压缩机压头-流量曲线系数
%b = [9.0719;2.2055;-0.0161];	%压缩机效率-流量曲线系数
Pressure_ini = zeros(SpaceSteps+1,1);		%沿线压力分布
MassFlux_ini = zeros(SpaceSteps+1,1);	%沿线质量流量分布

%边界条件
Qe = zeros(TimeSteps_Total,1);		%终点流量
for i = 1:Time_Secs			%根据时间点上的值设定整个时间段的流量
	%Qe(TimeSteps_Per_Sec*(i-1)+1:TimeSteps_Per_Sec*i) = Qbasic*Ff(i)*ones(TimeSteps_Per_Sec,1);
	for j = 1:TimeSteps_Per_Sec
		if i == 1
			Qe(TimeSteps_Per_Sec*(i-1)+j) = Qbasic*((Ff(i)-Ff(Time_Secs))*j/TimeSteps_Per_Sec + Ff(Time_Secs));
		else
			Qe(TimeSteps_Per_Sec*(i-1)+j) = Qbasic*((Ff(i)-Ff(i-1))*j/TimeSteps_Per_Sec + Ff(i-1));
		end
	end
end
Mse = Den_sta * Qe/Area;	%终点质量流量密度

%GA算法初始化-设置初始状态
%稳态模拟
tl = Len;			%管段长度
Ple = Pe_min;			%起点压力
Pressure_ini(SpaceSteps+1) = Ple;%沿线压力记录
i = SpaceSteps;
while tl>0 			%稳态模拟
	z = 1 + beta*Ple;	%压缩因子
	Pls = Ple^2 + lamda*z*Rel_Den*Temp*dx*Qe(TimeSteps_Total)^2/C0^2/Din^5;
	Pls = Pls^0.5;
	Pressure_ini(i) = Pls;
	i = i - 1;
	tl = tl - dx;
	Ple = Pls;
end
MassFlux_ini = (Den_sta*Qe(TimeSteps_Total)/Area)*ones(SpaceSteps+1,1);	%构造初始条件-质量流量密度
Pressure_ini = Pressure_ini';		%压力
%figure(1);
%plot(Pressure_ini);
%title('Pressure_ini');
%figure(2);
%plot(MassFlux_ini);
%title('MassFlux_ini');

%GA算法求解优化问题

%初始化参数
NP = 10;			%种群规模
Max_Gen = 10;			%最大遗传代数
PrePCoe = 1e-2;		%管段末段压力罚因子
%MsPCoe = 1e2;		%管段沿线流量罚因子
bits = 4;			%表示单个时段决策值需要的基因位数
dq = (Qs_Max - Qs_Min)/(power(2, bits) - 1);	%可行域离散精度
Total_bits = 4 * Time_Secs;	%基因长度
PC = 0.9;			%交叉概率
PM = 0.1;			%编译概率
gen = 1;			%种群代数
ades = 0.96;			%适应度函数标定参数
delta = 0;			%暂时设为零，以观察遗传算法的性状
MaxConsum = 1e7;		%前一种群最大函数值
GensPool = zeros(NP, Total_bits);	%基因池
FitsWheel = zeros(NP, 1);	%转轮法参数
MinObj = zeros(Max_Gen,1);	%记录每轮的最优目标函数值
OptQs = zeros(Max_Gen,Time_Secs);	%最优方案
OptRecs = zeros(Max_Gen,Total_bits);	%最优记录
AverageObjV = zeros(Max_Gen,1);	%记录每轮的平均适应度值

%开启并行计算池
startmatlabpool;

%产生初始种群
Chromes = round(rand(NP, Total_bits));

while gen <= Max_Gen
	%计算适应度函数值
	fprintf('Generation: %d\n', gen);
	tic;
	ObjV = zeros(NP, 1);
%	fprintf('%s\n', 'Computing Fitness...');
	parfor nn = 1:NP
%		fprintf('%d\n', nn);
		Qs_Temp = zeros(Time_Secs, 1);	%入口流量-整点处
		for kk = 1:Time_Secs
			dots = 0;			%离散点位置
			for mm = 1:bits
				dots = dots + Chromes(nn, bits*(kk - 1)+mm)*power(2, 4 - mm);
			end
			Qs_Temp(kk) = Qs_Min + dq*dots;
		end
		Qs = zeros(TimeSteps_Total,1);	%终点流量
		for oo = 1:Time_Secs			%根据时间点上的值设定整个时间段的流量
			for ee = 1:TimeSteps_Per_Sec
				if oo == 1
					Qs(TimeSteps_Per_Sec*(oo-1)+ee) = (Qs_Temp(oo)-Qs_Temp(Time_Secs))*ee/TimeSteps_Per_Sec + Qs_Temp(Time_Secs);
				else
					Qs(TimeSteps_Per_Sec*(oo-1)+ee) = (Qs_Temp(oo)-Qs_Temp(oo-1))*ee/TimeSteps_Per_Sec + Qs_Temp(oo-1);
				end
			end
		end
		Mss = (Den_sta/Area) * Qs;

		%计算适应度函数
		Pressure = Pressure_ini; MassFlux = MassFlux_ini;
		ComConsum = 0;			%压缩机功率
		for i = 1:TimeSteps_Total
			tf = @(x)transfun_re01(x,dt,dx,alpha,beta,lamda,Din,Pressure,MassFlux,Mss(i),Mse(i));	%构造方程
			x0 = zeros(2*SpaceSteps,1);	%准备初值
			x0(1) = Pressure(1);
			for j = 2:SpaceSteps
				x0(2*j-2) = Pressure(j);
				x0(2*j-1) = MassFlux(j);
			end
			x0(2*SpaceSteps) = Pressure(SpaceSteps+1);
			%options = optimset('MaxFunEvals',1E4,'MaxIter',1E4,'Display','iter-detailed');
			results = fsolve(tf,x0,options);		%计算
	%		results = (1-gama)*x0+gama*results;	%引入阻尼系数，防止计算结果震荡
			Pressure(1) = results(1);	%归档计算结果
			for j = 2: SpaceSteps
				Pressure(j) = results(2*j-2);
				MassFlux(j) = results(2*j-1);
			end
			Pressure(SpaceSteps+1) = results(2*SpaceSteps);
			Quan_Temp = (0.328*Area/Den_sta)*MassFlux(1);
			Sec_Com_Consum = dt*Quan_Temp*(2.682*(Pressure(1)/Pin)^0.217 - 2.658);			%该时步压缩机功率
			if Sec_Com_Consum < 0
				Sec_Com_Consum = 2e4;	%处理压缩机功率为负的情况
			end
			ComConsum = ComConsum + Sec_Com_Consum;	%计算压缩机功率
%			fprintf('%s%f\n', 'Compressor: ',dt*Quan_Temp*(2.682*(Pressure(1)/Pin)^0.217 - 2.658));
			if Pressure(SpaceSteps+1) < Pe_min
				ComConsum = ComConsum + PrePCoe*abs(Pressure(SpaceSteps+1) - Pe_min);	%引入罚函数部分
%				fprintf('%s%f\n', 'Pressure Punishment: ',PrePCoe*abs(Pressure(SpaceSteps+1) - Pe_min));
			end
%			if min(MassFlux) < 0
%				inc = MassFlux < 0;	%质量流量小与零的索引
%				ComConsum = ComConsum + MsPCoe * abs(inc' * MassFlux);
%				fprintf('%s%f\t', 'Quantity Punishment: ',MsPCoe * abs(inc' * MassFlux));
%				fprintf('%d\n', inc' * MassFlux);
%			end
%			fprintf('\n');
			if i < TimeSteps_Total
				MassFlux(1) = Mss(i+1);	%引入边界条件
				MassFlux(SpaceSteps+1) = Mse(i+1);
			end
		end
		ObjV(nn) = -ComConsum + MaxConsum + power(ades, gen)*delta;
		if ObjV(nn) < 0 
			ObjV(nn) = 0;				%剔除一些严重不可行方案
		end
	end
	MinObj(gen) = -1*(max(ObjV) - MaxConsum - power(ades, gen)*delta);	%最优目标函数值
	[Temp, MinRecNum] = max(ObjV);
	fprintf('%s%e\n', 'Minimum Object Function Value: ', MinObj(gen));
	OptRecs(gen,:) = Chromes(MinRecNum,:);	%最优记录
	Qs_Opt_Gen = zeros(Time_Secs, 1);	%入口流量-整点处
	for kk = 1:Time_Secs
		dots = 0;			%离散点位置
		for mm = 1:bits
			dots = dots + power(2, Chromes(MinRecNum, bits*(kk - 1)+mm));
		end
		Qs_Opt_Gen(kk) = Qs_Min + dq*dots;
	end
	OptQs(gen,:) = Qs_Opt_Gen;			%记录最优方案
	AverageObjV(gen) = mean(ObjV);		%平均目标函数值
	fprintf('%s%e\n', 'Average Fitness: ', AverageObjV(gen));
%	if mod(gen,10) == 0
%		plot(MinObj(1:gen));						%计算过程可视化
%		plot(AverageFit(1:gen));
%	end

	%生成基因池-转轮法
%	fprintf('%s\n', 'Generating Genpool...');
	TotalFits = sum(ObjV);
	FitsPro = ObjV / TotalFits;
	for pp = 1:NP
		FitsWheel(pp) = sum(FitsPro(1:pp));
	end
	for pp = 1:NP
		pro = rand;
		pos = ceil(NP / 2);
		pos_min = 0;
		pos_max = NP;
		while ~((pos == 1 && FitsWheel(pos)>= pro) || (FitsWheel(pos) >= pro && FitsWheel(pos - 1) <= pro))
			if FitsWheel(pos) > pro
				pos_max = pos;
				pos = ceil((pos_max + pos_min)/2);
			else
				pos_min = pos;
				pos = ceil((pos_max+pos_min)/2);
			end
		end
		GensPool(pp,:) = Chromes(pos,:);
	end

	%交叉
%	fprintf('%s\n', Implementing Crossover...');
	Chromes = zeros(NP, Total_bits);			%清空
	for pp = 1:NP/2
		gen1 = GensPool(round(1 + (NP - 1)*rand),:);	%抽取基因
		gen2 = GensPool(round(1 + (NP - 1)*rand),:);
		CPro = rand;					%确定是否进行交叉
		if CPro < PC
			pos = round(1 + (Total_bits - 2)*rand);	%确定交叉节点
			gen11 = gen1(1:pos);
			gen21 = gen2(1:pos);
			gen1 = [gen21 gen1(pos+1:end)];
			gen2 = [gen11 gen2(pos+1:end)];
			Chromes(2*pp - 1,:) = gen1;
			Chromes(2*pp,:) = gen2;
		else
			Chromes(2*pp - 1,:) = gen1;
			Chromes(2*pp,:) = gen2;
		end
	end
	
	%变易
%	fprintf('%s\n', 'Implementing Mutation...');
	for pp = 1:NP
		MPro = rand;
		if MPro < PM
			pos = round(1 + (Total_bits - 1)*rand);	%确定变易位置
			Chromes(pp,pos) = ~Chromes(pp,pos);
		end
	end

	%更新参数，准备下一轮计算
	MaxConsum_Temp = max(-1*(ObjV - MaxConsum - power(ades, gen)*delta));
	if MaxConsum_Temp < MaxConsum
		MaxConsum = (MaxConsum_Temp + MaxConsum) / 2;		%最大函数值，用于标定适应度函数
	end
	gen = gen + 1;			%增加遗传代数计数

	%计算过程可视化
	t1 = toc;
	fprintf('Time: %f s\n', t1);
end