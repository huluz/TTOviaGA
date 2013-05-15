%批量计算优化方案脚本

%上升情况
clc;clear;
load Increasing_InitialData;
TTOviaGA_re01(ff,time,dt,qsmin,qsmax,qbasic,fn,bits,np,maxgen,pc,pm,prepcoe);

%下降情况
load Decreasing_InitialData;
TTOviaGA_re01(ff,time,dt,qsmin,qsmax,qbasic,fn,bits,np,maxgen,pc,pm,prepcoe);

%波动情况
load Ramping_InitialData;
TTOviaGA_re01(ff,time,dt,qsmin,qsmax,qbasic,fn,bits,np,maxgen,pc,pm,prepcoe);

%长时间缓波动情况
load Slow_Ramping_24H_InitialData;
TTOviaGA_re01(ff,time,dt,qsmin,qsmax,qbasic,fn,bits,np,maxgen,pc,pm,prepcoe);

%长时间剧烈波动情况
load Severe_Ramping_24H_InitialData;
TTOviaGA_re01(ff,time,dt,qsmin,qsmax,qbasic,fn,bits,np,maxgen,pc,pm,prepcoe);