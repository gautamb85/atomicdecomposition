clear all;
clc;

a = 0.666*cos(2*pi*666*[0:8192*4-1]'/11000);
l = 8192*4;
a1 = cos(2*pi*333*[0:2047]'/11000);
p1 = zeros(l-2048,1);
a1 = [a1;p1];

a2 = 0.533*cos(2*pi*145*[0:4095]'/11000);
p2 = zeros(l-4096,1);
a2  = [a2;p2];


ip = a1 + a2 + a;
%do signal postpadding

lip = length(ip);

window = hanning(8192);

book =  MolecularMP(8192,ip,8192,window);
